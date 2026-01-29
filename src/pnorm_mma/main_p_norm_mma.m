%% =========================================================
%  main_p_norm_mma.m
%
%  Porosity-controlled topology optimization
%  using p-norm constraints and MMA
%
%  Author  : Musa Güngörürler
%  GitHub  : porosity-pnorm-topopt
% ==========================================================

function xPhys = main_p_norm_mma(rho_BR)
    %% PROBLEM SETUP
    nely     = size(rho_BR,1);
    nelx     = size(rho_BR,2);
    volfrac  = mean(rho_BR(:));
    nloop    = 2000;

    mdof     = [1 2];      % Active constraints
    rho_tol  = 0.025;
    
    penal = 3;             % SIMP penalty
    p     = 12;            % p-norm order
    
    rmin  = 2;             % Density filter radius
    move  = 0.025;         % MMA move limit
    n = nelx*nely;
    
    beta = 1;              % Heaviside continuation
    eta  = 0.5;            % Projection threshold
    
    %% MATERIAL PROPERTIES
    E0 = 1;
    Emin = 1e-9;
    nu = 0.3;
    
    %% FINITE ELEMENT PREPARATION
    A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
    A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
    B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
    B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
    
    KE = 1/(1-nu^2)/24 * ([A11 A12;A12' A11] + nu*[B11 B12;B12' B11]);
    
    nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
    edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
    edofMat = repmat(edofVec,1,8) + repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
    
    iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
    jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
    
    %%  LOADS AND BOUNDARY CONDITIONS (MBB-type)
    Fsparse = sparse(2*(nely+1)*(nelx+1),1);
    Fsparse(2*(nely+1)*nelx + nely + 2) = -1;
    
    fixeddofs = union(1:2*(nely+1),1);
    alldofs   = 1:2*(nely+1)*(nelx+1);
    freedofs = setdiff(alldofs,fixeddofs);
    
    U = zeros(2*(nely+1)*(nelx+1),1);
    
    %% DENSITY FILTER (Helmholtz-type)
    iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
    jH = ones(size(iH));
    sH = zeros(size(iH));
    k  = 0;
    
    for i1 = 1:nelx
        for j1 = 1:nely
            e1 = (i1-1)*nely+j1;
            for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
                for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
                    e2 = (i2-1)*nely+j2;
                    k = k+1;
                    iH(k) = e1;
                    jH(k) = e2;
                    sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
                end
            end
        end
    end
    
    H  = sparse(iH,jH,sH);
    Hs = sum(H,2);

    %%  REGIONAL PARTITIONING AND REFERENCE DENSITY DEFINITION

    Bx = 10;   % x yönünde blok genişliği
    By = 10;   % y yönünde blok yüksekliği
    
    nRx = nelx / Bx;
    nRy = nely / By;
    nReg = nRx * nRy;
    
    regionElems = cell(nReg,1);
    reg = 0;
    
    for ix = 1:nRx
        for iy = 1:nRy
            reg = reg + 1;
            ex = (ix-1)*Bx + (1:Bx);
            ey = (iy-1)*By + (1:By);
            [EX,EY] = meshgrid(ex,ey);
            regionElems{reg} = sub2ind([nely nelx],EY(:),EX(:));
        end
    end
    rhoRef_reg = zeros(nReg,1);
    for r = 1:nReg
        rhoRef_reg(r) = mean(rho_BR(regionElems{r}));
    end
    
    %%  INITIALIZE ITERATION
    x = repmat(volfrac,nely,nelx);
    xTilde = x;
    xPhys = (tanh(beta*eta) + tanh(beta*(xTilde-eta))) / (tanh(beta*eta) + tanh(beta*(1-eta)));
    xold1 = reshape(x,[nely*nelx,1]);
    xold2 = reshape(x,[nely*nelx,1]);
    low = 0;
    upp = 0;
    %rho_BR = (TF'*(LF'\(LF\(TF*rho_BR(:)))));
    %rho_BR = reshape(rho_BR,nely,nelx);
    loopbeta = 0;
    loop = 0;
    change = 1;
    c_hist = zeros(nloop,1);        % compliance
    
    %% START ITERATION
    while change > 0.0001 && loop < nloop
        loopbeta = loopbeta+1;
        loop = loop+1;
    
        %%  FEA
        sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
        K = sparse(iK,jK,sK); K = (K+K')/2;
        U(freedofs) = K(freedofs,freedofs)\Fsparse(freedofs);
    
        %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
        ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
        c = sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce));
        dc = -penal*(E0-Emin)*xPhys.^(penal-1).*ce;
    
        %% REGIONAL POROSITY CONSTRAINT (p-NORM AGGREGATION)
        rhoBar = zeros(nReg,1);
        for r = 1:nReg
            rhoBar(r) = mean(xPhys(regionElems{r}));
        end
        f_por = (sum(abs(rhoBar - rhoRef_reg).^p)/nReg)^(1/p) - rho_tol;
        dfdx_por = zeros(nelx*nely,1);
    
        normp = (sum(abs(rhoBar - rhoRef_reg).^p)/nReg)^(1/p);
        
        for r = 1:nReg
            elems = regionElems{r};
            coeff = (abs(rhoBar(r)-rhoRef_reg(r))^(p-1)) ...
                    / (nReg * normp^(p-1) + 1e-12);
            sgn = sign(rhoBar(r) - rhoRef_reg(r));
            dfdx_por(elems) = coeff * sgn / length(elems);
        end
        dfdx_por = H*(dfdx_por./Hs);
    
        %% FILTERING/MODIFICATION OF SENSITIVITIES
        dx = beta * (1-tanh(beta*(xTilde-eta)).*tanh(beta*(xTilde-eta))) / (tanh(beta*eta) + tanh(beta*(1-eta)));
        dc(:) = H*(dc(:).*dx(:)./Hs);
    
        %% METHOD OF MOVING ASYMPTOTES (MMA) UPDATE
        m = 2; 
        
        df0dx = reshape(dc,[nelx*nely,1]);
        dfdx = zeros(2,nelx*nely);
        dfdx(1,:) = ones(1,nelx*nely)/(nelx*nely*volfrac);
        dfdx(2,:) = dfdx_por';
    
        iter = loopbeta;
        xval = reshape(x,[nelx*nely,1]);
        xmin=max(0.0,xval-move);
        xmax=min(1,xval+move);
    
        f0val = c;
        fval = zeros(2,1);
        fval(1) = sum(xPhys(:))/(nelx*nely*volfrac) - 1;
        fval(2) = f_por;
    
        a0 = 1;
        a = zeros(m,1);     
        c_ = ones(m,1)*1000;
        d = zeros(m,1);
        [xmma,~,~,~,~,~,~,~,~,low,upp] = ...
            mmasub(m,n,iter,xval,xmin,xmax,xold1,xold2,...
            f0val,df0dx,fval(mdof),dfdx(mdof,:),low,upp,a0,a,c_,d);
    
        %% DESIGN UPDATE & VISUALIZATION
        xnew = reshape(xmma,[nely,nelx]);
        xold2 = xold1;
        xold1 = xval;
        
        xTilde(:) = (H*xnew(:))./Hs;
        xPhys = (tanh(beta*eta) + tanh(beta*(xTilde-eta))) / (tanh(beta*eta) + tanh(beta*(1-eta)));
    
        change = max(abs(xnew(:)-x(:)));
        x = xnew;
    
        if beta < 60 && (loopbeta >= 40 || change < 0.002)
            beta = 1.02 * beta;
            loopbeta = 0;
        end
        if loop < 200 && change < 0.002
            move = 0.05;
        end
    
        % Store current values
        c_hist(loop,1) = c;
    
        %% PLOT DENSITIES
        disp([' It.: ' sprintf('%4i',loop) ' Obj.: ' sprintf('%10.4f',c) ...
            ' Vol.: ' sprintf('%6.3f',sum(sum(xPhys))/(nelx*nely)) ...
            ' Ch.: ' sprintf('%6.3f',change) ...
            ' Cons.: ' sprintf('%6.3f',fval)]);
        figure(1);
        set(1, 'Position', [100, 450, 540, min(100+540*nely/nelx,540)]);
        colormap(gray); imagesc(1-xPhys); clim([0 1]); axis equal; axis off; drawnow;
    end
end