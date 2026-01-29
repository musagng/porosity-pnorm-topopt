%--------------------------------------------------------------------------
% BONE_REMODELING
%
% Bone-remodeling-inspired density adaptation algorithm for generating
% an initial density field with embedded stress-path information.
%
% The algorithm iteratively adapts element densities based on equivalent
% strain energy, targeting a prescribed global volume fraction.
%
% OUTPUT:
%   rho_BR : Final density field to be used as initialization for
%            p-norm constrained MMA optimization
%
% INPUT:
%   volfrac_target : Target volume fraction (e.g., 0.4)
%   nelx, nely     : Number of elements in x and y directions
%
%--------------------------------------------------------------------------
function rho_BR = bone_remodeling(volfrac_target, nelx, nely)

%% -------------------- PARAMETERS ----------------------------------------
penal = 3;
rmin  = 3;

E0   = 1;
Emin = 1e-9;
nu   = 0.3;

x_min = 0.001;
x_max = 1;

%% -------------------- FINITE ELEMENT SETUP -------------------------------
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];

KE = 1/(1-nu^2)/24 * ...
    ([A11 A12; A12' A11] + nu*[B11 B12; B12' B11]);

nodenrs = reshape(1:(1+nely)*(1+nelx),1+nelx,1+nely);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);

edofMat = repmat(edofVec,1,8) + ...
    repmat([0 1 2*nelx+[2 3 0 1] -2 -1],nelx*nely,1);

iK = reshape(kron(edofMat,ones(8,1))',[],1);
jK = reshape(kron(edofMat,ones(1,8))',[],1);

%% -------------------- LOADS AND SUPPORTS --------------------------------
F = zeros(2*(nely+1)*(nelx+1),1);
load_node = nodenrs(floor((nelx+2)/2), end);
F(load_node*2) = -1;                 % Vertical load
F = sparse(F);

left_nodes = nodenrs(:,1);
fixeddofs = unique([2*left_nodes(:)-1; 2*left_nodes(:)]);
alldofs   = 1:2*(nely+1)*(nelx+1);
freedofs = setdiff(alldofs, fixeddofs);

U = zeros(2*(nely+1)*(nelx+1),1);

%% -------------------- DENSITY FILTER ------------------------------------
iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));

k = 0;
for i1 = 1:nely
    for j1 = 1:nelx
        e1 = (i1-1)*nelx + j1;
        for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nely)
            for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nelx)
                e2 = (i2-1)*nelx + j2;
                k = k + 1;
                iH(k) = e1;
                jH(k) = e2;
                sH(k) = max(0, rmin - sqrt((i1-i2)^2 + (j1-j2)^2));
            end
        end
    end
end

H  = sparse(iH,jH,sH);
Hs = sum(H,2);
%% -------------------- INITIALIZATION ------------------------------------
x     = repmat(volfrac_target, nelx, nely);
xPhys = x;

delta1 = 0.02;
volfrac_d = 0.001;

%% Initial strain reference
sK = reshape(KE(:)*(Emin + xPhys(:)'.^penal*(E0-Emin)),[],1);
K  = sparse(iK,jK,sK); 
K  = (K+K')/2;

U(freedofs) = K(freedofs,freedofs)\F(freedofs);
ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nelx,nely);

strain_ref = mean(sqrt(2*ce(:))) * sqrt(volfrac_target);

%% -------------------- ITERATION LOOP ------------------------------------
loop = 0;
ch   = 0.001;

while ch < 99
      loop = loop + 1;

      % FEA
      sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
      K = sparse(iK,jK,sK); 
      K = (K+K')/2;
      U(freedofs) = K(freedofs,freedofs)\F(freedofs);

      % Equivalent strain
      ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nelx,nely);
      c = sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce));
      sigma_eq = reshape((sqrt(2.*ce)),nelx*nely,1);

      % Density update
      delta = (delta1 * (sigma_eq - strain_ref));
      delta = max(min(delta,(strain_ref)),-strain_ref);
      delta(abs(sigma_eq - strain_ref) < 0.15*strain_ref) = 0;
      x_new = x(:) + delta;
      x_new = max(min(x_new,x_max),x_min);
      x_new = reshape(x_new, nelx, nely);
      
      % Filtering
      xPhys(:) = (H*x_new(:))./Hs;
      xPhys = max(min(xPhys,x_max),x_min);
      xPhys = reshape(xPhys, nelx, nely);
      
      % Convergence check - Reeference Strain Update
      delta_x = abs(x_new(:)-x(:));
      ch = (length(find(delta_x == 0)) / (nelx*nely))*100;
      x = x_new;
      volfrac = mean(xPhys(:));
      strain_d = 2 * (abs(volfrac - volfrac_target));
      strain_d = max(min(strain_d,0.04),0.01);
      if loop < 100
          if ch > 50 && volfrac > volfrac_target * (1 + volfrac_d)
              strain_ref = strain_ref + strain_ref * strain_d;
          elseif ch > 60 && volfrac < volfrac_target * (1 - volfrac_d)
              strain_ref = strain_ref - strain_ref * strain_d;
          end
      else
          if ch > 50 && volfrac > volfrac_target * (1 + volfrac_d)
              strain_ref = strain_ref + strain_ref * strain_d;
          elseif ch > 90 && volfrac < volfrac_target * (1 - volfrac_d)
              strain_ref = strain_ref - strain_ref * strain_d;
          end
          if ch > 95
              delta1 = 0.01;
          elseif ch > 97
              delta1 = 0.005;
          end
      end

      fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,c, ...
      volfrac,ch);
      colormap(gray); imagesc(1-xPhys); clim([0 1]); axis equal; axis off; drawnow;
end
%% -------------------- OUTPUT --------------------------------------------
rho_BR = xPhys;
end