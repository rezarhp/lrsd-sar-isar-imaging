% MAIN_CONV_LRSD
% Conventional (CGS-based) LRSD SAR/ISAR imaging and decomposition demo.
%
% Please cite this article if you use this code:
% H. R. Hashempour, M. Moradikia, H. Bastami, A. Abdelhadi and M. Soltanalian,
% "Fast and Robust LRSD-Based SAR/ISAR Imaging and Decomposition,"
% IEEE Transactions on Geoscience and Remote Sensing, vol. 60, pp. 1-13, 2022,
% Art no. 5227413, doi: 10.1109/TGRS.2022.3172018.
%
% Dependencies:
%   - Extract_ROI.m
%   - Patch_Operator.m
%   - Inverse_Patch_Operator.m
%   - USFO.m
%   - sar_cmap.m
% Data:
%   - images/sub_image_from_SARPER.png
%
% (c) Hamid Reza Hashempour

clc;
clear; %#ok<CLSCR>
close all;

%% Settings
Nm = 256;
imgPath = fullfile("images", "sub_image_from_SARPER.png");

% Dependency checks (friendly errors)
assert(exist("Extract_ROI","file")==2, "Missing Extract_ROI.m (add repo to path).");
assert(exist("Patch_Operator","file")==2, "Missing Patch_Operator.m (add repo to path).");
assert(exist("Inverse_Patch_Operator","file")==2, "Missing Inverse_Patch_Operator.m (add repo to path).");
assert(exist("USFO","file")==2, "Missing USFO.m (add repo to path).");
assert(exist("sar_cmap","file")==2, "Missing sar_cmap.m (add repo to path).");
assert(isfile(imgPath), "Missing image file: %s", imgPath);

%% Load original scene
background1 = im2double(rgb2gray(imread(imgPath))); % stationary targets
background_cropped = Extract_ROI(background1, Nm, Nm);
original_scene = background_cropped;

fconv = original_scene(:);
Nx = size(original_scene, 1);
Ny = size(original_scene, 2);

%% Patch settings
Mx = 32; % patch width
My = 32; % patch height
Sx = 10; % stride in x
Sy = 10; % stride in y

rt = 1:Sy:Ny-My+1; if rt(end) ~= Ny-My+1, rt = [rt, Ny-My+1]; end %#ok<AGROW>
ct = 1:Sx:Nx-Mx+1; if ct(end) ~= Nx-Mx+1, ct = [ct, Nx-Mx+1]; end %#ok<AGROW>

Ns = numel(rt) * numel(ct);
M  = Mx * My;

%% Algorithm parameters
mu1      = 1;        % initial mu1
eta1     = 1;        % decreasing factor for mu1
mu1_min  = 1e-4;     % minimum mu1
rho1     = 0.95;     % decomposition parameter
lambda1  = rho1 * (max(M, Ns))^(-1/2);

delta1   = 1;        % penalty parameters
delta2   = 1;

maxIter  = 100;
tol      = 1e-6;
quiet    = false;

%% Sampling mask
undersampling = 1; % 1 => full sampling
if undersampling == 1
    ind1 = 1:Nx;
    ind2 = 1:Ny;
else
    ind1 = randi(Nx, round(undersampling * Nx), 1);
    ind2 = randi(Ny, round(undersampling * Ny), 1);
end

mask = zeros(Nx, Ny);
mask(ind1, ind2) = 1;
fprintf("Mask nnz = %d\n", nnz(mask));

%% Operators
R    = @(z) Patch_Operator(z, Nx, Ny, Mx, My, Ns, rt, ct);
Rinv = @(z) Inverse_Patch_Operator(z, Nx, Ny, Mx, My, Ns);

UF  = @(Z) USFO(Z, mask, 1);
UFh = @(Z) USFO(Z, mask, 2);

Betta  = ones(Nx * Ny, 1);
P_tild = ones(nnz(mask), 1);

E   = @(z) P_tild .* (UF((Betta) .* Rinv(z)));
Eh  = @(z) R(conj(Betta) .* UFh(P_tild .* z));
EhE = @(z) Eh(E(z));

%% Initialization
F = reshape(abs(R(fconv)), [M, Ns]);
y = E(F(:));

recidft = reshape(Rinv(Eh(y)), [Nx, Ny]);
plot_db_image(recidft, "2D-FT");

B  = F;
S  = zeros(M, Ns);
Z1 = zeros(M, Ns);
Z2 = zeros(M, Ns);

out.stopCriterion = zeros(maxIter, 1);
out.F = cell(maxIter, 1);
out.B = cell(maxIter, 1);

%% Main loop (conventional CGS)
tic;
converged = false;
iter = 0;

while ~converged && iter < maxIter

    Fp = F;

    lambdaB = mu1;
    lambdaS = lambda1 * mu1;

    % LRSD updates
    L = svt(F - S + Z1 / delta1, lambdaB / delta1);
    A = shrink(F - B + Z2 / delta2, lambdaS / delta1);

    bl = Eh(y) + delta1 * L(:) - Z1(:) - EhE(S);

    [temp, ~] = cgs(@(z) EhE(z(:)) + delta1 * z(:), bl, 1e-5, 100);
    B = reshape(temp, [M, Ns]);

    bs = Eh(y) + delta2 * A(:) - Z2(:) - EhE(B(:));

    [temp1, ~] = cgs(@(z) EhE(z(:)) + delta2 * z(:), bs, 1e-5, 100);
    S = reshape(temp1, [M, Ns]);

    Z1 = Z1 + delta1 * (B - L);
    Z2 = Z2 + delta2 * (S - A);

    F = B + S;

    % Update mu1
    mu1 = max(eta1 * mu1, mu1_min);

    % Bookkeeping
    iter = iter + 1;
    out.stopCriterion(iter) = norm(F(:) - Fp(:)) / max(norm(Fp(:)), eps);
    out.F{iter} = F;
    out.B{iter} = B;

    if ~quiet
        fprintf("iter=%3d | stop=%.3e\n", iter, out.stopCriterion(iter));
    end

    if out.stopCriterion(iter) <= tol && iter > 1
        converged = true;
        if ~quiet
            fprintf("Converged: tol=%.1e reached.\n", tol);
        end
    end
end

elapsedtime = toc;
fprintf("Total runtime: %.3f s\n", elapsedtime);

%% Reconstruction and plots
frec = Rinv(F);
brec = Rinv(B);
srec = Rinv(S);

fmat = abs(reshape(frec, [Ny, Nx]));
bmat = abs(reshape(brec, [Ny, Nx]));
smat = abs(reshape(srec, [Ny, Nx]));

plot_db_image(fmat, sprintf("f reconstructed (iter %d)", iter));
plot_db_image(bmat, sprintf("b reconstructed (iter %d)", iter));
plot_db_image(smat, sprintf("s reconstructed (iter %d)", iter));

figure(6);
plot(out.stopCriterion(1:iter), 'LineWidth', 1.2);
grid on;
title('Stop criterion vs iterations');

%% Save results (recommended for public repo)
if ~exist("results", "dir"), mkdir("results"); end
save(fullfile("results", "run_summary_conv.mat"), "out", "iter", "elapsedtime", ...
     "Mx", "My", "Sx", "Sy", "Nm", "lambda1", "mu1_min", "rho1", "delta1", "delta2");


%% ================= Local helper functions =================

function plot_db_image(x, ttl)
x = abs(x);
xdb = 20 * log10(max(x, 1e-12));
xdb(xdb < -30) = -30;

figure;
colormap(sar_cmap);
imagesc(xdb, [-30 0]);
axis image off;
colorbar;
title(ttl);
drawnow;
end

function B = svt(A, gamma)
[U, S, V] = svd(A, 0);
Z = diag(shrink(diag(S), gamma));
B = U * Z * V';
end

function A = shrink(B, gamma)
A = sign(B) .* max(abs(B) - gamma, 0);
end
