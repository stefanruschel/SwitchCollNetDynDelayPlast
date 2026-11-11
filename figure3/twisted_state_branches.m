%% twisted_state_branches.m
% Branches of m-twisted states of the ring of N = 11 nodes with constant delays (epsilon = 0) and their stability,
% drawn in grey in Fig. 3(a),(b) of the Letter.
%
% Requirements: MATLAB and DDE-Biftool 3.1 (https://github.com/DDE-BifTool/DDE-Biftool); set DDEBIFTOOL below.
% Run time: about one hour (stability of 6 branches x 250 points x 11 Fourier modes).
% Output:   data/twisted_state_branches.mat with the variables
%             theta     1 x 250, the parameter theta = Omega * tau along the branches, from 0 to 4 pi
%             Omega     cell array, Omega{m} = collective frequency along the branch of psi_m
%             tau1      cell array, tau1{m}  = delay per unit distance, tau = 1/v, along the branch of psi_m
%             Npos      cell array, Npos{m}  = number of unstable eigenvalues along the branch of psi_m
%           where m = 1..5 is the twist number and psi_0 is stored under m = N = 11 (MATLAB indices start at 1).
%           export_figure3.jl reads this file.
%
% Branches. The m-twisted state theta_j(t) = Omega t + 2 pi m j / N with all delays tau_j = d_j tau (d_j = ring
% distance of neighbour j) exists if
%     Omega = omega + sigma * sum_{j=1}^{N-1} W_j h(2 pi m j / N - d_j Omega tau),     h(x) = -sin(x - a) + r sin(2x).
% Writing theta = Omega tau, the right-hand side is a function of theta alone, which gives each branch in
% parametric form: Omega(theta) from the formula, tau(theta) = theta / Omega(theta).
%
% Stability. Perturbations of the twisted state decouple into Fourier modes l = 0..N-1; mode l obeys the scalar DDE
%     x'(t) = -kappa x(t) + sum_{j=1}^{N-1} sigma W_j h'(2 pi m j / N - d_j theta) exp(2 pi i l j / N) x(t - d_j tau),
%     kappa = sigma sum_j W_j h'(2 pi m j / N - d_j theta).
% Its rightmost eigenvalues are computed with DDE-Biftool's dde_stst_eig_cheb; an eigenvalue counts as unstable if
% its real part exceeds 1e-3 (mode l = 0 always has the trivial eigenvalue 0 from the phase shift symmetry).

clear; close all;
DDEBIFTOOL = fullfile('..', '..', 'DDE-Biftool');     % <-- folder that contains ddebiftool/, ddebiftool_utilities/, ...
addpath(fullfile(DDEBIFTOOL, 'ddebiftool'), fullfile(DDEBIFTOOL, 'ddebiftool_extra_psol'), ...
        fullfile(DDEBIFTOOL, 'ddebiftool_extra_nmfm'), fullfile(DDEBIFTOOL, 'ddebiftool_extra_symbolic'), ...
        fullfile(DDEBIFTOOL, 'ddebiftool_utilities'));
if exist('dde_stst_eig_cheb', 'file') ~= 2
    error('DDE-Biftool not found: set DDEBIFTOOL at the top of this script');
end

%% Parameters (caption of Fig. 3)
N     = 11;
omega = 1.0;
sigma = 0.1;
a     = 3./2*pi;      % h(x) = -sin(x - 3 pi/2) = sin(x - pi/2)
r     = 0.;           % no second harmonic

h      = @(x) -sin(x - a) + r.*sin(2.*x);
hprime = @(x) -cos(x - a) + 2.*r.*cos(2.*x);

theta = linspace(0, 4*pi, 250);                        % theta = Omega * tau

%% Ring coupling, seen from node 0: neighbours j = 1..N-1 at distance min(j, N - j), weights exp(-distance)
neighbour = 1:N-1;
distance  = min(neighbour, N - neighbour);
weight    = exp(-distance);
weight    = weight / sum(weight);                      % unit row sum (no self-coupling)
twist_phase = 2*pi*neighbour./N;                       % phase lag 2 pi j / N per unit twist

%% Branches
Omega = cell(N, 1);
tau1  = cell(N, 1);
for m = [N, 1:floor(N/2)]                              % m = N stands for psi_0
    Omega{m} = omega + sigma * sum( weight.' .* h( m*twist_phase.' - distance.' .* theta ), 1 );
    tau1{m}  = theta ./ Omega{m};
end

%% Stability: number of unstable eigenvalues at every branch point
Npos = cell(N, 1);
for m = [N, 1:floor(N/2)]
    Npos{m} = 0.*theta;
    first = 1 + (tau1{m}(1) == 0);                     % skip tau = 0 (no delay), filled in below
    for k = first:length(theta)
        slopes = sigma .* weight .* hprime(m*twist_phase - distance*theta(k));
        kappa  = sigma*sum(weight .* hprime(m*twist_phase - distance*theta(k)));
        delays = distance * tau1{m}(k);
        for l = 0:N-1
            coefficients = [-kappa, slopes .* exp(1i*l .* twist_phase)];     % coefficients of x(t), x(t - delays)
            A = reshape(coefficients, 1, 1, numel(coefficients));
            eigenvalues = dde_stst_eig_cheb(A, delays, 'minimal_real_part', -2, 'root_accuracy', 1e-8, ...
                                            'max_number_of_eigenvalues', 20, 'min_number_of_eigenvalues', 10);
            Npos{m}(k) = Npos{m}(k) + sum(real(eigenvalues.l0) > 1e-3);
        end
    end
    if tau1{m}(1) == 0
        Npos{m}(1) = Npos{m}(2);
    end
    fprintf('psi_%d done\n', mod(m, N));
end

save(fullfile('data', 'twisted_state_branches.mat'), 'theta', 'Omega', 'tau1', 'Npos');
fprintf('wrote data/twisted_state_branches.mat\n');
