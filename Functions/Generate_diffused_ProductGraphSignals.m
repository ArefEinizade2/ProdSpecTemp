function X = Generate_diffused_ProductGraphSignals(W, T, SNR, ProductType)
%% A function for Generating diffused product graph signals
%% In the case of finding its contents useful for your research work, kindly please cite the following paper. Thanks a lot for your attention.
% [1] Einizade, Aref, and Sepideh Hajipour Sardouie. "Learning Product Graphs from Spectral Templates." arXiv preprint arXiv:2211.02893 (2022).
%% Inputs: 
% >>> W: The true product graph

% >>> T: The number of graph signals

% >>> SNR: Signal to Noise Ratio (SNR) in db

% >>> ProductType: Type of graph product, 'Cartesian' or 'Strong'.

% >>> param.rho_init; % The initial values of rho

% >>> MaxIters = param.MaxIters; % Maximum ietartions of the algorithm

% >>> cnt = param.cnt; % A user-defined multiplying constant for increasing rho
%% Outpus: 
% >>> W [N x N]: The learned adjacency matrix
%%
N = size(W, 1); % The number of nodes of the product graph
%% Desiging the separable graph filters:
switch ProductType    
    case 'Cartesian'
        tau = 1/max(sum(W));
        H = expm(tau * W);
    case 'Strong'
        H = eye(N, N) + W;
    otherwise
        disp('Wrong content!')
end

%%  Diffusion process for generation of (possibly noisy) graph signals:
X = zeros(N, T); 

for t = 1 : T
    
  y = randn(N, 1); % the excitation vector y in the paper
      
  X(:, t) = H * y;
  
end


Noise = randn(size(X));

alpha = (10^(-SNR/20)) * norm(X, 'fro')/norm(Noise, 'fro');

X = X + alpha * Noise;


end