function W = SpecTemp_IALM(M, Gamma, DD, param)
%% The implementation of the SpecTemp-IALM method for accelerated learning graphs from spectral templates
%% In the case of finding its contents useful for your research work, kindly please cite the following paper. Thanks a lot for your attention.
% [1] Einizade, Aref, and Sepideh Hajipour Sardouie. "Learning Product Graphs from Spectral Templates." arXiv preprint arXiv:2211.02893 (2022).
%% Inputs: 
% >>> M = M_d * M_h, where M_d: duplication matrix and M_h: a matrix that vech(x) = M_h * vech(x), as explained in the paper

% >>> Gamma: Gamma matrix in paper

% >>> DD: Eigenvalue matrix of the Cx_i

% >>> param.Tol: % The user-defined toleronce for reaching convergence, e.g., 1e-2.

% >>> param.rho_init; % The initial values of rho

% >>> MaxIters = param.MaxIters; % Maximum ietartions of the algorithm

% >>> cnt = param.cnt; % A user-defined multiplying constant for increasing rho
%% Outpus: 
% >>> W [N x N]: The learned adjacency matrix

%% Initializations:

lambda_vec = diag(DD); % initializing the lambda with the eigenvalues of the Cx_i

Gamma_pinv = pinv(Gamma); % Moore-Penrose pseudo-inverse of Gamma

N = length(lambda_vec); % The number of nodes of the graph

m = N * (N - 1) / 2; % The number of non-diagonal upper diagonal entries of the graph

w_k_plus = ones(m,1); % initialization of w

s_k_plus = ones(m,1); % initialization of s

lambda_k_plus = lambda_vec(:); % initialization of lambda

gamma1_k_plus = w_k_plus - Gamma * lambda_k_plus; % initialization of gamma1

gamma2_k_plus = w_k_plus - s_k_plus; % initialization of gamma2

%% User-defined settings:

Tol = param.Tol; % The user-defined toleronce for reaching convergence

rho = param.rho_init; % The initial values of rho

MaxIters = param.MaxIters; % Maximum ietartions of the algorithm

cnt = param.cnt; % A user-defined multiplying constant for increasing rho

%% While loop for reaching convergence:

iter = 0; 

Start = true;

while ( (Start) ||  (norm(lambda_k_plus - lambda_k) / norm(lambda_k) > Tol && iter <= MaxIters) )
    %% initialization:
    Start = false;

    iter = iter + 1;
    
    s_k = s_k_plus; lambda_k = lambda_k_plus; 

    gamma1_k = gamma1_k_plus; gamma2_k = gamma2_k_plus;
    %% Update w:
    param.verbose = 0;
    
    x = (rho * (Gamma * lambda_k + s_k) + gamma1_k + gamma2_k) / (rho + rho);
        
    w_k_plus = prox_l1(x, 1 / (rho + rho), param);
        
    %% Update lambda:
    
    lambda_k_plus = Gamma_pinv * (w_k_plus - (gamma1_k / rho));
    
    %% Update s:
    
    s_k_plus = w_k_plus - (gamma2_k / rho);
    
    s_k_plus(s_k_plus < 0) = 0; % Euclidean projection on non-negative set
    
    % The additional constraint to avoid trival all-zero solution by enforing the first node to have unit degree:

%     s_k_plus(1:N-1) = s_k_plus(1:N-1)/sum(s_k_plus(1:N-1)); 
    
    % Another additional constraint to avoid trival all-zero solution by enforcing w to have unit norm, which gets better results:

    s_k_plus = s_k_plus / norm(s_k_plus, 'fro'); 
    
    %% Update gamma1 & gamma2:
    gamma1_k_plus = gamma1_k - rho * (w_k_plus - Gamma * lambda_k_plus);
    
    gamma2_k_plus = gamma2_k - rho * (w_k_plus - s_k_plus);
        
    %% Update rho:
    rho = rho * cnt;

end

if iter < MaxIters - 1
    
    disp('>>> Converged...')
    
else

    disp('>>> Not converged...')

end

%% Output:
W = reshape(M * s_k_plus, [N, N]);
W = W - diag(diag(W));
W = W / max(W(:));
end