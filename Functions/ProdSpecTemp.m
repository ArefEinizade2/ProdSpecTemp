function W_cell = ProdSpecTemp(X, dims, param)
%% The implementation of the ProdSpecTemp method for learning product graphs
%% In the case of finding its contents useful for your research work, kindly please cite the following paper. Thanks a lot for your attention.
% [1] Einizade, Aref, and Sepideh Hajipour Sardouie. "Learning Product Graphs from Spectral Templates." arXiv preprint arXiv:2211.02893 (2022).
%% Inputs: 
% >>> X: N x T, where N = P1 x P2 x ... x Pn, where Pi is the number of the i-th factor graphs

% >>> dims = [Pn,...,P2, P1], a list containing the number of nodes of the factor graphs

% >>> param: needed user-defined params for SpecTemp-IALM method
%% Outpus: 
% >>> W_cell = {W1, W2, ..., Wn}, a cell containing the learned adjacencies of the factor graphs
%%
T = size(X, 2);

n = length(dims); % The number of 
%%

for i = 1 : n % Loop on the different modes of the initial tensor product graph siganls
    %% Estimating the mode-i covariance matrix via sample covariance:
    
    Cx_i = 0;
    
    for t = 1 : T
                
        X_t = reshape(X(:, t), dims); % Build the initial tensor
        
        X_t_i = Unfold_a_tensor(X_t, dims, i); % Mode-i unfolding of the tensor
       
        Cx_i = Cx_i + (X_t_i * X_t_i' / T); 
    
    end
        
    %%
    
    [V, D] = eig(Cx_i); % EVD of Cx_i 

    [M, Gamma] = get_M_A(V); % Obtain Gamma and M = M_d * M_h, where M_d: duplication matrix and M_h: a matrix that vech(x) = M_h * vech(x), as explained in the paper

    W = SpecTemp_IALM(M, Gamma, D, param); % Lines 4-13 of the SpecTemp-IALM in the paper

%     W = W - diag(diag(W));

    W_cell{n - i + 1} = W;
    
    disp(['mode-', num2str(i), ' recovered'])
    
end


end