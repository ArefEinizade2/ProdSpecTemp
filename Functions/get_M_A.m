function [M, Gamma, t_Ours] = get_M_A(U)

%% A function to obtain Gamma and M = M_d * M_h, where M_d: duplication matrix and M_h: a matrix that vech(x) = M_h * vech(x), as explained in the paper
%% In the case of finding its contents useful for your research work, kindly please cite the following paper. Thanks a lot for your attention.
% [1] Einizade, Aref, and Sepideh Hajipour Sardouie. "Learning Product Graphs from Spectral Templates." arXiv preprint arXiv:2211.02893 (2022).

N = size(U, 1);

M_p = full(duplication_matrix(N));

vech_idx = pinv(M_p)*vec(eye(N));


t1 = tic;

Sel = eye(length(vech_idx));

Sel = Sel(vech_idx == 0, :);

M_N = pinv(Sel);

M = M_p * M_N;

pinv_M = pinv(M);
%
V_p = [];

for i = 1 : N
    
   V_p(:, i) = vec(U(:, i) * U(:, i)'); 
   
end

Gamma = pinv_M * V_p;

t_Ours = toc(t1);


end