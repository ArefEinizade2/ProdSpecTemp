function W_prod = CartProdTwoGraphs(W1, W2)
%% A function to generate the Cartesian graph product of two factor graph adjacencies:
%% In the case of finding its contents useful for your research work, kindly please cite the following paper. Thanks a lot for your attention.
% [1] Einizade, Aref, and Sepideh Hajipour Sardouie. "Learning Product Graphs from Spectral Templates." arXiv preprint arXiv:2211.02893 (2022).

N1 = size(W1, 1); N2 = size(W2, 1);

W_prod = kron(W1, eye(N2)) + kron(eye(N1), W2);


end