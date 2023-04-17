function W_prod = KronProdTwoGraphs(W1, W2)
%% A function to generate the Kronecker graph product of two factor graph adjacencies:
%% In the case of finding its contents useful for your research work, kindly please cite the following paper. Thanks a lot for your attention.
% [1] Einizade, Aref, and Sepideh Hajipour Sardouie. "Learning Product Graphs from Spectral Templates." arXiv preprint arXiv:2211.02893 (2022).

W_prod = kron(W1, W2);


end