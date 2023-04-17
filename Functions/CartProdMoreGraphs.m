function W_prod = CartProdMoreGraphs(W_cell)
%% A function to generate the Cartesian graph product of posssibly more than two factor graph adjacencies:
%% In the case of finding its contents useful for your research work, kindly please cite the following paper. Thanks a lot for your attention.
% [1] Einizade, Aref, and Sepideh Hajipour Sardouie. "Learning Product Graphs from Spectral Templates." arXiv preprint arXiv:2211.02893 (2022).

n = length(W_cell);

W_prod = CartProdTwoGraphs(W_cell{1}, W_cell{2});

for i = 3 : n
    
    W_prod = CartProdTwoGraphs(W_prod, W_cell{i});
    
end


end