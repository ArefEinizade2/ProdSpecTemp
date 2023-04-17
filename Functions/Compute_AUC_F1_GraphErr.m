function [AUC, F1, Graph_err] = Compute_AUC_F1_GraphErr(W, W_est)
%% A function to calculte the graph learning metrics AUC, F1, and Relative Graph Estiomation (RGE): 
%% In the case of finding its contents useful for your research work, kindly please cite the following paper. Thanks a lot for your attention.
% [1] Einizade, Aref, and Sepideh Hajipour Sardouie. "Learning Product Graphs from Spectral Templates." arXiv preprint arXiv:2211.02893 (2022).

[~, ~, ~, AUC] = perfcurve(vec(full(double(squareform_sp_Mine(W)))), vec(full(double(squareform_sp_Mine(W_est)))), 1);

[tpr, fpr] = roc(vec(full(double(squareform_sp_Mine(W))))', vec(full(double(squareform_sp_Mine(W_est))))');

F = (tpr) ./ (fpr + tpr); F = F(~isnan(F)); 

F1 = mean(F); 

Graph_err = (norm(W - W_est, 'fro')/norm(W, 'fro'))^2;

end
