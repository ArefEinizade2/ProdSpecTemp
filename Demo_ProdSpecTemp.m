%% This a demo for learning product graphs by recovering three factor graphs from a stream of product graph signals
%% In the case of finding its contents useful for your research work, kindly please cite the following paper. Thanks a lot for your attention.
% [1] Einizade, Aref, and Sepideh Hajipour Sardouie. "Learning Product Graphs from Spectral Templates." arXiv preprint arXiv:2211.02893 (2022).
clear; close all; clc;
%% The codes were written in MATLAB 2018b, and the GSPBOX and UNLocBoX toolboxes are also needed to run this demo.
%% Please note that in the case of getting an error named "Less than two classes are found in the array of true class labels.", just run the demo another time. In fact, in this case, the function gsp_erdos_renyi from GSPBOX has produced an all-zero adjacency for Erdos-Renyi graphs.
%% Besides, in the case of changing the current settings or hyperparameters, please tune the other parameters to get the desired results.
rng(8);
t0 = tic;
%% Add needed functions:
addpath('./Functions')
%% Settings for simulations:
N1 = 12; % The number of nodes of the first Erdos-Renyi factor graph
p1 = 0.3; % The edge probability of the first Erdos-Renyi factor graph
N2 = 10; % The number of nodes of the second Erdos-Renyi factor graph
p2 = 0.3; % The edge probability of the second Erdos-Renyi factor graph
N3 = 9; % The number of nodes of the third Erdos-Renyi factor graph
p3 = 0.3; % The edge probability of the third Erdos-Renyi factor graph
param_ER.connected = 1; % For generating connected Erdos-Renyi factor graphs
SNR = -10; % Signal to Noise Ratio (SNR) in db 
ProductIdx_vec = [1, 2]; % The product graph indices 
T_vec = 10.^[1, 2, 3, 4]; % The number of temporal samples, here 10^[1:4]
NumReal = 10; % The number of realizations
Prod_names = {'Cartesian', 'Strong'}; % Types of graph products, here Cartesian and Strong product graphs
%% needed user-defined params for SpecTemp-IALM method:
param.Tol = 1e-1; % The user-defined toleronce for reaching convergence
param.rho_init = 1; % The initial values of rho
param.MaxIters = 50; % Maximum ietartions of the algorithm
param.cnt = 1e3; % A user-defined multiplying constant for increasing rho

%%
for gen = 1 : NumReal

G1 = gsp_erdos_renyi(N1, p1, param_ER); % Generating the first Erdos-Renyi factor graph
G2 = gsp_erdos_renyi(N2, p2, param_ER); % Generating the second Erdos-Renyi factor graph
G3 = gsp_erdos_renyi(N3, p3, param_ER); % Generating the third Erdos-Renyi factor graph

W1 = full(double(G1.W)); % True adjacency matrix of the first Erdos-Renyi factor graph
W2 = full(double(G2.W)); % True adjacency matrix of the second Erdos-Renyi factor graph
W3 = full(double(G3.W)); % True adjacency matrix of the third Erdos-Renyi factor graph

W_1_true_array(:, :, gen) = W1; % Stacking the first Erdos-Renyi factor graphs
W_2_true_array(:, :, gen) = W2; % Stacking the second Erdos-Renyi factor graphs
W_3_true_array(:, :, gen) = W3; % Stacking the third Erdos-Renyi factor graphs


W_cart = CartProdMoreGraphs({W1, W2, W3}); % Generating Cartesian graph product of W1, W2, and W3
W_kron = KronProdMoreGraphs({W1, W2, W3}); % Generating Kronecker graph product of W1, W2, and W3
W_strong = StrongProdMoreGraphs({W1, W2, W3}); % Generating Strong graph product of W1, W2, and W3


W_cell = {W_cart, W_strong}; % A cell to consider the true product graphs


for ProductIdx = ProductIdx_vec % Loop on product graphs, here Cartesian and Strong product graphs

W = W_cell{ProductIdx}; % True product graph in the current loop

%% generate diffused product graph signals:
X_all = Generate_diffused_ProductGraphSignals(W, T_vec(end), SNR, Prod_names{ProductIdx});
%%
for T = T_vec % Loop on the number of samples (product graph signals)
    
disp(['>>>>>>>>>>>>>> Realization #', num2str(gen), ', Product type: ', Prod_names{ProductIdx}, ', Number of samples: ', num2str(T), ' >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'])
%% Learn graphs:
X = X_all(:, 1 : T); % The number of samples in the current loop (T)
W_cell_est = ProdSpecTemp(X, [N3, N2, N1], param); % Learning the underlying factor graphs via the proposed ProdSpecTemp
W1_est = W_cell_est{1}; % The recovered first factor graph 
W1_est_array(find(ProductIdx_vec==ProductIdx), find(T_vec==T), :, :, gen) = W1_est; % Stacking the recovered first factor graph
W2_est = W_cell_est{2}; % The recovered second factor graph
W2_est_array(find(ProductIdx_vec==ProductIdx), find(T_vec==T), :, :, gen) = W2_est; % Stacking the recovered second factor graph
W3_est = W_cell_est{3}; % The recovered third factor graph
W3_est_array(find(ProductIdx_vec==ProductIdx), find(T_vec==T), :, :, gen) = W3_est; % Stacking the recovered third factor graph
%% Evaluating graph recovery performance:

% Computing AUC, F1, and L2-error in recovering the first factorg graph:
[AUC(1, find(T_vec==T), find(ProductIdx_vec==ProductIdx), gen), ...
    F1(1, find(T_vec==T), find(ProductIdx_vec==ProductIdx), gen), ...
    Graph_err(1, find(T_vec==T), find(ProductIdx_vec==ProductIdx), gen)] = Compute_AUC_F1_GraphErr(W1, W1_est); 

% Computing AUC, F1, and L2-error in recovering the first factorg graph:
[AUC(2, find(T_vec==T), find(ProductIdx_vec==ProductIdx), gen), ...
    F1(2, find(T_vec==T), find(ProductIdx_vec==ProductIdx), gen), ...
    Graph_err(2, find(T_vec==T), find(ProductIdx_vec==ProductIdx), gen)] = Compute_AUC_F1_GraphErr(W2, W2_est);

% Computing AUC, F1, and L2-error in recovering the first factorg graph:
[AUC(3, find(T_vec==T), find(ProductIdx_vec==ProductIdx), gen), ...
    F1(3, find(T_vec==T), find(ProductIdx_vec==ProductIdx), gen), ...
    Graph_err(3, find(T_vec==T), find(ProductIdx_vec==ProductIdx), gen)] = Compute_AUC_F1_GraphErr(W3, W3_est);

end


end


end
%% Plotting the averaged AUC, F1, and L2-error metrics in recoverying all factror graphs:

f = figure; 

f.Position = [50 50 725 725];

% Graph #1, AUC metric:
subplot(3, 3, 1); AUC_temp = mean(squeeze(AUC(1, :, :, :)), 3); 
plot(1:length(T_vec), AUC_temp(:, 1)', 'b->', 1:length(T_vec), AUC_temp(:, 2)', 'r--o', 'LineWidth', 2); 
title('G_{1}', 'fontWeight', 'bold', 'fontsize',12); 
ylabel('AUC', 'fontWeight', 'bold', 'fontsize',10)
xticks(1:length(T_vec)); xticklabels({'10^1', '10^2', '10^3'}); xlabel('T', 'fontWeight', 'bold', 'fontsize',10)

% Graph #1, F1 metric:
subplot(3, 3, 2); F1_temp = mean(squeeze(F1(1, :, :, :)), 3); 
plot(1:length(T_vec), F1_temp(:, 1)', 'b->', 1:length(T_vec), F1_temp(:, 2)', 'r--o', 'LineWidth', 2); 
title('G_{1}', 'fontWeight', 'bold', 'fontsize',12); 
ylabel('F1', 'fontWeight', 'bold', 'fontsize',10)
xticks(1:length(T_vec)); xticklabels({'10^1', '10^2', '10^3'}); xlabel('T', 'fontWeight', 'bold', 'fontsize',10)

% Graph #1, L2-error metric:
subplot(3, 3, 3); Graph_err_temp = mean(squeeze(Graph_err(1, :, :, :)), 3); 
plot(1:length(T_vec), Graph_err_temp(:, 1)', 'b->', 1:length(T_vec), Graph_err_temp(:, 2)', 'r--o', 'LineWidth', 2); 
title('G_{1}', 'fontWeight', 'bold', 'fontsize',12); 
ylabel('edge L_2', 'fontWeight', 'bold', 'fontsize',10)
xticks(1:length(T_vec)); xticklabels({'10^1', '10^2', '10^3'}); xlabel('T', 'fontWeight', 'bold', 'fontsize',10)

% Graph #2, AUC metric:
subplot(3, 3, 4); AUC_temp = mean(squeeze(AUC(2, :, :, :)), 3); 
plot(1:length(T_vec), AUC_temp(:, 1)', 'b->', 1:length(T_vec), AUC_temp(:, 2)', 'r--o', 'LineWidth', 2); 
title('G_{2}', 'fontWeight', 'bold', 'fontsize',12); 
ylabel('AUC', 'fontWeight', 'bold', 'fontsize',10)
xticks(1:length(T_vec)); xticklabels({'10^1', '10^2', '10^3'}); xlabel('T', 'fontWeight', 'bold', 'fontsize',10)

% Graph #2, F1 metric:
subplot(3, 3, 5); F1_temp = mean(squeeze(F1(2, :, :, :)), 3); 
plot(1:length(T_vec), F1_temp(:, 1)', 'b->', 1:length(T_vec), F1_temp(:, 2)', 'r--o', 'LineWidth', 2); 
title('G_{2}', 'fontWeight', 'bold', 'fontsize',12); 
ylabel('F1', 'fontWeight', 'bold', 'fontsize',10)
xticks(1:length(T_vec)); xticklabels({'10^1', '10^2', '10^3'}); xlabel('T', 'fontWeight', 'bold', 'fontsize',10)

% Graph #2, L2-error metric:
subplot(3, 3, 6); Graph_err_temp = mean(squeeze(Graph_err(2, :, :, :)), 3); 
plot(1:length(T_vec), Graph_err_temp(:, 1)', 'b->', 1:length(T_vec), Graph_err_temp(:, 2)', 'r--o', 'LineWidth', 2); 
title('G_{2}', 'fontWeight', 'bold', 'fontsize',12); 
ylabel('edge L_2', 'fontWeight', 'bold', 'fontsize',10)
xticks(1:length(T_vec)); xticklabels({'10^1', '10^2', '10^3'}); xlabel('T', 'fontWeight', 'bold', 'fontsize',10)

% Graph #3, AUC metric:
subplot(3, 3, 7); AUC_temp = mean(squeeze(AUC(3, :, :, :)), 3); 
plot(1:length(T_vec), AUC_temp(:, 1)', 'b->', 1:length(T_vec), AUC_temp(:, 2)', 'r--o', 'LineWidth', 2); 
title('G_{3}', 'fontWeight', 'bold', 'fontsize',12); 
ylabel('AUC', 'fontWeight', 'bold', 'fontsize',10)
xticks(1:length(T_vec)); xticklabels({'10^1', '10^2', '10^3'}); xlabel('T', 'fontWeight', 'bold', 'fontsize',10)

% Graph #3, F1 metric:
subplot(3, 3, 8); F1_temp = mean(squeeze(F1(3, :, :, :)), 3); 
plot(1:length(T_vec), F1_temp(:, 1)', 'b->', 1:length(T_vec), F1_temp(:, 2)', 'r--o', 'LineWidth', 2); 
title('G_{3}', 'fontWeight', 'bold', 'fontsize',12); 
ylabel('F1', 'fontWeight', 'bold', 'fontsize',10)
xticks(1:length(T_vec)); xticklabels({'10^1', '10^2', '10^3'}); xlabel('T', 'fontWeight', 'bold', 'fontsize',10)

% Graph #2, L2-error metric:
subplot(3, 3, 9); Graph_err_temp = mean(squeeze(Graph_err(3, :, :, :)), 3); 
plot(1:length(T_vec), Graph_err_temp(:, 1)', 'b->', 1:length(T_vec), Graph_err_temp(:, 2)', 'r--o', 'LineWidth', 2); 
title('G_{3}', 'fontWeight', 'bold', 'fontsize',12); 
ylabel('edge L_2', 'fontWeight', 'bold', 'fontsize',10)
xticks(1:length(T_vec)); xticklabels({'10^1', '10^2', '10^3'}); xlabel('T', 'fontWeight', 'bold', 'fontsize',10)


lg = legend('applied on Cartesian prodcuts', 'applied on Strong prodcuts', 'Orientation', 'horizontal', 'Location', 'NE');

lg.Position = [0.25 0.03 0.5 0.01] ;

%% Plotting the true and learned third factor graph from Cartesian and Strong product graph signals: 
GraphCateg = length(ProductIdx_vec);
n = 3;
rows = GraphCateg;
columns = length(T_vec) + 1;
W1_true_mean = squeeze(mean(W_1_true_array, 3));    
W1_true_mean = W1_true_mean/max(W1_true_mean(:));
W2_true_mean = squeeze(mean(W_2_true_array, 3));    
W2_true_mean = W2_true_mean/max(W2_true_mean(:));
W3_true_mean = squeeze(mean(W_3_true_array, 3));    
W3_true_mean = W3_true_mean/max(W3_true_mean(:));

Counter = 0;

f = figure;

f.Position = [100 100 1100 400];

for ProductIdx = ProductIdx_vec
    
for T = [T_vec, 1]
    
    Counter = Counter + 1;
        
    if mod(Counter, columns) ~= 0
    
        W3_est_mean = squeeze(mean(W3_est_array(find(ProductIdx_vec==ProductIdx), find(T_vec==T), :, :, :), 5));
        W3_est_mean = W3_est_mean/max(W3_est_mean(:));
    
        subplot(rows, columns, Counter); imagesc(W3_est_mean); title(['G3, ', Prod_names{ProductIdx}, ', T: ', num2str(T)])
    else
        
        subplot(rows, columns, Counter); imagesc(W3_true_mean); title('G3, true')
        
    end
end

end
cb = colorbar(); 
%manually change the position of colorbar to what you want 
cb.Position = [0.95 0.3000 0.0100 0.5000] ;
%%
t1 = toc(t0);

disp(['>>>>> run-time: ', num2str(round(t1/60,2)), ' minutes']);

