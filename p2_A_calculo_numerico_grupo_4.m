                                      % Análise 1.1 %

%Foram escolhidos os pares MDVP:Fo(Hz) e MDVP:Flo(Hz) para o grupo sem Parkinson, pois ambos mostraram média e
%dispersão mais altas neste grupo, indicando maior estabilidade na frequência vocal. Para o grupo com Parkinson,
%foi selecionado MDVP:Shimmer(dB) e MDVP:APQ , que apresentaram médias e dispersão mais elevadas nos pacientes com a
%doença, refletindo a instabilidade vocal característica. Esses pares evidenciam diferenças importantes entre os grupos.

                                       % Análise 2.2 %

%Média e Valor Absoluto:

%Para métricas de frequência, como Fo(Hz), Fhi(Hz) e Flo(Hz), a média dos valores para o grupo de controle é
%consistentemente maior que para o grupo Parkinson, sugerindo uma maior capacidade de controle e variação da frequência
%vocal em indivíduos sem Parkinson. Em variáveis relacionadas à instabilidade e inconsistência vocal, como Jitter(%),
%Shimmer e APQ, as médias são geralmente mais altas no grupo Parkinson. Isso indica que a voz de indivíduos com Parkinson
%apresenta maior variação e instabilidade, características associadas ao tremor vocal e à rigidez muscular.

%Dispersão e Desvio Padrão:

%Em termos de dispersão, o grupo de Parkinson apresenta maior desvio padrão para quase todas as métricas, especialmente
%para Jitter e Shimmer, métricas associadas à instabilidade da frequência e amplitude vocal. Esse aumento na dispersão
%reflete a variabilidade de sintomas vocais entre os pacientes, como flutuações na capacidade de manter uma frequência
%estável. O maior desvio padrão no grupo Parkinson para as métricas de jitter e shimmer indica uma maior variabilidade
%da intensidade e frequência de voz, o que é esperado devido às dificuldades motoras que afetam o controle da fala.

%Interpretação Clínica:

%A média mais baixa em variáveis de frequência para o grupo Parkinson sugere que a voz desses pacientes pode ser menos modulada
%e com menor amplitude de variação. O aumento no desvio padrão indica uma maior dificuldade em manter consistência vocal, o que
%é característico dos sintomas motores de Parkinson, como a rigidez e o tremor que dificultam o controle preciso dos músculos vocais.

                                        % Anáise 3.3 %













% Load necessary packages and clear workspace
pkg load statistics;  % Carrega o pacote de estatísticas, necessário para análises estatísticas.
clc;  % Limpa o console, removendo qualquer saída anterior.
clearvars;  % Limpa todas as variáveis do espaço de trabalho, garantindo um ambiente limpo para a execução.

% ===================== Function Definitions =====================

% Function for Pearson correlation
function corr = pearson_corr(x, y)
    xm = x - mean(x);  % Centraliza os dados de x subtraindo a média.
    ym = y - mean(y);  % Centraliza os dados de y subtraindo a média.
    corr = sum(xm .* ym) / sqrt(sum(xm .^ 2) * sum(ym .^ 2));  % Calcula a correlação de Pearson entre x e y.
end


% Function for linear regression
function [slope, intercept] = linear_regression(x, y)
    n = length(x);  % Obtém o número de pontos de dados.
    X = [ones(n, 1), x];  % Cria a matriz de design X com uma coluna de 1s (para o intercepto).
    A = X' * X;  % Calcula a matriz A como o produto da transposta de X com X.
    b = X' * y;  % Calcula o vetor b como o produto da transposta de X com y.
    coeffs = solve_lu(A, b);  % Usa a função solve_lu para resolver o sistema linear Ax = b.
    intercept = coeffs(1);  % O primeiro coeficiente é o intercepto da regressão.
    slope = coeffs(2);  % O segundo coeficiente é a inclinação da linha de regressão.
end

% Function for performing double regression and obtaining metrics
function [coef1, coef2, Sr1, Sr2, r2_1, r2_2, syx1, syx2] = regressao_dupla(x1, y1, x2, y2)
    % Number of data points
    n = length(x1);  % Obtém o número de pontos de dados.

    % Model 1
    A1 = [x1, ones(n, 1)];  % Cria a matriz de design A1 para o primeiro modelo (x1).
    coef1 = lu_solve(A1, y1);  % Resolve o sistema para obter os coeficientes do modelo 1.
    y1_pred = A1 * coef1;  % Calcula as previsões para y1 usando os coeficientes obtidos.
    Sr1 = sum((y1 - y1_pred).^2);  % Calcula a soma dos resíduos quadráticos para o modelo 1.
    St1 = sum((y1 - mean(y1)).^2);  % Calcula a soma total dos quadrados para o modelo 1.
    r2_1 = 1 - (Sr1 / St1);  % Calcula o coeficiente de determinação R² para o modelo 1.
    syx1 = sqrt(Sr1 / (n - 2));  % Calcula o erro padrão da estimativa para o modelo 1.

    % Model 2
    A2 = [x2, ones(n, 1)];  % Cria a matriz de design A2 para o segundo modelo (x2).
    coef2 = lu_solve(A2, y2);  % Resolve o sistema para obter os coeficientes do modelo 2.
    y2_pred = A2 * coef2;  % Calcula as previsões para y2 usando os coeficientes obtidos.
    Sr2 = sum((y2 - y2_pred).^2);  % Calcula a soma dos resíduos quadráticos para o modelo 2.
    St2 = sum((y2 - mean(y2)).^2);  % Calcula a soma total dos quadrados para o modelo 2.
    r2_2 = 1 - (Sr2 / St2);  % Calcula o coeficiente de determinação R² para o modelo 2.
    syx2 = sqrt(Sr2 / (n - 2));  % Calcula o erro padrão da estimativa para o modelo 2.
end

% Função solve_lu para resolver Ax = b usando decomposição LU
function x = solve_lu(A, b)
    [L, U, P] = lu(A);  % Decomposição LU com pivotamento
    y = substituicao_progressiva(L, P * b);  % Resolve L * y = P * b
    x = substituicao_regressiva(U, y);       % Resolve U * x = y
endfunction

% Função para substituição progressiva (resolvendo L*y = b)
function y = substituicao_progressiva(L, b)
    n = length(b);
    y = zeros(n, 1);
    for i = 1:n
        y(i) = (b(i) - L(i, 1:i-1) * y(1:i-1)) / L(i, i);
    endfor
endfunction

% Função para substituição regressiva (resolvendo U*x = y)
function x = substituicao_regressiva(U, y)
    n = length(y);
    x = zeros(n, 1);
    for i = n:-1:1
        x(i) = (y(i) - U(i, i+1:n) * x(i+1:n)) / U(i, i);
    endfor
endfunction

% ===================== Main Script ===================== %

% Load the data from the CSV file
data = csvread('parkinson_mdvp.csv', 1, 0); % Ignora o cabeçalho ao ler o arquivo CSV.

% Definition of column names for future reference
columns_names = {'Fo(Hz)', 'Fhi(Hz)', 'Flo(Hz)', ...  % Nomes das colunas representando características da voz.
                 'Jitter(%)', 'Jitter(Abs)', 'RAP', ...  % Medidas de irregularidade e variação.
                 'PPQ', 'Shimmer', 'Shimmer(dB)', ...  % Medidas de amplitude e variação na amplitude.
                 'APQ'};  % Outros parâmetros da análise.

% Logarithmic normalization and data scaling
data(:, 1:end-1) = log(data(:, 1:end-1) + 1);  % Aplica a normalização logarítmica para evitar valores extremos.
data(:, 1:end-1) = (data(:, 1:end-1) - mean(data(:, 1:end-1))) ./ std(data(:, 1:end-1));  % Z-normaliza os dados para média 0 e desvio padrão 1.

% Modify the labels of the last column to 'Parkinson' and 'Control'
group_labels = cell(size(data, 1), 1);  % Inicializa uma célula para armazenar os rótulos dos grupos.
group_labels(data(:, end) == 1) = {'Parkinson'};  % Atribui 'Parkinson' aos dados com valor 1 na última coluna.
group_labels(data(:, end) == 0) = {'Control'};  % Atribui 'Control' aos dados com valor 0 na última coluna.

[unique_groups, ~, group_indices] = unique(group_labels);  % Obtém os grupos únicos e os índices para categorização.

% Create boxplots of the variables
figure;  % Cria uma nova figura para os boxplots.
for i = 1:length(columns_names)  % Loop através de cada coluna de dados para criar os boxplots.
    subplot(4, 3, i);  % Cria um subplot em uma grade 4x3 para acomodar todos os gráficos.
    boxplot(data(:, i), group_indices);  % Cria o boxplot para a i-ésima variável, categorizado por grupos.
    title(columns_names{i});  % Define o título do gráfico como o nome da coluna correspondente.

    % t-test between the two categories
    group1 = data(group_indices == 1, i);  % Extrai os dados do grupo 'Parkinson'.
    group2 = data(group_indices == 2, i);  % Extrai os dados do grupo 'Control'.

    if length(group1) > 1 && length(group2) > 1  % Verifica se ambos os grupos têm mais de um elemento.
        [~, p_value] = ttest2(group1, group2);  % Realiza um teste t de duas amostras entre os grupos.
        minus_log10_p = -log10(p_value);  % Calcula o valor negativo do logaritmo na base 10 do p-valor.
    else
        minus_log10_p = NaN;  % Se um dos grupos não tem dados suficientes, define o p-valor como NaN.
    end

    mean_group1 = mean(group1);  % Calcula a média do grupo 'Parkinson'.
    mean_group2 = mean(group2);  % Calcula a média do grupo 'Control'.
    std_group1 = std(group1);  % Calcula o desvio padrão do grupo 'Parkinson'.
    std_group2 = std(group2);  % Calcula o desvio padrão do grupo 'Control'.

    ylim_current = ylim();  % Obtém os limites atuais do eixo y do gráfico.
    y_position = ylim_current(2) + (ylim_current(2) - ylim_current(1)) * 0.15;  % Calcula a posição y para anotações.
    text(mean(xlim()), y_position, sprintf('-log_{10}(p) = %.2f', minus_log10_p), ...  % Adiciona texto ao gráfico com o valor -log10(p).
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontWeight', 'bold');

    ylim([ylim_current(1), y_position + (ylim_current(2) - ylim_current(1)) * 0.2]);  % Ajusta os limites do eixo y para acomodar o texto.
    hold on;  % Mantém o gráfico atual para adicionar mais elementos.
    positions = get(gca, 'XTick');  % Obtém as posições atuais das marcas no eixo x.
    y_offset = (ylim_current(2) - ylim_current(1)) * 0.05;  % Define um pequeno deslocamento para a anotação da média.

    text(positions(1), mean_group1 + y_offset, sprintf('%.2f ± %.2f', mean_group1, std_group1), ...  % Adiciona a média e desvio padrão do grupo 'Parkinson' ao gráfico.
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Color', 'blue', 'FontWeight', 'bold');

    text(positions(2), mean_group2 + y_offset, sprintf('%.2f ± %.2f', mean_group2, std_group2), ...  % Adiciona a média e desvio padrão do grupo 'Control' ao gráfico.
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Color', 'red', 'FontWeight', 'bold');
    hold off;  % Libera o gráfico atual para futuras alterações.
    set(gca, 'XTickLabel', unique_groups);  % Define os rótulos do eixo x com os nomes dos grupos únicos.
end

% Adiciona uma anotação ao topo da figura.
annotation('textbox', [0.5, 0.98, 0, 0], 'String', 'Boxplots of Variables by Category', ...
           'FitBoxToText', 'on', 'HorizontalAlignment', 'center', 'FontSize', 12, 'LineStyle', 'none');

% Separate data of patients with and without Parkinson
% Aqui, os dados são filtrados em duas variáveis separadas:
% 'data_parkinson' contém os dados de pacientes com Parkinson (onde a 11ª coluna é 1)
% e 'data_no_parkinson' contém os dados de pacientes sem a doença (onde a 11ª coluna é 0).
data_parkinson = data(data(:, 11) == 1, :);
data_no_parkinson = data(data(:, 11) == 0, :);

% Analysis for the pairs of columns
% O código a seguir realiza a análise de correlação e regressão linear
% para dois pares de variáveis selecionadas.

% Analisando as colunas 1 e 3 para pacientes com Parkinson
[x1_parkinson, y1_parkinson] = deal(data_parkinson(:, 1), data_parkinson(:, 3));
correlation_1_3_parkinson = pearson_corr(x1_parkinson, y1_parkinson);  % Calcula a correlação de Pearson
[slope_1_3_parkinson, intercept_1_3_parkinson] = linear_regression(x1_parkinson, y1_parkinson);  % Realiza a regressão linear

% Analisando as colunas 1 e 3 para pacientes sem Parkinson
[x1_no_parkinson, y1_no_parkinson] = deal(data_no_parkinson(:, 1), data_no_parkinson(:, 3));
correlation_1_3_no_parkinson = pearson_corr(x1_no_parkinson, y1_no_parkinson);  % Calcula a correlação de Pearson
[slope_1_3_no_parkinson, intercept_1_3_no_parkinson] = linear_regression(x1_no_parkinson, y1_no_parkinson);  % Realiza a regressão linear

% Analisando as colunas 9 e 10 para pacientes com Parkinson
[x2_parkinson, y2_parkinson] = deal(data_parkinson(:, 9), data_parkinson(:, 10));
correlation_9_10_parkinson = pearson_corr(x2_parkinson, y2_parkinson);  % Calcula a correlação de Pearson
[slope_9_10_parkinson, intercept_9_10_parkinson] = linear_regression(x2_parkinson, y2_parkinson);  % Realiza a regressão linear

% Analisando as colunas 9 e 10 para pacientes sem Parkinson
[x2_no_parkinson, y2_no_parkinson] = deal(data_no_parkinson(:, 9), data_no_parkinson(:, 10));
correlation_9_10_no_parkinson = pearson_corr(x2_no_parkinson, y2_no_parkinson);  % Calcula a correlação de Pearson
[slope_9_10_no_parkinson, intercept_9_10_no_parkinson] = linear_regression(x2_no_parkinson, y2_no_parkinson);  % Realiza a regressão linear

% Display analysis results
% Aqui, os resultados das análises de correlação e regressão linear são exibidos na tela.
disp('===========================================================');
disp('                         Analysis 1.3             ');  % Título da análise
disp('===========================================================');

disp('Pair MDVP:Fo(Hz) and MDVP:Flo(Hz):')  % Descrição do primeiro par de variáveis
disp(['With Parkinson: Correlation = ', num2str(correlation_1_3_parkinson), ', Slope = ', num2str(slope_1_3_parkinson), ', Intercept = ', num2str(intercept_1_3_parkinson)]);  % Resultados para pacientes com Parkinson
disp(['Without Parkinson: Correlation = ', num2str(correlation_1_3_no_parkinson), ', Slope = ', num2str(slope_1_3_no_parkinson), ', Intercept = ', num2str(intercept_1_3_no_parkinson)]);  % Resultados para pacientes sem Parkinson
disp('===========================================================');
disp('Pair MDVP:Shimmer(dB) and MDVP:APQ:')  % Descrição do segundo par de variáveis
disp(['With Parkinson: Correlation = ', num2str(correlation_9_10_parkinson), ', Slope = ', num2str(slope_9_10_parkinson), ', Intercept = ', num2str(intercept_9_10_parkinson)]);  % Resultados para pacientes com Parkinson
disp(['Without Parkinson: Correlation = ', num2str(correlation_9_10_no_parkinson), ', Slope = ', num2str(slope_9_10_no_parkinson), ', Intercept = ', num2str(intercept_9_10_no_parkinson)]);  % Resultados para pacientes sem Parkinson
disp('===========================================================');

% Indices of variables of interest
% Seleciona as colunas de interesse para análise adicional.
selected_columns = [9, 10];  % Índices das colunas 'Shimmer(dB)' e 'APQ'.

% Create a figure for boxplots of selected variables
% Cria uma figura para mostrar boxplots das variáveis selecionadas, facilitando a visualização das distribuições.
figure;
for i = 1:length(selected_columns)  % Itera sobre as colunas selecionadas
    col_index = selected_columns(i);  % Acessa o índice da coluna atual
    subplot(1, 2, i);  % Define o layout para apenas 2 gráficos.
    boxplot(data(:, col_index), group_indices);  % Cria um boxplot da variável selecionada, agrupando pelos índices numéricos.
    title(columns_names{col_index});  % Define o título do boxplot com o nome da variável.

    % Perform t-test between the two categories
    % Aqui, realiza um teste t para comparar as médias das duas categorias.
    group1 = data(group_indices == 1, col_index);  % Dados do grupo 'Controle'.
    group2 = data(group_indices == 2, col_index);  % Dados do grupo 'Parkinson'.

    % Check if groups have more than one observation
    if length(group1) > 1 && length(group2) > 1  % Verifica se ambos os grupos têm mais de uma observação
        [~, p_value] = ttest2(group1, group2);  % Realiza o teste t para comparar as médias.
        minus_log10_p = -log10(p_value);  % Calcula o logaritmo negativo na base 10 do p-valor para visualização.
    else
        minus_log10_p = NaN;  % Se o teste não puder ser realizado, define como NaN.
    end

    % Calculate means and standard deviations of each group
    mean_group1 = mean(group1);  % Média do grupo 'Controle'.
    mean_group2 = mean(group2);  % Média do grupo 'Parkinson'.
    std_group1 = std(group1);  % Desvio padrão do grupo 'Controle'.
    std_group2 = std(group2);  % Desvio padrão do grupo 'Parkinson'.
    % Add -log10(p-value) at the top of the boxplot
    ylim_current = ylim();  % Get current y-axis limits.
    y_position = ylim_current(2) + (ylim_current(2) - ylim_current(1)) * 0.15;  % Position for the text above the boxplot.
    text(mean(xlim()), y_position, sprintf('-log_{10}(p) = %.2f', minus_log10_p), ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontWeight', 'bold');

    % Adjust ylim to accommodate the text
    ylim([ylim_current(1), y_position + (ylim_current(2) - ylim_current(1)) * 0.2]);

    % Add means and standard deviations above each boxplot
    hold on;  % Keep current plot to overlay next elements.
    positions = get(gca, 'XTick');  % Get x-axis tick positions.
    y_offset = (ylim_current(2) - ylim_current(1)) * 0.05;  % Calculate offset to avoid overlap.

    % Add mean and standard deviation for 'Control' group (index 1)
    text(positions(1), mean_group1 + y_offset, sprintf('%.2f ± %.2f', mean_group1, std_group1), ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Color', 'blue', 'FontWeight', 'bold');

    % Add mean and standard deviation for 'Parkinson' group (index 2)
    text(positions(2), mean_group2 + y_offset, sprintf('%.2f ± %.2f', mean_group2, std_group2), ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Color', 'red', 'FontWeight', 'bold');
    hold off;  % Release current plot.

    % Adjust x-axis labels to show 'Control' and 'Parkinson'
    set(gca, 'XTickLabel', unique_groups);  % Set x-axis labels with group names.
end

% Add general title to the figure
annotation('textbox', [0.5, 0.98, 0, 0], 'String', 'Boxplots of Selected Variables by Category', ...
           'FitBoxToText', 'on', 'HorizontalAlignment', 'center', 'FontSize', 12, 'LineStyle', 'none');

% Clear variables and console
clearvars;
pkg load statistics;

% Load data from CSV
data = csvread('parkinson_mdvp.csv', 1, 0); % Ignore header

% Logarithmic normalization
data(:, 1:end-1) = log(data(:, 1:end-1) + 1);

% Scaling
data(:, 1:end-1) = (data(:, 1:end-1) - mean(data(:, 1:end-1))) ./ std(data(:, 1:end-1));

% Define columns for the two models
x1 = data(:, 9);  % MDVP:Shimmer(dB)
y1 = data(:, end);  % Status (1 for Parkinson, 0 for Control)

x2 = data(:, 10);  % MDVP:APQ
y2 = data(:, end);

% Call the double regression function
[coef1, coef2, Sr1, Sr2, r2_1, r2_2, syx1, syx2] = regressao_dupla(x1, y1, x2, y2);

% Display the results of Analysis 3.1
disp('===========================================================');
disp('                         Analysis 3.1             ');
disp('===========================================================');
disp('Coefficients of Model 1 (MDVP:Shimmer(dB)):');
disp(coef1);
disp('Sum of Residuals Model 1:'), disp(Sr1);
disp('Coefficient of Determination r² Model 1:'), disp(r2_1);
disp('Standard Deviation of Residuals Model 1:'), disp(syx1);

disp('-----------------------------------');
disp('Coefficients of Model 2 (MDVP:APQ):');
disp(coef2);
disp('Sum of Residuals Model 2:'), disp(Sr2);
disp('Coefficient of Determination r² Model 2:'), disp(r2_2);
disp('Standard Deviation of Residuals Model 2:'), disp(syx2);

% Comparison of models
disp('-----------------------------------');
if r2_1 > r2_2 % Verifica se o coeficiente de determinação do Modelo 1 é maior que o do Modelo 2
    disp('Model 1 is better based on the coefficient of determination r².');
else
    disp('Model 2 is better based on the coefficient of determination r².');
end
disp('===========================================================');

% Additional analysis (Analysis 3.2)
y3 = data(:, end);  % Status (1 for Parkinson, 0 for Control)
n = length(y3);  % Number of data points

% Define as combinações a serem testadas
combinations = {
    [x1, x2], 'MDVP:Shimmer(dB) e MDVP:APQ'; % Combinação de variáveis independentes com o nome do modelo
};


best_r2 = -Inf;  % Inicializa o melhor r² como um valor muito baixo
best_model = '';  % Inicializa a variável para o melhor modelo
best_coef = [];   % Inicializa a variável para os coeficientes do melhor modelo
best_Sr = [];     % Inicializa a variável para a soma dos resíduos do melhor modelo
best_syx = [];    % Inicializa a variável para o desvio padrão dos resíduos do melhor modelo


for i = 1:size(combinations, 1) % Loop para testar cada combinação de modelos
    X = combinations{i, 1};  % Variáveis independentes
    model_name = combinations{i, 2};  % Nome do modelo

    % Adiciona uma coluna de 1's para o intercepto
    X = [X, ones(size(X, 1), 1)];

    % Calcula os coeficientes usando decomposição LU
    coef = lu_solve(X, y3); % Chama a função para resolver o sistema usando decomposição LU
    y_pred = X * coef;  % Prediz os valores

    % Calcula a soma dos resíduos
    Sr = sum((y3 - y_pred).^2); % Soma dos quadrados dos resíduos

    % Calcula a soma total dos quadrados
    St = sum((y3 - mean(y3)).^2); % Soma total dos quadrados

    % Calcula o coeficiente de determinação r²
    r2 = 1 - (Sr / St); % Coeficiente de determinação

    % Calcula o desvio padrão dos resíduos
    syx = sqrt(Sr / (length(y3) - size(X, 2))); % Desvio padrão dos resíduos

    % Verifica se este modelo é o melhor até agora
    if r2 > best_r2
        best_r2 = r2; % Atualiza o melhor r²
        best_model = model_name; % Atualiza o melhor modelo
        best_coef = coef; % Atualiza os coeficientes do melhor modelo
        best_Sr = Sr; % Atualiza a soma dos resíduos do melhor modelo
        best_syx = syx; % Atualiza o desvio padrão dos resíduos do melhor modelo
    end
end

% Display the results of Analysis 3.2
disp('===========================================================');
disp('                         Analysis 3.2             ');
disp('===========================================================');
disp('===========================================================');
disp(['The best model is: ' best_model]);
disp('Coefficients of the Best Model (a0,3, a1,3, a2,3):');
disp(['a0,3: ', num2str(best_coef(end))]); % Intercept
disp(['a1,3: ', num2str(best_coef(1))]);   % Coefficient of x1
disp(['a2,3: ', num2str(best_coef(2))]);   % Coefficient of x2
disp('Sum of Residuals of the Best Model:'), disp(best_Sr);
disp('Coefficient of Determination r^2 of the Best Model:'), disp(best_r2);
disp('Standard Deviation of Residuals of the Best Model:'), disp(best_syx);
disp('===========================================================');

