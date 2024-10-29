Análise 1.1

Foram escolhidos os pares MDVP:Fo(Hz) e MDVP:Flo(Hz) para o grupo sem Parkinson, pois ambos mostraram média e dispersão mais altas neste grupo, indicando maior estabilidade na frequência vocal. Para o grupo com Parkinson, 
foi selecionado MDVP:Shimmer(dB) e MDVP:APQ , que apresentaram médias e dispersão mais elevadas nos pacientes com a doença, refletindo a instabilidade vocal característica. Esses pares evidenciam diferenças importantes entre os grupos.



Análise 2.2

Média e Valor Absoluto:

Para métricas de frequência, como Fo(Hz), Fhi(Hz) e Flo(Hz), a média dos valores para o grupo de controle é consistentemente maior que para o grupo Parkinson, sugerindo uma maior capacidade de controle e variação da frequência vocal em indivíduos sem Parkinson.
Em variáveis relacionadas à instabilidade e inconsistência vocal, como Jitter(%), Shimmer e APQ, as médias são geralmente mais altas no grupo Parkinson. Isso indica que a voz de indivíduos com Parkinson apresenta maior variação e instabilidade, características associadas ao tremor vocal e à rigidez muscular.

Dispersão e Desvio Padrão:

Em termos de dispersão, o grupo de Parkinson apresenta maior desvio padrão para quase todas as métricas, especialmente para Jitter e Shimmer, métricas associadas à instabilidade da frequência e amplitude vocal. Esse aumento na dispersão reflete a variabilidade de sintomas vocais entre os pacientes, como flutuações na capacidade de manter uma frequência estável.
O maior desvio padrão no grupo Parkinson para as métricas de jitter e shimmer indica uma maior variabilidade da intensidade e frequência de voz, o que é esperado devido às dificuldades motoras que afetam o controle da fala.

Interpretação Clínica:

A média mais baixa em variáveis de frequência para o grupo Parkinson sugere que a voz desses pacientes pode ser menos modulada e com menor amplitude de variação.
O aumento no desvio padrão indica uma maior dificuldade em manter consistência vocal, o que é característico dos sintomas motores de Parkinson, como a rigidez e o tremor que dificultam o controle preciso dos músculos vocais.

Anáise 3.3

Baixo Coeficiente de Determinação (R²): O valor de 𝑅2 = 0.1414 indica que o modelo explica apenas cerca de 14,14% da variabilidade dos
dados de resposta. Em modelos de regressão, um 𝑅2 mais próximo de 1 sugere que o modelo é eficaz para capturar a variação da variável dependente.
No entanto, com um valor de 0,1414, a maior parte da variabilidade no status de Parkinson não é explicada por essa combinação de variáveis.

Desvio Padrão dos Resíduos (syx): O desvio padrão dos resíduos, 𝑠𝑦𝑥=0.4023, sugere que as previsões do modelo têm uma variação relativamente alta
em torno da linha de regressão. Isso indica que o modelo apresenta um erro considerável ao prever novos dados, o que reduz sua confiabilidade para diagnósticos.


Codigo para calcular o EQM -
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


% Predictions for Model 1
y1_pred = [ones(length(x1), 1), x1] * coef1; % Predições para o Modelo 1

% Predictions for Model 2
y2_pred = [ones(length(x2), 1), x2] * coef2; % Predições para o Modelo 2

% Calculate MSE for both models
mse1 = calcular_mse(y1, y1_pred); % MSE para o Modelo 1
mse2 = calcular_mse(y2, y2_pred); % MSE para o Modelo 2

% Display the results of Analysis 3.1
disp('===========================================================');
disp('                         Analysis 3.1             ');
disp('===========================================================');
disp('Coefficients of Model 1 (MDVP:Shimmer(dB)):'), disp(coef1);
disp('Sum of Residuals Model 1:'), disp(Sr1);
disp('Coefficient of Determination r² Model 1:'), disp(r2_1);
disp('Standard Deviation of Residuals Model 1:'), disp(syx1);
disp('MSE for Model 1:'), disp(mse1); % MSE para o Modelo 1

disp('-----------------------------------');
disp('Coefficients of Model 2 (MDVP:APQ):'), disp(coef2);
disp('Sum of Residuals Model 2:'), disp(Sr2);
disp('Coefficient of Determination r² Model 2:'), disp(r2_2);
disp('Standard Deviation of Residuals Model 2:'), disp(syx2);
disp('MSE for Model 2:'), disp(mse2); % MSE para o Modelo 2

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


for i = 1
    X = combinations{1};  % Variáveis independentes
    model_name = combinations{i, 2};  % Nome do modelo

    % Adiciona uma coluna de 1's para o intercepto
    X = [X, ones(size(X, 1), 1)];

    % Verifica se X e y3 têm o mesmo número de linhas
    assert(size(X, 1) == length(y3), 'X e y3 devem ter o mesmo número de linhas');

    % Calcula os coeficientes usando decomposição LU
    coef = solve_lu(X' * X, X' * y3); % Usa solve_lu para resolver o sistema linear
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


% Clear variables and console
clearvars;
pkg load statistics;

% Load data from CSV
data = csvread('parkinson_mdvp.csv', 1, 0); % Ignore header

% Logarithmic normalization
data(:, 1:end-1) = log(data(:, 1:end-1) + 1);

% Scaling
data(:, 1:end-1) = (data(:, 1:end-1) - mean(data(:, 1:end-1))) ./ std(data(:, 1:end-1));

% Extrair as colunas relevantes e o 'status' real
X1 = data(:, 9);  % Coluna MDVP:Shimmer(dB)
X2 = data(:, 10); % Coluna MDVP:APQ
y_real = data(:, 11);  % Coluna status (valores reais)

% Coeficientes do melhor modelo
a0 = 0.75385;
a1 = 0.10122;
a2 = 0.063025;

% Calcular as previsões (y_pred) com base no modelo
y_pred = a0 + a1 * X1 + a2 * X2;

% Função para calcular o MSE
function mse_value = calcular_mse(y_real, y_pred)
    if length(y_real) != length(y_pred)
        error('Os vetores y_real e y_pred devem ter o mesmo tamanho.');
    end
    n = length(y_real);
    mse_value = sum((y_real - y_pred).^2) / n;
    fprintf('MSE: %.4f\n', mse_value);
end

% Calcular e exibir o MSE para o melhor modelo
mse = calcular_mse(y_real, y_pred);
