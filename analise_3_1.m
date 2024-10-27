% Limpar variáveis e console
clearvars;
clc;
pkg load statistics;


function [coef1, coef2, Sr1, Sr2, r2_1, r2_2, syx1, syx2] = regressao_dupla(x1, y1, x2, y2)
    % Número de dados
    n = length(x1);

    % Modelo 1: y1 = a0,1 + a1,1 * x1
    A1 = [x1, ones(n, 1)];
    coef1 = lu_solve(A1, y1);  % Coeficientes

    y1_pred = A1 * coef1;  % Valores preditos
    Sr1 = sum((y1 - y1_pred).^2);  % Soma dos resíduos para o modelo 1
    St1 = sum((y1 - mean(y1)).^2);  % Soma total dos quadrados para y1
    r2_1 = 1 - (Sr1 / St1);  % Coeficiente de determinação para o modelo 1
    syx1 = sqrt(Sr1 / (n - 2));  % Desvio padrão dos resíduos para o modelo 1

    % Modelo 2: y2 = a0,2 + a1,2 * x2
    A2 = [x2, ones(n, 1)];
    coef2 = lu_solve(A2, y2);  % Coeficientes

    y2_pred = A2 * coef2;  % Valores preditos
    Sr2 = sum((y2 - y2_pred).^2);  % Soma dos resíduos para o modelo 2
    St2 = sum((y2 - mean(y2)).^2);  % Soma total dos quadrados para y2
    r2_2 = 1 - (Sr2 / St2);  % Coeficiente de determinação para o modelo 2
    syx2 = sqrt(Sr2 / (n - 2));  % Desvio padrão dos resíduos para o modelo 2

    % Comparação dos modelos
    if r2_1 > r2_2
        disp('Modelo 1 é melhor com base no coeficiente de determinação r^2.');
    else
        disp('Modelo 2 é melhor com base no coeficiente de determinação r^2.');
    end
end

% Função auxiliar para resolver o sistema linear

function x = lu_solve(A, b)

    [L, U, P] = lu(A);  % Decomposição LU com pivoteamento
    y = L \ (P * b);    % Resolver L * y = P * b
    x = U \ y;          % Resolver U * x = y
end


% Carregar os dados do CSV


data = csvread('parkinson_mdvp.csv');
data(1, :) = [];  % Remover cabeçalho, se necessário

% Normalização logarítmica
data(:, 1:end-1) = log(data(:, 1:end-1) + 1);

% Escalonamento
data(:, 1:end-1) = (data(:, 1:end-1) - mean(data(:, 1:end-1))) ./ std(data(:, 1:end-1));

% Definir colunas para os dois modelos
x1 = data(:, 1);  % MDVP:Fo(Hz)
y1 = data(:, 3);  % MDVP:Flo(Hz)

x2 = data(:, 8);  % MDVP:Shimmer(dB)
y2 = data(:, 10); % MDVP:APQ

% Chamar a função de regressão dupla
[coef1, coef2, Sr1, Sr2, r2_1, r2_2, syx1, syx2] = regressao_dupla(x1, y1, x2, y2);

% Exibir os resultados
disp('Coeficientes do Modelo 1 (MDVP:Fo(Hz) vs MDVP:Flo(Hz)):'), disp(coef1);
disp('Soma dos Resíduos Modelo 1:'), disp(Sr1);
disp('Coeficiente de Determinação r^2 Modelo 1:'), disp(r2_1);
disp('Desvio Padrão dos Resíduos Modelo 1:'), disp(syx1);

disp('Coeficientes do Modelo 2 (MDVP:Shimmer(dB) vs MDVP:APQ):'), disp(coef2);
disp('Soma dos Resíduos Modelo 2:'), disp(Sr2);
disp('Coeficiente de Determinação r^2 Modelo 2:'), disp(r2_2);
disp('Desvio Padrão dos Resíduos Modelo 2:'), disp(syx2);

% Comparação dos modelos
if r2_1 > r2_2
    disp('Modelo 1 é melhor com base no coeficiente de determinação r^2.');
else
    disp('Modelo 2 é melhor com base no coeficiente de determinação r^2.');
end



