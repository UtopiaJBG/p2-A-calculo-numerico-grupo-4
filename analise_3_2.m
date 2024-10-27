% Limpar variáveis e console
clearvars;
clc;
pkg load statistics;

% Carregar o resultado do 3.1
run('analise_3_1.m');

% Definir y3 como a variável de saída (status de Parkinson)
y3 = data(:, end);  % Status (1 para Parkinson, 0 para Sem Parkinson)

% Número de dados
n = length(y3);

% Lista de combinações possíveis de pares de variáveis
combinacoes = {
    [x1, x2], 'MDVP:Fo(Hz) e MDVP:Shimmer(dB)';
    [x1, y1], 'MDVP:Fo(Hz) e MDVP:Flo(Hz)';
    [x1, y2], 'MDVP:Fo(Hz) e MDVP:APQ';
    [x2, y1], 'MDVP:Shimmer(dB) e MDVP:Flo(Hz)';
    [x2, y2], 'MDVP:Shimmer(dB) e MDVP:APQ';
    [y1, y2], 'MDVP:Flo(Hz) e MDVP:APQ'
};

% Inicializar variáveis para armazenar o melhor modelo
melhor_r2 = -Inf;
melhor_modelo = '';
melhor_coef = [];
melhor_Sr = 0;
melhor_syx = 0;

% Loop para cada combinação de pares de variáveis
for i = 1:length(combinacoes)
    % Obter o par de variáveis atual e o nome do modelo
    par_vars = combinacoes{i, 1};
    nome_modelo = combinacoes{i, 2};

    % Criar a matriz para a regressão múltipla
    A = [par_vars, ones(n, 1)];

    % Coeficientes do modelo de regressão múltipla
    coef = lu_solve(A, y3);

    % Valores preditos para o modelo de regressão múltipla
    y3_pred = A * coef;

    % Soma dos resíduos para o modelo
    Sr = sum((y3 - y3_pred).^2);

    % Soma total dos quadrados para y3
    St = sum((y3 - mean(y3)).^2);

    % Coeficiente de determinação para o modelo
    r2 = 1 - (Sr / St);

    % Desvio padrão dos resíduos para o modelo
    syx = sqrt(Sr / (n - 3));

    % Exibir os resultados da combinação atual
    disp(['Coeficientes do Modelo de Regressão Múltipla (' nome_modelo '):']), disp(coef');
    disp(['Soma dos Resíduos para ' nome_modelo ':']), disp(Sr);
    disp(['Coeficiente de Determinação r^2 para ' nome_modelo ':']), disp(r2);
    disp(['Desvio Padrão dos Resíduos para ' nome_modelo ':']), disp(syx);
    disp('------------------------------------------------------');

    % Verificar se este modelo é o melhor até agora
    if r2 > melhor_r2
        melhor_r2 = r2;
        melhor_modelo = nome_modelo;
        melhor_coef = coef;
        melhor_Sr = Sr;
        melhor_syx = syx;
    end
end

% Exibir o melhor modelo encontrado
disp('==========================================');
disp(['O melhor modelo é: ' melhor_modelo]);
disp('Coeficientes do Melhor Modelo:'), disp(melhor_coef');
disp('Soma dos Resíduos do Melhor Modelo:'), disp(melhor_Sr);
disp('Coeficiente de Determinação r^2 do Melhor Modelo:'), disp(melhor_r2);
disp('Desvio Padrão dos Resíduos do Melhor Modelo:'), disp(melhor_syx);

% Função auxiliar para resolver o sistema linear (mesma usada no 3.1)
function x = lu_solve(A, b)
    [L, U, P] = lu(A);  % Decomposição LU com pivoteamento
    y = L \ (P * b);    % Resolver L * y = P * b
    x = U \ y;          % Resolver U * x = y
end

