function x = lu_solve(A, b)
    % Decomposição LU
    [L, U, P] = lu(A); % P é a matriz de permutação

    % Resolver o sistema Ly = Pb
    y = L \ (P * b);

    % Resolver o sistema Ux = y
    x = U \ y;
end
