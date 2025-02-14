function Complex_Check(matrix)
    Name = inputname(1);
    if any(~isreal(matrix))
        [row, col] = find(~isreal(matrix));
        fprintf('Complex number found at %s(%d, %d)\n', Name, row(1), col(1));
    end
end
