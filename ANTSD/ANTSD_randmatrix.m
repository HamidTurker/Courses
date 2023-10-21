function [matrix] = ANTSD_randmatrix(n_rows,n_cols)

matrix = zeros(n_rows*n_cols,3);
random_matrix = rand(n_rows,n_cols);
counter = 0;
for r = 1:n_rows
    for c = 1:n_cols

        % Save current row and column values to output matrix
        counter = counter + 1;
        matrix(counter,1) = r; matrix(counter,2) = c;

        % Adjust ordinality string
        if (c == 1)
            ord_c = "st";
        elseif (c == 2)
            ord_c = "nd";
        elseif (c == 3)
            ord_c = "rd";
        else
            ord_c = "th";
        end

        if (r == 1)
            ord_r = "st";
        elseif (r == 2)
            ord_r = "nd";
        elseif (r == 3)
            ord_r = "rd";
        else
            ord_r = "th";
        end

        % Check if greater than 0.5 and save result to output matrix
        if (random_matrix(r,c) > .5) % If the given row-column is > 0.5
            disp(strcat("The ", num2str(r), ord_r, " row and ", num2str(c), ord_c," column has a value of ", num2str(random_matrix(r,c))," and is bigger than 0.5."))
            matrix(counter,3) = 1;
        else
            disp(strcat("The ", num2str(r), ord_r, " row and ", num2str(c), ord_c," column has a value of ", num2str(random_matrix(r,c))," and is not bigger than 0.5."))
        end
    end
end

end