function alpha_string = nz_vector_elements_string(alpha)
%NZ_VECTOR_ELEMENTS_STRING returns a string ending with a newline containing
%an entry "a{i}=val; " for all nonzero elements of a vector.
    alpha_string = "";
    for ii = 1 : length(alpha)
        if alpha(ii) ~= 0
           alpha_string = strcat(alpha_string, sprintf("a{%d} = %.2f; ", ii, alpha(ii)) );
        end
    end
    alpha_string = strcat(alpha_string, "\n");
end