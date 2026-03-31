function ret_val = ifelsef(if_val, then_val, else_val)
    if if_val
        ret_val = then_val;
    else
        ret_val = else_val;
    end
end