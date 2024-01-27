function gradf=standard_fun_grad(z,fun_grad,mu,v)

    % given a function gradient, this returns the gradient of the
    % function standardized wrt mu and v.

    % parse inputs
    parser=inputParser;
    parser.KeepUnmatched=true;
    addRequired(parser,'z',@iscolumn);
    addRequired(parser,'fun_grad',@(x) isa(x,'function_handle'));
    addRequired(parser,'mu',@isnumeric);
    addRequired(parser,'v',@isnumeric);

    parse(parser,z,fun_grad,mu,v);

    S=sqrtm(v);
    gradf=S'*fun_grad(S*z+mu);

end