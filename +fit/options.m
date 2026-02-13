function fit_options = options(varargin)

fit_options.display = true;
fit_options.imax = 1e9;
fit_options.gradeps = 0;
fit_options.lambda = eps;
fit_options.upFactor = 3;
fit_options.downFactor = 2;
fit_options.A = [];
fit_options.a = [];
fit_options.w = [];

for n = 1:2:length(varargin)
    fit_options.(varargin{n}) = varargin{n+1};
end