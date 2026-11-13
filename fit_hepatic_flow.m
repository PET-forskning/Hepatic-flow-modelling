function r = fit_hepatic_flow(model, p0, tac, data)
%FIT_HEPATIC_FLOW  Weighted least-squares fit of one hepatic flow model variant.
%
%   r = FIT_HEPATIC_FLOW(model, p0, tac, data) fits the model variant with
%   id MODEL (see HEPATIC_FLOW_MODELS) to the tissue time-activity curve TAC,
%   starting from the initial guess P0 (free parameters only). DATA is
%   described in HEPATIC_FLOW_MODEL.
%
%   Output struct r
%     p       estimated free parameters
%     p_sd    standard errors of p (from the 95% confidence intervals)
%     fit     fitted (unweighted) model curve, sampled at data.t_tac
%     row     1x16 vector with the full parameter set and standard errors,
%             [F SD PF SD VT SD kg SD delay SD VB SD Fp SD Fa SD];
%             parameters that are fixed or omitted have SD = 0
%     aic     Akaike information criterion
%     sc      Schwarz criterion
%
%   Requires the Optimization Toolbox (lsqcurvefit) and the Statistics and
%   Machine Learning Toolbox (nlparci).

models = hepatic_flow_models();
m  = models(strcmp({models.id}, model));
lb = zeros(size(m.ub));
ub = m.ub;

options = optimset(optimset('lsqcurvefit'), 'MaxFunEvals', 1E6, ...
    'TolFun', 1E-6, 'TolX', 1E-6, 'MaxIter', 1E6);

% Weighted fit: both model and data are multiplied by sqrt(weights)
model_fun = @(p, d) hepatic_flow_model(p, d, model);
[p, resnorm, residual, ~, ~, ~, jacobian] = lsqcurvefit(model_fun, p0, ...
    data, tac .* sqrt(data.weights), lb, ub, options);

% Standard errors from the 95% confidence intervals
ci   = nlparci(p, residual, 'jacobian', jacobian);
p_sd = abs(ci(:,2) - ci(:,1))' / (2*1.95996);

% Information criteria; k = number of fitted parameters, n = number of frames
n = length(data.t_tac);
k = length(p);
r.aic = n*log(resnorm) + 2*k;
r.sc  = n*log(resnorm) + log(n)*k;

% Fitted curve without weighting (for plotting)
data_unweighted = data;
data_unweighted.weights = ones(size(data.weights));
r.fit = hepatic_flow_model(p, data_unweighted, model);

r.p    = p;
r.p_sd = p_sd;

% Full parameter set [kg delay VB] with standard errors (0 if not fitted)
switch model
    case 'full'
        kg = [p(4) p_sd(4)];  delay = [p(5) p_sd(5)];  vb = [p(6) p_sd(6)];
    case 'kgfix'
        kg = [data.kg 0];  delay = [p(4) p_sd(4)];  vb = [p(5) p_sd(5)];
    case 'kgdelayfix'
        kg = [data.kg 0];  delay = [data.delay 0];  vb = [p(4) p_sd(4)];
    case 'kgdelayfix_novb'
        kg = [data.kg 0];  delay = [data.delay 0];  vb = [0 0];
    case 'coupled'
        kg = [p(4) p_sd(4)]*6;  delay = [p(4) p_sd(4)];  vb = [p(5) p_sd(5)];
    case 'coupled_novb'
        kg = [p(4) p_sd(4)]*6;  delay = [p(4) p_sd(4)];  vb = [0 0];
    case 'coupled_fixed'
        kg = [data.delaycoupled*6 0];  delay = [data.delaycoupled 0];
        vb = [p(4) p_sd(4)];
    case 'coupled_fixed_novb'
        kg = [data.delaycoupled_novb*6 0];  delay = [data.delaycoupled_novb 0];
        vb = [0 0];
end

% Portal (Fp = F*PF) and arterial (Fa = F*(1-PF)) flow, with error propagation
Fp    = p(1)*p(2);
Fp_sd = p(1)*p(2)*sqrt((p_sd(1)/p(1))^2 + (p_sd(2)/p(2))^2);
Fa    = p(1)*(1-p(2));
Fa_sd = p(1)*(1-p(2))*sqrt((p_sd(1)/p(1))^2 + (p_sd(2)/(1-p(2)))^2);

r.row = [p(1) p_sd(1) p(2) p_sd(2) p(3) p_sd(3) kg delay vb Fp Fp_sd Fa Fa_sd];

end
