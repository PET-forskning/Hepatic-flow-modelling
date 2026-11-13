function c_tis = hepatic_flow_model(p, data, model)
%HEPATIC_FLOW_MODEL  Dual-input one-tissue model for hepatic [15O]water PET.
%
%   c_tis = HEPATIC_FLOW_MODEL(p, data, model) returns the model liver
%   time-activity curve, multiplied by sqrt(data.weights) so that it can be
%   used directly in a weighted least-squares fit.
%
%   Model (Kudomi et al., Eur J Nucl Med Mol Imaging 2008):
%
%     Portal-vein input, derived from the arterial input Ca(t):
%       Cp(t) = kg * exp(-kg*t) (*) Ca(t - delay)
%
%     Liver tissue:
%       Ctis(t) = (1-VB) * [ Fa*Ca(t) + Fp*Cp(t) ] (*) exp(-F/VT * t)
%                 + VB * [ (1-PF)*Ca(t) + PF*Cp(t) ]
%
%     with F = Fa + Fp, PF = Fp/F, and (*) denoting convolution.
%
%   Inputs
%     p      vector of free parameters for the chosen model variant
%            (see HEPATIC_FLOW_MODELS)
%     data   struct with fields
%              t_tac    frame mid-times (min)
%              c_pl     arterial input function (kBq/ml), sampled at t_tac
%              weights  fit weights per frame
%              kg, delay, delaycoupled, delaycoupled_novb
%                       fixed values, used only by the variants that need them
%     model  model id (see HEPATIC_FLOW_MODELS)

[F, PF, VT, kg, delay, VB, has_vb] = expand_parameters(p, data, model);

% Portal-vein input: delayed and dispersed arterial input
cp = convexvarc(data.t_tac, [kg kg 0], data.t_tac + delay, data.c_pl);

% Tissue contributions from the hepatic artery and the portal vein
c_ha = convexvarc(data.t_tac, [F*(1-PF) F/VT 0], data.t_tac, data.c_pl);
c_pv = convexvarc(data.t_tac, [F*PF     F/VT 0], data.t_tac, cp);

if has_vb
    % Blood-volume contributions
    vb_a  = (1-PF) * data.c_pl;
    vb_p  = PF * cp;
    c_tis = (1-VB)*(c_ha + c_pv) + VB*(vb_a + vb_p);
else
    c_tis = c_ha + c_pv;
end

c_tis = c_tis .* sqrt(data.weights);

end


function [F, PF, VT, kg, delay, VB, has_vb] = expand_parameters(p, data, model)
% Map the free-parameter vector of a model variant to the full parameter set.
F  = p(1);
PF = p(2);
VT = p(3);
has_vb = true;
switch model
    case 'full'
        kg = p(4);  delay = p(5);  VB = p(6);
    case 'kgfix'
        kg = data.kg;  delay = p(4);  VB = p(5);
    case 'kgdelayfix'
        kg = data.kg;  delay = data.delay;  VB = p(4);
    case 'kgdelayfix_novb'
        kg = data.kg;  delay = data.delay;  VB = 0;  has_vb = false;
    case 'coupled'
        delay = p(4);  kg = delay*6;  VB = p(5);
    case 'coupled_novb'
        delay = p(4);  kg = delay*6;  VB = 0;  has_vb = false;
    case 'coupled_fixed'
        delay = data.delaycoupled;  kg = delay*6;  VB = p(4);
    case 'coupled_fixed_novb'
        delay = data.delaycoupled_novb;  kg = delay*6;  VB = 0;  has_vb = false;
    otherwise
        error('hepatic_flow_model:unknownModel', 'Unknown model "%s".', model);
end
end
