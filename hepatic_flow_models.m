function models = hepatic_flow_models()
%HEPATIC_FLOW_MODELS  Definitions of the dual-input liver models that are fitted.
%
%   models = HEPATIC_FLOW_MODELS() returns a struct array with one element
%   per model variant. All variants are based on the dual-input
%   one-tissue compartment model for [15O]water described by
%   Kudomi et al., Eur J Nucl Med Mol Imaging 2008 (see HEPATIC_FLOW_MODEL).
%
%   The full parameter set is
%       F      total hepatic blood flow (ml/cm3/min)
%       PF     portal fraction of the total flow, Fp/F
%       VT     distribution volume of water (ml/cm3)
%       kg     dispersion rate constant of the gut/spleen (1/min)
%       delay  delay of the portal-vein input relative to the arterial
%              input (min)
%       VB     blood volume fraction
%
%   The variants differ in which of kg, delay and VB are estimated, fixed,
%   coupled (kg = 6*delay) or omitted:
%
%   id                   free parameters          kg            delay              VB
%   kgfix                F PF VT delay VB         data.kg       fitted             fitted
%   kgdelayfix           F PF VT VB               data.kg       data.delay         fitted
%   full                 F PF VT kg delay VB      fitted        fitted             fitted
%   coupled              F PF VT delay VB         6*delay       fitted             fitted
%   coupled_fixed        F PF VT VB               6*delay       data.delaycoupled  fitted
%   coupled_novb         F PF VT delay            6*delay       fitted             0
%   coupled_fixed_novb   F PF VT                  6*delay       data.delaycoupled_novb  0
%   kgdelayfix_novb      F PF VT                  data.kg       data.delay         0
%
%   Fields
%     id     identifier passed to HEPATIC_FLOW_MODEL / FIT_HEPATIC_FLOW
%     title  section heading in the results (.csv) file
%     label  name used for the AIC-preferred model in the results file
%     ub     upper bounds of the free parameters (lower bounds are 0)
%
%   The order of the models determines the tie-breaking when the
%   preferred model (lowest AIC) is selected.

models = struct( ...
    'id',    {'kgfix', 'kgdelayfix', 'full', 'coupled', 'coupled_fixed', ...
              'coupled_novb', 'coupled_fixed_novb', 'kgdelayfix_novb'}, ...
    'title', {'Kudomi kg fixed', 'Kudomi kg delay fixed', 'Kudomi', ...
              'Kudomi kg delay coupled', 'Kudomi kg delay coupled and fixed', ...
              'Kudomi kg delay coupled no vb', ...
              'Kudomi kg delay coupled and fixed no vb', ...
              'Kudomi kg delay fixed no vb'}, ...
    'label', {'Kudomi kg fixed', 'Kudomi kg delay fixed', 'Kudomi full', ...
              'Kudomi kg delay coupled', 'Kudomi kg delay coupled and fixed', ...
              'Kudomi kg delay coupled no vb', ...
              'Kudomi kg delay coupled and fixed no vb', ...
              'Kudomi kg delay fixed no vb'}, ...
    'ub',    {[4 1 1.3 0.6 1], ...      % F PF VT delay VB
              [4 1 1.3 1], ...          % F PF VT VB
              [4 1 1.3 3 0.6 1], ...    % F PF VT kg delay VB
              [4 1 1.3 0.6 1], ...      % F PF VT delay VB
              [4 1 1.3 1], ...          % F PF VT VB
              [4 1 1.3 0.6], ...        % F PF VT delay
              [4 1 1.3], ...            % F PF VT
              [4 1 1.3]});              % F PF VT
end
