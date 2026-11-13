function analyze_water_hepatic_flow(aorta_file, voi_file, aorta_index)
%ANALYZE_WATER_HEPATIC_FLOW  Segmental hepatic blood flow from dynamic [15O]water PET.
%
%   ANALYZE_WATER_HEPATIC_FLOW() asks for the input files with file dialogs.
%   ANALYZE_WATER_HEPATIC_FLOW(aorta_file, voi_file, aorta_index) runs the
%   analysis without dialogs.
%
%   Inputs
%     aorta_file   VOIager statistics file (*_stat.txt) with the arterial
%                  (aorta) time-activity curve used as input function
%     voi_file     HERMES Excel export (*.xlsx) with the liver segment VOIs
%     aorta_index  index of the arterial curve in aorta_file (if omitted,
%                  a menu is shown)
%
%   Both files must contain the same dynamic frames.
%
%   Analysis
%     1. A whole-liver TAC is formed as the volume-weighted mean of all
%        VOIs in voi_file.
%     2. Eight variants of the dual-input one-tissue model (Kudomi et al.,
%        Eur J Nucl Med Mol Imaging 2008; see HEPATIC_FLOW_MODELS) are
%        fitted to the whole-liver TAC. The portal-vein delay from the
%        'kg fixed' fit (kg = 1.22 /min) is used as the fixed delay in the
%        variants with fixed delay, and the delays from the 'kg delay
%        coupled' fits are used in the corresponding 'coupled and fixed'
%        variants.
%     3. All variants are fitted to each VOI, with the whole-liver estimates
%        as initial guesses.
%     4. For each VOI the preferred model is the one with the lowest AIC.
%
%   Frames are weighted by frame duration times the decay factor of 15O
%   (half-life 2.03 min); frames with mid-time later than 6 min are
%   excluded (weight 0).
%
%   Output (written next to voi_file)
%     <voi_file>_hepaticflow.pdf   input curves and fitted TACs
%     <voi_file>_hepaticflow.csv   parameter estimates with standard errors,
%                                  AIC and Schwarz criterion for every
%                                  model and VOI, and the preferred model
%
%   Requires the Optimization Toolbox (lsqcurvefit) and the Statistics and
%   Machine Learning Toolbox (nlparci).

%% Settings
KG_FIXED        = 1.22;   % fixed gut dispersion rate constant kg (1/min)
T_HALF_O15      = 2.03;   % physical half-life of 15O (min)
T_MAX_FIT       = 6;      % frames with mid-time > T_MAX_FIT (min) get weight 0

% Initial guesses for the whole-liver fits
p0_kgfix      = [1 0.8 0.83 0.30 0.2];       % F, PF, VT, delay, VB
p0_kgdelayfix = [1 0.8 0.83 0.2];            % F, PF, VT, VB
p0_full       = [1 0.8 0.83 1.22 0.30 0.2];  % F, PF, VT, kg, delay, VB

%% Select input files
if nargin < 1 || isempty(aorta_file)
    [file, path] = uigetfile('*aorta*.txt', 'Open VOIager TAC file (*_stat.txt) with arterial input');
    if isequal(file, 0), return; end
    aorta_file = fullfile(path, file);
end
if nargin < 2 || isempty(voi_file)
    [aorta_path, aorta_name] = fileparts(aorta_file);
    subject = aorta_name(1:min(6, end));   % first characters of the file name
    [file, path] = uigetfile(fullfile(aorta_path, [subject '*.xlsx']), 'Open HERMES TAC file');
    if isequal(file, 0), return; end
    voi_file = fullfile(path, file);
end

%% Arterial input function
aorta = read_roi_voiager(aorta_file);
if nargin < 3 || isempty(aorta_index)
    aorta_index = menu('Choose AA TAC:', aorta.name);
end

%% Liver VOIs
data = struct();
[roi, data] = read_voi_HERMES(voi_file, data);
nvoi = length(roi.name);

% Whole liver: volume-weighted mean of all VOIs
livervoi = nvoi + 1;
sizes = cell2mat(roi.size);
roi.name{livervoi} = 'Whole liver';
roi.size{livervoi} = sum(sizes);
roi.tac{livervoi}  = cell2mat(roi.tac) * sizes' / sum(sizes);

%% Data for the model fits
data.c_pl  = aorta.tac{aorta_index};   % arterial input (kBq/ml)
data.t_tac = data.t_tac / 60;          % s -> min
data.kg    = KG_FIXED;
if numel(data.c_pl) ~= numel(data.t_tac)
    error('The arterial TAC (%d frames) and the VOI TACs (%d frames) differ in number of frames.', ...
        numel(data.c_pl), numel(data.t_tac));
end

% Weights: frame duration x decay factor, zero after T_MAX_FIT
data.weights = data.duration .* exp(-log(2)*data.t_tac/T_HALF_O15);
data.weights(data.t_tac > T_MAX_FIT) = 0;

%% Fit all model variants to the whole liver
models = hepatic_flow_models();
nmod   = length(models);
ids    = {models.id};
mi     = @(id) find(strcmp(ids, id));   % model index from id
res    = cell(livervoi, nmod);          % fit results per VOI and model
liver_tac = roi.tac{livervoi};

res{livervoi, mi('kgfix')} = fit_hepatic_flow('kgfix', p0_kgfix, liver_tac, data);
data.delay = res{livervoi, mi('kgfix')}.p(4);     % fixed delay for the 'delay fixed' variants

res{livervoi, mi('kgdelayfix')}      = fit_hepatic_flow('kgdelayfix', p0_kgdelayfix, liver_tac, data);
res{livervoi, mi('full')}            = fit_hepatic_flow('full', p0_full, liver_tac, data);
res{livervoi, mi('coupled')}         = fit_hepatic_flow('coupled', p0_kgfix, liver_tac, data);
res{livervoi, mi('coupled_novb')}    = fit_hepatic_flow('coupled_novb', p0_kgfix(1:4), liver_tac, data);
res{livervoi, mi('kgdelayfix_novb')} = fit_hepatic_flow('kgdelayfix_novb', p0_kgdelayfix(1:3), liver_tac, data);

% Delays for the 'coupled and fixed' variants
data.delaycoupled      = res{livervoi, mi('coupled')}.p(4);
data.delaycoupled_novb = res{livervoi, mi('coupled_novb')}.p(4);

% The 'coupled and fixed' variants are not fitted to the whole liver; they
% are identical to the 'coupled' fits, with kg and delay reported as fixed
% (SD = 0).
r = res{livervoi, mi('coupled')};
r.row([8 10]) = 0;
res{livervoi, mi('coupled_fixed')} = r;
r = res{livervoi, mi('coupled_novb')};
r.row([8 10]) = 0;
res{livervoi, mi('coupled_fixed_novb')} = r;

% Portal-vein input with fixed kg and the whole-liver delay (for display)
data.c_pv = convexvarc(data.t_tac, [data.kg data.kg 0], data.t_tac + data.delay, data.c_pl);

%% Figure: input curves and whole-liver fit
f = figure;
set(f, 'Position', [100 100 600 900], 'PaperOrientation', 'Portrait', ...
    'PaperType', 'A4', 'PaperPosition', [0 0 8.2*2.56 11.6*2.56]);
subplot(5,1,1); title('Segmental hepatic flow'); axis off;
text(0, 0.6, 'Blood file:');
text(0, 0.5, aorta_file, 'Interpreter', 'none');
text(0, 0.4, 'VOI file:');
text(0, 0.3, voi_file, 'Interpreter', 'none');
text(0, 0.1, date);

subplot(5,3,4);
plot(data.t_tac, data.c_pl, 'r-', data.t_tac, data.c_pv, 'b');
legend('Arterial', 'Portal vein'); legend boxoff;
title('Input curves');
ylabel('kBq/ml');
a = axis;
text(5, a(4)*0.5, ['Kg ' num2str(data.kg, '%4.2f') ' /min']);
text(5, a(4)*0.4, ['Delay ' num2str(data.delay*60, '%4.1f') ' s']);

% Values shown in the figure are from the 'kg delay fixed no vb' model
shown = mi('kgdelayfix_novb');
subplot(5,3,5);
plot_voi_fit(data.t_tac, roi.tac{livervoi}, {res{livervoi, shown}.fit}, {'r-'}, ...
    res{livervoi, shown}.row, roi.name{livervoi});

%% Fit all model variants to each VOI
% Initial guesses are the whole-liver estimates.
p_kgfix      = res{livervoi, mi('kgfix')}.p;
p_kgdelayfix = res{livervoi, mi('kgdelayfix')}.p;
p_full       = res{livervoi, mi('full')}.p;
p_coupled    = res{livervoi, mi('coupled')}.p;
p_coupled_nv = res{livervoi, mi('coupled_novb')}.p;

p0 = cell(1, nmod);
p0{mi('kgfix')}              = p_kgfix;
p0{mi('kgdelayfix')}         = p_kgdelayfix;
p0{mi('coupled')}            = p_coupled;
p0{mi('coupled_novb')}       = p_coupled_nv;
p0{mi('coupled_fixed_novb')} = p_coupled_nv(1:3);
p0{mi('kgdelayfix_novb')}    = p_kgdelayfix(1:3);
% NOTE: the two initial guesses below are kept as used for the published
% analysis. For 'full', the whole-liver [kg delay] estimates are placed in
% the [delay VB] positions (kg starts at 1.22); for 'coupled_fixed', the
% whole-liver delay is used as initial guess for VB. Initial guesses only
% affect where the optimisation starts.
p0{mi('full')}          = [p_full(1:3) 1.22 p_full(4:5)];
p0{mi('coupled_fixed')} = p_coupled(1:4);

line_styles = {'r-', 'b-', 'g-', 'c-', 'm-', 'y-', 'k-', 'k--'};   % in model order

for i = 1:nvoi
    for m = 1:nmod
        res{i, m} = fit_hepatic_flow(models(m).id, p0{m}, roi.tac{i}, data);
    end

    subplot(5,3,i+6);
    plot_voi_fit(data.t_tac, roi.tac{i}, cellfun(@(r) r.fit, res(i,:), 'UniformOutput', false), ...
        line_styles, res{i, shown}.row, roi.name{i});
    if i >= 6
        xlabel('Time (min)');
    end
    if any(i == [1 4 7])
        ylabel('kBq/ml');
    end
    drawnow;
end

% Legend for the model curves
subplot(5,3,6);
hold on;
for m = 1:nmod
    plot(NaN, NaN, line_styles{m});
end
hold off;
axis off;
legend({models.label}, 'FontSize', 6, 'Location', 'west'); legend boxoff;

[out_path, out_name] = fileparts(voi_file);
out_base = fullfile(out_path, [out_name '_hepaticflow']);
print(f, [out_base '.pdf'], '-dpdf');

%% Save results
% Columns: F, PF, VT, kg, delay and VB with standard errors (SD = 0 for
% fixed parameters), portal (Fp) and arterial (Fa) flow, AIC and Schwarz
% criterion. Flow in ml/cm3/min, volumes in ml/cm3, kg in 1/min, delay in min.
csv_order = [mi('full') mi('kgdelayfix') mi('kgfix') mi('coupled') ...
    mi('coupled_fixed') mi('coupled_novb') mi('coupled_fixed_novb') mi('kgdelayfix_novb')];
row_fmt = ['%s' repmat(', %4.2f', 1, 19) '\n'];

fid = fopen([out_base '.csv'], 'w');
fprintf(fid, 'Hepatic flow analysis\n');
for m = csv_order
    fprintf(fid, '\n%s\n', models(m).title);
    fprintf(fid, 'ROI name,volume,F,SD,PF,SD,VT,SD,kg,SD,Delay,SD,VB,SD,Fp,SD,Fa,SD,AIC,SCH\n');
    for i = 1:livervoi
        fprintf(fid, row_fmt, roi.name{i}, roi.size{i}, res{i,m}.row, res{i,m}.aic, res{i,m}.sc);
    end
end

% Preferred model: lowest AIC (first model in case of ties)
fprintf(fid, '\nROI name, Preferred model\n');
for i = 1:livervoi
    [~, best] = min(cellfun(@(r) r.aic, res(i,:)));
    fprintf(fid, '%s, %s\n', roi.name{i}, models(best).label);
end
fclose(fid);

fprintf('Results written to\n  %s.pdf\n  %s.csv\n', out_base, out_base);

end


function plot_voi_fit(t, tac, fits, styles, row, name)
% Plot a measured TAC with fitted curves and print F, PF and VT
% (row = [F SD PF SD VT SD ...], see FIT_HEPATIC_FLOW).
args = {t, tac, 'ko'};
for k = 1:length(fits)
    args = [args, {t, fits{k}, styles{k}}]; %#ok<AGROW>
end
plot(args{:});
title(name);
a = axis;
axis([0 6 a(3) a(4)]);
text(2, a(4)*0.3, ['F ' num2str(row(1), '%4.2f') ' ml/g/min']);
text(2, a(4)*0.2, ['PV fraction ' num2str(row(3)*100, '%3.0f') '%']);
text(2, a(4)*0.1, ['V_T ' num2str(row(5), '%4.2f') ' ml/g']);
end
