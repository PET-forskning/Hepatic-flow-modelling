# Segmental hepatic blood flow from dynamic [15O]water PET

MATLAB code used for the kinetic analysis in:

> [Runesson et al, Eur J Nucl Med Mol Imaging 2027, under review]

The code fits a dual-input (hepatic artery + portal vein) one-tissue compartment
model for [15O]water (Kudomi et al., *Eur J Nucl Med Mol Imaging* 2008) to
time-activity curves (TACs) of liver segments. The portal-vein input is derived
from the arterial input by a delay and a mono-exponential dispersion (rate
constant kg). Eight variants of the model, differing in which of kg, delay and
blood volume (VB) are estimated, fixed or omitted, are fitted to every segment,
and the preferred variant is chosen by the Akaike information criterion.

## Requirements

- MATLAB with the Optimization Toolbox (`lsqcurvefit`) and the Statistics and
  Machine Learning Toolbox (`nlparci`).

## Input data

1. **Arterial input function**: a GE VOIager statistics file (`*_stat.txt`)
   containing an aorta VOI.
2. **Liver segment VOIs**: an Excel export from HERMES (`*.xlsx`) with one block
   of rows per VOI (frame start time in column 2, mean activity concentration in
   column 6, VOI volume in column 16).

Both files must contain the same dynamic frames. Activity concentrations are
expected in Bq/ml and frame times in seconds.

## Usage

```matlab
% With file dialogs
analyze_water_hepatic_flow

% Without dialogs (third argument: index of the aorta VOI in the VOIager file)
analyze_water_hepatic_flow('subj01_aorta_stat.txt', 'subj01_liver.xlsx', 1)
```

## Output

Written next to the HERMES file:

- `<name>_hepaticflow.pdf`: input curves and measured/fitted TACs for the whole
  liver and every VOI.
- `<name>_hepaticflow.csv`: for every model variant and VOI, the estimates and
  standard errors of total flow F, portal fraction PF, distribution volume VT,
  kg, delay, VB, portal flow Fp and arterial flow Fa, together with AIC and
  Schwarz criterion; followed by the preferred model per VOI.

Units: flow in ml/cm3/min, volumes in ml/cm3, kg in 1/min, delay in min.

## Files

| File | Purpose |
| --- | --- |
| `analyze_water_hepatic_flow.m` | Main analysis: reads data, fits all models, writes figure and results |
| `hepatic_flow_models.m` | Definitions and parameter bounds of the eight model variants |
| `hepatic_flow_model.m` | Model equations |
| `fit_hepatic_flow.m` | Weighted least-squares fit, standard errors and information criteria |
| `convexvarc.m` | Numerical convolution with an exponential, averaged over PET frames |
| `read_roi_voiager.m` | Reader for VOIager `*_stat.txt` files |
| `read_voi_HERMES.m` | Reader for HERMES Excel exports |

## Analysis details

- Frames are weighted by frame duration × 15O decay factor (half-life 2.03 min);
  frames with mid-time after 6 min are excluded.
- kg is fixed at 1.22 /min in the variants where it is not estimated or coupled
  to the delay (kg = 6 × delay in the "coupled" variants).
- Fixed delays are taken from fits to the whole-liver TAC (volume-weighted mean
  of all VOIs), and whole-liver estimates are used as initial guesses for the
  segment fits.
- Standard errors are derived from the 95% confidence intervals given by
  `nlparci`.

## Licence

[Add licence, e.g. MIT]
