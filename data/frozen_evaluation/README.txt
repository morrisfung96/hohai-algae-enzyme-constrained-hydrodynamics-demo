Frozen evaluation inputs
========================

Run from the repository root:
  python scripts/reproduce_metrics.py
  python -m unittest discover -s tests -p test_evaluation.py

Only numpy and pandas are needed for this evaluation. No lake simulation,
optimization solver, FACHB905 model, metabolic LP, or new fit is invoked.

Evaluation identities
---------------------
* 164 MC station-date rows, 2019-2022: surface particulate+dissolved total MC,
  ELISA, 0.75-m observation support, nondetects represented by LOD/2. MC error
  is mean absolute log10(prediction+0.01)-log10(observation+0.01).
* 64 biomass rows, 2022: surface-expression comparison on exactly the same
  station-date rows. Observation proxy is measured Chl-a/10, gDW/m3; it does
  not contain a model-derived community fraction. Biomass errors use log10
  without an offset. Original state support/hour counts are retained and may
  be partial on a boundary day; these are not the 164 MC daily supports.
* 181 biomass rows: 141 development rows (2020-2022), 40 unchanged 2019
  conditional transfer-holdout rows. The x1.45 surface scale was frozen from
  development. This is distinct from the 64-row surface-expression comparison.
* 124 rows: environmental/stage stratification on the existing 2020-2022 MC
  subset. This subset must not replace the full 164-row performance set.
* S24: 12 archived LOYO OLS coefficient-interval rows. The release checks their
  point coefficient identity and distributes the archived intervals. It does
  not regenerate their bootstrap fits or a 492-point fitted-center envelope.
* S27: fpr10 two-stage predictions and nested selection records, evaluated
  separately from the formal background OLS and q0.80 predictions.

Source names and reproducibility
--------------------------------
Original archive column names and source values are retained. RevB is only
an internal archive identifier; model_dictionary.json supplies reader-facing
names. The formal prediction table includes historical fixed-mapping and
alternative q0.80 columns so that all 180 archived statistics can be checked.
They are not silently substituted for the formal seven-model comparison.

provenance.json records the source-relative path, source SHA256, released
SHA256, shape and any column-selection operation for each CSV. Absolute local
paths are not required. Hashes are checked before scoring; changed coverage
or values cause an error rather than a changed model-performance claim.
column_dictionary.csv lists fields used in reproduction and reference files.

Current outcomes are preserved, including weak OLS high-tail detection:
dynamic OLS log MAE 0.366237359737, maximum 4.6541237288 ug/L, >=5 TP 0/36;
dynamic q0.80 >=5 TP 19, FP 30; fpr10 >=5 TP 11, FP 20. These branches measure
different prediction tasks and must remain separately labelled.

The paired error intervals resample the frozen predictions in 32 pooled
station-year blocks, 2000 draws, numpy RNG seed 290222; year sensitivity
uses all 256 ordered four-year draws with replacement. No prediction or
mapping is refitted. Source-reference stage similarity and stratified
difference intervals are supplied as archived results; the command rescales
neither their inputs nor their uncertainty.

The upstream FACHB905 network and LP generation are not distributed. This
evaluation reproduces statistics from the supplied computed states and
fixed predictions; it is not an independent reconstruction of that network.
