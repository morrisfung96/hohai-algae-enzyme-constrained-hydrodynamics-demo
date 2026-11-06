"""Metrics used for the archived station-date evaluations.

MC uses log10(concentration + 0.01 ug/L); biomass uses log10(biomass)
with no offset. Correlations use average ranks for ties. Rows are never
silently removed: invalid or missing evaluation inputs raise ValueError.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def _paired(observed, predicted):
    observed = np.asarray(observed, dtype=float)
    predicted = np.asarray(predicted, dtype=float)
    if observed.ndim != 1 or observed.shape != predicted.shape or not len(observed):
        raise ValueError('Expected equal, nonempty one-dimensional evaluation arrays.')
    if not np.isfinite(observed).all() or not np.isfinite(predicted).all():
        raise ValueError('Missing or nonfinite values cannot be silently removed.')
    return observed, predicted


def log_errors(observed, predicted, offset=0.01):
    """Signed prediction minus observation log10 error on the exact rows."""
    observed, predicted = _paired(observed, predicted)
    if offset < 0 or (observed < 0).any() or (predicted < 0).any():
        raise ValueError('Concentrations, biomass and the log offset must be nonnegative.')
    if (observed + offset <= 0).any() or (predicted + offset <= 0).any():
        raise ValueError('Values plus offset must be strictly positive.')
    return np.log10(predicted + offset) - np.log10(observed + offset)


def spearman(observed, predicted):
    """Spearman correlation, using average tied ranks and no scipy dependency."""
    observed, predicted = _paired(observed, predicted)
    a = pd.Series(observed).rank(method='average').to_numpy()
    b = pd.Series(predicted).rank(method='average').to_numpy()
    if np.ptp(a) == 0 or np.ptp(b) == 0:
        return np.nan
    return float(np.corrcoef(a, b)[0, 1])


def rank_auc(positive, score):
    """Mann-Whitney AUROC; ties receive average ranks."""
    positive = np.asarray(positive, dtype=bool)
    score = np.asarray(score, dtype=float)
    if positive.shape != score.shape or positive.ndim != 1 or not np.isfinite(score).all():
        raise ValueError('Invalid classification labels or scores.')
    n_positive = int(positive.sum())
    n_negative = len(positive) - n_positive
    if not n_positive or not n_negative:
        return np.nan
    ranks = pd.Series(score).rank(method='average').to_numpy()
    return float((ranks[positive].sum() - n_positive * (n_positive + 1) / 2)
                 / (n_positive * n_negative))


def concentration_metrics(observed, predicted, offset=0.01):
    observed, predicted = _paired(observed, predicted)
    e = log_errors(observed, predicted, offset)
    return dict(n=len(observed), MAE_log10=float(np.abs(e).mean()), bias_log10=float(e.mean()),
                prediction_max=float(predicted.max()), prediction_p95=float(np.quantile(predicted, .95)),
                prediction_median=float(np.median(predicted)), spearman=spearman(observed, predicted))


def threshold_metrics(observed, predicted, threshold):
    observed, predicted = _paired(observed, predicted)
    truth = observed >= threshold
    call = predicted >= threshold
    tp = int((truth & call).sum()); fp = int((~truth & call).sum())
    fn = int((truth & ~call).sum()); tn = int((~truth & ~call).sum())
    return dict(threshold_ugL=float(threshold), observed_positive=tp+fn,
                predicted_positive=tp+fp, TP=tp, FP=fp, FN=fn, TN=tn,
                recall=tp/(tp+fn) if tp+fn else np.nan,
                PPV=tp/(tp+fp) if tp+fp else np.nan, AUROC=rank_auc(truth, predicted))
