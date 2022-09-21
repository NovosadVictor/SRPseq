import numpy as np
import pandas as pd
from scipy.stats import pearsonr, mannwhitneyu, zscore
from sklearn.metrics import r2_score

from src.consts import muted_columns
from src.utils.common import logit, inlogit


def prepare_model_data(df, is_numpy=False):
    cols = [c for c in df.columns if c not in muted_columns]
    x = df[cols]
    y = logit(df['fraq'])
    #
    return cols, x.to_numpy() if is_numpy else x, y.to_numpy().reshape((-1, 1)) if is_numpy else y


def get_accuracy(model, df, is_numpy=False):
    cols, df_np, _ = prepare_model_data(df, is_numpy=is_numpy)
    if len(cols) == 0:
        predictions = pd.Series([model.intercept_] * len(df_np), index=df.index)
        cor, r2 = 1, 1
    else:
        predictions = model.predict(df_np)
        cor = pearsonr(predictions, logit(df['fraq']))[0]
        r2 = r2_score(predictions, logit(df['fraq']))
    mann_w = mannwhitneyu(inlogit(predictions), df['fraq'])[1]
    #
    mean_pred, mean_true = np.median(inlogit(predictions)), np.median(df['fraq'])
    uplift = (mean_pred - mean_true) / mean_true
    mds = max(0, 1 - abs(uplift))
    #
    return {
        'cor': cor,
        'r2': r2,
        'mann-w': mann_w,
        'uplift': uplift,
        'mds': mds,
        'mean_pred': mean_pred,
        'mean_true': mean_true,
    }


def remove_outliers(df, tr=5):
    df = df[~df['fraq'].isna()]
    return df[np.abs(zscore(df['fraq'])) < tr]
