import os
import re
from functools import partial
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, mannwhitneyu
from sklearn.metrics import r2_score

from src.consts import muted_columns


def save_plt_fig(name, format):
    if format == 'tiff':
        kwargs = {'compression': 'tiff_lzw'} if format == 'tiff' else None
        plt.savefig(name, format=format, pil_kwargs=kwargs, dpi=350)
    else:
        plt.savefig(name, format=format, dpi=350)


def getattr_with_kwargs(module, method):
    if isinstance(method, dict):
        return partial(getattr(module, method['name']), **method.get('kwargs', {}))

    return getattr(module, method)


def intersect_dfs(dfs):
    common_index = set(dfs[0].index).intersection(*map(lambda df: df.index, dfs[1:]))
    return [df.loc[common_index] for df in dfs]


def make_sure_dir_exists(dir_path):
    if not os.path.isdir(dir_path):
        os.makedirs(dir_path)

    return dir_path


def logit(x):
    return np.log(x / (1 - x))


def inlogit(y):
    return np.exp(y) / (1 + np.exp(y))


def find_substring_occurrences(substring, string):
    return [m.start() for m in re.finditer(substring, string)]


def predict(df, coefs, logit=True):
    cols = set(df.columns) & set(coefs.index)
    df_pred = df[cols].mul(coefs.loc[cols, 'Estimate'], axis=1).sum(axis=1)
    if '(Intercept)' in coefs.index:
        df_pred += coefs.loc['(Intercept)']['Estimate']
    return df_pred if logit else inlogit(df_pred)


def get_accuracy(coefs, df):
    predictions = predict(df, coefs)
    return get_scores(predictions, df['fraq'])


def get_scores(pred, true):
    nans = pred.isna() | true.isna()
    pred, true = pred[~nans], true[~nans]
    if pred.empty:
        cor = 0
        r2 = 0
        mann_w = 0
        mean_pred = 0
        mean_true = 0
        uplift = 1
        mds = 0
    else:
        cor = pearsonr(pred, logit(true))[0]
        r2 = max(0, r2_score(pred, logit(true)))
        mann_w = mannwhitneyu(inlogit(pred), true)[1]
        #
        mean_pred, mean_true = np.median(inlogit(pred)), np.median(true)
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


def aggregated_score(transcript_accuracies, tissue):
    transcripts_variance_train = sum(transcript_accuracies[iso]['tissue'][tissue]['var.train'].dropna().var() for iso in transcript_accuracies)
    transcripts_variance_val = sum(transcript_accuracies[iso]['tissue'][tissue]['var.validation'].dropna().var() for iso in transcript_accuracies)
    return {
        'train': {
            'cor': sum(
                transcript_accuracies[iso]['tissue'][tissue]['train']['cor'] * transcript_accuracies[iso]['tissue'][tissue]['var.train'].dropna().var()
                for iso in transcript_accuracies
            ) / transcripts_variance_train,
            'mds': sum(
                transcript_accuracies[iso]['tissue'][tissue]['train']['mds'] * transcript_accuracies[iso]['tissue'][tissue]['var.train'].dropna().var()
                for iso in transcript_accuracies
            ) / transcripts_variance_train
        },
        'validation': {
            'cor': sum(
                transcript_accuracies[iso]['tissue'][tissue]['validation']['cor'] * transcript_accuracies[iso]['tissue'][tissue][
                    'var.validation'].dropna().var()
                for iso in transcript_accuracies
            ) / transcripts_variance_val,
            'mds': sum(
                transcript_accuracies[iso]['tissue'][tissue]['validation']['mds'] * transcript_accuracies[iso]['tissue'][tissue][
                    'var.validation'].dropna().var()
                for iso in transcript_accuracies
            ) / transcripts_variance_val
        }
    }


def add_freq_to_df(df, min_samples=10, min_expr=3):
    tissue_group = df.groupby(['Tissue'])
    df = tissue_group.filter(lambda g: len(g) >= min_samples)
    tissues = pd.DataFrame(df['Tissue'].value_counts())
    df = df[[c for c in df.columns if c not in muted_columns if (tissue_group[c].quantile(0.8) >= min_expr).any()] + list(set(muted_columns)&set(df.columns))]
    tissues.columns = ['Freq']
    tissues['Freq'] = tissues['Freq'].max() / tissues['Freq']
    tissues['Freq'] = np.log2(tissues['Freq'] + 1)
    tissues['Tissue'] = tissues.index
    df = df.reset_index().merge(tissues, on='Tissue').set_index('index')
    #
    return df
