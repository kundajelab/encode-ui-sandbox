# Some discriminative learning algorithms, building off scikit-learn.
# Author: Akshay Balsubramani

import base64, io, os, time, json, numpy as np, scipy as sp, pandas as pd
from sklearn.metrics import roc_auc_score, f1_score, precision_recall_curve, auc, log_loss
from sklearn.model_selection import KFold
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.cluster.bicluster import SpectralCoclustering




#================================================
#================== Algorithms ==================
#================================================

def discriminate(
    fg_data, bg_data, 
    balance_data=True, num_folds=5, model_type='lr', 
    max_num_pos=5000, 
    mean_center=False, std_scale=True, calc_loss='error', 
    control_features=False
):
    ssc = StandardScaler(with_mean=mean_center, with_std=std_scale)
    fg_data = ssc.fit_transform(fg_data)
    bg_data = ssc.fit_transform(bg_data)
    pos_ndces = np.arange(fg_data.shape[0])
    if (max_num_pos != 0 and max_num_pos < fg_data.shape[0]): # (Else leave pos_ndces as is.)
        pos_ndces = np.random.choice(pos_ndces, size=max_num_pos)
    if balance_data:     # Sample bg_data to be the same size as fg_data
        neg_ndces = np.random.choice(np.arange(bg_data.shape[0]), size=len(pos_ndces))
    datamat = np.row_stack((fg_data, bg_data))
    if control_features:
        datamat = np.column_stack((datamat, permute_columns(datamat)))
    
    labels = np.zeros(fg_data.shape[0] + bg_data.shape[0])
    labels[0:fg_data.shape[0]] = 1
    cv = KFold(n_splits=num_folds, shuffle=True)
    losses = []
    impts = []
    for tr_ndces, te_ndces in cv.split(datamat, labels):
        if model_type == 'tree':
            model_obj = GradientBoostingClassifier(
                min_samples_split=2, max_depth=None, max_features=12, n_estimators=250)
        elif model_type == 'lr':
            model_obj = LogisticRegression(penalty='l1', solver='liblinear')
        X_train = datamat[tr_ndces]
        X_test = datamat[te_ndces]
        y_train = labels[tr_ndces]
        y_test = labels[te_ndces]
        itime = time.time()
        model_obj.fit(X_train, y_train)
        # print("Time to fit {}: {}".format(model_type, str(time.time() - itime)))
        if model_type == 'lr':
            feat_impts = np.squeeze(model_obj.coef_)
        elif model_type == 'tree':       # If sklearn classifier...
            feat_impts = np.array(model_obj.feature_importances_)
        if calc_loss == 'error':
            losses.append(calc_err(y_test, model_obj.predict_proba(X_test)[:,1]))
        elif calc_loss == 'auc':
            losses.append(roc_auc_score(y_test, model_obj.predict_proba(X_test)[:,1]))
        impts.append(feat_impts/np.sum(feat_impts))
    return losses, np.row_stack(impts)



#===============================================
#================== Utilities ==================
#===============================================


# Permute each column of a matrix independently (instead of just reordering the rows)
def permute_columns(x):
    row_ndces = np.random.sample(x.shape).argsort(axis=0)
    col_ndces = np.tile(np.arange(x.shape[1]), (x.shape[0], 1))
    return x[row_ndces, col_ndces]



#===============================================
#=================== Metrics ===================
#===============================================

# Assumes y_true, y_score are in [0,1].
def calc_err(y_true, y_score):
    return 0.5*(1.0 - np.dot((2.0*y_true - 1.0), (2.0*y_score - 1.0))*(1.0/len(y_true)))


def calc_log_loss(y_true, y_score):
    return log_loss(y_true, y_score, labels=[0.0, 1.0])


def calc_recall_at_FDR(y_true, y_score, fdr_cutoff=0.50, weights=None):
    precision, recall, thresholds = precision_recall_curve(y_true, y_score, sample_weight=weights)
    fdr = 1-precision
    cutoff_index = next(i for i, x in enumerate(fdr) if x <= fdr_cutoff)
    return recall[cutoff_index]


def scikitlearn_calc_auPRC(y_true, y_score, weights=None):
    precision, recall, thresholds = precision_recall_curve(y_true, y_score, sample_weight=weights)
    return auc(recall, precision)


def calc_auPRC(y_true, y_score, weights=None, rmode=True):
    if rmode:
        ro.globalenv['pred'] = y_score
        ro.globalenv['labels'] = y_true
        return ro.r('library(PRROC); pr.curve(scores.class0=pred, weights.class0=labels)$auc.davis.goadrich')[0]
    else:  # return sklearn auPRC
        return scikitlearn_calc_auPRC(y_true, y_score, weights=weights)


def calc_auROC(true_labels, predictions, weights=None):
    return roc_auc_score(true_labels, predictions, sample_weight=weights)


def calc_metrics(y_true, y_pred, weights=None):
    au_roc = calc_auROC(y_true, y_pred, weights=weights)
    au_prc = calc_auPRC(y_true, y_pred, weights=weights)
    recall_at_50_fdr = calc_recall_at_FDR(y_true, y_pred, fdr_cutoff=0.50, weights=weights)
    recall_at_25_fdr = calc_recall_at_FDR(y_true, y_pred, fdr_cutoff=0.25, weights=weights)
    recall_at_10_fdr = calc_recall_at_FDR(y_true, y_pred, fdr_cutoff=0.10, weights=weights)
    recall_at_05_fdr = calc_recall_at_FDR(y_true, y_pred, fdr_cutoff=0.05, weights=weights)
    return [au_roc, au_prc, recall_at_50_fdr, recall_at_25_fdr, recall_at_10_fdr, recall_at_05_fdr]


def calc_fpfn_par(y_true, y_score, nonambs_only=False):
    numpos = np.sum(y_true == 1)
    if nonambs_only:
        nonamb_indices = np.where(y_true != -1)[0]
        top_ndces = np.argsort(y_score[nonamb_indices])[::-1]
        ytr_ord = y_true[nonamb_indices][top_ndces]
    else:
        top_ndces = np.argsort(y_score)[::-1]
        ytr_ord = y_true[top_ndces]
    # ypr_ord = np.zeros(len(ytr_ord))
    # ypr_ord[:numpos] = 1
    return {x: np.sum(ytr_ord[:numpos] == x) for x in np.unique(ytr_ord[:numpos])}


def calc_EPR(y_true, y_score, nonambs_only=False):
    m = calc_fpfn_par(y_true, y_score, nonambs_only=nonambs_only)
    p = m[1] if 1 in m else 0
    if (len(m) == 0) or (0 not in m):
        return None
    else:
        return 1.0*p/(m[0]+p)


def thresh_EPR(y_true, y_score, nonambs_only=False, ret_threshold=False):
    numpos = np.sum(y_true == 1)
    toret = np.zeros(len(y_score))
    if nonambs_only:
        nonamb_indices = np.where(y_true != -1)[0]
        natop_ndces = np.argsort(y_score[nonamb_indices])[::-1][0:numpos]
        top_ndces = nonamb_indices[natop_ndces]
    else:
        top_ndces = np.argsort(y_score)[::-1][0:numpos]
    toret[top_ndces] = 1.0
    return toret if not ret_threshold else (toret, y_score[top_ndces[-1]])


def binom_test(tot_fore, hits_fore, tot_bg, hits_bg): # Assumes tot_bg >> tot_fore. TBD: or Fisher's exact test
    return sp.stats.binom_test(hits_fore, n=tot_fore, p=(1.0*hits_bg/tot_bg), )