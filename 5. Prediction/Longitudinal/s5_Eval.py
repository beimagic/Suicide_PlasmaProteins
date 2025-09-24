
import numpy as np
import pandas as pd
from tqdm import tqdm
from Utility.Training_Utilities import *
from sklearn.metrics import mean_squared_error, mean_absolute_error
from sklearn.metrics import brier_score_loss, recall_score, roc_auc_score, average_precision_score
from sklearn.metrics import roc_curve
from tqdm import tqdm
pd.options.mode.chained_assignment = None  # default='warn'

def Find_Optimal_Cutoff(target, predicted):
    fpr, tpr, threshold = roc_curve(target, predicted)
    i = np.arange(len(tpr))
    roc = pd.DataFrame({'tf': pd.Series(tpr - (1 - fpr), index=i), 'threshold': pd.Series(threshold, index=i)})
    roc_t = roc.iloc[(roc.tf - 0).abs().argsort()[:1]]
    return list(roc_t['threshold'])

def get_eval(y_test, pred_prob, cutoff):
    pred_binary = threshold(pred_prob, cutoff)
    tn, fp, fn, tp = confusion_matrix(y_test, pred_binary).ravel()
    acc = (tp + tn) / (tp + tn + fp + fn)
    sens = tp / (tp + fn)
    spec = tn / (tn + fp)
    prec = tp / (tp + fp)
    Youden = sens + spec - 1
    f1 = 2 * prec * sens / (prec + sens)
    auc = roc_auc_score(y_test, pred_prob)
    apr = average_precision_score(y_test, pred_prob)
    brier = brier_score_loss(y_test, pred_prob)
    nnd = 1 / Youden
    evaluations = np.round((cutoff, acc, sens, spec, prec, Youden, f1, auc, apr, nnd, brier), 4)
    evaluations = pd.DataFrame(evaluations).T
    evaluations.columns = ['Cutoff', 'Acc', 'Sens', 'Spec', 'Prec', 'Youden', 'F1', 'AUC', 'APR', 'NND', 'BRIER']
    return evaluations

def get_avg_output(mydf, gt_col, pred_col, cutoff, nb_iters):
    idx_lst = [ele for ele in range(len(mydf))]
    out_df = pd.DataFrame()
    for i in range(nb_iters):
        random.seed(i)
        bt_idx = [random.choice(idx_lst) for _ in range(len(idx_lst))]
        mydf_bt = mydf.copy()
        mydf_bt = mydf_bt.iloc[bt_idx, :]
        tmpout_df = get_eval(mydf_bt[gt_col], mydf_bt[pred_col], cutoff)
        out_df = pd.concat([out_df, tmpout_df], axis = 0)
    result_df = out_df.T
    result_df['Median'] = result_df.median(axis=1)
    result_df['LBD'] = result_df.quantile(0.025, axis=1)
    result_df['UBD'] = result_df.quantile(0.975, axis=1)
    output_lst = []
    for i in range(11):
        output_lst.append('{:.3f}'.format(result_df['Median'][i]) + ' [' +
                          '{:.3f}'.format(result_df['LBD'][i]) + ' - ' +
                          '{:.3f}'.format(result_df['UBD'][i]) + ']')
    result_df['output'] = output_lst
    myout = result_df.T
    return myout.iloc[-4:,:]


type = 'Longitudinal'
imp_type = 'TotalGain'

dpath = '/Volumes/JasonWork/Projects/OtherProjects/ZhangBei/SuicideProtein1005/'
outfile = dpath + 'Results/'+type+'_HoldOut/s5_Eval.csv'
mydf = pd.read_csv(dpath + 'Results/'+type+'_HoldOut/s4_Predictions.csv')

ct_pro = Find_Optimal_Cutoff(mydf.target_y, mydf.y_pred_pro)[0]
ct_pro_cov = Find_Optimal_Cutoff(mydf.target_y, mydf.y_pred_pro_cov)[0]
eval_df_pro = get_avg_output(mydf, 'target_y', 'y_pred_pro', ct_pro, nb_iters = 1000)
eval_df_pro_cov = get_avg_output(mydf, 'target_y', 'y_pred_pro_cov', ct_pro_cov, nb_iters = 1000)

eval_df = pd.concat([eval_df_pro, eval_df_pro_cov], axis = 0)
eval_df.to_csv(outfile, index = True)


