import sys
import os
sys.path.append('src/')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.style as style
import logging
from sklearn.preprocessing import StandardScaler, MinMaxScaler, RobustScaler, PowerTransformer
from sklearn.metrics import confusion_matrix
from sklearn import metrics
import matplotlib.dates as mdates
import shap  # package used to calculate Shap values
from sklearn.metrics import r2_score, mean_squared_error, make_scorer
import train_test
import random

logger = logging.getLogger('Vital-E')

global model_names
model_names = ['Lasso', 'Random Forest', 'LGBM', 'CatBoost', 'XGBoost']

def bootstrapper(rt):
    resampled_rt = []
    for _ in range(len(rt)):
        x = np.random.choice(rt, size=len(rt), replace=True)
        resampled_rt.append(x.mean())
    return resampled_rt

def take_rt_stats(rt, stat = 0.5):
    rt_df = rt.groupby('collection_date')['Rt'].quantile(stat).to_frame().reset_index()
    rt_df['collection_date'] =  pd.to_datetime(rt_df['collection_date'])
    return rt_df

def make_output_df(dates, y):
    """
    Input: NumPy array of Dates and Output
    Output: A DataFrame of collection Date and output Rt
    """
    dates_df = pd.DataFrame(data=dates, columns=['collection_date'])
    y_df = pd.DataFrame(data=y, columns=['Rt'])
    dates_df.reset_index(inplace = True)
    y_df.reset_index(inplace=True)
    df = pd.concat([dates_df, y_df], axis = 1)
    df['collection_date'] = pd.to_datetime(df['collection_date'])
    return df

def plot_rt_estimate_all_models(y_test, models, X_test, dates, output_dir, filename = 'Rtestimate.png', title = 'Predicted Rt with True Rt', bootstrap_rt = False):
    plt.figure(figsize=(6,6))
    y_test_df = make_output_df(dates, y_test)
    if bootstrap_rt:
        # Remove Duplicates here
        y_test_df.drop_duplicates(subset=['collection_date'], inplace=True)
    plt.plot(y_test_df['collection_date'].dt.to_pydatetime(), y_test_df['Rt'], label='True Rt')
    for i in range(0,len(models)):
        y_pred_raw = models[i].predict(X_test)
        y_pred_df = make_output_df(dates, y_pred_raw)
        if bootstrap_rt:
            title = 'Predicted Rt with True Rt (bootstrapped)'
            y_pred_df = take_rt_stats(y_pred_df)
        plt.plot(y_pred_df['collection_date'].dt.to_pydatetime(), y_pred_df['Rt'],label=model_names[i])
    plt.title(title)
    plt.xlabel('Date')
    plt.ylabel('Rt')
    ax = plt.gca()
    ax.xaxis.set_major_locator(mdates.DayLocator(interval=10))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%m/%d/%Y'))
    ax.set_ylim([-1, 2])
    plt.gcf().autofmt_xdate()
    plt.legend()
    plt.savefig(os.path.join(output_dir,filename), bbox_inches='tight')


def MSE(y_true,y_pred):
    mse = mean_squared_error(y_true, y_pred)
    return mse

def summary_stats_models(models, X_test, y_test, output_dir):
    mse_list = []
    best_hyper_params = []
    for i in range(0, len(models)):
        mse = MSE(y_test, models[i].predict(X_test))
        logger.info("MSE of " + model_names[i] + ": " + str(mse))
        mse_list.append(mse)

        # Also store the Best Hyperparams for each model
        best_hyper_params.append(models[i].best_params_)
        logger.info(str(models[i].best_params_))

    # Write out the MSE in a file format
    summary_df = pd.DataFrame(data = [mse_list], columns = model_names)
    summary_df.to_csv(os.path.join(output_dir, 'summary_stats.csv'), index=False)
    hyperparams_df = pd.DataFrame(data=[best_hyper_params], columns = model_names)
    hyperparams_df.to_csv(os.path.join(output_dir, 'hyperparams.csv'), index=False)

def shap_summary_models(models, X_test, y_test, output_dir):
    feature_names = ['Mean', 'Median', 'Variance', 'Kurtosis', 'Skewness']
    for i in range(0,len(models)):
        plt.figure(figsize = (20,20))
        model = models[i].best_estimator_
        try:
            explainer = shap.TreeExplainer(model)
        except:
            explainer = shap.LinearExplainer(model, X_test)
        shap_values = explainer.shap_values(X_test)
        shap.summary_plot(shap_values, features = X_test, show = False, plot_size = 'auto', feature_names = feature_names)
        plt.title('SHAP summary plot of ' + model_names[i])
        plt.savefig(os.path.join(output_dir, 'SHAP_summary_' +  str(model_names[i]) + '.png'), bbox_inches='tight')
