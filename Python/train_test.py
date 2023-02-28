#!/usr/bin/env python
# Train and test
import sys
import os
sys.path.append('src/')
import train
import test
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse
from sklearn.preprocessing import StandardScaler, MinMaxScaler, RobustScaler, PowerTransformer
from sklearn.decomposition import PCA
import logging
import pickle
from sklearn import preprocessing
import math
import time
# Import the ML Models
from sklearn.ensemble import RandomForestRegressor
from catboost import CatBoostRegressor
import lightgbm as lgb
import xgboost

global model_names
model_names = ['Lasso', 'Random Forest', 'LGBM', 'CatBoost', 'XGBoost']

def impute(df):
    """
    Set to imputing a 0
    """
    df.fillna(0, inplace=True)
    return df

def config_logger(output_dir):
    logger = logging.getLogger("Vital-E")
    logger.setLevel(logging.DEBUG)
    # create handlers
    fh = logging.FileHandler(os.path.join(output_dir, 'train_test_log.txt'))
    fh.setLevel(logging.INFO)
    sh = logging.StreamHandler()
    sh.setLevel(logging.INFO)
    # create formatter
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    sh.setFormatter(formatter)
    logger.addHandler(fh)
    logger.addHandler(sh)
    return logger

def _parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output_dir', default='output', help='Output dir for results')
    parser.add_argument('-d', '--data', default='', help='Path to CSV file of data, see README for format.')
    parser.add_argument('-t', '--train_size', default=0.5, help='Control the training data size, see README for format.')
    parser.add_argument('-b', '--bootstrap', default=False, type=bool, help='Bootstrap predictions, see README for format.')
    return parser.parse_args()

def create_output_dirs(output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

def format_data(df):
    """
    Formats the Rt data to X, y
    """
    y = df['Rt']
    y = y.to_numpy()
    col_names = df.columns.tolist()
    col_names.remove('Rt')
    X = df[col_names]

    return X,y

def main():
    args = _parse_args()
    create_output_dirs(args.output_dir)
    logger = config_logger(args.output_dir)

    # preprocess and split into train and test
    logger.info("Loading and preprocessing data...")

    # Load the data (We ignore the date columns - Assume ORDERED)
    all_data = pd.read_csv(args.data, index_col=False)
    all_data = impute(all_data)


    # train test split (Time series version)
    split = math.floor(float(args.train_size)*len(all_data)) # 50% train
    train_data = all_data[:split]
    test_data  = all_data[split:]
    train_dates = train_data['collection_date']
    test_dates = test_data['collection_date']
    train_dates = pd.to_datetime(train_dates)
    test_dates = pd.to_datetime(test_dates)

    train_data.drop(columns=['collection_date'], inplace=True)
    test_data.drop(columns=['collection_date'], inplace=True)
    # print some stats on data
    logger.info("Total samples: {}".format(len(all_data)))
    logger.info("    Train samples: {}".format(len(train_data)))
    logger.info("    Test samples: {}".format(len(test_data)))

    # Split into X_train and y_train
    X_train, y_train = format_data(train_data)
    X_test, y_test = format_data(test_data)

    # preprocess
    higher_moments = ['Mean', 'Median', 'Variance', 'Kurtosis', 'Skewness']
    processor = train.fit_processor(X_train[higher_moments], args.output_dir)
    norm_X_train = X_train.copy()
    norm_X_train[higher_moments] = processor.transform(X_train[higher_moments])
    norm_X_test = X_test.copy()
    norm_X_test[higher_moments]  = processor.transform(X_test[higher_moments])

    # Train all models
    lasso, rf, lgbm, catboost, xgboost = train.train_all_models(norm_X_train, y_train, args.output_dir)

    # Predict
    all_models = [lasso, rf, lgbm, catboost, xgboost]

    # Print and save plots
    test.plot_rt_estimate_all_models(y_test, all_models, norm_X_test, test_dates, args.output_dir, bootstrap_rt = args.bootstrap)
    test.shap_summary_models(all_models, norm_X_test, y_test, args.output_dir)
    test.summary_stats_models(all_models, norm_X_test, y_test, args.output_dir)

    test.plot_rt_estimate_all_models(y_train, all_models, norm_X_train, train_dates, args.output_dir, filename = 'Rtestimate_train.png', title = 'Predicted Rt with True Rt on Training Set')
if __name__ == "__main__":
    main()
