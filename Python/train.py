SEED = 2022
#!/usr/bin/env python
# Train and test
import sys
import os
sys.path.append('src/')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse
from sklearn.preprocessing import StandardScaler, MinMaxScaler, RobustScaler, PowerTransformer
from sklearn.decomposition import PCA
import logging
import pickle
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import (
    TimeSeriesSplit,
    cross_val_score,
    cross_validate,
    train_test_split,
)
import time
# Import the ML Models
from sklearn import linear_model
from sklearn.ensemble import RandomForestRegressor
from catboost import CatBoostRegressor
import lightgbm as lgb
from xgboost import XGBRegressor

logger = logging.getLogger('Vital-E')

def fit_processor(X_train, output_dir):
    """
    Applies a standard scaler on X_train
    Returns Processor
    """
    logger.info("Preprocessing X_train")
    scaler = StandardScaler()
    norm_X_train = scaler.fit(X_train)
    with open(os.path.join(output_dir, 'processor.pkl'),'wb') as f:
        pickle.dump(scaler, f)
    return scaler

def time_series_cross_val(model, X_train, y_train, param_search):
    """
    Performs Time Series Split Cross Validation on Training Data
    Returns best model
    """

    tscv = TimeSeriesSplit(n_splits=10)
    gsearch = GridSearchCV(estimator=model, cv=tscv,
                            param_grid=param_search, n_jobs = -1)
    best_model = gsearch.fit(X_train, y_train)
    return best_model

def train_all_models(X_train, y_train, output_dir):
    """
    Train All the ML models on X_train and y_train
    Apply cross_validation in a select rolling fashion
    """
    logger.info("Training all models")
    # Initialise models

    #  Lasso
    logger.info("Training Lasso Model")
    lasso_params = {
        'alpha' : [0.001, 0.01, 0.1],
        'max_iter': [500, 1000, 2000]
    }
    lasso = time_series_cross_val(linear_model.Lasso(), X_train, y_train, lasso_params)
    with open(os.path.join(output_dir, 'lasso_model.pkl'), 'wb') as f:
        pickle.dump(lasso, f)

    # Random Forest
    logger.info("Training Random Forest Regressor")
    rf_params = {
        'max_depth': [2, 4, 8, 16, 32],
        'n_estimators': [4, 16, 64, 256],
        'min_samples_split': [2, 4, 8, 16, 32]
    }
    rf = time_series_cross_val(RandomForestRegressor(), X_train, y_train, rf_params)
    with open(os.path.join(output_dir, 'rf_model.pkl'), 'wb') as f:
        pickle.dump(rf, f)

    # LGBM
    logger.info("Training LGBM Regressor")
    lgbm_params = {
        'num_leaves': [4,8,16,32,64,128],
         'max_depth': [2,4,8],
         'min_data_in_leaf': [2,4,8,16,32]
    }
    lgbm = time_series_cross_val(lgb.LGBMRegressor(), X_train, y_train, lgbm_params)
    with open(os.path.join(output_dir, 'lgbm_model.pkl'), 'wb') as f:
        pickle.dump(lgbm, f)

    # Catboost
    logger.info("Training Catboost Regressor")
    catboost_params = {'depth':[1,5,10],
            'iterations':[250, 500, 1000],
            'learning_rate':[0.001, 0.01, 0.1],
            'l2_leaf_reg':[1,5, 10],
            'border_count':[10, 50, 100],
            'thread_count':[4]}
    catboost = time_series_cross_val(CatBoostRegressor(), X_train, y_train, catboost_params)
    with open(os.path.join(output_dir, 'catboost_model.pkl'), 'wb') as f:
        pickle.dump(catboost, f)

    # XGBoost
    logger.info("Training XGBoost Regressor")
    xgboost_params = {
        # Parameters that we are going to tune.
        'max_depth':[6, 10],
        'min_child_weight': [1, 3],
        'eta':[.3, .7],
        'subsample': [1],
        'colsample_bytree': [1],
        'alpha' : [0.01, 0.03],
        'lambda' : [2, 4],
        # Other parameters
        'objective':['reg:squarederror']
    }
    xgboost = time_series_cross_val(XGBRegressor(), X_train, y_train, xgboost_params)
    with open(os.path.join(output_dir, 'xgboost_model.pkl'), 'wb') as f:
        pickle.dump(xgboost, f)

    return lasso, rf, lgbm, catboost, xgboost
