"""
regressorMetaLearner.py
=================================================================================
Version 2.0: A collection of regression models for meta-learner regression models

-----
Major changes from V1.0:
    - use mlens package instead of manually building a new learner object
    - separate sub-learners by clusters of metabolites (variables)
    - add a new testing/debugging module
"""
# Essential packages
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys
#sys.path.append('../GeneralMethods/')
#from AnalysisKits import *
from SCOOTI.MatplotProp import CanvasStyle, PltProps
from SCOOTI.regressionAnalyzer import *
import networkx as nx
import pingouin as pg
PltProps()
import warnings; warnings.simplefilter('ignore')
from statsmodels.stats.multitest import fdrcorrection
import scipy.stats as ss
import os
from tqdm import tqdm


# regression functions for metabolic modeling results
from sklearn.linear_model import Lasso, LassoLars, ElasticNet, LarsCV, LassoCV, LassoLarsIC, ElasticNetCV, QuantileRegressor, LinearRegression
from sklearn.feature_selection import RFECV
from sklearn.svm import SVR


# packages for super learner models for regression
from math import sqrt
from numpy import hstack
from numpy import vstack
from numpy import asarray
from sklearn.datasets import make_regression
from sklearn.model_selection import KFold
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import ElasticNet


# packages referred to "regressors" that is used to calculate pvalues for sklearn models
#import SCOOBI.regression_stats as regss


from sklearn.metrics import mean_squared_error
from sklearn.exceptions import ConvergenceWarning
from tqdm import tqdm

warnings.simplefilter('ignore', ConvergenceWarning)
sns.set_style("whitegrid")


class Adaline:
    def __init__(self, learning_rate=0.001, epochs=10000):
        self.learning_rate = learning_rate
        self.epochs = epochs

    def fit(self, X, y):
        # inittiate the weights with random array
        self.weights = np.random.random(X.shape[1]+1)

        stop = 0.001
        Error=[stop +1]
        step = 0
        # check the stop condition for the network
        while (Error[-1] > stop or Error[-1]-Error[-2] > 0.0001) and step<self.epochs:
            step += 1
            output = self.net_input(X)
            errors = y - output
            self.weights[1:] += self.learning_rate * X.T.dot(errors)
            self.weights[0] += self.learning_rate * errors.sum()
            # Store sum of square errors
            Error.append((errors**2).sum())    
        #print('Error :',Error[-1], Error[1])

    def net_input(self, X):
        return np.dot(X, self.weights[1:]) + self.weights[0]

    def predict(self, X):
        return self.net_input(X)



class regression_methods:
    
    def __init__(
            self,
            df_variables,
            df_response,
            cluster_label_path='',
            learner='L',
            learning_rate=0.001,
            epo=10000
            ):
        self.df_var = df_variables
        self.df_res = df_response
        if len(cluster_label_path)>0:
            # read cluster labels
            cluster_labels = pd.read_csv(cluster_label_path)
            # remove unmatched rows
            self.df_var = self.df_var.T[
                    self.df_var.columns.isin(cluster_labels.iloc[:, 0])
                    ].T
            # get cluster names and labels
            cluster_label_names = cluster_labels.iloc[:, 0].to_numpy()
            cluster_labels = cluster_labels['cluster'].to_numpy()
            self.cluster_labels = cluster_labels
            self.cluster_label_names = cluster_label_names
        else:
            self.cluster_labels = []
            self.cluster_label_names = []
        self.models = None
        self.meta_model = None
        self.learner = learner
        self.learning_rate = learning_rate
        self.epo = epo
    
    @staticmethod
    def selector_models(X, y):
        """Integrate feature selections into regression models
        
        This function is a static method that supports building regression models
        with feature selection pipelines. Current version only supports recursive
        feature elimination (RFE) which is a function from scikit-learn. RFE was
        employed to mimic the idea of step-wise regression.

        Parameters
        ----------
        X : {array-like},
            data of variables (predictors)
        
        y : {array-like},
            data of outcomes (responses)

        
        Returns
        -------
        coef_dict : dictionary,
            the dictionary contains selected features as keys 
            that comes with corresponding coefficients as values.

        model : class object,
            a regression model that was fitted with selected features.
            the class object is a scikit-learn model instance.
        
        feature_bool : boolean array,
            an boolean array that indicates the feature which is actually columns
            selected by a selector if True; otherwise, 0 means a column (feature)
            is removed.

        Notes
        -----
        The selector fits data with five-fold cross validation that is used to choose
        a best combination of features.
        """
        # define a linear regression model
        estimator = LinearRegression(positive=True, fit_intercept=False)
        # RFE feature selection
        selector = RFECV(estimator, cv=5, step=1)
        # fit a model
        selector = selector.fit(X, y)
        # coefficient
        Xtmp = X.copy()
        # make unselected features 0
        Xtmp[Xtmp.columns[selector.support_]] = 0
        # fit the model again with modified variables
        model = LinearRegression(
                positive=True, # non-negative coefficients
                fit_intercept=False
                ).fit(Xtmp, y)
        
        return model

    
    @staticmethod
    def sklearn_models(X, y, model_select='LinearRegression'):
        """Linear regressors with non-negative coefficients without intercept
        
        This function is a static method that supports building regression models
        with different loss functions. Current version only supports
        linear regression, LASSO, ElasticNet, and LASSO-LARS which are functions
        from scikit-learn.

        Parameters
        ----------
        X : {array-like},
            data of variables (predictors)
        
        y : {array-like},
            data of outcomes (responses)
        
        Returns
        -------
        coef_dict : dictionary,
            the dictionary contains selected features as keys 
            that comes with corresponding coefficients as values.

        model : class object,
            a regression model that was fitted with selected features.
            the class object is a scikit-learn model instance.
        
        Notes
        -----
        The selector fits data with five-fold cross validation that is used to choose
        a best combination of features.
        """
        # define models
        if model_select=='LassoCV':
            # LASSO
            model = LassoCV(
                cv=5, positive=True, fit_intercept=False, random_state=0
            ).fit(X, y)
        elif model_select=='ElasticNetCV':
            # ElasticNet
            model = ElasticNetCV(
                    cv=5, positive=True, fit_intercept=False, random_state=0
                    ).fit(X, y)
        elif model_select=='LassoLarsIC':
            # LASSO-Lars
            model = LassoLarsIC(
                criterion='aic', positive=True, normalize=False, fit_intercept=False,
            ).fit(X, y)
        else:
            # typical linear regression without feature selections
            model = LinearRegression(positive=True, fit_intercept=False).fit(X, y)

        return model

        
    # realize models
    def get_models(self, X, y, return_dict=False):
        if return_dict:
            # make an empty list
            models = {}
            # RFE linear regressor
            models['rfelm'] = self.selector_models(X, y)
            # linear models
            for model_sel, key in zip(['LassoCV', 'ElasticNetCV', 'LassoLarsIC'], ['lasso', 'EN', 'LL']):
                models[key] = self.sklearn_models(X, y, model_select=model_sel)
            return models
        else:
            # make an empty list
            models = list()
            # RFE linear regressor
            models.append(self.selector_models(X, y))
            # linear models
            for model_sel in ['LassoCV', 'ElasticNetCV', 'LassoLarsIC']:
                models.append(self.sklearn_models(X, y, model_select=model_sel))
            return models


    # collect out of fold predictions form k-fold cross validation
    def get_out_of_fold_predictions(self, X, y):
        meta_X, meta_y = list(), list()
        # define split of data
        kfold = KFold(n_splits=5, shuffle=True, random_state=0)
        # enumerate splits
        for train_ix, test_ix in kfold.split(X):
            fold_yhats = list()
            # train-test split
            train_X, test_X = X.iloc[train_ix,:], X.iloc[test_ix,:]
            train_y, test_y = y.iloc[train_ix], y.iloc[test_ix]
            meta_y.extend(test_y)
            # fit and make predictions with each sub-model
            models = self.get_models(train_X, train_y)
            for model in models:
                # make predictions based on training data
                yhat = model.predict(test_X)
                # store columns
                fold_yhats.append(yhat.reshape(len(yhat),1))

            # store fold yhats as columns
            meta_X.append(hstack(fold_yhats))

        return vstack(meta_X), asarray(meta_y)

    ## fit a meta model
    #@staticmethod
    #def fit_meta_model(X, y):
    #    print('fitting meta model...')
    #    model = LinearRegression(positive=True, fit_intercept=False)
    #    model.fit(X, y)
    #    #from sklearn.model_selection import cross_val_score
    #    #scores = cross_val_score(model, X, y, cv=5)
    #    #print(scores)

    #    # 5 fold CV
    #    kfold = KFold(n_splits=5, shuffle=True, random_state=1000)
    #    scores = []
    #    for train_ix, test_ix in kfold.split(X):
    #        train_X, test_X = X[train_ix, :], X[test_ix, :]
    #        train_y, test_y = y[train_ix], y[test_ix]
    #        y_pred = model.predict(test_X)
    #        mse = mean_squared_error(test_y, y_pred)
    #        scores.append(mse)  # Just to keep consistent with the structure
    #    # print(model.coef_)
    #    return model, scores


    # fit a meta model
    @staticmethod
    def fit_meta_model(X, y, learner, learning_rate=0.001, epo=10000):
        if learner=='L':
            print('fitting linear meta model...')
            meta_model = LinearRegression(positive=True, fit_intercept=False)
            meta_model.fit(X, y)

        else: # ADALINE
            print('Fitting Adaline meta model...')
            meta_model = Adaline(learning_rate=learning_rate, epochs=epo)
            meta_model.fit(X, y)

        # 5 fold CV
        kfold = KFold(n_splits=5, shuffle=True, random_state=1000)
        scores = []
        submodels = []
        for train_ix, test_ix in kfold.split(X):
            train_X, test_X = X[train_ix, :], X[test_ix, :]
            train_y, test_y = y[train_ix], y[test_ix]
            y_pred = meta_model.predict(test_X)
            mse = mean_squared_error(test_y, y_pred)
            scores.append(mse)  # Just to keep consistent with the structure
            # calculate robustness
            if learner=='L':
                submodel = LinearRegression(positive=True, fit_intercept=False)
            else:
                submodel = Adaline(learning_rate=learning_rate, epochs=epo)
            # fit fold data
            submodel.fit(train_X, train_y)
            submodels.append(submodel)

        
        #print(f'Mean Squared Error: {mse}')
        return meta_model, scores, submodels


    # fit all base models on the entire dataset
    def fit_base_models(self):
        # fit and make predictions with each sub-model
        res_collect = []
        df_res_tmp = pd.DataFrame(self.df_res)
        #for i in tqdm(range(len(df_res_tmp.columns))):
        for i in range(len(df_res_tmp.columns)):
            # data
            X, y = self.df_var, df_res_tmp.iloc[:, i]
            # fit the meta model
            model_dict = self.get_models(X, y)
            res_collect.append(model_dict)

        # make coefficients into tables
        model_coef_collect = []
        for j in range(len(model_dict)):
            model_coef = {}
            for i in range(len(res_collect)):
                model_coef[i] = res_collect[i][j].coef_
            model_coef_collect.append(pd.DataFrame(model_coef))

        self.models = res_collect

        return model_coef_collect

    # evaluate a list of model on a dataset
    def evaluate_models(self, models):
        for model in models:
            yhat = model.predict(self.df_var)
            mse = mean_squared_error(self.df_res, yhat)
            print('%s: RMSE %.3f' % (model.__class__.__name__, sqrt(mse)))
     
    # make predictions with stacked model
    def super_learner_predictions(self, X, models, meta_model):
        meta_X = list()
        for model in models:
            yhat = model.predict(X)
            meta_X.append(yhat.reshape(len(yhat),1))
        meta_X = hstack(meta_X)
        # predict
        return meta_model.predict(meta_X)

    def clusterSuperLearner(self):
        # save the original data
        ori_var = self.df_var.copy()
        ori_res = self.df_res.copy()
        # list to collect cluster superlearners
        cluster_model_res_list = []
        for cluster in np.unique(self.cluster_labels):
            # make predictions based on training data
            self.df_var = ori_var[ori_var.columns[self.cluster_labels==cluster]]
            # fit superlearner model
            model_res_list = self.SuperLearner()
            # collect coefficients
            cluster_model_res_list.append(model_res_list[-1])
            # replace the variables with the original table
            self.df_var = ori_var.copy()
            self.df_res = ori_res.copy()

        # fit meta-learner models
        meta_coef = {}
        scores_collect = {}
        # fit and make predictions with each sub-model
        for i in tqdm(range(len(self.df_res.columns))):
            # empty lists
            meta_X, meta_y = list(), list()
            # define split of data
            kfold = KFold(n_splits=5, shuffle=True, random_state=0)
            # data
            X, y = self.df_var.copy(), self.df_res.iloc[:, i].copy()
            # enumerate splits
            for train_ix, test_ix in kfold.split(X):
                fold_yhats = list()
                # train-test split
                train_X, test_X = X.iloc[train_ix,:], X.iloc[test_ix,:]
                train_y, test_y = y.iloc[train_ix], y.iloc[test_ix]
                meta_y.extend(test_y)
                for cluster in np.unique(self.cluster_labels):
                    print('cluster:', cluster)
                    print(train_X.shape)
                    # make predictions based on training data
                    self.df_var = train_X[
                            train_X.columns[self.cluster_labels==cluster]
                            ].copy()
                    self.df_res = train_y
                    # fit base models
                    model_res_list = self.fit_base_models()
                    # fit superlearner model
                    model_res_list = self.SuperLearner()
                    # make predictions based on training data
                    cluster_test_X = test_X[
                            test_X.columns[self.cluster_labels==cluster]
                            ].copy()
                    yhat = self.super_learner_predictions(
                            cluster_test_X, self.models[0], self.meta_model
                            )
                    # store columns
                    fold_yhats.append(yhat.reshape(len(yhat),1))

                # store fold yhats as columns
                meta_X.append(hstack(fold_yhats))
            # change format
            meta_X, meta_y = vstack(meta_X), asarray(meta_y)
            # replace the variables with the original table
            self.df_var = ori_var.copy()
            self.df_res = ori_res.copy()
            print('Meta ', meta_X.shape, meta_y.shape)
            # fit the meta model
            meta_model, scores = self.fit_meta_model(
                    meta_X, meta_y, self.learner, learning_rate=self.learning_rate, epo=self.epo
                    )
            if self.learner=='L':
                meta_coef[i] = meta_model.coef_
            else:
                meta_coef[i] = meta_model.weights[1:]  # Adaline does not have a coef_ attribute
            # collect CV scores
            scores_collect[i] = scores
        
        # calculate the final coefficients by the new weight
        integrate_res = cluster_model_res_list[0].copy()
        #for col in cluster_model_res_list[0].columns:
        #print(np.concatenate([model_w[col].mul(meta_coef[col][i]).to_numpy() for i, model_w in enumerate(cluster_model_res_list)], axis=0))
        for i in range(len(cluster_model_res_list)):
            print('debug')
            print(cluster_model_res_list[i].shape)
            
        #return cluster_model_res_list, meta_coef
        # calculate the final coefficients by the new weight
        meta_df = pd.DataFrame(meta_coef)
        integrate_res = []
        for i in range(len(cluster_model_res_list)):
            integrate_res.append(cluster_model_res_list[i].mul(meta_df.iloc[i, :]))
        integrate_res = pd.concat(integrate_res, axis=0)
        cluster_model_res_list.append(integrate_res)
        cluster_model_res_list.append(meta_df)

        return cluster_model_res_list


    def SuperLearner(self):
        # get models
        model_res_list = self.fit_base_models()
        
        # fit meta-learner models
        meta_coef = {}
        scores_collect = {}
        subWeights_collect = {}
        meta_model_collect = []
        df_res_tmp = pd.DataFrame(self.df_res)
        #for i in tqdm(range(len(df_res_tmp.columns))):
        for i in range(len(df_res_tmp.columns)):
            # data
            X, y = self.df_var, df_res_tmp.iloc[:, i]
            # get out of fold predictions
            meta_X, meta_y = self.get_out_of_fold_predictions(X, y)
            #print('Meta ', meta_X.shape, meta_y.shape)
            # fit the meta model
            meta_model, scores, submodels = self.fit_meta_model(
                    meta_X, meta_y, self.learner, learning_rate=self.learning_rate, epo=self.epo
                    )
            if self.learner=='L':
                meta_coef[i] = meta_model.coef_
            else:
                meta_coef[i] = meta_model.weights[1:]  # Adaline does not have a coef_ attribute
            meta_model_collect.append(meta_model)
            scores_collect[i] = scores
            subWeights_collect[i] = submodels

        # calculate correlations of predicted coefficients based on different folds
        mean_corrs = []
        for col in model_res_list[0].columns:
            submodels = subWeights_collect[col]
            subWeights = []
            for submodel in submodels:
                if self.learner=='L':
                    subcoef = submodel.coef_
                else:
                    subcoef = submodel.weights[1:]  # Adaline does not have a coef_ attribute
                subWeight = np.sum([model_w[col].mul(subcoef[i]).to_numpy() for i, model_w in enumerate(model_res_list)], axis=0)
                subWeights.append(pd.DataFrame(subWeight))
            corr = pd.concat(subWeights, axis=1).corr()
            mean_corr = pd.DataFrame(corr.mean(axis=1))
            mean_corr.columns = [col]
            mean_corrs.append(mean_corr)

        # merge correlations
        mean_corrs = pd.concat(mean_corrs, axis=1)

        # only works for one-column response
        self.meta_model = meta_model
        
        # calculate the final coefficients by the new weight
        scores_df = pd.DataFrame(scores_collect)
        print(scores_df.shape)
        print(model_res_list[0].shape)
        scores_df.columns = model_res_list[0].columns
        integrate_res = model_res_list[0].copy()
        for col in model_res_list[0].columns:
            integrate_res[col] = np.sum([model_w[col].mul(meta_coef[col][i]).to_numpy() for i, model_w in enumerate(model_res_list)], axis=0)
        model_res_list.append(integrate_res)
        # replace index with the metabolite name
        for i in range(len(model_res_list)):
           model_res_list[i].index = self.df_var.columns

        # append scores
        model_res_list.append(scores_df)

        # append correlations
        model_res_list.append(mean_corrs)

        return model_res_list


#"""Testing module
#"""
#def testing():
#    # import SCOOTI
#    from SCOOTI.metObjAnalyzer import metObjAnalyzer
#    from SCOOTI.regressionAnalyzer import *
#    #def testing_func():
#    # get unconstrained models
#    root_path = '/nfs/turbo/umms-csriram/daweilin/fluxPrediction/unconstrained_models/pfba/all_mets/DMEMF12/'
#    uncon_res = unconstrained_models(root_path, norm=True, medium='DMEMF12')
#    # get constrained models
#    root_paths = f'/nfs/turbo/umms-csriram/daweilin/fluxPrediction/prolif_qui/'
#    con_res = constrained_models(
#            root_paths+'/',
#            CFR_paraScan=True,
#            norm=False,
#            CFR_k=[0.1,],
#            CFR_r=[10,],# input_path_pattern='NCI60'
#            )   
#    # loading datasets
#    uncon_models, con_models = remove_all_zeros(uncon_res, con_res)
#    con_models = con_models.iloc[:, :3]
#    uncon_models = uncon_models.iloc[:, :]
#    pg_uncon_models, pg_con_models = uncon_models, con_models
#    
#    # cluster labels should be integers
#    #cluster_labels = pd.read_csv('/nfs/turbo/umms-csriram/daweilin/fluxPrediction/unconstrained_models/pfba/cluster_info/recon1_all_mets.csv')
#    #cluster_label_names = cluster_labels.iloc[:, 0].to_numpy()
#    #cluster_labels = cluster_labels['cluster'].to_numpy()
#    reg = regression_methods(
#            pg_uncon_models,
#            pg_con_models,
#            cluster_label_path='/nfs/turbo/umms-csriram/daweilin/fluxPrediction/unconstrained_models/pfba/cluster_info/recon1_all_mets.csv'
#            )
#    #model_res = reg.fit_base_models()
#    #resList = reg.SuperLearner()
#    clusterLearners = reg.clusterSuperLearner()
#
#    return clusterLearners
#
#
#
#
#
#    
