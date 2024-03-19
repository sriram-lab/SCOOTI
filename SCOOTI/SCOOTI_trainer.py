"""
SCOOTI_trainer.py
========================================================
executable python script to initiate regression training
"""
# packages
import numpy as np
import pandas as pd
import os, sys
from tqdm import tqdm
from SCOOTI.regressorTraining import *
import argparse

# parser
parser = argparse.ArgumentParser()
parser.add_argument("--unconModel", help="Fluxes of unconstrained models")
parser.add_argument("--conModel", help="Fluxes of constrained models")
parser.add_argument("--savePath", help="Root path to save regression models")

# optional
parser.add_argument(
        "--kappaArr",
        default='10,1,0.1,0.01,0.001',
        help="Parameter kappa for CFR; 0 for non-CFR models")
parser.add_argument(
        "--rhoArr",
        default='10,1,0.1,0.01,0.001',
        help="Parameter rho for CFR; 0 for non-CFR models"
        )
parser.add_argument("--expName", default='regression', help="Name of the experiment")
parser.add_argument(
        "--unconNorm", choices=['T', 'F'], default='T', help="Normalization"
        )
parser.add_argument(
        "--conNorm", choices=['T', 'F'], default='F', help="Normalization"
        )
parser.add_argument("--medium", default='DNENF12', help="Medium name")
parser.add_argument("--method", default='cfr', help="Method name; CFR, DFA, or Compass")
parser.add_argument("--model", default='recon1', help="Name of the metabolic model")
parser.add_argument(
        "--inputType", default='flux', help="Type of inputs; flux, pageRank, or stack"
        )
parser.add_argument(
        "--clusterPath", default='', help="path to access the .csv file that clusters metabolites"
        )
parser.add_argument(
        "--objListPath", default='', help="path to access the .csv file that select metabolites as objectives"
        )
parser.add_argument(
        "--rank",
        choices=['T', 'F'],
        default='F',
        help="Rank data if set true"
        )

# get arguments
args = parser.parse_args()
if args.kappaArr:
    kappaArr = [float(s) for s in args.kappaArr.split(',')]

if args.rhoArr:
    rhoArr = [float(s) for s in args.rhoArr.split(',')]
print('debug')
print(args.conNorm)
print(args.unconNorm)
print(args.rank)
args.unconNorm = True if args.unconNorm=='T' else False
args.conNorm = True if args.conNorm=='T' else False
args.rank = True if args.rank=='T' else False
print(args.conNorm)
print(args.unconNorm)
print(args.rank)


# run regression
regression_exe = regressorTraining(
    args.unconModel,
    args.conModel,
    args.savePath,
    kappa_arr=kappaArr,
    rho_arr=rhoArr,
    expName=args.expName,
    uncon_norm=args.unconNorm,
    con_norm=args.conNorm,
    medium=args.medium,
    method=args.method,
    model=args.model,
    input_type=args.inputType,
    cluster_path=args.clusterPath,
    objList_path=args.objListPath,
    rank=args.rank
)


#regression_exe = regressorTraining(
#    unconModel,
#    conModel,
#    savePath,
#    kappa_arr=kappaArr,
#    rho_arr=rhoArr,
#    expName=expName,
#    uncon_norm=True,
#    con_norm=False,
#    medium=medium,
#    method=method,
#    model=model,
#    input_type=inputType
#)
