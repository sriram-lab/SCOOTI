COBRA_path='/home/daweilin/cobratoolbox/'
GEM_path='/home/daweilin/StemCell/Project_mESC_JinZhang/SCOOTI/SCOOTI/metabolicModel/GEMs/Shen2019.mat'
model_name='Recon1' 
obj_candidate_list_file='/home/daweilin/StemCell/Project_mESC_JinZhang/SCOOTI/SCOOTI/metabolicModel/obj52_metabolites_recon1.csv'
input_obj_tb=''

DFA_kappa=-1
CFR_kappa=-1
CFR_rho=-1
paraLen=5 
random_para=0 
init_objective=1 
genekoflag=1
rxnkoflag=0
FVAflag=0
pfba=1 
medium_perturbation=0 
pairwise_CFR_model=0
algorithm='MOOMIN'

data_dir='/nfs/turbo/umms-csriram/daweilin/data/Tzelepis2016/sigMat/'
prefix_name='ESCMOOMIN'
medium='DMEMF12'
late_stage='upgenes'
early_stage='dwgenes'
simulation='MOOMIN'
constraint=1
save_root_path='/nfs/turbo/umms-csriram/daweilin/fluxPrediction/Tzelepis2016/Recon1_MOOMIN/'
CFR_model_path=''
extraWeight=0;

START_NUM=1
END_NUM=paraLen*paraLen

%addpath('/home/daweilin/StemCell/Project_mESC_JinZhang/SCOOTI/SCOOTI/metabolicModel/')
for run=START_NUM:END_NUM,
  multiObj_CBM(run, DFA_kappa, CFR_kappa, CFR_rho, COBRA_path,GEM_path, model_name, obj_candidate_list_file, input_obj_tb, paraLen, random_para, init_objective, genekoflag, rxnkoflag, FVAflag, pfba, medium_perturbation, data_dir, prefix_name, medium, late_stage, early_stage, simulation, constraint, save_root_path, CFR_model_path, pairwise_CFR_model, extraWeight, algorithm)
end



