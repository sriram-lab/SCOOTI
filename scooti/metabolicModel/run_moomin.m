function [contFlux, objFlux]=run_moomin(model, df),

% A MILP-solver needs to be set for COBRA TB
% e.g. changeCobraSolver('ibm_cplex','milp');
  
  % process expression data
  geneid = [];
  ppde = [];
  fc = [];
  disp(df)
  disp(size(df.GeneID))
  disp(size(df.PPDE))
  disp(size(df.FC))
  for i=1:length(df.GeneID),
    ppde = [ppde; df.PPDE(i)];
    fc = [fc; df.FC(i)];
  end
  df.GeneID(cellfun(@(x)any(isnan(x)),df.GeneID))={''}
  %model.genes(cellfun(@(x)any(isnan(x)),model.genes))={''}
  model.genes = upper(model.genes);
  df.GeneID = upper(df.GeneID);
  expression.GeneID = char(df.GeneID);
  expression.PPDE = ppde;
  expression.FC = fc;
  % solve for one solution using stoichiometric contraints
  addpath('/home/daweilin/StemCell/moomin/')
  % run moomin
  [model_stoich, MILPsolution_stoich, MILPproblem_stoich] = moomin(model, expression);
  disp('Successfully computed WT fluxes via MOOMIN...')
  contFlux = MILPsolution_stoich{1}.cont;
  objFlux = MILPsolution_stoich{1}.obj;
  %% write the standard output file
  %writeMoominOutput(model_stoich, 'example/output_stoich.txt');

  
end

%% initialize COBRA toolbox
%addpath('/home/daweilin/StemCell/cobratoolbox/')
%run initCobraToolbox%(false)
%changeCobraSolver('gurobi')
%% Load metabolic network model
%model = load('/home/daweilin/StemCell/Project_mESC_JinZhang/SCOOTI/SCOOTI/metabolicModel/GEMs/Shen2019.mat');
%fn = fieldnames(model);
%model = getfield(model, fn{1});
%
%
%% fix rules
%addpath('/home/daweilin/StemCell/Project_mESC_JinZhang/SCOOTI/SCOOTI/metabolicModel/')
%initial_model = model;
%initial_model.b = model.b(:, 1);
%initial_model = fixModelRules(initial_model);
%initial_model.c = initial_model.c*0
%% load expression
%df = load('/nfs/turbo/umms-csriram/daweilin/data/dawei_prolif_qui_datasets/sigMat/Sharma_21_EGF_expression.mat')

%geneid = []
%ppde = []
%fc = []
%for i=1:length(df.GeneID),
%  %geneid = geneid; df.GeneID{i}]
%  ppde = [ppde; df.PPDE(i)];
%  fc = [fc; df.FC(i)];
%end
%expression
%expression.GeneID = char(df.GeneID)
%expression.PPDE = ppde
%expression.FC = fc
%% example format for MOOMIN
%path = '~/StemCell/moomin/example/expression.mat'
%example = load(path);
%example = example.expression
%example
%
%
%function sharma_21_cells_sigGenes_mat(SC_p)
%
%    % single cell cycle path
%    SC_paths = '/nfs/turbo/umms-csriram/daweilin/data/dawei_prolif_qui_datasets/sharma_21_prolif_pmid_34274479.csv'
%    % read data
%    expData = readtable(SC_paths);%, index_col=0)
%    
%    contains(expData.Properties.VariableNames(2:end), prefix_str)
%
%    %expData.columns = pd.Series(expData.columns).apply(lambda x: int(x.split('c')[1]))
%    %expData = expData[sorted(expData.columns)]
%    for i=2:size(expData, 2),
%      ref = sum(expData{:, i})
%      expData{:, i} = expData{:, i}/ref
%    end
%    % get samples
%    cols = expData.Properties.VariableNames(2:end)
%    ctrl1 = {'c1','c5','c9','c15'}
%    ctrl2 = {'c7','c11','c17'}
%    exp1 = {'c2','c6','c10'}
%    exp2 = {'c13','c16'}
%    exp3 = {'c8','c12'}
%    exp4 = {'c14','c18'}
%    % get sig genes
%    df_dict = {
%            'EGF':[exp1, ctrl1],
%            'TPA':[exp2, ctrl1],
%            'H89-EGF':[exp3, ctrl2],
%            'H89-TPA':[exp4, ctrl2]
%            }
%    % calculate pvalues
%    for k in df_dict.keys():
%        genedf = expData[df_dict[k][0]+df_dict[k][1]].T
%        print(genedf.shape)
%        # calculate pvalues
%        pvalues = []
%        foldChanges = []
%        for i in range(genedf.shape[1]):
%            _, p = ss.ttest_ind(
%                    genedf.iloc[:len(df_dict[k][0]), i],
%                    genedf.iloc[len(df_dict[k][0]):, i]
%                    )
%            #print(genedf.iloc[:len(df_dict[k][0]), i].shape)
%            #print(genedf.iloc[len(df_dict[k][0]):, i].shape)
%            p = p if np.isnan(p)==0 else 1
%
%            pvalues.append(p) 
%            fc = genedf.iloc[:len(df_dict[k][0]), i].mean()/genedf.iloc[len(df_dict[k][0]):, i].mean()
%            foldChanges.append(fc)
%        # get up- or down-regulated genes for each cells
%        dwgenedf = (np.array(pvalues)<0.05) & (np.array(foldChanges)<1)
%        upgenedf = (np.array(pvalues)<0.05) & (np.array(foldChanges)>1)
%        
%        print(upgenedf.sum())
%        cell_dir = '/nfs/turbo/umms-csriram/daweilin/data/dawei_prolif_qui_datasets/sigGenes/prolif_qui/'
%        isExist = os.path.exists(cell_dir)
%        if not isExist:
%            # create a new dir if not existing
%            os.makedirs(cell_dir)
%            print(f'Create a folder for {cell_dir}')
%        
%        # save up-/down-regulated genes for each cells
%        pd.DataFrame({
%            'upgenes':genedf.columns[upgenedf].to_numpy()
%            }).to_csv(cell_dir+'Sharma_21'+f'_{k}_upgenes.csv')
%
%        pd.DataFrame({
%            'dwgenes':genedf.columns[dwgenedf].to_numpy()
%            }).to_csv(cell_dir+'Sharma_21'+f'_{k}_dwgenes.csv')
%
%        # save regulators for reverse conditions
%        pd.DataFrame({
%            'dwgenes':genedf.columns[upgenedf].to_numpy()
%            }).to_csv(cell_dir+'Sharma_21_P'+f'_{k}_dwgenes.csv')
%
%        pd.DataFrame({
%            'upgenes':genedf.columns[dwgenedf].to_numpy()
%            }).to_csv(cell_dir+'Sharma_21_P'+f'_{k}_upgenes.csv')
%
%
%def min_18_cells_sigGenes(
%        datasets_repo_path='/nfs/turbo/umms-csriram/daweilin/data/'
%        ):
%
%    # single cell cycle path
%    SC_paths = '/nfs/turbo/umms-csriram/daweilin/data/dawei_prolif_qui_datasets/min_18_GSE122927_ReadCount.csv'
%    # read data
%    expData = pd.read_csv(SC_paths)
%    expData.index = expData.gene
%    # get sample names
%    samples = pd.Series(expData.columns[2:]).apply(
%            lambda x: '_'.join(x.split('_')[:-1])
%            )
%    df_dict = {}
%    # p21_2E2
%    cols = expData.columns[pd.Series(expData.columns).str.contains(r'(?=.*p21_high)(?=.*2E2)|(?=.*p21_low)(?=.*2E2)', regex=True)]
%    cols = sorted(cols)
%    genedf = expData[cols].T
%    df_dict['p21_2E2'] = genedf
%    # p21_3B6
%    cols = expData.columns[pd.Series(expData.columns).str.contains(r'(?=.*p21)(?=.*3B6)', regex=True)]
%    cols = sorted(cols)
%    genedf = expData[cols].T
%    df_dict['p21_3B6'] = genedf
%    # Serum Starvation 2E2
%    cols = expData.columns[pd.Series(expData.columns).str.contains(r'(?=.*SerumStarvation)(?=.*2E2)|(?=.*2E2)(?=.*Control)', regex=True)]
%    genedf = expData[cols].T
%    df_dict['SerumStarvation_2E2'] = genedf
%    # Serum Starvation 3B6
%    cols = expData.columns[pd.Series(expData.columns).str.contains(r'(?=.*SerumStarvation)(?=.*3B6)|(?=.*3B6)(?=.*Control)', regex=True)]
%    genedf = expData[cols].T
%    df_dict['SerumStarvation_3B6'] = genedf
%    # Meki 2E2
%    cols = expData.columns[pd.Series(expData.columns).str.contains(r'(?=.*Meki)(?=.*2E2)|(?=.*2E2)(?=.*Control)', regex=True)]
%    genedf = expData[cols].T
%    df_dict['Meki_2E2'] = genedf
%    # Meki 3B6
%    cols = expData.columns[pd.Series(expData.columns).str.contains(r'(?=.*Meki)(?=.*3B6)|(?=.*3B6)(?=.*Control)', regex=True)]
%    genedf = expData[cols].T
%    df_dict['Meki_3B6'] = genedf
%    # CDK46i 2E2
%    cols = expData.columns[pd.Series(expData.columns).str.contains(r'(?=.*CDK46i)(?=.*2E2)|(?=.*2E2)(?=.*Control)', regex=True)]
%    genedf = expData[cols].T
%    df_dict['CDK46i_2E2'] = genedf
%    # CDK46i 3B6
%    cols = expData.columns[pd.Series(expData.columns).str.contains(r'(?=.*CDK46i)(?=.*3B6)|(?=.*3B6)(?=.*Control)', regex=True)]
%    genedf = expData[cols].T
%    df_dict['CDK46i_3B6'] = genedf
%    # ContactIn 2E2
%    cols = expData.columns[pd.Series(expData.columns).str.contains(r'(?=.*ContactIn)(?=.*2E2)|(?=.*2E2)(?=.*Control)', regex=True)]
%    genedf = expData[cols].T
%    df_dict['ContactIn_2E2'] = genedf
%    # ContactIn 3B6
%    cols = expData.columns[pd.Series(expData.columns).str.contains(r'(?=.*ContactIn)(?=.*3B6)|(?=.*3B6)(?=.*Control)', regex=True)]
%    genedf = expData[cols].T
%    df_dict['ContactIn_3B6'] = genedf
%
%    for k in df_dict.keys():
%        k='ContactIn_2E2'
%        genedf = df_dict[k]
%        # calculate pvalues
%        pvalues = []
%        foldChanges = []
%        for i in range(genedf.shape[1]):
%            _, p = ss.ttest_ind(genedf.iloc[len(genedf)//2:, i], genedf.iloc[:len(genedf)//2, i])
%            p = p if np.isnan(p)==0 else 1
%
%            pvalues.append(p) 
%            fc = genedf.iloc[len(genedf)//2:, i].mean()/genedf.iloc[:len(genedf)//2, i].mean()
%            foldChanges.append(fc)
%        # get up- or down-regulated genes for each cells
%        dwgenedf = (np.array(pvalues)<0.05) & (np.array(foldChanges)<1)
%        upgenedf = (np.array(pvalues)<0.05) & (np.array(foldChanges)>1)
%        
%        print(upgenedf.sum())
%        cell_dir = '/nfs/turbo/umms-csriram/daweilin/data/dawei_prolif_qui_datasets/sigGenes/prolif_qui/'
%        isExist = os.path.exists(cell_dir)
%        if not isExist:
%            # create a new dir if not existing
%            os.makedirs(cell_dir)
%            print(f'Create a folder for {cell_dir}')
%        
%        # save up-/down-regulated genes for each cells
%        pd.DataFrame({
%            'upgenes':genedf.columns[upgenedf].to_numpy()
%            }).to_csv(cell_dir+'min_18_GSE122927'+f'_{k}_upgenes.csv')
%
%        pd.DataFrame({
%            'dwgenes':genedf.columns[dwgenedf].to_numpy()
%            }).to_csv(cell_dir+'min_18_GSE122927'+f'_{k}_dwgenes.csv')
%
%        # save regulators for reverse conditions
%        pd.DataFrame({
%            'dwgenes':genedf.columns[upgenedf].to_numpy()
%            }).to_csv(cell_dir+'min_18_GSE122927_P'+f'_{k}_dwgenes.csv')
%
%        pd.DataFrame({
%            'upgenes':genedf.columns[dwgenedf].to_numpy()
%            }).to_csv(cell_dir+'min_18_GSE122927_P'+f'_{k}_upgenes.csv')
%
%
%def QPcells_sigGenes(
%        datasets_repo_path='/nfs/turbo/umms-csriram/daweilin/data/'
%        ):
%
%    # single cell cycle path
%    SC_paths = '/nfs/turbo/umms-csriram/daweilin/data/dawei_prolif_qui_datasets/johnson_18_GSE117444_prolif_qui_count.csv'
%
%
%    expData = pd.read_csv(SC_paths)
%    expData.index = expData.gene
%    expData = expData.iloc[:, 1:7]
%    
%
%    # assign the expression table to the objective
%    genedf = expData.T
%
%    # calculate pvalues
%    pvalues = []
%    foldChanges = []
%    for i in range(genedf.shape[1]):
%        _, p = ss.ttest_ind(genedf.iloc[:3, i], genedf.iloc[3:, i])
%        pvalues.append(p)
%        fc = genedf.iloc[3:, i].mean()/genedf.iloc[:3, i].mean()
%        foldChanges.append(fc)
%
%    # get up- or down-regulated genes for each cells
%    dwgenedf = (np.array(pvalues)<0.05) & (np.array(foldChanges)<1)
%    upgenedf = (np.array(pvalues)<0.05) & (np.array(foldChanges)>1)
%    cell_dir = '/nfs/turbo/umms-csriram/daweilin/data/dawei_prolif_qui_datasets/sigGenes/prolif_qui/'
%    isExist = os.path.exists(cell_dir)
%    if not isExist:
%        # create a new dir if not existing
%        os.makedirs(cell_dir)
%        print(f'Create a folder for {cell_dir}')
%    
%    # save up-/down-regulated genes for each cells
%    pd.DataFrame({
%        'upgenes':genedf.columns[upgenedf].to_numpy()
%        }).to_csv(cell_dir+'johnson_18_GSE117444'+'_upgenes.csv')
%
%    pd.DataFrame({
%        'dwgenes':genedf.columns[dwgenedf].to_numpy()
%        }).to_csv(cell_dir+'johnson_18_GSE117444'+'_dwgenes.csv')
%
%    # save reverse conditions
%    pd.DataFrame({
%        'dwgenes':genedf.columns[upgenedf].to_numpy()
%        }).to_csv(cell_dir+'johnson_18_GSE117444_P'+'_dwgenes.csv')
%
%    pd.DataFrame({
%        'upgenes':genedf.columns[dwgenedf].to_numpy()
%        }).to_csv(cell_dir+'johnson_18_GSE117444_P'+'_upgenes.csv')
