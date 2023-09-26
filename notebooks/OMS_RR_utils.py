import pandas as pd
import numpy as np


def has_fill(df):
    """
    Input: dataframe
    
    Returns: dataframe 
    
    This function returns a dataframe where only entries that have a non-null "fill_number" attribute are kept.
    """
    dfcopy=df[df.fill_number.isna() == False].copy()
    return dfcopy

def add_loc_wrt_fill(df):
    import numpy as np
    dfcopy= df.copy()
    dfcopy["Fill location"] = np.empty(len(dfcopy))
    temp=dfcopy.groupby("fill_number")["run_number"].transform('count')    
    x=1
    for count,i in enumerate(temp):
        if x==i :
            dfcopy['Fill location'].iloc[count]="({}/{})".format(x,i)
            x=1
            continue
        dfcopy['Fill location'].iloc[count]="({}/{})".format(x,i)
        x+=1
        
    return dfcopy



def makeDF(json):
    datadict=json['data'][0]['attributes']
    keys=datadict.keys()
    
    datasetlist=[]
    
    for i in range(len(json['data'])):
        values=json['data'][i]['attributes'].values()
        datasetlist.append(values)
    return pd.DataFrame(datasetlist,columns=keys)#\
#     .set_index(['fill_number','run_number'])\
#     .sort_index()

def convert_check_addFillLoc(json):
    """
    Expects a json from with the attribute "fill_number" in the query
    
    """
    df=makeDF(json)
    #now filter out runs that don't  have fill number
    df=has_fill(df)
    #now add run location wrt fill
    DF_withloc=add_loc_wrt_fill(df)
    return DF_withloc

def get_collisions(rundf,lsdf):
    rundf_coll=rundf[rundf['l1_hlt_mode'].str.contains("collisions")]
    
    lsdf_coll=lsdf.set_index("run_number")\
    .loc[[i for i in rundf_coll.run_number.tolist() if i in lsdf.run_number.unique().tolist()]]\
    .reset_index()
    return rundf_coll,lsdf_coll
    
    
def get_runs_in_ls_df(lsdf):
    """
    Gets the list of run numbers in the lumi df
    """
    return lsdf['run_number'].unique()

def missing_runs(runsdf,lsdf,fromlumi=True):
    '''
    Function takes both run and lumi df and find missing runs.
    Returns:  list of missing runs
    
    ----
    fromlumi :  defaults to True
                This will look for missing run numbers from the Lumi df
                if False : returns missing run numbers from the Run df
    '''
    
    miss_runs=[]
    runsINls=get_runs_in_ls_df(lsdf)
    if fromlumi:
        for i in runsdf['run_number'].values:
            if i not in runsINls:
                miss_runs.append(i)
    else:
        for i in runsINls:
            if i not in runsdf['run_number'].values:
                miss_runs.append(i)
    return miss_runs



def get_pileup(lsdf):
    """
    Gets the pileup from the lumi df and returns a Pandas Series with mean PU per run.
    """
    return lsdf.groupby(['run_number'])['pileup'].mean()

def get_avg_initLumi(lsdf):
    """
    Gets average init luminosity for each run from the LS dataframe
    """
    
    return lsdf.groupby(["run_number"])["init_lumi"].mean()

def relative_diff(df,target,key=None):
    if key==None:
        return 100*(df - target)/(target)
    else:
        return 100*(df[key]- target)/(target)
    
def make_grading(rundf,target,version="V3"):
    if version == "V1":
        # calculate the weights of a run as if it were a grade 
        luminosity_weight = .5*rundf["inst_lumi_delta %"]/100 #(rundf["ave_inst_lumi"]*100)
        pileup_weight = .25*rundf["pileup_delta %"]/100 #rundf["ave_pileup"]
        run_number_weight = .25*rundf["run_number_delta %"]/100# rundf["run_number"]
    elif version == "V2":
        run_instlumi = target[0]
        run_pu = target[1]
        run_num = target[2]
        
        # calculate the weights of a run as if it were a grade 
        luminosity_weight = abs(.5*(rundf["inst_lumi_delta %"]*run_instlumi))/(100*rundf["ave_inst_lumi"])
        pileup_weight = abs(.25*rundf["pileup_delta %"])/run_pu
        run_number_weight = abs(.25*rundf["run_number_delta %"])/run_num

    elif version == "V3":
        run_instlumi = target[0]
        run_pu = target[1]
        run_num = target[2]
        last_lumi = target[3]
#         ZB_rate = target[4]
        
        # calculate the weights of a run as if it were a grade 
        luminosity_weight = abs(rundf["inst_lumi_delta %"])/100
        pileup_weight = abs(rundf["pileup_delta %"])/100
        run_number_weight = abs(rundf["run_number_delta"])/100
        Nlumisections_weight = abs(rundf["num_lumi_delta %"])/100
#         Zerobias_weight= abs(rundf['ZB_rate_delta %'])/100
    
    grade = luminosity_weight + pileup_weight +Nlumisections_weight  + run_number_weight
    
    
    return grade

def drop_short_runs(rundf):
    last_lumi=rundf.last_lumisection_number
    return rundf[last_lumi >= 500]


def ref_rank(DF,runnbs,Trim = False,**kwargs):
    """
    Takes a run number and looks at all the numbers before it in the given pandas dataframe.
    If run number is the first in the frame retruns Null
    
    Usage:
    ref_rank(dictionary,runnbs,Trim=False,**kwargs)
    
    Returns:
    dataframe of runs with ranks per criteria w.r.t given run number
    
    Best rank of the runnb given is 1
    Best rank for the pileup is 0
    Best rank in general is as close to 0 as possible
    """
    if isinstance(DF,(pd.DataFrame,dict)):
        pass
    else:
        raise TypeError("Expecting a pandas dataframe or dictionary but got something else")
        
    if runnbs not in DF['lumisections']['run_number'].values:
        for i in DF['runs']['run_number']:
            if abs(runnbs - i) < 10:
                print("You can try {}".format(i))
        raise KeyError("{} not in the the dataframe\nPlease use a dataframe that includes this run or select any from the suggested runs above".format(runnbs))
        
    ################## Experimental ###################
    if Trim:
            #### RR input ####

        # Idea is to drop runs if not found in GoldenJsons

        golden_UL_2017_runs=list(DF["golden_UL_2017"].keys())
        golden_UL_2018_runs=list(DF["golden_UL_2018"].keys())
        q=[]
        for i in DF['runs'].run_number:
            if str(i) not in golden_UL_2018_runs and str(i) not in golden_UL_2017_runs :
#                 print(i,"is dropped")
                q.append(i)
        runs_skim=DF["runs"].set_index("run_number").drop(q).reset_index()

        # Drop short runs (<500 LS)
        runs_skim = drop_short_runs(runs_skim)
        
        # filter out missing runs
        miss_fromRuns=missing_runs(runs_skim,DF['lumisections'],fromlumi=False)
        miss_fromLumi=missing_runs(runs_skim,DF['lumisections'],fromlumi=True)
        runs_skim= runs_skim.set_index("run_number").drop(miss_fromLumi).reset_index()
        ls_skim = DF["lumisections"].set_index("run_number").drop(miss_fromRuns).reset_index()
#         ZB_missing =missing_runs(runs_skim,DF['Zerobias'],fromlumi=False)
#         ZB_skim = DF["Zerobias"].set_index("run_number").drop(ZB_missing).reset_index()
        
        
        
        # only work with runs before given run number
        runs_tocheck=runs_skim[runs_skim['run_number']< runnbs ]
        ls_tocheck = ls_skim[ls_skim["run_number"]< runnbs ]
        print("Trimming dataframes\n")
    ###################################################
    else: 
        # only work with runs before the given run number
        runs_tocheck= DF['runs'][DF['runs']['run_number']< runnbs ]
        ls_tocheck = DF['lumisections'][DF['lumisections']["run_number"]< runnbs ]

    # get run_number ranks based on distance
    runs_tocheck['run_number_delta']= runs_tocheck.index
    # get a signed percentage of difference
    runs_tocheck["run_number_delta %"]=relative_diff(runs_tocheck,runnbs,"run_number")
    
    # get avg run_number pileups
    
    ls_pileup= get_pileup(ls_tocheck)  
    # ls_skim.groupby(['run_number'])['pileup'].mean()
    runs_tocheck["ave_pileup"] = ls_pileup.values
    runnbspileup= get_pileup(DF['lumisections']).loc[runnbs]  
    # DF['lumisections'].groupby(['run_number'])['pileup'].mean().loc[runnbs]
    runs_tocheck['pileup_delta']= ls_pileup.values - runnbspileup 
    runs_tocheck["pileup_delta %"] = relative_diff(ls_pileup.values,runnbspileup,None)
    
    #### Get number of lumisections rank ####
    last_lumi_runnbs = DF["runs"].set_index("run_number").loc[runnbs].last_lumisection_number
    runs_tocheck["number_of_lumisections_delta"]= runs_tocheck.last_lumisection_number - last_lumi_runnbs 
    runs_tocheck["num_lumi_delta %"] = relative_diff(runs_tocheck,last_lumi_runnbs,"last_lumisection_number")
    
    #### Get luminosity rank #### (average of init lumi per lumisection)
    runAvgInstlumi = get_avg_initLumi(DF["lumisections"]).loc[runnbs]
    runs_tocheck['ave_inst_lumi'] = get_avg_initLumi(ls_tocheck).tolist()
    runs_tocheck["inst_lumi_delta"] =  (get_avg_initLumi(ls_tocheck) - runAvgInstlumi).tolist()
    runs_tocheck["inst_lumi_delta %"] = relative_diff(get_avg_initLumi(ls_tocheck).tolist(),runAvgInstlumi,None )
    
    
    #### Dataset Rate ####
    
    # some code goes here # 
#     ZB = DF["ZeroBias"]
#     runnbZB_rate=ZB.groupby("run_number").mean().rate.loc[runnbs]
#     runs_tocheck['ZB_rate_delta'] = ZB.groupby("run_number").mean().rate.
    
    
    
    
    
    #### Get Ranking #####
    targ= [runAvgInstlumi,runnbspileup,runnbs,last_lumi_runnbs]
    runs_tocheck["Run_Rank"] = make_grading(runs_tocheck,targ,version = "V3")
    
    
    
    
    
#     print(DF["runs"].set_index('run_number').loc[runnbs])
    print("Target run : {}".format(runnbs))
    print("Fill location : {}".format(DF['runs'].set_index('run_number')["Fill location"].loc[runnbs]))
    print("Fill number : {}".format(DF['runs'].set_index('run_number')["fill_number"].loc[runnbs]))
    print("average pileup : {}".format(runnbspileup))
    print("last lumi : {}".format(last_lumi_runnbs))
    print("avg inst lumi : {}".format(runAvgInstlumi))
    print("L1 HLT mode : {}".format(DF['runs'].set_index('run_number')["l1_hlt_mode"].loc[runnbs]))
    
    return runs_tocheck