###############################################################
###################### Gadget3 - Vendace ######################
###############################################################
#Author: Willian T.A.F. Silva.
#E-mail: willian.silva@evobiolab.com

###############################################################
########################### R setup ###########################
###############################################################

#Clear R environment.
rm(list=ls())
#history(inf)
#sessionInfo()

#Call libraries.
library("rstudioapi")
library("mfdb")
library("tidyverse")
library("gadget3")
library("Rgadget")
library("gadgetplots")
library("gadgetutils")

#Set the working directory to folder of the current script file.
#getwd()
scriptfiledir<-dirname(getSourceEditorContext()$path)
setwd(scriptfiledir)
#getwd()

###############################################################
#################### Create model folders #####################
###############################################################

#Create gadget model directory.
modelversion<-format(Sys.time(),"%Y%m%d-%H%M%S")
dirName<-paste0("g3model-vendace-",modelversion)

#Create directories for data, model files and output plots.
modelfolders<-file.path(dirName)
sapply(modelfolders,dir.create)

###############################################################
###################### Import area data #######################
###############################################################
#Load ICES rectangles and subdivisions.
squares<-read.table("../datafiles/ICES_RECT.csv",header=T,sep=";")
squares<-squares[,c("SD","ICES_Rectangle","Area")]
tmp<-aggregate(squares$Area,list("rect"=squares$ICES_Rectangle),sum)
squares<-squares[order(squares[,"ICES_Rectangle"],squares[,"Area"],decreasing=T),]
squares<-squares[!duplicated(squares[,"ICES_Rectangle"]),]
squares$Area<-tmp$x[match(squares$ICES_Rectangle,tmp$rect)]
squares<-squares[squares$SD>=22,]
squares %>% filter(SD %in% c(30,31))
#Add subdivisions as areacell (to avoid assigning dummy rectangles to data at the subdivision scale).
#subdiv<-aggregate(squares$Area,list(sub=squares$SD),sum)
#system("say ICES rectangles loaded!")

###############################################################
##################### Basic model setup #######################
###############################################################

#Define time period. E.g., year_range<-1991:2018.
year_range<-1991:2018 #TBD from data.

#Define species defaults.
ven_defaults<-list(
  area=mfdb_group(#"rect30"=as.character(squares$ICES_Rectangle[squares$SD==30]),
                  "rect31"=as.character(squares$ICES_Rectangle[squares$SD==31])),
  timestep=mfdb_timestep_quarterly,
  year=year_range,
  species="FVE")

#Map area names to integer area numbers (in this case only "rect30"==>1, but could be anything).
ven_areas<-structure(
  seq_along(ven_defaults$area),
  names=names(ven_defaults$area))

#Define time scale.
#An action (i.e. list of formula objects) that will:
#1. Define cur_* variables: cur_time, cur_step, cur_step_size, cur_year, cur_step_final, cur_year_projection, total_steps, total_years.
#If we've reached the end of the model, return nll (negative log likelihood).
ven_time_actions<-list(
  g3a_time( #Add timekeeping to a g3 model
    start_year=min(ven_defaults$year), #Year model run will start.
    end_year=max(ven_defaults$year), #After this year, model run will stop.
    step_lengths=ven_defaults$timestep, #Either an MFDB time grouping, e.g., mfdb::mfdb_timestep_quarterly, or a vector of step lengths which should should sum to 12, e.g., c(3,3,3,3) for quarterly steps within a year.
    final_year_steps=length(ven_defaults$timestep), #Number of steps of final year to include. 
    #Either as an integer or quoted code, in which case it will be calculated when the model runs. 
    #For example: 0 - Model stops before the start of end_year (it is exclusive); 
    #length(step_lengths) - Model stops at the end of end_year (it is inclusive); 
    #2 - Model stops at the second step of end_year, mid-year if step_lengths is quarterly.
    project_years=~g3_param("project_years",default=0,optimise=FALSE), #Number of years to continue running after the "end" of the model. Must be >= 0. Defaults to an unoptimized project_years parameter, set to 0 (i.e. no projection). Generally, you would change this parameter in the parameter template, rather than changing here.
    retro_years=~g3_param("retro_years",default=0,optimise=FALSE), #Adjust end_year to finish model early. Must be >= 0. 
    #Can be used in conjunction with project_years to project instead. The true end year of the model will be end_year-retro_years+project_years.
    #Defaults to an unoptimized retro_years parameter, set to 0. Generally, you would change this parameter in the parameter template, rather than changing here.
    run_at=0) #Integer order that actions will be run within model.
)

#Define age groups.
#ven_agegroups<-c(imm=0:1,mat=2:10) #Not used anywhere. This can be deleted.

###############################################################
########################### Stocks ############################
###############################################################

#Define stock locations.
ven_imm_areas<-ven_areas[c(#"rect30",
                           "rect31")] 
ven_mat_areas<-ven_areas[c(#"rect30",
                           "rect31")] 

#Define stock length groups.
ven_minlength<-3.5 #TBD from data.
ven_maxlength<-20.5 #TBD from data.
ven_dl<-0.5 #TBD from data.
ven_imm_lengroups<-seq(from=ven_minlength,to=17.5,by=ven_dl) #TBD from data.
ven_mat_lengroups<-seq(from=ven_minlength,to=ven_maxlength,by=ven_dl) #TBD from data.
#Maximum number of length groups for each stock within a time step (maxlengthgroupgrowth).
ven_imm_mlgg<-5 #TBD.
ven_mat_mlgg<-5 #TBD.

#Define stock age range.
ven_minage<-0
ven_maxage<-10
ven_imm_minage<-ven_minage #TBD from data.
ven_imm_maxage<-1 #TBD from data.
ven_mat_minage<-ven_minage #TBD from data.
ven_mat_maxage<-ven_maxage #TBD from data.

ven_catchlen<-mfdb_interval("len",seq(ven_minlength,ven_maxlength,by=ven_dl))
ven_catchage<-mfdb_interval("all",c(ven_minage,ven_maxage),open_ended=c("upper","lower"))

#Define gadget3 stocks and stock actions.
ven_imm<-
  g3_stock("ven_imm",ven_imm_lengroups,open_ended=FALSE) %>%
  g3s_livesonareas(ven_imm_areas) %>% #Add area dimensions to g3_stock classes.
  #When iterating over the stock, iterate over each area in turn, area will be set to the current integer area.
  #When intersecting with another stock, only do anything if area is also part of our list of areas.
  #Each area will be defined as a variable in your model as area_x, allowing you to use names in formulas, e.g., run_f=quote(area==area_x).
  g3s_age(ven_imm_minage,ven_imm_maxage) #%>% #Add age dimensions to g3_stock classes.
  #When iterating over the stock, iterate over each age in turn, age will be set to the current integer age.
  #When intersecting with another stock, only do anything if age is between minage and maxage.
  #If an age dimension already exists, it is redefined with new parameters.

ven_mat<-
  g3_stock("ven_mat",ven_mat_lengroups,open_ended=FALSE) %>%
  g3s_livesonareas(ven_mat_areas) %>% #Add area dimensions to g3_stock classes.
  #When iterating over the stock, iterate over each area in turn, area will be set to the current integer area.
  #When intersecting with another stock, only do anything if area is also part of our list of areas.
  #Each area will be defined as a variable in your model as area_x, allowing you to use names in formulas, e.g., run_f=quote(area==area_x).
  g3s_age(ven_mat_minage,ven_mat_maxage) #%>% #Add age dimensions to g3_stock classes.
  #When iterating over the stock, iterate over each age in turn, age will be set to the current integer age.
  #When intersecting with another stock, only do anything if age is between minage and maxage.
  #If an age dimension already exists, it is redefined with new parameters.

#List of stocks.
ven_stocks<-list(ven_imm,ven_mat)

#names(ven_imm$env) #Show which objects can be retrieved from the stock.

###############################################################
###################### Stock life cycle #######################
###############################################################

#Here we simply state which model (equations) we will use to calculate the initial conditions, mortality, ageing, growth, etc. 
#Each model has a set of parameters, whose initial values needs to be given in the parameter template table (created further down).
#The parameters can be specific to a particular stock (by_stock=TRUE) or generalized to all (or a few) stocks (by_stock=list(stock1)).

#Order of calculations.
growmature_order<-1
growmature_transition_order<-2
spawning_order<-3
recruitment_order<-4
renewal_order<-5
natural_mortality_order<-6
ven_fleet_order<-7
ageing_order<-8
ageing_transition_order<-9
migration_order<-10
random_order<-11
likelihood_order<-12

###############################################################
##################### Initial conditions ######################
###############################################################

#Stock initial conditions.
#Add initialconditions to a g3 model.
#An action (i.e., list of formula objects) that will, for the given stock, iterate over each area/age/etc combination, and generate a lengthgroup vector of new individuals and weights using num_f and wgt_f.
#renewal will add fish to the existing collection, whereas initialconditions will assume the stock is currently empty.
ven_initial_conditions<-list(
  #ven_imm.
  g3a_initialconditions_normalparam( #Equations: n=exp(-(((L-mu)/sigma)^2)); N=F*10000*n/sum(n); W=alpha*L^beta.
    ven_imm, #The g3_stock to apply to.
    factor_f= #F. 
      g3a_renewal_initabund( #Formula: scalar∗init∗e−1∗(M+init.F)∗age.
        scalar=g3_parameterized("init.scalar",by_stock=TRUE), #scalar.
        init=g3_parameterized("init",by_stock=TRUE,by_age=TRUE), #init.
        M=g3_parameterized("init.M",by_stock=TRUE), #M.
        init_F=g3_parameterized("init.F",by_stock=ven_stocks), #init.F
        by_stock=TRUE,
        by_stock_f=FALSE),
    mean_f= #mu.
      g3a_renewal_vonb( #Formula: Linf*(1-exp(-1*K*(a-(1+(log(1-r/Linf)/K))))).
        Linf=g3_parameterized("init.Linf",by_stock=ven_stocks), #Linf.
        K=g3_parameterized("init.K",by_stock=ven_stocks,scale=0.001), #K.
        recl=g3_parameterized("init.recl",by_stock=ven_stocks), #r.
        by_stock=ven_stocks),
    stddev_f= #sigma.
      g3_parameterized("init.sd",by_stock=ven_stocks,by_age=TRUE),
    alpha_f= #alpha.
      g3_parameterized("init.alpha",by_stock=TRUE),
    beta_f= #beta.
      g3_parameterized("init.beta",by_stock=TRUE),
    by_age=TRUE, 
    wgt_by_stock=TRUE),
  
  #ven_mat.
  g3a_initialconditions_normalparam( #Equations: n=exp(-(((L-mu)/sigma)^2)); N=F*10000*n/sum(n); W=alpha*L^beta.
    ven_mat, #The g3_stock to apply to.
    factor_f= #F. 
      g3a_renewal_initabund( #Formula: scalar∗init∗e−1∗(M+init.F)∗age.
        scalar=g3_parameterized("init.scalar",by_stock=TRUE), #scalar.
        init=g3_parameterized("init",by_stock=TRUE,by_age=TRUE), #init.
        M=g3_parameterized("init.M",by_stock=TRUE), #M.
        init_F=g3_parameterized("init.F",by_stock=ven_stocks), #init.F
        by_stock=TRUE,
        by_stock_f=FALSE),
    mean_f= #mu.
      g3a_renewal_vonb( #Formula: Linf*(1-exp(-1*K*(a-(1+(log(1-r/Linf)/K))))).
        Linf=g3_parameterized("init.Linf",by_stock=ven_stocks), #Linf.
        K=g3_parameterized("init.K",by_stock=ven_stocks,scale=0.001), #K.
        recl=g3_parameterized("init.recl",by_stock=ven_stocks), #r.
        by_stock=ven_stocks),
    stddev_f= #sigma.
      g3_parameterized("init.sd",by_stock=ven_stocks,by_age=TRUE),
    alpha_f= #alpha.
      g3_parameterized("init.alpha",by_stock=TRUE),
    beta_f= #beta.
      g3_parameterized("init.beta",by_stock=TRUE),
    by_age=TRUE, 
    wgt_by_stock=TRUE))

###############################################################
##################### Natural mortality #######################
###############################################################

#Stock natural mortality.
#Add natural mortality to a g3 model.
#A model can have any number of g3a_naturalmortality actions, so long as the calling arguments are different.
#For instance, run_f=~age==5 and run_f=~age==7.
ven_natural_mortality<-list(
  #ven_imm.
  g3a_naturalmortality( #Action that will, for the given stock, remove a proportion of each stock group as calculated by the mortality formula mortality_f.
    ven_imm, #g3_stock mortality applies to.
    mortality_f= #A mortality formula, as defined by g3a_naturalmortality_exp.
      g3a_naturalmortality_exp( #Formula: exp(-M*deltat).
        param_f= #M.
          g3_parameterized("natmort.M",by_stock=TRUE,by_age=FALSE),
        by_stock=TRUE,
        by_age=FALSE,
        action_step_size_f=~cur_step_size), #action_step_size_f=1 if action only runs yearly.
    run_f=TRUE, #Formula specifying a condition for running this action, default always runs.
    run_at=natural_mortality_order), #Integer order that actions will be run within model.
  
  #ven_mat.
  g3a_naturalmortality( #Action that will, for the given stock, remove a proportion of each stock group as calculated by the mortality formula mortality_f.
    ven_mat, #g3_stock mortality applies to.
    mortality_f= #A mortality formula, as defined by g3a_naturalmortality_exp.
      g3a_naturalmortality_exp( #Formula: exp(-M*deltat).
        param_f= #M.
          g3_parameterized("natmort.M",by_stock=TRUE,by_age=FALSE),
        by_stock=TRUE,
        by_age=FALSE,
        action_step_size_f=~cur_step_size), #action_step_size_f=1 if action only runs yearly.
    run_f=TRUE, #Formula specifying a condition for running this action, default always runs.
    run_at=natural_mortality_order)) #Integer order that actions will be run within model. 

###############################################################
########################### Ageing ############################
###############################################################

#Stock ageing.
#An action (i.e., list of formula objects) that will, for the given stock:
#1. Move the final age group into temporary storage, stock__transitioning_num / stock__transitioning_wgt.
#2. Move the contents of all other age groups into the age group above.
#3. Move the contents of the temporary storage into output_stocks.
ven_ageing<-list(
  #ven_imm.
  g3a_age(ven_imm, #g3_stock to age.
          output_stocks=list(ven_mat), #List of g3_stocks that oldest specimens in stock should move into.
          output_ratios=rep(1/length(list(ven_mat)),times=length(list(ven_mat))), #Vector of proportions for how to distribute into output_stocks, default evenly spread.
          run_f=~cur_step_final, #Formula specifying a condition for running this action, default is end of model year.
          run_at=ageing_order, #Integer order that actions will be run within model.
          transition_at=ageing_transition_order), #Integer order that transition actions will be run within model.
  
  #ven_mat.
  g3a_age(ven_mat, #g3_stock to age.
          output_stocks=list(), #List of g3_stocks that oldest specimens in stock should move into.
          output_ratios=rep(1/length(list()),times=length(list())), #Vector of proportions for how to distribute into output_stocks, default evenly spread.
          run_f=~cur_step_final, #Formula specifying a condition for running this action, default is end of model year.
          run_at=ageing_order, #Integer order that actions will be run within model.
          transition_at=ageing_transition_order)) #Integer order that transition actions will be run within model.

###############################################################
################### Growth and maturation #####################
###############################################################

#Stock maturation.
ven_growmature<-list(
  #ven_imm.
  g3a_growmature( #Add growth/maturity actions to a g3 model.
    #An action (i.e. list of formula objects) that will, for the given stock:
    #1. Move any maturing individuals into temporary storage, stock__transitioning_num / stock__transitioning_wgt.
    #2. Calculate increase in length/weight using growth_f and impl_f.
    #3. Move the contents of the temporary storage into output_stocks
    #A model can have any number of g3a_growmature actions, so long as the calling arguments are different. 
    #For instance, run_f=~age==5 and run_f=~age==7.
    #impl_f's dependent variables are analysed to see what will affect growth. If nothing but cur_step_size will affect growth, then growth will only be recalculated when the step size changes.
    ven_imm, #g3_stock to grow.
    impl_f= #A pair of formula objects, as defined by g3a_grow_impl_bbinom. Both define a matrix of length groups i to length group deltas j (0...maxlengthgroupgrowth), the values in the first indicate the proportion of individuals moving from i to i+j, the values in the second indicate the corresponding weight increase of individuals moving from i to i+j.
      g3a_grow_impl_bbinom( #Formula object converting mean growths using beta-binomial distribution.
        delta_len_f= #A formula defining a non-negative vector for mean increase in length for stock for each lengthgroup, as defined by g3a_grow_lengthvbsimple.
          g3a_grow_lengthvbsimple( #Formula: delta_len_f=(Linf-Li)*(1-exp(-K*deltat)).
            linf_f=g3_parameterized("growth.Linf",by_stock=ven_stocks), #Linf.
            kappa_f=g3_parameterized("growth.K",by_stock=ven_stocks,scale=0.001), #K.
            by_stock=ven_stocks), 
        delta_wgt_f= #A formula defining the corresponding weight increase as a matrix of lengthgroup to lengthgroup delta for stock, as defined by g3a_grow_weightsimple.
          g3a_grow_weightsimple( #Formula: delta_wgt_f=alpha*((Li+delta_lenj)^beta-Li^beta). deltat is the length of the current timestep, delta_len is all possible length group increases, i.e, 0...maxlengthgroupgrowth.
            alpha_f=g3_parameterized("growth.alpha",by_stock=ven_stocks), #alpha.
            beta_f=g3_parameterized("growth.beta",by_stock=ven_stocks), #beta.
            by_stock=ven_stocks), 
        beta_f=g3_parameterized("growth.bbin",by_stock=ven_stocks,scale=10), #beta (binomial).
        maxlengthgroupgrowth=ven_imm_mlgg, #An integer with the maximum length groups an individual can jump in one step.
        by_stock=ven_stocks), 
    maturity_f= #A maturity formula, as defined by g3a_mature_constant.
      g3a_mature_continuous( #Formula: m0*(alpha*delta_len+beta*deltat)^T. m0=1/(1+exp(-alpha*(l-l50)-beta*(a-a50)-gamma*(k-k50))). l=length of stock; l50=length of stock when 50% are mature; a=age of stock; a50=age of stock when 50% are mature; k=weight of stock; k50=weight of stock when 50% are mature.
        alpha=g3_parameterized("mat.alpha",by_stock=ven_stocks,scale=0.001), #alpha.
        l50=g3_parameterized("mat.l50",by_stock=ven_stocks), #A formula to substitute for l50. Must be defined if alpha is defined.
        beta=0, #beta.
        a50=0, #A formula to substitute for a50. Must be defined if beta is defined.
        by_stock=ven_stocks),
    output_stocks=list(ven_mat), #List of g3_stocks that maturing stock should move into.
    output_ratios=rep(1/length(list(ven_mat)),times=length(list(ven_mat))), #Vector of proportions for how to distribute into output_stocks, summing to 1, default evenly spread.
    transition_f=~cur_step_final, #Formula specifying a contition for running maturation steps as well as growth, default final step of year.
    run_f=~TRUE, #Formula specifying a condition for running this action, default always runs.
    run_at=growmature_order, #Integer order that actions will be run within model.
    transition_at=growmature_transition_order),
  
  #ven_mat.
  g3a_growmature( #Add growth/maturity actions to a g3 model.
    #An action (i.e. list of formula objects) that will, for the given stock:
    #1. Move any maturing individuals into temporary storage, stock__transitioning_num / stock__transitioning_wgt.
    #2. Calculate increase in length/weight using growth_f and impl_f.
    #3. Move the contents of the temporary storage into output_stocks
    #A model can have any number of g3a_growmature actions, so long as the calling arguments are different. 
    #For instance, run_f=~age==5 and run_f=~age==7.
    #impl_f's dependent variables are analysed to see what will affect growth. If nothing but cur_step_size will affect growth, then growth will only be recalculated when the step size changes.
    ven_mat, #g3_stock to grow.
    impl_f= #A pair of formula objects, as defined by g3a_grow_impl_bbinom. Both define a matrix of length groups i to length group deltas j (0...maxlengthgroupgrowth), the values in the first indicate the proportion of individuals moving from i to i+j, the values in the second indicate the corresponding weight increase of individuals moving from i to i+j.
      g3a_grow_impl_bbinom( #Formula object converting mean growths using beta-binomial distribution.
        delta_len_f= #A formula defining a non-negative vector for mean increase in length for stock for each lengthgroup, as defined by g3a_grow_lengthvbsimple.
          g3a_grow_lengthvbsimple( #Formula: delta_len_f=(Linf-Li)*(1-exp(-K*deltat)).
            linf_f=g3_parameterized("growth.Linf",by_stock=ven_stocks), #Linf.
            kappa_f=g3_parameterized("growth.K",by_stock=ven_stocks,scale=0.001), #K.
            by_stock=ven_stocks), 
        delta_wgt_f= #A formula defining the corresponding weight increase as a matrix of lengthgroup to lengthgroup delta for stock, as defined by g3a_grow_weightsimple.
          g3a_grow_weightsimple( #Formula: delta_wgt_f=alpha*((Li+delta_lenj)^beta-Li^beta). deltat is the length of the current timestep, delta_len is all possible length group increases, i.e, 0...maxlengthgroupgrowth.
            alpha_f=g3_parameterized("growth.alpha",by_stock=ven_stocks), #alpha.
            beta_f=g3_parameterized("growth.beta",by_stock=ven_stocks), #beta.
            by_stock=ven_stocks), 
        beta_f=g3_parameterized("growth.bbin",by_stock=ven_stocks,scale=10), #beta (binomial).
        maxlengthgroupgrowth=ven_imm_mlgg, #An integer with the maximum length groups an individual can jump in one step.
        by_stock=ven_stocks), 
    maturity_f= #A maturity formula, as defined by g3a_mature_constant.
      g3a_mature_continuous( #Formula: m0*(alpha*delta_len+beta*deltat)^T. m0=1/(1+exp(-alpha*(l-l50)-beta*(a-a50)-gamma*(k-k50))). l=length of stock; l50=length of stock when 50% are mature; a=age of stock; a50=age of stock when 50% are mature; k=weight of stock; k50=weight of stock when 50% are mature.
        alpha=g3_parameterized("mat.alpha",by_stock=ven_stocks,scale=0.001), #alpha.
        l50=g3_parameterized("mat.l50",by_stock=ven_stocks), #A formula to substitute for l50. Must be defined if alpha is defined.
        beta=0, #beta.
        a50=0, #A formula to substitute for a50. Must be defined if beta is defined.
        by_stock=ven_stocks),
    output_stocks=list(ven_mat), #List of g3_stocks that maturing stock should move into.
    output_ratios=rep(1/length(list(ven_mat)),times=length(list(ven_mat))), #Vector of proportions for how to distribute into output_stocks, summing to 1, default evenly spread.
    transition_f=~cur_step_final, #Formula specifying a contition for running maturation steps as well as growth, default final step of year.
    run_f=~TRUE, #Formula specifying a condition for running this action, default always runs.
    run_at=growmature_order, #Integer order that actions will be run within model.
    transition_at=growmature_transition_order)) #Integer order that transition actions will be run within model.

###############################################################
################# Spawning and recruitment ####################
###############################################################

#Reproduction (spawning).
ven_reproduction<-list(
  g3a_spawn( #Add spawning to a g3 model.
    #An action (i.e. list of formula objects) that will, for the given stock:
    #1. Use proportion_f to calculate the total parent stock that will spawn.
    #2. Use recruitment_f to derive the total newly spawned stock.
    #3. Apply weightloss_f and mortality_f to the parent stock.
    #Then, at recruitment stage:
    #1. Recruit evenly into output_stocks, using mean_f, stddev_f, alpha_f, beta_f as-per g3a_renewal_normalparam.
    ven_mat, #The mature g3_stock that will spawn in this action.
    recruitment_f= #A list of formula generated by one of the g3a_spawn_recruitment_* functions, containing S (formula run for each subset of stock) and R (final formula for calculating number of recruits for spawning action).
      g3a_spawn_recruitment_simplessb(mu=g3_parameterized("spawn.rec.mu",by_stock=ven_stocks)), #Equations: S=Nal*p*Wal; R=mu*S. Nal=Number of parent stock; p=Proportion of parent stock spawning, from proportion_f; Wal=Weight of parent stock.
    #g3a_spawn_recruitment_fecundity(p0,p1,p2,p3,p4), 
    #g3a_spawn_recruitment_ricker(mu,lambda),
    #g3a_spawn_recruitment_bevertonholt(mu,lambda),
    #g3a_spawn_recruitment_hockeystick(r0,blim),
    proportion_f= #Formula generated by one of the g3_suitability_* functions, describing the proportion of stock that will spawn at this timestep.
      g3_suitability_exponentiall50( #Formula: 1/(1+exp(-alpha*(l-l50))). l=Vector of stock midlength for each lengthgroup; l50=Length of the stock with a 50% probability of predation.
        alpha=g3_parameterized("spawn.prop.alpha",by_stock=ven_stocks,scale=-1),
        l50=g3_parameterized("spawn.prop.l50",by_stock=ven_stocks)), 
    #g3_suitability_andersen(p0,p1,p2,p3=p4,p4,p5=~pred_stock__midlen),
    #g3_suitability_gamma(alpha,beta,gamma),
    #g3_suitability_exponential(alpha,beta,gamma,delta),
    #g3_suitability_straightline(alpha,beta),
    #g3_suitability_constant(alpha),
    #g3_suitability_richards(p0,p1,p2,p3,p4),
    mortality_f= #Formula generated by one of the g3_suitability_* functions, describing the proportion of spawning stock that will die during spawning.
      g3_suitability_straightline( #Formula: alpha+beta*l. l=Vector of stock midlength for each lengthgroup.
        alpha=g3_parameterized("spawn.mort.alpha"),
        beta=g3_parameterized("spawn.mort.beta")), 
    weightloss_f= #Formula generated by one of the g3_suitability_* functions, describing the overall weight loss during spawning.
      g3_suitability_constant( #Formula: alpha.
        alpha=g3_parameterized("spawn.wgtloss.alpha")), 
    output_stocks=list(ven_imm), #List of g3_stocks that will be spawned into.
    output_ratios=rep(1/length(list(ven_imm)),times=length(list(ven_imm))), #Vector of proportions for how to distribute into output_stocks, summing to 1, default evenly spread.
    #Recruit evenly into output_stocks, using mean_f, stddev_f, alpha_f, beta_f as-per g3a_renewal_normalparam.
    mean_f= #mu. Formula substituted into normalparam calculations.
      g3a_renewal_vonb( #Formula: Linf*(1-exp(-1*K*(a-(1+(log(1-r/Linf)/K))))).
        Linf=g3_parameterized("renewal.Linf",by_stock=ven_stocks), #Linf.
        K=g3_parameterized("renewal.K",by_stock=ven_stocks,scale=0.001), #K.
        recl=g3_parameterized("renewal.recl",by_stock=ven_stocks), #r.
        by_stock=ven_stocks), 
    stddev_f= #sigma. Formula substituted into normalparam calculations.
      g3_parameterized("renewal.rec.sd",by_stock=ven_stocks,by_age=FALSE), 
    alpha_f= #alpha. Formula substituted into normalparam calculations.
      g3_parameterized("renewal.alpha",by_stock=ven_imm), 
    beta_f= #beta. Formula substituted into normalparam calculations.
      g3_parameterized("renewal.beta",by_stock=ven_imm), 
    run_f=~TRUE, #Formula specifying a condition for running this action, default always runs (~TRUE).
    run_at=spawning_order, #Integer order that spawning actions will be run within model.
    recruit_at=recruitment_order)) #Integer order that recruitment from spawning will be run within model.


###############################################################
########################## Renewal ############################
###############################################################

#Stock renewal.
#Add renewal to a g3 model.
#A model can have any number of g3a_renewal_* actions, so long as the calling arguments are different.
#For instance, run_f=~age==5 and run_f=~age==7.
#The g3a_renewal_* actions will define the following stock instance variables for stock:
#stock__renewalnum: Extra individuals added to the stock.
#stock__renewalwgt: Mean weight of added individuals.
#An action (i.e., list of formula objects) that will, for the given stock, iterate over each area/age/etc combination, and generate a lengthgroup vector of new individuals and weights using num_f and wgt_f.
#renewal will add fish to the existing collection, whereas initialconditions will assume the stock is currently empty.
ven_renewal<-list(
  g3a_renewal_normalparam( #Equations: n=exp(-(((L-mu)/sigma)^2)); N=F*10000*n/sum(n); W=alpha*L^beta.
    ven_imm, #The g3_stock to apply to.
    factor_f= #F. Formula substituted into normalparam calculations.
      g3_parameterized("renewal.rec",
                       by_stock=ven_stocks, 
                       by_year=TRUE,
                       scale=g3_parameterized(name="renewal.rec.scalar",
                                              by_stock=ven_stocks, 
                                              exponentiate=FALSE),
                       ifmissing=NaN), 
    mean_f= #mu. Formula substituted into normalparam calculations.
      g3a_renewal_vonb( #Formula: Linf*(1-exp(-1*K*(a-(1+(log(1-r/Linf)/K))))).
        Linf=g3_parameterized("renewal.Linf",by_stock=ven_stocks), #Linf.
        K=g3_parameterized("renewal.K",by_stock=ven_stocks,scale=0.001), #K.
        recl=g3_parameterized("renewal.recl",by_stock=ven_stocks), #r.
        by_stock=TRUE), 
    stddev_f= #sigma. Formula substituted into normalparam calculations.
      g3_parameterized("renewal.rec.sd",by_stock=ven_stocks,by_age=FALSE), 
    alpha_f= #alpha. Formula substituted into normalparam calculations.
      g3_parameterized("renewal.alpha",by_stock=ven_stocks), 
    beta_f= #beta. Formula substituted into normalparam calculations.
      g3_parameterized("renewal.beta",by_stock=ven_stocks), 
    by_stock=ven_stocks, 
    by_age=FALSE, 
    wgt_by_stock=ven_stocks, 
    run_f=~TRUE, #Formula specifying a condition for running this action, default always runs for renewal, first timestep for initialconditions.
    run_at=renewal_order)) #Integer order that actions will be run within model.

###############################################################
######################### Migration ###########################
###############################################################

#Migration.
#Add migration to a g3 model.
#To restrict movement to a particular step in a year, or a particular area, use run_f. For example:
#cur_step == 1: Migration will happen on first step of every year.
#cur_step == 1 && cur_year >= 1990: Migration will happen on first step of every year after 1990.
#cur_step == 2 && area = 1: Migration will happen on second step of every year, in the first area.
#Multiple migration actions can be added, for a separate spring and autumn migration, for instance.
#The action will define the following stock instance variables for each given stock: 
#stock__migratematrix: a × a array, containing proportion of (stock) moved from one area to another. If NaN, no movement has occurred.
ven_migration<-list(
  #ven_imm.
  g3a_migrate(
    #An action (i.e. list of formula objects) that will, for the given stock:
    #1. Fill in stock__migratematrix using migrate_f and normalize_f.
    #2. Apply movement to stock.
    ven_imm, #The g3_stock that will migrate in this action.
    migrate_f=~if(area==1 && dest_area==2){0}else{0}, #A formula describing the migration in terms of (source) area and dest_area.
    normalize_f= #Function to normalize a vector of possible destinations, to make sure fish aren't added or destroyed.
      g3a_migrate_normalize(
        #A formula transforming stock__migratematrix[,stock__area_idx] (i.e. all possible destinations from a given area) by:
        #1. Squaring so values are all positive.
        #2. Altering the proportion of static individuals so a row sums to row_total.
        #3. Dividing by row_total so a row sums to 1.
        row_total=1), #When calculating the proportion of individuals that will stay in place, use this total for what rows are expected to sum to.
    run_f=~TRUE, #Formula specifying a condition for running this action, default always runs.
    run_at=migration_order), #Integer order that spawning actions will be run within model.
  
  #ven_mat.
  g3a_migrate(
    ven_mat, #The g3_stock that will migrate in this action.
    migrate_f=~if(area==1 && dest_area==2){0}else{0}, #A formula describing the migration in terms of (source) area and dest_area.
    normalize_f= #Function to normalize a vector of possible destinations, to make sure fish aren't added or destroyed.
      g3a_migrate_normalize( 
        #A formula transforming stock__migratematrix[,stock__area_idx] (i.e. all possible destinations from a given area) by:
        #1. Squaring so values are all positive.
        #2. Altering the proportion of static individuals so a row sums to row_total.
        #3. Dividing by row_total so a row sums to 1.
        row_total=1), #When calculating the proportion of individuals that will stay in place, use this total for what rows are expected to sum to.
    run_f=~TRUE, #Formula specifying a condition for running this action, default always runs. 
    run_at=migration_order)) #Integer order that spawning actions will be run within model.

###############################################################
################### Import empirical data #####################
###############################################################

source("gadget3-vendace-data-v4.R")

###############################################################
########################### Fleets ############################
###############################################################

#Define fleets.
ven_fleet1<-
  g3_fleet("ven_fleet1") %>%
  g3s_livesonareas(ven_areas)

ven_fleet2<-
  g3_fleet("ven_fleet2") %>%
  g3s_livesonareas(ven_areas)

ven_fleet3<-
  g3_fleet("ven_fleet3") %>%
  g3s_livesonareas(ven_areas)

#Define fleet actions.
#Add predation to a g3 model.
#g3a_predate_fleet:
#An action (i.e. list of formula objects) that will:
#1. Zero fleet and prey catch counters.
#2. For each prey, collect maximum desired by fleet for all prey, into a prey_stock__predby_fleet_stock variable.
#3. After all fleet consumption is done, scale consumption using catchability_f, sum into prey_stock__totalpredate.
#4. After all consumption is done, temporarily convert prey_stock__predby_fleet_stock to a proprotion of prey_stock__totalpredate.
#5. Calculate prey_stock__consratio (ratio of consumed to available), capping using overconsumption_f. Update prey_stock__num.
#6. Recalculate prey_stock__predby_fleet_stock, fleet_stock__catch, post-overconsumption.
ven_fleet_actions<-list(
  #ven_fleet1: com.catch.ven.bench2021 (comven1).
  g3a_predate_fleet( 
    fleet_stock=ven_fleet1, #g3_stock that describes the harvesting fleet.
    prey_stocks=ven_stocks, #List of g3_stocks that maturing stock should move into.
    suitabilities=list(#Either a list of stock names to formula objects, with an optional unnamed default option, or a formula object (which is always used). Each formula should define suitability of a stock, for example by using g3_suitability_exponentiall50.
      g3_suitability_exponentiall50( #Formula: 1/(1+exp(-alpha*(l-l50))). l=Vector of stock midlength for each lengthgroup; l50=Length of the stock with a 50% probability of predation.
        alpha=g3_parameterized("fleet1.alpha",by_stock=ven_stocks), #alpha.
        l50=g3_parameterized("fleet1.l50",by_stock=ven_stocks)) #l50.
        #g3_suitability_andersen(p0,p1,p2,p3=p4,p4,p5=~pred_stock__midlen),
        #g3_suitability_gamma(alpha,beta,gamma),
        #g3_suitability_exponential(alpha,beta,gamma,delta),
        #g3_suitability_straightline(alpha,beta),
        #g3_suitability_constant(alpha),
        #g3_suitability_richards(p0,p1,p2,p3,p4),
    ),
    catchability_f= #A formula generated by e.g. g3a_predate_catchability_totalfleet expressing the catchability (i.e. total biomass) for that fleet.
      g3a_predate_catchability_totalfleet( #Formula: ?????. E=Biomass caught by fleet. A formula that returns the total biomass a stock can harvest in the current time/area, generally defined by a g3_timeareadata table.
        E=g3_timeareadata("ven_fleet1_landings",ven_fleet1_landings,value_field="total_weight")), 
      #g3a_predate_catchability_numberfleet(E), #E=Numbers caught by fleet.
      #g3a_predate_catchability_linearfleet(E),
      #g3a_predate_catchability_effortfleet(catchability_fs,E),
      #g3a_predate_catchability_quotafleet(quota_table,E,sum_stocks=list(),recalc_f=NULL),
    overconsumption_f= #Overconsumption rule, a formula that should cap all values in stock__consratio to <= 95.
      quote(logspace_add_vec(stock__consratio*-1e3,0.95*-1e3)/-1e3),
    run_f=~TRUE, #Formula specifying a condition for running this action, default always runs.
    run_at=ven_fleet_order), #Integer order that actions will be run within model.
  
  #ven_fleet2: com.catch.ven.bench2021 (comven2).
  g3a_predate_fleet( 
    fleet_stock=ven_fleet2, #g3_stock that describes the harvesting fleet.
    prey_stocks=ven_stocks, #List of g3_stocks that maturing stock should move into.
    suitabilities=list(#Either a list of stock names to formula objects, with an optional unnamed default option, or a formula object (which is always used). Each formula should define suitability of a stock, for example by using g3_suitability_exponentiall50.
      g3_suitability_exponentiall50( #Formula: 1/(1+exp(-alpha*(l-l50))). l=Vector of stock midlength for each lengthgroup; l50=Length of the stock with a 50% probability of predation.
        alpha=g3_parameterized("fleet2.alpha",by_stock=ven_stocks), #alpha.
        l50=g3_parameterized("fleet2.l50",by_stock=ven_stocks)) #l50.
        #g3_suitability_andersen(p0,p1,p2,p3=p4,p4,p5=~pred_stock__midlen),
        #g3_suitability_gamma(alpha,beta,gamma),
        #g3_suitability_exponential(alpha,beta,gamma,delta),
        #g3_suitability_straightline(alpha,beta),
        #g3_suitability_constant(alpha),
        #g3_suitability_richards(p0,p1,p2,p3,p4),
    ),
    catchability_f= #A formula generated by e.g. g3a_predate_catchability_totalfleet expressing the catchability (i.e. total biomass) for that fleet.
      g3a_predate_catchability_totalfleet( #Formula: ?????. E=Biomass caught by fleet. A formula that returns the total biomass a stock can harvest in the current time/area, generally defined by a g3_timeareadata table.
        E=g3_timeareadata("ven_fleet2_landings",ven_fleet2_landings,value_field="total_weight")), 
      #g3a_predate_catchability_numberfleet(E), #E=Numbers caught by fleet.
      #g3a_predate_catchability_linearfleet(E),
      #g3a_predate_catchability_effortfleet(catchability_fs,E),
      #g3a_predate_catchability_quotafleet(quota_table,E,sum_stocks=list(),recalc_f=NULL),
    overconsumption_f= #Overconsumption rule, a formula that should cap all values in stock__consratio to <= 95.
      quote(logspace_add_vec(stock__consratio*-1e3,0.95*-1e3)/-1e3),
    run_f=~TRUE, #Formula specifying a condition for running this action, default always runs.
    run_at=ven_fleet_order), #Integer order that actions will be run within model.
  
  #ven_fleet3: cpue_comm_all_bench2021.
  g3a_predate_fleet( 
    fleet_stock=ven_fleet3, #g3_stock that describes the harvesting fleet.
    prey_stocks=ven_stocks, #List of g3_stocks that maturing stock should move into.
    suitabilities=list(#Either a list of stock names to formula objects, with an optional unnamed default option, or a formula object (which is always used). Each formula should define suitability of a stock, for example by using g3_suitability_exponentiall50.
      g3_suitability_exponentiall50( #Formula: 1/(1+exp(-alpha*(l-l50))). l=Vector of stock midlength for each lengthgroup; l50=Length of the stock with a 50% probability of predation.
        alpha=g3_parameterized("fleet3.alpha",by_stock=ven_stocks), #alpha.
        l50=g3_parameterized("fleet3.l50",by_stock=ven_stocks)) #l50.
        #g3_suitability_andersen(p0,p1,p2,p3=p4,p4,p5=~pred_stock__midlen),
        #g3_suitability_gamma(alpha,beta,gamma),
        #g3_suitability_exponential(alpha,beta,gamma,delta),
        #g3_suitability_straightline(alpha,beta),
        #g3_suitability_constant(alpha),
        #g3_suitability_richards(p0,p1,p2,p3,p4),
    ),
    catchability_f= #A formula generated by e.g. g3a_predate_catchability_totalfleet expressing the catchability (i.e. total biomass) for that fleet.
      #g3a_predate_catchability_totalfleet(E), #E=Biomass caught by fleet. A formula that returns the total biomass a stock can harvest in the current time/area, generally defined by a g3_timeareadata table.
      g3a_predate_catchability_numberfleet(
        E=g3_timeareadata("ven_fleet3_landings",ven_fleet3_landings,value_field="number")), #Formula: ?????. E=Numbers caught by fleet.
      #g3a_predate_catchability_linearfleet(E),
      #g3a_predate_catchability_effortfleet(catchability_fs,E),
      #g3a_predate_catchability_quotafleet(quota_table,E,sum_stocks=list(),recalc_f=NULL),
    overconsumption_f= #Overconsumption rule, a formula that should cap all values in stock__consratio to <= 95.
      quote(logspace_add_vec(stock__consratio*-1e3,0.95*-1e3)/-1e3),
    run_f=~TRUE, #Formula specifying a condition for running this action, default always runs.
    run_at=ven_fleet_order) #Integer order that actions will be run within model.
)

###############################################################
######################### Likelihood ##########################
###############################################################

#Define likelihood actions.
ven_likelihood_actions<-list(
  g3l_understocking( #Add rates of understocking in a g3 model to nll.
    #The model report will contain nll_understocking__wgt, the results of the formula below.
    #If nll_breakdown is TRUE, this will be an array with one entry per timestep.
    #An action (i.e. list of formula objects) that will: 
    #1. Sum the total biomass adjustment due to overstocking for each prey according to the formula: l=sum_{time} sum_{areas} (sum_{prey_stocks} U_trs)^p. 
    #Where p is the power coefficient from power_f, U_trs is the total biomass adjustment to predator consumption due to overconsumption.
    prey_stocks=ven_stocks, #A list of g3_stock objects to collect catch data for.
    power_f=~2, #A formula representing power coefficient p to use.
    nll_breakdown=TRUE, #Should the nll report be broken down by time? TRUE/FALSE.
    weight=10, #Weighting applied to this likelihood component.
    run_at=likelihood_order), #Integer order that actions will be run within model.
  
  #Use g3l_catchdistribution for fleet catches and survey indices, and g3l_abundancedistribution for survey indices or stock distribution (with fleets=list(), empty).
  #Fleet catch distribution.
  g3l_catchdistribution( #Gather nll in a g3 model.
    #An action (i.e. list of formula objects) that will:
    #1. For all fleets and stocks combinations, collect catch data into modelstock__num or modelstock__wgt, depending on the columns provided in obs_data.
    #2. Compare modelstock__num/wgt with obsstock__num/wgt, using function_f.
    #The output of function_f is summed over all stock dimensions (age/area) and time and added to nll.
    nll_name="ven_catchdist1", #Character string, used to define the variable name for obsstock and modelstock.
    obs_data=ven_catchdist1, #Data.frame of observation data, for example the results of mfdb_sample_count. Should at least have a year column, and a length or weight column. For more information, see "obs_data and data aggregation" below.
    fleets=list(ven_fleet1,ven_fleet2), #A list of g3_stock objects to collect catch data for. If empty, will collect abundance data for stocks instead.
    stocks=ven_stocks, #A list of g3_stock objects to collect catch or abundance data for, depending if stocks were provided.
    function_f= #A formula to compare obsstock__x to modelstock__x and generate nll, defined by one of the g3l_distribution_* functions. This will be adapted to compare either number (modelstock__num) or weight (modelstock__wgt) depending on what columns obs_data has.
      g3l_distribution_sumofsquares( #Formula: sum_{lengths}(Ntral/Ntr-vtral/vtr)^2. Ntral=Observation sample size for current time/area/age/length combination; vtral=Model sample size for current time/area/age/length combination; Ntr=Total observation sample size for current time/area (or dimensions set in over); vtr=Total model sample size for current time/area (or dimensions set in over).
        over=c("area")), #When comparing proportions of lengthgroups, specifies the dimensions that define the total. For example the default "area" means the proprtion of the current lengthgroup to all individuals in that area.
    transform_fs=list(), #A list of dimension name to formula to apply to model data before collating.
    missing_val=0, #Where there are missing values in the incoming data, value to replace them with.
    area_group=NULL, #mfdb_group or list mapping area names used in obs_data to integer model areas.
    report=TRUE, #If TRUE, add model and observation arrays to the model report, called cdist_nll_name_model__num/wgt and cdist_nll_name_obs__num/wgt respectively.
    nll_breakdown=TRUE, #Should the nll report be broken down by time? TRUE/FALSE.
    weight= #Weighting applied to this likelihood component. Default is a g3_param that defaults to 1, allowing weights to be altered without recompiling.
      substitute(g3_param(n,optimise=FALSE,value=1),
                 list(n=paste0("ven_catchdist1","_weight"))),
    run_at=likelihood_order), #Integer order that actions will be run within model.
  
  #Fleet catch distribution.
  g3l_catchdistribution( #Gather nll in a g3 model.
    #An action (i.e. list of formula objects) that will:
    #1. For all fleets and stocks combinations, collect catch data into modelstock__num or modelstock__wgt, depending on the columns provided in obs_data.
    #2. Compare modelstock__num/wgt with obsstock__num/wgt, using function_f.
    #The output of function_f is summed over all stock dimensions (age/area) and time and added to nll.
    nll_name="ven_catchdist2", #Character string, used to define the variable name for obsstock and modelstock.
    obs_data=ven_catchdist2, #Data.frame of observation data, for example the results of mfdb_sample_count. Should at least have a year column, and a length or weight column. For more information, see "obs_data and data aggregation" below.
    fleets=list(ven_fleet1,ven_fleet2), #A list of g3_stock objects to collect catch data for. If empty, will collect abundance data for stocks instead.
    stocks=ven_stocks, #A list of g3_stock objects to collect catch or abundance data for, depending if stocks were provided.
    function_f= #A formula to compare obsstock__x to modelstock__x and generate nll, defined by one of the g3l_distribution_* functions. This will be adapted to compare either number (modelstock__num) or weight (modelstock__wgt) depending on what columns obs_data has.
      g3l_distribution_sumofsquares( #Formula: sum_{lengths}(Ntral/Ntr-vtral/vtr)^2. Ntral=Observation sample size for current time/area/age/length combination; vtral=Model sample size for current time/area/age/length combination; Ntr=Total observation sample size for current time/area (or dimensions set in over); vtr=Total model sample size for current time/area (or dimensions set in over).
        over=c("area")), #When comparing proportions of lengthgroups, specifies the dimensions that define the total. For example the default "area" means the proprtion of the current lengthgroup to all individuals in that area.
    transform_fs=list(), #A list of dimension name to formula to apply to model data before collating.
    missing_val=0, #Where there are missing values in the incoming data, value to replace them with.
    area_group=NULL, #mfdb_group or list mapping area names used in obs_data to integer model areas.
    report=TRUE, #If TRUE, add model and observation arrays to the model report, called cdist_nll_name_model__num/wgt and cdist_nll_name_obs__num/wgt respectively.
    nll_breakdown=TRUE, #Should the nll report be broken down by time? TRUE/FALSE.
    weight= #Weighting applied to this likelihood component. Default is a g3_param that defaults to 1, allowing weights to be altered without recompiling.
      substitute(g3_param(n,optimise=FALSE,value=1),
                 list(n=paste0("ven_catchdist2","_weight"))),
    run_at=likelihood_order), #Integer order that actions will be run within model.
  
  #Fleet catch distribution.
  g3l_catchdistribution( #Gather nll in a g3 model.
    #An action (i.e. list of formula objects) that will:
    #1. For all fleets and stocks combinations, collect catch data into modelstock__num or modelstock__wgt, depending on the columns provided in obs_data.
    #2. Compare modelstock__num/wgt with obsstock__num/wgt, using function_f.
    #The output of function_f is summed over all stock dimensions (age/area) and time and added to nll.
    nll_name="ven_catchdist3", #Character string, used to define the variable name for obsstock and modelstock.
    obs_data=ven_catchdist3, #Data.frame of observation data, for example the results of mfdb_sample_count. Should at least have a year column, and a length or weight column. For more information, see "obs_data and data aggregation" below.
    fleets=list(ven_fleet3), #A list of g3_stock objects to collect catch data for. If empty, will collect abundance data for stocks instead.
    stocks=ven_stocks, #A list of g3_stock objects to collect catch or abundance data for, depending if stocks were provided.
    function_f= #A formula to compare obsstock__x to modelstock__x and generate nll, defined by one of the g3l_distribution_* functions. This will be adapted to compare either number (modelstock__num) or weight (modelstock__wgt) depending on what columns obs_data has.
      g3l_distribution_sumofsquares( #Formula: sum_{lengths}(Ntral/Ntr-vtral/vtr)^2. Ntral=Observation sample size for current time/area/age/length combination; vtral=Model sample size for current time/area/age/length combination; Ntr=Total observation sample size for current time/area (or dimensions set in over); vtr=Total model sample size for current time/area (or dimensions set in over).
        over=c("area")), #When comparing proportions of lengthgroups, specifies the dimensions that define the total. For example the default "area" means the proprtion of the current lengthgroup to all individuals in that area.
    transform_fs=list(), #A list of dimension name to formula to apply to model data before collating.
    missing_val=0, #Where there are missing values in the incoming data, value to replace them with.
    area_group=NULL, #mfdb_group or list mapping area names used in obs_data to integer model areas.
    report=TRUE, #If TRUE, add model and observation arrays to the model report, called cdist_nll_name_model__num/wgt and cdist_nll_name_obs__num/wgt respectively.
    nll_breakdown=TRUE, #Should the nll report be broken down by time? TRUE/FALSE.
    weight= #Weighting applied to this likelihood component. Default is a g3_param that defaults to 1, allowing weights to be altered without recompiling.
      substitute(g3_param(n,optimise=FALSE,value=1),
                 list(n=paste0("ven_catchdist3","_weight"))),
    run_at=likelihood_order), #Integer order that actions will be run within model.
  
  #Fleet catch distribution.
  g3l_catchdistribution( #Gather nll in a g3 model.
    #An action (i.e. list of formula objects) that will:
    #1. For all fleets and stocks combinations, collect catch data into modelstock__num or modelstock__wgt, depending on the columns provided in obs_data.
    #2. Compare modelstock__num/wgt with obsstock__num/wgt, using function_f.
    #The output of function_f is summed over all stock dimensions (age/area) and time and added to nll.
    nll_name="ven_catchdist4", #Character string, used to define the variable name for obsstock and modelstock.
    obs_data=ven_catchdist4, #Data.frame of observation data, for example the results of mfdb_sample_count. Should at least have a year column, and a length or weight column. For more information, see "obs_data and data aggregation" below.
    fleets=list(ven_fleet3), #A list of g3_stock objects to collect catch data for. If empty, will collect abundance data for stocks instead.
    stocks=ven_stocks, #A list of g3_stock objects to collect catch or abundance data for, depending if stocks were provided.
    function_f= #A formula to compare obsstock__x to modelstock__x and generate nll, defined by one of the g3l_distribution_* functions. This will be adapted to compare either number (modelstock__num) or weight (modelstock__wgt) depending on what columns obs_data has.
      g3l_distribution_sumofsquares( #Formula: sum_{lengths}(Ntral/Ntr-vtral/vtr)^2. Ntral=Observation sample size for current time/area/age/length combination; vtral=Model sample size for current time/area/age/length combination; Ntr=Total observation sample size for current time/area (or dimensions set in over); vtr=Total model sample size for current time/area (or dimensions set in over).
        over=c("area")), #When comparing proportions of lengthgroups, specifies the dimensions that define the total. For example the default "area" means the proprtion of the current lengthgroup to all individuals in that area.
    transform_fs=list(), #A list of dimension name to formula to apply to model data before collating.
    missing_val=0, #Where there are missing values in the incoming data, value to replace them with.
    area_group=NULL, #mfdb_group or list mapping area names used in obs_data to integer model areas.
    report=TRUE, #If TRUE, add model and observation arrays to the model report, called cdist_nll_name_model__num/wgt and cdist_nll_name_obs__num/wgt respectively.
    nll_breakdown=TRUE, #Should the nll report be broken down by time? TRUE/FALSE.
    weight= #Weighting applied to this likelihood component. Default is a g3_param that defaults to 1, allowing weights to be altered without recompiling.
      substitute(g3_param(n,optimise=FALSE,value=1),
                 list(n=paste0("ven_catchdist4","_weight"))),
    run_at=likelihood_order), #Integer order that actions will be run within model.
  
  #Survey indices distribution.
  g3l_abundancedistribution( #Gather nll in a g3 model.
    #An action (i.e. list of formula objects) that will:
    #1. For all fleets and stocks combinations, collect catch data into modelstock__num or modelstock__wgt, depending on the columns provided in obs_data.
    #2. Compare modelstock__num/wgt with obsstock__num/wgt, using function_f.
    #The output of function_f is summed over all stock dimensions (age/area) and time and added to nll.
    nll_name="ven_surveyindices1", #Character string, used to define the variable name for obsstock and modelstock.
    obs_data=ven_surveyindices1, #Data.frame of observation data, for example the results of mfdb_sample_count. Should at least have a year column, and a length or weight column. For more information, see "obs_data and data aggregation" below.
    fleets=list(), #A list of g3_stock objects to collect catch data for. If empty, will collect abundance data for stocks instead.
    stocks=ven_stocks, #A list of g3_stock objects to collect catch or abundance data for, depending if stocks were provided.
    function_f= #A formula to compare obsstock__x to modelstock__x and generate nll, defined by one of the g3l_distribution_* functions. This will be adapted to compare either number (modelstock__num) or weight (modelstock__wgt) depending on what columns obs_data has.
      g3l_distribution_surveyindices_log( #Formula: sum_{lengths}(alpha+beta*Ntral-vtral)^2. Ntral=Observation sample size for current time/area/age/length combination; vtral=Model sample size for current time/area/age/length combination.
        alpha=NULL, #Formula substituted into surveyindices calcuations to fix intercept of linear regression, or NULL if not fixed.
        beta=NULL), #Formula substituted into surveyindices calcuations to fix slope of linear regression, or NULL if not fixed.
    transform_fs=list(), #A list of dimension name to formula to apply to model data before collating.
    missing_val=0, #Where there are missing values in the incoming data, value to replace them with.
    area_group=NULL, #mfdb_group or list mapping area names used in obs_data to integer model areas.
    report=TRUE, #If TRUE, add model and observation arrays to the model report, called cdist_nll_name_model__num/wgt and cdist_nll_name_obs__num/wgt respectively.
    nll_breakdown=TRUE, #Should the nll report be broken down by time? TRUE/FALSE.
    weight= #Weighting applied to this likelihood component. Default is a g3_param that defaults to 1, allowing weights to be altered without recompiling.
      substitute(g3_param(n,optimise=FALSE,value=1),
                 list(n=paste0("ven_surveyindices1","_weight"))),
    run_at=likelihood_order), #Integer order that actions will be run within model.
  
  #Survey indices distribution.
  g3l_abundancedistribution( #Gather nll in a g3 model.
    #An action (i.e. list of formula objects) that will:
    #1. For all fleets and stocks combinations, collect catch data into modelstock__num or modelstock__wgt, depending on the columns provided in obs_data.
    #2. Compare modelstock__num/wgt with obsstock__num/wgt, using function_f.
    #The output of function_f is summed over all stock dimensions (age/area) and time and added to nll.
    nll_name="ven_surveyindices2", #Character string, used to define the variable name for obsstock and modelstock.
    obs_data=ven_surveyindices2, #Data.frame of observation data, for example the results of mfdb_sample_count. Should at least have a year column, and a length or weight column. For more information, see "obs_data and data aggregation" below.
    fleets=list(), #A list of g3_stock objects to collect catch data for. If empty, will collect abundance data for stocks instead.
    stocks=ven_stocks, #A list of g3_stock objects to collect catch or abundance data for, depending if stocks were provided.
    function_f= #A formula to compare obsstock__x to modelstock__x and generate nll, defined by one of the g3l_distribution_* functions. This will be adapted to compare either number (modelstock__num) or weight (modelstock__wgt) depending on what columns obs_data has.
      g3l_distribution_surveyindices_log( #Formula: sum_{lengths}(alpha+beta*log(Ntral)-log(vtral))^2. Ntral=Observation sample size for current time/area/age/length combination; vtral=Model sample size for current time/area/age/length combination.
        alpha=NULL, #Formula substituted into surveyindices calcuations to fix intercept of linear regression, or NULL if not fixed.
        beta=NULL), #Formula substituted into surveyindices calcuations to fix slope of linear regression, or NULL if not fixed.
    transform_fs=list(), #A list of dimension name to formula to apply to model data before collating.
    missing_val=0, #Where there are missing values in the incoming data, value to replace them with.
    area_group=NULL, #mfdb_group or list mapping area names used in obs_data to integer model areas.
    report=TRUE, #If TRUE, add model and observation arrays to the model report, called cdist_nll_name_model__num/wgt and cdist_nll_name_obs__num/wgt respectively.
    nll_breakdown=TRUE, #Should the nll report be broken down by time? TRUE/FALSE.
    weight= #Weighting applied to this likelihood component. Default is a g3_param that defaults to 1, allowing weights to be altered without recompiling.
      substitute(g3_param(n,optimise=FALSE,value=1),
                 list(n=paste0("ven_surveyindices2","_weight"))),
    run_at=likelihood_order), #Integer order that actions will be run within model.
  
  #Stock distribution.
  g3l_abundancedistribution(
    #Assuming fleets is empty, an action (i.e. list of formula objects) that will:
    #1. For all stocks, collect catch data into modelstock__num or modelstock__wgt, depending on the columns provided in obs_data.
    #2. Compare modelstock__num/wgt with obsstock__num/wgt, using function_f.
    #The output of function_f is summed over all stock dimensions (age/area) and time and added to nll.
    nll_name="ven_stockdist1", #Character string, used to define the variable name for obsstock and modelstock.
    obs_data=ven_stockdist1, #Data.frame of observation data, for example the results of mfdb_sample_count. Should at least have a year column, and a length or weight column. For more information, see "obs_data and data aggregation" below.
    fleets=list(), #A list of g3_stock objects to collect catch data for. If empty, will collect abundance data for stocks instead.
    stocks=ven_stocks, #A list of g3_stock objects to collect catch or abundance data for, depending if stocks were provided.
    function_f= #A formula to compare obsstock__x to modelstock__x and generate nll, defined by one of the g3l_distribution_* functions. This will be adapted to compare either number (modelstock__num) or weight (modelstock__wgt) depending on what columns obs_data has.
      g3l_distribution_sumofsquares( #Formula: sum_{lengths}(Ntral/Ntr-vtral/vtr)^2. Ntral=Observation sample size for current time/area/age/length combination; vtral=Model sample size for current time/area/age/length combination; Ntr=Total observation sample size for current time/area (or dimensions set in over); vtr=Total model sample size for current time/area (or dimensions set in over).
        over=c("area")), #When comparing proportions of lengthgroups, specifies the dimensions that define the total. For example the default "area" means the proprtion of the current lengthgroup to all individuals in that area.
    transform_fs=list(), #A list of dimension name to formula to apply to model data before collating.
    missing_val=0, #Where there are missing values in the incoming data, value to replace them with.
    area_group=NULL, #mfdb_group or list mapping area names used in obs_data to integer model areas.
    report=TRUE, #If TRUE, add model and observation arrays to the model report, called cdist_nll_name_model__num/wgt and cdist_nll_name_obs__num/wgt respectively.
    nll_breakdown=TRUE, #Should the nll report be broken down by time? TRUE/FALSE.
    weight= #Weighting applied to this likelihood component. Default is a g3_param that defaults to 1, allowing weights to be altered without recompiling.
      substitute(g3_param(n,optimise=FALSE,value=1),
                 list(n=paste0("ven_stockdist1","_weight"))),
    run_at=likelihood_order), #Integer order that actions will be run within model.
  
  #Stock distribution.
  g3l_abundancedistribution(
    #Assuming fleets is empty, an action (i.e. list of formula objects) that will:
    #1. For all stocks, collect catch data into modelstock__num or modelstock__wgt, depending on the columns provided in obs_data.
    #2. Compare modelstock__num/wgt with obsstock__num/wgt, using function_f.
    #The output of function_f is summed over all stock dimensions (age/area) and time and added to nll.
    nll_name="ven_stockdist2", #Character string, used to define the variable name for obsstock and modelstock.
    obs_data=ven_stockdist2, #Data.frame of observation data, for example the results of mfdb_sample_count. Should at least have a year column, and a length or weight column. For more information, see "obs_data and data aggregation" below.
    fleets=list(), #A list of g3_stock objects to collect catch data for. If empty, will collect abundance data for stocks instead.
    stocks=ven_stocks, #A list of g3_stock objects to collect catch or abundance data for, depending if stocks were provided.
    function_f= #A formula to compare obsstock__x to modelstock__x and generate nll, defined by one of the g3l_distribution_* functions. This will be adapted to compare either number (modelstock__num) or weight (modelstock__wgt) depending on what columns obs_data has.
      g3l_distribution_sumofsquares( #Formula: sum_{lengths}(Ntral/Ntr-vtral/vtr)^2. Ntral=Observation sample size for current time/area/age/length combination; vtral=Model sample size for current time/area/age/length combination; Ntr=Total observation sample size for current time/area (or dimensions set in over); vtr=Total model sample size for current time/area (or dimensions set in over).
        over=c("area")), #When comparing proportions of lengthgroups, specifies the dimensions that define the total. For example the default "area" means the proprtion of the current lengthgroup to all individuals in that area.
    transform_fs=list(), #A list of dimension name to formula to apply to model data before collating.
    missing_val=0, #Where there are missing values in the incoming data, value to replace them with.
    area_group=NULL, #mfdb_group or list mapping area names used in obs_data to integer model areas.
    report=TRUE, #If TRUE, add model and observation arrays to the model report, called cdist_nll_name_model__num/wgt and cdist_nll_name_obs__num/wgt respectively.
    nll_breakdown=TRUE, #Should the nll report be broken down by time? TRUE/FALSE.
    weight= #Weighting applied to this likelihood component. Default is a g3_param that defaults to 1, allowing weights to be altered without recompiling.
      substitute(g3_param(n,optimise=FALSE,value=1),
                 list(n=paste0("ven_stockdist2","_weight"))),
    run_at=likelihood_order)#, #Integer order that actions will be run within model.
  
  # #Stock distribution.
  # g3l_abundancedistribution(
  #   #Assuming fleets is empty, an action (i.e. list of formula objects) that will:
  #   #1. For all stocks, collect catch data into modelstock__num or modelstock__wgt, depending on the columns provided in obs_data.
  #   #2. Compare modelstock__num/wgt with obsstock__num/wgt, using function_f.
  #   #The output of function_f is summed over all stock dimensions (age/area) and time and added to nll.
  #   nll_name="ven_imm_stockdist1", #Character string, used to define the variable name for obsstock and modelstock.
  #   obs_data=ven_imm_stockdist1, #Data.frame of observation data, for example the results of mfdb_sample_count. Should at least have a year column, and a length or weight column. For more information, see "obs_data and data aggregation" below.
  #   fleets=list(), #A list of g3_stock objects to collect catch data for. If empty, will collect abundance data for stocks instead.
  #   stocks=list(ven_imm), #A list of g3_stock objects to collect catch or abundance data for, depending if stocks were provided.
  #   function_f= #A formula to compare obsstock__x to modelstock__x and generate nll, defined by one of the g3l_distribution_* functions. This will be adapted to compare either number (modelstock__num) or weight (modelstock__wgt) depending on what columns obs_data has.
  #     g3l_distribution_sumofsquares( #Formula: sum_{lengths}(Ntral/Ntr-vtral/vtr)^2. Ntral=Observation sample size for current time/area/age/length combination; vtral=Model sample size for current time/area/age/length combination; Ntr=Total observation sample size for current time/area (or dimensions set in over); vtr=Total model sample size for current time/area (or dimensions set in over).
  #       over=c("area")), #When comparing proportions of lengthgroups, specifies the dimensions that define the total. For example the default "area" means the proprtion of the current lengthgroup to all individuals in that area.
  #   transform_fs=list(), #A list of dimension name to formula to apply to model data before collating.
  #   missing_val=0, #Where there are missing values in the incoming data, value to replace them with.
  #   area_group=NULL, #mfdb_group or list mapping area names used in obs_data to integer model areas.
  #   report=TRUE, #If TRUE, add model and observation arrays to the model report, called cdist_nll_name_model__num/wgt and cdist_nll_name_obs__num/wgt respectively.
  #   nll_breakdown=TRUE, #Should the nll report be broken down by time? TRUE/FALSE.
  #   weight= #Weighting applied to this likelihood component. Default is a g3_param that defaults to 1, allowing weights to be altered without recompiling.
  #     substitute(g3_param(n,optimise=FALSE,value=1),
  #                list(n=paste0("ven_imm_stockdist1","_weight"))),
  #   run_at=likelihood_order), #Integer order that actions will be run within model.
  # 
  # #Stock distribution.
  # g3l_abundancedistribution(
  #   #Assuming fleets is empty, an action (i.e. list of formula objects) that will:
  #   #1. For all stocks, collect catch data into modelstock__num or modelstock__wgt, depending on the columns provided in obs_data.
  #   #2. Compare modelstock__num/wgt with obsstock__num/wgt, using function_f.
  #   #The output of function_f is summed over all stock dimensions (age/area) and time and added to nll.
  #   nll_name="ven_mat_stockdist1", #Character string, used to define the variable name for obsstock and modelstock.
  #   obs_data=ven_mat_stockdist1, #Data.frame of observation data, for example the results of mfdb_sample_count. Should at least have a year column, and a length or weight column. For more information, see "obs_data and data aggregation" below.
  #   fleets=list(), #A list of g3_stock objects to collect catch data for. If empty, will collect abundance data for stocks instead.
  #   stocks=list(ven_mat), #A list of g3_stock objects to collect catch or abundance data for, depending if stocks were provided.
  #   function_f= #A formula to compare obsstock__x to modelstock__x and generate nll, defined by one of the g3l_distribution_* functions. This will be adapted to compare either number (modelstock__num) or weight (modelstock__wgt) depending on what columns obs_data has.
  #     g3l_distribution_sumofsquares( #Formula: sum_{lengths}(Ntral/Ntr-vtral/vtr)^2. Ntral=Observation sample size for current time/area/age/length combination; vtral=Model sample size for current time/area/age/length combination; Ntr=Total observation sample size for current time/area (or dimensions set in over); vtr=Total model sample size for current time/area (or dimensions set in over).
  #       over=c("area")), #When comparing proportions of lengthgroups, specifies the dimensions that define the total. For example the default "area" means the proprtion of the current lengthgroup to all individuals in that area.
  #   transform_fs=list(), #A list of dimension name to formula to apply to model data before collating.
  #   missing_val=0, #Where there are missing values in the incoming data, value to replace them with.
  #   area_group=NULL, #mfdb_group or list mapping area names used in obs_data to integer model areas.
  #   report=TRUE, #If TRUE, add model and observation arrays to the model report, called cdist_nll_name_model__num/wgt and cdist_nll_name_obs__num/wgt respectively.
  #   nll_breakdown=TRUE, #Should the nll report be broken down by time? TRUE/FALSE.
  #   weight= #Weighting applied to this likelihood component. Default is a g3_param that defaults to 1, allowing weights to be altered without recompiling.
  #     substitute(g3_param(n,optimise=FALSE,value=1),
  #                list(n=paste0("ven_mat_stockdist1","_weight"))),
  #   run_at=likelihood_order), #Integer order that actions will be run within model.
  # 
  # #Stock distribution.
  # g3l_abundancedistribution(
  #   #Assuming fleets is empty, an action (i.e. list of formula objects) that will:
  #   #1. For all stocks, collect catch data into modelstock__num or modelstock__wgt, depending on the columns provided in obs_data.
  #   #2. Compare modelstock__num/wgt with obsstock__num/wgt, using function_f.
  #   #The output of function_f is summed over all stock dimensions (age/area) and time and added to nll.
  #   nll_name="ven_imm_stockdist2", #Character string, used to define the variable name for obsstock and modelstock.
  #   obs_data=ven_imm_stockdist2, #Data.frame of observation data, for example the results of mfdb_sample_count. Should at least have a year column, and a length or weight column. For more information, see "obs_data and data aggregation" below.
  #   fleets=list(), #A list of g3_stock objects to collect catch data for. If empty, will collect abundance data for stocks instead.
  #   stocks=list(ven_imm), #A list of g3_stock objects to collect catch or abundance data for, depending if stocks were provided.
  #   function_f= #A formula to compare obsstock__x to modelstock__x and generate nll, defined by one of the g3l_distribution_* functions. This will be adapted to compare either number (modelstock__num) or weight (modelstock__wgt) depending on what columns obs_data has.
  #     g3l_distribution_sumofsquares( #Formula: sum_{lengths}(Ntral/Ntr-vtral/vtr)^2. Ntral=Observation sample size for current time/area/age/length combination; vtral=Model sample size for current time/area/age/length combination; Ntr=Total observation sample size for current time/area (or dimensions set in over); vtr=Total model sample size for current time/area (or dimensions set in over).
  #       over=c("area")), #When comparing proportions of lengthgroups, specifies the dimensions that define the total. For example the default "area" means the proprtion of the current lengthgroup to all individuals in that area.
  #   transform_fs=list(), #A list of dimension name to formula to apply to model data before collating.
  #   missing_val=0, #Where there are missing values in the incoming data, value to replace them with.
  #   area_group=NULL, #mfdb_group or list mapping area names used in obs_data to integer model areas.
  #   report=TRUE, #If TRUE, add model and observation arrays to the model report, called cdist_nll_name_model__num/wgt and cdist_nll_name_obs__num/wgt respectively.
  #   nll_breakdown=TRUE, #Should the nll report be broken down by time? TRUE/FALSE.
  #   weight= #Weighting applied to this likelihood component. Default is a g3_param that defaults to 1, allowing weights to be altered without recompiling.
  #     substitute(g3_param(n,optimise=FALSE,value=1),
  #                list(n=paste0("ven_imm_stockdist2","_weight"))),
  #   run_at=likelihood_order), #Integer order that actions will be run within model.
  # 
  # #Stock distribution.
  # g3l_abundancedistribution(
  #   #Assuming fleets is empty, an action (i.e. list of formula objects) that will:
  #   #1. For all stocks, collect catch data into modelstock__num or modelstock__wgt, depending on the columns provided in obs_data.
  #   #2. Compare modelstock__num/wgt with obsstock__num/wgt, using function_f.
  #   #The output of function_f is summed over all stock dimensions (age/area) and time and added to nll.
  #   nll_name="ven_mat_stockdist2", #Character string, used to define the variable name for obsstock and modelstock.
  #   obs_data=ven_mat_stockdist2, #Data.frame of observation data, for example the results of mfdb_sample_count. Should at least have a year column, and a length or weight column. For more information, see "obs_data and data aggregation" below.
  #   fleets=list(), #A list of g3_stock objects to collect catch data for. If empty, will collect abundance data for stocks instead.
  #   stocks=list(ven_mat), #A list of g3_stock objects to collect catch or abundance data for, depending if stocks were provided.
  #   function_f= #A formula to compare obsstock__x to modelstock__x and generate nll, defined by one of the g3l_distribution_* functions. This will be adapted to compare either number (modelstock__num) or weight (modelstock__wgt) depending on what columns obs_data has.
  #     g3l_distribution_sumofsquares( #Formula: sum_{lengths}(Ntral/Ntr-vtral/vtr)^2. Ntral=Observation sample size for current time/area/age/length combination; vtral=Model sample size for current time/area/age/length combination; Ntr=Total observation sample size for current time/area (or dimensions set in over); vtr=Total model sample size for current time/area (or dimensions set in over).
  #       over=c("area")), #When comparing proportions of lengthgroups, specifies the dimensions that define the total. For example the default "area" means the proprtion of the current lengthgroup to all individuals in that area.
  #   transform_fs=list(), #A list of dimension name to formula to apply to model data before collating.
  #   missing_val=0, #Where there are missing values in the incoming data, value to replace them with.
  #   area_group=NULL, #mfdb_group or list mapping area names used in obs_data to integer model areas.
  #   report=TRUE, #If TRUE, add model and observation arrays to the model report, called cdist_nll_name_model__num/wgt and cdist_nll_name_obs__num/wgt respectively.
  #   nll_breakdown=TRUE, #Should the nll report be broken down by time? TRUE/FALSE.
  #   weight= #Weighting applied to this likelihood component. Default is a g3_param that defaults to 1, allowing weights to be altered without recompiling.
  #     substitute(g3_param(n,optimise=FALSE,value=1),
  #                list(n=paste0("ven_mat_stockdist2","_weight"))),
  #   run_at=likelihood_order) #Integer order that actions will be run within model.
)

###############################################################
####################### Random effects ########################
###############################################################

#Add likelihood components for random effects.
#The model report will contain nll_random_dnorm_dnorm_lin__dnorm, the results of applying dnorm. If nll_breakdown is TRUE, this will be an array with one entry per timestep.
# ven_random_actions<-list(
#   g3l_random_dnorm( 
#     #An action (i.e. list of formula objects) that will: 
#     #1. On the final model step, calculate dnorm(param_f, mean_f, sigma_f) & add to nll.
#     nll_name=, #Character string, used to define the variable name for dnorm output.
#     param_f=, #A formula representing the value to apply dnorm to. Invariably a g3_param for g3l_random_dnorm, a g3_param_table with cur_year for g3l_random_walk.
#     mean_f=0, #A formula representing mean in dnorm.
#     sigma_f=1, #A formula representing sigma in dnorm.
#     log_f=TRUE, #A formula representing log in dnorm.
#     weight=if(log_f){-1.0}else{1.0}, #Weighting applied to this likelihood component.
#     run_at=random_order), #Integer order that actions will be run within model.
#   
#   g3l_random_walk(
#     #An action (i.e. list of formula objects) that will:
#     #1. Calculate dnorm(param_f, previous param_f, sigma_f) (at final year if period = year).
#     #2. Add to nll.
#     nll_name=, #Character string, used to define the variable name for dnorm output.
#     param_f=, #A formula representing the value to apply dnorm to. Invariably a g3_param for g3l_random_dnorm, a g3_param_table with cur_year for g3l_random_walk.
#     sigma_f=1, #A formula representing sigma in dnorm.
#     log_f=TRUE, #A formula representing log in dnorm.
#     period="year", #When dnorm should be recalculated. Once per year or every step.
#     nll_breakdown=FALSE, #Should the nll report be broken down by time? TRUE/FALSE.
#     weight=if(log_f){-1.0}else{1.0}, #Weighting applied to this likelihood component.
#     run_at=random_order) #Integer order that actions will be run within model.
# )

###############################################################
################## Generate R and TMB codes ###################
###############################################################

#Stock actions.
ven_stock_actions<-c(ven_initial_conditions,
                     ven_natural_mortality,
                     ven_ageing,
                     ven_renewal,
                     ven_growmature,
                     ven_reproduction,
                     ven_migration)

#Collate actions.
ven_actions<-c(
  ven_time_actions,
  ven_stock_actions,
  ven_fleet_actions,
  #ven_random_actions,
  ven_likelihood_actions)

ven_actions<-c(ven_actions,
               list(g3a_report_history(ven_actions)))

#Save model actions.
save(ven_actions,
     file=file.path(dirName,"model_actions.Rdata"))

#Turn actions into an R function.
ven_r_model<-g3_to_r(ven_actions,
                     trace=TRUE,
                     strict=TRUE)

#Turn actions into TMB function.
ven_tmb_model<-g3_to_tmb(ven_actions,
                         trace=TRUE,
                         strict=TRUE)

#Save model.
save(ven_r_model,
     ven_tmb_model,
     file=file.path(dirName,"model_R_TMP_code.Rdata"))

###############################################################
##################### Fill parameter table ####################
###############################################################

#Get the parameter template.
ven_tmb_param<-attr(ven_tmb_model,"parameter_template")
ven_tmb_param

#Fill it with initial conditions.
ven_tmb_param<- 
  ven_tmb_param %>% 
  g3_init_guess("ven_imm.init.scalar",value=1e3,lower=1,upper=1e8,optimise=0) %>%
  g3_init_guess("ven_mat.init.scalar",value=1e3,lower=1,upper=1e8,optimise=0) %>%
  g3_init_guess("ven_imm.init.[0-9]",value=1,lower=0.001,upper=100,optimise=1) %>%
  g3_init_guess("ven_mat.init.[0-9]",value=1,lower=0.001,upper=100,optimise=1) %>%
  g3_init_guess("ven_imm.init.M",value=mortrate.constant,lower=0.001,upper=1,optimise=1) %>% 
  g3_init_guess("ven_mat.init.M",value=mortrate.constant,lower=0.001,upper=1,optimise=1) %>% 
  g3_init_guess("ven_imm.init.alpha",value=ven_lw.constants$a,lower=1e-10,upper=1,optimise=0) %>%
  g3_init_guess("ven_mat.init.alpha",value=ven_lw.constants$a,lower=1e-10,upper=1,optimise=0) %>%
  g3_init_guess("ven_imm.init.beta",value=ven_lw.constants$b,lower=2,upper=4,optimise=0) %>%
  g3_init_guess("ven_mat.init.beta",value=ven_lw.constants$b,lower=2,upper=4,optimise=0) %>%
  # g3_init_guess("ven_imm.growth.alpha",value=ven_lw.constants$a,lower=1e-10,upper=1,optimise=0) %>%
  # g3_init_guess("ven_mat.growth.alpha",value=ven_lw.constants$a,lower=1e-10,upper=1,optimise=0) %>%
  # g3_init_guess("ven_imm.growth.beta",value=ven_lw.constants$b,lower=2,upper=4,optimise=0) %>%
  # g3_init_guess("ven_mat.growth.beta",value=ven_lw.constants$b,lower=2,upper=4,optimise=0) %>%
  g3_init_guess("growth.alpha",value=ven_lw.constants$a,lower=1e-10,upper=1,optimise=0) %>%
  g3_init_guess("growth.beta",value=ven_lw.constants$b,lower=2,upper=4,optimise=0) %>%
  g3_init_guess("renewal.alpha",value=ven_lw.constants$a,lower=1e-10,upper=1,optimise=0) %>%
  g3_init_guess("renewal.beta",value=ven_lw.constants$b,lower=2,upper=4,optimise=0) %>%
  # g3_init_guess("ven_imm.renewal.alpha",value=ven_lw.constants$a,lower=1e-10,upper=1,optimise=0) %>%
  # g3_init_guess("ven_imm.renewal.beta",value=ven_lw.constants$b,lower=2,upper=4,optimise=0) %>%
  g3_init_guess("ven_imm.natmort.M",value=mortrate.constant,lower=0.001,upper=1,optimise=1) %>% 
  g3_init_guess("ven_mat.natmort.M",value=mortrate.constant,lower=0.001,upper=1,optimise=1) %>% 
  g3_init_guess("init.F",value=0.3,lower=0.1,upper=1,optimise=1) %>%
  g3_init_guess("init.Linf",value=ven_al.constants$Linf,lower=12,upper=25,optimise=0) %>%
  g3_init_guess("init.K",value=ven_al.constants$k,lower=0.1,upper=10,optimise=0) %>%
  g3_init_guess("init.recl",value=9,lower=5,upper=15,optimise=1) %>%
  g3_init_guess("init.sd.[0-9]",value=mean(ven_init.sigma$lengthsd),lower=0.01,upper=10,optimise=0) %>%
  g3_init_guess("growth.Linf",value=ven_al.constants$Linf,lower=12,upper=25,optimise=0) %>%
  g3_init_guess("growth.K",value=ven_al.constants$k,lower=0.1,upper=10,optimise=0) %>%
  g3_init_guess("growth.bbin",value=0.9,lower=0.001,upper=50,optimise=1) %>%
  g3_init_guess("mat.alpha",value=ven_mat.constants$a,lower=0.1,upper=10,optimise=0) %>%
  g3_init_guess("mat.l50",value=ven_mat.constants$l50,lower=1,upper=30,optimise=0) %>%
  g3_init_guess("spawn.prop.alpha",value=ven_mat.constants$a,lower=0.1,upper=10,optimise=0) %>%
  g3_init_guess("spawn.prop.l50",value=ven_mat.constants$l50,lower=1,upper=30,optimise=0) %>%
  g3_init_guess("spawn.rec.mu",value=0.4,lower=0.01,upper=1,optimise=1) %>%
  g3_init_guess("spawn.wgtloss.alpha",value=0.02,lower=0.001,upper=1,optimise=1) %>%
  g3_init_guess("spawn.mort.alpha",value=0.56,lower=0.001,upper=1,optimise=1) %>% #Karjalainen and Marjomäki (2017).
  g3_init_guess("spawn.mort.beta",value=0.0,lower=0.0,upper=1,optimise=1) %>% #Karjalainen and Marjomäki (2017).
  g3_init_guess("renewal.Linf",value=ven_al.constants$Linf,lower=12,upper=25,optimise=0) %>%
  g3_init_guess("renewal.K",value=ven_al.constants$k,lower=0.1,upper=10,optimise=0) %>%
  g3_init_guess("renewal.recl",value=9,lower=5,upper=15,optimise=1) %>%
  g3_init_guess("renewal.rec.sd",value=1.1,lower=0.01,upper=15,optimise=1) %>%
  g3_init_guess("renewal.rec.scalar",value=1e3,lower=1,upper=1e8,optimise=1) %>%
  g3_init_guess("renewal.rec.[0-9]",value=1,lower=0.001,upper=100,optimise=1) %>%
  g3_init_guess("fleet1.alpha",value=0.9,lower=0.01,upper=3,optimise=1) %>%
  g3_init_guess("fleet1.l50",value=12,lower=2,upper=20,optimise=1) %>%
  g3_init_guess("fleet2.alpha",value=0.9,lower=0.01,upper=3,optimise=1) %>%
  g3_init_guess("fleet2.l50",value=12,lower=2,upper=20,optimise=1) %>%
  g3_init_guess("fleet3.alpha",value=0.9,lower=0.01,upper=3,optimise=1) %>%
  g3_init_guess("fleet3.l50",value=12,lower=2,upper=20,optimise=1) #%>%
  #g3_init_guess("ven_imm_stockdist1_weight",value=,lower=,upper=,optimise=1) %>%
  #g3_init_guess("ven_imm_stockdist2_weight",value=,lower=,upper=,optimise=1) %>%
  #g3_init_guess("ven_mat_stockdist1_weight",value=,lower=,upper=,optimise=1) %>%
  #g3_init_guess("ven_mat_stockdist2_weight",value=,lower=,upper=,optimise=1) %>%
  #g3_init_guess("ven_surveyindices1_weight",value=,lower=,upper=,optimise=1) %>%
  #g3_init_guess("ven_surveyindices2_weight",value=,lower=,upper=,optimise=1) %>%
  #g3_init_guess("ven_catchdist1_weight",value=,lower=,upper=,optimise=1) %>%
  #g3_init_guess("ven_catchdist2_weight",value=,lower=,upper=,optimise=1) %>%
  #g3_init_guess("ven_catchdist3_weight",value=,lower=,upper=,optimise=1) %>%
  #g3_init_guess("ven_catchdist4_weight",value=,lower=,upper=,optimise=1) %>%
  #g3_init_guess("retro_years",value=0,lower=0,upper=10,optimise=0) %>%
  #g3_init_guess("project_years",value=0,lower=0,upper=10,optimise=0) 
ven_tmb_param

save(ven_tmb_param,
     file=file.path(dirName,"model_paramvalues.Rdata"))

###############################################################
###################### Run and optimize #######################
###############################################################

#Run R model with provided parameters.
ven_r_model_output<-ven_r_model(ven_tmb_param$value)
ven_r_model_output
#ven_r_model_output[[1]]
#ven_tmb_param$value$project_years<-10
#ven_tmb_param$value$tac.hr.high<-0.25
#ven_tmb_param$value$tac.trigger<-9e6
#ven_r_model_output<-ven_r_model(ven_tmb_param$value)

#Get the report from the last model run.
ven_r_model_report<-attributes(ven_r_model_output)

#List all available reports.
print(names(attributes(ven_r_model_output)))

#Compile and generate TMB ADFun (see ?TMB::MakeADFun).
ven_tmb_obj.fun<-g3_tmb_adfun(ven_tmb_model,parameters=ven_tmb_param,
                              compile_flags=
                                if(.Platform$OS.type=="windows"){c("-O1","-march=native")}
                              else{c("-O3","-flto","-g")}) #Inform Bjarki about adding -g.

#Run model once, using g3_tmb_par to reshape tmb_param into param vector.
#Get nll (negative log likelihood). We maximize the probability of choosing the correct category by minimizing the negative log likelihood.
ven_tmb_obj.fun$fn(g3_tmb_par(ven_tmb_param))

#Run model once, returning model report.
examplereport<-ven_tmb_obj.fun$report(g3_tmb_par(ven_tmb_param))
examplereport

#Run model through R optimiser, using bounds set in ven_tmb_param.
ven_r_model_fit.opt<-optim(ven_tmb_obj.fun$par,
                           ven_tmb_obj.fun$fn,
                           ven_tmb_obj.fun$gr,
                           lower=gadget3:::g3_tmb_bound(ven_tmb_param,"lower",include_random=T),
                           upper=gadget3:::g3_tmb_bound(ven_tmb_param,"upper",include_random=T),
                           method="L-BFGS-B",
                           control=list(trace=6,
                                        maxit=100,
                                        factr=.Machine$double.eps^2)) 
ven_r_model_fit.opt$value

#Save model optimization output.
save(ven_r_model_fit.opt,
     file=file.path(dirName,"model_optim_output.Rdata"))

#Fit model with optimized values.
ven_r_model_fit<-gadgetutils::g3_fit(ven_r_model,
                                     g3_tmb_relist(ven_tmb_param,ven_r_model_fit.opt$par))

#Save model fit output.
save(ven_r_model_fit,
     file=file.path(dirName,"model_fit_output.Rdata"))

#system("say Job finished!")
