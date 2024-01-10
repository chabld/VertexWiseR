## FUNCTION FOR MIXED-EFFECT VERTEX-WISE TFCE ANALYSIS
## FOR USE IN THE COGNITIVE AND BRAIN HEALTH LABORATORY

############################################################################################################################
############################################################################################################################
##Main function

TFCE.vertex_analysis.mixed=function(model,contrast, CT_data, random, nperm=100, tail=2, nthread=10, smooth_FWHM, perm_within_between=F)
{
  ##load other TFCE and vertex-wise functions
  source("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/vertTFCE.r?raw=TRUE")
  source("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/otherfunc.r?raw=TRUE")
  
  ##checks
    #check random variable and recode to numeric
    if(missing("random"))  {stop("random variable is missing")}
    else 
    { #recoding subject variable
      random=dat_beh$SUB_ID
      random=match(random,unique(random))
    }
    
    #check if required packages are installed
    packages=c("foreach","doParallel","parallel","doSNOW","reticulate")
    new.packages = packages[!(packages %in% installed.packages()[,"Package"])]
    if(length(new.packages)) 
    {
      cat(paste("The following package(s) are required and will be installed:\n",new.packages,"\n"))
      install.packages(new.packages)
    }  
    #check if nrow is consistent for model and CT_data
    if(NROW(CT_data)!=NROW(model))  {stop(paste("The number of rows for CT_data (",NROW(CT_data),") and model (",NROW(model),") are not the same",sep=""))}
    
    #incomplete data check
    idxF=which(complete.cases(model)==F)
    if(length(idxF)>0)
    {
      cat(paste("model contains",length(idxF),"subjects with incomplete data. Subjects with incomplete data will be excluded in the current analysis"))
      model=model[-idxF,]
      contrast=contrast[-idxF]
      CT_data=CT_data[-idxF,]
    }
    
    #check contrast
    for(colno in 1:(NCOL(model)+1))
    {
      if(colno==(NCOL(model)+1))  {stop("contrast is not contained within model")}
      
      if(class(contrast) != "integer" & class(contrast) != "numeric") 
      {
        if(identical(contrast,model[,colno]))  {break} 
      } else 
      {
        if(identical(as.numeric(contrast),as.numeric(model[,colno])))  {break}
      }
    }
    
    #check and recode categorical variables
    if(NCOL(model)>1)
    {
      for (column in 1:NCOL(model))
      {
        if(class(model[,column]) != "integer" & class(model[,column]) != "numeric")
        {
          if(length(unique(model[,column]))==2)
          {
            cat(paste("The binary variable '",colnames(model)[column],"' will be recoded such that ",unique(model[,column])[1],"=0 and ",unique(model[,column])[2],"=1 for the analysis\n",sep=""))
            
            recode=rep(0,NROW(model))
            recode[model[,column]==unique(model[,column])[2]]=1
            model[,column]=recode
            IV_of_interest=model[,colno]
          } else if(length(unique(model[,column]))>2)    {stop(paste("The categorical variable '",colnames(model)[column],"' contains more than 2 levels, please code it into binarized dummy variables",sep=""))}
        }      
      }
    } else
    {
      for (column in 1:NCOL(model))
      {
        if(class(model[column]) != "integer" & class(model[column]) != "numeric")
        {
          if(length(unique(model[column]))==2)
          {
            cat(paste("The binary variable '",colnames(model)[column],"' will be recoded such that ",unique(model[column])[1],"=0 and ",unique(model[column])[2],"=1 for the analysis\n",sep=""))
            
            recode=rep(0,NROW(model))
            recode[model[column]==unique(model[column])[2]]=1
            model[,column]=recode
            IV_of_interest=model[,colno]
          } else if(length(unique(model[column]))>2)    {stop(paste("The categorical variable '",colnames(model)[column],"' contains more than 2 levels, please code it into binarized dummy variables",sep=""))}
        }      
      }
    }
    
    #check length of CT data and load the appropriate fsaverage files
    n_vert=ncol(CT_data)
    if(n_vert==20484)
    {
      template="fsaverage5"
      load(file = url("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/data/ROImap_fs5.rdata?raw=TRUE"))
      load(file = url("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/data/edgelistfs5.rdata?raw=TRUE"),envir = globalenv())
    }
    else if (n_vert==81924)
    {
      template="fsaverage6"
      load(file = url("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/data/ROImap_fs6.rdata?raw=TRUE"))
      load(file = url("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/data/edgelistfs6.rdata?raw=TRUE"),envir = globalenv())
    }
    else {stop("CT_data should only contain 20484 (fsaverage5) or 81924 (fsaverage6) columns")}
    
    #check for collinearity
    if(NCOL(model)>1)
    {
      cormat=cor(model,use = "pairwise.complete.obs")
      cormat.0=cormat
      cormat.0[cormat.0==1]=NA
      if(max(abs(cormat.0),na.rm = T) >0.5)
      {
        warning(paste("correlations among variables in model are observed to be as high as ",round(max(abs(cormat.0),na.rm = T),2),", suggesting potential collinearity among predictors.\nAnalysis will continue...\n",sep=""))
      }
    }
  
  ##smoothing
  n_vert=NCOL(CT_data)
  if(missing("smooth_FWHM"))
  {
    if(n_vert==20484) 
    {
      reticulate::source_python("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/smooth.py?raw=TRUE")
      cat("CT_data will be smoothed using the default 10mm FWHM kernel for fsaverage5 images\n")
      CT_data=mesh_smooth(CT_data, FWHM=10)
    }
    else if(n_vert==81924) 
    {
      reticulate::source_python("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/smooth.py?raw=TRUE")
      cat("CT_data will be smoothed using the default 5mm FWHM kernel for fsaverage6 images")
      CT_data=mesh_smoothsmooth(CT_data, FWHM=5)
    }
  } else if(smooth_FWHM>0) 
  {
    reticulate::source_python("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/smooth.py?raw=TRUE")
    cat(paste("CT_data will be smoothed using a ", smooth,"mm FWHM kernel", sep=""))
    CT_data=mesh_smooth(CT_data, FWHM=smooth_FWHM)
  }
  
  ##unpermuted model
    #preparing mask for model
    mask=array(rep(T,NCOL(CT_data)))
    maskNA=which(colSums(CT_data != 0) == 0)
    mask[which(colSums(CT_data != 0) == 0)]=F
    
    #construct model
    start=Sys.time()
    cat("Estimating unpermuted TFCE image...")
    brainstat.stats.terms=reticulate::import("brainstat.stats.terms")
    brainstat.stats.SLM=reticulate::import("brainstat.stats.SLM")
    terms=brainstat.stats.terms$MixedEffect(ran = random,fix = model,"_check_categorical" = F)
    model.fit=brainstat.stats.SLM$SLM(model = terms,
                                      contrast=contrast,
                                      surf = template, 
                                      mask=mask,
                                      correction="None",
                                      cluster_threshold=1)
    #fit model
    model.fit$fit(CT_data)
    
    #save output from model
    tmap.orig=as.numeric(model.fit$t)
    TFCE.orig=TFCE.multicore(tmap.orig,tail=2,nthread=10)
    
    end=Sys.time()
    
    cat(paste("Completed in",round(difftime(end,start, units="secs"),1),"secs\nEstimating permuted TFCE images...\n",sep=" "))
  
  ##permuted model
    #generating permutation sequences  
    set.seed(123)
    permseq=matrix(NA, nrow=NROW(model), ncol=nperm)
    
    if(perm_within_between==T) {for (perm in 1:nperm)  {permseq[,perm]=perm_within_between(random)}} 
    else if(perm_within_between==F) {for (perm in 1:nperm)  {permseq[,perm]=sample.int(NROW(model))}}
    
    #activate parallel processing
    unregister_dopar = function() {
      env = foreach:::.foreachGlobals
      rm(list=ls(name=env), pos=env)
    }
    unregister_dopar()
    
    cl=parallel::makeCluster(nthread)
    doParallel::registerDoParallel(cl)
    `%dopar%` = foreach::`%dopar%`
    
    #progress bar
    doSNOW::registerDoSNOW(cl)
    pb=txtProgressBar(max = nperm, style = 3)
    progress=function(n) setTxtProgressBar(pb, n)
    opts=list(progress = progress)
    
    #fitting permuted model and extracting max-TFCE values in parallel streams
    start=Sys.time()
    TFCE.max=foreach::foreach(perm=1:nperm, .combine="rbind",.export=c("TFCE","edgelist","getClusters"), .options.snow = opts)  %dopar%
      {
        brainstat.stats.terms=reticulate::import("brainstat.stats.terms")
        brainstat.stats.SLM=reticulate::import("brainstat.stats.SLM")

      ##commented out alternative method of permutation
        #model.permuted=model
        #model.permuted[,colno]=model.permuted[permseq[,perm],colno] ##permute only the contrast
        #terms=brainstat.stats.terms$MixedEffect(ran = random,fix = model.permuted,"_check_categorical" = F)
        #model.fit=brainstat.stats.SLM$SLM(model = terms,
        #                                   contrast=contrast,
        #                                   surf = template, 
        #                                   mask=mask,
        #                                   correction="None",
        #                                   cluster_threshold=1)
        #model.fit$fit(CT_data)

        terms=brainstat.stats.terms$MixedEffect(ran = random,fix = model,"_check_categorical" = F)
        model.fit=brainstat.stats.SLM$SLM(model = terms,
                                          contrast=contrast,
                                          surf = template, 
                                          mask=mask,
                                          correction="None",
                                          cluster_threshold=1)
        model.fit$fit(CT_data[permseq[,perm],])
        
        return(max(abs(suppressWarnings(TFCE(data = as.numeric(model.fit$t),tail = 2)))))
      }
    end=Sys.time()
    cat(paste("\nCompleted in ",round(difftime(end, start, units='mins'),1)," minutes \n",sep=""))
    
  ##saving list objects
  returnobj=list(tmap.orig,TFCE.orig, TFCE.max,tail)
  names(returnobj)=c("t_stat","TFCE.orig","TFCE.max","tail")
  
  return(returnobj)
}

############################################################################################################################
############################################################################################################################
##Example

#model=TFCE.vertex_analysis.mixed(model = dat_beh[,c(2:4)],contrast = dat_beh$AGE_AT_SCAN,random = dat_beh$SUB_ID,CT_data = dat_CT,nperm = 100,tail = 2, nthread = 10, smooth_FWHM = 0)

#results=TFCE.threshold(model)
#results$cluster_level_results
