## FUNCTION FOR MIXED-EFFECT VERTEX-WISE TFCE ANALYSIS
## FOR USE IN THE COGNITIVE AND BRAIN HEALTH LABORATORY

############################################################################################################################
############################################################################################################################
#' @title Vertex-wise analysis with TFCE (mixed effect)
#'
#' @description Fits a model with the whole-brain and hippocampal surface data in template space. The data is smoothed and fit to a linear model with mixed effects, and returns a brain-wide or hippocampal t-value maps, as well as cluster-corrected maps with threshold-free cluster enhancement.
#' 
#' @details The TFCE method for estimating unpermuted TFCE statistics is adapted from the \href{https://github.com/nilearn/nilearn/blob/main/nilearn/mass_univariate/_utils.py#L7C8-L7C8}{nilearn python library}. 
#' 
#' @param model A data.frame object containing the variables to include in the model at each column, and rows of values assigned to each participant.
#' @param contrast An object containing the values of the independent variable of interest for which to fit a contrast
#' @param surf_data A matrix object containing the surface data, see SURFvextract() output format. 
#' @param random An object containing the values of the random variable 
#' @param nperm A numeric integer object stating the number of permutations wanted for the cluster-correction (default = 100)
#' @param tail A numeric integer object stating whether to test a one-sided (1,-1) or two-sided (2) model
#' @param nthread Maximum number of cpu cores to allocate 
#' @param smooth_FWHM A numeric vector object containing the desired smoothing width in mm 
#' @param perm_within_between A logical object stating whether to implement a permutation function for random subject effects 
#'
#'
#' @returns A list object containing  the threshold t-test, the TFCE cluster output, and permuted TFCE cluster maps. 
#' @examples
#'model=TFCE.vertex_analysis.mixed(model = dat_beh[,c(2:4)],contrast = dat_beh$AGE_AT_SCAN,random = dat_beh$SUB_ID,surf_data = dat_CT,nperm = 100,tail = 2, nthread = 10, smooth_FWHM = 0)
#'
#'results=TFCE.threshold(model)
#'results$cluster_level_results
#'
#' @importFrom reticulate import r_to_py
#' @importFrom foreach foreach 
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom doSNOW registerDoSNOW
#' @export


##Main function

TFCE.vertex_analysis.mixed=function(model,contrast, surf_data, random, nperm=100, tail=2, nthread=10, smooth_FWHM, perm_within_between=F)
{
  
  #If the contrast/model is a tibble (e.g., taken from a read_csv output)
  #converts the columns to regular data.frame column types
  if ('tbl_df' %in% class(contrast)) {
    if (class(contrast[[1]])=="character") {contrast = contrast[[1]]}
    if(class(contrast[[1]])=="integer") {contrast = as.numeric(contrast[[1]])}
  } 
  if ('tbl_df' %in% class(model)) {
    model=as.data.frame(model)
    if (NCOL(model)==1) {model = model[[1]]
    } else { for (c in 1:NCOL(model)) { 
      if(class(model[,c])=="double") {model[,c] = as.numeric(model[,c])}
    }  }
  }
  
  if(class(contrast)=="integer") {contrast=as.numeric(contrast)}
  
  ##load other TFCE and vertex-wise functions
  source("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/R/vertTFCE.r?raw=TRUE")
  source("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/R/otherfunc.r?raw=TRUE")
  
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
    #check if nrow is consistent for model and surf_data
    if(NROW(surf_data)!=NROW(model))  {stop(paste("The number of rows for surf_data (",NROW(surf_data),") and model (",NROW(model),") are not the same",sep=""))}
    if(length(random)!=NROW(model))  {stop(paste("The number of rows for random (",NROW(surf_data),") and model (",NROW(model),") are not the same",sep=""))}
  
    #incomplete data check
    idxF=which(complete.cases(model)==F)
    if(length(idxF)>0)
    {
      cat(paste("model contains",length(idxF),"subjects with incomplete data. Subjects with incomplete data will be excluded in the current analysis\n"))
      model=model[-idxF,]
      contrast=contrast[-idxF]
      surf_data=surf_data[-idxF,]
    }
    
    #check contrast
    if(NCOL(model)>1)
    {
      for(colno in 1:(NCOL(model)+1))
      {
        if(colno==(NCOL(model)+1))  {warning("contrast is not contained within model")}
        
        if (class(contrast)=="character") 
        {
          if(identical(data.matrix(contrast),data.matrix(model)[,colno]))  {break} 
        } else 
        {
          if(identical(as.numeric(contrast),as.numeric(model[,colno])))  {break}
        }
      }
    }  else
    {
      if (class(contrast)=="character") 
      {
        if(identical(contrast,model))  {colno=1} 
        else  {warning("contrast is not contained within model")}
      } else
      {
        if(identical(as.numeric(contrast),as.numeric(model)))  {colno=1}
        else  {warning("contrast is not contained within model")}
      }
    }
    
    #check categorical and recode variable
    if(NCOL(model)>1)
    {
      for (column in 1:NCOL(model))
      {
        if(class(model[,column])=="character") 
        {
          if(length(unique(model[,column]))==2)
          {
            cat(paste("The binary variable '",colnames(model)[column],"' will be recoded with ",unique(model[,column])[1],"=0 and ",unique(model[,column])[2],"=1 for the analysis\n",sep=""))
            
            recode=rep(0,NROW(model))
            recode[model[,column]==unique(model[,column])[2]]=1
            model[,column]=recode
            contrast=model[,colno]
          } else if(length(unique(model[,column]))>2)    {stop(paste("The categorical variable '",colnames(model)[column],"' contains more than 2 levels, please code it into binarized dummy variables",sep=""))}
        }      
      }
    } else
    {
      if (class(model)=="character") 
      {
        if(length(unique(model))==2)
        {
          cat(paste("The binary variable '",colnames(model),"' will be recoded such that ",unique(model)[1],"=0 and ",unique(model)[2],"=1 for the analysis\n",sep=""))
          
          recode=rep(0,NROW(model))
          recode[model==unique(model)[2]]=1
          model=recode
          model=model
        } else if(length(unique(model))>2)    {stop(paste("The categorical variable '",colnames(model),"' contains more than 2 levels, please code it into binarized dummy variables",sep=""))}
      }      
    }
    
    #check length of CT data and load the appropriate fsaverage files
    n_vert=ncol(surf_data)
    if(n_vert==20484)
    {
      template="fsaverage5"
      load(file = url("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/data/edgelistfs5.rdata?raw=TRUE"),envir = globalenv())
    }
    else if (n_vert==81924)
    {
      template="fsaverage6"
      load(file = url("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/data/edgelistfs6.rdata?raw=TRUE"),envir = globalenv())
    }
    else if (n_vert==14524)
    {
      if(file.exists("hip_template.fs")==F)
      {
        cat("\nhip_template.fs is not detected in the current working directory. The hippocampus surface template will be downloaded\n")
        download.file(url="https://raw.githubusercontent.com/CogBrainHealthLab/VertexWiseR/main/data/hip_template.fs",destfile ="hip_template.fs",mode = "wb")
      } 
      brainspace.mesh.mesh_io=reticulate::import("brainspace.mesh.mesh_io")
      template=brainspace.mesh.mesh_io$read_surface("hip_template.fs")
  
      load(file = url("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/data/edgelistHIP.rdata?raw=TRUE"),envir = globalenv())
    }
    else {stop("data vector should only contain 20484 (fsaverage5), 81924 (fsaverage6) or 14524 (hippocampal vertices) columns")}
    
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
    n_vert=NCOL(surf_data)
    if(missing("smooth_FWHM"))
    {
      if(n_vert==20484) 
      {
        cat("surf_data will be smoothed using the default 10mm FWHM kernel for fsaverage5 images\n")
        surf_data=smooth(surf_data, FWHM=10)
      }
      else if(n_vert==81924) 
      {
        cat("surf_data will be smoothed using the default 5mm FWHM kernel for fsaverage6 images\n")
        surf_data=smooth(surf_data, FWHM=5)
      }
      else if(n_vert==14524) 
      {
        cat("surf_data will be smoothed using the default 5mm FWHM kernel for hippocampal maps\n")
        surf_data=smooth(surf_data, FWHM=5)
      }
    } else if(smooth_FWHM>0) 
    {
      cat(paste("surf_data will be smoothed using a ", smooth_FWHM,"mm FWHM kernel\n", sep=""))
      surf_data=smooth(surf_data, FWHM=smooth_FWHM)
    }
    surf_data[is.na(surf_data)]=0
  
  ##unpermuted model
    #preparing mask for model
    mask=array(rep(T,NCOL(surf_data)))
    maskNA=which(colSums(surf_data != 0) == 0)
    mask[which(colSums(surf_data != 0) == 0)]=F
    
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
    model.fit$fit(surf_data)
    
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

      ##commented out alternative method of permutation— permuting only the contrast variable
        #model.permuted=model
        #model.permuted[,colno]=model.permuted[permseq[,perm],colno] ##permute only the contrast
        #terms=brainstat.stats.terms$MixedEffect(ran = random,fix = model.permuted,"_check_categorical" = F)
        #model.fit=brainstat.stats.SLM$SLM(model = terms,
        #                                   contrast=contrast,
        #                                   surf = template, 
        #                                   mask=mask,
        #                                   correction="None",
        #                                   cluster_threshold=1)
        #model.fit$fit(surf_data)

        terms=brainstat.stats.terms$MixedEffect(ran = random,fix = model,"_check_categorical" = F)
        model.fit=brainstat.stats.SLM$SLM(model = terms,
                                          contrast=contrast,
                                          surf = template, 
                                          mask=mask,
                                          correction="None",
                                          cluster_threshold=1)
        model.fit$fit(surf_data[permseq[,perm],])
        
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

#model=TFCE.vertex_analysis.mixed(model = dat_beh[,c(2:4)],contrast = dat_beh$AGE_AT_SCAN,random = dat_beh$SUB_ID,surf_data = dat_CT,nperm = 100,tail = 2, nthread = 10, smooth_FWHM = 0)

#results=TFCE.threshold(model)
#results$cluster_level_results
