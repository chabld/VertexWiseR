## FUNCTIONS FOR VERTEX-WISE TFCE ANALYSIS
## FOR USE IN THE COGNITIVE AND BRAIN HEALTH LABORATORY

############################################################################################################################
############################################################################################################################
loadRData <- function(fileName){
    #loads an RData file, and returns it
    load(fileName)
    get(ls()[ls() != "fileName"])
  }
##Main function

  
TFCE.vertex_analysis=function(model,contrast, CT_data, nperm=100, tail=2, nthread=10, smooth_FWHM)
{
  if(class(contrast)=="integer") {contrast=as.numeric(contrast)}
  ##load other vertex-wise functions
  source("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/otherfunc.r?raw=TRUE")
  
  ##checks
    #check if required packages are installed
    packages=c("foreach","doParallel","parallel","doSNOW")
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
      cat(paste("model contains",length(idxF),"subjects with incomplete data. Subjects with incomplete data will be excluded in the current analysis\n"))
      model=model[-idxF,]
      contrast=contrast[-idxF]
      CT_data=CT_data[-idxF,]
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
          if(identical(suppressWarnings(as.numeric(contrast)),suppressWarnings(as.numeric(model[,colno]))))  {break}
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
    n_vert=ncol(CT_data)
      if(n_vert==81924) 
    {
      edgelist=loadRData(url("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/data/edgelistfs6.rdata?raw=TRUE"))
      assign("edgelist", edgelist, envir = globalenv())
    } 
    else if(n_vert==14524) 
    {
      edgelist=loadRData(url("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/data/edgelistHIP.rdata?raw=TRUE"))
      assign("edgelist", edgelist, envir = globalenv())
    }
    else if (n_vert==14524)  
    {
      edgelist=loadRData(url("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/data/edgelistHIP.rdata?raw=TRUE"))
      assign("edgelist", edgelist, envir = globalenv())
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
  n_vert=NCOL(CT_data)
  if(missing("smooth_FWHM"))
  {
    if(n_vert==20484) 
    {
      cat("CT_data will be smoothed using the default 10mm FWHM kernel for fsaverage5 images\n")
      CT_data=smooth(CT_data, FWHM=10)
    }
    else if(n_vert==81924) 
    {
      cat("CT_data will be smoothed using the default 5mm FWHM kernel for fsaverage6 images\n")
      CT_data=smooth(CT_data, FWHM=5)
    }
    else if(n_vert==14524) 
    {
      cat("CT_data will be smoothed using the default 5mm FWHM kernel for hippocampal maps\n")
      CT_data=smooth(CT_data, FWHM=5)
    }
  } else if(smooth_FWHM>0) 
  {
    cat(paste("CT_data will be smoothed using a ", smooth_FWHM,"mm FWHM kernel\n", sep=""))
    CT_data=smooth(CT_data, FWHM=smooth_FWHM)
  }
  
  ##unpermuted model
  model=data.matrix(model)
  mod=.lm.fit(y=CT_data,x=data.matrix(cbind(1,model)))
  
  #extract tstat and calculate tfce image
  start=Sys.time()
  cat("Estimating unpermuted TFCE image...")
  
  tmap.orig=extract.t(mod,colno+1)
  TFCE.orig=suppressWarnings(TFCE.multicore(data = tmap.orig,tail = tail,nthread=nthread))
  remove(mod)
  
  end=Sys.time()
  cat(paste("Completed in",round(difftime(end,start, units="secs"),1),"secs\nEstimating permuted TFCE images...\n",sep=" "))
  
  ##permuted models
  #generating permutation sequences  
  set.seed(123)
  permseq=matrix(NA, nrow=NROW(model), ncol=nperm)
  for (perm in 1:nperm)  {permseq[,perm]=sample.int(NROW(model))}
  
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
  
  ##fitting permuted regression model and extracting t-stats in parallel streams
  start=Sys.time()
  
  TFCE.max=foreach::foreach(perm=1:nperm, .combine="rbind",.export=c("TFCE","extract.t","getClusters","edgelist"), .options.snow = opts)  %dopar%
    {
      ##commented out alternative method of permutation— permuting only the contrast variable
      #model.permuted=model
      #model.permuted[,colno]=model.permuted[permseq[,perm],colno] ##permute only the contrast
      #mod.permuted=lm(CT_data~data.matrix(model.permuted))
      
      mod.permuted=.lm.fit(y=CT_data[permseq[,perm],],x=data.matrix(cbind(1,model)))
      tmap=extract.t(mod.permuted,colno+1)
      
      remove(mod.permuted,model.permuted)
      return(max(abs(suppressWarnings(TFCE(data = tmap,tail = tail)))))
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
##TFCE single core— for estimating permuted TFCE statistics
##adapted from nilearn python library: https://github.com/nilearn/nilearn/blob/main/nilearn/mass_univariate/_utils.py#L7C8-L7C8
TFCE=function(data,tail=tail)
{
  #selecting tail type
  if (tail==2) 
  {
    signs = c(-1, 1)
    max_score = max(abs(data),na.rm = T)
  } else if(tail==1)
  {
    signs = 1
    max_score = max(data,na.rm = T)
  } else if(tail==-1)
  {
    signs = -1
    max_score = max(-data,na.rm = T)
  }
  
  #define TFCE parameters
  step=max_score / 100 #calculating number of steps for TFCE estimation
  score_threshs = seq(step, max_score, by = step) #Set based on determined step size
  
  #loop across different signs(i.e., for two tailed test)
  for (sign.idx in 1:length(signs)) 
  {
    temp_data = data * signs[sign.idx]
    tfce=rep(0,length(temp_data))
    
    #loop across different score_threshs values for TFCE estimation
    for(thresh.no in 1:length(score_threshs))
    {
      temp_data[temp_data < score_threshs[thresh.no]] = 0
      if(length(which(temp_data>0))>1) #if less than 2 vertices, skip the following steps
      {
        clust.dat=getClusters(temp_data)
        if (clust.dat[[2]][1]!="noclusters") #if no clusters, skip the following steps
        {
          non_zero_inds = which(clust.dat[[1]] >0)
          labeled_non_zero = clust.dat[[1]][non_zero_inds]
          cluster_tfces = signs[sign.idx] * clust.dat[[2]] * (score_threshs[thresh.no] ^ 2) #using the E=1 , H=2 paramters for 2D (vertex-wise data)
          tfce_step_values = rep(0, length(clust.dat[[1]]))
          tfce_step_values[non_zero_inds] = cluster_tfces[labeled_non_zero]
          tfce[non_zero_inds]=tfce_step_values[non_zero_inds]+tfce[non_zero_inds] #cumulatively add up all TFCE values at each vertex
          remove(clust.dat,non_zero_inds,cluster_tfces,tfce_step_values, labeled_non_zero)
        }
      }
    }
    #combine results from positive and negative tails if necessary 
    if(sign.idx==1){tfce_step_values.all=tfce}
    else if (sign.idx==2){tfce_step_values.all=tfce_step_values.all+tfce}
    remove(tfce)
  }
  return(tfce_step_values.all)
}
############################################################################################################################
############################################################################################################################
##TFCE multicore— for estimating unpermuted TFCE statistics
##adapted from nilearn python library: https://github.com/nilearn/nilearn/blob/main/nilearn/mass_univariate/_utils.py#L7C8-L7C8
TFCE.multicore=function(data,tail=tail,nthread)
{
  #selecting tail type
  if (tail==2) 
  {
    signs = c(-1, 1)
    max_score = max(abs(data),na.rm = T)
  } else if(tail==1)
  {
    signs = 1
    max_score = max(data,na.rm = T)
  } else if(tail==-1)
  {
    signs = -1
    max_score = max(-data,na.rm = T)
  }
  
  #define TFCE parameters
  step=max_score / 100 #calculating number of steps for TFCE estimation
  score_threshs = seq(step, max_score, by = step) #Set based on determined step size
  
  #loop across different signs(i.e., for two tailed test)
  for (sign.idx in 1:length(signs)) 
  {
    temp_data = data * signs[sign.idx]
    tfce=rep(0,length(temp_data))
    
    #activate parallel processing
    unregister_dopar = function() {
      env = foreach:::.foreachGlobals
      rm(list=ls(name=env), pos=env)
    }
    unregister_dopar()
    
    cl=parallel::makeCluster(nthread)
    doParallel::registerDoParallel(cl)
    `%dopar%` = foreach::`%dopar%`
    
    #parallel loop across different score_threshs values for TFCE estimation
    tfce=foreach::foreach(thresh.no=1:length(score_threshs), .combine="rbind", .export=c("getClusters","edgelist"))  %dopar%
      {
        temp_data[temp_data < score_threshs[thresh.no]] = 0
        if(length(which(temp_data>0))>1) #if less than 2 vertices, skip the following steps
        {
          clust.dat=getClusters(temp_data)
          if (clust.dat[[2]][1]!="noclusters") #if no clusters, skip the following steps
          {
            non_zero_inds = which(clust.dat[[1]] >0)
            labeled_non_zero = clust.dat[[1]][non_zero_inds]
            cluster_tfces = signs[sign.idx] * clust.dat[[2]] * (score_threshs[thresh.no] ^ 2)
            tfce_step_values = rep(0, length(clust.dat[[1]]))
            tfce_step_values[non_zero_inds] = cluster_tfces[labeled_non_zero]
            return(tfce_step_values)
          }
        }
      }
    #suppressWarnings(closeAllConnections())
    
    #combine results from positive and negative tails if necessary 
    if(length(tfce)>length(temp_data))
    {
      tfce=colSums(tfce)
    } else if(length(tfce)==0)
    {
      tfce=0
    }
    if(sign.idx==1)
    {
      tfce_step_values.all=tfce
    } else if (sign.idx==2)
    {
      tfce_step_values.all=tfce_step_values.all+tfce
    }
  }
  parallel::stopCluster(cl)
  return(tfce_step_values.all)
}
############################################################################################################################
############################################################################################################################
TFCE.threshold=function(TFCE.output, p=0.05, atlas=1, k=20)
{
  nperm=length(TFCE.output$TFCE.max)
  
  #check if number of permutations is adequate
  if(nperm<1/p)  {warning(paste("Not enough permutations were carried out to estimate the p<",p," threshold precisely\nConsider setting an nperm to at least ",ceiling(1/p),sep=""))}
  
  #check which template is used and load appropriate tempalte files
  n_vert=length(TFCE.output$t_stat)
  if(n_vert==20484) 
  {    
    load(file = url("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/data/edgelistfs5.rdata?raw=TRUE"))
    load(file = url("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/data/ROImap_fs5.rdata?raw=TRUE"))
    load(file = url("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/data/MNImap_fs5.rdata?raw=TRUE"))

    edgelist=edgelistfs5    
    ROImap=ROImap_fs5
    MNImap=MNImap_fs5
    remove(edgelistfs5,ROImap_fs5,MNImap_fs5)
  }
  else if (n_vert==81924) 
  {
    load(file = url("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/data/edgelistfs6.rdata?raw=TRUE"))
    load(file = url("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/data/ROImap_fs6.rdata?raw=TRUE"))
    load(file = url("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/data/MNImap_fs6.rdata?raw=TRUE"))

    edgelist=edgelistfs6   
    ROImap=ROImap_fs6
    MNImap=MNImap_fs6
    remove(edgelistfs6,ROImap_fs6,MNImap_fs6)
  } 
  else if (n_vert==14524) 
  {
    load(file = url("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/data/ROImap_hip.rdata?raw=TRUE"))
    ROImap=list(data.matrix(ROImap_HIP[[1]]),ROImap_HIP[[2]])
    load(file = url("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/data/edgelistHIP.rdata?raw=TRUE"))
    load(file = url("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/data/MNImap_hip.rdata?raw=TRUE"))
    
    edgelist=edgelistHIP  
    MNImap=MNImap_hip
    remove(edgelistHIP,ROImap_HIP,MNImap_hip)
  } 
  ##generating p map
  tfce.p=rep(NA,n_vert)
  TFCE.output$t_stat[is.na(TFCE.output$t_stat)]=0
  for (vert in 1:n_vert)  {tfce.p[vert]=length(which(TFCE.output$TFCE.max>abs(TFCE.output$TFCE.orig[vert])))/nperm}
  
  ##generating thresholded t-stat map
  TFCE.output$t_stat[is.na(TFCE.output$t_stat)]=0
  
  t_stat.thresholdedP=TFCE.output$t_stat
  t_stat.thresholdedP[tfce.p>p]=0
  
  ##Cluster level results
  ##positive cluster
  if(TFCE.output$tail==1 |TFCE.output$tail==2)
  {
    #applying p thresholding
    pos.t_stat.thresholdedP=t_stat.thresholdedP
    pos.t_stat.thresholdedP[pos.t_stat.thresholdedP<0]=0
    
    if(length(which(pos.t_stat.thresholdedP!=0))>1) #skip if no clusters detected
    {
      pos.clusters0=getClusters(pos.t_stat.thresholdedP) ## 1st getCluster() to identify all clusters with no. vertices > 1
      #applying k thresholding
      pos.clustID.remove=which(pos.clusters0[[2]]<k)
      pos.clusters0[[1]][which(!is.na(match(pos.clusters0[[1]],pos.clustID.remove)))]=NA
      
      #generating mask
      pos.clusters=getClusters(pos.clusters0[[1]]) ## 2nd getCluster() to identify all clusters from the k-thresholded clustermap
      pos.clusters[[1]][is.na(pos.clusters[[1]])]=0
      pos.mask=rep(0,n_vert)
      ROImap[[1]][, atlas]
      ROImap=list(data.matrix(ROImap[[1]]),ROImap[[2]])

      #results table
      pos.clustermap=rep(NA,n_vert)
      
      if(pos.clusters[[2]][1]!="noclusters")
      {
        pos.mask[pos.clusters[[1]]>0]=1
        pos.clust.results=data.frame(matrix(NA,nrow=length(pos.clusters[[2]]), ncol=8))
        colnames(pos.clust.results)=c("clusid","nverts","P","X","Y","Z","tstat","region")
        clust.idx=1
        
        for(clust.no in order(pos.clusters[[2]],decreasing = T))
        {
          clust.vert.idx=which(pos.clusters[[1]]==clust.no)
          pos.clustermap[clust.vert.idx]=clust.idx
          
          pos.clust.results[clust.idx,1]=clust.idx
          pos.clust.results[clust.idx,2]=length(clust.vert.idx)
          max.vert.idx=clust.vert.idx[which(abs(TFCE.output$t_stat[clust.vert.idx])==max(abs(TFCE.output$t_stat[clust.vert.idx]),na.rm = T))[1]]
          pos.clust.results[clust.idx,3]=round(tfce.p[max.vert.idx],3)
          if(pos.clust.results[clust.idx,3]==0) {pos.clust.results[clust.idx,3]=paste("<",1/nperm,sep="")}
          pos.clust.results[clust.idx,c(4,5,6)]=round(MNImap[,max.vert.idx],1)
          pos.clust.results[clust.idx,7]=round(abs(TFCE.output$t_stat[max.vert.idx]),2)
          
          atlas.idx=ROImap[[1]][,atlas][max.vert.idx]
          if(atlas.idx>0) {pos.clust.results[clust.idx,8]=ROImap[[2]][,atlas][atlas.idx] } ##to deal with desikan atlas missing vertex mappings
          else {pos.clust.results[clust.idx,8]="unknown (use another atlas)"}
          
          clust.idx=clust.idx+1
        }
      } else if(pos.clusters[[2]][1]=="noclusters")
      { 
        pos.clust.results="No significant clusters"
        pos.clustermap="No significant clusters"
      }
    } else ## 2nd getClusters()
    {
      pos.clust.results="No significant clusters"
      pos.clustermap="No significant clusters"
      pos.mask=rep(0,n_vert)
    }
  } else if(TFCE.output$tail==-1)
  {
    pos.clust.results="Positive contrast not analyzed, only negative one-tailed TFCE statistics were estimated)"
    pos.clustermap="No significant clusters"
    pos.mask=rep(0,n_vert)
  } 
  
  ##negative cluster
  if(TFCE.output$tail==-1 |TFCE.output$tail==2)
  {
    #applying p thresholding
    neg.t_stat.thresholdedP=t_stat.thresholdedP
    neg.t_stat.thresholdedP[neg.t_stat.thresholdedP>0]=0
    
    if(length(which(neg.t_stat.thresholdedP!=0))>1) #skip if no clusters detected
    {
      neg.clusters0=getClusters(neg.t_stat.thresholdedP) ## 1st getCluster() to identify all clusters with no. vertices > 1
      #applying k thresholding
      neg.clustID.remove=which(neg.clusters0[[2]]<k)
      neg.clusters0[[1]][which(!is.na(match(neg.clusters0[[1]],neg.clustID.remove)))]=NA
      
      #generating mask
      neg.clusters=getClusters(neg.clusters0[[1]]) ## 2nd getCluster() to identify all clusters from the k-thresholded clustermap
      neg.clusters[[1]][is.na(neg.clusters[[1]])]=0
      neg.mask=rep(0,n_vert)
      
      #results table
      neg.clustermap=rep(NA,n_vert)
      
      if(neg.clusters[[2]][1]!="noclusters")
      {
        neg.mask[neg.clusters[[1]]>0]=1
        neg.clust.results=data.frame(matrix(NA,nrow=length(neg.clusters[[2]]), ncol=8))
        colnames(neg.clust.results)=c("clusid","nverts","P","X","Y","Z","tstat","region")
        clust.idx=1
        
        for(clust.no in order(neg.clusters[[2]],decreasing = T))
        {
          clust.vert.idx=which(neg.clusters[[1]]==clust.no)
          neg.clustermap[clust.vert.idx]=clust.idx
          
          neg.clust.results[clust.idx,1]=clust.idx
          neg.clust.results[clust.idx,2]=length(clust.vert.idx)
          max.vert.idx=clust.vert.idx[which(abs(TFCE.output$t_stat[clust.vert.idx])==max(abs(TFCE.output$t_stat[clust.vert.idx]),na.rm = T))[1]]
          neg.clust.results[clust.idx,3]=round(tfce.p[max.vert.idx],3)
          if(neg.clust.results[clust.idx,3]==0) {neg.clust.results[clust.idx,3]=paste("<",1/nperm,sep="")}
          neg.clust.results[clust.idx,c(4,5,6)]=round(MNImap[,max.vert.idx],1)
          neg.clust.results[clust.idx,7]=round(abs(TFCE.output$t_stat[max.vert.idx]),2)
          
          atlas.idx=ROImap[[1]][,atlas][max.vert.idx]
          if(atlas.idx>0) {neg.clust.results[clust.idx,8]=ROImap[[2]][,atlas][atlas.idx] } ##to deal with desikan atlas missing vertex mappings
          else {neg.clust.results[clust.idx,8]="unknown (use another atlas)"}
          
          clust.idx=clust.idx+1
        }
      } else if(neg.clusters[[2]][1]=="noclusters")
      { 
        neg.clust.results="No significant clusters"
        neg.clustermap="No significant clusters"
      } 
    } else  ## 2nd getClusters()
    {
      neg.clust.results="No significant clusters"
      neg.clustermap="No significant clusters"
      neg.mask=rep(0,n_vert)
    }
  }
  else if(TFCE.output$tail==1)
  {
    neg.clust.results="Negative contrast not analyzed, only negative one-tailed TFCE statistics were estimated"
    neg.clustermap="No significant clusters"
    neg.mask=rep(0,n_vert)
  } 
  ##saving list objects
  cluster_level_results=list(pos.clust.results,neg.clust.results)
  names(cluster_level_results)=c("Positive contrast", "Negative contrasts")
  
  t_stat.thresholdedPK=TFCE.output$t_stat*(pos.mask+neg.mask)
  
  returnobj=list(cluster_level_results, t_stat.thresholdedPK, pos.clustermap, neg.clustermap)  
  names(returnobj)=c("cluster_level_results","thresholded_tstat_map","pos_clusterIDmap","neg_clusterIDmap")
  returnobj$cluster_level_results
  
  return(returnobj)
}  
############################################################################################################################
############################################################################################################################
##example

#pos=TFCE.vertex_analysis(model =all_pred, contrast = dat_beh$age, CT_data = dat_CT, tail=1, nperm=100, nthread = 10)
#neg=TFCE.vertex_analysis(model =all_pred, contrast = dat_beh$age, CT_data = dat_CT, tail=-1 ,nperm=100, nthread = 10)
#two=TFCE.vertex_analysis(model =all_pred, contrast = dat_beh$age, CT_data = dat_CT, tail=2 ,nperm=100, nthread = 10)

#pos.results=TFCE.threshold(pos)
#pos.results$cluster_level_results
#plotCT(pos.results$thresholded_tstat_map, filename="pos.png")

#neg.results=TFCE.threshold(neg)
#neg.results$cluster_level_results
#plotCT(neg.results$thresholded_tstat_map, filename="neg.png")

#two.results=TFCE.threshold(two)
#two.results$cluster_level_results
#plotCT(two.results$thresholded_tstat_map, filename="two.png")
