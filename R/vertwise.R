## FUNCTION FOR VERTEX-WISE ANALYSIS WITH MIXED EFFECTS
## ADAPTED FROM brainstat python library (https://brainstat.readthedocs.io/en/master/_modules/brainstat/stats/SLM.html#SLM)
## FOR USE IN THE COGNITIVE AND BRAIN HEALTH LABORATORY

############################################################################################################################
############################################################################################################################
#' @title Vertex-wise analysis
#'
#' @description Fits a model with the whole-brain and hippocampal surface data in template space. The data is smoothed and fit to a linear model with fixed or mixed effects, and returns a brain-wide or hippocampal t-value maps, as well as random field theory cluster maps. 
#'
#' @details The function imports and adapts the \href{https://brainstat.readthedocs.io/en/master/_modules/brainstat/stats/SLM.html#SLM)}{brainstat python library}. 
#'
#' @param model A data.frame object containing the variables to include in the model at each column, and rows of values assigned to each participant.
#' @param contrast An object containing the values of the independent variable of interest for which to fit a contrast
#' @param random An object containing the values of the random variable (optional)
#' @param surf_data A matrix object containing the surface data, see SURFvextract() output format. 
#' @param p A numeric object stating the p-value threshold for the linear model and cluster-correction
#' @param atlas A numeric integer object corresponding to the atlas of interest. 1=Desikan, 2=Schaefer-100, 3=Schaefer-200, 4=Glasser-360, 5=Destrieux-148.
#' @param smooth_FWHM A numeric vector object containing the desired smoothing width in mm 
#'
#'
#' @return A list object containing the results at cluster level, the threshold t-test map, positive and negative results maps, positive and negative  cluster maps. 
#' @examples
#' model=vertex_analysis(model = dat_beh[,2:3], contrast = dat_beh$Age, random = dat_beh$SUB_ID, surf_data = dat_CT,p = 0.01, atlas=1)
#' model$cluster_level_results
#' @importFrom reticulate import r_to_py
#' @export

##vertex wise analysis with mixed effects
vertex_analysis=function(model,contrast, random, surf_data, p=0.05, atlas=1, smooth_FWHM)  ## atlas: 1=Desikan, 2=Schaefer-100, 3=Schaefer-200, 4=Glasser-360, 5=Destrieux-148; ignored for hippocampal surfaces
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
  
  ##load other vertex-wise functions (not needed when package is in library)
  #source("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/R/otherfunc.r?raw=TRUE")
  
  ##checks
    #check contrast
    if(NCOL(model)>1)
    {
      for(colno in 1:(NCOL(model)+1))
      {
        if(colno==(NCOL(model)+1))  {warning("contrast is not contained within model")}
        
        if(class(contrast)=="character") 
        {
          if(identical(contrast,data.matrix(model)[,colno]))  {break} 
        } else 
        {
          if(identical(suppressWarnings(as.numeric(contrast)),suppressWarnings(as.numeric(model[,colno]))))  {break}
        }
      }
    }  else
    {
      if(class(contrast)=="character") 
      {
        if(identical(contrast,model))  {colno=1} 
        else  {warning("contrast is not contained within model")}
      } else
      {
        if(identical(as.numeric(contrast),as.numeric(model)))  {colno=1}
        else  {warning("contrast is not contained within model")}
      }
    }
    
    #check if nrow is consistent for model and surf_data
    if(NROW(surf_data)!=NROW(model))  {stop(paste("The number of rows for surf_data (",NROW(surf_data),") and model (",NROW(model),") are not the same",sep=""))}
    
    #incomplete data check
    idxF=which(complete.cases(model)==F)
    if(length(idxF)>0)
    {
      cat(paste("model contains",length(idxF),"subjects with incomplete data. Subjects with incomplete data will be excluded in the current analysis\n"))
      model=model[-idxF,]
      contrast=contrast[-idxF]
      surf_data=surf_data[-idxF,]
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
      if(class(model)=="character") 
      {
        if(length(unique(model))==2)
        {
          cat(paste("The binary variable '",colnames(model),"' will be recoded such that ",unique(model)[1],"=0 and ",unique(model)[2],"=1 for the analysis\n",sep=""))
          
          recode=rep(0,NROW(model))
          recode[model==unique(model)[2]]=1
          model=recode
          contrast=model
        } else if(length(unique(model))>2)    {stop(paste("The categorical variable '",colnames(model),"' contains more than 2 levels, please code it into binarized dummy variables",sep=""))}
      }      
    }
    
    #check length of CT data and load the appropriate fsaverage files
    n_vert=ncol(surf_data)
    if(n_vert==20484)
    {
      template="fsaverage5"
      load(file = url("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/data/ROImap_fs5.rdata?raw=TRUE"))
    }
    else if (n_vert==81924)
    {
      template="fsaverage6"
      load(file = url("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/data/ROImap_fs6.rdata?raw=TRUE"))
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
      load(file = url("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/data/ROImap_hip.rdata?raw=TRUE"))
    }
    else {stop("data vector should only contain 20484 (fsaverage5), 81924 (fsaverage6) or 14524 (hippocampal vertices) columns")}
  
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
        cat("surf_data will be smoothed using the default 5mm FWHM kernel for fsaverage6 images")
        surf_data=smooth(surf_data, FWHM=5)
      }
      else if(n_vert==14524) 
      {
        cat("surf_data will be smoothed using the default 5mm FWHM kernel for hippocampal maps\n")
        surf_data=smooth(surf_data, FWHM=5)
      }
    } else if(smooth_FWHM>0) 
    {
      cat(paste("surf_data will be smoothed using a ",smooth_FWHM,"mm FWHM kernel", sep=""))
      surf_data=smooth(surf_data, FWHM=smooth_FWHM)
    }
    surf_data[is.na(surf_data)]=0
  
  ##import python libaries
  brainstat.stats.terms=reticulate::import("brainstat.stats.terms")
  brainstat.stats.SLM=reticulate::import("brainstat.stats.SLM")
  brainstat.datasets=reticulate::import("brainstat.datasets")  
  
  ##fitting model
  #preparing mask for model
  mask=array(rep(T,NCOL(surf_data)))
  maskNA=which(colSums(surf_data != 0) == 0)
  mask[which(colSums(surf_data != 0) == 0)]=F
  
  #fit model
  if(missing("random")) {model0=brainstat.stats.terms$FixedEffect(model, "_check_categorical" = F)}
  else {model0=brainstat.stats.terms$MixedEffect(ran = random,fix = model,"_check_categorical" = F)}
  model=brainstat.stats.SLM$SLM(model = model0,
                                contrast=contrast,
                                surf = template,
                                mask=mask,
                                correction=c("fdr", "rft"),
                                cluster_threshold=p)
  model$fit(surf_data)
  
  #extracting tstats
  tstat=model$t
  
  ##extracting positive results
  cluster_pos=reticulate::py_to_r(model$P[["clus"]][[1]]) #pulling out results from brainstat's output
  cluster_pos=cluster_pos[cluster_pos$P<p,] #removing clusters that are not significant
  
  #extracting positive cluster map
  pos_clusterIDmap=model$P$clusid[[1]]
  
  if(NROW(cluster_pos)==0) #if no sig clusters emerged
  {
    cluster_pos="No significant clusters"
    pos_clusterIDmap=rep(0, NCOL(surf_data))
  } else
  {
    #creating new result variables in the cluster_pos objects
    cluster_pos$P=round(cluster_pos$P,3)
    cluster_pos$P[cluster_pos$P==0]="<0.001"
    cluster_pos=cluster_pos[ , !(names(cluster_pos) %in% "resels")] #removing the 'resels' column from the original brainstat output
    cluster_pos$X=NA
    cluster_pos$Y=NA
    cluster_pos$Z=NA
    cluster_pos$tstat=NA
    cluster_pos$region=NA
    
    #entering results for each cluster
    for (clusno in cluster_pos$clusid)
    {
      clus_tstat=tstat
      clus_tstat[pos_clusterIDmap!=clusno]=0
      cluster_pos$tstat[clusno]=round(clus_tstat[which.max(clus_tstat)],2)
      cluster_pos[clusno,4:6]=round(model$coord[,which.max(abs(clus_tstat))],1)
      
      #identifying region by matching the indices
      idx_pos=ROImap[[1]][,atlas][which.max(clus_tstat)]
      if(idx_pos>0){cluster_pos$region[clusno]=ROImap[[2]][,atlas][idx_pos] } ##to deal with desikan atlas missing vertex mappings
      else {cluster_pos$region[clusno]="unknown (use another atlas)"}
      
      remove(clus_tstat,idx_pos)
    }
    #thresholding positive cluster map
    pos_clusterIDmap[pos_clusterIDmap>max(cluster_pos$clusid)]=0
  }
  
  ##extracting negative results
  cluster_neg=reticulate::py_to_r(model$P[["clus"]][[2]]) #pulling out results from brainstat's output
  cluster_neg=cluster_neg[cluster_neg$P<p,] #removing clusters that are not significant
  
  #extracting negative cluster map
  neg_clusterIDmap=model$P$clusid[[2]]
  if(NROW(cluster_neg)==0) #if no sig clusters emerged
  {
    cluster_neg="No significant clusters"
    neg_clusterIDmap=rep(0, NCOL(surf_data))
  } else
  { #creating new result variables in the cluster_pos objects
    cluster_neg$P=round(cluster_neg$P,3)
    cluster_neg$P[cluster_neg$P==0]="<0.001"
    cluster_neg=cluster_neg[ , !(names(cluster_neg) %in% "resels")] #removing the 'resels' column from the original brainstat output
    cluster_neg$X=NA
    cluster_neg$Y=NA
    cluster_neg$Z=NA
    cluster_neg$tstat=NA
    cluster_neg$region=NA
    
    #entering results for each cluster
    for (clusno in cluster_neg$clusid)
    {
      clus_tstat=tstat
      clus_tstat[neg_clusterIDmap!=clusno]=0
      cluster_neg$tstat[clusno]=round(clus_tstat[which.min(clus_tstat)],2)
      cluster_neg[clusno,4:6]=round(model$coord[,which.max(abs(clus_tstat))],1)
      
      #identifying region by matching the indices
      idx_neg=ROImap[[1]][,atlas][which.min(clus_tstat)]
      if(idx_neg>0){cluster_neg$region[clusno]=ROImap[[2]][,atlas][idx_neg] } ##to deal with desikan atlas missing vertex mappings
      else {cluster_neg$region[clusno]="unknown (use another atlas)"}
      
      remove(clus_tstat,idx_neg)
    }
    #thresholding negative cluster map
    neg_clusterIDmap[neg_clusterIDmap>max(cluster_neg$clusid)]=0
  }
  ##combining results from both clusters into a list object
  cluster_results=list(cluster_pos,cluster_neg)
  names(cluster_results)=c("Positive contrast","Negative contrast")
  
  ##generating thresholded t-stat vector (for plotting)
  tstat[intersect(which(neg_clusterIDmap==0),which(pos_clusterIDmap==0))]=NA
  tstat[is.na(tstat)]=0
  tstat[maskNA]=NA
  
  ##generating positive and negative masks
  posmask=array(rep(0,NCOL(surf_data)))
  posmask[which(tstat>0)]=1
  
  negmask=array(rep(0,NCOL(surf_data)))
  negmask[which(tstat<0)]=1
  
  #listing objects to return
  returnobj=list(cluster_results,tstat,posmask,negmask,pos_clusterIDmap,neg_clusterIDmap)
  names(returnobj)=c("cluster_level_results","thresholded_tstat_map","pos_mask","neg_mask","pos_clusterIDmap","neg_clusterIDmap")
  return(returnobj)
}