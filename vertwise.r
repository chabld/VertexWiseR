## FUNCTION FOR VERTEX-WISE ANALYSIS
## ADAPTED FROM brainstat python library (https://brainstat.readthedocs.io/en/master/_modules/brainstat/stats/SLM.html#SLM)
## FOR USE IN THE COGNITIVE AND BRAIN HEALTH LABORATORY

############################################################################################################################
############################################################################################################################
##vertex wise analysis
vertex_analysis=function(all_predictors,IV_of_interest, CT_data, p=0.05, atlas=1)  ## atlas: 1=Desikan, 2=Schaefer-100, 3=Schaefer-200, 4=Glasser-360, 5=Destrieux-148
{
    ##checks
    # check required packages
    list.of.packages ="reticulate"
    new.packages = list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
    if(length(new.packages)) 
    {
      cat(paste("The following package(s) are required and will be installed:\n",new.packages,"\n"))
      install.packages(new.packages)
    }
    #check if nrow is consistent for all_predictors and FC_data
    if(NROW(CT_data)!=NROW(all_predictors))  {stop(paste("The number of rows for CT_data (",NROW(CT_data),") and all_predictors (",NROW(all_predictors),") are not the same",sep=""))}
    #check categorical variable
    for (column in 1:NCOL(all_predictors))
    {
      if(class(all_predictors[,column])  != "integer" & class(all_predictors[,column])  != "numeric")  {stop(paste(colnames(all_predictors)[column],"is not a numeric variable, please recode it into a numeric variable"))}
    }
  
    #incomplete data check
    idxF=which(complete.cases(all_predictors)==F)
    if(length(idxF)>0)
    {
      cat(paste("all_predictors contains",length(idxF),"subjects with incomplete data. Subjects with incomplete data will be excluded in the current analysis"))
      all_predictors=all_predictors[-idxF,]
      IV_of_interest=IV_of_interest[-idxF]
      CT_data=CT_data[-idxF,]
    }
    remove(colno)

  ##import python libaries
  brainstat.stats.terms=reticulate::import("brainstat.stats.terms")
  brainstat.stats.SLM=reticulate::import("brainstat.stats.SLM")
  brainstat.datasets=reticulate::import("brainstat.datasets")  
  
  ##fitting model
  mask=array(rep(T,NCOL(CT_data)))
  maskNA=which(colSums(CT_data != 0) == 0)
  mask[which(colSums(CT_data != 0) == 0)]=F
  model0 = brainstat.stats.terms$FixedEffect(all_predictors, "_check_categorical" = F)
  model=brainstat.stats.SLM$SLM(model = model0,
                                contrast=IV_of_interest,
                                surf = "fsaverage5", 
                                mask=mask,
                                correction=c("fdr", "rft"),
                                cluster_threshold=p)
  model$fit(CT_data)
  
  tstat=model$t
  ##extracting positive results
  cluster_pos=reticulate::py_to_r(model$P[["clus"]][[1]])
  cluster_pos=cluster_pos[cluster_pos$P<p,]
  if(NROW(cluster_pos)==0)
  {
    cluster_pos="No significant clusters"
    pos_clusterIDmap=rep(0, NCOL(CT_data))
  } else
  {
    cluster_pos$P=round(cluster_pos$P,3)
    cluster_pos$P[cluster_pos$P==0]="<0.001"
    cluster_pos=cluster_pos[ , !(names(cluster_pos) %in% "resels")]
    cluster_pos$X=NA
    cluster_pos$Y=NA
    cluster_pos$Z=NA
    cluster_pos$tstat=NA
    cluster_pos$region=NA
    
    pos_clusterIDmap=model$P$clusid[[1]]
    for (clusno in cluster_pos$clusid)
    {
      clus_tstat=tstat
      clus_tstat[pos_clusterIDmap!=clusno]=0
      cluster_pos$tstat[clusno]=round(clus_tstat[which.max(clus_tstat)],2)
      cluster_pos[clusno,4:6]=round(model$coord[,which.max(abs(clus_tstat))],1)
  
      load(file = url("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/data/ROImap.rdata?raw=TRUE"))
      idx_pos=ROImap[[1]][,atlas][which.max(clus_tstat)]
      if(idx_pos>0){cluster_pos$region[clusno]=ROImap[[2]][,atlas][idx_pos] } ##to deal with desikan atlas missing vertex mappings
      else {cluster_pos$region[clusno]="unknown (use another atlas)"}
      
      remove(clus_tstat,idx_pos)
    }
    pos_clusterIDmap=model$P$clusid[[1]]
    pos_clusterIDmap[pos_clusterIDmap>max(cluster_pos$clusid)]=0
  }
  ##extracting negative results
  cluster_neg=reticulate::py_to_r(model$P[["clus"]][[2]])
  cluster_neg=cluster_neg[cluster_neg$P<p,]
  if(NROW(cluster_neg)==0)
  {
    cluster_neg="No significant clusters"
    neg_clusterIDmap=rep(0, NCOL(CT_data))
  } else
  {
    cluster_neg$P=round(cluster_neg$P,3)
    cluster_neg$P[cluster_neg$P==0]="<0.001"
    cluster_neg=cluster_neg[ , !(names(cluster_neg) %in% "resels")]
    cluster_neg$X=NA
    cluster_neg$Y=NA
    cluster_neg$Z=NA
    cluster_neg$tstat=NA
    cluster_neg$region=NA
    
    neg_clusterIDmap=model$P$clusid[[2]]
    for (clusno in cluster_neg$clusid)
    {
      clus_tstat=tstat
      clus_tstat[neg_clusterIDmap!=clusno]=0
      cluster_neg$tstat[clusno]=round(clus_tstat[which.min(clus_tstat)],2)
      cluster_neg[clusno,4:6]=round(model$coord[,which.max(abs(clus_tstat))],1)
      if(!exists("ROImap", inherit=F))
        {
        load(file = url("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/data/ROImap.rdata?raw=TRUE"))
        } 
      idx_neg=ROImap[[1]][,atlas][which.min(clus_tstat)]
      if(idx_neg>0){cluster_neg$region[clusno]=ROImap[[2]][,atlas][idx_neg] } ##to deal with desikan atlas missing vertex mappings
      else {cluster_neg$region[clusno]="unknown (use another atlas)"}

      remove(clus_tstat,idx_neg)
    }
    neg_clusterIDmap=model$P$clusid[[2]]
    neg_clusterIDmap[neg_clusterIDmap>max(cluster_neg$clusid)]=0
  }
  cluster_results=list(cluster_pos,cluster_neg)
  names(cluster_results)=c("Positive contrast","Negative contrast")
  
  tstat[intersect(which(neg_clusterIDmap==0),which(pos_clusterIDmap==0))]=NA
  tstat[is.na(tstat)]=0
  tstat[maskNA]=NA
  
  posmask=array(rep(0,NCOL(CT_data)))
  posmask[which(tstat>0)]=1
  
  negmask=array(rep(0,NCOL(CT_data)))
  negmask[which(tstat<0)]=1
  
  returnobj=list(cluster_results,tstat,posmask,negmask,pos_clusterIDmap,neg_clusterIDmap)
  
  names(returnobj)=c("cluster_level_results","thresholded_tstat_map","pos_mask","neg_mask","pos_clusterIDmap","neg_clusterIDmap")
  
  return(returnobj)
}
############################################################################################################################
############################################################################################################################
##load other vertex-wise functions
source("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/otherfunc.r?raw=TRUE")

##example
# source("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/vertwise.r?raw=TRUE")
# 
# model=vertex_analysis(all_predictors =all_pred, IV_of_interest = dat_beh$Sex, CT_data = dat_CT)
# model$cluster_level_results

