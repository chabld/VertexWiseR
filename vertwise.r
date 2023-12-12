## WRAPPERS FOR VERTEX-WISE FUNCTIONS
## FOR USE IN THE COGNITIVE AND BRAIN HEALTH LABORATORY

############################################################################################################################
############################################################################################################################
## To extract atlas ROI values from fsaverage5 vertex-wise data
fs5_to_atlas=function(data,atlas) ## atlas: 1=Desikan, 2=Schaefer-100, 3=Schaefer-200, 4=Glasser-360, 5=Destrieux-148
{
  load(file = url("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/ROImap.rdata?raw=TRUE"))
  nregions=max(ROImap[[1]][,atlas])
  if(length(data)%%20484!=0)
  {
    stop("Length of data is not a multiple of 20484")
  }
  data[is.na(data)]=0
  if(length(data)==20484)
  {
  data=matrix(data,ncol=20484,nrow=1)  
  ROI=rep(NA,nregions)
    for (region in 1:nregions)
    {
      ROI[region]=mean(data[which(ROImap[[1]][,atlas]==region)])
    }
  } else 
  {
    ROI=matrix(NA, nrow=NROW(data), ncol=nregions)
    for (region in 1:nregions)
    {
      ROI[,region]=rowMeans(data[,which(ROImap[[1]][,atlas]==region)])
    }
  }
  return(ROI)
}

atlas_to_fs5=function(data,atlas) ## atlas: 1=Desikan, 2=Schaefer-100, 3=Schaefer-200, 4=Glasser-360
  {
  load(file = url("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/ROImap.rdata?raw=TRUE"))
  nregions=max(ROImap[[1]][,atlas])
  fs5_dat=rep(NA,20484)
  for (region in 1:nregions)
    {
    fs5_dat[which(ROImap[[1]][,atlas]==region)]=data[region]
    }
  return(as.numeric(fs5_dat))
  }
############################################################################################################################
############################################################################################################################
##smoothing fsaverage5 and fsaverage6 data
smooth=function(data,FWHM=10)
{
  brainstat.mesh.data=reticulate::import("brainstat.mesh.data")
  brainstat.datasets=reticulate::import("brainstat.datasets")  
  col0=which(colSums(CT_dat==0) == nrow(CT_dat))
  if(ncol(CT_dat)==20484)
    {
    surftemp=brainstat.datasets$fetch_template_surface("fsaverage5", join=T)
    vert_mm=3.5
    }
  if(ncol(CT_dat)==81924)
    {
    surftemp=brainstat.datasets$fetch_template_surface("fsaverage6", join=T)
    vert_mm=2
    }
  smooth=brainstat.mesh.data$mesh_smooth(Y=CT_dat,surf=surftemp, FWHM = FWHM/vert_mm)
  smooth[,col0]=0
  return(smooth)  
}

############################################################################################################################
############################################################################################################################
##vertex wise analysis

vertex_analysis=function(all_predictors,IV_of_interest, CT_data, p=0.05, atlas=1)  ## atlas: 1=Desikan, 2=Schaefer-100, 3=Schaefer-200, 4=Glasser-360, 5=Destrieux-148
{
  for (column in 1:NCOL(all_predictors))
  {
    if(class(all_predictors[,column])  != "integer" & class(all_predictors[,column])  != "numeric")
    {
      stop(paste(colnames(all_predictors)[column],"is not a numeric variable, please recode it into a numeric variable"))
    }
  }
  
  brainstat.stats.terms=reticulate::import("brainstat.stats.terms")
  brainstat.stats.SLM=reticulate::import("brainstat.stats.SLM")
  brainstat.datasets=reticulate::import("brainstat.datasets")  
  
  ##incomplete data check
  idxF=which(complete.cases(all_predictors)==F)

  if(length(idxF)>0)
  {
    cat(paste("all_predictors contains",length(idxF),"subjects with incomplete data. Subjects with incomplete data will be excluded in the current analysis"))
    all_predictors=all_predictors[-idxF,]
    IV_of_interest=IV_of_interest[-idxF]
    CT_data=CT_data[-idxF,]
  }
  
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
  #cluster level
  cluster_pos=reticulate::py_to_r(model$P[["clus"]][[1]])
  cluster_pos=cluster_pos[cluster_pos$P<p,]
  if(NROW(cluster_pos)==0)
  {
    cluster_pos="No significant clusters"
    pos_clusterIDmap=rep(0, NCOL(CT_data))
  } else
  {
    cluster_pos$P=round(cluster_pos$P,3)
    cluster_pos$P[cluster_pos$P==0]="<.001"
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
  
      load(file = url("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/ROImap.rdata?raw=TRUE"))
      idx_pos=ROImap[[1]][,atlas][which.max(clus_tstat)]
      cluster_pos$region[clusno]=ROImap[[2]][,atlas][idx_pos]
      
      remove(clus_tstat,idx_pos)
    }
    pos_clusterIDmap=model$P$clusid[[1]]
    pos_clusterIDmap[pos_clusterIDmap>max(cluster_pos$clusid)]=0
  }
  ##extracting negative results
  #cluster level
  cluster_neg=reticulate::py_to_r(model$P[["clus"]][[2]])
  cluster_neg=cluster_neg[cluster_neg$P<p,]
  if(NROW(cluster_neg)==0)
  {
    cluster_neg="No significant clusters"
    neg_clusterIDmap=rep(0, NCOL(CT_data))
  } else
  {
    
    cluster_neg$P=round(cluster_neg$P,3)
    cluster_neg$P[cluster_neg$P==0]="<.001"
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
        load(file = url("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/ROImap.rdata?raw=TRUE"))
        } 
      idx_neg=ROImap[[1]][,atlas][which.min(clus_tstat)]
      cluster_neg$region[clusno]=ROImap[[2]][,atlas][idx_neg]

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
##CT surface plots

plotCT=function(data, filename,title="",surface="inflated",cmap,fs_path, range=NULL)
{
  if(length(data) != 20484)
  {
    stop("Data has to be a numeric vector with 20484 values")
  } 
  if(!missing("fs_path")){cat("The fs_path parameter and the fsaverage5 files are no longer needed in the updated plotCT function\n")}
  
  brainstat.datasets=reticulate::import("brainstat.datasets")  
  brainspace.plotting=reticulate::import("brainspace.plotting")  

  left=brainstat.datasets$fetch_template_surface("fsaverage5", join=F, layer=surface)[1]
  right=brainstat.datasets$fetch_template_surface("fsaverage5", join=F, layer=surface)[2]
  
  if(missing("cmap"))
  {
    if(range(data,na.rm = T)[1]>=0)
      {
      cmap="Reds"
      range=NULL
      }
    else if (range(data,na.rm = T)[2]<=0)
      {cmap="Blues"
      range=NULL
      }
    else
      {
      cmap="RdBu_r"
      range="sym"
      }  
  }
  
  CTplot=brainspace.plotting$plot_hemispheres(left[[1]], right[[1]],  array_name=reticulate::np_array(data),cmap=cmap, 
                                              size=reticulate::tuple(as.integer(c(1920,400))),nan_color=reticulate::tuple(c(0.7, 0.7, 0.7, 1)),
                                              return_plotter=T,background=reticulate::tuple(as.integer(c(1,1,1))),zoom=1.25,color_range=range,
                                              label_text=list('left'=list(title)),interactive=F, color_bar=T,  transparent_bg=FALSE)
  CTplot$screenshot(filename=filename,transparent_bg = F)

}

############################################################################################################################
############################################################################################################################
##CT image decoding
decode_img=function(img,contrast="positive")
{
  ##input checks
  if(length(img) != 20484)
  {
    stop("img object has to be a numeric vector with 20484 values")
  } 
  if(contrast != "positive" & contrast != "negative" )
  {
    stop("contrast has to be either positive or negative")
  } 
  tflow=reticulate::import("templateflow.api")
  np=reticulate::import("numpy")
  interpolate=reticulate::import("brainstat.mesh.interpolate")
  discrete=reticulate::import("nimare.decode")
  nimare.dataset=reticulate::import("nimare.dataset")
  
  ##selecting contrasts
  if(contrast=="positive")
  {
    img[is.na(img)]=0
    img[img<0]=0
    img[img>0]=1
    stat_labels=r_to_py(img)
    stat_nii = interpolate$`_surf2vol`("fsaverage5", stat_labels$flatten())
  } else if (contrast=="negative")
  {
    img[is.na(img)]=0
    img[img>0]=0
    img[img<0]=1
    stat_labels=r_to_py(img)
    stat_nii = interpolate$`_surf2vol`("fsaverage5", stat_labels$flatten())
  }
  
  ##download neurosynth database if necessary 
  if(file.exists("neurosynth_dataset.pkl.gz")==F)
  {
    cat("\nneurosynth_dataset.pkl.gz is not detected in the current working directory. The neurosynth database will be downloaded\n")
    download.file(url="https://blogs.ntu.edu.sg/cogbrainhealthlab/files/2023/10/neurosynth_dataset.pkl_.gz",destfile = "neurosynth_dataset.pkl.gz")
  } 
  ##running the decoding process
  neurosynth_dset = nimare.dataset$Dataset$load("neurosynth_dataset.pkl.gz")
  cat("Correlating input image with images in the neurosynth database. This may take a while\n")
  decoder = discrete$ROIAssociationDecoder(stat_nii)
  decoder$fit(neurosynth_dset)
  decoder_df = data.matrix(decoder$transform())
  
  ##compiling the results
  row.names(decoder_df)=gsub(pattern = "terms_abstract_tfidf__",x=row.names(decoder_df), replacement = "")
  result=data.frame(row.names(decoder_df),round(as.numeric(decoder_df),3))
  colnames(result)=c("keyword","r")
  result=result[order(-result$r),]
  
  return(result)
}  
  
