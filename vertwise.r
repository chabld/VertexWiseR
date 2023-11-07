## WRAPPERS FOR VERTEX-WISE FUNCTIONS
## FOR USE IN THE COGNITIVE AND BRAIN HEALTH LABORATORY

############################################################################################################################
############################################################################################################################
##vertex wise analysis

vertex_analysis=function(all_predictors,IV_of_interest, CT_data, p=0.05)
{
  list.of.packages <- c("label4MRI", "reticulate")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) 
  {
    cat(paste("The following package(s) are required and will be installed:\n",new.packages,"\n"))
    install.packages(new.packages)
  }  
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
  cluster_pos=py_to_r(model$P[["clus"]][[1]])
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
      cluster_pos$region[clusno]=mni_to_region_name(cluster_pos[clusno,4],cluster_pos[clusno,5],cluster_pos[clusno,6])$aal.label
      remove(clus_tstat)
    }
    pos_clusterIDmap=model$P$clusid[[1]]
    pos_clusterIDmap[pos_clusterIDmap>max(cluster_pos$clusid)]=0
  }
  ##extracting negative results
  #cluster level
  cluster_neg=py_to_r(model$P[["clus"]][[2]])
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
      cluster_neg$region[clusno]=mni_to_region_name(cluster_neg[clusno,4],cluster_neg[clusno,5],cluster_neg[clusno,6])$aal.label
      remove(clus_tstat)
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
plotCT=function(data, fs_path, filename, surface="inflated", hot="#F8766D", cold="#00BFC4", limits)
{
  list.of.packages <- "fsbrain"
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) 
  {
    cat(paste("The following package are required and will be installed:\n",new.packages,"\n"))
    install.packages(new.packages)
  }  
  if(fs_path== "fsaverage5")
  {
    fs_path=paste(getwd(),"/fsaverage5", sep="")
  } else if(fs_path== "fsaverage5/")
  {
    fs_path=paste(getwd(),"/fsaverage5", sep="")
  }
  if(length(data) != 20484)
  {
    stop("Data has to be a numeric vector with 20484 values")
  } 
  
  if(file.exists(fs_path)== F)
  {
    stop("fs_path does not exist")
  }   
    
  if(range(data,na.rm = T)[1]>=0)
    {
      limits=c(0,max(data,na.rm = T))
      symm=F
      colfunc=colorRampPalette(c("white",hot))
    } 
    else if (range(data,na.rm = T)[2]<=0)
    {
      limits=c(-max(abs(data),na.rm = T),0)
      colfunc=colorRampPalette(c(cold,"white"))
      symm=F
    } 
    else 
    {
      limits=c(-max(abs(data),na.rm = T), max(abs(data),na.rm = T))
      colfunc=colorRampPalette(c(cold,"white",hot))
      symm=T
    }

  plotCT=vis.data.on.subject(gsub("fsaverage5","",fs_path), "fsaverage5", morph_data_both = data, surface=surface, 
                             views=NULL, makecmap_options = list('colFn'=colfunc, range=limits,symm=symm,col.na="gray80"))
  img=suppressWarnings(export(plotCT,output_img = filename, grid=F, silent=T))
}

plotCT2=function(data, filename)
{
  list.of.packages <- "reticulate"
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) 
  {
    cat(paste("The following package are required and will be installed:\n",new.packages,"\n"))
    install.packages(new.packages)
  }  
  if(length(data) != 20484)
  {
    stop("Data has to be a numeric vector with 20484 values")
  } 
  
  data[data==0]=NaN
  brainstat.datasets=reticulate::import("brainstat.datasets")  
  brainspace.plotting=reticulate::import("brainspace.plotting")  
  
  left=brainstat.datasets$fetch_template_surface("fsaverage5", join=F)[1]
  right=brainstat.datasets$fetch_template_surface("fsaverage5", join=F)[2]
  
  if(range(data,na.rm = T)[1]>=0)
  {
    CTplot=brainspace.plotting$plot_hemispheres(left[[1]], right[[1]],  array_name=np_array(data),cmap="Reds", nan_color=tuple(as.integer(c(1,1,1,1))),
                                            size=tuple(as.integer(c(1920,500))),return_plotter=T,background=tuple(as.integer(c(1,1,1))),
                                            interactive=F, color_bar=T,  transparent_bg=FALSE)
    CTplot$screenshot(filename=filename,transparent_bg = F)
  } 
  else if (range(data,na.rm = T)[2]<=0)
  {
    CTplot=brainspace.plotting$plot_hemispheres(left[[1]], right[[1]],  array_name=np_array(data),cmap="Blues", nan_color=tuple(as.integer(c(1,1,1,1))),
                                                size=tuple(as.integer(c(1920,500))),return_plotter=T,background=tuple(as.integer(c(1,1,1))),
                                                interactive=F, color_bar=T,  transparent_bg=FALSE)
    CTplot$screenshot(filename=filename,transparent_bg = F)
  } 
  else 
  {
    CTplot=brainspace.plotting$plot_hemispheres(left[[1]], right[[1]],  array_name=np_array(data),cmap="RdBu", nan_color=tuple(as.integer(c(1,1,1,1))),color_range="sym",
                                                size=tuple(as.integer(c(1920,500))),return_plotter=T,background=tuple(as.integer(c(1,1,1))),
                                                interactive=F, color_bar=T,  transparent_bg=FALSE)
    CTplot$screenshot(filename=filename,transparent_bg = F)
  }
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
  tflow=import("templateflow.api")
  np=import("numpy")
  interpolate=import("brainstat.mesh.interpolate")
  discrete=import("nimare.decode")
  nimare.dataset=import("nimare.dataset")
  
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
    cat("neurosynth_dataset.pkl.gz is not detected in the current working directory. The neurosynth database will be downloaded\n")
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
  
