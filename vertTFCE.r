## FUNCTIONS FOR VERTEX-WISE TFCE ANALYSIS
## FOR USE IN THE COGNITIVE AND BRAIN HEALTH LABORATORY

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
        }
      }
    }
    #combine results from positive and negative tails if necessary 
    if(sign.idx==1){tfce_step_values.all=tfce}
    else if (sign.idx==2){tfce_step_values.all=tfce_step_values.all+tfce}
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
    cl=parallel::makeCluster(nthread)
    doParallel::registerDoParallel(cl)
    `%dopar%` = foreach::`%dopar%`

    #parallel loop across different score_threshs values for TFCE estimation
    tfce=foreach::foreach(thresh.no=1:length(score_threshs), .combine="rbind", .export=c("getClusters","fs5_edgelist"))  %dopar%
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
    suppressWarnings(closeAllConnections())
    #combine results from positive and negative tails if necessary 
    if(sign.idx==1){tfce_step_values.all=colSums(tfce)}
    else if (sign.idx==2){tfce_step_values.all=tfce_step_values.all+colSums(tfce)}
  }
  return(tfce_step_values.all)
  parallel::stopCluster(cl)
}
############################################################################################################################
############################################################################################################################
##Main function

TFCE.vertex_analysis=function(all_predictors,IV_of_interest, CT_data, nperm=5, tail=2, nthread=10)
{
  ##checks
    # check required packages
    list.of.packages = c("parallel", "doParallel","igraph","doSNOW","foreach")
    new.packages = list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
    if(length(new.packages)) 
    {
      cat(paste("The following package(s) are required and will be installed:\n",new.packages,"\n"))
      install.packages(new.packages)
    }  
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
  
  ##load edgelist data
  if(!exists("fs5_edgelist"))  {load(file = url("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/data/fs5edgelist.rdata?raw=TRUE"),envir = globalenv())} 
  
  ##unpermuted model
    all_predictors=data.matrix(all_predictors)
    mod=lm(CT_data~data.matrix(all_predictors))
  
    #identify contrast
    for(colno in 1:(NCOL(all_predictors)+1))
    {
      if(colno==(NCOL(all_predictors)+1))
      {stop("IV_of_interest is not contained within all_predictors")}
      if(identical(IV_of_interest,all_predictors[,colno]))
      {break}
    }
  
    #extract tstat and calculate tfce image
    start=Sys.time()
    cat("\nEstimating unpermuted TFCE image...")
    
    tmap.orig=extract.t(mod,colno+1)
    TFCE.orig=TFCE.multicore(data = tmap.orig,tail = tail,nthread=nthread)
    
    end=Sys.time()
    cat(paste("Completed in",round(difftime(end,start, units="secs"),1),"secs\nEstimating permuted TFCE images...\n",sep=" "))
  
  ##permuted models
    ##generating permutation sequences  
    permseq=matrix(NA, nrow=NROW(all_predictors), ncol=nperm)
    for (perm in 1:nperm)  {permseq[,perm]=sample.int(NROW(all_predictors))}
  
    #activate parallel processing
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
  
    TFCE.max=foreach::foreach(perm=1:nperm, .combine="rbind",.export=c("TFCE","extract.t","getClusters","fs5_edgelist"), .options.snow = opts)  %dopar%
      {
        all_predictors.permuted=all_predictors
        all_predictors.permuted[,colno]=all_predictors.permuted[permseq[,perm],colno]
        mod.permuted=lm(CT_data~data.matrix(all_predictors.permuted))
        tmap=extract.t(mod.permuted,colno+1)
        
        return(max(abs(TFCE(data = tmap,tail = tail))))
      }
  
    end=Sys.time()
    cat(paste("\nCompleted in :",round(difftime(end, start, units='mins'),1)," minutes \n",sep=""))
    suppressWarnings(closeAllConnections())
    
    ##saving list objects
    returnobj=list(tmap.orig,TFCE.orig, TFCE.max,tail)
    names(returnobj)=c("t_stat","TFCE.orig","TFCE.max","tail")
  
    return(returnobj)
}
############################################################################################################################
############################################################################################################################
TFCE.threshold=function(TFCE.output, p=0.05, atlas=1, k=20)
{
  nperm=length(TFCE.output$TFCE.max)

  #check if number of permutations is adequate
  if(nperm<1/p)
  {
    if(TFCE.output$tail==2){warning(paste("Not enough permutations were carried out to estimate the two-tailed p<",p*2," threshold precisely\nConsider setting an nperm to at least ",ceiling(1/p),sep=""))} 
    else{warning(paste("Not enough permutations were carried out to estimate the p<",p," threshold precisely\nConsider setting nperm to at least ",ceiling(1/p),sep=""))}
  }
  
  ##loading vertex mapping data
  if(!exists("ROImap", inherit=F))  {load(file = url("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/data/ROImap.rdata?raw=TRUE"))} 
  if(!exists("MNImap", inherit=F))  {load(file = url("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/data/MNImap.rdata?raw=TRUE"))} 
  if(!exists("fs5_edgelist", inherit=F))  {load(file = url("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/data/fs5edgelist.rdata?raw=TRUE"))} 
  
  ##generating p map
  tfce.p=rep(NA,20484)
  TFCE.output$t_stat[is.na(TFCE.output$t_stat)]=0
  for (vert in 1:20484)
  {
  tfce.p[vert]=length(which(TFCE.output$TFCE.max>abs(TFCE.output$TFCE.orig[vert])))/nperm
  }
  
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
    
    pos.clusters0=getClusters(pos.t_stat.thresholdedP)
    if(pos.clusters0[[2]][1]!="noclusters") #skip if no clusters detected
    {
      #applying k thresholding
      pos.clustID.remove=which(pos.clusters0[[2]]<k)
      pos.clusters0[[1]][which(!is.na(match(pos.clusters0[[1]],pos.clustID.remove)))]=NA
      
      #generating mask
      pos.clusters=getClusters(pos.clusters0[[1]])
      pos.clusters[[1]][is.na(pos.clusters[[1]])]=0
      pos.mask=rep(0,20484)
      
      if(pos.clusters[[2]][1]!="noclusters") {pos.mask[pos.clusters[[1]]>0]=1}
      
      #results table
      pos.clustermap=rep(NA,20484)
      
      if(pos.clusters[[2]][1]!="noclusters")
      {
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
          if(atlas.idx>0){pos.clust.results[clust.idx,8]=ROImap[[2]][,atlas][atlas.idx] } ##to deal with desikan atlas missing vertex mappings
          else {pos.clust.results[clust.idx,8]="unknown (use another atlas)"}
          
          clust.idx=clust.idx+1
        }
      } else if(pos.clusters[[2]][1]=="noclusters")
      { 
        pos.clust.results="No significant clusters"
        pos.clustermap="No significant clusters"
      }
    } else 
    {
      pos.clust.results="No significant clusters"
      pos.clustermap="No significant clusters"
      pos.mask=rep(0,20484)
    }
  } else if(TFCE.output$tail==-1)
  {
    pos.clust.results="Positive contrast not analyzed, only negative one-tailed TFCE statistics were estimated)"
    pos.clustermap="No significant clusters"
    pos.mask=rep(0,20484)
  } 
  
  ##negative cluster
  if(TFCE.output$tail==-1 |TFCE.output$tail==2)
  {
    #applying p thresholding
    neg.t_stat.thresholdedP=t_stat.thresholdedP
    neg.t_stat.thresholdedP[neg.t_stat.thresholdedP>0]=0
    
    neg.clusters0=getClusters(neg.t_stat.thresholdedP)
    if(neg.clusters0[[2]][1]!="noclusters") #skip if no clusters detected
    {
      #applying k thresholding
      neg.clustID.remove=which(neg.clusters0[[2]]<k)
      neg.clusters0[[1]][which(!is.na(match(neg.clusters0[[1]],neg.clustID.remove)))]=NA
      
      #generating mask
      neg.clusters=getClusters(neg.clusters0[[1]])
      neg.clusters[[1]][is.na(neg.clusters[[1]])]=0
      neg.mask=rep(0,20484)
      
      if(neg.clusters[[2]][1]!="noclusters") {neg.mask[neg.clusters[[1]]>0]=1}
      
      #results table
      neg.clustermap=rep(NA,20484)
      
      if(neg.clusters[[2]][1]!="noclusters")
      {
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
          if(neg.clust.results[clust.idx,3]==0) {neg.clust.results[clust.idx,3]="<0.001"}
          neg.clust.results[clust.idx,c(4,5,6)]=round(MNImap[,max.vert.idx],1)
          neg.clust.results[clust.idx,7]=round(abs(TFCE.output$t_stat[max.vert.idx]),2)
          
          atlas.idx=ROImap[[1]][,atlas][max.vert.idx]
          if(atlas.idx>0){neg.clust.results[clust.idx,8]=ROImap[[2]][,atlas][atlas.idx] } ##to deal with desikan atlas missing vertex mappings
          else {neg.clust.results[clust.idx,8]="unknown (use another atlas)"}
          
          clust.idx=clust.idx+1
        }
      } else if(neg.clusters[[2]][1]=="noclusters")
      { 
        neg.clust.results="No significant clusters"
        neg.clustermap="No significant clusters"
      }
    } else 
    {
      neg.clust.results="No significant clusters"
      neg.clustermap="No significant clusters"
      neg.mask=rep(0,20484)
    }
  } else if(TFCE.output$tail==-1)
  {
    neg.clust.results="Negative contrast not analyzed, only negative one-tailed TFCE statistics were estimated)"
    neg.clustermap="No significant clusters"
    neg.mask=rep(0,20484)
  } 
  ##saving list objects
  
  cluster_level_results=list(pos.clust.results,neg.clust.results)
  names(cluster_level_results)=c("Positive contrast", "Negative contrasts")
  
  t_stat.thresholdedPK=TFCE.output$t_stat*(pos.mask+neg.mask)
  
  returnobj=list(cluster_level_results, t_stat.thresholdedPK, pos.clustermap, neg.clustermap)  
  names(returnobj)=c("cluster_level_results","thresholded_tstat_map","pos_clustermap","neg_clustermap")
  returnobj$cluster_level_results
  
  return(returnobj)
}  
############################################################################################################################
############################################################################################################################
##load other vertex-wise functions
source("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/otherfunc.r?raw=TRUE")

##example

#pos=TFCE.vertex_analysis(all_predictors =all_pred, IV_of_interest = dat_beh$age, CT_data = dat_CT, tail=1, nperm=100, nthread = 10)
#neg=TFCE.vertex_analysis(all_predictors =all_pred, IV_of_interest = dat_beh$age, CT_data = dat_CT, tail=-1 ,nperm=100, nthread = 10)
#two=TFCE.vertex_analysis(all_predictors =all_pred, IV_of_interest = dat_beh$age, CT_data = dat_CT, tail=2 ,nperm=100, nthread = 10)

#pos.results=TFCE.threshold(pos)
#pos.results$cluster_level_results
#plotCT(pos.results$thresholded_tstat_map, filename="pos.png")

#neg.results=TFCE.threshold(neg)
#neg.results$cluster_level_results
#plotCT(neg.results$thresholded_tstat_map, filename="neg.png")

#two.results=TFCE.threshold(two)
#two.results$cluster_level_results
#plotCT(two.results$thresholded_tstat_map, filename="two.png")
