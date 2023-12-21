## Find clusters from fsaverage5 template
## FOR USE IN THE COGNITIVE AND BRAIN HEALTH LABORATORY
############################################################################################################################
############################################################################################################################
getClusters=function(data)
{
  vert.all=which(abs(data)>0)
  lh.vert.all=vert.all[vert.all<10243]
  rh.vert.all=vert.all[vert.all>10242]
  edgelist.all=matrix(NA,nrow=0,ncol=2)
  
  ##LH
  if(length(lh.vert.all)>1)
  {
    for (vert.no in 1:length(lh.vert.all))
    {
      vert.connected=lh.vert.all[match(fs5_adj[[lh.vert.all[vert.no]]],lh.vert.all)]
      if(anyNA(vert.connected))
      {
        vert.connected=vert.connected[-which(is.na(vert.connected))] 
      }
      if(length(vert.connected)>0)
      {
        edgelist=matrix(NA,nrow=length(vert.connected),ncol=2)
        for (edge.idx in 1:length(vert.connected))
        {
          edgelist[edge.idx,1]=lh.vert.all[[vert.no]]
          edgelist[edge.idx,2]=vert.connected[edge.idx]
        }
        edgelist.all=rbind(edgelist.all,edgelist)
      }
    }
    if(NROW(edgelist.all)>0)
    {
      idx= !duplicated(t(apply(edgelist.all,  1, sort)))
      edgelist.all=edgelist.all[idx,]
      names(edgelist.all)=c("N1","N2")
      if(length(edgelist.all)==2)
      {
        LH.clust.map=rep(NA,20484)
        LH.clust.map[edgelist.all[1]]=1
        LH.clust.map[edgelist.all[2]]=1
        LH.clust.size=2
      } else
      {
        LH.com=igraph::components(igraph::graph.data.frame(edgelist.all, directed = F))
        LH.clust.size=LH.com$csize
        LH.clust.map=rep(NA,20484)
        for(clust.no in 1:LH.com$no)
        {
          LH.clust.map[as.numeric(names(which(LH.com$membership==clust.no)))]=clust.no
        } 
      }
    } else
    {
      LH.clust.map="noclusters"
      LH.clust.size=NA
    }
  } else 
  {
    LH.clust.map="noclusters"
    LH.clust.size=NA
  }
  remove(edgelist)
  edgelist.all=matrix(NA,nrow=0,ncol=2)
  ##RH
  if(length(rh.vert.all)>1)
  {
    for (vert.no in 1:length(rh.vert.all))
    {
      vert.connected=rh.vert.all[match(fs5_adj[[rh.vert.all[vert.no]]],rh.vert.all)]
      if(anyNA(vert.connected))
      {
        vert.connected=vert.connected[-which(is.na(vert.connected))] 
      }
      if(length(vert.connected)>0)
      {
        edgelist=matrix(NA,nrow=length(vert.connected),ncol=2)
        for (edge.idx in 1:length(vert.connected))
        {
          edgelist[edge.idx,1]=rh.vert.all[[vert.no]]
          edgelist[edge.idx,2]=vert.connected[edge.idx]
        }
        edgelist.all=rbind(edgelist.all,edgelist)
      }
    }
    if(NROW(edgelist.all)>0)
    {
      idx= !duplicated(t(apply(edgelist.all,  1, sort)))
      edgelist.all=edgelist.all[idx,]
      names(edgelist.all)=c("N1","N2")
      if(length(edgelist.all)==2)
      {
        RH.clust.map=rep(NA,20484)
        RH.clust.map[edgelist.all[1]]=1
        RH.clust.map[edgelist.all[2]]=1
        RH.clust.size=2
      } else
      {
        RH.com=igraph::components(igraph::graph.data.frame(edgelist.all, directed = F))
        RH.clust.size=RH.com$csize
        RH.clust.map=rep(NA,20484)
        for(clust.no in 1:RH.com$no)
        {
          RH.clust.map[as.numeric(names(which(RH.com$membership==clust.no)))]=clust.no
        } 
      }
    } else
    {
      RH.clust.map="noclusters"
      RH.clust.size=NA
    } 
  } else
  {
    RH.clust.map="noclusters"
    RH.clust.size=NA
  }
  if(!is.na(LH.clust.size[1]) & is.na(RH.clust.size[1]))
  {
    clust.map=LH.clust.map
    clust.size=LH.clust.size
  } else if(is.na(LH.clust.size[1]) & !is.na(RH.clust.size[1]))
  {
    clust.map=RH.clust.map
    clust.size=RH.clust.size
  } else if(!is.na(LH.clust.size[1]) & !is.na(RH.clust.size[1]))
  {
    RH.clust.map[which(RH.clust.map>0)]=RH.clust.map[which(RH.clust.map>0)]+max(LH.clust.map,na.rm = T)
    RH.clust.map[is.na(RH.clust.map)]=0
    LH.clust.map[is.na(LH.clust.map)]=0
    clust.map=LH.clust.map+RH.clust.map
    clust.map[clust.map==0]=NA
    clust.size=c(LH.clust.size,RH.clust.size)
  } else if(is.na(LH.clust.size[1]) & is.na(RH.clust.size[1]))
  {
    clust.map="noclusters"
    clust.size="noclusters"
  }
  return(list(clust.map,clust.size))
}
############################################################################################################################
############################################################################################################################
getClustersOLD=function(data)
{
  vert.all=which(abs(data)>0)
  lh.vert.all=vert.all[vert.all<10243]
  rh.vert.all=vert.all[vert.all>10242]
  edgelist.all=matrix(NA,nrow=0,ncol=2)
  
  ##LH
  if(length(lh.vert.all)>1)
  {
    for (vert.no in 1:length(lh.vert.all))
    {
      vert.connected=lh.vert.all[match(fs5_adj[[lh.vert.all[vert.no]]],lh.vert.all)]
      if(anyNA(vert.connected))
      {
        vert.connected=vert.connected[-which(is.na(vert.connected))] 
      }
      if(length(vert.connected)>0)
      {
        edgelist=matrix(NA,nrow=length(vert.connected),ncol=2)
        for (edge.idx in 1:length(vert.connected))
        {
          edgelist[edge.idx,1]=lh.vert.all[[vert.no]]
          edgelist[edge.idx,2]=vert.connected[edge.idx]
        }
        edgelist.all=rbind(edgelist.all,edgelist)
      }
    }
  }
  ##RH
  if(length(rh.vert.all)>1)
  {
    for (vert.no in 1:length(rh.vert.all))
    {
      vert.connected=rh.vert.all[match(fs5_adj[[rh.vert.all[vert.no]]],rh.vert.all)]
      if(anyNA(vert.connected))
      {
        vert.connected=vert.connected[-which(is.na(vert.connected))] 
      }
      if(length(vert.connected)>0)
      {
        edgelist=matrix(NA,nrow=length(vert.connected),ncol=2)
        for (edge.idx in 1:length(vert.connected))
        {
          edgelist[edge.idx,1]=rh.vert.all[[vert.no]]
          edgelist[edge.idx,2]=vert.connected[edge.idx]
        }
        edgelist.all=rbind(edgelist.all,edgelist)
      }
    }
  }
  if(NROW(edgelist.all)>0)
  {
    idx= !duplicated(t(apply(edgelist.all,  1, sort)))
    edgelist.all=edgelist.all[idx,]
    names(edgelist.all)=c("N1","N2")
    if(length(edgelist.all)==2)
    {
      clust.map=clust.map=rep(NA,20484)
      clust.map[edgelist.all[1]]=1
      clust.map[edgelist.all[2]]=1
      clust.size=2
    } else
    {
      com=igraph::components(igraph::graph.data.frame(edgelist.all, directed = F))
      clust.size=com$csize
      clust.map=rep(NA,20484)
      for(clust.no in 1:com$no)
      {
        clust.map[as.numeric(names(which(com$membership==clust.no)))]=clust.no
      } 
    }
  } else
  {
    clust.map="noclusters"
    clust.size="noclusters"
  }
  return(list(clust.map,clust.size))
}
