tab<- read.table("details.txt")
       library("phangorn")
       library("ape")
       i<- 1
    while(i<nrow(tab)+1){
                          main<-read.tree("supermatrix_partition_iqtree_reroot_no_bootstraps.tree")
    a<- strsplit(as.matrix(tab[i,]),"\\|")
        genes<-a[[1]][1]
    tfl<-"tfl"
          j<-2
    while(j<length(a[[1]])+1){
                               main<-drop.tip(main,a[[1]][j])
    j<-j+1
           }
           rfoutput=paste(genes,tfl,sep=".")
           write.tree(file=rfoutput,main)
           i<-i+1
    }
