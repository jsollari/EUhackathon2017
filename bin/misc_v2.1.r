#autor:      Joao Sollari Lopes
#local:      INE, Lisboa
#criado:     13.03.2017
#modificado: 17.05.2017

## 1. FUNCTIONS
{
# Key:
# a - auxiliar
# d - dataframe
# f - file name
# i - iterator
# l - list
# n - vector size
# s - selected samples
# t - title
# v - vector
# w - weights

# Calculate correlations
# dat1 - data frame of variables to be analysed
# method - type of correlation to use {"pearson","spearman"}
# f1 - output filename
calc_cor = function(dat1,method="spearman",f1){
    library("Hmisc")
    lab1 = colnames(dat1)
    m1 = as.matrix(dat1)
    rres = rcorr(m1,type=method)
    for(i1 in names(rres)){
        diag(rres[[i1]]) = NA
        rres[[i1]][lower.tri(rres[[i1]])] = NA
    }
    rres$padj = matrix(p.adjust(rres$P,method="BH"),ncol=ncol(rres$P),
      byrow=FALSE,dimnames=list(lab1,lab1))
    f2 = paste(f1,".RDS",sep="")
    saveRDS(rres,f2)
    return(rres)
}

# Calculate correlations (dat1 vs dat2)
# dat1 - first data frame of variables to be analysed
# dat2 - second data frame of variables to be analysed
# method - type of correlation to use {"pearson","spearman"}
# f1 - output filename
calc_cor_v2 = function(dat1,dat2,method="spearman",f1){
    library("Hmisc")
    n1 = ncol(dat1)
    n2 = ncol(dat2)
    lab1 = colnames(dat1)
    lab2 = colnames(dat2)
    m1 = as.matrix(cbind(dat1,dat2))
    rres1 = rcorr(m1,type=method)
    rres2 = vector("list",length=4)
    names(rres2) = c("r","n","P","padj")
    for(i1 in names(rres1)){
        rres2[[i1]] = rres1[[i1]][1:n1,(1+n1):(n1+n2),drop=FALSE]
    }
    rres2$padj = matrix(p.adjust(rres2$P,method="BH"),ncol=ncol(rres2$P),
      byrow=FALSE,dimnames=list(lab1,lab2))
    f2 = paste(f1,".RDS",sep="")
    saveRDS(rres2,f2)
    return(rres2)
}

# Plot levelplots of distances (wide format)
# cor1 - correlation results between variables to be analysed
# stat - statistic to report {"pval","adjp"}
# reord - reorder classes according to hierarchical clustering
# mirror - rotate results to visualize them as in the matrix
# rm_xlab - remove x label information
# rm_ylab - remove y label information
# f1 - output filename
# main - title for plots
plot_levelplots_cor_v2 = function(cor1,stat="padj",reord=FALSE,mirror=TRUE,rm_xlab=FALSE,
  rm_ylab=FALSE,f1,main=NULL){
    library("lattice")
    library("grid")
    plotpval = function(...){
        panel.levelplot(...)
        tmp = which(pnt1 < 0.10 & pnt1 >= 0.050,arr.ind=TRUE)
        if(length(tmp) > 0)
            grid.points(tmp[,1],tmp[,2],pch="+",gp=gpar(cex=0.6))
        tmp = which(pnt1 < 0.05 & pnt1 >= 0.010,arr.ind=TRUE)
        if(length(tmp) > 0)
            grid.points(tmp[,1],tmp[,2],pch="*",gp=gpar(cex=0.9))
        tmp = which(pnt1 < 0.01 & pnt1 >= 0.001,arr.ind=TRUE)
        if(length(tmp) > 0){
            grid.points(tmp[,1]-0.15,tmp[,2],pch="*",gp=gpar(cex=0.9))
            grid.points(tmp[,1]+0.15,tmp[,2],pch="*",gp=gpar(cex=0.9))
        }
        tmp = which(pnt1 < 0.001,arr.ind=TRUE)
        if(length(tmp) > 0){
            grid.points(tmp[,1]-0.3,tmp[,2],pch="*",gp=gpar(cex=0.9))
            grid.points(tmp[,1],tmp[,2],pch="*",gp=gpar(cex=0.9))
            grid.points(tmp[,1]+0.3,tmp[,2],pch="*",gp=gpar(cex=0.9))
        }
    }
    m1 = as.matrix(cor1$r)
    if(rm_xlab){
        lab1 = paste(rep("v",nrow(m1)),1:nrow(m1),sep="")
    }else{
        lab1 = rownames(m1)
    }
    if(rm_ylab){
        lab2 = paste(rep("v",ncol(m1)),1:ncol(m1),sep="")
    }else{
        lab2 = colnames(m1)
    }
    dimnames(m1) = list(lab1,lab2)
    if(reord){
        cl1 = hclust(x=cor1$r,method="average")
        m1 = m1[cl1$order,cl1$order]
    }
    pnt1 = as.matrix(cor1[[stat]])
    if(mirror){
        m1 = m1[,nrow(m1):1]
        pnt1 = pnt1[,nrow(pnt1):1]
    }
    n1 = 10
    col1 = colorRampPalette(c("Blue","White","Red"))(n1)
    for(i1 in c("eps","tiff")){
        if(i1 == "eps"){
            postscript(file=paste(f1,".eps",sep=""),width=10,height=10,
              colormodel="rgb",horizontal=FALSE,onefile=FALSE,paper="special",
              family="Helvetica")
        }else{
            tiff(file=paste(f1,".tif",sep=""),units="in",width=10,height=10,
              res=300,compression="lzw",family="sans")
        }
        trellis.par.set(list(fontsize=list(text=9,points=8)))
        opt = list(
          layout.heights=list(top.padding=1,bottom.padding=1),
          layout.widths=list(left.padding=1,right.padding=1),
          axis.text=list(cex=0.7),
          par.xlab.text=list(cex=0.7),
          par.ylab.text=list(cex=0.7),
          par.main.text=list(cex=0.8))
        print(levelplot(m1,at=seq(-1,1,len=n1),col.regions=col1,xlab="",
          ylab="",panel=plotpval,scales=list(y=list(tck=0,rot=0),x=list(tck=0,
          rot=90)),main=list(main),colorkey=list(tck=0,width=1.0),
          par.settings=opt))
        invisible(dev.off())
    }
}

# Perform Partition around medoids clustering (wide format) [allow for specific number of partitions]
# dst1 - distance matrix of variables to be analysed
# kbest - specify number of partitions [if is.null(kbest) then calculate best number of partitions]
# col1 - colors for points
# border - add border to columns
# f1 - output filename
# main - title for plots
plot_pam_v2 = function(dst1,kbest=NULL,col1=NULL,border=1,f1,main=NULL){
    library("cluster")
    kmax <- dim(as.matrix(dst1))[1]-1
    asw <- rep(0,kmax)
    for(k in 2:kmax)
        asw[k] = pam(dst1,k)$silinfo$avg.width
    if(is.null(kbest)){
        kbest = which.max(asw)
    }
    pres1 = pam(dst1,kbest)
    d1 = cbind(1:kmax,asw)
    colnames(d1) = c("K","Average silhouette width")
    d2 = silhouette(pres1)
    f2 = paste(f1,"_1.csv",sep="")
    f3 = paste(f1,"_2.csv",sep="")
    f4 = paste(f1,"_3.csv",sep="")
    write.table(d1,f2,quote=TRUE,sep=",",eol="\n",na="NA",dec=".",
      col.names=TRUE,row.names=FALSE)
    write.table(pres1$clusinfo,f3,quote=TRUE,sep=",",eol="\n",na="NA",dec=".",
      col.names=TRUE,row.names=FALSE)
    write.table(d2,f4,quote=TRUE,sep=",",eol="\n",na="NA",dec=".",col.names=NA)
    for(i1 in c("eps","tiff")){
        if(i1 == "eps"){
            postscript(file=paste(f1,".eps",sep=""),width=7.0,height=3.5,
              colormodel="rgb",horizontal=FALSE,onefile=FALSE,paper="special",
              family="Helvetica")
        }else{
            tiff(file=paste(f1,".tif",sep=""),units="in",width=7.0,height=3.5,
              res=300,compression="lzw",family="sans")
        }
        oldpar = par(mfcol=c(1,2),mar=c(2.2,2.2,1.2,1.2),mgp=c(1.0,0.2,0.0),
          tcl=-0.2,cex=0.8,cex.lab=0.8,cex.axis=0.8,cex.main=0.9)
        plot(d1,type="h",main=main)
        points(kbest,asw[kbest],col="red",type= "h")
        plot(d2,main=main,sub="",xlab="Silhouette width",do.n.k=FALSE,
          border=border,do.clus.stat=FALSE,mgp=c(1.0,0.2,0.0),
          col=col1[rownames(d2)])
        par(oldpar)
        invisible(dev.off())
    }
}

# min-max transformation
# y - vector
# miny - minimum value to scale y
# maxy - maximum value to scale y
minmax = function(y,miny=NULL,maxy=NULL){
    if(is.null(miny)){miny = min(y)}
    if(is.null(maxy)){maxy = max(y)}
    return((y - miny)/(maxy - miny))
}

# Inverse min-max transformation
# y - vector
# miny - minimum value to scale y
# maxy - maximum value to scale y
minmax_inv = function(y,miny,maxy){
    return(y*(maxy - miny) + miny)
}

# Construct network (wide format) [allow for very many nodes]
# dst1 - distance matrix of variables to be analysed
# col1 - colors for points
# vmode - network layout
# edge_thr - threshold for edges to plot
# edge_lwd - scale for edges
# label_cex - size of node labels
# f1 - output filename
# main - title for plots
plot_sna_v2 = function(dst1,col1=NULL,vmode="fruchtermanreingold",edge_thr=0.80,
  edge_lwd=10,label_cex=0.8,f1,main=NULL){
    library("sna")
    m1 = as.matrix(dst1)
    m1 = -m1                            #invert distances
    mind = min(m1)
    maxd = max(m1)
    m1 = minmax(m1,mind,maxd)           #min-max transformation
    m1[m1 < edge_thr] = NA              #remove dist < thresh
    m1 = ((m1 - edge_thr))/(1 - edge_thr) #rescale dist
    edge_col = "lightblue"
    vertex_cex=1.5
    for(i1 in c("eps","tiff")){
        if(i1 == "eps"){
            postscript(file=paste(f1,".eps",sep=""),width=7,height=10,
              colormodel="rgb",horizontal=FALSE,onefile=FALSE,paper="special",
              family="Helvetica")
        }else{
            tiff(file=paste(f1,".tif",sep=""),units="in",width=7,height=10,
              res=300,compression="lzw",family="sans")
        }
        oldpar = par(mar=c(1.1,1.1,2.1,1.1),mgp=c(0,0,0),xpd=TRUE,cex.main=0.9)
        gplot(m1,gmode="graph",displaylabels=TRUE,label.cex=label_cex,
          label.pos=5,mode=vmode,jitter=FALSE,main=paste(main,
          " [minmax(-dist) > ",edge_thr,"]",sep=""),vertex.cex=vertex_cex,
          vertex.col=col1[colnames(m1)],vertex.border=1,vertex.lty=1,
          vertices.last=TRUE,edge.lwd=edge_lwd,edge.lty=1,edge.lty.neg=1,
          edge.col=edge_col)
        par(oldpar)
        invisible(dev.off())
    }
}

# Summary statistics for network
# dst1 - distance matrix of variables to be analysed
# f1 - output filename
analyse_net = function(dst1,f1){
    library("sna")
    m1 = as.matrix(dst1)
    m1 = -m1                            #invert distances
    mind = min(m1)
    maxd = max(m1)
    m1 = (m1 - mind)/(maxd - mind)      #min-max transformation
    f2 = paste(f1,"_degree.csv",sep="")
    f3 = paste(f1,"_eigen.csv",sep="")
    f4 = paste(f1,"_dens.csv",sep="")
    #Degree - number of links incident upon a node
    net_deg = degree(m1,gmode="graph")
    s1 = order(net_deg,decreasing=TRUE)
    d1 = cbind(rownames(m1)[s1],net_deg[s1])
    colnames(d1) = c("Node","Degree")
    write.table(d1,f2,row.names=FALSE,col.names=TRUE,sep=",",dec=".",eol="\n",
      na="NA")
    #Eigenvector centrality - Influence of a node in a network
    net_evc = evcent(m1,gmode="graph",use.eigen=TRUE)
    s1 = order(abs(net_evc),decreasing=TRUE)
    d1 = cbind(rownames(m1)[s1],net_evc[s1])
    colnames(d1) = c("Node","Eigenvector")
    write.table(d1,f3,row.names=FALSE,col.names=TRUE,sep=",",dec=".",eol="\n",
      na="NA")
    #Density - the proportion of direct ties in a network relative to the total number possible
    nvert = nrow(m1)
    m2 = m1; diag(m2) = NA; m2[upper.tri(m2)] = NA
    nedge = sum(m2 != 0 & !is.na(m2))
    tedge = nties(m1,mode="graph")
    dens = gden(m1,mode="graph")
    d1 = data.frame(nvert,nedge,tedge,dens)
    write.table(d1,f4,row.names=FALSE,col.names=TRUE,sep=",",dec=".",eol="\n",
      na="NA")
}

# Calculate descriptive statistics (wide format)
# dat1 - data frame of variables to be analysed
# f1 - file name
describe_data_2 = function(dat1,f1){
    a1 = apply(dat1,2,function(y) mean(y,na.rm=TRUE))
    a2 = apply(dat1,2,function(y) length(y[!is.na(y)]))
    a3 = apply(dat1,2,function(y) sd(y,na.rm=TRUE))
    d1 = cbind(a1,a2,a3)
    colnames(d1) = c("mean","n","stdev")
    write.table(d1,f1,quote=TRUE,sep=",",eol="\n",na="NA",dec=".",col.names=NA)
}

# Calculate between groups variance [allow for only-NAs groups]
# dat1 - data frame of variables to be analysed
# grps - groups
calc_between_var_v2 = function(dat1,grps){
    vb = c()
    for(i1 in 1:ncol(dat1)){
         gmean = mean(dat1[,i1],na.rm=TRUE)
         ntot = 0
         ngrps = 0
         for(g1 in levels(grps)){
            v1 = dat1[grps==g1,i1]
            mean1 = mean(v1,na.rm=TRUE)
            if(!is.na(mean1)){
                len1 = length(v1[!is.na(v1)])
                ntot = ntot + len1*((mean1 - gmean)^2)
                ngrps = ngrps + 1
            }
         }
         n = max(ngrps - 1,0)
         vb = c(vb,ntot/n)
     }
     names(vb) = colnames(dat1)
     return(vb)
}

# Calculate within-group variance [allow for 1-element groups]
# dat1 - data frame of variables to be analysed
# grps - groups
calc_within_var_v2 = function(dat1,grps){
    vw = c()
    for(i1 in 1:ncol(dat1)){
        ntot = 0
        dtot = 0
        ngrps = 0
        for(g1 in levels(grps)){
            v1 = dat1[grps==g1,i1]
            sd1 = sd(v1,na.rm=TRUE)
            if(!is.na(sd1)){
                len1 = length(v1[!is.na(v1)])
                ntot = ntot + (len1 - 1)*(sd1*sd1)
                dtot = dtot + len1
                ngrps = ngrps + 1
            }
        }
        vw = c(vw,ntot/(dtot - ngrps))
    }
    names(vw) = colnames(dat1)
    return(vw)
}

# Calculate between groups covariance [allow for only-NAs groups]
# dat1 - data frame of variables to be analysed
# grps - groups
calc_between_cov_v2 = function(dat1,grps){
# Note: Relates two variables for different groups.
    vars = colnames(dat1)
    nvars = length(vars)
    m1 = matrix(NA,ncol=nvars,nrow=nvars,dimnames=list(vars,vars))
    d1 = 0
    for(i1 in 1:nvars){
        d1 = d1 + 1
        for(i2 in 1:nvars){
            d2 = d1 + (i2-i1)
            if(i1 > i2){
                v1 = dat1[,i1]
                v2 = dat1[,i2]
                gmean1 = mean(v1,na.rm=TRUE)
                gmean2 = mean(v2,na.rm=TRUE)
                covb = 0
                ngrps = 0
                if(!is.na(gmean1) & !is.na(gmean2)){ 
                    for(g1 in levels(grps)){
                        mean1 = mean(v1[grps==g1],na.rm=TRUE)
                        mean2 = mean(v2[grps==g1],na.rm=TRUE)
                        if(!is.na(mean1) & !is.na(mean2)){
                            len1 = length(v1[grps==g1 & !is.na(v1)])
                            covb = covb + (mean1-gmean1)*(mean2-gmean2)*(len1)
                            ngrps = ngrps + 1
                        }
                    }
                }
                n = max(ngrps - 1,0)
                m1[d1,d2] = covb/n
            }
        }
    }
    return(m1)
}

# Calculate within-group covariance [allow for only-NAs groups]
# dat1 - data frame of variables to be analysed
# grps - groups
calc_within_cov_v2 = function(dat1,grps){
# Note: Relates two variables for same-group individuals.
    vars = colnames(dat1)
    nvars = length(vars)
    m1 = matrix(NA,ncol=nvars,nrow=nvars,dimnames=list(vars,vars))
    d1 = 0
    for(i1 in 1:nvars){
        d1 = d1 + 1
        for(i2 in 1:nvars){
            d2 = d1 + (i2-i1)
            if(i1 > i2){
                v1 = dat1[,i1]
                v2 = dat1[,i2]
                glen1 = sum(!is.na(v1) & !is.na(v2))
                covw = 0
                ngrps = 0
                for(g1 in levels(grps)){
                    mean1 = mean(v1[grps==g1],na.rm=TRUE)
                    len1 = length(v1[grps==g1])
                    mean2 = mean(v2[grps==g1],na.rm=TRUE)
                    term1 = 0
                    if(!is.na(mean1) & !is.na(mean2)){
                        for(i3 in 1:len1){
                            if(!is.na(v1[grps==g1][i3])
                              & !is.na(v2[grps==g1][i3]))
                                term1 = term1 + ((v1[grps==g1][i3] - 
                                  mean1)*(v2[grps==g1][i3] - mean2))
                        }
                        ngrps = ngrps + 1
                    }
                    covw = covw + term1
                }
                m1[d1,d2] = covw/(glen1 - ngrps)
            }
        }
    }
    return(m1)
}

# Calculates top correlated variables
# dat1 - data frame of variables to be analysed
# method - type of correlation to use {"pearson","spearman"}
# k1 - number of top correlations to show
calc_highly_cor = function(dat1,method="spearman",k1){
    library("Hmisc")
    rres = rcorr(as.matrix(dat1),type=method)$r
    diag(rres) = NA
    rres[lower.tri(rres)] = NA
    a1 = as.data.frame(as.table(round(rres,3)))
    names(a1) <- c("v1", "v2","cor")
    return(head(a1[order(abs(a1$cor),decreasing=T),],n=k1))
}

# Calculate descriptive statistics for variables in grouped data (wide format) [allow for only-NAs and 1-element groups]
# dat1 - data frame of variables to be analysed
# grps - groups
describe_vars_v2 = function(dat1,grps,f1){
# Note: "separation" calculates how well a variable is able to discriminate groups (the greater the better).
#       "between-groups covariance" relates two variables for different groups.
#       "within-groups covariance" relates two variables for same-group individuals.
#       "correlation" measures association between variables across all groups.
    dat1 = dat1[,!apply(dat1,2,function(x) all(is.na(x)))] #remove variables without information
    vb = calc_between_var_v2(dat1,grps) #between groups variance
    vw = calc_within_var_v2(dat1,grps)  #within group variance
    sep = vb/vw                         #separation
    cb = calc_between_cov_v2(dat1,grps) #between groups covariance
    cw = calc_within_cov_v2(dat1,grps)  #within group covariance
    a1 = calc_highly_cor(dat1,k1=50)    #correlation between variables accross groups
    f2 = paste(f1,"_vb.csv",sep="")
    f3 = paste(f1,"_vw.csv",sep="")
    f4 = paste(f1,"_sep.csv",sep="")
    f5 = paste(f1,"_cb.csv",sep="")
    f6 = paste(f1,"_cw.csv",sep="")
    f7 = paste(f1,"_cor.csv",sep="")
    write.table(vb,f2,quote=TRUE,sep=",",eol="\n",na="NA",dec=".",col.names=NA)
    write.table(vw,f3,quote=TRUE,sep=",",eol="\n",na="NA",dec=".",col.names=NA)
    write.table(sep,f4,quote=TRUE,sep=",",eol="\n",na="NA",dec=".",col.names=NA)
    write.table(cb,f5,quote=TRUE,sep=",",eol="\n",na="NA",dec=".",col.names=NA)
    write.table(cw,f6,quote=TRUE,sep=",",eol="\n",na="NA",dec=".",col.names=NA)
    write.table(a1,f7,quote=TRUE,sep=",",eol="\n",na="NA",dec=".",col.names=NA)
}

# Perform pairwise t-test in grouped data (wide format) [select significant tests]
# dat1 - data frame of variables to be analysed
# grps - groups
# alternative - alternative hypothesis {"two.sided","greater","less"}
# paired - paired t-test
# var.equl - pooled variance
# p.adjust - correction method {"holm","hochberg","hommel","bonferroni","BH","BY","fdr","none"}
# stat - statistic to report {"stat","pval","adjp"}
# thr - threshold of significance to report
# f1 - output filename
calc_pairwise_ttest_v2 = function(dat1,grps,alternative="two.sided",
  paired=FALSE,var.equal=TRUE,p.adjust="BH",stat="adjp",thr=0.05,f1){
    vars = colnames(dat1)
    nvars = length(vars)
    lgrps = levels(grps)
    ngrps = length(lgrps)
    ntests = choose(ngrps,2)
    tests = c()
    for(i1 in 1:ngrps){
        for(i2 in 1:ngrps){
            if(i1 < i2){
                tests = c(tests,paste("X",i1,"vsX",i2,sep=""))
            }
        }
    }
    res1 = matrix(NA,ncol=ntests,nrow=nvars,dimnames=list(vars,tests))
    res2 = matrix(NA,ncol=ntests,nrow=nvars,dimnames=list(vars,tests))
    res3 = matrix(NA,ncol=ntests,nrow=nvars,dimnames=list(vars,tests))
    res4 = matrix(NA,ncol=ntests,nrow=nvars,dimnames=list(vars,tests))
    k1 = 0
    for(i1 in 1:nvars){
        k1 = k1 + 1
        k2 = 0
        for(j1 in 1:ngrps){
            for(j2 in 1:ngrps){
                if(j1 < j2){
                    k2 = k2 + 1
                    v1 = dat1[grps==lgrps[j1],i1]
                    v2 = dat1[grps==lgrps[j2],i1]
                    a1 = v1[!is.na(v1)]
                    a2 = v2[!is.na(v2)]
                    if(paired){
                        a1 = nrow(unique(na.omit(cbind(v1,v2))))
                        skip = a1 <= 1
                    }else{
                        a1 = length(unique(na.omit(v1)))
                        a2 = length(unique(na.omit(v2)))
                        skip = a1 == 0 | a2 == 0 | (a1 == 1 & a2 == 1) 
                    }
                    if(!skip){
                        tres = t.test(v1,v2,alternative=alternative,paired=paired,
                          var.equal=var.equal,na.action="na.omit")
                        res1[k1,k2] = tres$statistic
                        res2[k1,k2] = tres$parameter
                        res3[k1,k2] = tres$p.value
                    }
                }
            }
        }
        res4[k1,] = p.adjust(res3[k1,],method=p.adjust)
    }
    if(stat == "stat"){
        res5 = res1
    }else if(stat == "pval"){
        res5 = res3
    }else if(stat == "adjp"){
        res5 = res4
    }else{
        res5 = matrix(NA,nrow=nrow(res1),ncol=ncol(res1),
          dimnames=list(rownames(res1),colnames(res1)))
    }
    s1 = apply(res5,1,function(x){any(!is.na(x) & x < 0.05)})
    res5 = res5[s1,]
    f2 = paste(f1,"_stat.csv",sep="")
    f3 = paste(f1,"_df.csv",sep="")
    f4 = paste(f1,"_pval.csv",sep="")
    f5 = paste(f1,"_adjp.csv",sep="")
    f6 = paste(f1,"_sign.csv",sep="")
    write.table(res1,file=f2,sep=",",row.names=TRUE,col.names=NA)
    write.table(res2,file=f3,sep=",",row.names=TRUE,col.names=NA)
    write.table(res3,file=f4,sep=",",row.names=TRUE,col.names=NA)
    write.table(res4,file=f5,sep=",",row.names=TRUE,col.names=NA)
    write.table(res5,file=f6,sep=",",row.names=TRUE,col.names=NA)
}

# Plot value~grps (wide format) [specify labels for y-axis]
# dat1 - data frame of variables to be analysed
# grps - groups
# col1 - colors associated to subjects
# units1 - labels for y-axis
# k1 - number of plots per row to print
# k2 - number of plots per columns to print
# rm_xlab - remove x label information
# rm_ylab - remove y label information
# rm_main - remove main label information
# horiz - plot as landscape
# f1 - output filename
boxplot_data_v2 = function(dat1,grps,col1=NULL,units1,k1=8,k2=4,rm_xlab=TRUE,
  rm_ylab=TRUE,rm_main=TRUE,horiz=FALSE,f1){
    nplot = k1*k2
    ngrps = length(levels(grps))
    if(is.null(col1)){
        col1 = rainbow(ngrps)
    }
    if(horiz){
        iwidth = 10
        iheight = 7
    }else{
        iwidth = 7
        iheight = 10
    }
    for(i1 in c("eps","tiff")){
        f2 = paste(f1,"_1",sep="")
        if(i1 == "eps"){
            postscript(file=paste(f2,".eps",sep=""),width=iwidth,height=iheight,
              colormodel="rgb",horizontal=FALSE,onefile=FALSE,paper="special",
              family="Helvetica")
        }else{
            tiff(file=paste(f2,".tif",sep=""),units="in",width=iwidth,
              height=iheight,res=300,compression="lzw",family="sans")
        }
        if(horiz){
            oldpar = par(mfrow=c(k2,k1),mar=c(1.7,2.3,1.0,0.2),mgp=c(1.0,0.2,0),
              tcl=-0.2,cex=0.8,cex.lab=0.8,cex.axis=0.8,cex.main=0.9)
        }else{
            oldpar = par(mfrow=c(k1,k2),mar=c(1.7,2.3,1.0,0.2),mgp=c(1.0,0.2,0),
              tcl=-0.2,cex=0.8,cex.lab=0.8,cex.axis=0.8,cex.main=0.9)
        }
        iplot = 0
        for(i2 in 1:ncol(dat1)){
            iplot = iplot + 1
            if(rm_xlab){
                xlab1 = ""
            }else{
                xlab1 = colnames(dat1)[i2]
            }
            if(rm_ylab){
                ylab1 = ""
                yaxt1 = "n"
            }else{
                ylab1 = units1[i2]
                yaxt1 = "s"
            }
            if(rm_main){
                main1 = paste("v",i2,sep="")
            }else{
                main1 = colnames(dat1)[i2]
            }
            boxplot(dat1[,i2]~grps,main=main1,xlab=xlab1,ylab=ylab1,xaxt="n",
              yaxt=yaxt1,col=col1)
            if(iplot%%nplot==0 & iplot < ncol(dat1)){
                invisible(dev.off())
                f2 = paste(f1,"_",1+iplot/nplot,sep="")
                if(i1 == "eps"){
                    postscript(file=paste(f2,".eps",sep=""),width=iwidth,
                      height=iheight,colormodel="rgb",horizontal=FALSE,
                      onefile=FALSE,paper="special",family="Helvetica")
                }else{
                    tiff(file=paste(f2,".tif",sep=""),units="in",width=iwidth,
                      height=iheight,res=300,compression = "lzw",family="sans")
                }
                if(horiz){
                    oldpar = par(mfrow=c(k2,k1),mar=c(1.7,2.3,1.0,0.2),
                      mgp=c(1.0,0.2,0),tcl=-0.2,cex=0.8,cex.lab=0.8,
                      cex.axis=0.8,cex.main=0.9)
                }else{
                    oldpar = par(mfrow=c(k1,k2),mar=c(1.7,2.3,1.0,0.2),
                      mgp=c(1.0,0.2,0),tcl=-0.2,cex=0.8,cex.lab=0.8,
                      cex.axis=0.8,cex.main=0.9)
                }
            }
        }
        invisible(dev.off())
    }
}

# Remove missing data
# dat1 - data frame of variables to be analysed
# init - parameter for initial threshold 2^-i where i={init, init+1, ...}
# col1st - analyse columns firts
# plotit - plot analysis information
# verbose - print analysis information
remove_na = function(dat1,init,col1st=TRUE,plotit=TRUE,verbose=FALSE){
    res = dat1
    n1 = nrow(res)
    n2 = ncol(res)
    srow = apply(res,1,function(x){sum(is.na(x))})
    scol = apply(res,2,function(x){sum(is.na(x))})
    if(plotit){
        par(mfrow=c(2,1),ask=TRUE)
    }
    while(sum(c(srow,scol)) != 0){
        thrs = 2^-init
        #by columns
        if(col1st){
            scol = apply(res,2,function(x){sum(is.na(x))})
            res = res[,scol/n1 < thrs] #remove columns with is.na() > 2^-init
            if(verbose | plotit){
                tcol = table(scol)
                names(tcol) = round(100*as.numeric(names(tcol))/n1,1)
            }
            if(verbose){
                print(tcol)
            }
            if(plotit){
                main1 = paste("thrs = ",round(100*thrs,1),"%",sep="")
                plot(tcol,main=main1,ylab="Freq(cols)")
                abline(v=100*thrs,col="red",lty=2)
            }
        }
        #by row
        srow = apply(res,1,function(x){sum(is.na(x))})
        res = res[srow/n2 < thrs,]
        if(verbose | plotit){
            trow = table(srow)
            names(trow) = round(100*as.numeric(names(trow))/n2,1)
        }
        if(verbose){
            print(trow)
        }
        if(plotit){
            main1 = paste("thrs = ",round(100*thrs,1),"%",sep="")
            plot(trow,main=main1,ylab="Freq(rows)")
            abline(v=100*thrs,col="red",lty=2)
        }
        #by columns
        if(!col1st){
            thrs = 2^-init
            scol = apply(res,2,function(x){sum(is.na(x))})
            res = res[,scol/n1 < thrs]
            if(verbose | plotit){
                tcol = table(scol)
                names(tcol) = round(100*as.numeric(names(tcol))/n1,1)
            }
            if(verbose){
                print(tcol)
            }
            if(plotit){
                main1 = paste("thrs = ",round(100*thrs,1),"%",sep="")
                plot(tcol,main=main1,ylab="Freq(cols)")
                abline(v=100*thrs,col="red",lty=2)
            }
        }
        init = init + 1
    }
    return(res)
}

# Reduce number of predictors
# dat1 - data frame of variables to be analysed (only predictors)
# y - response variable
# thr1 - upper threshold for correlation between predictors 
# thr2 - lower threshold for correlation between predictors and response variables [if thr2==NULL then reduces till maxsize]
# thr3 - lower threshold for coefficient of variation [if thr3==NULL then reduces till maxsize]
# thr4 - upper threshold for Variance Inflation Factor (VIF)
# maxsize - maximum number of predictors [if maxsize==NULL then reduces till nsubj - 1]
# f1 - output filename
reduce_predictors_v2 = function(dat1,y,thr1=0.95,method="spearman",thr2=NULL,
  thr3=NULL,thr4=10,maxsize=NULL,f1){
# Note: Sqrt(R^2) for continuous variables is the correlation between variables. Sqrt(R^2)
# for continuous and categorical variables is the correlation between observed and predicted
# values.
    if(is.null(maxsize)){
        maxsize = nrow(dat1) - 1
    }
    if(is.null(thr2)){
        doStep2 = TRUE
    }else if(thr2 > 0){
        doStep2 = TRUE
    }else{
        doStep2 = FALSE
    }
    if(is.null(thr3)){
        doStep3 = TRUE
    }else if(thr3 > 0){
        doStep3 = TRUE
    }else{
        doStep3 = FALSE
    }
    write("Step1: Highly correlated between themselves",file=f1,append=FALSE)
    if(thr1 < 1){                       #remove variables highly correlated between themselves
        vars = colnames(dat1)
        rm_vars = c()
        for(i1 in vars){
            for(i2 in vars){
                if(i1 < i2 & !i1 %in% rm_vars){
                    rres = cor(dat1[,i1],dat1[,i2],method=method,use="complete")
                    if(abs(rres) > thr1){
                        write(paste("  Remove ",i2," (r_",i1," = ",round(rres,2),
                          ")",sep=""),file=f1,append=TRUE)
                        rm_vars = c(rm_vars,i2)
                    }
                }
            }
        }
        dat1 = dat1[,!vars %in% rm_vars]
    }else{
        write("<Skip step>",file=f1,append=TRUE)
    }
    write("Step2: Low correlated with response variable",file=f1,append=TRUE)
    if(doStep2){                        #remove variables with low correlation with dependent variable
        vars = colnames(dat1)
        if(is.factor(y)){
            y = as.character(y)
        }
        cors = apply(dat1,2,function(x){sqrt(summary(lm(y~x))$r.squared)})
        if(is.null(thr2)){
            n_rm = ncol(dat1) - maxsize #n_rm, so that, nvars = nsubj - 1
            if(n_rm <= 0){
                thr2 = 0.00
            }else{
                thr2 = sort(cors)[n_rm]
            }
        }
        rm_vars = cors <= thr2
        if(sum(rm_vars) == 0){
            write("<None removed>",file=f1,append=TRUE)
        }else{
            dat1 = dat1[,!rm_vars]
            tmp = cbind(vars[rm_vars],cors[rm_vars])
            tmp = apply(tmp,1,function(x){paste("  Remove ",x[1]," (r = ",
              round(as.numeric(x[2]),2),")",sep="")})
            write(tmp,file=f1,append=TRUE)
        }
    }else{
        write("<Skip step>",file=f1,append=TRUE)
    }
    write("Step3: Variables with low variation",file=f1,append=TRUE)
    if(doStep3){                        #remove variables with low variation
        vars = colnames(dat1)
        sds = apply(dat1,2,sd,na.rm=TRUE)
        means = apply(dat1,2,mean,na.rm=TRUE)
        if(is.null(thr3)){
            n_rm = ncol(dat1) - maxsize #n_rm, so that, nvars = nsubj - 1
            if(n_rm <= 0){
                thr3 = 0.00
            }else{
                thr3 = sort(sds/means)[n_rm]
            }
        }
        rm_vars = sds/means <= thr3
        if(sum(rm_vars) == 0){
            write("<None removed>",file=f1,append=TRUE)
        }else{
            dat1 = dat1[,!rm_vars]
            tmp = cbind(vars[rm_vars],sds/means[rm_vars])
            tmp = apply(tmp,1,function(x){paste("  Remove ",x[1]," (r = ",
              round(as.numeric(x[2]),2),")",sep="")})
            write(tmp,file=f1,append=TRUE)
        }
    }else{
        write("<Skip step>",file=f1,append=TRUE)
    }
    write("Step4: Variables with low VIF",file=f1,append=TRUE)
    if(thr4 < Inf){                     #remove variables using VIF (reduce multicollinearity)
        vars = colnames(dat1)
        keep_vars = stepwise_vif(dat1,thresh=thr4,trace=FALSE)
        if(sum(!keep_vars) == 0){
            write("<None removed>",file=f1,append=TRUE)
        }else{
            dat1 = dat1[,keep_vars]
            tmp = paste("  Remove",vars[!vars %in% keep_vars])
            write(tmp,file=f1,append=TRUE)
        }
    }else{
        write("<Skip step>",file=f1,append=TRUE)
    }
    return(dat1)
}
 
# Stepwise Variance Inflation Factors (VIF) selection to reduce collinearity
# [From "Marcus W Beck" in https://gist.github.com/fawda123/4717702#file-vif_fun-r]
# in_frame - data frame of variables to be analyzed
# thresh - threshold for VIF
# trace - print output of each iteration
stepwise_vif = function(in_frame,thresh=10,trace=TRUE){
    if(class(in_frame) != 'data.frame'){
        in_frame = data.frame(in_frame)
    }
    #get initial vif value for all comparisons of variables
    var_names = names(in_frame)
    vif_init = vif(in_frame)
    vif_max = max(as.numeric(vif_init),na.rm=TRUE)
    if(vif_max < thresh){
        if(trace==TRUE){
            print(data.frame(vif_init))
            cat('\n')
            cat(paste('All variables have VIF < ', thresh,', max VIF ',
              round(vif_max,2), sep=''),'\n\n')
        }
        return(var_names)
    }else{
        in_dat = in_frame
        #backwards selection of explanatory variables, stops when all VIF values are below 'thresh'
        while(vif_max >= thresh){
            var_names = names(in_dat)
            vif_vals = vif(in_dat)
            imax = which(vif_vals == max(as.numeric(vif_vals),na.rm=TRUE))[1]
            vif_max = as.numeric(vif_vals[imax])
            if(vif_max < thresh){
                break
            }
            if(trace==TRUE){
                print(data.frame(vif_vals))
                cat('\n')
                cat('removed: ',names(vif_vals)[imax],vif_max,'\n\n')
                flush.console()
            }
            in_dat = in_dat[,!names(in_dat) %in% names(vif_vals)[imax]]
        }
        return(names(in_dat))
    }
}

# Calculate Variance Inflation Factors (VIF) at a model-level
# [based on HH::vif()]
# dat1 - data frame of variables to be analysed (only predictors)
vif = function(dat1){
#see https://en.wikipedia.org/wiki/Variance_inflation_factor
#see car::vif() for calculation of VIF at variable-level
#see HH::vif() for calculation of VIF at model-level
#see faraway::vif() for calculation of VIF at model-level
#see fmsb::vif() for calculation of VIF at model-level
# Note: VIF can be defined at variable-level (using the partial R^2) and at the model
# level (using R^2, i.e. coefficient of determination).
    vars = colnames(dat1)
    nvars = length(vars)
    r2 = vector("numeric",length=nvars)
    names(r2) = vars
    for(i1 in 1:nvars){
        tmp = lm(dat1[,i1] ~ data.matrix(dat1[,-i1]),na.action="na.omit")
        r2[i1] = 1/(1 - summary(tmp)$r.squared)
    }
    return(r2)
}

#Register "glmulti::getfit()" for "multinom" objects
#[adapted from ("Using glmulti...",Vincent Calcagno,February 18 2012)]
library("nnet")
library("glmulti")
setOldClass("multinom")
setMethod('getfit','multinom',function(object,...){
    coefs = summary(object)$coefficients
    namez = dimnames(coefs)
    neonamez = unlist(lapply(namez[[2]],function(x){lapply(namez[[1]],
      function(y){neonamez = paste(x,y,sep="/")})}))
    ses = summary(object)$standard.errors
    neocoefs = as.vector(coefs)
    neoses = as.vector(ses)
    neodf = rep(summary(object)$edf,length(coefs))
    return(data.frame(estimate=neocoefs,se=neoses,df=neodf,row.names=neonamez))
})

# Automatic model selection for binomial logistic regression
# dat1 - data frame of variables to be analysed
# grps - groups
# col1 - colors for each group
# thr1 - upper threshold for correlation between predictors 
# thr2 - lower threshold for correlation between predictors and response variables [if thr2==NULL then reduces till maxsize]
# thr3 - lower threshold for coefficient of variation [if thr3==NULL then reduces till maxsize]
# thr4 - upper threshold for Variance Inflation Factor (VIF)
# maxsize - maximum number of candidate predictors [if maxsize==NULL then reduces till nsubj - 1]
# delrows - obtain complete cases by eliminating subjects
# maxs1 -  maximum number of terms in candidate models {"hard","soft1","soft2","soft3","soft4","soft5","soft6","exhaustive"}
# confs - size of the returned confidence set 
# level - level models' complexity {1: only main effects; 2: pairwise interactions}
# crit - information criteria to use {"aic","aicc","bic","qaic","qaicc"}
# thrs1 - threshold for number of candidate predictors to calculate nmods exactly
# thrs2 - threshold for nmods to perform exhaustive search
# nreps - number of replicas for genetic algorithm
# f1 - output filename
mod_select_binom = function(dat1,grps,col1=NULL,thr1=0.95,thr2=0.50,
  thr3=NULL,thr4=Inf,maxsize=30,delrows=FALSE,maxs1="exhaustive",confs=100,
  level=1,crit="aicc",thrs1=20,thrs2=100000,nreps=2,f1){
    library("glmulti")
    library("MASS")
    dat1 = dat1[,!apply(dat1,2,function(x) all(is.na(x)))] #remove variables with all NA
    if(delrows){
        ina = apply(dat1,1,function(x){any(is.na(x))})
        dat1 = dat1[!ina,]                  #remove all subjects with NAs
        grps = grps[!ina]                   #remove all subjects with NAs
        grps = factor(grps,levels=sort(unique(as.character(grps)))) #in case some level was lost
    }else{
        ina = apply(dat1,2,function(x){any(is.na(x))})
        dat1 = dat1[,!ina]                  #remove all variables with NAs
    }
    if(thr1 != 1.00 | thr2 != 0.00 | thr3 != 0.00 | thr4 != Inf){
        f2 = paste(f1,"_datred.log",sep="")
        dat1 = reduce_predictors_v2(dat1,grps,thr1=thr1,thr2=thr2,thr3=thr3,
          thr4=thr4,maxsize=maxsize,f1=f2)  #reduce number of candidate predictors
    }
    d1 = data.frame(grps,dat1)
    colnames(d1) = sapply(colnames(d1),function(x){gsub("[-|.]","_",x)}) #remove symbols used in "formulas"
    nvars = ncol(d1)
    vars = colnames(d1)
    a1 = ifelse(crit=="aicc",3,2)
    hmax = nrow(d1) - 1 - a1            #hard max (i.e. nparams = nsubj - [2 or 3])
    if(is.numeric(maxs1)){
        maxs1 = maxs1
    }else if(maxs1 == "hard"){          #hard max nterms
        maxs1 = hmax
    }else if(maxs1 == "soft1"){         #soft max nterms
        maxs1 = nrow(d1) - 10
    }else if(maxs1 == "soft2"){         #soft max nterms
        maxs1 = nrow(d1) - 20
    }else if(maxs1 == "soft3"){         #soft max nterms
        maxs1 = nrow(d1) - 50
    }else if(maxs1 == "soft4"){         #soft max nterms
        maxs1 = nrow(d1) - 100
    }else if(maxs1 == "soft5"){         #soft max nterms
        maxs1 = nrow(d1) - 200
    }else if(maxs1 == "soft6"){         #soft max nterms
        maxs1 = nrow(d1) - 400
    }else if(maxs1 == "exhaustive"){    #max nterms so that nmods < thrs1 (i.e. force exhaustive search)
        maxs1 = 0
        tmp = 0
        while(tmp < thrs2){
            maxs1 = maxs1 + 1
            tmp = sum(sapply(1:maxs1,function(x){choose(nvars-1,x)})) + 1
        }
        maxs1 = maxs1 - 1
    }
    maxs1 = min(maxs1,hmax)             #prevent maxs1 > hard max
    maxs1 = max(maxs1,1)                #prevent maxs1 <= 0
    if(nvars - 1 <= thrs1){             #accounts for possible filters (slow)
        nmods = glmulti(y=vars[1],xr=vars[2:nvars],data=d1,level=1,method="d",
          maxsize=maxs1,report=FALSE,fitfunction="multinom",trace=FALSE)
    }else{                              #does not account for possible filters
        nmods = sum(sapply(1:maxs1,function(x){choose(nvars-1,x)})) + 1
    }
    if(nmods < thrs2){                  #exhaustive search
        res3 = glmulti(y=vars[1],xr=vars[2:nvars],data=d1,level=1,maxsize=maxs1,
          method="h",crit=crit,confsetsize=confs,plotty=FALSE,report=FALSE,
          includeobjects=TRUE,marginality=FALSE,model=TRUE,
          fitfunction="glm",family=binomial)
    }else{                              #genetic algorithm
        res1 = vector("list",length=nreps)
        res2 = vector("list",length=nreps)
        for(i1 in 1:nreps){
            f2 = paste(f1,"_rep",i1,"_1",sep="")
            f3 = paste(f1,"_rep",i1,"_2",sep="")
            res1[[i1]] = glmulti(y=vars[1],xr=vars[2:nvars],data=d1,level=1,
              maxsize=maxs1,method="g",crit=crit,confsetsize=confs,plotty=FALSE,
              report=FALSE,name=f2,includeobjects=TRUE,marginality=FALSE,
              model=TRUE,fitfunction="glm",family=binomial)
            popsize = 200               #Population size (default = 100)
            mutrate = 0.01              #Per locus (i.e. per term) mutation rate, [0,1], (default = 10^-3)
            sexrate = 0.2               #Sexual reproduction rate, [0,1], (default = 0.1)
            imm = 0.6                   #Immigration rate, [0,1] (default = 0.3)
            deltaM = 0.005              #Stop Rule: change in mean IC (default = 0.05)
            deltaB = 0.005              #Stop Rule: change in best IC (default = 0.05)
            conseq = 10                 #Stop Rule: times with no improvement (default = 5)
            res2[[i1]] = glmulti(y=vars[1],xr=vars[2:nvars],data=d1,level=1,
              maxsize=maxs1,method="r",crit=crit,confsetsize=confs,plotty=FALSE,
              report=FALSE,name=f3,includeobjects=TRUE,marginality=FALSE,
              model=TRUE,popsize=popsize,mutrate=mutrate,sexrate=sexrate,
              imm=imm,deltaM=deltaM,deltaB=deltaB,conseq=conseq,resumefile=f2,
              fitfunction="glm",family=binomial)
        }
        if(nreps > 1){
            res3 = glmulti::consensus(xs=res2,confsetsize=confs)
        }else{
            res3 = res2
        }
    }
    t1 = weightable(res3)
    t2 = t1[t1[,crit] <= min(t1[,crit]) + 2,]
    t3 = coef(res3)
    f2 = paste(f1,".log",sep="")
    f3 = paste(f1,"_1.csv",sep="")
    f4 = paste(f1,"_2.csv",sep="")
    capture.output(print(res3),file=f2,append=FALSE,type="output")
    write.table(t2,file=f3,sep=",",row.names=TRUE,col.names=NA)
    write.table(t3,file=f4,sep=",",row.names=TRUE,col.names=NA)
    for(i1 in c("eps","tiff")){
        if(i1 == "eps"){
            postscript(file=paste(f1,".eps",sep=""),width=7.0,height=10,
              colormodel="rgb",horizontal=FALSE,onefile=FALSE,paper="special",
              family="Helvetica")
        }else{
            tiff(file=paste(f1,".tif",sep=""),units="in",width=7.0,height=10,
              res=300,compression="lzw",family="sans")
        }
        oldpar = par(mfcol=c(3,1),mar=c(2.2,2.2,1.2,1.2),mgp=c(1.0,0.2,0.0),
          oma=c(0,0,0,0),tcl=-0.2,cex=0.8,cex.lab=0.8,cex.axis=0.8,cex.main=0.9)
        par(mar=c(2.2,10.2,1.2,1.2))
        plot(res3,type="s",cex.names=0.6)
        par(mar=c(2.2,2.2,1.2,1.2))
        plot(res3,type="p")
        plot(res3,type="w")
        par(oldpar)
        invisible(dev.off())
    }
    a1 = formula(res3@objects[[1]]$terms) #formula
    a2 = res3@objects[[1]]$model        #data
    l1 = list(a1,a2)
    names(l1) = c("formula","data")
    return(l1)
}

# Fit binomial logistic regression using GLM [account for very many levels]
# form1 - formula of the model to fit
# dat1 - data frame of variables to fit 
# col1 - colors for each group
# f1 - output filename
# main - title for plots
fit_binomial_v2 = function(form1,dat1,col1=NULL,f1,main){
#see https://en.wikipedia.org/wiki/Logistic_regression
# Note: Logistic regression assumptions:
# 1. The error terms need to be independent.
# 2. The independent variables should be independent from each other (i.e. little or no multicollinearity).
# 3. Linearity between independent variables and log-odds.
# 4. Large sample sizes (i.e. 10-30 cases per parameter).
    library("MASS")
    library("gtools")
    library("car")
    fit1 = glm(form1,family=binomial,data=dat1,model=TRUE)
    vars = all.vars(form1)
    grps = dat1[,vars[1]]               #response variable
    levs = levels(grps)                 #levels of response variable
    nlev = length(levs)                 #number of levels of response variable
    sfit1 = summary(fit1)
    form0 = formula(paste(vars[1],"~1",sep=""))
    fit0 = glm(form0,family=binomial,data=dat1,model=TRUE)
    r1 = 1 - logLik(fit1)/logLik(fit0)
    pval1 = anova(fit1,fit0)$"Pr(Chi)"[2]
    pred1 = predict(fit1,dat1)
    r2 = sum(pred1 == grps)/length(grps)
    t1 = sfit1$coefficients
    t2 = cbind(grps,1-fit1$fit,fit1$fit)
    colnames(t2) = c("Cluster",paste(rep("Prob(",nlev),gsub("$",")",levs),sep=""))
    t3 = Anova(fit1,type="II")
    f2 = paste(f1,".log",sep="")
    f3 = paste(f1,"_1.csv",sep="")
    f4 = paste(f1,"_2.csv",sep="")
    f5 = paste(f1,"_3.csv",sep="")
    capture.output(print(sfit1),file=f2,append=FALSE,type="output")
    write(paste("R McFadden = ",r1,"(",pval1,")",sep=""),file=f2,append=TRUE)
    write(paste("R count = ",r2,sep=""),file=f2,append=TRUE)
    write.table(t1,file=f3,sep=",",row.names=TRUE,col.names=NA)
    write.table(t2,file=f4,sep=",",row.names=TRUE,col.names=NA)
    write.table(t3,file=f5,sep=",",row.names=TRUE,col.names=NA)
    minf = 0-1e-12/2
    maxf = 1+1e-12/2
    logitfit = (round(fit1$fit,12) - minf)/(maxf - minf) #avoid exact 0s and 1s
    if(is.null(col1)){
        col1 = rainbow(nlev)
    }
    col2 = col1[as.numeric(grps)]       #observed group
    if(nrow(dat1) > 100){
        width1 = 7.0
        height1 = 10.0
        mfrow1 = c(2,1)
    }else{
        width1 = 7.0
        height1 = 3.5
        mfrow1 = c(1,2)
    }
    if(nlev > 15){
        ncol1 = ceiling(nlev/15)
    }else{
        ncol1 = 1
    }
    f2 = paste(f1,"_1",sep="")
    for(i1 in c("eps","tiff")){
        if(i1 == "eps"){
            postscript(file=paste(f2,".eps",sep=""),width=width1,height=height1,
              colormodel="rgb",horizontal=FALSE,onefile=FALSE,paper="special",
              family="Helvetica")
        }else{
            tiff(file=paste(f2,".tif",sep=""),units="in",width=width1,
              height=height1,res=300,compression="lzw",family="sans")
        }
        oldpar = par(mfrow=mfrow1,mar=c(2.2,2.2,1.2,1.2),mgp=c(1.0,0.2,0.0),
          tcl=-0.2,cex=0.8,cex.lab=0.8,cex.axis=0.8,cex.main=0.9)
        plot(fit1$residuals~logitfit,xlab="Predicted values",ylab="Residuals",
          main="Residuals vs Fitted",col=replicate(3,col2),cex=0.6,pch=23)
        legend("bottomleft",inset=c(0.02,0.02),fill=col1,legend=levs,cex=0.5,
          ncol=ncol1)
        abline(h=0,col=8,lty=3)
        barplot(rbind(1-fit1$fit,fit1$fit),col=col1,las=2,ylab="Probability",
          main=main,cex.names=0.6)
        legend("topright",inset=c(0.01,0.01),fill=col1,legend=levs,bg="white",
          cex=0.5,ncol=ncol1)
        par(oldpar)
        invisible(dev.off())
    }
    vars = all.vars(fit1$terms)         #predictors (including "intercept")
    nvar = length(vars)                 #number of predictors (including "intercept")
    y1 = fit1$fit                       #probability
    y2 = fit1$residuals                 #residuals
    pos1 = rep(1,nrow(fit1$model))
    pos1[y1 < 0.5] = 3
    pos2 = rep(1,nrow(fit1$model))
    pos2[y2 == min(y2)] = 3
    col3 = col1[apply(cbind(1-fit1$fit,fit1$fit),1,which.max)]  #predicted group
    nplot = 3
    for(i1 in c("eps","tiff")){
        f2 = paste(f1,"_2",sep="")
        if(i1 == "eps"){
            postscript(file=paste(f2,".eps",sep=""),width=7.0,height=10,
              colormodel="rgb",horizontal=FALSE,onefile=FALSE,paper="special",
              family="Helvetica")
        }else{
            tiff(file=paste(f2,".tif",sep=""),units="in",width=7.0,height=10,
              res=300,compression="lzw",family="sans")
        }
        oldpar = par(mfrow=c(3,2),mar=c(2.2,2.2,1.2,1.2),mgp=c(1.0,0.2,0.0),
          tcl=-0.2,cex=0.8,cex.lab=0.8,cex.axis=0.8,cex.main=0.9)
        iplot = 0
        for(i2 in 2:nvar){
            iplot = iplot + 1
            xlab1 = vars[i2]
            x = fit1$model[,xlab1]
            plot(y1~x,col=col3,ylab="Fitted value",ylim=c(0,1),xlab=xlab1)
            text(y1~x,labels=rownames(fit1$model),cex=0.6,pos=pos1,col=col3)
            abline(h=1/nlev,lty=2)
            legend("bottomright",inset=c(0.02,0.02),fill=col1,legend=levs,
              cex=0.8)
            plot(y2~x,col=col3,ylab="Residuals",xlab=xlab1)
            text(y2~x,labels=rownames(fit1$model),cex=0.6,pos=pos2,col=col3)
            abline(h=0,lty=2)
            legend("bottomright",inset=c(0.02,0.02),fill=col1,legend=levs,
              cex=0.5,ncol=ncol1)
            if(iplot%%nplot==0 & iplot < nvar - 1){
                par(oldpar)
                invisible(dev.off())
                f2 = paste(f1,2+iplot/nplot,sep="_")
                if(i1 == "eps"){
                    postscript(file=paste(f2,".eps",sep=""),width=7.0,
                      height=10.0,colormodel="rgb",horizontal=FALSE,
                      onefile=FALSE,paper="special",family="Helvetica")
                }else{
                    tiff(file=paste(f2,".tif",sep=""),units="in",width=7.0,
                      height=10.0,res=300,compression = "lzw",family="sans")
                }
                oldpar = par(mfrow=c(3,2),mar=c(2.2,2.2,1.2,1.2),mgp=c(1.0,0.2,
                  0.0),tcl=-0.2,cex=0.8,cex.lab=0.8,cex.axis=0.8,cex.main=0.9)
            }
        }
        par(oldpar)
        invisible(dev.off())
    }
}

# Automatic model selection for multinomial logistic regression
# dat1 - data frame of variables to be analysed
# grps - groups
# thr1 - upper threshold for correlation between predictors 
# thr2 - lower threshold for correlation between predictors and response variables [if thr2==NULL then reduces till maxsize]
# thr3 - lower threshold for coefficient of variation [if thr3==NULL then reduces till maxsize]
# thr4 - upper threshold for Variance Inflation Factor (VIF)
# maxsize - maximum number of candidate predictors [if maxsize==NULL then reduces till nsubj - 1]
# delrows - obtain complete cases by eliminating subjects
# maxs1 -  maximum number of terms in candidate models {"hard","soft1","soft2","soft3","soft4","soft5","soft6","exhaustive"}
# confs - size of the returned confidence set 
# level - level models' complexity {1: only main effects; 2: pairwise interactions}
# crit - information criteria to use {"aic","aicc","bic","qaic","qaicc"}
# thrs1 - threshold for number of candidate predictors to calculate nmods exactly
# thrs2 - threshold for nmods to perform exhaustive search
# nreps - number of replicas for genetic algorithm
# f1 - output filename
mod_select_multinom = function(dat1,grps,thr1=0.95,thr2=0.50,thr3=NULL,thr4=Inf,
  maxsize=30,delrows=FALSE,maxs1="exhaustive",confs=100,level=1,crit="aic",
  thrs1=20,thrs2=100000,nreps=2,f1){
# Note: AICC improves AIC by taking into account sample sizes (but "nnet" can only use AIC)
    library("nnet")
    library("glmulti")
    dat1 = dat1[,!apply(dat1,2,function(x) all(is.na(x)))] #remove variables with all NA
    if(delrows){
        ina = apply(dat1,1,function(x){any(is.na(x))})
        dat1 = dat1[!ina,]                  #remove all subjects with NAs
        grps = grps[!ina]                   #remove all subjects with NAs
        grps = factor(grps,levels=sort(unique(as.character(grps)))) #in case some level was lost
    }else{
        ina = apply(dat1,2,function(x){any(is.na(x))})
        dat1 = dat1[,!ina]                  #remove all variables with NAs
    }
    if(thr1 != 1.00 | thr2 != 0.00 | thr3 != 0.00 | thr4 != Inf){
        f2 = paste(f1,"_datred.log",sep="")
        dat1 = reduce_predictors_v2(dat1,grps,thr1=thr1,thr2=thr2,thr3=thr3,
          thr4=thr4,maxsize=maxsize,f1=f2)  #reduce number of candidate predictors
    }
    d1 = data.frame(grps,dat1)
    colnames(d1) = sapply(colnames(d1),function(x){gsub("[-|.]","_",x)}) #remove symbols used in "formulas"
    nvars = ncol(d1)
    vars = colnames(d1)
    hmax = 0
    tmp = 0
    while(tmp <= (nrow(d1) - 1)){       #hard max (i.e. nparam = nsubj - 1)
        hmax = hmax + 1
        tmp = (length(levels(grps)) - 1)*(hmax + 1)
    }
    hmax = hmax - 1
    if(is.numeric(maxs1)){
        maxs1 = maxs1
    }else if(maxs1 == "hard"){          #hard max nterms
        maxs1 = hmax
    }else if(maxs1 == "soft1"){         #soft max nterms
        maxs1 = nrow(d1) - 10
    }else if(maxs1 == "soft2"){         #soft max nterms
        maxs1 = nrow(d1) - 20
    }else if(maxs1 == "soft3"){         #soft max nterms
        maxs1 = nrow(d1) - 50
    }else if(maxs1 == "soft4"){         #soft max nterms
        maxs1 = nrow(d1) - 100
    }else if(maxs1 == "soft5"){         #soft max nterms
        maxs1 = nrow(d1) - 200
    }else if(maxs1 == "soft6"){         #soft max nterms
        maxs1 = nrow(d1) - 400
    }else if(maxs1 == "exhaustive"){    #max nterms so that nmods < thrs1 (i.e. force exhaustive search)
        maxs1 = 0
        tmp = 0
        while(tmp < thrs2){
            maxs1 = maxs1 + 1
            tmp = sum(sapply(1:maxs1,function(x){choose(nvars-1,x)})) + 1
        }
        maxs1 = maxs1 - 1
    }
    maxs1 = min(maxs1,hmax)             #prevent maxs1 > hard max
    maxs1 = max(maxs1,1)                #prevent maxs1 <= 0
    if(nvars - 1 <= thrs1){             #accounts for possible filters (slow)
        nmods = glmulti(y=vars[1],xr=vars[2:nvars],data=d1,level=1,method="d",
          maxsize=maxs1,report=FALSE,fitfunction="multinom",trace=FALSE)
    }else{                              #does not account for possible filters
        nmods = sum(sapply(1:maxs1,function(x){choose(nvars-1,x)})) + 1
    }
    if(nmods < thrs2){                  #exhaustive search
        res3 = glmulti(y=vars[1],xr=vars[2:nvars],data=d1,level=1,maxsize=maxs1,
          method="h",crit=crit,confsetsize=confs,plotty=FALSE,report=FALSE,
          includeobjects=TRUE,marginality=FALSE,model=TRUE,
          fitfunction="multinom",trace=FALSE)
    }else{                              #genetic algorithm
        res1 = vector("list",length=nreps)
        res2 = vector("list",length=nreps)
        for(i1 in 1:nreps){
            f2 = paste(f1,"_rep",i1,"_1",sep="")
            f3 = paste(f1,"_rep",i1,"_2",sep="")
            res1[[i1]] = glmulti(y=vars[1],xr=vars[2:nvars],data=d1,level=1,
              maxsize=maxs1,method="g",crit=crit,confsetsize=confs,plotty=FALSE,
              report=FALSE,name=f2,includeobjects=TRUE,marginality=FALSE,
              model=TRUE,fitfunction="multinom",trace=FALSE)
            popsize = 200               #Population size (default = 100)
            mutrate = 0.01              #Per locus (i.e. per term) mutation rate, [0,1], (default = 10^-3)
            sexrate = 0.2               #Sexual reproduction rate, [0,1], (default = 0.1)
            imm = 0.6                   #Immigration rate, [0,1] (default = 0.3)
            deltaM = 0.005              #Stop Rule: change in mean IC (default = 0.05)
            deltaB = 0.005              #Stop Rule: change in best IC (default = 0.05)
            conseq = 10                 #Stop Rule: times with no improvement (default = 5)
            res2[[i1]] = glmulti(y=vars[1],xr=vars[2:nvars],data=d1,level=1,
              maxsize=maxs1,method="r",crit=crit,confsetsize=confs,plotty=FALSE,
              report=FALSE,name=f3,includeobjects=TRUE,marginality=FALSE,
              model=TRUE,popsize=popsize,mutrate=mutrate,sexrate=sexrate,
              imm=imm,deltaM=deltaM,deltaB=deltaB,conseq=conseq,resumefile=f2,
              fitfunction="multinom",trace=FALSE)
        }
        if(nreps > 1){
            res3 = glmulti::consensus(xs=res2,confsetsize=confs)
        }else{
            res3 = res2
        }
    }
    t1 = weightable(res3)
    t2 = t1[t1[,crit] <= min(t1[,crit]) + 2,]
    t3 = coef(res3)
    f2 = paste(f1,".log",sep="")
    f3 = paste(f1,"_1.csv",sep="")
    f4 = paste(f1,"_2.csv",sep="")
    capture.output(print(res3),file=f2,append=FALSE,type="output")
    write.table(t2,file=f3,sep=",",row.names=TRUE,col.names=NA)
    write.table(t3,file=f4,sep=",",row.names=TRUE,col.names=NA)
    for(i1 in c("eps","tiff")){
        if(i1 == "eps"){
            postscript(file=paste(f1,".eps",sep=""),width=7.0,height=10,
              colormodel="rgb",horizontal=FALSE,onefile=FALSE,paper="special",
              family="Helvetica")
        }else{
            tiff(file=paste(f1,".tif",sep=""),units="in",width=7.0,height=10,
              res=300,compression="lzw",family="sans")
        }
        oldpar = par(mfcol=c(3,1),mar=c(2.2,2.2,1.2,1.2),mgp=c(1.0,0.2,0.0),
          oma=c(0,0,0,0),tcl=-0.2,cex=0.8,cex.lab=0.8,cex.axis=0.8,cex.main=0.9)
        par(mar=c(2.2,10.2,1.2,1.2))
        plot(res3,type="s",cex.names=0.6)
        par(mar=c(2.2,2.2,1.2,1.2))
        plot(res3,type="p")
        plot(res3,type="w")
        par(oldpar)
        invisible(dev.off())
    }
    a1 = formula(res3@objects[[1]]$terms) #formula
    a2 = res3@objects[[1]]$model        #data
    l1 = list(a1,a2)
    names(l1) = c("formula","data")
    return(l1)
}

# Fit multinomial logistic regression using ANN [account for very many levels]
# form1 - formula of the model to fit
# dat1 - data frame of variables to fit 
# col1 - colors for each group
# eqsamp - use equal-sized samples per level to avoid identification bias
# f1 - output filename
# main - title for plots
fit_multinomial_v2 = function(form1,dat1,col1=NULL,eqsamp=FALSE,f1,main){
#see https://en.wikipedia.org/wiki/Multinomial_logistic_regression
#see https://en.wikipedia.org/wiki/Artificial_neural_network
# Note: Multivariate logistic regression assumptions:
# 1. The error terms need to be independent.
# 2. The independent variables should be independent from each other (i.e. little or no multicollinearity).
# 3. Linearity between independent variables and log-odds.
# 4. Large sample sizes (i.e. 10-30 cases per parameter).
    library("nnet")
    library("gtools")
    library("car")
    if(eqsamp){
        grps = dat1[,all.vars(form1)[1]] #response variable
        levs = levels(grps)             #levels of response variable
        fmin = min(table(grps))         #minimum level frequency
        idx =  1:nrow(dat1)             #indexes of response variable
        s1 = unlist(lapply(levs,function(x){sample(idx[grps==x],fmin,
          replace=FALSE)}))
        dat1 = dat1[s1,]                #equal-sized samples per level
    }
    fit1 = multinom(form1,data=dat1,model=TRUE,trace=FALSE)
    vars = all.vars(form1)
    grps = dat1[,vars[1]]               #response variable
    levs = levels(grps)                 #levels of response variable
    nlev = length(levs)                 #number of levels of response variable
    sfit1 = summary(fit1)
    form0 = formula(paste(vars[1],"~1",sep=""))
    fit0 = multinom(form0,data=dat1,model=TRUE,trace=FALSE)
    r1 = 1 - logLik(fit1)/logLik(fit0)
    pval1 = anova(fit1,fit0)$"Pr(Chi)"[2]
    pred1 = predict(fit1,dat1)
    r2 = sum(pred1 == grps)/length(grps)
    a1 = getfit(fit1)[,1,drop=FALSE]    #estimates and 
    a2 = getfit(fit1)[,2,drop=FALSE]    #standard errors
    a3 = a1[,1]/a2[,1]                  #z-values (i.e. Wald-ratios)
    a4 = pnorm(abs(a3),0,1,low=FALSE)*2 #p-values
    t1 = cbind(estimate=a1,se=a2,zval=a3,pval=a4)
    t2 = cbind(grps,fit1$fit)
    colnames(t2) = c("Cluster",paste(rep("Prob(",nlev),gsub("$",")",levs),sep=""))
    t3 = Anova(fit1,type="II")
    f2 = paste(f1,".log",sep="")
    f3 = paste(f1,"_1.csv",sep="")
    f4 = paste(f1,"_2.csv",sep="")
    f5 = paste(f1,"_3.csv",sep="")
    capture.output(print(sfit1),file=f2,append=FALSE,type="output")
    write(paste("R McFadden = ",r1,"(",pval1,")",sep=""),file=f2,append=TRUE)
    write(paste("R count = ",r2,sep=""),file=f2,append=TRUE)
    write.table(t1,file=f3,sep=",",row.names=TRUE,col.names=NA)
    write.table(t2,file=f4,sep=",",row.names=TRUE,col.names=NA)
    write.table(t3,file=f5,sep=",",row.names=TRUE,col.names=NA)
    logitfit = gtools::logit(round(fit1$fit,12),min=-1e-12/2,max=1+1e-12/2)
    if(is.null(col1)){
        col1 = rainbow(nlev)
    }
    col2 = col1[as.numeric(grps)]       #observed group
    if(nrow(dat1) > 100){
        width1 = 7.0
        height1 = 10.0
        mfrow1 = c(2,1)
    }else{
        width1 = 7.0
        height1 = 3.5
        mfrow1 = c(1,2)
    }
    if(nlev > 15){
        ncol1 = ceiling(nlev/15)
    }else{
        ncol1 = 1
    }
    f2 = paste(f1,"_1",sep="")
    for(i1 in c("eps","tiff")){
        if(i1 == "eps"){
            postscript(file=paste(f2,".eps",sep=""),width=width1,height=height1,
              colormodel="rgb",horizontal=FALSE,onefile=FALSE,paper="special",
              family="Helvetica")
        }else{
            tiff(file=paste(f2,".tif",sep=""),units="in",width=width1,
              height=height1,res=300,compression="lzw",family="sans")
        }
        oldpar = par(mfrow=mfrow1,mar=c(2.2,2.2,1.2,1.2),mgp=c(1.0,0.2,0.0),
          tcl=-0.2,cex=0.8,cex.lab=0.8,cex.axis=0.8,cex.main=0.9)
        plot(fit1$residuals~logitfit,xlab="Predicted values",ylab="Residuals",
          main="Residuals vs Fitted",col=replicate(3,col2),cex=0.6,pch=23)
        legend("bottomleft",inset=c(0.02,0.02),fill=col1,legend=levs,cex=0.5,
          ncol=ncol1)
        abline(h=0,col=8,lty=3)
        barplot(t(fit1$fit),col=col1,las=2,ylab="Probability",main=main,
          cex.names=0.6)
        legend("topright",inset=c(0.01,0.01),fill=col1,legend=levs,bg="white",
          cex=0.5,ncol=ncol1)
        par(oldpar)
        invisible(dev.off())
    }
    vars = fit1$coefnames               #predictors (including "intercept")
    nvar = length(vars)                 #number of predictors (including "intercept")
    y1 = apply(fit1$fit,1,max)          #maximum probability
    y2 = fit1$residuals                 #residuals
    col3 = col1[apply(fit1$fit,1,which.max)]  #predicted group
    col3[col2 != col3] = 1              #uncolor when predicted!=observed
    nplot = 3                           #Three 1-by-2 plots
    for(i1 in c("eps","tiff")){
        f2 = paste(f1,"_2",sep="")
        if(i1 == "eps"){
            postscript(file=paste(f2,".eps",sep=""),width=7.0,height=10,
              colormodel="rgb",horizontal=FALSE,onefile=FALSE,paper="special",
              family="Helvetica")
        }else{
            tiff(file=paste(f2,".tif",sep=""),units="in",width=7.0,height=10,
              res=300,compression="lzw",family="sans")
        }
        oldpar = par(mfrow=c(3,2),mar=c(2.2,2.2,1.2,1.2),mgp=c(1.0,0.2,0.0),
          tcl=-0.2,cex=0.8,cex.lab=0.8,cex.axis=0.8,cex.main=0.9)
        iplot = 0
        for(i2 in 2:nvar){
            iplot = iplot + 1
            xlab1 = vars[i2]
            x = fit1$model[,xlab1]
            plot(y1~x,col=col3,ylab="Fitted value",ylim=c(0,1),xlab=xlab1)
            text(y1~x,labels=rownames(fit1$model),cex=0.6,pos=1,col=col3)
            abline(h=1/nlev,lty=2)
            legend("bottomright",inset=c(0.02,0.02),fill=col1,legend=levs,
              cex=0.5,ncol=ncol1)
            x = replicate(nlev,fit1$model[,xlab1])
            plot(y2~x,col=col3,ylab="Residuals",xlab=xlab1)
            text(y2~x,labels=rownames(fit1$model),cex=0.6,pos=1,col=col3)
            abline(h=0,lty=2)
            legend("bottomright",inset=c(0.02,0.02),fill=col1,legend=levs,
              cex=0.5,ncol=ncol1)
            if(iplot%%nplot==0 & iplot < nvar - 1){
                par(oldpar)
                invisible(dev.off())
                f2 = paste(f1,2+iplot/nplot,sep="_")
                if(i1 == "eps"){
                    postscript(file=paste(f2,".eps",sep=""),width=7.0,
                      height=10.0,colormodel="rgb",horizontal=FALSE,
                      onefile=FALSE,paper="special",family="Helvetica")
                }else{
                    tiff(file=paste(f2,".tif",sep=""),units="in",width=7.0,
                      height=10.0,res=300,compression = "lzw",family="sans")
                }
                oldpar = par(mfrow=c(3,2),mar=c(2.2,2.2,1.2,1.2),mgp=c(1.0,0.2,
                  0.0),tcl=-0.2,cex=0.8,cex.lab=0.8,cex.axis=0.8,cex.main=0.9)
            }
        }
        par(oldpar)
        invisible(dev.off())
    }
}

# Fits multiple regression via neural networks
# [adapted from nnet::nnet.formula()]
# formula - formula of the regression
# data - data frame of variables to fit
# weights - weights for each subject
# ... - arguments for "nnet" {e.g. size, decay, maxit, abstol, reltol, Hess, trace}
# na.action - usual model parameter
# contrast - usual model parameter
multivar = function(formula,data,weights,...,subset,na.action,contrasts=NULL){
    m = match.call(expand.dots=FALSE)
    if(is.matrix(eval.parent(m$data))){
        m$data <- as.data.frame(data)
    }
    m$... = m$contrasts = NULL
    m[[1L]] = quote(stats::model.frame)
    m = eval.parent(m)
    Terms = attr(m, "terms")
    x = model.matrix(Terms,m,contrasts)
    cons = attr(x,"contrast")
    xint = match("(Intercept)",colnames(x),nomatch = 0L)
    if(xint > 0L){
        x = x[,-xint,drop=FALSE]
    }
    if(attr(Terms,"intercept") == 0){
        Qr = qr(x)
    }else{
        Qr = qr(cbind(1,x))
        colnames(Qr$qr) = c("(Intercept)",colnames(Qr$qr)[-1])
    }
    w = model.weights(m)
    if(length(w) == 0L){
        w = rep(1,nrow(x))
    }
    y = model.response(m)
    miny = min(y)
    maxy = max(y)
    y_tr = minmax(y,miny,maxy)
    res = nnet.default(x,y_tr,w,linout=TRUE,...)
    res$terms = Terms
    res$coefnames = colnames(x)
    res$call = match.call()
    res$na.action = attr(m,"na.action")
    res$contrasts = cons
    res$call$formula = formula
    res$xlevels = .getXlevels(Terms, m)
    if(0){
        res$tr.fitted.values = matrix(minmax_inv(res$fitted.values[,1],
          miny,maxy))
        res$tr.residuals = matrix(y - res$tr.fitted.values)
    }else{
        res$tr.fitted.values = res$fitted.values
        res$tr.residuals = res$residuals
        res$fitted.values = matrix(minmax_inv(res$fitted.values[,1],miny,maxy))
        res$residuals = matrix(y - res$fitted.values)
    }
    res$deviance = sum(res$residuals^2)
    res$qr = Qr
    res$rank = res$qr$rank
    res$df.residual = length(y) - res$rank 
    res$model = m
    class(res) = c("multivar","nnet.formula","nnet")
    res
}

# Create logLik() function for "multivar" objects
# [adapted from stats:::logLik.lm()]
# object - "multivar" object
logLik.multivar = function(object,...){
    res = object$residuals
    p = object$rank
    N = length(res)
    if(is.null(w <- object$weights)){
        w = rep.int(1,N)
    }else{
        excl = w == 0
        if(any(excl)){
            res = res[!excl]
            N = length(res)
            w = w[!excl]
        }
    }
    val = 0.5*(sum(log(w)) - N*(log(2*pi) + 1 - log(N) + log(sum(w*res^2))))
    attr(val, "nobs") = N
    attr(val, "df") = p + 1
    class(val) = "logLik"
    val
}

# Create nobs() function for "multivar" objects
# [adapted from stats:::nobs.lm()]
# object - "multivar" object
nobs.multivar = function(object){
    if(!is.null(w <- object$weights)){
        sum(w != 0)
    }else{
        nrow(object$residuals)
    }
}

# Automatic model selection for multiple regression using ANN [add finer tune for ANN]
# dat1 - data frame of independet variables (predictors)
# resp - dependent variable (response variable)
# thr1 - upper threshold for correlation between predictors 
# thr2 - lower threshold for correlation between predictors and response variables [if thr2==NULL then reduces till maxsize]
# thr3 - lower threshold for coefficient of variation [if thr3==NULL then reduces till maxsize]
# thr4 - upper threshold for Variance Inflation Factor (VIF)
# maxsize - maximum number of candidate predictors [if maxsize==NULL then reduces till nsubj - 1]
# delrows - obtain complete cases by eliminating subjects
# maxs1 -  maximum number of terms in candidate models {"hard","soft1","soft2","soft3","soft4","soft5","soft6","exhaustive"}
# nn_size - number of units in the hidden layer [ANN]
# nn_decay - parameter for weight decay [ANN]
# nn_maxit - maximum number of iterations [ANN]
# nn_abstol - absolute fit criterion [ANN]
# nn_reltol - relative fit criterion [ANN]
# confs - size of the returned confidence set 
# level - level models' complexity {1: only main effects; 2: pairwise interactions}
# crit - information criteria to use {"aic","aicc","bic","qaic","qaicc"}
# thrs1 - threshold for number of candidate predictors to calculate nmods exactly
# thrs2 - threshold for nmods to perform exhaustive search
# nreps - number of replicas for genetic algorithm
# f1 - output filename
mod_select_nnet_v2 = function(dat1,resp,thr1=0.95,thr2=0.50,thr3=NULL,thr4=Inf,
  maxsize=30,delrows=FALSE,maxs1="exhaustive",nn_size=10,nn_decay=0.01,
  nn_maxit=1000,nn_abstol=1.0e-4,nn_reltol=1.0e-8,confs=100,level=1,crit="aicc",
  thrs1=20,thrs2=100000,nreps=2,f1){
    library("nnet")
    library("glmulti")
    dat1 = dat1[,!apply(dat1,2,function(x) all(is.na(x)))] #remove variables with all NA
    if(delrows){
        ina = apply(dat1,1,function(x){any(is.na(x))})
        dat1 = dat1[!ina,]                  #remove all subjects with NAs
        resp = resp[!ina]                   #remove all subjects with NAs
    }else{
        ina = apply(dat1,2,function(x){any(is.na(x))})
        dat1 = dat1[,!ina]                  #remove all variables with NAs
    }
    if(thr1 != 1.00 | thr2 != 0.00 | thr3 != 0.00 | thr4 != Inf){
        f2 = paste(f1,"_datred.log",sep="")
        dat1 = reduce_predictors_v2(dat1,grps,thr1=thr1,thr2=thr2,thr3=thr3,
          thr4=thr4,maxsize=maxsize,f1=f2)  #reduce number of candidate predictors
    }
    d1 = data.frame(resp,dat1)
    colnames(d1) = sapply(colnames(d1),function(x){gsub("[-|.]","_",x)}) #remove symbols used in "formulas"
    nvars = ncol(d1)
    vars = colnames(d1)
    a1 = ifelse(crit=="aicc",3,2)
    hmax = nrow(d1) - 1 - a1            #hard max (i.e. nparam = nsubj - [2 or 3])
    if(is.numeric(maxs1)){
        maxs1 = maxs1
    }else if(maxs1 == "hard"){          #hard max nterms
        maxs1 = hmax
    }else if(maxs1 == "soft1"){         #soft max nterms
        maxs1 = nrow(d1) - 10
    }else if(maxs1 == "soft2"){         #soft max nterms
        maxs1 = nrow(d1) - 20
    }else if(maxs1 == "soft3"){         #soft max nterms
        maxs1 = nrow(d1) - 50
    }else if(maxs1 == "soft4"){         #soft max nterms
        maxs1 = nrow(d1) - 100
    }else if(maxs1 == "soft5"){         #soft max nterms
        maxs1 = nrow(d1) - 200
    }else if(maxs1 == "soft6"){         #soft max nterms
        maxs1 = nrow(d1) - 400
    }else if(maxs1 == "exhaustive"){    #max nterms so that nmods < thrs1 (i.e. force exhaustive search)
        maxs1 = 0
        tmp = 0
        while(tmp < thrs2){
            maxs1 = maxs1 + 1
            tmp = sum(sapply(1:maxs1,function(x){choose(nvars-1,x)})) + 1
        }
        maxs1 = maxs1 - 1
    }
    maxs1 = min(maxs1,hmax)             #prevent maxs1 > hard max
    maxs1 = max(maxs1,1)                #prevent maxs1 <= 0
    if(nvars - 1 <= thrs1){             #accounts for possible filters (slow)
        nmods = glmulti(y=vars[1],xr=vars[2:nvars],data=d1,level=1,method="d",
          maxsize=maxs1,report=FALSE,fitfunction="multivar",size=nn_size,
          decay=nn_decay,maxit=nn_maxit,abstol=nn_abstol,reltol=nn_reltol,
          trace=FALSE)
    }else{                              #does not account for possible filters
        nmods = sum(sapply(1:maxs1,function(x){choose(nvars-1,x)})) + 1
    }
    if(nmods < thrs2){                  #exhaustive search
        res3 = glmulti(y=vars[1],xr=vars[2:nvars],data=d1,level=1,maxsize=maxs1,
          method="h",crit=crit,confsetsize=confs,plotty=FALSE,report=FALSE,
          includeobjects=TRUE,marginality=FALSE,model=TRUE,
          fitfunction="multivar",size=nn_size,decay=nn_decay,maxit=nn_maxit,
          abstol=nn_abstol,reltol=nn_reltol,trace=FALSE)
    }else{                              #genetic algorithm
        res1 = vector("list",length=nreps)
        res2 = vector("list",length=nreps)
        for(i1 in 1:nreps){
            f2 = paste(f1,"_rep",i1,"_1",sep="")
            f3 = paste(f1,"_rep",i1,"_2",sep="")
            res1[[i1]] = glmulti(y=vars[1],xr=vars[2:nvars],data=d1,level=1,
              maxsize=maxs1,method="g",crit=crit,confsetsize=confs,plotty=FALSE,
              report=FALSE,name=f2,includeobjects=TRUE,marginality=FALSE,
              model=TRUE,fitfunction="multivar",size=nn_size,decay=nn_decay,
              maxit=nn_maxit,abstol=nn_abstol,reltol=nn_reltol,trace=FALSE)
            popsize = 200               #Population size (default = 100)
            mutrate = 0.01              #Per locus (i.e. per term) mutation rate, [0,1], (default = 10^-3)
            sexrate = 0.2               #Sexual reproduction rate, [0,1], (default = 0.1)
            imm = 0.6                   #Immigration rate, [0,1] (default = 0.3)
            deltaM = 0.005              #Stop Rule: change in mean IC (default = 0.05)
            deltaB = 0.005              #Stop Rule: change in best IC (default = 0.05)
            conseq = 10                 #Stop Rule: times with no improvement (default = 5)
            res2[[i1]] = glmulti(y=vars[1],xr=vars[2:nvars],data=d1,level=1,
              maxsize=maxs1,method="r",crit=crit,confsetsize=confs,plotty=FALSE,
              report=FALSE,name=f3,includeobjects=TRUE,marginality=FALSE,
              model=TRUE,popsize=popsize,mutrate=mutrate,sexrate=sexrate,
              imm=imm,deltaM=deltaM,deltaB=deltaB,conseq=conseq,resumefile=f2,
              fitfunction="multivar",size=nn_size,decay=nn_decay,
              maxit=nn_maxit,abstol=nn_abstol,reltol=nn_reltol,trace=FALSE)
        }
        if(nreps > 1){
            res3 = glmulti::consensus(xs=res2,confsetsize=confs)
        }else{
            res3 = res2
        }
    }
    t1 = weightable(res3)
    t2 = t1[t1[,crit] <= min(t1[,crit]) + 2,]
    f2 = paste(f1,".log",sep="")
    f3 = paste(f1,".csv",sep="")
    capture.output(print(res3),file=f2,append=FALSE,type="output")
    write.table(t2,file=f3,sep=",",row.names=TRUE,col.names=NA)
    for(i1 in c("eps","tiff")){
        if(i1 == "eps"){
            postscript(file=paste(f1,".eps",sep=""),width=7.0,height=10,
              colormodel="rgb",horizontal=FALSE,onefile=FALSE,paper="special",
              family="Helvetica")
        }else{
            tiff(file=paste(f1,".tif",sep=""),units="in",width=7.0,height=10,
              res=300,compression="lzw",family="sans")
        }
        oldpar = par(mfcol=c(3,1),mar=c(2.2,2.2,1.2,1.2),mgp=c(1.0,0.2,0.0),
          oma=c(0,0,0,0),tcl=-0.2,cex=0.8,cex.lab=0.8,cex.axis=0.8,cex.main=0.9)
        par(mar=c(2.2,10.2,1.2,1.2))
        plot(res3,type="s",cex.names=0.6)
        par(mar=c(2.2,2.2,1.2,1.2))
        plot(res3,type="p")
        plot(res3,type="w")
        par(oldpar)
        invisible(dev.off())
    }
    a1 = formula(res3@objects[[1]]$terms) #formula
    a2 = res3@objects[[1]]$model           #data
    l1 = list(a1,a2)
    names(l1) = c("formula","data")
    return(l1)
}

# Create summary() function for "multivar" objects
# [adapted from stats:::summary.lm()]
# object - "multivar" object
summary_multivar = function(object){
    z = object
    p = z$rank
    rdf = z$df.residual
    if(p == 0){
        r = z$residuals
        n = length(r)
        w = z$weights
        if(is.null(w)){
            rss = sum(r^2)
        }else{
            rss = sum(w*r^2)
            r = sqrt(w)*r
        }
        resvar = rss/rdf
        ans = z[c("call", "terms", if (!is.null(z$weights)) "weights")]
        ans$residuals = r
        ans$df = n
        ans$sigma = sqrt(resvar)
        ans$r.squared = ans$adj.r.squared <- 0
        return(ans)
    }
    Qr = object$qr
    n = nrow(Qr$qr)
    r = z$residuals
    f = z$fitted.values
    w = z$weights
    if(is.null(w)){
        mss = if(attr(z$terms, "intercept")){sum((f - mean(f))^2)}else{sum(f^2)}
        rss = sum(r^2)
    }else{
        mss = if(attr(z$terms, "intercept")){m = sum(w * f/sum(w));
          sum(w*(f - m)^2)}else{sum(w*f^2)}
        rss = sum(w*r^2)
        r = sqrt(w)*r
    }
    resvar = rss/rdf
    if(is.finite(resvar) && resvar < (mean(f)^2 + var(f))*1e-30) 
        warning("essentially perfect fit: summary may be unreliable")
    p1 = 1L:p
    ans = z[c("call", "terms",if(!is.null(z$weights)){"weights"})]
    ans$residuals = r
    ans$sigma = sqrt(resvar)
    ans$df = c(p,rdf,ncol(Qr$qr))
    if(p != attr(z$terms,"intercept")){
        df.int = if(attr(z$terms, "intercept")){1L}else{0L}
        ans$r.squared = mss/(mss + rss)
        ans$adj.r.squared = 1 - (1 - ans$r.squared)*((n - df.int)/rdf)
        ans$fstatistic = c(value=(mss/(p - df.int))/resvar,numdf=p - df.int,
          dendf=rdf)
    }else{
        ans$r.squared = ans$adj.r.squared = 0
    }
    if(!is.null(z$na.action)){
        ans$na.action = z$na.action
    }
    ans
}

# Fit multiple regression using ANN [add finer tune for ANN]
# form1 - formula of the model to fit
# dat1 - data frame of variables to fit 
# nn_size - number of units in the hidden layer [ANN]
# nn_decay - parameter for weight decay [ANN]
# nn_maxit - maximum number of iterations [ANN]
# nn_abstol - absolute fit criterion [ANN]
# nn_reltol - relative fit criterion [ANN]
# f1 - output filename
# main - title for plots
fit_nnet_v2 = function(form1,dat1,nn_size=10,nn_decay=0.01,nn_maxit=1000,
  nn_abstol=1e-4,nn_reltol=1e-8,f1,main){
#see https://en.wikipedia.org/wiki/Artificial_neural_network
# Note: Multiple regression assumptions:
# 1. The error terms need to be independent.
# 2. The independent variables should be independent from each other (i.e. little or no multicollinearity).
# 3. Large sample sizes (i.e. 10-30 cases per parameter).
# Note: Neuronal network caveats
# 1. optimal number of hidden units is problem specific (risk of overfitting)
    library("nnet")
    library("caret")
    library("NeuralNetTools")
    vars = all.vars(form1)
    if(is.null(nn_size) || is.null(nn_decay)){
        if(!is.null(nn_size)){
            gr_size = nn_size
        }else{
            gr_size = c(1,3*(1:10))
        }
        if(!is.null(nn_decay)){
            gr_decay = nn_decay
        }else{
            gr_decay = 1*10^(-(1:9))
        }
        y = as.matrix(dat1[,vars[1]])
        miny = min(y)
        maxy = max(y)
        dat2 = dat1
        dat2[,vars[1]] = minmax(y,miny,maxy)
        grid1 = expand.grid(size=gr_size,decay=gr_decay)
        train1 = train(form1,data=dat2,method="nnet",maxit=nn_maxit,
          abstol=nn_abstol,reltol=nn_reltol,linout=TRUE,trace=FALSE,
          trControl=trainControl(method="boot",number=100),tuneGrid=grid1,
          metric="RMSE")
        nn_size = train1$bestTune[1,1]
        nn_decay = train1$bestTune[1,2]
    }
    fit1 = multivar(form1,dat1,size=nn_size,decay=nn_decay,maxit=nn_maxit,
      abstol=nn_abstol,reltol=nn_reltol,trace=FALSE)
    sfit1 = summary(fit1)
    sfit2 = summary_multivar(fit1)
    if(length(vars) > 1){
        t1 = olden(fit1,bar_plot=FALSE)
    }
    f2 = paste(f1,"_1.log",sep="")
    f3 = paste(f1,"_2.log",sep="")
    f4 = paste(f1,"_3.log",sep="")
    f5 = paste(f1,".csv",sep="")
    capture.output(print(sfit1),file=f2,append=FALSE,type="output")
    capture.output(stats:::print.summary.lm(sfit2),file=f3,append=FALSE,
      type="output")
    if(exists("train1")){
        capture.output(train1,file=f4,append=FALSE,type="output")
    }
    if(length(vars) > 1){
        write.table(t1,file=f5,sep=",",row.names=TRUE,col.names=NA)
    }
    f2 = paste(f1,"_1",sep="")
    y = dat1[,vars[1]]
    f = fit1$fitted.values
    r = fit1$residuals
    for(i1 in c("eps","tiff")){
        if(i1 == "eps"){
            postscript(file=paste(f2,".eps",sep=""),width=7.0,height=7.0,
              colormodel="rgb",horizontal=FALSE,onefile=FALSE,paper="special",
              family="Helvetica")
        }else{
            tiff(file=paste(f2,".tif",sep=""),units="in",width=7.0,height=7.0,
              res=300,compression="lzw",family="sans")
        }
        oldpar = par(mfcol=c(2,2),mar=c(2.2,2.2,1.2,1.2),mgp=c(1.0,0.2,0.0),
          tcl=-0.2,cex=0.8,cex.lab=0.8,cex.axis=0.8,cex.main=0.9)

        qqnorm(r,col=4); qqline(r)                        #wg-err_norm: qqplot 
        plot(f,sqrt(abs(r)),xlab="Fitted",ylab="sqrt(abs(Residuals))",col=4)
        abline(h=0,lty=2)                                 #wg-err_var: scatterplot
        hist(r,freq=FALSE,main="",col=4,xlab="Residuals") #wg-err_norm: histogram
        plot(f,y,xlab="Fitted values",ylab="Values",main=main,col=4)
        abline(c(0,1),lty=2)                              #goodness-of-fit: scatterplot
        par(oldpar)
        invisible(dev.off())
    }
    if(length(vars) > 1){
        f2 = paste(f1,"_2",sep="")
        for(i1 in c("eps","tiff")){
            if(i1 == "eps"){
                postscript(file=paste(f2,".eps",sep=""),width=7.0,height=7.0,
                  colormodel="rgb",horizontal=FALSE,onefile=FALSE,
                  paper="special",family="Helvetica")
            }else{
                tiff(file=paste(f2,".tif",sep=""),units="in",width=7.0,
                  height=7.0,res=300,compression="lzw",family="sans")
            }
            oldpar = par(mar=c(1.1,1.1,2.1,1.1),mgp=c(0,0,0),xpd=TRUE,
              cex.main=0.9)
            plotnet(fit1,cex_val=0.7,circle_cex=3)
            par(oldpar)
            invisible(dev.off())
        }
    }
}

# Automatic model selection for multiple linear regression using General LM
# dat1 - data frame of independet variables (predictors)
# resp - dependent variable (response variable)
# thr1 - upper threshold for correlation between predictors 
# thr2 - lower threshold for correlation between predictors and response variables [if thr2==NULL then reduces till maxsize]
# thr3 - lower threshold for coefficient of variation [if thr3==NULL then reduces till maxsize]
# thr4 - upper threshold for Variance Inflation Factor (VIF)
# maxsize - maximum number of candidate predictors [if maxsize==NULL then reduces till nsubj - 1]
# delrows - obtain complete cases by eliminating subjects
# maxs1 -  maximum number of terms in candidate models {"hard","soft1","soft2","soft3","soft4","soft5","soft6","exhaustive"}
# confs - size of the returned confidence set 
# level - level models' complexity {1: only main effects; 2: pairwise interactions}
# crit - information criteria to use {"aic","aicc","bic","qaic","qaicc"}
# thrs1 - threshold for number of candidate predictors to calculate nmods exactly
# thrs2 - threshold for nmods to perform exhaustive search
# nreps - number of replicas for genetic algorithm
# f1 - output filename
mod_select_lm = function(dat1,resp,thr1=0.95,thr2=0.50,thr3=NULL,thr4=Inf,
  maxsize=30,delrows=FALSE,maxs1="exhaustive",confs=100,level=1,crit="aicc",
  thrs1=20,thrs2=100000,nreps=2,f1){
    library("glmulti")
    dat1 = dat1[,!apply(dat1,2,function(x) all(is.na(x)))] #remove variables with all NA
    if(delrows){
        ina = apply(dat1,1,function(x){any(is.na(x))})
        dat1 = dat1[!ina,]                  #remove all subjects with NAs
        resp = resp[!ina]                   #remove all subjects with NAs
    }else{
        ina = apply(dat1,2,function(x){any(is.na(x))})
        dat1 = dat1[,!ina]                  #remove all variables with NAs
    }
    if(thr1 != 1.00 | thr2 != 0.00 | thr3 != 0.00 | thr4 != Inf){
        f2 = paste(f1,"_datred.log",sep="")
        dat1 = reduce_predictors_v2(dat1,grps,thr1=thr1,thr2=thr2,thr3=thr3,
          thr4=thr4,maxsize=maxsize,f1=f2)  #reduce number of candidate predictors
    }
    d1 = data.frame(resp,dat1)
    colnames(d1) = sapply(colnames(d1),function(x){gsub("[-|.]","_",x)}) #remove symbols used in "formulas"
    nvars = ncol(d1)
    vars = colnames(d1)
    a1 = ifelse(crit=="aicc",3,2)
    hmax = nrow(d1) - 1 - a1            #hard max (i.e. nparam = nsubj - [2 or 3])
    if(is.numeric(maxs1)){
        maxs1 = maxs1
    }else if(maxs1 == "hard"){          #hard max nterms
        maxs1 = hmax
    }else if(maxs1 == "soft1"){         #soft max nterms
        maxs1 = nrow(d1) - 10
    }else if(maxs1 == "soft2"){         #soft max nterms
        maxs1 = nrow(d1) - 20
    }else if(maxs1 == "soft3"){         #soft max nterms
        maxs1 = nrow(d1) - 50
    }else if(maxs1 == "soft4"){         #soft max nterms
        maxs1 = nrow(d1) - 100
    }else if(maxs1 == "soft5"){         #soft max nterms
        maxs1 = nrow(d1) - 200
    }else if(maxs1 == "soft6"){         #soft max nterms
        maxs1 = nrow(d1) - 400
    }else if(maxs1 == "exhaustive"){    #max nterms so that nmods < thrs1 (i.e. force exhaustive search)
        maxs1 = 0
        tmp = 0
        while(tmp < thrs2){
            maxs1 = maxs1 + 1
            tmp = sum(sapply(1:maxs1,function(x){choose(nvars-1,x)})) + 1
        }
        maxs1 = maxs1 - 1
    }
    maxs1 = min(maxs1,hmax)             #prevent maxs1 > hard max
    maxs1 = max(maxs1,1)                #prevent maxs1 <= 0
    if(nvars - 1 <= thrs1){             #accounts for possible filters (slow)
        nmods = glmulti(y=vars[1],xr=vars[2:nvars],data=d1,level=1,method="d",
          maxsize=maxs1,report=FALSE,fitfunction="lm")
    }else{                              #does not account for possible filters
        nmods = sum(sapply(1:maxs1,function(x){choose(nvars-1,x)})) + 1
    }
    if(nmods < thrs2){                  #exhaustive search
        res3 = glmulti(y=vars[1],xr=vars[2:nvars],data=d1,level=1,maxsize=maxs1,
          method="h",crit=crit,confsetsize=confs,plotty=FALSE,report=FALSE,
          includeobjects=TRUE,marginality=FALSE,model=TRUE,fitfunction="lm")
    }else{                              #genetic algorithm
        res1 = vector("list",length=nreps)
        res2 = vector("list",length=nreps)
        for(i1 in 1:nreps){
            f2 = paste(f1,"_rep",i1,"_1",sep="")
            f3 = paste(f1,"_rep",i1,"_2",sep="")
            res1[[i1]] = glmulti(y=vars[1],xr=vars[2:nvars],data=d1,level=1,
              maxsize=maxs1,method="g",crit=crit,confsetsize=confs,plotty=FALSE,
              report=FALSE,name=f2,includeobjects=TRUE,marginality=FALSE,
              model=TRUE,fitfunction="lm")
            popsize = 200               #Population size (default = 100)
            mutrate = 0.01              #Per locus (i.e. per term) mutation rate, [0,1], (default = 10^-3)
            sexrate = 0.2               #Sexual reproduction rate, [0,1], (default = 0.1)
            imm = 0.6                   #Immigration rate, [0,1] (default = 0.3)
            deltaM = 0.005              #Stop Rule: change in mean IC (default = 0.05)
            deltaB = 0.005              #Stop Rule: change in best IC (default = 0.05)
            conseq = 10                 #Stop Rule: times with no improvement (default = 5)
            res2[[i1]] = glmulti(y=vars[1],xr=vars[2:nvars],data=d1,level=1,
              maxsize=maxs1,method="r",crit=crit,confsetsize=confs,plotty=FALSE,
              report=FALSE,name=f3,includeobjects=TRUE,marginality=FALSE,
              model=TRUE,popsize=popsize,mutrate=mutrate,sexrate=sexrate,
              imm=imm,deltaM=deltaM,deltaB=deltaB,conseq=conseq,resumefile=f2,
              fitfunction="lm")
        }
        if(nreps > 1){
            res3 = glmulti::consensus(xs=res2,confsetsize=confs)
        }else{
            res3 = res2
        }
    }
    t1 = weightable(res3)
    t2 = t1[t1[,crit] <= min(t1[,crit]) + 2,]
    t3 = coef(res3)
    f2 = paste(f1,".log",sep="")
    f3 = paste(f1,"_1.csv",sep="")
    f4 = paste(f1,"_2.csv",sep="")
    capture.output(print(res3),file=f2,append=FALSE,type="output")
    write.table(t2,file=f3,sep=",",row.names=TRUE,col.names=NA)
    write.table(t3,file=f4,sep=",",row.names=TRUE,col.names=NA)
    for(i1 in c("eps","tiff")){
        if(i1 == "eps"){
            postscript(file=paste(f1,".eps",sep=""),width=7.0,height=10,
              colormodel="rgb",horizontal=FALSE,onefile=FALSE,paper="special",
              family="Helvetica")
        }else{
            tiff(file=paste(f1,".tif",sep=""),units="in",width=7.0,height=10,
              res=300,compression="lzw",family="sans")
        }
        oldpar = par(mfcol=c(3,1),mar=c(2.2,2.2,1.2,1.2),mgp=c(1.0,0.2,0.0),
          oma=c(0,0,0,0),tcl=-0.2,cex=0.8,cex.lab=0.8,cex.axis=0.8,cex.main=0.9)
        par(mar=c(2.2,10.2,1.2,1.2))
        plot(res3,type="s",cex.names=0.6)
        par(mar=c(2.2,2.2,1.2,1.2))
        plot(res3,type="p")
        plot(res3,type="w")
        par(oldpar)
        invisible(dev.off())
    }
    a1 = formula(res3@objects[[1]]$terms) #formula
    a2 = res3@objects[[1]]$model        #data
    l1 = list(a1,a2)
    names(l1) = c("formula","data")
    return(l1)
}

# Fit multiple linear regression using General LS
# form1 - formula of the model to fit
# dat1 - data frame of variables to fit 
# f1 - output filename
# main - title for plots
fit_lm = function(form1,dat1,f1,main){
#see https://en.wikipedia.org/wiki/General_linear_model
# Note: Multiple linear regression assumptions:
# 1. The error terms need to be independent.
# 2. The error variance is homoskedastic.
# 3. The errors are normally distributed with mean = 0
# 4. The independent variables should be independent from each other (i.e. little or no multicollinearity).
# 5. Linearity between independent variables and dependent variable.
# 6. Large sample sizes (i.e. 10-30 cases per parameter).
    library("car")
    fit1 = lm(form1,data=dat1)
    sfit1 = summary(fit1)
    t1 = Anova(fit1,type="II")
    f2 = paste(f1,".log",sep="")
    f3 = paste(f1,".csv",sep="")
    capture.output(print(sfit1),file=f2,append=FALSE,type="output")
    write.table(t1,file=f3,sep=",",row.names=TRUE,col.names=NA)
    vars = all.vars(form1)
    y = dat1[,vars[1]]
    f = fit1$fitted.values
    r = fit1$residuals
    for(i1 in c("eps","tiff")){
        if(i1 == "eps"){
            postscript(file=paste(f1,".eps",sep=""),width=7.0,height=7.0,
              colormodel="rgb",horizontal=FALSE,onefile=FALSE,paper="special",
              family="Helvetica")
        }else{
            tiff(file=paste(f1,".tif",sep=""),units="in",width=7.0,height=7.0,
              res=300,compression="lzw",family="sans")
        }
        oldpar = par(mfcol=c(2,2),mar=c(2.2,2.2,1.2,1.2),mgp=c(1.0,0.2,0.0),
          tcl=-0.2,cex=0.8,cex.lab=0.8,cex.axis=0.8,cex.main=0.9)
        qqnorm(r,col=4); qqline(r)                        #wg-err_norm: qqplot 
        plot(f,sqrt(abs(r)),xlab="Fitted",ylab="sqrt(abs(Residuals))",col=4)
        abline(h=0,lty=2)                                 #wg-err_var: scatterplot
        hist(r,freq=FALSE,main="",col=4,xlab="Residuals") #wg-err_norm: histogram
        plot(f,y,xlab="Fitted values",ylab="Values",main=main,col=4)
        abline(c(0,1),lty=2)                              #goodness-of-fit: scatterplot
        par(oldpar)
        invisible(dev.off())
    }
}

# Perform Over-representation (a.k.a. Enrichment) analysis [allow NAs]
# samp - sample
# univ - universe
# dat1 - data frame of variables to be analysed for enrichment
# prob - quantile to transform numeric to binomial variable
# dirct - direction to transform numeric to binomial variable {"ge","le"}
# stat - statistic to report {"pval","adjp"}
# cut1 - threshold for significance
# f1 - output filename
performORA_v2 = function(samp,univ,dat1,prob=0.75,dirct="ge",stat="adjp",
  cut1=0.05,f1){
    vars = colnames(dat1)
    nvars = ncol(dat1)
    subj = rownames(dat1)
    db1 = vector("list",length=nvars)
    names(db1) = vars
    nsamp = vector("numeric",length=nvars)   #sample size
    nuniv = vector("numeric",length=nvars)   #universe size
    for(i1 in 1:nvars){
        ina = is.na(dat1[,i1])
        v1 = dat1[!ina,i1]
        subj1 = subj[!ina]
        nsamp[i1] = length(samp %in% subj1)
        nuniv[i1] = length(univ %in% subj1)
        #from numeric to binary
        if(is.null(levels(v1))){
            thr1 = quantile(v1,prob)
            if(dirct == "ge"){
                db1[[i1]] = subj1[v1 >= thr1]
            }else{
                db1[[i1]] = subj1[v1 <= thr1]
            }
        #assumes base is first level of categorical variable
        }else{
            db1[[i1]] = subj1[v1 == levels(v1)[1]]
        }
    }
    insamp = vector("numeric",length=nvars)  #hit in sample per test
    inuniv = vector("numeric",length=nvars)  #hit in universe per test
    ntest = vector("numeric",length=nvars)   #number of elements per test
    for(i1 in 1:nvars){
        nam = vars[i1]
        insamp[i1] = sum(samp %in% db1[[nam]])
        inuniv[i1] = sum(univ %in% db1[[nam]])
        ntest[i1] = length(db1[[nam]])
    }
    outsamp = nsamp - insamp                 #fail in sample
    outuniv = nuniv - inuniv                 #fail in universe
    pval_up = phyper(insamp-1,nsamp,nuniv-nsamp,inuniv,lower.tail=FALSE)
    pval_dn = phyper(insamp,nsamp,nuniv-nsamp,inuniv,lower.tail=TRUE)
    adjp_up = p.adjust(pval_up,method="BH")
    adjp_dn = p.adjust(pval_dn,method="BH")
    oddsratio = (insamp/outsamp)/((inuniv - insamp)/(outuniv - outsamp))
    expected = (inuniv/nuniv)*nsamp
    df_up = data.frame(vars=vars,pval=pval_up,adjp=adjp_up,odds=oddsratio,
      exp=expected,obs=insamp,size=ntest,stringsAsFactors=FALSE)
    df_up = df_up[order(df_up[,"pval"]),]
    rownames(df_up) = 1:nrow(df_up)
    df_dn = data.frame(vars=vars,pval=pval_dn,adjp=adjp_dn,odds=oddsratio,
      exp=expected,obs=insamp,size=ntest,stringsAsFactors=FALSE)
    df_dn = df_dn[order(df_dn[,"pval"]),]
    rownames(df_dn) = 1:nrow(df_dn)
    s1 = df_up[,stat] < cut1
    s2 = df_dn[,stat] < cut1
    f2 = paste(f1,"_up.log",sep="")
    f3 = paste(f1,"_up.csv",sep="")
    f4 = paste(f1,"_dn.log",sep="")
    f5 = paste(f1,"_dn.csv",sep="")
    capture.output(
      cat("Sets to variables test for over-representation\n"),
      cat(paste(nvars,"variables tested,",sum(s1),"have p <",cut1,"\n")),
      cat(paste("Selected set size:",length(samp),"\n")),
      cat(paste("    Universe size:",length(univ),"\n")),
      file=f2,append=FALSE,type="output")
    if(sum(s1) > 0){
        capture.output(print(df_up[s1,]),file=f2,append=TRUE,type="output")
    }    
    write.table(df_up,f3,quote=TRUE,sep=",",eol="\n",na="NA",dec=".",
      row.names=FALSE,col.names=TRUE)
    capture.output(
      cat("Sets to variables test for under-representation\n"),
      cat(paste(nvars,"variables tested,",sum(s2),"have p <",cut1,"\n")),
      cat(paste("Selected set size:",length(samp),"\n")),
      cat(paste("    Universe size:",length(univ),"\n")),
      file=f4,append=FALSE,type="output")
    if(sum(s2) > 0){
        capture.output(print(df_dn[s2,]),file=f4,append=TRUE,type="output")
    }    
    write.table(df_dn,f5,quote=TRUE,sep=",",eol="\n",na="NA",dec=".",
      row.names=FALSE,col.names=TRUE)
}

# Perform Weighted Correlation Network Analysis [allow NAs]
# dat1 - data frame of variables to be analysed
# tra1 - data frame of traits (response variables) to be analysed
# mins1 - minimum module size [WCNA]
# spl1 - cluster splitting [WCNA]
# meth1 - dynamic Tree Cut {"hybrid","tree"} [WCNA]
# cut1 - merge cut height [if is.null(cut) then do not merge][WCNA]
# pwr - soft-threshold power [if is.null(pwr) then calculate best pwr] [WCNA]
# pwrrg - range for flexible selection of soft-threshold power
# thr1 - thr. for vars to consider [Standard Cor]
# ntop - no. of top vars to choose from all modules
# thr2 - thr. for choosing significant Modules
# maxs1 - max no. vars for labeling vars
# thr3 - thr. for labeling vars [Module Membership]
# thr4 - thr. for labeling vars [Trait Significance]
# thr5 - thr. (quantile) for labeling vars [Intramod Connect]
# netscreen - calculate module-based weighted vars-trait significance
# labeltop - label top vars-trait significance variables
# f1 - file name
analyse_WCNA_v2 = function(dat1,tra1,mins1=20,spl1=FALSE,meth1="hybrid",
  cut1=0.20,pwr=NULL,pwrrg=0.95,thr1=1.00,ntop=30,thr2=1.00,maxs1=50,
  thr3=0.70,thr4=0.70,thr5=0.75,netscreen=FALSE,labeltop=FALSE,f1){
#see https://en.wikipedia.org/wiki/Weighted_correlation_network_analysis
# Note: Authors advise not to use datasets with fewer than 15 samples, and recommend at least 20
# samples. Soft-threshold for scale-free topology index (samples < 20): a)hybrid_networks = 10;
# b)signed_networks = 20.
#       Grey (or 0) labels variables outside of all modules. And traits are usually continuous.
    library("WGCNA")
    library("nnet")
    allowWGCNAThreads()
    #choose power to use for adjacency calculation
    if(is.null(pwr)){
        pwrs = c(c(1:5),seq(from=6,to=10,by=2),seq(from=15,to=40,by=5))
        sfts = pickSoftThreshold(dat1,powerVector=pwrs,verbose=0)
        rsqs = -sign(sfts$fitIndices[,"slope"])*sfts$fitIndices[,"SFT.R.sq"]
        pwr = pwrs[which(rsqs > max(rsqs)*pwrrg)][1]
    }
    #network construction and module detection
    adjs = abs(WGCNA::cor(dat1,method="s",use="p"))^pwr #co-expression and adjacency
    diss = 1 - TOMsimilarity(adjs,verbose=0)    #turn adjacency into topological overlap
    cl1 = hclust(as.dist(diss),method="average")
    dynM = cutreeDynamic(dendro=cl1,minClusterSize=mins1,method=meth1,distM=diss,
      deepSplit=spl1,pamRespectsDendro=FALSE,verbose=0) #module identification using dynamic tree cut
    col1 = labels2colors(dynM)
    MEs0 = moduleEigengenes(dat1,colors=col1)$eigengenes #calculate eigengenes
    ina = apply(MEs0,1,function(x){any(is.na(x))})
    MEs0 = MEs0[!ina,]
    dat1 = dat1[!ina,]
    tra1 = tra1[!ina,]
    dist1 = 1 - cor(MEs0,method="s",use="a") #calculate dissimilarity of module eigengenes
    MEclst = hclust(as.dist(dist1),method="average") #cluster module eigengenes
    #merge modules
    if(!is.null(cut1)){
        modCol = mergeCloseModules(dat1,col1,cutHeight=cut1,verbose=0)$colors
    }else{
        modCol = col1
    }
    #calculate module eigengenes (i.e. first Principal Component)
    MEs = moduleEigengenes(dat1,modCol)$eigengenes #recalculate MEs with color labels
    MEs = orderMEs(MEs)                            #put close eigengenes close to each other
    modNam = substring(names(MEs),3)
    #calculate module membership and intramodular connectivity
    varIC = intramodularConnectivity(adjs,modCol)   #var Intramodular Connectivity
    varMM = as.data.frame(cor(dat1,MEs,method="s",use="pairwise.complete.obs")) #var Module Membership
    varMMPval = as.data.frame(corPvalueStudent(as.matrix(varMM),nrow(dat1)))
    colnames(varMM) = modNam
    colnames(varMMPval) = modNam
    varInfo = data.frame(var=colnames(dat1),module=modCol,kTotal=varIC$kTotal,
      kIn=varIC$kWithin)
    for(i1 in 1:ncol(varMM)){
      oldnam = colnames(varInfo)
      varInfo = data.frame(varInfo,varMM[,i1],varMMPval[,i1])
      colnames(varInfo) = c(oldnam,paste("MM.",modNam[i1],sep=""),paste("p.MM.",
        modNam[i1],sep=""))
    }
    modOrd <- rep(NA,length(modCol))
    for(i1 in 1:length(modNam)){
        modOrd[modCol == modNam[i1]] = i1
    }
    s1 = order(modOrd)                          #Order var by module Color
    varInfo = varInfo[s1,]
    #quantifying standard var-trait significance
    traits = colnames(tra1)
    varTS = data.frame(matrix(vector(),nrow=ncol(dat1),ncol=ncol(tra1),
      dimnames=list(colnames(dat1),traits)))
    varTSPval = data.frame(matrix(vector(),nrow=ncol(dat1),ncol=ncol(tra1),
      dimnames=list(colnames(dat1),traits)))
    for(i1 in traits){
        if(is.null(levels(tra1[,i1]))){
            varTS[,i1] = WGCNA::cor(dat1,tra1[,i1],use="p",method="s")
            varTSPval[,i1] <- corPvalueStudent(varTS[,i1],nrow(dat1))
        }else{
            for(i2 in 1:ncol(dat1)){
                ina = is.na(dat1[,i2])
                if(sum(!ina) > 2){
                    d1 = cbind(tra1[!ina,i1,drop=FALSE],dat1[!ina,i2])
                    d1[,1] = factor(d1[,1],
                      levels=sort(unique(as.numeric(as.character(d1[,1]))))) #in case some level was lost
                    fit0 = multinom(paste(i1,"~ 1"),data=d1,model=FALSE,
                      trace=FALSE)
                    fit1 = multinom(paste(i1,"~ ."),data=d1,model=FALSE,
                      trace=FALSE)
                    RMcFadden = max(0,1 - deviance(fit1)/deviance(fit0))
                    varTS[i2,i1] = sqrt(RMcFadden)
                    varTSPval[i2,i1] = anova(fit1,fit0)$"Pr(Chi)"[2]
                }
            }
        }
    }
    if(netscreen){
        #quantifying module-based weighted var-trait significance [experimental]
        varNS = vector("list",length=ncol(tra1))
        names(varNS) = traits
        for(i1 in traits){
            if(is.null(levels(tra1[,i1]))){
                varNS[[i1]] = networkScreening(y=as.integer(tra1[,i1]),
                  datME=MEs,datExpr=dat1,oddPower=3,blockSize=ncol(dat1),
                  corOptions="use='p',method='s'",minimumSampleSize=4,
                  addMEy=TRUE,removeDiag=FALSE,weightESy=0.5,getQValues=FALSE)
            }else{
                varNS[[i1]] = networkScreeningGS(datME=MEs,datExpr=dat1,
                  GS=varTS[,i1],oddPower=3,blockSize=ncol(dat1),
                  minimumSampleSize=4,addGS=FALSE)
            }
        }
    }
    #quantifying "modules' eigengenes"-based module-trait association
    modTS = data.frame(matrix(vector(),nrow=ncol(MEs),ncol=ncol(tra1),
      dimnames=list(colnames(MEs),traits)))
    modTSPval = data.frame(matrix(vector(),nrow=ncol(MEs),ncol=ncol(tra1),
      dimnames=list(colnames(MEs),traits)))
    for(i1 in traits){
        if(is.null(levels(tra1[,i1]))){
            modTS[,i1] = cor(MEs,tra1[,i1],use="p",method="spearman")
            modTSPval[,i1] = corPvalueStudent(modTS[,i1],nrow(MEs))
        }else{
            for(i2 in 1:ncol(MEs)){
                ina = is.na(MEs[,i2])
                if(sum(!ina) > 2){
                    d1 = cbind(tra1[!ina,i1,drop=FALSE],MEs[!ina,i2])
                    d1[,1] = factor(d1[,1],
                      levels=sort(unique(as.numeric(as.character(d1[,1]))))) #in case some level was lost
                    fit0 = multinom(paste(i1,"~ 1"),data=d1,model=FALSE,
                      trace=FALSE)
                    fit1 = multinom(paste(i1,"~ ."),data=d1,model=FALSE,
                      trace=FALSE)
                    RMcFadden = 1 - deviance(fit1)/deviance(fit0)
                    modTS[i2,i1] = sqrt(RMcFadden)
                    modTSPval[i2,i1] = anova(fit1,fit0)$"Pr(Chi)"[2]
                }
            }
        }
    }
    #quantifying "modules' average var-trait significance"-based module-trait association
    modTS_2 = data.frame(matrix(vector(),nrow=ncol(MEs),ncol=ncol(tra1),
      dimnames=list(colnames(MEs),traits)))
    modTSPval_2 = data.frame(matrix(vector(),nrow=ncol(MEs),ncol=ncol(tra1),
      dimnames=list(colnames(MEs),traits)))
    for(i1 in traits){
        modTS_2[,i1] = tapply(varTS[,i1],modCol,mean,na.rm=TRUE)[modNam]  #combine geneTS
        modTSPval_2[,i1] = tapply(p.adjust(varTSPval[,i1],method="BH"),modCol,
          function(x){1 - pchisq(-2*sum(log(x)),df=2*length(x))})[modNam] #combine geneTSPval
    }
    varTSInfo = vector("list",length=ncol(tra1))
    modTSInfo = vector("list",length=ncol(tra1))
    names(varTSInfo) = traits
    for(i1 in traits){
        modTSInfo[[i1]] = data.frame(module=modNam,eigenvar=modTS[,i1],
          p.eigenvar=modTSPval[,i1],meanVarTS=modTS_2[,i1],
          p.meanVarTS=modTSPval_2[,i1])
        modOrder = order(-modTS[,i1])   #Order modules by module Trait Significance
        modTSInfo[[i1]] = modTSInfo[[i1]][modOrder,]
        sel = varTSPval[,i1] <= thr1    #Standard Trait Significance p-value
        varTSInfo[[i1]] = data.frame(var=rownames(varTS[sel,]),
          module=modCol[sel],VS=varTS[sel,i1],p.VS=varTSPval[sel,i1])
        modOrd = rep(NA,length(modCol))
        for(mod in modNam){
            modOrd[modCol == mod] = modTS_2[mod,i1]
        }
        varOrd = order(-modOrd[sel],-varTSInfo[[i1]][,"VS"]) #Order var by module and var-trait significance
        varTSInfo[[i1]] = varTSInfo[[i1]][varOrd,]
        rownames(varTSInfo[[i1]]) = 1:nrow(varTSInfo[[i1]])
    }
    #save data
    t1 = table(modCol)
    t2 = cor(MEs,method="s",use="pairwise.complete.obs")
    t3 = MEs
    rownames(t3) = rownames(dat1)
    f2 = paste(f1,"_1.csv",sep="")
    f3 = paste(f1,"_2.csv",sep="")
    f4 = paste(f1,"_subjInfo.csv",sep="")
    f5 = paste(f1,"_varInfo.csv",sep="")
    f6 = paste(f1,"_modTSInfo.csv",sep="")
    f7 = paste(f1,"_varTSInfo.csv",sep="")
    f8 = paste(f1,"_modTSInfo",sep="")
    f9 = paste(f1,"_varTSInfo",sep="")
    write.table(t1,f2,quote=TRUE,sep=",",eol="\n",na="NA",dec=".",
      row.names=FALSE)
    write.table(t2,f3,quote=TRUE,sep=",",eol="\n",na="NA",dec=".",col.names=NA)
    write.table(t3,f4,quote=TRUE,sep=",",eol="\n",na="NA",dec=".",col.names=NA)
    write.table(varInfo,f5,quote=TRUE,sep=",",eol="\n",na="NA",dec=".",
      col.names=NA)
    write.table(cbind(t(modTS),t(modTSPval)),f6,quote=TRUE,sep=",",eol="\n",
      na="NA",dec=".",col.names=NA)
    write.table(cbind(t(varTS),t(varTSPval)),f7,quote=TRUE,sep=",",eol="\n",
      na="NA",dec=".",col.names=NA)
    for(i1 in traits){
        f10 = paste(f8,"_",i1,".csv",sep="")
        f11 = paste(f9,"_",i1,".csv",sep="")
        write.table(modTSInfo[[i1]],f10,quote=TRUE,sep=",",eol="\n",na="NA",
          dec=".",col.names=TRUE,row.names=FALSE)
        write.table(varTSInfo[[i1]],f11,quote=TRUE,sep=",",eol="\n",na="NA",
          dec=".",col.names=TRUE,row.names=FALSE)
    }
    #plotting
    if(exists("sfts")){
        f2 = paste(f1,"_net_pwr",sep="")
        mcon = sfts$fitIndices[,"mean.k."]
        xlab1 = "Soft Threshold (power)"
        ylab1 = "signed R^2"
        ylab2 = "Mean Connectivity"
        main1 = "Scale independence"
        main2 = "Mean connectivity"
        for(i1 in c("eps","tiff")){
            if(i1 == "eps"){
                postscript(file=paste(f2,".eps",sep=""),width=7.0,height=3.5,
                  colormodel="rgb",horizontal=FALSE,onefile=FALSE,
                  paper="special",family="Helvetica")
            }else{
                tiff(file=paste(f2,".tif",sep=""),units="in",width=7.0,
                  height=3.5,res=300,compression="lzw",family="sans")
            }
            oldpar = par(mfcol=c(1,2),mar=c(2.2,2.2,1.2,1.2),mgp=c(1.0,0.2,0.0),
              tcl=-0.2,cex=0.8,cex.lab=0.8,cex.axis=0.8,cex.main=0.9)
            plot(pwrs,rsqs,xlab=xlab1,ylab=ylab1,type="c",main=main1)
            text(pwrs,rsqs,labels=pwrs,cex=0.9,col=2)
            abline(h=c(0.5,0.7,0.8,0.9),col=2,lty=c(3,2,4,1))
            abline(h=max(rsqs)*pwrrg,col=1,lty=3)
            plot(pwrs,mcon,xlab=xlab1,ylab=ylab1,type="c",main=main2)
            text(pwrs,mcon,labels=pwrs,cex=0.9,col=2)
            par(oldpar)
            invisible(dev.off())
        }
    }
    f2 = paste(f1,"_net_mrg",sep="")
    main1 = "Clustering of module eigengenes"
    for(i1 in c("eps","tiff")){
        if(i1 == "eps"){
            postscript(file=paste(f2,".eps",sep=""),width=7.0,height=6.0,
              colormodel="rgb",horizontal=FALSE,onefile=FALSE,paper="special",
              family="Helvetica")
        }else{
            tiff(file=paste(f2,".tif",sep=""),units="in",width=7.0,height=6.0,
              res=300,compression="lzw",family="sans")
        }
        oldpar = par(mar=c(2.2,2.2,1.2,1.2),mgp=c(1.0,0.2,0.0),tcl=-0.2,cex=0.8,
          cex.lab=0.8,cex.axis=0.8,cex.main=0.9)
        plot(MEclst,main=main1,xlab="",sub="")
        if(!is.null(cut1)){
            abline(h=cut1,col="red",lty=2)
        }
        par(oldpar)
        invisible(dev.off())
    }
    f2 = paste(f1,"_net_dendr",sep="")
    if(is.null(cut1)){
        lcol = modCol
        llab = "Dynamic Cut"
    }else{
        lcol = cbind(col1,modCol)
        llab = c("Dynamic Cut",paste("Merged (",cut1,")",sep=""))
    }
    main1 = "Variables dendrogram and modules"
    for(i1 in c("eps","tiff")){
        if(i1 == "eps"){
            postscript(file=paste(f2,".eps",sep=""),width=7.0,height=5.0,
              colormodel="rgb",horizontal=FALSE,onefile=FALSE,paper="special",
              family="Helvetica")
        }else{
            tiff(file=paste(f2,".tif",sep=""),units="in",width=7.0,height=5.0,
              res=300,compression="lzw",family="sans")
        }
        oldpar = par(mar=c(2.2,2.2,1.2,1.2),mgp=c(1.0,0.2,0.0),tcl=-0.2,cex=0.8,
          cex.lab=0.8,cex.axis=0.8,cex.main=0.9)
        plotDendroAndColors(cl1,colors=lcol,groupLabels=llab,dendroLabels=FALSE,
          hang=0.03,addGuide=TRUE,guideHang=0.05,main=main1)
        par(oldpar)
        invisible(dev.off())
    }
    f2 = paste(f1,"_net_mds",sep="")
    cmd1 = cmdscale(as.dist(diss),2)
    main1 = "Multi-Dimensional Scaling plot"
    xlab1 = "Scaling 1"
    ylab2 = "Scaling 2"
    for(i1 in c("eps","tiff")){
        if(i1 == "eps"){
            postscript(file=paste(f2,".eps",sep=""),width=7.0,height=6.0,
              colormodel="rgb",horizontal=FALSE,onefile=FALSE,paper="special",
              family="Helvetica")
        }else{
            tiff(file=paste(f2,".tif",sep=""),units="in",width=7.0,height=6.0,
              res=300,compression="lzw",family="sans")
        }
        oldpar = par(mar=c(2.2,2.2,1.2,1.2),mgp=c(1.0,0.2,0.0),tcl=-0.2,cex=0.8,
          cex.lab=0.8,cex.axis=0.8,cex.main=0.9)
        plot(cmd1,col=modCol,main=main1,xlab=xlab1,ylab=ylab1)
        par(oldpar)
        invisible(dev.off())
    }
    if(ncol(MEs) > 3){
        f2 = paste(f1,"_mod_clst",sep="")
        for(i1 in c("eps","tiff")){
            if(i1 == "eps"){
                postscript(file=paste(f2,".eps",sep=""),width=5.0,height=7.5,
                  colormodel="rgb",horizontal=FALSE,onefile=FALSE,
                  paper="special",family="Helvetica")
            }else{
                tiff(file=paste(f2,".tif",sep=""),units="in",width=5.0,
                  height=7.5,res=300,compression="lzw",family="sans")
            }
            oldpar = par(mar=c(2.2,2.2,1.2,1.2),mgp=c(1.0,0.2,0.0),tcl=-0.2,
              cex=0.7,cex.lab=0.7,cex.axis=0.7,cex.main=0.8)
            plotEigengeneNetworks(MEs,setLabels="",marDendro=c(0,3.2,1,2.8),
              marHeatmap=c(2,4,1,0),xLabelsAngle=90)
            par(oldpar)
            invisible(dev.off())
        }
    }
    f2 = paste(f1,"_mod_pair",sep="")
    title1 = "Relationship between module eigengenes"
    for(i1 in c("eps","tiff")){
        if(i1 == "eps"){
            postscript(file=paste(f2,".eps",sep=""),width=7.0,height=7.0,
              colormodel="rgb",horizontal=FALSE,onefile=FALSE,paper="special",
              family="Helvetica")
        }else{
            tiff(file=paste(f2,".tif",sep=""),units="in",width=7.0,height=7.0,
              res=300,compression="lzw",family="sans")
        }
        oldpar = par(mar=c(0.5,0.5,0.5,0.5),mgp=c(3.0,0.5,0.0),tcl=-0.2,cex=0.6)
        plotMEpairs(MEs,cex.main=1.0,cex.labels=0.8,oma=c(1.2,1.2,3.1,1.2),
          main="")
        title(title1,line=-0.5)
        par(oldpar)
        invisible(dev.off())
    }
    f2 = paste(f1,"_vars_ICvsMM",sep="")
    xlab1 = "Connectivity"
    ylab1 = paste("(Module Membership)^",pwr,sep="")
    title1 = "Relationship between module eigengenes"
    for(i1 in c("eps","tiff")){
        if(i1 == "eps"){
            postscript(file=paste(f2,".eps",sep=""),width=7.0,height=10.0,
              colormodel="rgb",horizontal=FALSE,onefile=FALSE,paper="special",
              family="Helvetica")
        }else{
            tiff(file=paste(f2,".tif",sep=""),units="in",width=7.0,height=10.0,
              res=300,compression="lzw",family="sans")
        }
        oldpar = par(mfrow=c(ceiling((length(modNam)-1)/2),2),mar=c(2.0,2.3,1.0,0.2),mgp=c(1.0,0.2,0),
              tcl=-0.2,cex=0.8,cex.lab=0.8,cex.axis=0.8,cex.main=0.9)
        for(i2 in modNam){
            if(i2 != "grey"){
                s1 = modCol == i2
                plot(varIC$kWithin[s1],varMM[s1,i2]^pwr,col=i2,main=i2,
                    xlab=xlab1,ylab=ylab1,cex=1.0)
            }
        }
        par(oldpar)
        invisible(dev.off())
    }
    for(i1 in traits){
        f2 = paste(f1,"_trait_TOM_",i1,sep="")
        main1 = paste("TOM for signifcant vars (",i1,")",
          sep="")
        s1 = rank(-abs(varTS[,i1]),ties.method="first") <= ntop
        labs = rownames(varTS)[s1]; labs = gsub("[-:|]",".",labs)
        d1 = dat1[,s1]
        for(i2 in c("eps","tiff")){
            if(i2 == "eps"){
                postscript(file=paste(f2,".eps",sep=""),width=7.0,height=7.0,
                  colormodel="rgb",horizontal=FALSE,onefile=FALSE,
                  paper="special",family="Helvetica")
            }else{
                tiff(file=paste(f2,".tif",sep=""),units="in",width=7.0,
                  height=7.0,res=300,compression="lzw",family="sans")
            }
            oldpar = par(mar=c(0.5,0.5,0.5,0.5),mgp=c(3.0,0.5,0.0),tcl=-0.2,
              cex=0.6)
            plotNetworkHeatmap(d1,plotGenes=labs,networkType="unsigned",
              useTOM=TRUE,power=pwr,main=main1)
            par(oldpar)
            invisible(dev.off())
        }
    }
    if(netscreen){
        main1 = "Network-based weighted versus standard GS\n"
        main2 = "Network-based weighted vs standard sqrt(pseudo-r^2)\n"
        xlab1 = "Weighted.GS"
        xlab2 = "Weighted.sqrt(pseudo-r^2)"
        ylab1 = "Standard.GS"
        ylab2 = "Standard.sqrt(pseudo-r^2)"
        for(i1 in traits){
            f2 = paste(f1,"_trait_WSvsSS_",i1,sep="")
            for(i2 in c("eps","tiff")){
                if(i2 == "eps"){
                    postscript(file=paste(f2,".eps",sep=""),width=5,height=5,
                      colormodel="rgb",horizontal=FALSE,onefile=FALSE,
                      paper="special",family="Helvetica")
                }else{
                    tiff(file=paste(f2,".tif",sep=""),units="in",width=5,
                      height=5,res=300,compression="lzw",family="sans")
                }
                oldpar = par(mar=c(3.8,3.8,2.8,1),tck=-0.01,mgp=c(2.5,1.0,0),
                  cex=0.5,cex.main=0.5,cex.lab=0.5,cex.axis=0.5)
                if(is.null(levels(tra1[,i1]))){
                    verboseScatterplot(varNS[[i1]][,"cor.Weighted"],
                      varNS[[i1]][,"cor.Standard"],xlab=xlab1,ylab=ylab1,
                      col=modCol,main=main1)
                    abline(0,1)
                }else{
                    verboseScatterplot(varNS[[i1]][,"GS.Weighted"],varNS[[i1]][,
                      "GS"],xlab=xlab2,ylab=ylab2,col=modCol,main=main2)
                    abline(0,-1)
                    abline(0,1)
                }
                par(oldpar)
                invisible(dev.off())
            }
        }
    }
    f2 = paste(f1,"_trait_MM1",sep="")
    text1 = paste(signif(as.matrix(modTS),2),"\n(",signif(as.matrix(modTSPval),
      1),")",sep="")
    dim(text1) = dim(modTS)
    col1 = colorRampPalette(c("red","white","green"))(n=10)
    main1 = "Module-trait relationships"
    for(i1 in c("eps","tiff")){
        if(i1 == "eps"){
            postscript(file=paste(f2,".eps",sep=""),width=7.0,height=8.0,
              colormodel="rgb",horizontal=FALSE,onefile=FALSE,paper="special",
              family="Helvetica")
        }else{
            tiff(file=paste(f2,".tif",sep=""),units="in",width=7.0,height=8.0,
              res=300,compression="lzw",family="sans")
        }
        oldpar = par(mar=c(4.0,3.8,3,0.2),cex.lab=0.8,cex.main=1.0)
        labeledHeatmap(Matrix=modTS,xLabels=traits,yLabels=colnames(MEs),
          ySymbols=colnames(MEs),colorLabels=FALSE,yColorWidth=0.01,colors=col1,
          textMatrix=text1,setStdMargins=FALSE,cex.text=0.5,cex.lab.x=0.6,
          cex.lab.y=0.6,zlim=c(-1,1),main=main1,plotLegend=FALSE)
        par(oldpar)
        invisible(dev.off())
    }
    f2 = paste(f1,"_trait_MM2",sep="")
    text1 = paste(signif(as.matrix(modTS_2),2),"\n(",
      signif(as.matrix(modTSPval_2),1),")",sep="")
    dim(text1) = dim(modTS_2)
    col1 = colorRampPalette(c("red","white","green"))(n=10)
    main1 = "Module-trait relationships"
    for(i1 in c("eps","tiff")){
        if(i1 == "eps"){
            postscript(file=paste(f2,".eps",sep=""),width=7.0,height=8.0,
              colormodel="rgb",horizontal=FALSE,onefile=FALSE,paper="special",
              family="Helvetica")
        }else{
            tiff(file=paste(f2,".tif",sep=""),units="in",width=7.0,height=8.0,
              res=300,compression="lzw",family="sans")
        }
        oldpar = par(mar=c(4.0,3.8,3,0.2),cex.lab=0.8,cex.main=1.0)
        labeledHeatmap(Matrix=modTS_2,xLabels=traits,yLabels=colnames(MEs),
          ySymbols=colnames(MEs),colorLabels=FALSE,yColorWidth=0.01,colors=col1,
          textMatrix=text1,setStdMargins=FALSE,cex.text=0.5,cex.lab.x=0.6,
          cex.lab.y=0.6,zlim=c(-1,1),main=main1,plotLegend=FALSE)
        par(oldpar)
        invisible(dev.off())
    }
    f2 = paste(f1,"_trait_MM3",sep="")
    s1 = vector("logical",length=ncol(tra1))
    for(i1 in 1:ncol(tra1)){
        s1[i1] = is.null(levels(tra1[,i1]))
    }
    d1 = cbind(MEs,tra1[,s1])           #remove categorical variables
    d1 = remove_na(d1,init=1,plotit=FALSE,verbose=FALSE)
    MET = orderMEs(d1)
    for(i1 in c("eps","tiff")){
        if(i1 == "eps"){
            postscript(file=paste(f2,".eps",sep=""),width=7.0,height=10.0,
              colormodel="rgb",horizontal=FALSE,onefile=FALSE,paper="special",
              family="Helvetica")
        }else{
            tiff(file=paste(f2,".tif",sep=""),units="in",width=7.0,height=10.0,
              res=300,compression="lzw",family="sans")
        }
        oldpar = par(mar=c(2.2,2.2,1.2,1.2),mgp=c(1.0,0.2,0.0),tcl=-0.2,cex=0.7,
          cex.lab=0.7,cex.axis=0.7,cex.main=0.8)
        plotEigengeneNetworks(MET,setLabels="",marDendro=c(0,4.8,1,4.6),
          marHeatmap=c(5,6,1,1),xLabelsAngle=90)
        par(oldpar)
        invisible(dev.off())
    }
    f2 = paste(f1,"_trait_MM4",sep="")
    s1 = vector("logical",length=ncol(tra1))
    for(i1 in 1:ncol(tra1)){
        s1[i1] = is.null(levels(tra1[,i1]))
    }
    d1 = cbind(MEs,tra1[,s1])           #remove categorical variables
    d1 = remove_na(d1,init=1,plotit=FALSE,verbose=FALSE)
    MET = orderMEs(d1)
    for(i1 in c("eps","tiff")){
        if(i1 == "eps"){
            postscript(file=paste(f2,".eps",sep=""),width=7.0,height=7.0,
              colormodel="rgb",horizontal=FALSE,onefile=FALSE,paper="special",
              family="Helvetica")
        }else{
            tiff(file=paste(f2,".tif",sep=""),units="in",width=7.0,height=7.0,
              res=300,compression="lzw",family="sans")
        }
        oldpar = par(mar=c(0.5,0.5,0.5,0.5),mgp=c(3,0.5,0),tcl=-0.2,cex=0.6)
        plotMEpairs(MET,cex.main=1.0,cex.labels=0.8,oma=c(1.2,1.2,3.1,1.2),
          main="")
        par(oldpar)
        invisible(dev.off())
    }
    for(i1 in traits){
        f2 = paste(f1,"_trait_TS_",i1,sep="")
        main1 = paste("\"",i1,"\",",sep="")
        for(i2 in c("eps","tiff")){
            if(i2 == "eps"){
                postscript(file=paste(f2,".eps",sep=""),width=7.0,height=7.0,
                  colormodel="rgb",horizontal=FALSE,onefile=FALSE,
                  paper="special",family="Helvetica")
            }else{
                tiff(file=paste(f2,".tif",sep=""),units="in",width=7.0,
                  height=7.0,res=300,compression="lzw",family="sans")
            }
            oldpar = par(mar=c(6.2,3.1,2.1,1.1),mgp=c(2,0.5,0),tcl=-0.01,las=2,
              cex=0.8,cex.lab=0.8,cex.axis=0.8,cex.main=0.9)
            plotModuleSignificance(varTS[,i1],modCol,boxplot=TRUE,main=main1)
            abline(h=c(-0.9,-0.8,-0.7,-0.5,0.5,0.7,0.8,0.9),col="red",lty=c(1,4,
              2,3,3,2,4,1))
            par(oldpar)
            invisible(dev.off())
        }
    }
    modules = apply(modTSPval,2,function(x) modNam[x <= thr2])
    if(class(modules)=="list"){
        modules = lapply(modules,function(x) x[x!="grey"]) #remove "grey"
    }else{
        modules = rbind(apply(modules,2,function(x) x[x!="grey"])) #remove "grey"
        modules = as.list(as.data.frame(modules,stringsAsFactors=FALSE))
    }
    for(i1 in unique(unlist(modules))){
        f2 = paste(f1,"_trait_mod1_",i1,sep="")
        title1 = paste("ME",i1,sep="")
        xlab1 = ""
        ylab1 = ""
        modVars = modCol == i1
        clabs = rownames(dat1)
        if(sum(modVars) <= maxs1){
            rlabs = colnames(dat1)[modVars]
        }else{
            rlabs = FALSE
        }
        ME=MEs[,paste("ME",i1,sep="")]
        for(i2 in c("eps","tiff")){
            if(i2 == "eps"){
                postscript(file=paste(f2,".eps",sep=""),width=7.0,height=10.0,
                  colormodel="rgb",horizontal=FALSE,onefile=FALSE,
                  paper="special",family="Helvetica")
            }else{
                tiff(file=paste(f2,".tif",sep=""),units="in",width=7.0,
                  height=10.0,res=300,compression="lzw",family="sans")
            }
            oldpar = par(mfrow=c(2,1),mgp=c(1.0,0.2,0),mar=c(0.0,4.5,1.2,1.3),
              oma=c(0.0,0.0,0.0,0.0),tck=0,cex.main=1.2,cex=0.8)
            plotMat(t(scale(dat1[,modVars])),nrgcols=30,rlabels=rlabs,
              clabels=FALSE)
            title(title1,line=0.25)
            par(mar=c(2.0,2.5,0.5,0.1),tck=-0.01,cex.lab=0.6,cex.axis=0.6,cex=1)
            barplot(ME,col=i1,main="",cex.main=2,ylab=ylab1,xlab=xlab1,
              names.arg=clabs,cex.names=0.5,las=2)
            par(oldpar)
            invisible(dev.off())
        }
    }
    for(i1 in traits){
        main1 = paste("Module membership vs. vars-",i1," significance\n",sep="")
        main2 = paste("Intramodular Connectivity vs. vars-",i1,
          " significance\n",sep="")
        for(i2 in modules[[i1]]){
            f2 = paste(f1,"_trait_mod2_",i1,"_",i2,sep="")
            modVars = modCol == i2
            d1 = abs(varMM[modVars,i2])
            d2 = varIC[modVars,1]
            d3 = varTS[modVars,i1]
            xlab1 = paste("MM in",i2,"module")
            xlab2 = paste("IC in",i2,"module")
            ylab1 = "Significance"
            if(labeltop){       
                n = min(ntop,length(y))
                dst1 = rowMeans(cbind(scale(d1),scale(d3)))
                s1 = order(dst1,decreasing=TRUE)
                s2 = s1[d1[s1] > thr3 & d3[s1] > thr4][1:n]
                labs1 = colnames(dat1)[modVars][s2]
                dst2 = rowMeans(cbind(scale(d2),scale(d3)))
                thr5 = quantile(d2,thr5)
                s1 = order(dst2,decreasing=TRUE)
                s2 = s1[d2[s1] >= thr5 & d3[s1] > thr4][1:n]
                labs2 = colnames(dat1)[modVars][s2]
            }
            for(i3 in c("eps","tiff")){
                if(i3 == "eps"){
                    postscript(file=paste(f2,".eps",sep=""),width=7,height=10,
                      colormodel="rgb",horizontal=FALSE,onefile=FALSE,
                      paper="special",family="Helvetica")
                }else{
                    tiff(file=paste(f2,".tif",sep=""),units="in",width=7,
                      height=10,res=300,compression="lzw",family="sans")
                }
                oldpar = par(mfrow=c(2,1),mar=c(3.8,3.8,2.8,1.0),mgp=c(2.5,1,0),
                  tcl=-0.2,cex=0.8,cex.lab=0.8,cex.axis=0.8,cex.main=0.9)
                plot(d1,d3,col=i2,xlab=xlab1,ylab=ylab1,main=main1,cex=1.0)
                if(labeltop){
                    abline(v=thr3,lty=2)
                    abline(h=thr4,lty=2)
                    if(length(s2) > 0){
                        text(d1[s2],d3[s2],labels=labs1,cex=0.8,col="black",
                          pos=2,offset=0.3)
                    }
                }
                plot(d2,d3,col=i2,xlab=xlab2,ylab=ylab1,main=main1,cex=1.0)
                if(labeltop){
                    abline(v=thr5,lty=2)
                    abline(h=thr4,lty=2)
                    if(length(s2) > 0){
                        text(d2[s2],d3[s2],labels=labs2,cex=0.8,col="black",
                          pos=2,offset=0.3)
                    }
                }
                par(oldpar)
                invisible(dev.off())
            }
        }
    }
}  

}
