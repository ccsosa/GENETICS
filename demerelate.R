wrequire(Demerelate);require(adegenet)

setwd("E:/Kathe_2015_02")
demerelpop<-read.table("./.txt",header=T,sep="\t")
demerelpop2<-read.table("./.txt",header=T,sep="\t",row.names = 1)

  Palma_ind<-genind(tab=demerelpop2[,2:21],pop=demerelpop2[,1],ploidy = as.integer(2))

grp <- find.clusters(Palma_ind, max.n.clust=10)

DPAL<-dapc(Palma_ind)

#,n.pca = 4,n.da = 2,scale = T)
oD<-optim.a.score(DPAL)
pred.sup <- predict.dapc(DPAL)

temp <- summary(DPAL)
 par(mar=c(4.5,7.5,1,1))

#  barplot(temp, xlab="% of reassignment to actual breed", horiz=TRUE)

demerelpop.sp <- split(demerelpop,demerelpop[,2])
compoplot(DPAL, posi=list(x=0,y=1.2),lab=row.names(DPAL))#

scatter.dapc(DPAL,scree.da=FALSE)
, bg="white", pch=20, cell=0, cstar=0, solid=.4,
              cex=3,clab=0, leg=TRUE)



########DEMERELATE######
dem.results <- Demerelate(demerelpop, value="rxy",iteration=1000,Fis=T,
                          file.output=T, object=TRUE, pairs=1000,NA.rm=T)




#empirical.result <- emp.calc(demerelpop.sp[[1]], value="Mxy",
                             ref.pop="NA")
fstat.results <- F.stat(demerelpop, iteration = 1000,
                        directory.name = "./",object=T,
                        out.name = "Summary_Fstat")

  Loci.results <- Loci.test(tab.pop=demerelpop.sp[[1]], object = TRUE,
                          value = "rxy", bt = 1000,file.output=T)
