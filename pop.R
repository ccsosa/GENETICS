require(polysat);require(ade4);require("scatterplot3d");require(adegenet);require(dendextend);require(ape);require(poppr);require(pegas);require(PopGenKit)




mati2<-read.table("")
#mati2_pop<-mati2[,4:47]
mati2_pop<-mati2[,4:39] ##COlUMNS WITH MARKERS TO ANALYZZE

mati3_pop<-mati2_pop

#for(i in 1:ncol(mati3_pop)){
  
  
 # mati3_pop[,i]<-factor(mati3_pop[,i])  
#}



locname<-colnames(mati3_pop)
locname<-sub(".1","",locname,ignore.case = F,fixed=T)
locname<-sub(".2","",locname,ignore.case = F,fixed=T)
locname<-sub(".3","",locname,ignore.case = F,fixed=T)
locname<-sub(".4","",locname,ignore.case = F,fixed=T)
locname<-sub(".5","",locname,ignore.case = F,fixed=T)
locname<-sub(".6","",locname,ignore.case = F,fixed=T)
locname<-unique(locname)


mati4<-as.matrix(mati3_pop)
PAdata <- new("genbinary", samples=row.names(mati4),loci=locname)
Genotypes(PAdata) <- mati4

Ploidies(PAdata)<-3

PopInfo(PAdata)<-mati2$Pop_number
PopNames(PAdata)<-as.character(unique(mati2$Pop_name))
Usatnts(PAdata)<-rep(x = 2,length(Usatnts(PAdata)))


#PAdata_P<-reformatPloidies(PAdata,output = "sample",na.rm = T)

PAdata_am<-genbinary.to.genambig(PAdata)
PAdata_am<- reformatPloidies(PAdata_am, output="sample")
PAdataFreq <- simpleFreq(PAdata,samples = Samples(PAdata),loci = Loci(PAdata))
PAdataFreq_Am<-simpleFreq(PAdata_am,samples = Samples(PAdata_am),loci = Loci(PAdata_am))

write.csv(PAdataFreq,"./Simple_Allele_freq.csv")
write.csv(PAdataFreq_Am,"./Simple_Allele_freq_am.csv")
                         #,loci=Loci(PAdata_am),samples = Samples(PAdata_am))


meanDist<-meandistance.matrix2(PAdata_am,loci = Loci(PAdata_am),freq = PAdataFreq,all.distances=F,distmetric = Bruvo.distance)
write.csv(meanDist,"E:/poliplo/Ind_distance.csv")
gM<-gendata.to.genind(PAdata_am,samples = Samples(PAdata_am),loci = Loci(PAdata_am))
gM2<-gendata.to.genind(PAdata_am)
gMpop<-genind2loci(gM)
gMpop2<-as.loci(gM)
gMpop2<-gMpop2[complete.cases(gMpop2),]

require(gdata)
s<-Loci(PAdata)

s2<-matrix(nrow = length(unique(mati2$Pop_name)),ncol=length(s))

colnames(s2)<-s
row.names(s2)<-unique(mati2$Pop_name)

for(i in 1:length(s)){
  for(j in 1:nrow(s2)){
x<-  matchcols(PAdataFreq_Am,with=s[i])
x2<-PAdataFreq_Am[j,x]^2
x3<-1-sum(x2)

s2[j,i]<-x3
  };rm(j)
};rm(i)

s2<-t(s2)
write.csv(s2,"./HExp_marker.csv",quote=F)

l_t<-locus_table(gM2, index = "simpson", lev = "allele", population = "ALL",information = TRUE)

ggppop<-genind2genpop(gM)

require(diveRsity)
basicStats(ggppop)
           

write.csv(l_t,"./locus_data.csv")
sLD<-LD(gMpop2,locus=c(as.integer(1),as.integer(14)),details = T)
require(LDcorSV)

LD_gM<-LD.Measures(gM$tab,data="G",na.presence = T,supinfo = T)


ddd<-dapc(gM,n.pca = 2,n.da = 2,var.contrib = T,scale = F)
lab<-pop(gM)
var<-mati2$Variedad
compoplot(ddd, posi=list(x=0,y=1.2),lab=row.names(ddd))#
#find.clusters(gM)
scatter.dapc(ddd,scree.da=T)

contrib <- loadingplot(ddd$var.contr, axis=1, thres=.05, lab.jitter=1)

PAdiversity<-alleleDiversity(PAdata_am)



#tr <- nj(as.dist(meanDi2st),X = )
hcD<-hclust(as.dist(meanDist))
plot(hcD)
M<-as.dendrogram(hcD)
colors_to_use <- as.numeric(mati2$Pop_name)
colors_to_use <- colors_to_use[order.dendrogram(M)]
labels_colors(M) <- colors_to_use
dend1 <- color_branches(M, k = 2)
M<-set(M, "labels_cex", 1)
plot_horiz.dendrogram(M,type="rectangle",side=F,center=T,text_pos=3,xlab="Bruvo distance",main="#############3",dLeaf= 0)
legend("topleft",legend=c("POP1","POP2"),col=c("red","black"),lty=1,cex = 0.65)   


pca <- cmdscale(meanDist, eig=TRUE)
pca$eig/sum(pca$eig)
mycol <- c("red", "green")
plot(pca$points[,1], pca$points[,2], col=mycol[PopInfo(PAdata_am)], main = "PCA with Bruvo distance")

locN<-unique(as.matrix(as.data.frame(strsplit(names(PAdataFreq), split = ".", fixed = TRUE),stringsAsFactors = FALSE))[1, ])


  SimFst <- calcFst(PAdataFreq)
SimFst_Am<- calcFst(PAdataFreq_Am,colnames(PAdataFreq_Am))
loci = c(MAOCEN08.1 MAOCEN08.3 MAOCEN04.1 MAOCEN04.2 MAOCEN04.3 MAOCEN17.1))
PAdataFreq_Am
#, pops = row.names(PAdataFreq), loci = unique(as.matrix(
    #as.data.frame(strsplit(names(PAdataFreq), split = ".", fixed = TRUE),
     #             stringsAsFactors = FALSE))[1, ]))




genDiv<-genotypeDiversity(PAdata_am,index=Shannon)
genDiv2<-genotypeDiversity(PAdata_am,index=Simpson)
write.csv(genDiv,"./Shannon_diversity_p_locus.csv")
write.csv(genDiv2,"./Simpson_diversity_p_locus.csv")


Pr<-as.genclone(gM)

Pinf_rc <- recode_polyploids(Pr)



summary(Pinf_rc)


poppr_table<-poppr(Pinf_rc,sample=999,hist=T,legend=T)
colnames(poppr_table)<-c("Population","n_individues","Multilocus genotypes_n","Multilocus genotypes_exp","Standard_error","Shannon-Weiner","Stoddard","Expected heterozygosity","Evenness"," Index of Association","p-value_IA","Standardized Index of Association","p-value for rbarD","File")

write.csv(poppr_table,"./summmary_table.csv")
inf<-info_table(Pinf_rc, plot = TRUE)

#fstat(gM, pop=NULL, fstonly=FALSE)


###########################################
     ####MULTIVARIATE APPROACH####
###########################################
library(rgl);require(ggplot2)

mat2<-read.table("./##########.txt",header = T,row.names=1,sep="\t")
mat2_pop<-mat2[,4:39]
mat3_pop<-mat2_pop

for(i in 1:ncol(mat3_pop)){
  
  
  mat3_pop
  mat3_pop[,i]<-factor(mat3_pop[,i])  
}

#xyz = MCA(mat_a,quali.sup = 1:2, graph = F,na.method = "NA")

xyz<-dudi.acm(mat3_pop,nf = 5,scannf=F)
xyz_in<-inertia.dudi(xyz,row.inertia = T,col.inertia = T)
s.label(xyz$li,cpoint = 0.1,csub = 0.2,clabel = 0.5,boxes = F)
summary(xyz)
boxplot(xyz)

xyz_in$col.rel
Loc<-rownames(xyz$co);loc1<-Loc
Loc<-sub(".0","",Loc,ignore.case = F,fixed=T)
Loc<-sub(".1","",Loc,ignore.case = F,fixed=T)
Loc<-sub(".2","",Loc,ignore.case = F,fixed=T)
Loc<-sub(".3","",Loc,ignore.case = F,fixed=T)
Loc<-sub(".4","",Loc,ignore.case = F,fixed=T)
Loc<-sub(".5","",Loc,ignore.case = F,fixed=T)
Loc<-sub(".6","",Loc,ignore.case = F,fixed=T)
Loc<-sub(".1","",Loc,ignore.case = F,fixed=T)

nam<-unique(Loc)
nam2<-c("A","B","C","D","E","F","G","H","I","J","K","L")
Loc2<-loc1

for(i in 1:length(nam)){
Loc2<-sub(nam[[i]],nam2[[i]],Loc2,ignore.case =F,fixed=F,useBytes=T)
}

row.names(xyz$co)<-Loc2



Loc_df<-cbind(loc1,Loc2,Loc)
colnames(Loc_df)<-c("original_name","code","Loci")
write.csv(Loc_df,"./code_graph_MCA.csv",row.names=F)





write.csv(xyz_in$col.abs/100,"./contribuciones_relativas_MCA.csv",row.names=T)
write.csv(xyz$li[,1:2],"./coords_ind_MCA.csv",row.names=T)
write.csv(xyz$co[,1:2],"./coords_var_contr_MCA.csv",row.names=T)

#####plot variables###

mca3_vars_df = data.frame(xyz$co, Loci = Loc)


ggplot(data = mca3_vars_df, 
       aes(x = Comp1, y = Comp2, label = rownames(mca3_vars_df))) +
  geom_hline(yintercept = 0, colour = "gray70") +
  geom_vline(xintercept = 0, colour = "gray70") +
  geom_text(aes(colour = Loci), size=2.4) +
  ggtitle("MCA plot of variables")


#####scatter variables

#xlab = "Dim 1 27.09%",ylab = "Dim 2 15.43%",zlab = "Dim 3 10.14%"
rR<-scatterplot3d(mca3_vars_df$Comp1,mca3_vars_df$Comp2,mca3_vars_df$Comp3,type = "p", box = T,xlab = "Dim 1 31.04%",ylab = "Dim 14.66%",zlab = "Dim 3 10.8%",highlight.3d = F,pch=20
                  , angle = 100)#10
rR.coords <- rR$xyz.convert(mca3_vars_df$Comp1,mca3_vars_df$Comp2,mca3_vars_df$Comp3)

rR$points3d(mca3_vars_df[,1:3][which(mca3_vars_df$Loci==nam[[1]]),],col = "gray9", type = "h", pch = 19)
rR$points3d(mca3_vars_df[,1:3][which(mca3_vars_df$Loci==nam[[2]]),],col = "forestgreen", type = "h", pch = 20)
rR$points3d(mca3_vars_df[,1:3][which(mca3_vars_df$Loci==nam[[3]]),],col = "burlywood4", type = "h", pch = 20)
rR$points3d(mca3_vars_df[,1:3][which(mca3_vars_df$Loci==nam[[4]]),],col = "red", type = "h", pch = 20)
rR$points3d(mca3_vars_df[,1:3][which(mca3_vars_df$Loci==nam[[5]]),],col = "green", type = "h", pch = 20)
rR$points3d(mca3_vars_df[,1:3][which(mca3_vars_df$Loci==nam[[6]]),],col = "gold", type = "h", pch = 20)
rR$points3d(mca3_vars_df[,1:3][which(mca3_vars_df$Loci==nam[[7]]),],col = "cornflowerblue", type = "h", pch = 20)
rR$points3d(mca3_vars_df[,1:3][which(mca3_vars_df$Loci==nam[[8]]),],col = "cyan2", type = "h", pch = 20)
rR$points3d(mca3_vars_df[,1:3][which(mca3_vars_df$Loci==nam[[9]]),],col = "darkorange", type = "h", pch = 20)
rR$points3d(mca3_vars_df[,1:3][which(mca3_vars_df$Loci==nam[[10]]),],col = "mediumspringgreen", type = "h", pch = 20)
#rR$points3d(mca3_vars_df[,1:3][which(mca3_vars_df$Loci==nam[[11]]),],col = "lightgoldenrod3", type = "h", pch = 20)
#rR$points3d(mca3_vars_df[,1:3][which(mca3_vars_df$Loci==nam[[12]]),],col = "dodgerblue3", type = "h", pch = 20)

text(rR.coords$x,rR.coords$y,labels = row.names(mca3_vars_df),pos=2,offset = 0.5,cex = 0.5)
#legend(rR$xyz.convert(-2.5, 0.7, 0.4), col =c("gray9","forestgreen","burlywood4","red","green","gold","cornflowerblue","cyan2","darkorange","mediumspringgreen","lightgoldenrod3","dodgerblue3"), yjust =0,legend = nam, cex = 0.5,lwd=1)

legend(rR$xyz.convert(-2.5, 0.5, 0.1), col =c("gray9","forestgreen","burlywood4","red","green","gold","cornflowerblue","cyan2","darkorange","mediumspringgreen"), yjust =0,legend = nam, cex = 0.4,lwd=1)


#####scatter ind


xpop<-cbind(mati2$Pop_name,xyz$li[,1:3])

plot(xyz$li[,1],xyz$li[,2],col=mati2$Pop,pch=15:16)
text(xyz$li[,1],xyz$li[,2],labels = row.names(mati2),offset = 0.2,cex = 0.5,adj = 1,pos=3)

#xlab = "Dim 1 27.09%",ylab = "Dim 2 15.43%",zlab = "Dim 3 10.14%"

#
rR<-scatterplot3d(xyz$li[,1],xyz$li[,2],xyz$li[,3],type = "p", box = T,xlab = "Dim 1 31.04%",ylab = "Dim 2 14.66%",zlab = "Dim 3 10.8%",highlight.3d = F,pch=20
                  , angle = 250)#250, 120,39
rR.coords <- rR$xyz.convert(xyz$li[,1],xyz$li[,2],xyz$li[,3])

#rR$points3d(mca3_vars_df[,1:3][which(mca3_vars_df$Loci=="MAOCEN08"),],col = "red", type = "h", pch = 19)
rR$points3d(xpop[,2:4][which(xpop[]=="Buenavista"),],col = "blue", type = "h", pch = 20)
rR$points3d(xpop[,2:4][which(xpop[]=="Montenegro"),],col = "red", type = "h", pch = 20)

text(rR.coords$x,rR.coords$y,labels = 1:48,pos=2,offset = 0.5,cex = 0.5)
legend(rR$xyz.convert(0.2, -0.2, 0.4), col = c("red","blue"), yjust =0,legend = c("Montenegro", "Buenavista"), cex = 0.6,lwd=1)
#legend(rR$xyz.convert(-0.9, 0.4, 0.5), col = c("red","blue"), yjust =0,legend = c("Montenegro", "Buenavista"), cex = 0.6,lwd=1)


s.class(xyz$li[,1:2],mat2$Pop_name,cellipse = 01,cstar = 1,csub = 0.5,col=c("red","blue"),xax = 1,yax = 2)


# rR$points3d


Populations<-mat2$Pop_name

xxx<-xyz$li
attach(xxx)
ind<-1:length(row.names(mat2))


sp<-ggplot(xxx,aes(x=Axis1, y=Axis2,colour=Populations,shape=Populations))#+ geom_point()
#scale_shape_manual(values=c(1, 2)) 
#scale_colour_brewer(palette="Set1" )
sp + geom_text(aes(label=ind),  size=3,hjust=0)
#sp + geom_point() + stat_density2d()
#,geom_point() +
 #        scale_shape_manual(values=c(1, 2)) +
  #       scale_colour_brewer(palette="Set1" ))
#p=qplot(xyz$li[,1],xyz$li[,2],color=Populations,xlab="27.09% Inertia",ylab="15.4% Inertia",main="Multiple Correspondence Analysis for plantain population in Quindío")

#p=geom_text(aes(label=row.names(mat2)),  vjust=1.5,  colour="black" )
#p + theme_bw() 

#co<-xyz$co[,1:2]
#sp<-ggplot(co,aes(x=Comp1, y=Comp2))#+ geom_point()
#scale_shape_manual(values=c(1, 2)) 
#scale_colour_brewer(palette="Set1" )
#sp + geom_text(aes(label=row.names(co)),  size=2,hjust=0)
  



  #xlab = "Dim 1 31.04%",ylab = "Dim 14.66%",zlab = "Dim 3 10.8%"
  interleave <- function(v1,  v2) as.vector(rbind(v1, v2))
  
  plot3d(xyz$li[,1],xyz$li[,2],xyz$li[,3],  axes=F,type="p" ,  size=0.75,  lit=FALSE,xlab = "" ,  ylab = "" ,  zlab = "" )
  segments3d(interleave(xyz$li[,1],  xyz$li[,1]),
             interleave(xyz$li[,2],  xyz$li[,2]),
             interleave(xyz$li[,3],  min(xyz$li[,3])),
             alpha = 0.4,  col = "blue" )
  axes3d(edges=c("x--" , "y+-" , "z--" ),
         ntick=6,  # Attempt 6 tick marks on each side
         cex=.75)
  mtext3d("Dim 1 31.04%" ,  edge="x--" ,  line=2)#27.09
  mtext3d("Dim 2 14.66%" ,  edge="y+-" ,  line=3)#15.43
  mtext3d("Dim 3 10.8%" ,  edge="z--" ,  line=3)#10.14
  rgl.bbox(color="grey50" ,  # grey60 surface and black text
           emission="grey50" ,  # emission color is grey50
           xlen=0,  ylen=0,  zlen=0) #
  rgl.material(color="black" )
  
  text3d(xyz$li[,1],xyz$li[,2],xyz$li[,3],texts = ind,font=1,cex=0.7,adj=c(1,1),color="red")
  
  rgl.snapshot('./3dplot.png' ,  fmt='png' )



