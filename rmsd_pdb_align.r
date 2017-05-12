######################################
  #####RMSD cluster script#####
######################################
###written by CCSOSA 2015_01_27

require(bio3d);require(dendextend);library(seqinr);require(adegenet)



###list files###
  s_names_models<-list.files("E:/pdi_2016_09_23",pattern = ".pdb$",full.names = T,recursive = F)
snames<-list.files("E:/pdi_2016_09_23",pattern = ".pdb$",full.names = F,recursive = F)
snames<-sub(".pdb","",snames)
#s_names_models<-paste0("E:/pdi_2016_09_23",snames,"/","model","/","01","/","model.pdb")
#snames_2<-substr(snames, 1, nchar(snames)-11) # could run if you use the date into the folder names 
s_model<-lapply(s_names_models,read.pdb)

#alignment using MUSCLE Algorithm, it is a good idea filter the alignment in order to get a better one
s_aln_seq<-pdbaln(s_names_models,exefile="E:/SOFTWARE/Seaview/seaview4/muscle.exe")##exefile parameter must be changed accordly the muscle path
#s_aln_seq$id<-s_names_models




fit<-pdbfit(s_aln_seq,outpath ="E:/pdi_2016_09_23" )
#ss<-pdbaln(s_names_models)


###get distances matrix####
inds <- gap.inspect(s_aln_seq$xyz)

rd <- rmsd(a=s_aln_seq,fit=T,a.inds = inds$f.inds)
row.names(rd)<-snames
colnames(rd)<-snames

hc.rd <- hclust(as.dist(rd),method = "ward.D") 
#hc.rd$labels<-snames
plot(hc.rd)


#hclustplot(hc.rd, k=3,cex=0.8, ylab="RMSD (Å)", 
 #          main="RMSD Cluster Dendrogram", fillbox=FALSE)



###PLOT DENDOGRAM####
dend <- as.dendrogram(hc.rd)
plot_horiz.dendrogram(dend,type="rectangle",side=F,center=T,text_pos=3,xlab="RMSD (Å)",main="RMSD PDI Proteins Cluster Dendrogram",nodePar = list(cex=0.07),dLeaf= 0,xlim = c(30,-2))
abline(v=7,col="red")
legend("topleft",legend="RMSD: 7Å",col="red",lty=1,cex = 0.6)   


#################

mase.res <- read.alignment(file="E:/pdi_2016_09_23/result.mase", format = "mase")
x <- alignment2genind(mase.res,
                      pop = c("Eimeria tenella",
                                       "Toxoplasma gondii",
                                       "Sarcocystis neurona",
                                       "Sarcocystis neurona",
                                       "Neospora caninum",
                                       "Toxoplasma gondii",
                                       "Hammondia hammondi",
                                       "Sarcocystis neurona",
                                       "Toxoplasma gondii",
                                       "Neospora caninum",
                                       "Hammondia hammondi",
                                       "Sarcocystis neurona",
                                       "Neospora caninum",
                                       "Toxoplasma gondii",
                                       "Hammondia hammondi",
                                       "Neospora caninum",
                                       "Toxoplasma gondii",
                                       "Hammondia hammondi",
                                       "Neospora caninum",
                                       "Eimeria tenella",
                                       "Toxoplasma gondii",
                                       "Sarcocystis neurona",
                                       "Sarcocystis neurona",
                                       "Eimeria tenella",
                                       "Sarcocystis neurona",
                                       "Neospora caninum",
                                       "Toxoplasma gondii",
                                       "Hammondia hammondi",
                                       "Toxoplasma gondii",
                                       "Eimeria tenella",
                                       "Eimeria tenella",
                                       "Eimeria tenella",
                                       "Sarcocystis neurona",
                                       "Neospora caninum",
                                       "Toxoplasma gondii",
                                       "Hammondia hammondi",
                                       "Eimeria tenella",
                                       "Sarcocystis neurona",
                                       "Neospora caninum",
                                       "Toxoplasma gondii",
                                       "Hammondia hammondi",
                                       "Eimeria tenella",
                                       "Eimeria tenella",
                                       "Sarcocystis neurona",
                                       "Neospora caninum",
                                       "Toxoplasma gondii",
                                       "Hammondia hammondi",
                                       "Eimeria tenella",
                                       "Sarcocystis neurona",
                                       "Eimeria tenella",
                                       "Neospora caninum",
                                       "Toxoplasma gondii",
                                       "Hammondia hammondi",
                                       "Sarcocystis neurona",
                                       "Neospora caninum",
                                       "Toxoplasma gondii",
                                       "Hammondia hammondi",
                                       "Neospora caninum",
                                       "Toxoplasma gondii",
                                       "Hammondia hammondi",
                                       "Sarcocystis neurona",
                                       "Neospora caninum",
                                       "Toxoplasma gondii",
                                       "Hammondia hammondi",
                                       "Neospora caninum",
                                       "Toxoplasma gondii",
                                       "Hammondia hammondi",
                                       "Sarcocystis neurona",
                                       "Neospora caninum",
                                       "Toxoplasma gondii",
                                       "Hammondia hammondi",
                                       "Neospora caninum",
                                       "Toxoplasma gondii",
                                       "Hammondia hammondi",
                                       "Sarcocystis neurona",
                                       "Toxoplasma gondii",
                                       "Hammondia hammondi",
                                       "Neospora caninum",
                                       "Eimeria tenella",
                                        "Eimeria tenella"))
tabAA <- genind2df(x)


D <- as.matrix(dist(tab(x)),"euclidean")
D<-D[1:74,1:74]
# library(ape)
# tre <- nj(D)
# par(xpd=TRUE)
# plot(tre, type="unrooted", edge.w=2)
# edgelabels(tex=round(tre$edge.length,1), bg=rgb(.8,.8,1,.8))
pco1 <- dudi.pco(as.dist(D), scannf=FALSE,nf=2)
# s.label(pco1$li*1.1, clab=0, pch="")
# textplot(pco1$li[,1], pco1$li[,2], words=rownames(pco1$li),
#          cex=1.4, new=FALSE, xpd=TRUE)
# title("Principal Coordinate Analysis\n-based on proteic distances-")

D

