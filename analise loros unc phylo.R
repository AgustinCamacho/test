#packages
require(geiger)
require(phytools)
require(nlme)
require(nortest)
require(plyr)

####PREPARE RECEIVING DFS
setwd("C:/Users/Agus/Dropbox/Agus_Erica/dados loros")
#setwd("~/My Dropbox/Agus_Erica/dados loros/")
# read data and tree, and prepare their names
# bayesian tree, generate the parsymony too.
# to save list of trees from mesquite: taxa a&trees save trees from file...save stored trees give extension .newick
mparams=data.frame(matrix(nrow=0,ncol=10))
sizeparams=data.frame(matrix(nrow=0,ncol=10))
shapeparams=data.frame(matrix(nrow=0,ncol=10))
evofit=data.frame(matrix(nrow=0,ncol=3))
pPCA=data.frame(matrix(nrow=0,ncol=6))

colnames(mparams)=c("Value","Std.Error","t-value","p-value","model","psignal","psignalvar", "normp","factors","tree")
write.table(mparams,"mparams.csv",append = F, sep=",")
colnames(shapeparams)=c("Value","Std.Error","t-value","p-value","model","psignal","psignalvar","normp","factors","tree")
write.table(shapeparams,"shapeparams.csv",append = F, sep=",")
colnames(pPCA)=c("PC1","PC2","lambda","tree","variables")
write.table(pPCA,"pPCA.csv",append = F, sep=",")
colnames(evofit)=c("Brfit","OUfit","tree")
write.table(evofit,"evofit.csv",append = F, sep=",")




for (i in 1:10000)
{
  system.time(
print(i)
n=read.csv("loro morpho geo.csv")
file<-"parrotstrees.newick"
text<-scan(file,what="character",skip=as.character(i),nlines = 1)
tree<-read.tree(text=text)
tree$tip.label=gsub("_"," ",tree$tip.label)
tree$tip.label=tolower(tree$tip.label)
tree=drop.tip(tree,setdiff(tree$tip.label,n$Species))
rownames(n)=n$Species
n=n[tree$tip.label,]

# transform latitude to abs value
n$cent_Latitude=abs(n$cent_Latitude)
n$cent_Longitude=abs(n$cent_Longitude)
n=na.omit(n)# necessary for PCA takes out two species:leptosittaca branickii,prioniturus leuconensis,psittacula columboides
tree=drop.tip(tree,setdiff(tree$tip.label,n$Species))# adjust tree after exceprting spp
str(n)
###################################################################
# phylogenetic PCA reduction
m=n[,13:(ncol(n)-4)]
rownames(m)=n$Species
pp=phyl.pca(tree,m,method="lambda", mode="corr")
rownames(pp$L)
pp$L# reasonable loadings with centered data and correlations, very low loadings with uncentered and covariance
pp$lambda# high lambda, thus conservated morphology, good for ancestral reconstruction
n$size=pp$S[,"PC1"]
n$pwingbilltar=pp$S[,"PC2"]
n$Pptail=pp$S[,"PC3"]# very low loadings most of times
pPCA=rbind(pPCA,cbind(pp$L[,1:3],paste("tree",i)))  
n$Region=as.factor(n$Region)
levels(n$Region)=c("island","continent")
  
#################################################
# test evolutionary model

x=n$Area_Km2
          names(x)=tree$tip.label
          BMfit <- geiger::fitContinuous(tree, x, model="BM");BMfit
          OUfit <- fitContinuous(tree, x, model="OU");OUfit

tfit=cbind(BMfit$opt$aicc,OUfit$opt$aicc,paste("tree",i))
x=n$cent_Latitude
          names(x)=tree$tip.label
          BMfit <- fitContinuous(tree, x, model="BM");BMfit
          OUfit <- fitContinuous(tree, x, model="OU");OUfit
tfit=rbind(tfit,cbind(BMfit$opt$aicc,OUfit$opt$aicc,paste("tree",i)))
x=n$pwingbilltar
          names(x)=tree$tip.label
          BMfit <- fitContinuous(tree, x, model="BM");BMfit
          OUfit <- fitContinuous(tree, x, model="OU");OUfit
tfit=rbind(tfit,cbind(BMfit$opt$aicc,OUfit$opt$aicc,paste("tree",i)))
x=n$ptail
          names(x)=tree$tip.label
          BMfit <- fitContinuous(tree, x, model="BM");BMfit
          OUfit <- fitContinuous(tree, x, model="OU");OUfit
tfit=rbind(tfit,cbind(BMfit$opt$aicc,OUfit$opt$aicc,paste("tree",i)))
BROU=sum(as.numeric(as.character(tfit[,1])))>sum(as.numeric(as.character(tfit[,2])))
evofit=rbind(evofit,tfit)
write.table(evofit,"evofit.csv",append = T,col.names=F)
  
# ou always bvetter, even with different transformations(sc,log) 
#################################################
#tests

# keep filogenies fixed because they often give error when fitting also the alpha or gamma value
if
  (BROU==FALSE)
  {
  m=gls(sqrt(Area_Km2)~Region+ptail+pwingbilltar+size,method="ML", corMartins(1, tree, fixed=F),data = n)
  summ <- summary(m)
  alpha=attributes(m$apVar)$Pars[1];names(alpha)=NULL
  sigma=attributes(m$apVar)$Pars[2];names(sigma)=NULL
  mparams=rbind(mparams,cbind(summ$tTable,"OU",alpha,sigma,lillie.test(chol(solve(vcv(tree)))%*%residuals(m))$p.value,rownames(summ$tTable),i)) # store results
  rm(summ,m)
}

if
  (BROU==TRUE)
  {
  m=gls(sqrt(Area_Km2)~minlat+maxlat+Region+ptail+pwingbilltar+size,method="ML", corPagel(1, tree, fixed=T),data = n)
  summ <- summary(m)
  lambda=attributes(m$apVar)$Pars[1];names(lambda)=NULL
  sigma=attributes(m$apVar)$Pars[2];names(sigma)=NULL
  mparams=rbind(mparams,cbind(summ$tTable,"OU",lambda,sigma,lillie.test(chol(solve(vcv(tree)))%*%residuals(m))$p.value,rownames(summ$tTable),i)) # store results
# lambda and sigma do not have sense when phylo signal is fixed because
  # store results
  rm(summ,m)
}

)
}#end of for
colnames(mparams)=c("Value","Std.Error","t-value","p-value","model","norm.p")
colnames(evofit)=c("Brfit","OUfit","turn")
head(mparams)
mparams$parameter=rep(c("(Intercept)","cent_Latitude","Regioncontinent","ptail","pwingbilltar","size"),100)

write.table(evofit,"evofit.csv",row.names=F,sep=",")


# total variance estimated by

# variancia total= variancia do coeficiente + a media da variancia
varBeta<-var(mparams[mparams$parameter=="Regioncontinent","Value"])+mean(mparams[mparams$parameter=="Regioncontinent","Std.Error"])
# t value
t.beta<-mean(mparams[mparams$parameter=="Regioncontinent","Value"])/sqrt(varBeta)

# p value
P.beta<-2*pt(abs(t.beta),df=length(trees[[1]]$tip)-2, lower.tail=FALSE)
hist(P.beta)


mparams[mparams$parameter=="Regioncontinent","Std.Error"]

par(mfrow=c(2,3))
hist(mparams[mparams$parameter=="Regioncontinent","p-value"],main="Region", xlab="p-value")
hist(mparams[mparams$parameter=="cent_Latitude","p-value"],main="latitude", xlab="p-value")
hist(as.numeric(mparams[mparams$parameter=="ptail","p-value"]),main="RTL", xlab="p-value")
hist(mparams[mparams$parameter=="Regioncontinent","Value"],main="Region", xlab="slope")
hist(mparams[mparams$parameter=="cent_Latitude","Value"],main="Latitude", xlab="slope")
hist(mparams[mparams$parameter=="ptail","Value"],main="RTL", xlab="Value")
# one more histogram for normality assumption.

# normality test controlling for phylogeny

#  passed the normality assumption

#########################################
#visualize


pdf(file='C:/Users/Agus/Dropbox/Agus_Erica/dados loros/PANELA.pdf',
width=3, height=3)

## panel A
n$color=n$Region
levels(n$color)=c("black","blue")
plot(n$cent_Latitude,n$Area_Km2/1000000,
            xlab="Latitude (decimal)",ylab="GRS (Km2e06)",
            pch=19, cex=0.8,
            col=n$color)
dev.off()

#panel B
pdf(file='C:/Users/Agus/Dropbox/Agus_Erica/dados loros/PANELB.pdf',
    width=3, height=3)
plot(n$Region,n$Area_Km2/1000000,
          ylab = NA,
         col=n$color)
     box(col="black",lwd=1.5)
     mtext(side = 1, "", line = 2)
     mtext(side = 2, "GRS KM2e06", line = 2)
     text(0.6,7,"A")
box(col="black",lwd=1.5)
dev.off()
# Panel C
pdf(file='C:/Users/Agus/Dropbox/Agus_Erica/dados loros/PANELC.pdf',
    width=3, height=3)
plot(n$ptail,n$Area_Km2/1000000,
     xlab="RTL",ylab="GRS (Km2e06)",
     pch=19, cex=0.8,col=n$color)
dev.off()

# Prediction of tail proportion
m=gls(ptail~minlat+maxlat+Region,method="ML", corPagel(1, tree),data = n)
summary(m)
lillie.test(chol(solve(vcv(tree)))%*%residuals(m))


# summarizing PCA results across trees
p=read.csv("PCA.csv")
s=ddply(p,.(variables), summarize,  
          quantile(PC1,probs=c(0.025)),
          quantile(PC1,probs=c(0.975)),
          quantile(PC2,probs=c(0.025)),
          quantile(PC2,probs=c(0.975)),
          quantile(PC3,probs=c(0.025)),
          quantile(PC3,probs=c(0.975))
        )
s
write.table(s,"table PCA loro.csv",row.names=F,sep=",")
##summarizing evofit results
#summarizing GLS results

g=read.csv("mparamsloro.csv")
colnames(g)
class(g$p.value)
s=ddply(p,.(Parameters), summarize,  
        quantile(value,probs=c(0.025)),
        quantile(PC1,probs=c(0.975)),
        quantile(PC2,probs=c(0.025)),
        quantile(PC2,probs=c(0.975)),
        quantile(PC3,probs=c(0.025)),
        quantile(PC3,probs=c(0.975))
)
s








g=read.csv("mparamsloro.csv")
str(g)
colnames(g)
# liams, way
s=ddply(g,.(Parameters), summarize,
        val=mean(Value),
        varval=var(Value)+mean(Std.Error),
        tval=mean(t.value),
        vartval=var(t.value)+mean(Std.Error),
        quantile(p.value,probs=c(0.025)),
        quantile(p.value,probs=c(0.975)),
        quantile(norm.p,probs=c(0.025)),
        quantile(norm.p,probs=c(0.975))
        
)
write.csv(s,"table mparamsloro.csv",row.names=F,sep=",")

mean(g$norm.p>=0.05)

