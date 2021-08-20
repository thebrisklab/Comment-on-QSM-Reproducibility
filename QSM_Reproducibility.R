####################
# Ben Risk
# Simulations to examine the single-visit variance
# components models from Cogswell et al 2021
#
# This code supports the commentary:
# On the reproducibility of Quantitative Susceptibility Mapping and its potential as a clinical biomarker: a comment on Cogswell et al. 2021


library(lme4)
library(MASS)
library(ggplot2)
library(gridExtra)
set.seed(1234)


# QSM intrascanner at 3T:
# 

# Here, we use the variances reported in Tables 2 and 5 of
# Feng, X., Deistung, A., & Reichenbach, J. R. (2018). Quantitative susceptibility mapping (QSM) and R2* in the human brain at 3 T: Evaluation of intra-scanner repeatability. Zeitschrift für Medizinische Physik, 28(1), 36-48.

# Use fWM normalized sds:

# 1. Accumbens
# 2. Caudate
# 3. Globus Pallidus
# 4. Hippocampus
# 5. Putamen
# 6. Thalamus
sigmae1 = 6.4397
sigmae2 = 2.9023
sigmae3 = 5.1835
sigmae4 = 3.5113
sigmae5 = 2.7242
sigmae6 = 2.2921


sigmab1 = 17.5220
sigmab2 = 6.2399
sigmab3 = 31.8223
sigmab4 = 8.6990
sigmab5 = 12.1866
sigmab6 = 8.6665

#true ratio of sd:
sigmab1/sigmae1 #reported in manuscript
sigmab2/sigmae2
sigmab3/sigmae3
sigmab4/sigmae4
sigmab5/sigmae5
sigmab6/sigmae6 #reported in manuscript

# true ICC: 
sigmab1^2/(sigmab1^2+sigmae1^2)
sigmab2^2/(sigmab2^2+sigmae2^2)
sigmab3^2/(sigmab3^2+sigmae3^2)
sigmab4^2/(sigmab4^2+sigmae4^2)
sigmab5^2/(sigmab5^2+sigmae5^2)
sigmab6^2/(sigmab6^2+sigmae6^2)


# means are from the first scan
mu1 = 39.42
mu2 = 56.25
mu3 = 131.08
mu4 = 14.58
mu5 = 52.65
mu6 = 19.94

# Simulate data:
var(c(mu1,mu2,mu3,mu4,mu5,mu6))

n = 1000
#n = 10000

R = 3
id = c(1:n)%x%rep(1,R)

temp = diag(c(sigmab1,sigmab2,sigmab3,sigmab4,sigmab5,sigmab6))

# Correlation between random effects in different regions varies.
# In practice, this will be different for each pair of brain regions.
# For 6 brain regions, this involves 6*5/2 = 15 correlations.
# We simplify this setting by setting all correlations
# 15 correlations equal to the same value. 
# We consider two correlations: 1 for all regions; 0.5 for all regions:
rho=c(1,0.5)
varcomp.mismodel.list=list(rho1=NULL,rho0.5=NULL)
varcomp.modelcorrect.list.r1 = list(rho1=NULL,rho0.5=NULL)
varcomp.modelcorrect.list.r6 = list(rho1=NULL,rho0.5=NULL)


for (i in 1:length(rho)) {
  corrMat = matrix(rho[i],6,6)
  diag(corrMat) = 1
  Sigma = temp%*%corrMat%*%temp
  
  b12 = mvrnorm(n=n, mu=c(0,0,0,0,0,0),Sigma=Sigma)
  
  yi1r= b12[,1]%x%rep(1,R)+mu1+rnorm(n=n*R,sd=sigmae1)
  yi2r= b12[,2]%x%rep(1,R)+mu2+rnorm(n=n*R,sd=sigmae2)
  yi3r= b12[,3]%x%rep(1,R)+mu3+rnorm(n=n*R,sd=sigmae3)
  yi4r= b12[,4]%x%rep(1,R)+mu4+rnorm(n=n*R,sd=sigmae4)
  yi5r= b12[,5]%x%rep(1,R)+mu5+rnorm(n=n*R,sd=sigmae5)
  yi6r= b12[,6]%x%rep(1,R)+mu6+rnorm(n=n*R,sd=sigmae6)
  
  
  sim.data = data.frame(y=c(yi1r,yi2r,yi3r,yi4r,yi5r,yi6r),id=c(id,id,id,id,id,id),region=c(c(1:6)%x%rep(1,n*R)),replicate = rep(1,n*R)%x%c(1:R))
  
  #mean.data = sim.data %>% group_by(id,region) %>% summarize_at(vars(y),list(ybar=mean))
  
  # check how well we estimate the true variance components.
  # These results were added to the revised manuscript:
  varcomp.modelcorrect.list.r1[[i]] = as.data.frame(VarCorr(lmer(yi1r~(1|id))))
  varcomp.modelcorrect.list.r6[[i]] = as.data.frame(VarCorr(lmer(yi6r~(1|id))))
 
  sub.data = sim.data[sim.data$replicate==1,]
  mismodel = lmer(y~(1|id)+(1|region),data=sim.data[sim.data$replicate==1,])
  summary(mismodel)
  
  varcomp.mismodel = as.data.frame(VarCorr(mismodel))
  
  varcomp.mismodel.list[[i]] = varcomp.mismodel
}


# ratio of "between" to "within":
varcomp.mismodel.list


# these numbers appear in the manuscript:
# ratio of inter-participant to "error" variance in mis-specified model:
# rho=1
varcomp.mismodel.list[[1]][1,5]/varcomp.mismodel.list[[1]][3,5]

# rho=0.5
varcomp.mismodel.list[[2]][1,5]/varcomp.mismodel.list[[2]][3,5]




## Using the correct model:
varcomp.modelcorrect.list.r1[[1]][1,5]/varcomp.modelcorrect.list.r1[[1]][2,5]
varcomp.modelcorrect.list.r1[[2]][1,5]/varcomp.modelcorrect.list.r1[[2]][2,5]

varcomp.modelcorrect.list.r6[[1]][1,5]/varcomp.modelcorrect.list.r6[[1]][2,5]
varcomp.modelcorrect.list.r6[[2]][1,5]/varcomp.modelcorrect.list.r6[[2]][2,5]




########### Create plots:
# rho=1: 1.33; sigmab = 13.6; sigmae = 10.3

# rho=0.75: 0.96;
# rho=0.5 0.66; 9.4, 14.2
# rho=0: can't distinguish from error

df <- data.frame(SD=factor(rep(c("Participant", "Error"), each=10000)), Random.effect=c(rnorm(10000, mean=0, sd=sigmab1),rnorm(10000, mean=0, sd=sigmae1)))
p0<-ggplot(df, aes(x=Random.effect, fill=SD)) + geom_density(alpha=0.4,adjust=3)+theme(legend.position=c(0.8,0.7))+ggtitle('a) Truth: Region 1')+scale_fill_manual(values=c(Participant='yellow',Error='salmon'))


###############
df <- data.frame(SD=factor(rep(c("Participant", "Error"), each=10000)),Random.effect=c(rnorm(10000, mean=0, sd=sigmab6),rnorm(10000, mean=0, sd=sigmae6)))
# Use semi-transparent fill
p1<-ggplot(df, aes(x=Random.effect, fill=SD)) +
  geom_density(alpha=0.3,adjust=3)+theme(legend.position=c(0.8,0.7))+ggtitle('b) Truth: Region 6')+scale_fill_manual(values=c(Participant='yellow',Error='salmon'))


#######################
df <- data.frame(SD=factor(rep(c("Participant", "Error"), each=10000)), Random.effect=c(rnorm(10000, mean=0, sd=varcomp.mismodel.list$rho1[1,5]),rnorm(10000, mean=0, sd=varcomp.mismodel.list$rho1[3,5])))

p2<-ggplot(df, aes(x=Random.effect, fill=SD))+geom_density(alpha=0.4,adjust=3)+theme(legend.position=c(0.8,0.7))+ggtitle('c) Scenario 1: Misspecified')+scale_fill_manual(values=c(Participant='yellow',Error='salmon'))

#######################
df <- data.frame(SD=factor(rep(c("Participant", "Error"), each=10000)), Random.effect=c(rnorm(10000, mean=0, sd=varcomp.mismodel.list$rho0.5[1,5]),rnorm(10000, mean=0, sd=varcomp.mismodel.list$rho0.5[3,5])))

p3<-ggplot(df, aes(x=Random.effect, fill=SD))+geom_density(alpha=0.4,adjust=3)+theme(legend.position=c(0.8,0.7))+ggtitle('d) Scenario 2: Misspecified')+scale_fill_manual(values=c(Participant='yellow',Error='salmon'))


#pdf(file='~/Dropbox/Apps/Overleaf/QSM Reliability/Figure1.pdf',width=12,height=5)
grid.arrange(p0,p1,p2,p3,nrow=1)
#dev.off()
