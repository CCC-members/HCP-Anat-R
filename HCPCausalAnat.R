## This is implementation for the Causal Analysis on HCP data for paper
##Causal Effects of Cingulate Morphology on Executive functions in Healthy Young Adults
##Code is based on book and related code by
##Robins M. James, Miguel A. Hernán. 2020. 
##Foundations of Agnostic Statistics Causal Inference - What If. 
##Chapman & Hall/CRC.

###---Author: Fuleah A. Razzaq


library(ggplot2)
library(ggpubr)
library(ggbiplot)
library(WeightIt)
library(geepack)
library(EValue)
library(ipw)
library(ggpmisc)
rm(list = ls(all.names = TRUE)) 
#load data
data=("path")

#####  association ########
C=c("Gender","Age","Hand","Race","GMVol","Fid","Mid")
dataarr=data[,!names(data) %in% C]
colnames(dataarr)=c("ID",
                    "CardSort","Flanker","ListSort",
                    "ACC Area","ACC Thickness","ACC Volume",
                    "PCC Area","PCC Thickness","PCC Volume")
n=nrow(dataarr)
score1=subset(dataarr,select=c(ID,Flanker))
score2=subset(dataarr,select=c(ID,CardSort))
score3=subset(dataarr,select=c(ID,ListSort))

type=rep("Flanker",n)
score1=cbind.data.frame(score1,type)

type=rep("Card Sort",n)
score2=cbind.data.frame(score2,type)

type=rep("List Sort",n)
score3=cbind.data.frame(score3,type)

colnames(score1)=colnames(score2)=colnames(score3)=c("ID","TestScore","TestType")

plotLM=function(dataarr,score,ytitle)
{
  dataplot=NULL
  annText=NULL
  for (i in 1:ncol(dataarr)) 
  {
    col=colnames(dataarr)[i]
    Size=(dataarr[,i])
    type=rep(col,n)
    temp=cbind.data.frame(score,Size,type)
    #temp=rbind.data.frame(temp,temp) #use when both
    dataplot=rbind.data.frame(dataplot,temp[,-c(4:6)])
    dataplot$Size=as.numeric(scale(dataarr[,i]))
    ss=summary(lm(TestScore~Size+Size:a+Size:as.factor(r),data=temp))
    ss=paste("p = ", signif(ss$coefficients[2,4], digits = 4), sep = "")
    
    ann_textTMP = data.frame(TestScore = 75,Size = 0,lab = ss,
                             TestType=score$TestType[1],
                           type=factor(col,
                                       levels=c("ACC Area","ACC Thickness","ACC Volume",
                                                "PCC Area", "PCC Thickness","PCC Volume"
                                       )))
    annText=rbind.data.frame(annText,ann_textTMP)
    }
  
  dataplot$typef=as.factor(dataplot$type)
  my.formula = y ~ x
  
  p=ggplot(transform(dataplot,
                   type=factor(type,
                               levels=c("ACC Area","ACC Thickness","ACC Volume",
                                       "PCC Area", "PCC Thickness","PCC Volume"
                                        ))),
         aes(y= TestScore, x=Size, color = "blue"),show.legend = FALSE) +
    #geom_line(color="blue",alpha = 0.3) +
    geom_point(color="blue",alpha = 0.3) +
    facet_wrap(type~., scales = "fixed") +
    
    geom_smooth(method = "lm", formula = my.formula, se = F,color="red",alpha=.9)+
   
     ylim(70,150)+

     labs(x="Anatomical measurements",y=ytitle)+
    
    theme(panel.background = element_rect(fill = "white", 
                                          colour = "grey50"))
  
  p=p + geom_text(
    data    = annText,
    mapping = aes(x = -Inf, y = -Inf,label = lab),
    hjust   = -0.2,
    vjust   = -1,show.legend = FALSE,color="black")
  p
  return (p)
}

dataarr=dataarr[,-(1:4)]#exlude id and Y

score=score1#for flanker only
ytitle="Flanker"
plotLM(dataarr,score,ytitle)

score=score2 #for card sort only
ytitle="Card Sort"
plotLM(dataarr,score,ytitle)

score=score3 #for list sort only
ytitle="List Sort"
plotLM(dataarr,score,ytitle)

######### causal analysis ########
X=c("ACC_Area","PCC_Area","ACC_Thick","PCC_Thick")
C=c("Gender","Age","Hand","Race","GMVol","Fid","Mid")
Y=c("CardSort","Flanker","ListSort")

f=1
ressummc=NULL
Wtc=NULL

for (x in X)
{ 
  wdata=data[,c(x,C)]
  weightFormulaD=as.formula(paste0(x,"~", paste(C, collapse="+")))
  weightFormulaN=as.formula(paste0(x,"~1"))
  
  # estimation of denominator of ip weights
  denom = lm(weightFormulaD,data = wdata)
  pdenom = predict(denom, type = "response")
  dens.den = dnorm(wdata[,x],pdenom,summary(denom)$sigma)
  
  # estimation of numerator of ip weights
  numer=lm(weightFormulaN, data = wdata)
  pnumer = predict(numer, type = "response")
  dens.num = dnorm(wdata[,x],pnumer, 
                   summary(numer)$sigma)
  
  #weights
  W= dens.num / dens.den
  W[W>26]=26
  Wtc[[x]]=cbind.data.frame(id=data$ID,W)
  for (y in Y)
  {
    fdata=data[,c("ID",x,y,C)]
    fdata$wt=W
    H1=as.formula(paste0(y,"~",x,"+",x,":Age+",x,":as.factor(Race)"))
    
    msm = geeglm(H1,data = fdata,weights = wt,id = ID,
                 corstr = "independence")
    yout=y
    n=row.names(summary(msm)$coefficients)
    xout=x
    p=summary(msm)$coefficients[,4] 
    
    beta = coef(msm)
    SE = coef(summary(msm))[, 2]
    lcl = beta - qnorm(0.975) * SE
    ucl = beta + qnorm(0.975) * SE
    
    tres=cbind.data.frame(xout,yout,p,beta, lcl, ucl)
    
   #evalue from msm

    e=evalues.OLS(est = tres[2,"beta"],
            se = SE[2],
            sd = sd(fdata[,y]),
            delta = mean(fdata[,x])
            )
    
    tresen=cbind.data.frame(tres[2,],
                            e['RR',1], e['E-values',1],mean(fdata[,x]))
    ressummc=rbind.data.frame(ressummc,tresen)
    
    f=f+1
  }
}
sigres=ressummc
colnames(sigres)=c("X","Y","p-value","Point Est.",
                   "95%CI-Lower","95%CI-Upper","RR","E-Value","delta")
sigres=sigres[order(sigres$Y),]#arrange by Y
sigres

#####IPWT plot#####
dplotdata=NULL
W=Wtc
n=nrow(W[[1]])
for (i in 1:length(W))
{
  w=W[[i]][2]
  pdata=cbind.data.frame(log(w),rep(names(W)[i],n),rep(0,n))
  dplotdata=rbind.data.frame(dplotdata,pdata)
}
colnames(dplotdata)=c("Weights","Area","y")
mu = ddply(dplotdata, "Area", summarise, grp.mean=mean(Weights))

ggplot(dplotdata, aes(x=Weights, color=Area)) +
  geom_density()+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=Area),
             linetype="dashed")+
  labs(x="log (IPT Weights) ", y="Density") +
  theme(axis.text=element_text(size=6),
        axis.title=element_text(size=10,face="bold"),
        panel.background = element_rect(fill = "white", 
                                        colour = "grey50"))





