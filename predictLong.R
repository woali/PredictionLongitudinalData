####mailto: rphdstatistics@gmail.com
#Code by - Stanislaw Wanat, Alicja Wolny-Dominiak

library(copula)
library(CopulaRegression)
library(tweedie)
library(cplm)
library(gplots)
library(ggplot2)
library(gridExtra)

panel1 <- read.csv2(file="panel_2007.csv")
panel2 <- read.csv2(file="panel_2008.csv")
panel3 <- read.csv2(file="panel_2009.csv")
panel4 <- read.csv2(file="panel_2010.csv")

##Dealing with data
panel1$SEX<-ifelse(panel1$SEX==1,"M","F")
panel2$SEX<-ifelse(panel2$SEX==1,"M","F")
panel3$SEX<-ifelse(panel3$SEX==1,"M","F")
panel4$SEX<-ifelse(panel4$SEX==1,"M","F")

panel1$KIND_OF_PAYMENT<-ifelse(panel1$KIND_OF_PAYMENT=="cash","C","T")
panel2$KIND_OF_PAYMENT<-ifelse(panel2$KIND_OF_PAYMENT=="cash","C","T")
panel3$KIND_OF_PAYMENT<-ifelse(panel3$KIND_OF_PAYMENT=="cash","C","T")
panel4$KIND_OF_PAYMENT<-ifelse(panel4$KIND_OF_PAYMENT=="cash","C","T")

attach(panel4)
group <-interaction(SEX,ENGINE,KIND_OF_PAYMENT)
g.1 <-g.2<-group
for (j in 1:8) {g.2<-ifelse(g.1==levels(group)[j],j,g.2)}
group.nb <-as.factor(g.2)
legend <-data.frame(Number=levels(group.nb),Group=levels(group))
legend

##Exposure
Y1 <- panel1$CLAIM_AMOUNT/panel1$EXPOSURE
Y2 <- panel2$CLAIM_AMOUNT/panel2$EXPOSURE
Y3 <- panel3$CLAIM_AMOUNT/panel3$EXPOSURE
Y4 <- panel4$CLAIM_AMOUNT/panel4$EXPOSURE

##Fig1
m <-mean(c(mean(Y1), mean(Y2), mean(Y3), mean(Y4)))
par(mfrow=c(2,2))
hist(Y1[Y1>0], breaks=30, col = "lightgray", main="", xlab="Total claim amount >0, t=2007", cex.axis=1.3, cex.lab=1.3)
abline(v=m, lty=2)
hist(Y2[Y2>0], breaks=30, col = "lightgray", main="", xlab="Total claim amount >0, t=2008", cex.axis=1.3, cex.lab=1.3)
abline(v=m, lty=2)
hist(Y3[Y3>0], breaks=30, col = "lightgray", main="", xlab="Total claim amount >0, t=2009", cex.axis=1.3, cex.lab=1.3)
abline(v=m, lty=2)
hist(Y4[Y4>0], breaks=30, col = "lightgray", main="", xlab="Total claim amount >0, t=2010", cex.axis=1.3, cex.lab=1.3)
abline(v=m, lty=2)

##Fig2
par(mfrow=c(1,3))
plotmeans(Y4~panel4$SEX, xlab="Gender", ylab="Total claim amount, t=2010", cex.axis=1.3, cex.lab=1.3)
plotmeans(Y4~panel4$ENGINE, xlab="Engine", ylab="", cex.axis=1.3, cex.lab=1.3)
plotmeans(Y4~panel4$KIND_OF_PAYMENT, xlab="Kind of payment", ylab="", cex.axis=1.3, cex.lab=1.3)

##Estimation of Tweedie parameters
model.Y1=cpglm(Y1~SEX+ENGINE+KIND_OF_PAYMENT, link="log", data=panel1)
model.Y1
model.Y1@p   
model.Y1@phi

model.Y2=cpglm(Y2~SEX+ENGINE+KIND_OF_PAYMENT, link="log", data=panel2)
model.Y2
model.Y2@p   
model.Y2@phi

model.Y3=cpglm(Y3~SEX+ENGINE+KIND_OF_PAYMENT, link="log", data=panel3)
model.Y3
model.Y3@p   
model.Y3@phi

model.Y4=cpglm(Y4~SEX+ENGINE+KIND_OF_PAYMENT, link="log", data=panel4)
model.Y4
model.Y4@p   
model.Y4@phi

##Fitted values of total claim amounts
mu.1<-fitted.values(model.Y1)
mu.2<-fitted.values(model.Y2)
mu.3<-fitted.values(model.Y3)
mu.4<-fitted.values(model.Y4)

##Fig3
par(mfrow=c(2,2))
boxplot(mu.1~panel1$SEX+panel1$ENGINE+panel1$KIND_OF_PAYMENT, main="2007")
boxplot(mu.2~panel2$SEX+panel2$ENGINE+panel2$KIND_OF_PAYMENT, main="2008")
boxplot(mu.3~panel3$SEX+panel3$ENGINE+panel3$KIND_OF_PAYMENT, main="2009")
boxplot(mu.4~panel4$SEX+panel4$ENGINE+panel4$KIND_OF_PAYMENT, main="2010")

level <-interaction(panel1$SEX,panel1$ENGINE,panel1$KIND_OF_PAYMENT)
df1 <-data.frame(level, mu.1)
df2 <-data.frame(level, mu.2)
df3 <-data.frame(level, mu.3)
df4 <-data.frame(level, mu.4)

p1 <- ggplot(df1, aes(level, mu.1))+geom_boxplot()+theme_bw()+theme(panel.grid.major=element_blank()) 
p2 <- ggplot(df2, aes(level, mu.2))+geom_boxplot()+theme_bw()+theme(panel.grid.major=element_blank()) 
p3 <- ggplot(df3, aes(level, mu.3))+geom_boxplot()+theme_bw()+theme(panel.grid.major=element_blank()) 
p4 <- ggplot(df4, aes(level, mu.4))+geom_boxplot()+theme_bw()+theme(panel.grid.major=element_blank()) 

grid.arrange(p1, p2, p3, p4,
ncol=2, nrow=2)

##Estimation of copula parameter rho
U.Y1 <- ptweedie(Y1, power=model.Y1@p, mu=mu.1, phi=model.Y1@phi)
U.Y2 <- ptweedie(Y2, power=model.Y2@p, mu=mu.2, phi=model.Y2@phi)
U.Y3 <- ptweedie(Y3, power=model.Y3@p, mu=mu.3, phi=model.Y3@phi)
U.Y4 <- ptweedie(Y4, power=model.Y4@p, mu=mu.4, phi=model.Y4@phi)

U<-cbind(U.Y1,U.Y2,U.Y3,U.Y4)

fit.normal.ar1<- fitCopula(normalCopula(c(0),dim = 4, dispstr = "ar1"),U, start=0, method="ml")
fit.normal.ar1
rho.ar1<-coef(fit.normal.ar1)
summary(fit.normal.ar1)
rho.ar1

##Prediction (AR1)
zeta.t1<-qnorm(U.Y1)
zeta.t2<-qnorm(U.Y2)
zeta.t3<-qnorm(U.Y3)
zeta.t4<-qnorm(U.Y4)
zeta<-cbind(zeta.t1,zeta.t2,zeta.t3,zeta.t4)
zeta

sigma<-matrix(c(
1,rho.ar1,rho.ar1^2,rho.ar1^3,
rho.ar1,1,rho.ar1,rho.ar1^2,
rho.ar1^2,rho.ar1,1,rho.ar1,
rho.ar1^3,rho.ar1^2,rho.ar1,1), nrow=4,ncol=4)
sigma

sigma.4na5<-as.vector(c(rho.ar1^4,rho.ar1^3,rho.ar1^2,rho.ar1))
sigma.4na5

mi.z.5<-rep(1,245)
for (i in 1:245){mi.z.5[i]<-sigma.4na5%*%solve(sigma)%*%zeta[i,]}
mi.z.5

sigma.z.5<-1-sigma.4na5%*%solve(sigma)%*%sigma.4na5
sigma.z.5

############
###Simulation (takes time)
l.sym <- 500 
pred <- rep(1,245)
for (i in 1:245){
	sym.z <- rnorm(l.sym)
	fi <- pnorm(mi.z.5[i]+sigma.z.5*sym.z)
	F_odw <- qtweedie(fi, power=model.Y4@p, mu=fitted(model.Y4)[i], phi=model.Y4@phi)
	pred[i] <- mean(F_odw)
		}

#### Results data frame
results <- data.frame(SEX=SEX,ENGINE=ENGINE,KIND_OF_PAYMENT=KIND_OF_PAYMENT,Y.1=Y1,Y.2=Y2, Y.3=Y3, Y.4=Y4,Risk.Group=group, Pre.Tot.Claim=pred)

##Fig4
results.ogr1 <-results[results$Pre.Tot.Claim>500,]
results.ogr <-results[results$Pre.Tot.Claim<500,]

par(mfrow=c(1,2))
boxplot(Pre.Tot.Claim~Risk.Group, data=results,main="a) Predicted total claim amount")  
boxplot(Pre.Tot.Claim~Risk.Group, data=results.ogr,main="b) Predicted total claim amount <500")

############
##Fig5
mean.Y1<-rbind(mean(results[results$Risk.Group=="F.BEN.C", ]$Y.1),
mean(results[results$Risk.Group=="M.BEN.C", ]$Y.1),
mean(results[results$Risk.Group=="F.DIE.C", ]$Y.1),
mean(results[results$Risk.Group=="M.DIE.C", ]$Y.1),
mean(results[results$Risk.Group=="F.BEN.T", ]$Y.1),
mean(results[results$Risk.Group=="M.BEN.T", ]$Y.1),
mean(results[results$Risk.Group=="F.DIE.T", ]$Y.1),
mean(results[results$Risk.Group=="M.DIE.T", ]$Y.1))

mean.Y2<-rbind(mean(results[results$Risk.Group=="F.BEN.C", ]$Y.1),
mean(results[results$Risk.Group=="M.BEN.C", ]$Y.2),
mean(results[results$Risk.Group=="F.DIE.C", ]$Y.2),
mean(results[results$Risk.Group=="M.DIE.C", ]$Y.2),
mean(results[results$Risk.Group=="F.BEN.T", ]$Y.2),
mean(results[results$Risk.Group=="M.BEN.T", ]$Y.2),
mean(results[results$Risk.Group=="F.DIE.T", ]$Y.2),
mean(results[results$Risk.Group=="M.DIE.T", ]$Y.2))

mean.Y3<-rbind(mean(results[results$Risk.Group=="F.BEN.C", ]$Y.1),
mean(results[results$Risk.Group=="M.BEN.C", ]$Y.3),
mean(results[results$Risk.Group=="F.DIE.C", ]$Y.3),
mean(results[results$Risk.Group=="M.DIE.C", ]$Y.3),
mean(results[results$Risk.Group=="F.BEN.T", ]$Y.3),
mean(results[results$Risk.Group=="M.BEN.T", ]$Y.3),
mean(results[results$Risk.Group=="F.DIE.T", ]$Y.3),
mean(results[results$Risk.Group=="M.DIE.T", ]$Y.3))


mean.Y4<-rbind(mean(results[results$Risk.Group=="F.BEN.C", ]$Y.1),
mean(results[results$Risk.Group=="M.BEN.C", ]$Y.4),
mean(results[results$Risk.Group=="F.DIE.C", ]$Y.4),
mean(results[results$Risk.Group=="M.DIE.C", ]$Y.4),
mean(results[results$Risk.Group=="F.BEN.T", ]$Y.4),
mean(results[results$Risk.Group=="M.BEN.T", ]$Y.4),
mean(results[results$Risk.Group=="F.DIE.T", ]$Y.4),
mean(results[results$Risk.Group=="M.DIE.T", ]$Y.4))


tot.claim.g1<-mean(results[results$Risk.Group=="F.BEN.C", ]$Pre.Tot.Claim)
tot.claim.g2<-mean(results[results$Risk.Group=="M.BEN.C", ]$Pre.Tot.Claim)
tot.claim.g3<-mean(results[results$Risk.Group=="F.DIE.C", ]$Pre.Tot.Claim)
tot.claim.g4<-mean(results[results$Risk.Group=="M.DIE.C", ]$Pre.Tot.Claim)
tot.claim.g5<-mean(results[results$Risk.Group=="F.BEN.T", ]$Pre.Tot.Claim)
tot.claim.g6<-mean(results[results$Risk.Group=="M.BEN.T", ]$Pre.Tot.Claim)
tot.claim.g7<-mean(results[results$Risk.Group=="F.DIE.T", ]$Pre.Tot.Claim)
tot.claim.g8<-mean(results[results$Risk.Group=="M.DIE.T", ]$Pre.Tot.Claim)
tot.claim<-rbind(tot.claim.g1,tot.claim.g2,tot.claim.g3,tot.claim.g4, tot.claim.g5, tot.claim.g6, tot.claim.g7,tot.claim.g8)
rownames(tot.claim)<-legend$Group

par(mfrow=c(1,2))
plot(density(results.ogr$Pre.Tot.Claim[results.ogr$Risk.Group=="F.BEN.C"]), xlim=c(0,400), ylim=c(0,0.06), main="a) Female* ", xlab=" ")
lines(density(results.ogr$Pre.Tot.Claim[results.ogr$Risk.Group=="F.BEN.T"]), lty=2)
lines(density(results.ogr$Pre.Tot.Claim[results.ogr$Risk.Group=="F.DIE.C"]), lty=3)
legend("topright", c("F.BEN.C","F.BEN.T","F.DIE.C"), lty=1:3)


plot(density(results.ogr$Pre.Tot.Claim[results.ogr$Risk.Group=="M.BEN.C"]), xlim=c(0,400), ylim=c(0,0.06), main="b) Male ", xlab=" ")
lines(density(results.ogr$Pre.Tot.Claim[results.ogr$Risk.Group=="M.BEN.T"]), lty=2)
lines(density(results.ogr$Pre.Tot.Claim[results.ogr$Risk.Group=="M.DIE.C"]), lty=3)
lines(density(results.ogr$Pre.Tot.Claim[results.ogr$Risk.Group=="M.DIE.T"]), lty=4)
legend("topright", c("M.BEN.C","M.BEN.T","M.DIE.C","M.DIE.T"), lty=1:4)


##########
##Fig6

mean.claim.pred<-cbind(mean.Y1,mean.Y2,mean.Y3,mean.Y4, tot.claim)
colnames(mean.claim.pred)<-c("2007", "2008", "2009", "2010", "2011 (pred)")
mean.claim.pred

barplot(tot.claim[,1], axis.lty=1)
abline(h=0)
abline(h=mean(pred), lty=2)
mean<-expression(paste("Mean of predictors, t=2011"))
text(7, 700, pos=2, mean)
