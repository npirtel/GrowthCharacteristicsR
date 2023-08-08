#packages##########################################################################################
install.packages("ggplot2")
install.packages("gridExtra")
install.packages("ggpubr")
install.packages("dplyr")
install.packages("tidyr")
install.packages("ggforce")
install.packages("plotly")
install.packages("metR")
install.packages("DataCombine")
install.packages("vctrs")
install.packages("ggnewscale")
install.packages("rstatix")

library(ggplot2)
library(gridExtra)
library(ggpubr)
library(dplyr)
library(tidyr)
library(cowplot)
library(ggforce)
library(plotly)
library(grid)
library(metR)
library(ggnewscale)
library(rstatix)






##rooting depth vs diameter############################################################################# 

roots <- read.csv(file = "rootdiameter.csv")


#fir################################################################################
firroot <- filter(roots, species == 'Fir')
firroot$hclass <- factor(pinroot$hclass, levels = c("<150", "150-300", "300+"))
firroot150 <- subset(ponroot, hclass == "<150")
firroot150300 <- subset(ponroot, hclass == "150-300")
firroot300 <- subset(ponroot, hclass == "300+")

wide_firroot <- firroot %>% spread(sample, diameter)
write.csv(wide_firroot, file = "fir/fir rooting depth.csv")
wide_firroot1 <- read.csv(file = "fir/fir rooting depth1.csv")


firint <- 0
firslope <- 0
firp <- 0
firrsq <- 0
firrmse <- 0

for (i in (2:length(wide_firroot1[1,]))) {
  z <- wide_firroot1[!is.na(wide_firroot1[,i]), c(1,i)]
  fit <- lm(z$depth~(z[, 2]))
  firint[c(1,i)] <- coef(fit)["(Intercept)"]
  firslope[c(1,i)] <- coef(fit)["z[, 2]"]
  firp[c(1,i)] <- summary(fit)$coefficients[,4]
  firrsq[c(1,i)] <- summary(fit)$adj.r.squared
  firrmse[c(1,i)] <- sqrt(mean(fit$residuals^2))
}

dffir1 <- data.frame(rbind(firint, firslope, firp, firrsq, firrmse))


firint2 <- 0
firslope2 <- 0
firp2 <- 0
firrsq2 <- 0
firrmse2 <- 0

for (ii in (2:length(wide_firroot1[1,]))) {
  z <- wide_firroot1[!is.na(wide_firroot1[,ii]), c(1,ii)]
  fit2 <- lm(z$depth~log(z[, 2]))
             firint2[c(1,ii)] <- coef(fit2)["(Intercept)"]
             firslope2[c(1,ii)] <- coef(fit2)["log(z[, 2])"]
             firp2[c(1,ii)] <- summary(fit2)$coefficients[,4]
             firrsq2[c(1,ii)] <- summary(fit2)$adj.r.squared
             firrmse2[c(1,ii)] <- sqrt(mean(fit2$residuals^2))
}
dffir2 <- data.frame(rbind(firint2, firslope2, firp2, firrsq2, firrmse2))


firint3 <- 0
firslope3 <- 0
firp3 <- 0
firrsq3 <- 0
firslopea3 <- 0
firslopeb3 <- 0
firrmse3 <- 0

#ax^2+bx+c
for (iii in (2:length(wide_firroot1[1,]))) {
  z <- wide_firroot1[!is.na(wide_firroot1[,iii]), c(1,iii)]
  fit3 <- lm(z$depth~poly((z[, 2]), 2, raw=TRUE))
  firint3[c(1,iii)] <- coef(fit3)["(Intercept)"]
  firslopea3[c(1,iii)] <- coef(fit3)["poly((z[, 2]), 2, raw = TRUE)1"]
  firslopeb3[c(1,iii)] <- coef(fit3)["poly((z[, 2]), 2, raw = TRUE)2"]
  firp3[c(1,iii)] <- summary(fit3)$coefficients[,4]
  firrsq3[c(1,iii)] <- summary(fit3)$adj.r.squared
  firrmse3[c(1,iii)] <- sqrt(mean(fit3$residuals^2))
}
dffir3 <- data.frame(rbind(firint3, firslopea3, firslopeb3, firp3, firrsq3, firrmse3))


firint4 <- 0
firslope4 <- 0
firp4 <- 0
firrsq4 <- 0
firslopea4 <- 0
firslopeb4 <- 0
firslopec4 <- 0
firrmse4 <- 0

#ax^3+bx^2+cx+d
for (iv in (2:length(wide_firroot1[1,]))) {
  z <- wide_firroot1[!is.na(wide_firroot1[,iv]), c(1,iv)]
  fit4 <- lm(z$depth~poly((z[, 2]), 3, raw=TRUE))
  firint4[c(1,iv)] <- coef(fit4)["(Intercept)"]
  firslopea4[c(1,iv)] <- coef(fit4)["poly((z[, 2]), 3, raw = TRUE)1"]
  firslopeb4[c(1,iv)] <- coef(fit4)["poly((z[, 2]), 3, raw = TRUE)2"]
  firslopec4[c(1,iv)] <- coef(fit4)["poly((z[, 2]), 3, raw = TRUE)3"]
  firp4[c(1,iv)] <- summary(fit4)$coefficients[,4]
  firrsq4[c(1,iv)] <- summary(fit4)$adj.r.squared
  firrmse4[c(1,iv)] <- sqrt(mean(fit4$residuals^2))
}
dffir4 <- data.frame(rbind(firint4, firslopea4, firslopeb4, firslopec4, firp4, firrsq4,firrmse4))

allfircoef <- data.frame(rbind(dffir1,dffir2,dffir3,dffir4))
fircompare <- data.frame(rbind(firp,firp2,firp3,firp4,firrsq,firrsq2,firrsq3,firrsq4))
#compare sig p-values of different lines excluding first column/intercept
#keep linear and log
rowSums(fircompare[-1] < 0.05)

#compare avg fit of lines to distribution points
#log has higher R2 than linear
mean(firrsq)
mean(firrsq2)
mean(firrsq3)
mean(firrsq4)

write.csv(fircompare, file = "fir/fir prsq comparison.csv")
write.csv(allfircoef, file = "fir/fir reg coefficients.csv")

allfircoef1 <- read.csv(file = "fir/fir max root depths.csv")


fircoef150 <- subset(allfircoef1, hclass == "<150")
fircoef150300 <- subset(allfircoef1, hclass == "150-300")




#juniper####################################################################################
junroot <- filter(roots, species == 'Juniper')
wide_junroot <- junroot %>% spread(sample, diameter)
write.csv(wide_junroot, file = "juniper/juniper rooting depth.csv")
wide_junroot1 <- read.csv(file = "juniper/juniper rooting depth1.csv")


junint <- 0
junslope <- 0
junp <- 0
junrsq <- 0
junrmse <- 0

for (i in (2:length(wide_junroot1[1,]))) {
  x <- wide_junroot1[!is.na(wide_junroot1[,i]), c(1,i)]
  fit <- lm(x$depth~(x[, 2]))
  junint[c(1,i)] <- coef(fit)["(Intercept)"]
  junslope[c(1,i)] <- coef(fit)["x[, 2]"]
  junp[c(1,i)] <- summary(fit)$coefficients[,4]
  junrsq[c(1,i)] <- summary(fit)$adj.r.squared
  junrmse[c(1,i)] <- sqrt(mean(fit$residuals^2))
}
dfjun1 <- data.frame(rbind(junint, junslope, junp, junrsq, junrmse))


junint2 <- 0
junslope2 <- 0
junp2 <- 0
junrsq2 <- 0
junrmse2 <- 0

for (ii in (2:length(wide_junroot1[1,]))) {
  x <- wide_junroot1[!is.na(wide_junroot1[,ii]), c(1,ii)]
  fit2 <- lm(x$depth~log(x[, 2]))
  junint2[c(1,ii)] <- coef(fit2)["(Intercept)"]
  junslope2[c(1,ii)] <- coef(fit2)["log(x[, 2])"]
  junp2[c(1,ii)] <- summary(fit2)$coefficients[,4]
  junrsq2[c(1,ii)] <- summary(fit2)$adj.r.squared
  junrmse2[c(1,ii)] <- sqrt(mean(fit2$residuals^2))
}
dfjun2 <- data.frame(rbind(junint2, junslope2, junp2, junrsq2, junrmse2))


junint3 <- 0
junslopea3 <- 0
junslopeb3 <- 0
junp3 <- 0
junrsq3 <- 0
junrmse3 <- 0

#ax^2+bx+c
for (iii in (2:length(wide_junroot1[1,]))) {
  x <- wide_junroot1[!is.na(wide_junroot1[,iii]), c(1,iii)]
  fit3 <- lm(x$depth~poly((x[, 2]), 2, raw=TRUE))
  junint3[c(1,iii)] <- coef(fit3)["(Intercept)"]
  junslopea3[c(1,iii)] <- coef(fit3)["poly((x[, 2]), 2, raw = TRUE)1"]
  junslopeb3[c(1,iii)] <- coef(fit3)["poly((x[, 2]), 2, raw = TRUE)2"]
  junp3[c(1,iii)] <- summary(fit3)$coefficients[,4]
  junrsq3[c(1,iii)] <- summary(fit3)$adj.r.squared
  junrmse3[c(1,iii)] <- sqrt(mean(fit3$residuals^2))
}
dfjun3 <- data.frame(rbind(junint3, junslopea3, junslopeb3, junp3, junrsq3, junrmse3))


junint4 <- 0
junslopea4 <- 0
junslopeb4 <- 0
junslopec4 <- 0
junp4 <- 0
junrsq4 <- 0
junrmse4 <- 0

#ax^3+bx^2+cx+d
for (iv in (2:length(wide_junroot1[1,]))) {
  x <- wide_junroot1[!is.na(wide_junroot1[,iv]), c(1,iv)]
  fit4 <- lm(x$depth~poly((x[, 2]), 3, raw=TRUE))
  junint4[c(1,iv)] <- coef(fit4)["(Intercept)"]
  junslopea4[c(1,iv)] <- coef(fit4)["poly((x[, 2]), 3, raw = TRUE)1"]
  junslopeb4[c(1,iv)] <- coef(fit4)["poly((x[, 2]), 3, raw = TRUE)2"]
  junslopec4[c(1,iv)] <- coef(fit4)["poly((x[, 2]), 3, raw = TRUE)3"]
  junp4[c(1,iv)] <- summary(fit4)$coefficients[,4]
  junrsq4[c(1,iv)] <- summary(fit4)$adj.r.squared
  junrmse4[c(1,iv)] <- sqrt(mean(fit4$residuals^2))
}
dfjun4 <- data.frame(rbind(junint4, junslopea4, junslopeb4, junslopec4, junp4, junrsq4, junrmse4))

#for regression eq coefficients
summary(fit)
summary(fit2)
summary(fit3)
summary(fit4)


alljuncoef <- data.frame(rbind(dfjun1,dfjun2,dfjun3,dfjun4))
juncompare <- data.frame(rbind(junp,junp2,junp3,junp4,junrsq,junrsq2,junrsq3,junrsq4))
#compare sig p-values of different lines excluding first column/intercept
#keep linear and log
rowSums(juncompare[-1] < 0.05)

#compare avg fit of lines to distribution points
#log has higher average R2 than linear
mean(junrsq)
mean(junrsq2)
mean(junrsq3)
mean(junrsq4)

write.csv(juncompare, file = "juniper/jun prsq comparison.csv")
write.csv(alljuncoef, file = "juniper/jun reg coefficients.csv")


alljuncoef1 <- read.csv(file = "juniper/jun max root depths.csv")


juncoef150 <- subset(alljuncoef1, hclass == "<150")
juncoef150300 <- subset(alljuncoef1, hclass =! "<150")
juncoef300 <- subset(alljuncoef1, hclass == "300+")
junroot150 <- subset(junroot, hclass == "<150")
junroot150300 <- subset(junroot, hclass == "150-300")
junroot300 <- subset(junroot, hclass == "300+")




#pinon#######################################################################################
pinroot <- filter(roots, species == 'Pinon')
wide_pinroot <- pinroot %>% spread(sample, diameter)
write.csv(wide_pinroot, file = "pinon/pinon rooting depth.csv")
wide_pinroot1 <- read.csv(file = "pinon/pinon rooting depth1.csv")


pinint <- 0
pinslope <- 0
pinp <- 0
pinrsq <- 0
pinrmse <- 0

for (i in (2:length(wide_pinroot1[1,]))) {
  y <- wide_pinroot1[!is.na(wide_pinroot1[,i]), c(1,i)]
  fit <- lm(y$depth~(y[, 2]))
  pinint[c(1,i)] <- coef(fit)["(Intercept)"]
  pinslope[c(1,i)] <- coef(fit)["y[, 2]"]
  pinp[c(1,i)] <- summary(fit)$coefficients[,4]
  pinrsq[c(1,i)] <- summary(fit)$adj.r.squared
  pinrmse[c(1,i)] <- sqrt(mean(fit$residuals^2))
}
dfpin1 <- data.frame(rbind(pinint, pinslope, pinp, pinrsq, pinrmse))


pinint2 <- 0
pinslope2 <- 0
pinp2 <- 0
pinrsq2 <- 0
pinrmse2 <- 0

for (ii in (2:length(wide_pinroot1[1,]))) {
  y <- wide_pinroot1[!is.na(wide_pinroot1[,ii]), c(1,ii)]
  fit2 <- lm(y$depth~log(y[, 2]))
  pinint2[c(1,ii)] <- coef(fit2)["(Intercept)"]
  pinslope2[c(1,ii)] <- coef(fit2)["log(y[, 2])"]
  pinp2[c(1,ii)] <- summary(fit2)$coefficients[,4]
  pinrsq2[c(1,ii)] <- summary(fit2)$adj.r.squared
  pinrmse2[c(1,ii)] <- sqrt(mean(fit2$residuals^2))
}
dfpin2 <- data.frame(rbind(pinint2, pinslope2, pinp2, pinrsq2, pinrmse2))


pinint3 <- 0
pinslopea3 <- 0
pinslopeb3 <- 0
pinp3 <- 0
pinrsq3 <- 0
pinrmse3 <- 0

#ax^2+bx+c
for (iii in (2:length(wide_pinroot1[1,]))) {
  y <- wide_pinroot1[!is.na(wide_pinroot1[,iii]), c(1,iii)]
  fit3 <- lm(y$depth~poly((y[, 2]), 2, raw=TRUE))
  pinint3[c(1,iii)] <- coef(fit3)["(Intercept)"]
  pinslopea3[c(1,iii)] <- coef(fit3)["poly((y[, 2]), 2, raw = TRUE)1"]
  pinslopeb3[c(1,iii)] <- coef(fit3)["poly((y[, 2]), 2, raw = TRUE)2"]
  pinp3[c(1,iii)] <- summary(fit3)$coefficients[,4]
  pinrsq3[c(1,iii)] <- summary(fit3)$adj.r.squared
  pinrmse3[c(1,iii)] <- sqrt(mean(fit3$residuals^2))
}
dfpin3 <- data.frame(rbind(pinint3, pinslopea3, pinslopeb3, pinp3, pinrsq3, pinrmse3))


pinint4 <- 0
pinslopea4 <- 0
pinslopeb4 <- 0
pinslopec4 <- 0
pinp4 <- 0
pinrsq4 <- 0
pinrmse4 <- 0

#ax^3+bx^2+cx+d
for (iv in (2:length(wide_pinroot1[1,]))) {
  y <- wide_pinroot1[!is.na(wide_pinroot1[,iv]), c(1,iv)]
  fit4 <- lm(y$depth~poly((y[, 2]), 3, raw=TRUE))
  pinint4[c(1,iv)] <- coef(fit4)["(Intercept)"]
  pinslopea4[c(1,iv)] <- coef(fit4)["poly((y[, 2]), 3, raw = TRUE)1"]
  pinslopeb4[c(1,iv)] <- coef(fit4)["poly((y[, 2]), 3, raw = TRUE)2"]
  pinslopec4[c(1,iv)] <- coef(fit4)["poly((y[, 2]), 3, raw = TRUE)3"]
  pinp4[c(1,iv)] <- summary(fit4)$coefficients[,4]
  pinrsq4[c(1,iv)] <- summary(fit4)$adj.r.squared
  pinrmse4[c(1,iv)] <- sqrt(mean(fit4$residuals^2))
}
dfpin4 <- data.frame(rbind(pinint4, pinslopea4, pinslopeb4, pinslopec4, pinp4, pinrsq4, pinrmse4))

#for regression eq coefficients
summary(fit)
summary(fit2)
summary(fit3)
summary(fit4)


allpincoef <- data.frame(rbind(dfpin1,dfpin2,dfpin3,dfpin4))
pincompare <- data.frame(rbind(pinp,pinp2,pinp3,pinp4,pinrsq,pinrsq2,pinrsq3,pinrsq4))
#compare sig p-values of different lines excluding first column/intercept
#keep linear and log
rowSums(pincompare[-1] < 0.05)

#compare avg fit of lines to distribution points
#log has higher average R2 than linear
mean(pinrsq)
mean(pinrsq2)
mean(pinrsq3)
mean(pinrsq4)

write.csv(pincompare, file = "pinon/pin prsq comparison.csv")
write.csv(allpincoef, file = "pinon/pin reg coefficients.csv")


allpincoef1 <- read.csv(file = "pinon/pin max root depths.csv")


pincoef150 <- subset(allpincoef1, hclass == "<150")
pincoef150300 <- subset(allpincoef1, hclass =! "<150")
pincoef300 <- subset(allpincoef1, hclass == "300+") 
pinroot150 <- subset(pinroot, hclass == "<150")
pinroot150300 <- subset(pinroot, hclass == "150-300")
pinroot300 <- subset(pinroot, hclass == "300+")




#ponderosa#####################################################################################
ponroot <- filter(roots, species == 'Ponderosa')
ponroot$hclass <- factor(ponroot$hclass, levels = c("<150", "150-300", "300+"))
ponroot150 <- subset(ponroot, hclass == "<150")
ponroot150300 <- subset(ponroot, hclass == "150-300")
ponroot300 <- subset(ponroot, hclass == "300+")

wide_ponroot <- ponroot %>% spread(sample, diameter)
write.csv(wide_ponroot, file = "ponderosa/ponderosa rooting depth.csv")
wide_ponroot1 <- read.csv(file = "ponderosa/ponderosa rooting depth1.csv")


ponint <- 0
ponslope <- 0
ponp <- 0
ponrsq <- 0
ponrmse <- 0

for (i in (2:length(wide_ponroot1[1,]))) {
  w <- wide_ponroot1[!is.na(wide_ponroot1[,i]), c(1,i)]
  fit <- lm(w$depth~(w[, 2]))
  ponint[c(1,i)] <- coef(fit)["(Intercept)"]
  ponslope[c(1,i)] <- coef(fit)["w[, 2]"]
  ponp[c(1,i)] <- summary(fit)$coefficients[,4]
  ponrsq[c(1,i)] <- summary(fit)$adj.r.squared
  ponrmse[c(1,i)] <- sqrt(mean(fit$residuals^2))
}
dfpon1 <- data.frame(rbind(ponint, ponslope, ponp, ponrsq, ponrmse))


ponint2 <- 0
ponslope2 <- 0
ponp2 <- 0
ponrsq2 <- 0
ponrmse2 <- 0

for (ii in (2:length(wide_ponroot1[1,]))) {
  w <- wide_ponroot1[!is.na(wide_ponroot1[,ii]), c(1,ii)]
  fit2 <- lm(w$depth~log(w[, 2]))
  ponint2[c(1,ii)] <- coef(fit2)["(Intercept)"]
  ponslope2[c(1,ii)] <- coef(fit2)["log(w[, 2])"]
  ponp2[c(1,ii)] <- summary(fit2)$coefficients[,4]
  ponrsq2[c(1,ii)] <- summary(fit2)$adj.r.squared
  ponrmse2[c(1,ii)] <- sqrt(mean(fit2$residuals^2))
}
dfpon2 <- data.frame(rbind(ponint2, ponslope2, ponp2, ponrsq2, ponrmse2))


ponint3 <- 0
ponslopea3 <- 0
ponslopeb3 <- 0
ponp3 <- 0
ponrsq3 <- 0
ponrmse3 <- 0

#ax^2+bx+c
for (iii in (2:length(wide_ponroot1[1,]))) {
  w <- wide_ponroot1[!is.na(wide_ponroot1[,iii]), c(1,iii)]
  fit3 <- lm(w$depth~poly((w[, 2]), 2, raw=TRUE))
  ponint3[c(1,iii)] <- coef(fit3)["(Intercept)"]
  ponslopea3[c(1,iii)] <- coef(fit3)["poly((w[, 2]), 2, raw = TRUE)1"]
  ponslopeb3[c(1,iii)] <- coef(fit3)["poly((w[, 2]), 2, raw = TRUE)2"]
  ponp3[c(1,iii)] <- summary(fit3)$coefficients[,4]
  ponrsq3[c(1,iii)] <- summary(fit3)$adj.r.squared
  ponrmse3[c(1,iii)] <- sqrt(mean(fit3$residuals^2))
}
dfpon3 <- data.frame(rbind(ponint3, ponslopea3, ponslopeb3, ponp3, ponrsq3, ponrmse3))


ponint4 <- 0
ponslopea4 <- 0
ponslopeb4 <- 0
ponslopec4 <- 0
ponp4 <- 0
ponrsq4 <- 0
ponrmse4 <- 0

#ax^3+bx^2+cx+d
for (iv in (2:length(wide_ponroot1[1,]))) {
  w <- wide_ponroot1[!is.na(wide_ponroot1[,iv]), c(1,iv)]
  fit4 <- lm(w$depth~poly((w[, 2]), 3, raw=TRUE))
  ponint4[c(1,iv)] <- coef(fit4)["(Intercept)"]
  ponslopea4[c(1,iv)] <- coef(fit4)["poly((w[, 2]), 3, raw = TRUE)1"]
  ponslopeb4[c(1,iv)] <- coef(fit4)["poly((w[, 2]), 3, raw = TRUE)2"]
  ponslopec4[c(1,iv)] <- coef(fit4)["poly((w[, 2]), 3, raw = TRUE)3"]
  ponp4[c(1,iv)] <- summary(fit4)$coefficients[,4]
  ponrsq4[c(1,iv)] <- summary(fit4)$adj.r.squared
  ponrmse4[c(1,iv)] <- sqrt(mean(fit4$residuals^2))
}
dfpon4 <- data.frame(rbind(ponint4, ponslopea4, ponslopeb4, ponslopec4, ponp4, ponrsq4, ponrmse4))

#for regression eq coefficients
summary(fit)
summary(fit2)
summary(fit3)
summary(fit4)


allponcoef <- data.frame(rbind(dfpon1,dfpon2,dfpon3,dfpon4))
poncompare <- data.frame(rbind(ponp,ponp2,ponp3,ponp4,ponrsq,ponrsq2,ponrsq3,ponrsq4))
#compare sig p-values of different lines excluding first column/intercept
#keep linear, log, second degree polynomial
rowSums(poncompare[-1] < 0.05)

#compare avg fit of lines to distribution points
#2nd deg poly highest average R2, this and log much higher than linear
mean(ponrsq)
mean(ponrsq2)
mean(ponrsq3)
mean(ponrsq4)

write.csv(poncompare, file = "ponderosa/pon prsq comparison.csv")
write.csv(allponcoef, file = "ponderosa/pon reg coefficients.csv")


allponcoef1 <- read.csv(file = "ponderosa/pon max root depths.csv")


poncoef150 <- subset(allponcoef1, hclass == "<150")
poncoef150300 <- subset(allponcoef1, hclass == "150-300")
poncoef300 <- subset(allponcoef1, hclass == "300+")








#max rooting depth individual calculations##############################################################
diameter <- 0.3
diameter5 <- 0.5
diameter75 <- 0.75
diameter1 <- 1

#fir
rootdep1 <- (firslope*diameter)+firint
rootdep15 <- (firslope*(diameter5))+firint
rootdep175 <- (firslope*(diameter75))+firint
rootdep11 <- (firslope*(diameter1))+firint

rootdep2 <- (firslope2*log(diameter))+firint2
rootdep25 <- (firslope2*log(diameter5))+firint2
rootdep275 <- (firslope2*log(diameter75))+firint2
rootdep21 <- (firslope2*log(diameter1))+firint2

rootdep3 <- (firslopeb3*(diameter^2))+(firslopea3*diameter)+firint3
rootdep35 <- (firslopeb3*(diameter5^2))+(firslopea3*diameter5)+firint3
rootdep375 <- (firslopeb3*(diameter75^2))+(firslopea3*diameter75)+firint3
rootdep31 <- (firslopeb3*(diameter1^2))+(firslopea3*diameter1)+firint3

rootdep4 <- (firslopec4*(diameter^3))+(firslopeb4*(diameter^2))+(firslopea4*diameter)+firint4
rootdep45 <- (firslopec4*(diameter5^3))+(firslopeb4*(diameter5^2))+(firslopea4*diameter5)+firint4
rootdep475 <- (firslopec4*(diameter75^3))+(firslopeb4*(diameter75^2))+(firslopea4*diameter75)+firint4
rootdep41 <- (firslopec4*(diameter1^3))+(firslopeb4*(diameter1^2))+(firslopea4*diameter1)+firint4


firrootdeps <- data.frame(rbind(rootdep1,rootdep2,rootdep3,rootdep4))
firrootdeps5 <- data.frame(rbind(rootdep15,rootdep25,rootdep35,rootdep45))
firrootdeps75 <- data.frame(rbind(rootdep175,rootdep275,rootdep375,rootdep475))
firrootdeps1 <- data.frame(cbind(rootdep11,rootdep21,rootdep31,rootdep41))

write.csv(firrootdeps, file = "fir/fir max root depths.csv")


#jun
rootdep1 <- (junslope*diameter)+junint
rootdep15 <- (junslope*(diameter5))+junint
rootdep175 <- (junslope*(diameter75))+junint
rootdep11 <- (junslope*(diameter1))+junint

rootdep2 <- (junslope2*log(diameter))+junint2
rootdep25 <- (junslope2*log(diameter5))+junint2
rootdep275 <- (junslope2*log(diameter75))+junint2
rootdep21 <- (junslope2*log(diameter1))+junint2

rootdep3 <- (junslopeb3*(diameter^2))+(junslopea3*diameter)+junint3
rootdep35 <- (junslopeb3*(diameter5^2))+(junslopea3*diameter5)+junint3
rootdep375 <- (junslopeb3*(diameter75^2))+(junslopea3*diameter75)+junint3
rootdep31 <- (junslopeb3*(diameter1^2))+(junslopea3*diameter1)+junint3

rootdep4 <- (junslopec4*(diameter^3))+(junslopeb4*(diameter^2))+(junslopea4*diameter)+junint4
rootdep45 <- (junslopec4*(diameter5^3))+(junslopeb4*(diameter5^2))+(junslopea4*diameter5)+junint4
rootdep475 <- (junslopec4*(diameter75^3))+(junslopeb4*(diameter75^2))+(junslopea4*diameter75)+junint4
rootdep41 <- (junslopec4*(diameter1^3))+(junslopeb4*(diameter1^2))+(junslopea4*diameter1)+junint4

junrootdeps <- data.frame(cbind(rootdep1,rootdep2,rootdep3,rootdep4))
junrootdeps5 <- data.frame(cbind(rootdep15,rootdep25,rootdep35,rootdep45))
junrootdeps75 <- data.frame(cbind(rootdep175,rootdep275,rootdep375,rootdep475))
junrootdeps1 <- data.frame(cbind(rootdep11,rootdep21,rootdep31,rootdep41))

write.csv(junrootdeps, file = "juniper/jun max root depths.csv")


#pin
rootdep1 <- (pinslope*diameter)+pinint
rootdep15 <- (pinslope*(diameter5))+pinint
rootdep175 <- (pinslope*(diameter75))+pinint
rootdep11 <- (pinslope*(diameter1))+pinint

rootdep2 <- (pinslope2*log(diameter))+pinint2
rootdep25 <- (pinslope2*log(diameter5))+pinint2
rootdep275 <- (pinslope2*log(diameter75))+pinint2
rootdep21 <- (pinslope2*log(diameter1))+pinint2

rootdep3 <- (pinslopeb3*(diameter^2))+(pinslopea3*diameter)+pinint3
rootdep35 <- (pinslopeb3*(diameter5^2))+(pinslopea3*diameter5)+pinint3
rootdep375 <- (pinslopeb3*(diameter75^2))+(pinslopea3*diameter75)+pinint3
rootdep31 <- (pinslopeb3*(diameter1^2))+(pinslopea3*diameter1)+pinint3

rootdep4 <- (pinslopec4*(diameter^3))+(pinslopeb4*(diameter^2))+(pinslopea4*diameter)+pinint4
rootdep45 <- (pinslopec4*(diameter5^3))+(pinslopeb4*(diameter5^2))+(pinslopea4*diameter5)+pinint4
rootdep475 <- (pinslopec4*(diameter75^3))+(pinslopeb4*(diameter75^2))+(pinslopea4*diameter75)+pinint4
rootdep41 <- (pinslopec4*(diameter1^3))+(pinslopeb4*(diameter1^2))+(pinslopea4*diameter1)+pinint4

pinrootdeps <- data.frame(cbind(rootdep1,rootdep2,rootdep3,rootdep4))
pinrootdeps5 <- data.frame(cbind(rootdep15,rootdep25,rootdep35,rootdep45))
pinrootdeps75 <- data.frame(cbind(rootdep175,rootdep275,rootdep375,rootdep475))
pinrootdeps1 <- data.frame(cbind(rootdep11,rootdep21,rootdep31,rootdep41))

write.csv(pinrootdeps, file = "pinon/pin max root depths.csv")


#pon
rootdep1 <- (ponslope*diameter)+ponint
rootdep15 <- (ponslope*(diameter5))+ponint
rootdep175 <- (ponslope*(diameter75))+ponint
rootdep11 <- (ponslope*(diameter1))+ponint

rootdep2 <- (ponslope2*log(diameter))+ponint2
rootdep25 <- (ponslope2*log(diameter5))+ponint2
rootdep275 <- (ponslope2*log(diameter75))+ponint2
rootdep21 <- (ponslope2*log(diameter1))+ponint2

rootdep3 <- (ponslopeb3*(diameter^2))+(ponslopea3*diameter)+ponint3
rootdep35 <- (ponslopeb3*(diameter5^2))+(ponslopea3*diameter5)+ponint3
rootdep375 <- (ponslopeb3*(diameter75^2))+(ponslopea3*diameter75)+ponint3
rootdep31 <- (ponslopeb3*(diameter1^2))+(ponslopea3*diameter1)+ponint3

rootdep4 <- (ponslopec4*(diameter^3))+(ponslopeb4*(diameter^2))+(ponslopea4*diameter)+ponint4
rootdep45 <- (ponslopec4*(diameter5^3))+(ponslopeb4*(diameter5^2))+(ponslopea4*diameter5)+ponint4
rootdep475 <- (ponslopec4*(diameter75^3))+(ponslopeb4*(diameter75^2))+(ponslopea4*diameter75)+ponint4
rootdep41 <- (ponslopec4*(diameter1^3))+(ponslopeb4*(diameter1^2))+(ponslopea4*diameter1)+ponint4

ponrootdeps <- data.frame(cbind(rootdep1,rootdep2,rootdep3,rootdep4))
ponrootdeps5 <- data.frame(cbind(rootdep15,rootdep25,rootdep35,rootdep45))
ponrootdeps75 <- data.frame(cbind(rootdep175,rootdep275,rootdep375,rootdep475))
ponrootdeps1 <- data.frame(cbind(rootdep11,rootdep21,rootdep31,rootdep41))

write.csv(ponrootdeps, file = "ponderosa/pon max root depths.csv")







#curve fitting of 1 sample####################################################################################
x <- wide_ponroot1$NMVC.VCPP.1.2
y <- wide_ponroot1$depth

plot(x,y, pch = 19)
fit <- lm(y~x)
fit2 <- lm(y~log(x))
fit3 <- lm(y~poly(x,2,raw=TRUE))
fit4 <- lm(y~poly(x,3,raw=TRUE))

xx <- seq(0,4.1, length=32)
plot(x,y, pch = 19, xlim=c(0,4.1),ylim = c(-200,0),xlab = "Diameter (mm)", 
     ylab = "Rooting depth (mm)",xaxs="i",yaxs="i")
lines(xx, predict(fit, data.frame(x=xx)), col = 'red', lwd=2)
lines(xx, predict(fit2, data.frame(x=xx)), col = 'green', lwd=2)
lines(xx, predict(fit3, data.frame(x=xx)), col = 'blue', lwd=2)
lines(xx, predict(fit4, data.frame(x=xx)), col = 'orange', lwd=2)
text(1.5,10,expression(paste("Ponderosa pine: NMVC-PP-1")), col = "black", cex = 2, xpd=TRUE)
text(3,-110,expression(paste("RD = a*D+b")), col = "red", cex = 1.1)
text(3,-130,expression(paste("RD = a*log(D)+b")), col = "green", cex = 1.1)
text(3,-150,expression(paste("RD = a*D"^"2"*"+b*D+c")), col = "blue", cex = 1.1)
text(3,-170,expression(paste("RD = a*D"^"3"*"+b*D"^"2"*"+c*D+d")), col = "orange", cex = 1.1)

summary(fit) #through fit4 if wanted
coef(fit) #through fit 4 if wanted











#leaf, stem, root areas####################################################################################
leafarea <- read.csv(file = "areas/leaf_areas.csv")
stemarea <- read.csv(file = "areas/stem_areas.csv")
rootarea <- read.csv(file = "areas/root_areas.csv")


leafarea1 <- aggregate(area~location+forest+plotnum+samplenum+spp,data=leafarea,sum)
leafarea2 <- aggregate(area~location+spp,data=leafarea1,mean)
leafarea2 <- leafarea2[order(leafarea2$location),]
leafarea3 <- cbind(leafarea2,stddevl1)

stddevl <- leafarea1 %>% 
  group_by(location, spp) %>% 
  summarise(std = sd(area))
stddevl1 <- stddevl$std


stemarea1 <- aggregate(area~location+forest+plotnum+samplenum+spp,data=stemarea,sum)
stemarea2 <- aggregate(area~location+spp,data=stemarea1,mean)
stemarea2 <- stemarea2[order(stemarea2$location),]
stemarea3 <- cbind(leafarea2,stddevs1)

stddevs <- stemarea1 %>% 
  group_by(location, spp) %>% 
  summarise(std = sd(area))
stddevs1 <- stddevs$std


rootarea1 <- aggregate(area~location+forest+plotnum+samplenum+spp,data=rootarea,sum)
rootarea2 <- aggregate(area~location+spp,data=rootarea1,mean)
rootarea2 <- rootarea2[order(rootarea2$location),]
rootarea3 <- cbind(rootarea2,stddevr1)

stddevr <- rootarea1 %>% 
  group_by(location, spp) %>% 
  summarise(std = sd(area))

stddevr1 <- stddevr$std


write.csv(leafarea1, file = "areas/indleafarea.csv")
write.csv(stemarea1, file = "areas/indstemarea.csv")
write.csv(rootarea1, file = "areas/indrootarea.csv")








#tree rings###################################################################################################

rings <- read.csv(file = "treerings.csv")
rings <- subset(rings, Species != "Juniper") #if excluding juniper rings
rings$Species <- gsub("Pon", "Ponderosa pine", rings$Species)
rings$Species <- gsub("Pin", "Piñon pine", rings$Species)
rings$Species <- factor(rings$Species, levels = c("Juniper", "Piñon pine", "Ponderosa pine", "Fir"))


pinrings <- subset(rings, Species == c("Piñon pine"))

ponrings <- subset(rings, Species == "Ponderosa pine")

firrings <- subset(rings, Species == "Fir")

junrings <- subset(rings, Species == "Juniper")




#figure 2####################################################################################
#linear correlations of ring count to height and diameter and height to diameter


plot.new()
fig2 <- paste("figure 2.pdf",sep='')
pdf(fig2,width = 20,height = 60) 
layout(matrix(c(1,2,3),3,1, byrow = TRUE))
par(cex=1.5,cex.lab=1.5,mar = c(5, 10, 1.5, 2) + 0.3)


color2a <- c("orangered2","orange3","orchid3","palegreen4")[rings$Species]
reg2apjun <-lm(height ~ finalage, data = junrings)
reg2appin <-lm(height ~ finalage, data = pinrings)
reg2appon <-lm(height ~ finalage, data = ponrings)
reg2apfir <-lm(height ~ finalage, data = firrings)

plot(rings$finalage, rings$height, col=color2a,
     pch = c(16, 17, 15, 3, 7)[as.numeric(rings$core)], cex = 2,
     xlab = "", ylab = "", cex.lab=3,
     axes=FALSE, xaxs = "i", yaxs="i", xlim = c(0,43.5), ylim = c(0,900))
title(xlab = "Ring count", line = 4, cex.lab = 3)
axis(side = 1, at=c(0,10,20,30,40), labels = c(0,10,20,30,40), cex.axis=3, lwd = 2,
     mgp = c(3, 1.75, 0))
axis(side = 2, at=c(0,150,300,450,600), labels = c(0,150,300,450,600), cex.axis=3, las=2, lwd = 2,
     mgp = c(3, 0.75, 0))
#text(par("usr")[2] - 12, 0, srt=0, adj = c(1,1), labels = "Ring count", xpd = TRUE)
text(par("usr")[1] - 5.5, 15, srt=90, adj = -0.8, labels = "Height (mm)", xpd = TRUE, cex = 3)
abline(reg2apjun, col="orangered2", lwd=4)
abline(reg2appin, col="orange3", lwd=4)
abline(reg2appon, col="orchid3", lwd=4)
abline(reg2apfir, col="palegreen4", lwd=4)
rect(32, 650, 43.5, 920, border = FALSE, col = "white")
rect(40.5, 0, 43.5, 650, border = FALSE, col = "white")
text(11,860, labels = "Juniper: H = 6.6*RC + 49.5", cex = 3, xpd = TRUE)
text(34,860, labels = expression(paste("R"^"2"*" = 0.47")), cex = 3)
text(11,800, labels = "Piñon pine: H = 13*RC + 55", cex = 3, xpd = TRUE)
text(34,800, labels = expression(paste("R"^"2"*" = 0.72")), cex = 3)
text(11,740, labels = "Ponderosa pine: H = 16*RC + 56", cex = 3, xpd = TRUE)
text(34,740, labels = expression(paste("R"^"2"*" = 0.28")), cex = 3)
text(11,680, labels = "Fir: H = 5.2*RC + 52", cex = 3, xpd = TRUE)
text(34,680, labels = expression(paste("R"^"2"*" = 0.48")), cex = 3)
text(43.5,860, labels = "a", font = 2, cex = 5, xpd = TRUE)



reg2bpjun <-lm(diameter ~ 0 + finalage, data = junrings)
reg2bppin <-lm(diameter ~ finalage, data = pinrings)
reg2bppon <-lm(diameter ~ finalage, data = ponrings)
reg2bpfir <-lm(diameter ~ finalage, data = firrings)

plot(rings$finalage, rings$diameter, col=color2a,
     pch = c(16, 17, 15, 3, 7)[as.numeric(rings$core)], cex = 2,
     xlab = "", ylab = "",
     axes=FALSE, xaxs = "i", yaxs="i", xlim = c(0,43.5), ylim = c(0,36),
     cex.lab=3)
title(xlab = "Ring count", line = 4, cex.lab = 3)
axis(side = 1, at=c(0,10,20,30,40), labels = c(0,10,20,30,40), cex.axis=3, lwd=2,
     mgp = c(3, 1.75, 0))
axis(side = 2, at=c(0,5,10,15,20,25), labels = c(0,5,10,15,20,25), cex.axis=3, las=2, lwd=2,
     mgp = c(3, 0.75, 0))
text(par("usr")[1] - 5.5, 0, srt=90, adj = -0.6, labels = "Diameter (mm)", xpd = TRUE, cex = 3)
abline(reg2bpjun, col="orangered2", lwd=4)
abline(reg2bppin, col="orange3", lwd=4)
abline(reg2bppon, col="orchid3", lwd=4)
abline(reg2bpfir, col="palegreen4", lwd=4)
rect(40.5, 0, 43.5, 27, border = FALSE, col = "white")
text(11,34.4, labels = "Juniper: D = 0.25*RC", cex=3, xpd = TRUE)
text(34,34.4, labels = expression(paste("R"^"2"*" = 0.94")), cex=3)
text(11,32, labels = "Piñon pine: D = 0.41*RC + 2.1", cex=3, xpd = TRUE)
text(34,32, labels = expression(paste("R"^"2"*" = 0.85")), cex=3)
text(11,29.65, labels = "Ponderosa pine: D = 0.48*RC + 2.2", cex=3, xpd = TRUE)
text(34,29.65, labels = expression(paste("R"^"2"*" = 0.31")), cex=3)
text(11,27.25, labels = "Fir: D = 0.2*RC + 1.3", cex=3, xpd = TRUE)
text(34,27.25, labels = expression(paste("R"^"2"*" = 0.77")), cex=3)
text(43.5,34.5, labels = "b", font = 2, cex = 5, xpd = TRUE)




reg2cpjun <-lm(height ~ diameter, data = junrings)
reg2cppin <-lm(height ~ 0 + diameter, data = pinrings)
reg2cppon <-lm(height ~ 0 + diameter, data = ponrings)
reg2cpfir <-lm(height ~ diameter, data = firrings)

plot(rings$diameter, rings$height, col=color2a,
     pch = c(16, 17, 15, 3, 7)[as.numeric(rings$core)], cex = 2,
     xlab = "", ylab = "",
     axes=FALSE, xaxs = "i", yaxs="i", xlim = c(0,26.5), ylim = c(0,900),
     cex.lab=3)
title(xlab = "Diameter (mm)", line = 4, cex.lab = 3)
axis(side = 1, at=c(0,5,10,15,20,25), labels = c(0,5,10,15,20,25), cex.axis=3, lwd=2, 
     mgp = c(3, 1.75, 0))
axis(side = 2, at=c(0,150,300,450,600), labels = c(0,150,300,450,600), cex.axis=3, las=2, lwd=2,
     mgp = c(3, 0.75, 0))
text(par("usr")[1] - 3.5, 15, srt=90, adj = -0.8, labels = "Height (mm)", xpd = TRUE, cex = 3)
abline(reg2cpjun, col="orangered2", lwd=4)
abline(reg2cppin, col="orange3", lwd=4)
abline(reg2cppon, col="orchid3", lwd=4)
abline(reg2cpfir, col="palegreen4", lwd=4)
rect(17, 650, 28, 920, border = FALSE, col = "white")
rect(25.5, 0, 26.5, 650, border = FALSE, col = "white")
text(7,860, labels = "Juniper: H = 22*D + 68", cex=3, xpd = TRUE)
text(20.7,860, labels = expression(paste("R"^"2"*" = 0.74")), cex=3)
text(7,800, labels = "Piñon pine: H = 29.2*D", cex=3, xpd = TRUE)
text(20.7,800, labels = expression(paste("R"^"2"*" = 0.95")), cex=3)
text(7,740, labels = "Ponderosa pine: H = 31*D", cex=3, xpd = TRUE)
text(20.7,740, labels = expression(paste("R"^"2"*" = 0.91")), cex=3)
text(7,680, labels = "Fir: H = 26*D + 22", cex=3, xpd = TRUE)
text(20.7,680, labels = expression(paste("R"^"2"*" = 0.56")), cex=3)
text(26.25,860, labels = "c", font = 2, cex = 5, xpd = TRUE)
legend(14.2, 220, legend=levels(rings$Species), col = c("orangered2","orange3","orchid3","palegreen4"), cex = 2.7, pch = 16,
       lty=1, lwd=3.5, bg='white', box.lwd = 2)
legend(20, 470, legend=levels(rings$core), col = "black", cex = 2.7, pch = c(16, 17, 15, 3, 7),bg='white',
       box.lwd = 2)

dev.off()











#figure 3###############################################################################################
#difference in height and diameter between different microsite locations
install.packages('gapminder')
install.packages('forcats')
library(gapminder)
library(forcats)

alltrees <- read.csv(file = "measuredexcavated3.csv")
alltrees$leafarea <- alltrees$leafarea + (alltrees$leafarea*0.31) #add 31% to dry leaf area to account for wet
alltrees1 <- read.csv(file = "measuredexcavated1.csv")
measured <- read.csv(file = "measuredtrees.csv")


#juniper##############################################################################
classjun <- subset(measured, Species == "Juniper")
classjun$hclass <- factor(classjun$hclass, levels = c("<150", "150-300", "300+"))
classjun$covered <- gsub("Y", "Sheltered", classjun$covered)
classjun$covered <- gsub("N", "Unsheltered", classjun$covered)

classjun1 <- subset(alltrees1, Species == "Juniper")


classjun150c <- subset(classjun, hclass == "<150" & covered == "Covered")
classjun150u <- subset(classjun, hclass == "<150" & covered == "Uncovered")
classjun150a <- subset(classjun1, hclass == "<150")


classjun150300c <- subset(classjun, hclass == "150-300" & covered == "Covered")
classjun150300u <- subset(classjun, hclass == "150-300" & covered == "Uncovered")
classjun150300a <- subset(classjun1, hclass == "150-300")


classjun300c <- subset(classjun, hclass == "300+" & covered == "Covered")
classjun300u <- subset(classjun, hclass == "300+" & covered == "Uncovered")
classjun300a <- subset(classjun1, hclass == "300+")

classjun$hclass <- gsub("<150", "< 150", classjun$hclass)
classjun1$hclass <- gsub("<150", "< 150", classjun1$hclass)

anova_test(height ~ hclass*Species, data = alltrees1)
model  <- lm(height ~ hclass*Species, data = alltrees1)
alltrees1 %>%
  group_by(hclass) %>%
  anova_test(height ~ Species, error = model)

#height
fig3a <- ggplot() +
  geom_boxplot(classjun, mapping= aes(x=hclass, y=height, fill=covered), width=0.4,
               position = position_dodge(1), outlier.shape = NA) +
  geom_boxplot(classjun, mapping= aes(x=hclass, y=height, fill=covered), width=0.4,
               position = position_dodge(1), outlier.shape = NA) +
  geom_boxplot(classjun1, mapping= aes(x=hclass, y=height, fill="All"), width=0.2, outlier.shape = NA) + 
  stat_boxplot(classjun, mapping= aes(x=hclass, y=height, fill=covered), geom = "errorbar", width=0.2,
               position = position_dodge(1)) +
  stat_boxplot(classjun1, mapping= aes(x=hclass, y=height), geom = "errorbar", width=0.1) +
  labs(x="", y="Height (mm)", title = "Juniper") +
  scale_fill_manual(values = c("orangered2", "red3", "violetred"), labels = c("Sheltered", "All", "Unsheltered")) +
  scale_y_continuous(expand = c(0,0), breaks = c(0,250,500,750,1000), limits = c(0,1000)) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size=26), legend.position = c(0.2,0.75),
        legend.background = element_rect(linetype = "solid", color = "black"), 
        axis.title.x = element_blank(), legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 28),
        axis.text = element_text(size=26), 
        axis.title = element_text(size=26)) +
  annotate(geom = "text", x=0.75, y=240, label="abc", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=1, y=242, label="a", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=1.25, y=240, label="n=4", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=1.75, y=400, label="efh", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=2, y=400, label="f", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=2.25, y=325, label="fl", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=2.75, y=680, label="mnpq", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=3, y=900, label="mq", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=3.25, y=770, label="mnqr", fontface='italic', color="black", size=10) 

#no sig diffs
t.test(classjun150u$height, classjun150c$height, alternative = "greater")
t.test(classjun150u$height, classjun150a$height, alternative = "greater")
t.test(classjun150a$height, classjun150c$height, alternative = "greater")

#no sig diffs
t.test(classjun150300c$height, classjun150300u$height, alternative = "greater")
t.test(classjun150300a$height, classjun150300u$height, alternative = "greater")
t.test(classjun150300c$height, classjun150300a$height, alternative = "greater")

#no sig diffs
t.test(classjun300u$height, classjun300c$height, alternative = "greater")
t.test(classjun300u$height, classjun300a$height, alternative = "greater")
t.test(classjun300a$height, classjun300c$height, alternative = "greater")




#diameter
fig3b <- ggplot() +
  geom_boxplot(classjun, mapping= aes(x=hclass, y=diameter, fill=covered), width=0.4,
               position = position_dodge(1), outlier.shape = NA) +
  geom_boxplot(classjun1, mapping= aes(x=hclass, y=diameter, fill="All"), width=0.2, outlier.shape = NA) + 
  stat_boxplot(classjun, mapping= aes(x=hclass, y=diameter, fill=covered), geom = "errorbar", width=0.2,
               position = position_dodge(1)) +
  stat_boxplot(classjun1, mapping= aes(x=hclass, y=diameter), geom = "errorbar", width=0.1) +
  labs(x="Height class (mm)", y="Diameter (mm)", title="Juniper") +
  scale_fill_manual(values = c("orangered2", "red3", "violetred"), labels = c("Sheltered", "All", "Unsheltered")) +
  scale_y_continuous(expand = c(0,0), breaks = c(0,10,20,30,40,50), limits = c(0,50)) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.position = "none", 
        axis.title.x = element_blank(), axis.text = element_text(size = 26), axis.title = element_text(size = 26),
        plot.title = element_text(size=28)) +
  annotate(geom = "text", x=0.75, y=9, label="ade", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=1, y=9.5, label="acde", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=1.25, y=15, label="n=4", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=1.75, y=15, label="ijlm", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=2, y=20, label="jko", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=2.25, y=17, label="op", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=2.75, y=37.5, label="qst", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=3, y=33.5, label="qst", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=3.25, y=31.5, label="qstv", fontface='italic', color="black", size=10) 


t.test(classjun150u$diameter, classjun150c$diameter, alternative = "greater") #p=0.54
t.test(classjun150u$diameter, classjun150a$diameter, alternative = "greater") #not sig
t.test(classjun150u$diameter, classjun150a$diameter, alternative = "greater") #not sig


t.test(classjun150300u$diameter, classjun150300c$diameter, alternative = "greater") #sig greater p=0.01
t.test(classjun150300u$diameter, classjun150300a$diameter, alternative = "greater") #p=0.054
t.test(classjun150300a$diameter, classjun150300c$diameter, alternative = "greater") #not sig

#not sig
t.test(classjun300u$diameter, classjun300c$diameter, alternative = "greater")
t.test(classjun300u$diameter, classjun300a$diameter, alternative = "greater")
t.test(classjun300a$diameter, classjun300c$diameter, alternative = "greater")





#pinon#####################################################################################
classpin <- subset(measured, Species == "Pinon")
classpin$hclass <- factor(classpin$hclass, levels = c("<150", "150-300", "300+"))
classpin$covered <- gsub("Y", "Sheltered", classpin$covered)
classpin$covered <- gsub("N", "Unsheltered", classpin$covered)


classpin1 <- subset(alltrees1, Species == "Piñon")

classpin150c <- subset(classpin, hclass == "<150" & covered == "Covered")
classpin150u <- subset(classpin, hclass == "<150" & covered == "Uncovered")
classpin150a <- subset(classpin1, hclass == "<150")


classpin150300c <- subset(classpin, hclass == "150-300" & covered == "Covered")
classpin150300u <- subset(classpin, hclass == "150-300" & covered == "Uncovered")
classpin150300a <- subset(classpin1, hclass == "150-300")


classpin300c <- subset(classpin, hclass == "300+" & covered == "Covered")
classpin300u <- subset(classpin, hclass == "300+" & covered == "Uncovered")
classpin300a <- subset(classpin1, hclass == "300+")

classpin$hclass <- gsub("<150", "< 150", classpin$hclass)
classpin1$hclass <- gsub("<150", "< 150", classpin1$hclass)


#height
fig3c <- ggplot() +
  geom_boxplot(classpin, mapping= aes(x=hclass, y=height, fill=covered), width=0.4,
               position = position_dodge(1), outlier.shape = NA) +
  geom_boxplot(classpin1, mapping= aes(x=hclass, y=height, fill="All"), width=0.2, outlier.shape = NA) + 
  stat_boxplot(classpin, mapping= aes(x=hclass, y=height, fill=covered), geom = "errorbar", width=0.2,
               position = position_dodge(1)) +
  stat_boxplot(classpin1, mapping= aes(x=hclass, y=height), geom = "errorbar", width=0.1) +
  labs(x="", y="Height (mm)", title = "Piñon pine") +
  scale_fill_manual(values = c("tan1", "orange3", "wheat2"), labels = c("Sheltered", "All", "Unsheltered")) +
  scale_y_continuous(expand = c(0,0), breaks = c(0,250,500,750,1000), limits = c(0,1000)) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size=26), legend.position = c(0.2,0.75),
        legend.background = element_rect(linetype = "solid", color = "black"), 
        axis.title.x = element_blank(), axis.text = element_text(size = 26), axis.title = element_text(size = 26),
        plot.title = element_text(size = 28)) +
  annotate(geom = "text", x=0.75, y=250, label="abc", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=1, y=250, label="abd", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=1.25, y=250, label="a", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=1.75, y=390, label="efg", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=2, y=390, label="ehil", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=2.25, y=390, label="ij", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=2.75, y=950, label="mnq", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=3, y=875, label="mnq", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=3.25, y=700, label="mnq", fontface='italic', color="black", size=10) 

#no sig diffs
t.test(classpin150u$height, classpin150c$height, alternative = "greater")
t.test(classpin150u$height, classpin150a$height, alternative = "greater")
t.test(classpin150a$height, classpin150c$height, alternative = "greater")


t.test(classpin150300u$height, classpin150300c$height, alternative = "greater") # sig p=0.04
t.test(classpin150300u$height, classpin150300a$height, alternative = "greater") #p=0.056
t.test(classpin150300c$height, classpin150300a$height, alternative = "greater") #not sig

#no sig diffs
t.test(classpin300c$height, classpin300u$height, alternative = "greater")
t.test(classpin300a$height, classpin300u$height, alternative = "greater")
t.test(classpin300c$height, classpin300a$height, alternative = "greater")




#diameter
fig3d <- ggplot() +
  geom_boxplot(classpin, mapping= aes(x=hclass, y=diameter, fill=covered), width=0.4,
               position = position_dodge(1), outlier.shape = NA) +
  geom_boxplot(classpin1, mapping= aes(x=hclass, y=diameter, fill="All"), width=0.2, outlier.shape = NA) +
  stat_boxplot(classpin, mapping= aes(x=hclass, y=diameter, fill=covered), geom = "errorbar", width=0.2,
               position = position_dodge(1)) +
  stat_boxplot(classpin1, mapping= aes(x=hclass, y=diameter), geom = "errorbar", width=0.1) +
  labs(x="Height class (mm)", y="Diameter (mm)", title = "Piñon pine") +
  scale_fill_manual(values = c("tan1", "orange3", "wheat2"), labels = c("Sheltered", "All", "Unsheltered")) +
  scale_y_continuous(expand = c(0,0), breaks = c(0,10,20,30,40,50), limits = c(0,50)) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.position = "none", plot.title = element_text(size = 28),
        axis.title.x = element_blank(), axis.text = element_text(size = 26), axis.title = element_text(size = 26)) +
  annotate(geom = "text", x=0.75, y=12.5, label="abce", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=1, y=12.5, label="bc", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=1.25, y=14, label="bfi", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=1.75, y=18, label="jm", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=2, y=20.5, label="jo", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=2.25, y=25, label="p", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=2.75, y=30, label="qtu", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=3, y=33, label="qtu", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=3.25, y=33, label="sv", fontface='italic', color="black", size=10) 

#not sig
t.test(classpin150u$diameter, classpin150c$diameter, alternative = "greater") 
t.test(classpin150u$diameter, classpin150a$diameter, alternative = "greater") 
t.test(classpin150u$diameter, classpin150a$diameter, alternative = "greater") 


t.test(classpin150300u$diameter, classpin150300c$diameter, alternative = "greater") #sig greater p=0.004
t.test(classpin150300u$diameter, classpin150300a$diameter, alternative = "greater") #sig greater p=0.007
t.test(classpin150300a$diameter, classpin150300c$diameter, alternative = "greater") #not sig


t.test(classpin300u$diameter, classpin300c$diameter, alternative = "greater") #p=0.008
t.test(classpin300u$diameter, classpin300a$diameter, alternative = "greater") #p=0.02
t.test(classpin300a$diameter, classpin300c$diameter, alternative = "greater") #not sig







#ponderosa###########################################################################################
classpon <- subset(measured, Species == "Ponderosa")
classpon$hclass <- factor(classpon$hclass, levels = c("<150", "150-300", "300+"))
classpon$covered <- gsub("Y", "Sheltered", classpon$covered)
classpon$covered <- gsub("N", "Unsheltered", classpon$covered)

classpon1 <- subset(alltrees1, Species == "Ponderosa")


classpon150c <- subset(classpon, hclass == "<150" & covered == "Covered")
classpon150u <- subset(classpon, hclass == "<150" & covered == "Uncovered")
classpon150a <- subset(classpon1, hclass == "<150")


classpon150300c <- subset(classpon, hclass == "150-300" & covered == "Covered")
classpon150300u <- subset(classpon, hclass == "150-300" & covered == "Uncovered")
classpon150300a <- subset(classpon1, hclass == "150-300")


classpon300c <- subset(classpon, hclass == "300+" & covered == "Covered")
classpon300u <- subset(classpon, hclass == "300+" & covered == "Uncovered")
classpon300a <- subset(classpon1, hclass == "300+")

classpon$hclass <- gsub("<150", "< 150", classpon$hclass)
classpon1$hclass <- gsub("<150", "< 150", classpon1$hclass)


#height
fig3e <- ggplot() +
  geom_boxplot(classpon, mapping= aes(x=hclass, y=height, fill=covered), width=0.4,
               position = position_dodge(1), outlier.shape = NA) +
  geom_boxplot(classpon1, mapping= aes(x=hclass, y=height, fill="All"), width=0.2, outlier.shape = NA) + 
  stat_boxplot(classpon, mapping= aes(x=hclass, y=height, fill=covered), geom = "errorbar", width=0.2,
               position = position_dodge(1)) +
  stat_boxplot(classpon1, mapping= aes(x=hclass, y=height), geom = "errorbar", width=0.1) +
  labs(x="", y="Height (mm)", title = "Ponderosa pine") +
  scale_fill_manual(values = c("orchid3", "darkviolet", "mediumpurple"), labels = c("Sheltered", "All", "Unsheltered")) +
  scale_y_continuous(expand = c(0,0), breaks = c(0,250,500,750,1000), limits = c(0,1000)) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size=26), legend.position = c(0.2,0.75),
        legend.background = element_rect(linetype = "solid", color = "black"),
        axis.title.x = element_blank(), axis.text = element_text(size = 26), axis.title = element_text(size = 26),
        plot.title = element_text(size = 28)) +
  annotate(geom = "text", x=0.75, y=250, label="n=3", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=1, y=250, label="c", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=1.25, y=250, label="cd", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=1.75, y=450, label="n=6", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=2, y=420, label="fk", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=2.25, y=420, label="fk", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=2.75, y=750, label="n=6", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=3, y=840, label="or", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=3.25, y=840, label="oq", fontface='italic', color="black", size=10)

#no sig diffs
t.test(classpon150c$height, classpon150u$height, alternative = "greater")
t.test(classpon150u$height, classpon150a$height, alternative = "greater")
t.test(classpon150c$height, classpon150a$height, alternative = "greater")

#no sig diffs
t.test(classpon150300c$height, classpon150300u$height, alternative = "greater")
t.test(classpon150300a$height, classpon150300u$height, alternative = "greater")
t.test(classpon150300c$height, classpon150300a$height, alternative = "greater")

#no sig diffs
t.test(classpon300c$height, classpon300u$height, alternative = "greater")
t.test(classpon300a$height, classpon300u$height, alternative = "greater")
t.test(classpon300c$height, classpon300a$height, alternative = "greater")



#diameter
fig3f <- ggplot() +
  geom_boxplot(classpon, mapping= aes(x=hclass, y=diameter, fill=covered), width=0.4,
               position = position_dodge(1), outlier.shape = NA) +
  geom_boxplot(classpon1, mapping= aes(x=hclass, y=diameter, fill="All"), width=0.2, outlier.shape = NA) + 
  stat_boxplot(classpon, mapping= aes(x=hclass, y=diameter, fill=covered), geom = "errorbar", width=0.2,
               position = position_dodge(1)) +
  stat_boxplot(classpon1, mapping= aes(x=hclass, y=diameter), geom = "errorbar", width=0.1) +
  labs(x="Height class (mm)", y="Diameter (mm)", title = "Ponderosa pine") +
  scale_fill_manual(values = c("orchid3", "darkviolet", "mediumpurple"), labels = c("Sheltered", "All", "Unsheltered")) +
  scale_y_continuous(expand = c(0,0), breaks = c(0,10,20,30,40,50), limits = c(0,50)) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.position = "none", plot.title = element_text(size = 28),
        axis.title.x = element_blank(), axis.text = element_text(size = 26), axis.title = element_text(size = 26)) +
  annotate(geom = "text", x=0.75, y=10, label="n=3", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=1, y=11.5, label="eh", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=1.25, y=11.5, label="efg", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=1.75, y=14, label="n=6", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=2, y=17, label="ln", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=2.25, y=17, label="ln", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=2.75, y=46, label="n=6", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=3, y=36, label="sw", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=3.25, y=36, label="stuv", fontface='italic', color="black", size=10)

#not sig
t.test(classpon150c$diameter, classpon150u$diameter, alternative = "greater") 
t.test(classpon150u$diameter, classpon150a$diameter, alternative = "greater") 
t.test(classpon150c$diameter, classpon150a$diameter, alternative = "greater") 

#not sig
t.test(classpon150300c$diameter, classpon150300u$diameter, alternative = "greater") 
t.test(classpon150300u$diameter, classpon150300a$diameter, alternative = "greater") 
t.test(classpon150300c$diameter, classpon150300a$diameter, alternative = "greater") 

#not sig
t.test(classpon300c$diameter, classpon300u$diameter, alternative = "greater") 
t.test(classpon300a$diameter, classpon300u$diameter, alternative = "greater") 
t.test(classpon300c$diameter, classpon300a$diameter, alternative = "greater") 




#fir########################################################################################
classfir <- subset(measured, Species == "Fir")
classfir$hclass <- factor(classfir$hclass, levels = c("<150", "150-300", "300+"))
classfir$covered <- gsub("Y", "Sheltered", classfir$covered)
classfir$covered <- gsub("N", "Unsheltered", classfir$covered)


classfir1 <- subset(alltrees1, Species == "Fir")


classfir150c <- subset(classfir, hclass == "<150" & covered == "Covered")
classfir150u <- subset(classfir, hclass == "<150" & covered == "Uncovered")
classfir150a <- subset(classfir1, hclass == "<150")


classfir150300c <- subset(classfir, hclass == "150-300" & covered == "Covered")
classfir150300u <- subset(classfir, hclass == "150-300" & covered == "Uncovered")
classfir150300a <- subset(classfir1, hclass == "150-300")


classfir300c <- subset(classfir, hclass == "300+" & covered == "Covered")
classfir300u <- subset(classfir, hclass == "300+" & covered == "Uncovered")
classfir300a <- subset(classfir1, hclass == "300+")

classfir$hclass <- gsub("<150", "< 150", classfir$hclass)
classfir1$hclass <- gsub("<150", "< 150", classfir1$hclass)

#height
fig3g <- ggplot() +
  geom_boxplot(classfir, mapping= aes(x=hclass, y=height, fill=covered), width=0.4,
               position = position_dodge(1), outlier.shape = NA) +
  geom_boxplot(classfir1, mapping= aes(x=hclass, y=height, fill="All"), width=0.2, outlier.shape = NA) +
  stat_boxplot(classfir, mapping= aes(x=hclass, y=height, fill=covered), geom = "errorbar", width=0.2,
               position = position_dodge(1)) +
  stat_boxplot(classfir1, mapping= aes(x=hclass, y=height), geom = "errorbar", width=0.1) +
  labs(x="Height class (mm)", y="Height (mm)", title = "Fir") +
  scale_fill_manual(values = c("palegreen3", "palegreen4", "gray60"), labels = c("Sheltered", "All", "Unsheltered")) +
  scale_y_continuous(expand = c(0,0), breaks = c(0,250,500,750,1000), limits = c(0,1000)) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size=26), legend.position = c(0.2,0.8),
        legend.background = element_rect(linetype = "solid", color = "black"), axis.text = element_text(size = 26), 
        axis.title = element_text(size = 26), plot.title = element_text(size = 28)) +
  annotate(geom = "text", x=0.75, y=240, label="a", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=1, y=230, label="bc", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=1.25, y=230, label="bc", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=1.75, y=390, label="ehj", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=2, y=400, label="ehjl", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=2.25, y=400, label="ghikl", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=2.75, y=480, label="nop", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=3, y=590, label="p", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=3.25, y=590, label="p", fontface='italic', color="black", size=10)
  


t.test(classfir150c$height, classfir150u$height, alternative = "greater") #p=0.02
t.test(classfir150a$height, classfir150u$height, alternative = "greater") #not sig
t.test(classfir150c$height, classfir150a$height, alternative = "greater") #p=0.03

#no sig diffs
t.test(classfir150300c$height, classfir150300u$height, alternative = "greater")
t.test(classfir150300a$height, classfir150300u$height, alternative = "greater")
t.test(classfir150300c$height, classfir150300a$height, alternative = "greater")

#no sig diffs
t.test(classfir300c$height, classfir300u$height, alternative = "greater")
t.test(classfir300a$height, classfir300u$height, alternative = "greater")
t.test(classfir300c$height, classfir300a$height, alternative = "greater")



#diameter
fig3h <- ggplot() +
  geom_boxplot(classfir, mapping= aes(x=hclass, y=diameter, fill=covered), width=0.4,
               position = position_dodge(1), outlier.shape = NA) +
  geom_boxplot(classfir1, mapping= aes(x=hclass, y=diameter, fill="All"), width=0.2, outlier.shape = NA) +
  stat_boxplot(classfir, mapping= aes(x=hclass, y=diameter, fill=covered), geom = "errorbar", width=0.2,
               position = position_dodge(1)) +
  stat_boxplot(classfir1, mapping= aes(x=hclass, y=diameter), geom = "errorbar", width=0.1) +
  labs(x="Height class (mm)", y="Diameter (mm)", title = "Fir") +
  scale_fill_manual(values = c("palegreen3", "palegreen4", "gray60"), labels = c("Sheltered", "All", "Unsheltered")) +
  scale_y_continuous(expand = c(0,0), breaks = c(0,10,20,30,40,50), limits = c(0,50)) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.position = "none", axis.text = element_text(size = 26), 
        axis.title = element_text(size = 26), plot.title = element_text(size = 28)) +
  annotate(geom = "text", x=0.75, y=12, label="dh", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=1, y=11, label="d", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=1.25, y=12.5, label="dgh", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=1.75, y=19, label="jkn", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=2, y=20, label="km", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=2.25, y=20, label="jk", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=2.75, y=19, label="r", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=3, y=37, label="tw", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=3.25, y=37, label="tvw", fontface='italic', color="black", size=10)

#not sig
t.test(classfir150u$diameter, classfir150c$diameter, alternative = "greater") 
t.test(classfir150u$diameter, classfir150a$diameter, alternative = "greater") 
t.test(classfir150c$diameter, classfir150a$diameter, alternative = "greater") 

#not sig
t.test(classfir150300u$diameter, classfir150300c$diameter, alternative = "greater") 
t.test(classfir150300u$diameter, classfir150300a$diameter, alternative = "greater") 
t.test(classfir150300a$diameter, classfir150300c$diameter, alternative = "greater") 


t.test(classfir300u$diameter, classfir300c$diameter, alternative = "greater") #p=0.002
t.test(classfir300u$diameter, classfir300a$diameter, alternative = "greater") 
t.test(classfir300a$diameter, classfir300c$diameter, alternative = "greater") #p=0.006




plot_grid(fig3a, fig3b, fig3c, fig3d, fig3e, fig3f, fig3g, fig3h, nrow = 4,
          labels = c("a", "b", "c", "d", "e", "f", "g", "h"),
          label_x = 0.95, label_y = 1, label_size = 40)





#comparison between species
#height
t.test(classjun150c$height, classpin150c$height, alternative = "greater") 
t.test(classjun150a$height, classpin150a$height, alternative = "greater") 
t.test(classjun150c$height, classpin150u$height, alternative = "less") 
t.test(classjun150c$height, classpin150a$height, alternative = "greater") 
t.test(classjun150a$height, classpin150u$height, alternative = "less")
t.test(classjun150a$height, classpin150c$height, alternative = "greater") 

t.test(classjun150a$height, classpon150a$height, alternative = "greater") #p < 0.05
t.test(classjun150c$height, classpon150u$height, alternative = "greater") 
t.test(classjun150c$height, classpon150a$height, alternative = "greater") 
t.test(classjun150a$height, classpon150u$height, alternative = "greater") #p < 0.05

t.test(classjun150c$height, classfir150c$height, alternative = "less") 
t.test(classjun150a$height, classfir150a$height, alternative = "greater") #p < 0.05
t.test(classjun150c$height, classfir150u$height, alternative = "greater") 
t.test(classjun150c$height, classfir150a$height, alternative = "greater") 
t.test(classjun150a$height, classfir150u$height, alternative = "greater") #p < 0.05
t.test(classjun150a$height, classfir150c$height, alternative = "less") 

t.test(classpin150a$height, classpon150a$height, alternative = "greater") #p < 0.05
t.test(classpin150u$height, classpon150u$height, alternative = "greater") #p < 0.05
t.test(classpin150c$height, classpon150u$height, alternative = "greater") 
t.test(classpin150c$height, classpon150a$height, alternative = "greater") 
t.test(classpin150a$height, classpon150u$height, alternative = "greater")
t.test(classpin150u$height, classpon150a$height, alternative = "greater") #p < 0.05

t.test(classpin150c$height, classfir150c$height, alternative = "less") 
t.test(classpin150a$height, classfir150a$height, alternative = "greater") 
t.test(classpin150u$height, classfir150u$height, alternative = "greater") #p < 0.05
t.test(classpin150c$height, classfir150u$height, alternative = "greater") 
t.test(classpin150c$height, classfir150a$height, alternative = "greater") 
t.test(classpin150a$height, classfir150u$height, alternative = "greater")
t.test(classpin150u$height, classfir150c$height, alternative = "greater") 
t.test(classpin150a$height, classfir150c$height, alternative = "less") 
t.test(classpin150u$height, classfir150a$height, alternative = "greater") #p < 0.05

t.test(classpon150a$height, classfir150a$height, alternative = "less") 
t.test(classpon150u$height, classfir150u$height, alternative = "greater")
t.test(classpon150a$height, classfir150u$height, alternative = "greater")
t.test(classpon150u$height, classfir150a$height, alternative = "greater")
t.test(classpon150u$height, classfir150c$height, alternative = "less") #p < 0.05
t.test(classpon150a$height, classfir150c$height, alternative = "less") #p < 0.05



t.test(classjun150300c$height, classpin150300c$height, alternative = "less") 
t.test(classjun150300a$height, classpin150300a$height, alternative = "less") #p < 0.05
t.test(classjun150300u$height, classpin150300u$height, alternative = "less") #p < 0.05
t.test(classjun150300c$height, classpin150300u$height, alternative = "less") #p < 0.05
t.test(classjun150300c$height, classpin150300a$height, alternative = "less") 
t.test(classjun150300a$height, classpin150300u$height, alternative = "less") #p < 0.05
t.test(classjun150300u$height, classpin150300c$height, alternative = "less") 
t.test(classjun150300a$height, classpin150300c$height, alternative = "less") 
t.test(classjun150300u$height, classpin150300a$height, alternative = "less")

t.test(classjun150300a$height, classpon150300a$height, alternative = "less") 
t.test(classjun150300u$height, classpon150300u$height, alternative = "less") 
t.test(classjun150300c$height, classpon150300u$height, alternative = "greater") 
t.test(classjun150300c$height, classpon150300a$height, alternative = "less") 
t.test(classjun150300a$height, classpon150300u$height, alternative = "less")
t.test(classjun150300u$height, classpon150300a$height, alternative = "less")

t.test(classjun150300c$height, classfir150300c$height, alternative = "less") 
t.test(classjun150300a$height, classfir150300a$height, alternative = "less") #p < 0.05
t.test(classjun150300u$height, classfir150300u$height, alternative = "less") 
t.test(classjun150300c$height, classfir150300u$height, alternative = "less") 
t.test(classjun150300c$height, classfir150300a$height, alternative = "less") 
t.test(classjun150300a$height, classfir150300u$height, alternative = "less") #p < 0.05
t.test(classjun150300u$height, classfir150300c$height, alternative = "less") #p < 0.05 
t.test(classjun150300a$height, classfir150300c$height, alternative = "less") #p < 0.05
t.test(classjun150300u$height, classfir150300a$height, alternative = "less")

t.test(classpin150300a$height, classpon150300a$height, alternative = "greater") #p < 0.05
t.test(classpin150300u$height, classpon150300u$height, alternative = "greater") #p < 0.05
t.test(classpin150300c$height, classpon150300u$height, alternative = "greater") 
t.test(classpin150300c$height, classpon150300a$height, alternative = "greater") 
t.test(classpin150300a$height, classpon150300u$height, alternative = "greater") #p < 0.05
t.test(classpin150300u$height, classpon150300a$height, alternative = "greater") #p < 0.05

t.test(classpin150300c$height, classfir150300c$height, alternative = "less") 
t.test(classpin150300a$height, classfir150300a$height, alternative = "less") 
t.test(classpin150300u$height, classfir150300u$height, alternative = "greater") 
t.test(classpin150300c$height, classfir150300u$height, alternative = "less") 
t.test(classpin150300c$height, classfir150300a$height, alternative = "less") 
t.test(classpin150300a$height, classfir150300u$height, alternative = "less")
t.test(classpin150300u$height, classfir150300c$height, alternative = "greater") 
t.test(classpin150300a$height, classfir150300c$height, alternative = "less") 
t.test(classpin150300u$height, classfir150300a$height, alternative = "greater")

t.test(classpon150300a$height, classfir150300a$height, alternative = "less") #p < 0.05
t.test(classpon150300u$height, classfir150300u$height, alternative = "less")
t.test(classpon150300u$height, classfir150300c$height, alternative = "less") #p < 0.05
t.test(classpon150300a$height, classfir150300c$height, alternative = "less") #p < 0.05
t.test(classpon150300a$height, classfir150300u$height, alternative = "less")
t.test(classpon150300u$height, classfir150300a$height, alternative = "less") #p < 0.05



t.test(classjun300c$height, classpin300c$height, alternative = "less") 
t.test(classjun300a$height, classpin300a$height, alternative = "less") 
t.test(classjun300u$height, classpin300u$height, alternative = "greater")
t.test(classjun300c$height, classpin300u$height, alternative = "less") 
t.test(classjun300c$height, classpin300a$height, alternative = "less") 
t.test(classjun300a$height, classpin300u$height, alternative = "less")
t.test(classjun300u$height, classpin300c$height, alternative = "greater") 
t.test(classjun300a$height, classpin300c$height, alternative = "less") 
t.test(classjun300u$height, classpin300a$height, alternative = "greater")

t.test(classjun300a$height, classpon300a$height, alternative = "less") #p < 0.05
t.test(classjun300u$height, classpon300u$height, alternative = "greater") 
t.test(classjun300c$height, classpon300u$height, alternative = "less") 
t.test(classjun300c$height, classpon300a$height, alternative = "less") #p < 0.05
t.test(classjun300a$height, classpon300u$height, alternative = "less")
t.test(classjun300u$height, classpon300a$height, alternative = "less")

t.test(classjun300c$height, classfir300c$height, alternative = "less") 
t.test(classjun300a$height, classfir300a$height, alternative = "greater") #p < 0.05
t.test(classjun300u$height, classfir300u$height, alternative = "greater") #p < 0.05
t.test(classjun300c$height, classfir300u$height, alternative = "greater") 
t.test(classjun300c$height, classfir300a$height, alternative = "greater") 
t.test(classjun300a$height, classfir300u$height, alternative = "greater") #p < 0.05
t.test(classjun300u$height, classfir300c$height, alternative = "greater") 
t.test(classjun300a$height, classfir300c$height, alternative = "greater") #p < 0.05
t.test(classjun300u$height, classfir300a$height, alternative = "greater") #p < 0.05

t.test(classpin300a$height, classpon300a$height, alternative = "less") #p < 0.05
t.test(classpin300u$height, classpon300u$height, alternative = "less") 
t.test(classpin300c$height, classpon300u$height, alternative = "less") 
t.test(classpin300c$height, classpon300a$height, alternative = "less") #p < 0.05
t.test(classpin300a$height, classpon300u$height, alternative = "less")
t.test(classpin300u$height, classpon300a$height, alternative = "less") #p < 0.05

t.test(classpin300c$height, classfir300c$height, alternative = "greater") 
t.test(classpin300a$height, classfir300a$height, alternative = "greater") #p < 0.05
t.test(classpin300u$height, classfir300u$height, alternative = "greater") #p < 0.05
t.test(classpin300c$height, classfir300u$height, alternative = "greater") #p < 0.05
t.test(classpin300c$height, classfir300a$height, alternative = "greater") #p < 0.05
t.test(classpin300a$height, classfir300u$height, alternative = "greater") #p < 0.05
t.test(classpin300u$height, classfir300c$height, alternative = "greater") 
t.test(classpin300a$height, classfir300c$height, alternative = "greater") 
t.test(classpin300u$height, classfir300a$height, alternative = "greater") #p < 0.05

t.test(classpon300a$height, classfir300a$height, alternative = "greater") #p < 0.05
t.test(classpon300u$height, classfir300u$height, alternative = "greater") #p < 0.05
t.test(classpon300u$height, classfir300c$height, alternative = "greater") 
t.test(classpon300a$height, classfir300c$height, alternative = "greater") 
t.test(classpon300a$height, classfir300u$height, alternative = "greater") #p < 0.05
t.test(classpon300u$height, classfir300a$height, alternative = "greater") #p < 0.05



t.test(classjun150300c$height, classpin150c$height, alternative = "greater") #p < 0.05 all
t.test(classjun150300a$height, classpin150a$height, alternative = "greater") 
t.test(classjun150300u$height, classpin150u$height, alternative = "greater")
t.test(classjun150300c$height, classpin150u$height, alternative = "greater") 
t.test(classjun150300c$height, classpin150a$height, alternative = "greater") 
t.test(classjun150300a$height, classpin150u$height, alternative = "greater")
t.test(classjun150300u$height, classpin150c$height, alternative = "greater") 
t.test(classjun150300a$height, classpin150c$height, alternative = "greater") 
t.test(classjun150300u$height, classpin150a$height, alternative = "greater")

t.test(classjun150300a$height, classpon150a$height, alternative = "greater") 
t.test(classjun150300u$height, classpon150u$height, alternative = "greater") 
t.test(classjun150300c$height, classpon150u$height, alternative = "greater") 
t.test(classjun150300c$height, classpon150a$height, alternative = "greater") 
t.test(classjun150300a$height, classpon150u$height, alternative = "greater")
t.test(classjun150300u$height, classpon150a$height, alternative = "greater")

t.test(classjun150300c$height, classfir150c$height, alternative = "greater") 
t.test(classjun150300a$height, classfir150a$height, alternative = "greater") 
t.test(classjun150300u$height, classfir150u$height, alternative = "greater") 
t.test(classjun150300c$height, classfir150u$height, alternative = "greater") 
t.test(classjun150300c$height, classfir150a$height, alternative = "greater") 
t.test(classjun150300a$height, classfir150u$height, alternative = "greater")
t.test(classjun150300u$height, classfir150c$height, alternative = "greater") 
t.test(classjun150300a$height, classfir150c$height, alternative = "greater") 
t.test(classjun150300u$height, classfir150a$height, alternative = "greater")
 
t.test(classpin150300a$height, classpon150a$height, alternative = "greater") 
t.test(classpin150300u$height, classpon150u$height, alternative = "greater") 
t.test(classpin150300c$height, classpon150u$height, alternative = "greater") 
t.test(classpin150300c$height, classpon150a$height, alternative = "greater") 
t.test(classpin150300a$height, classpon150u$height, alternative = "greater")
t.test(classpin150300u$height, classpon150a$height, alternative = "greater")

t.test(classpin150300c$height, classfir150c$height, alternative = "greater") 
t.test(classpin150300a$height, classfir150a$height, alternative = "greater") 
t.test(classpin150300u$height, classfir150u$height, alternative = "greater") 
t.test(classpin150300c$height, classfir150u$height, alternative = "greater") 
t.test(classpin150300c$height, classfir150a$height, alternative = "greater") 
t.test(classpin150300a$height, classfir150u$height, alternative = "greater")
t.test(classpin150300u$height, classfir150c$height, alternative = "greater") 
t.test(classpin150300a$height, classfir150c$height, alternative = "greater") 
t.test(classpin150300u$height, classfir150a$height, alternative = "greater")

t.test(classpon150300a$height, classfir150a$height, alternative = "greater") 
t.test(classpon150300u$height, classfir150u$height, alternative = "greater") 
t.test(classpon150300c$height, classfir150u$height, alternative = "greater") 
t.test(classpon150300c$height, classfir150a$height, alternative = "greater") 
t.test(classpon150300a$height, classfir150u$height, alternative = "greater")
t.test(classpon150300u$height, classfir150a$height, alternative = "greater")



t.test(classjun300c$height, classpin150c$height, alternative = "greater") #p < 0.05 all
t.test(classjun300a$height, classpin150a$height, alternative = "greater") 
t.test(classjun300u$height, classpin150u$height, alternative = "greater")
t.test(classjun300c$height, classpin150u$height, alternative = "greater") 
t.test(classjun300c$height, classpin150a$height, alternative = "greater") 
t.test(classjun300a$height, classpin150u$height, alternative = "greater")
t.test(classjun300u$height, classpin150c$height, alternative = "greater") 
t.test(classjun300a$height, classpin150c$height, alternative = "greater") 
t.test(classjun300u$height, classpin150a$height, alternative = "greater")

t.test(classjun300a$height, classpon150a$height, alternative = "greater") 
t.test(classjun300u$height, classpon150u$height, alternative = "greater") 
t.test(classjun300c$height, classpon150u$height, alternative = "greater") 
t.test(classjun300c$height, classpon150a$height, alternative = "greater") 
t.test(classjun300a$height, classpon150u$height, alternative = "greater")
t.test(classjun300u$height, classpon150a$height, alternative = "greater")

t.test(classjun300c$height, classfir150c$height, alternative = "greater") 
t.test(classjun300a$height, classfir150a$height, alternative = "greater") 
t.test(classjun300u$height, classfir150u$height, alternative = "greater") 
t.test(classjun300c$height, classfir150u$height, alternative = "greater") 
t.test(classjun300c$height, classfir150a$height, alternative = "greater") 
t.test(classjun300a$height, classfir150u$height, alternative = "greater")
t.test(classjun300u$height, classfir150c$height, alternative = "greater") 
t.test(classjun300a$height, classfir150c$height, alternative = "greater") 
t.test(classjun300u$height, classfir150a$height, alternative = "greater")

t.test(classpin300a$height, classpon150a$height, alternative = "greater") 
t.test(classpin300u$height, classpon150u$height, alternative = "greater") 
t.test(classpin300c$height, classpon150u$height, alternative = "greater") 
t.test(classpin300c$height, classpon150a$height, alternative = "greater") 
t.test(classpin300a$height, classpon150u$height, alternative = "greater") 
t.test(classpin300u$height, classpon150a$height, alternative = "greater")

t.test(classpin300c$height, classfir150c$height, alternative = "greater") 
t.test(classpin300a$height, classfir150a$height, alternative = "greater") 
t.test(classpin300u$height, classfir150u$height, alternative = "greater") 
t.test(classpin300c$height, classfir150u$height, alternative = "greater") 
t.test(classpin300c$height, classfir150a$height, alternative = "greater") 
t.test(classpin300a$height, classfir150u$height, alternative = "greater")
t.test(classpin300u$height, classfir150c$height, alternative = "greater") 
t.test(classpin300a$height, classfir150c$height, alternative = "greater") 
t.test(classpin300u$height, classfir150a$height, alternative = "greater")

t.test(classpon300a$height, classfir150a$height, alternative = "greater") 
t.test(classpon300u$height, classfir150u$height, alternative = "greater")
t.test(classpon300a$height, classfir150u$height, alternative = "greater")
t.test(classpon300u$height, classfir150c$height, alternative = "greater") 
t.test(classpon300a$height, classfir150c$height, alternative = "greater") 
t.test(classpon300u$height, classfir150a$height, alternative = "greater")



t.test(classjun300c$height, classpin150300c$height, alternative = "greater") #p < 0.05 all
t.test(classjun300a$height, classpin150300a$height, alternative = "greater") 
t.test(classjun300u$height, classpin150300u$height, alternative = "greater")
t.test(classjun300c$height, classpin150300u$height, alternative = "greater") 
t.test(classjun300c$height, classpin150300a$height, alternative = "greater") 
t.test(classjun300a$height, classpin150300u$height, alternative = "greater")
t.test(classjun300u$height, classpin150300c$height, alternative = "greater") 
t.test(classjun300a$height, classpin150300c$height, alternative = "greater") 
t.test(classjun300u$height, classpin150300a$height, alternative = "greater")
 
t.test(classjun300a$height, classpon150300a$height, alternative = "greater") 
t.test(classjun300u$height, classpon150300u$height, alternative = "greater")
t.test(classjun300a$height, classpon150300u$height, alternative = "greater")
t.test(classjun300u$height, classpon150300c$height, alternative = "greater") 
t.test(classjun300a$height, classpon150300c$height, alternative = "greater") 
t.test(classjun300u$height, classpon150300a$height, alternative = "greater")

t.test(classjun300c$height, classfir150300c$height, alternative = "greater") 
t.test(classjun300a$height, classfir150300a$height, alternative = "greater") 
t.test(classjun300u$height, classfir150300u$height, alternative = "greater") 
t.test(classjun300c$height, classfir150300u$height, alternative = "greater") 
t.test(classjun300c$height, classfir150300a$height, alternative = "greater") 
t.test(classjun300a$height, classfir150300u$height, alternative = "greater")
t.test(classjun300u$height, classfir150300c$height, alternative = "greater") 
t.test(classjun300a$height, classfir150300c$height, alternative = "greater") 
t.test(classjun300u$height, classfir150300a$height, alternative = "greater")
 
t.test(classpin300a$height, classpon150300a$height, alternative = "greater") 
t.test(classpin300u$height, classpon150300u$height, alternative = "greater") 
t.test(classpin300c$height, classpon150300u$height, alternative = "greater") 
t.test(classpin300c$height, classpon150300a$height, alternative = "greater") 
t.test(classpin300a$height, classpon150300u$height, alternative = "greater") 
t.test(classpin300u$height, classpon150300a$height, alternative = "greater")

t.test(classpin300c$height, classfir150300c$height, alternative = "greater") 
t.test(classpin300a$height, classfir150300a$height, alternative = "greater") 
t.test(classpin300u$height, classfir150300u$height, alternative = "greater") 
t.test(classpin300c$height, classfir150300u$height, alternative = "greater") 
t.test(classpin300c$height, classfir150300a$height, alternative = "greater") 
t.test(classpin300a$height, classfir150300u$height, alternative = "greater")
t.test(classpin300u$height, classfir150300c$height, alternative = "greater") 
t.test(classpin300a$height, classfir150300c$height, alternative = "greater") 
t.test(classpin300u$height, classfir150300a$height, alternative = "greater")

t.test(classpon300a$height, classfir150300a$height, alternative = "greater") 
t.test(classpon300u$height, classfir150300u$height, alternative = "greater")
t.test(classpon300a$height, classfir150300u$height, alternative = "greater")
t.test(classpon300u$height, classfir150300c$height, alternative = "greater") 
t.test(classpon300a$height, classfir150300c$height, alternative = "greater") 
t.test(classpon300u$height, classfir150300a$height, alternative = "greater")





#diameter
t.test(classjun150c$diameter, classpin150c$diameter, alternative = "less") 
t.test(classjun150a$diameter, classpin150a$diameter, alternative = "less") 
t.test(classjun150c$diameter, classpin150u$diameter, alternative = "less") #p < 0.05
t.test(classjun150c$diameter, classpin150a$diameter, alternative = "less") #p < 0.05
t.test(classjun150a$diameter, classpin150u$diameter, alternative = "less") #p < 0.05
t.test(classjun150a$diameter, classpin150c$diameter, alternative = "less") 

t.test(classjun150a$diameter, classpon150a$diameter, alternative = "less") 
t.test(classjun150c$diameter, classpon150u$diameter, alternative = "less") 
t.test(classjun150c$diameter, classpon150a$diameter, alternative = "less") 
t.test(classjun150a$diameter, classpon150u$diameter, alternative = "less")
t.test(classjun150u$diameter, classpon150a$diameter, alternative = "greater")

t.test(classjun150c$diameter, classfir150c$diameter, alternative = "greater") 
t.test(classjun150a$diameter, classfir150a$diameter, alternative = "greater") 
t.test(classjun150c$diameter, classfir150u$diameter, alternative = "less") 
t.test(classjun150c$diameter, classfir150a$diameter, alternative = "greater") 
t.test(classjun150a$diameter, classfir150u$diameter, alternative = "greater")
t.test(classjun150u$diameter, classfir150c$diameter, alternative = "greater") 
t.test(classjun150a$diameter, classfir150c$diameter, alternative = "greater") 
t.test(classjun150u$diameter, classfir150a$diameter, alternative = "greater")

t.test(classpin150a$diameter, classpon150a$diameter, alternative = "greater") #p < 0.05
t.test(classpin150u$diameter, classpon150u$diameter, alternative = "greater") 
t.test(classpin150c$diameter, classpon150u$diameter, alternative = "greater") 
t.test(classpin150c$diameter, classpon150a$diameter, alternative = "greater") 
t.test(classpin150a$diameter, classpon150u$diameter, alternative = "greater") #p < 0.05
t.test(classpin150u$diameter, classpon150a$diameter, alternative = "greater") #p < 0.05

t.test(classpin150c$diameter, classfir150c$diameter, alternative = "greater") #p < 0.05
t.test(classpin150a$diameter, classfir150a$diameter, alternative = "greater") #p < 0.05
t.test(classpin150u$diameter, classfir150u$diameter, alternative = "greater") #p < 0.05
t.test(classpin150c$diameter, classfir150u$diameter, alternative = "greater") #p < 0.05
t.test(classpin150c$diameter, classfir150a$diameter, alternative = "greater") #p < 0.05
t.test(classpin150a$diameter, classfir150u$diameter, alternative = "greater") #p < 0.05
t.test(classpin150u$diameter, classfir150c$diameter, alternative = "greater") #p < 0.05
t.test(classpin150a$diameter, classfir150c$diameter, alternative = "greater") #p < 0.05
t.test(classpin150u$diameter, classfir150a$diameter, alternative = "greater") #p < 0.05

t.test(classpon150a$diameter, classfir150a$diameter, alternative = "greater") #p < 0.05
t.test(classpon150u$diameter, classfir150u$diameter, alternative = "greater") 
t.test(classpon150u$diameter, classfir150c$diameter, alternative = "greater") #p < 0.05
t.test(classpon150a$diameter, classfir150c$diameter, alternative = "greater") 
t.test(classpon150a$diameter, classfir150u$diameter, alternative = "greater")
t.test(classpon150u$diameter, classfir150a$diameter, alternative = "greater") #p < 0.05



t.test(classjun150300c$diameter, classpin150300c$diameter, alternative = "less") 
t.test(classjun150300a$diameter, classpin150300a$diameter, alternative = "less")
t.test(classjun150300u$diameter, classpin150300u$diameter, alternative = "greater") 
t.test(classjun150300c$diameter, classpin150300u$diameter, alternative = "less") #p < 0.05
t.test(classjun150300c$diameter, classpin150300a$diameter, alternative = "less") 
t.test(classjun150300a$diameter, classpin150300u$diameter, alternative = "less") #p < 0.05
t.test(classjun150300u$diameter, classpin150300c$diameter, alternative = "greater") #p < 0.05
t.test(classjun150300a$diameter, classpin150300c$diameter, alternative = "greater") 
t.test(classjun150300u$diameter, classpin150300a$diameter, alternative = "greater")

t.test(classjun150300a$diameter, classpon150300a$diameter, alternative = "greater") #p < 0.05
t.test(classjun150300u$diameter, classpon150300u$diameter, alternative = "greater") #p < 0.05
t.test(classjun150300c$diameter, classpon150300u$diameter, alternative = "greater") 
t.test(classjun150300c$diameter, classpon150300a$diameter, alternative = "greater") 
t.test(classjun150300a$diameter, classpon150300u$diameter, alternative = "greater") #p < 0.05
t.test(classjun150300u$diameter, classpon150300a$diameter, alternative = "greater") #p < 0.05

t.test(classjun150300c$diameter, classfir150300c$diameter, alternative = "less") 
t.test(classjun150300a$diameter, classfir150300a$diameter, alternative = "greater") 
t.test(classjun150300u$diameter, classfir150300u$diameter, alternative = "greater") #p < 0.05
t.test(classjun150300c$diameter, classfir150300u$diameter, alternative = "less") 
t.test(classjun150300c$diameter, classfir150300a$diameter, alternative = "less") 
t.test(classjun150300a$diameter, classfir150300u$diameter, alternative = "greater")
t.test(classjun150300u$diameter, classfir150300c$diameter, alternative = "greater") #p < 0.05
t.test(classjun150300a$diameter, classfir150300c$diameter, alternative = "greater") 
t.test(classjun150300u$diameter, classfir150300a$diameter, alternative = "greater") #p < 0.05

t.test(classpin150300a$diameter, classpon150300a$diameter, alternative = "greater") #p < 0.05
t.test(classpin150300u$diameter, classpon150300u$diameter, alternative = "greater") #p < 0.05
t.test(classpin150300c$diameter, classpon150300u$diameter, alternative = "greater") #p < 0.05
t.test(classpin150300c$diameter, classpon150300a$diameter, alternative = "greater") #p < 0.05
t.test(classpin150300a$diameter, classpon150300u$diameter, alternative = "greater") #p < 0.05
t.test(classpin150300u$diameter, classpon150300a$diameter, alternative = "greater") #p < 0.05

t.test(classpin150300c$diameter, classfir150300c$diameter, alternative = "greater") 
t.test(classpin150300a$diameter, classfir150300a$diameter, alternative = "greater") #p < 0.05
t.test(classpin150300u$diameter, classfir150300u$diameter, alternative = "greater") #p < 0.05
t.test(classpin150300c$diameter, classfir150300u$diameter, alternative = "greater") 
t.test(classpin150300c$diameter, classfir150300a$diameter, alternative = "greater") 
t.test(classpin150300a$diameter, classfir150300u$diameter, alternative = "greater")
t.test(classpin150300u$diameter, classfir150300c$diameter, alternative = "greater") #p < 0.05
t.test(classpin150300a$diameter, classfir150300c$diameter, alternative = "greater") 
t.test(classpin150300u$diameter, classfir150300a$diameter, alternative = "greater") #p < 0.05

t.test(classpon150300a$diameter, classfir150300a$diameter, alternative = "less") #p < 0.05
t.test(classpon150300u$diameter, classfir150300u$diameter, alternative = "less") #p < 0.05
t.test(classpon150300u$diameter, classfir150300c$diameter, alternative = "less") 
t.test(classpon150300a$diameter, classfir150300c$diameter, alternative = "less") 
t.test(classpon150300a$diameter, classfir150300u$diameter, alternative = "less") #p < 0.05
t.test(classpon150300u$diameter, classfir150300a$diameter, alternative = "less") #p < 0.05



t.test(classjun300c$diameter, classpin300c$diameter, alternative = "greater") 
t.test(classjun300a$diameter, classpin300a$diameter, alternative = "greater") 
t.test(classjun300u$diameter, classpin300u$diameter, alternative = "less")
t.test(classjun300c$diameter, classpin300u$diameter, alternative = "less") 
t.test(classjun300c$diameter, classpin300a$diameter, alternative = "greater") 
t.test(classjun300a$diameter, classpin300u$diameter, alternative = "less")
t.test(classjun300u$diameter, classpin300c$diameter, alternative = "greater") 
t.test(classjun300a$diameter, classpin300c$diameter, alternative = "greater") 
t.test(classjun300u$diameter, classpin300a$diameter, alternative = "greater")

t.test(classjun300a$diameter, classpon300a$diameter, alternative = "less") 
t.test(classjun300u$diameter, classpon300u$diameter, alternative = "greater") 
t.test(classjun300c$diameter, classpon300u$diameter, alternative = "less") 
t.test(classjun300c$diameter, classpon300a$diameter, alternative = "less") 
t.test(classjun300a$diameter, classpon300u$diameter, alternative = "less")
t.test(classjun300u$diameter, classpon300a$diameter, alternative = "greater")

t.test(classjun300c$diameter, classfir300c$diameter, alternative = "greater") #p < 0.05
t.test(classjun300a$diameter, classfir300a$diameter, alternative = "greater") 
t.test(classjun300u$diameter, classfir300u$diameter, alternative = "greater") 
t.test(classjun300c$diameter, classfir300u$diameter, alternative = "less") 
t.test(classjun300c$diameter, classfir300a$diameter, alternative = "greater") 
t.test(classjun300a$diameter, classfir300u$diameter, alternative = "less")
t.test(classjun300u$diameter, classfir300c$diameter, alternative = "greater") #p < 0.05
t.test(classjun300a$diameter, classfir300c$diameter, alternative = "greater") #p < 0.05
t.test(classjun300u$diameter, classfir300a$diameter, alternative = "greater")

t.test(classpin300a$diameter, classpon300a$diameter, alternative = "less") #p < 0.05
t.test(classpin300u$diameter, classpon300u$diameter, alternative = "greater") 
t.test(classpin300c$diameter, classpon300u$diameter, alternative = "less") 
t.test(classpin300c$diameter, classpon300a$diameter, alternative = "less") #p < 0.05
t.test(classpin300a$diameter, classpon300u$diameter, alternative = "less")
t.test(classpin300u$diameter, classpon300a$diameter, alternative = "greater")

t.test(classpin300c$diameter, classfir300c$diameter, alternative = "greater") #p < 0.05
t.test(classpin300a$diameter, classfir300a$diameter, alternative = "less") 
t.test(classpin300u$diameter, classfir300u$diameter, alternative = "greater") 
t.test(classpin300c$diameter, classfir300u$diameter, alternative = "less") 
t.test(classpin300c$diameter, classfir300a$diameter, alternative = "less") 
t.test(classpin300a$diameter, classfir300u$diameter, alternative = "less")
t.test(classpin300u$diameter, classfir300c$diameter, alternative = "greater") #p < 0.05
t.test(classpin300a$diameter, classfir300c$diameter, alternative = "greater") #p < 0.05
t.test(classpin300u$diameter, classfir300a$diameter, alternative = "greater") #p < 0.05

t.test(classpon300a$diameter, classfir300a$diameter, alternative = "greater") 
t.test(classpon300u$diameter, classfir300u$diameter, alternative = "greater") 
t.test(classpon300u$diameter, classfir300c$diameter, alternative = "greater") #p < 0.05
t.test(classpon300a$diameter, classfir300c$diameter, alternative = "greater") #p < 0.05
t.test(classpon300a$diameter, classfir300u$diameter, alternative = "greater")
t.test(classpon300u$diameter, classfir300a$diameter, alternative = "greater")



t.test(classjun150300c$diameter, classpin150c$diameter, alternative = "greater") #p < 0.05 all
t.test(classjun150300a$diameter, classpin150a$diameter, alternative = "greater") 
t.test(classjun150300u$diameter, classpin150u$diameter, alternative = "greater")
t.test(classjun150300c$diameter, classpin150u$diameter, alternative = "greater") #insig 
t.test(classjun150300c$diameter, classpin150a$diameter, alternative = "greater") 
t.test(classjun150300a$diameter, classpin150u$diameter, alternative = "greater")
t.test(classjun150300u$diameter, classpin150c$diameter, alternative = "greater") 
t.test(classjun150300a$diameter, classpin150c$diameter, alternative = "greater") 
t.test(classjun150300u$diameter, classpin150a$diameter, alternative = "greater")

t.test(classjun150300a$diameter, classpon150a$diameter, alternative = "greater") 
t.test(classjun150300u$diameter, classpon150u$diameter, alternative = "greater") 
t.test(classjun150300c$diameter, classpon150u$diameter, alternative = "greater") 
t.test(classjun150300c$diameter, classpon150a$diameter, alternative = "greater") 
t.test(classjun150300a$diameter, classpon150u$diameter, alternative = "greater")
t.test(classjun150300u$diameter, classpon150a$diameter, alternative = "greater")

t.test(classjun150300c$diameter, classfir150c$diameter, alternative = "greater") 
t.test(classjun150300a$diameter, classfir150a$diameter, alternative = "greater") 
t.test(classjun150300u$diameter, classfir150u$diameter, alternative = "greater") 
t.test(classjun150300c$diameter, classfir150u$diameter, alternative = "greater") 
t.test(classjun150300c$diameter, classfir150a$diameter, alternative = "greater") 
t.test(classjun150300a$diameter, classfir150u$diameter, alternative = "greater")
t.test(classjun150300u$diameter, classfir150c$diameter, alternative = "greater") 
t.test(classjun150300a$diameter, classfir150c$diameter, alternative = "greater") 
t.test(classjun150300u$diameter, classfir150a$diameter, alternative = "greater")

t.test(classpin150300a$diameter, classpon150a$diameter, alternative = "greater") 
t.test(classpin150300u$diameter, classpon150u$diameter, alternative = "greater") 
t.test(classpin150300c$diameter, classpon150u$diameter, alternative = "greater") 
t.test(classpin150300c$diameter, classpon150a$diameter, alternative = "greater") 
t.test(classpin150300a$diameter, classpon150u$diameter, alternative = "greater")
t.test(classpin150300u$diameter, classpon150a$diameter, alternative = "greater")

t.test(classpin150300c$diameter, classfir150c$diameter, alternative = "greater") 
t.test(classpin150300a$diameter, classfir150a$diameter, alternative = "greater") 
t.test(classpin150300u$diameter, classfir150u$diameter, alternative = "greater") 
t.test(classpin150300c$diameter, classfir150u$diameter, alternative = "greater") 
t.test(classpin150300c$diameter, classfir150a$diameter, alternative = "greater") 
t.test(classpin150300a$diameter, classfir150u$diameter, alternative = "greater")
t.test(classpin150300u$diameter, classfir150c$diameter, alternative = "greater") 
t.test(classpin150300a$diameter, classfir150c$diameter, alternative = "greater") 
t.test(classpin150300u$diameter, classfir150a$diameter, alternative = "greater")

t.test(classpon150300a$diameter, classfir150a$diameter, alternative = "greater") 
t.test(classpon150300u$diameter, classfir150u$diameter, alternative = "greater") 
t.test(classpon150300c$diameter, classfir150u$diameter, alternative = "greater") 
t.test(classpon150300c$diameter, classfir150a$diameter, alternative = "greater") 
t.test(classpon150300a$diameter, classfir150u$diameter, alternative = "greater")
t.test(classpon150300u$diameter, classfir150a$diameter, alternative = "greater")



t.test(classjun300c$diameter, classpin150c$diameter, alternative = "greater") #p < 0.05 all
t.test(classjun300a$diameter, classpin150a$diameter, alternative = "greater") 
t.test(classjun300u$diameter, classpin150u$diameter, alternative = "greater")
t.test(classjun300c$diameter, classpin150u$diameter, alternative = "greater") 
t.test(classjun300c$diameter, classpin150a$diameter, alternative = "greater") 
t.test(classjun300a$diameter, classpin150u$diameter, alternative = "greater")
t.test(classjun300u$diameter, classpin150c$diameter, alternative = "greater") 
t.test(classjun300a$diameter, classpin150c$diameter, alternative = "greater") 
t.test(classjun300u$diameter, classpin150a$diameter, alternative = "greater")

t.test(classjun300a$diameter, classpon150a$diameter, alternative = "greater") 
t.test(classjun300u$diameter, classpon150u$diameter, alternative = "greater") 
t.test(classjun300c$diameter, classpon150u$diameter, alternative = "greater") 
t.test(classjun300c$diameter, classpon150a$diameter, alternative = "greater") 
t.test(classjun300a$diameter, classpon150u$diameter, alternative = "greater")
t.test(classjun300u$diameter, classpon150a$diameter, alternative = "greater")

t.test(classjun300c$diameter, classfir150c$diameter, alternative = "greater") 
t.test(classjun300a$diameter, classfir150a$diameter, alternative = "greater") 
t.test(classjun300u$diameter, classfir150u$diameter, alternative = "greater") 
t.test(classjun300c$diameter, classfir150u$diameter, alternative = "greater") 
t.test(classjun300c$diameter, classfir150a$diameter, alternative = "greater") 
t.test(classjun300a$diameter, classfir150u$diameter, alternative = "greater")
t.test(classjun300u$diameter, classfir150c$diameter, alternative = "greater") 
t.test(classjun300a$diameter, classfir150c$diameter, alternative = "greater") 
t.test(classjun300u$diameter, classfir150a$diameter, alternative = "greater")

t.test(classpin300a$diameter, classpon150a$diameter, alternative = "greater") 
t.test(classpin300u$diameter, classpon150u$diameter, alternative = "greater") 
t.test(classpin300c$diameter, classpon150u$diameter, alternative = "greater") 
t.test(classpin300c$diameter, classpon150a$diameter, alternative = "greater") 
t.test(classpin300a$diameter, classpon150u$diameter, alternative = "greater") 
t.test(classpin300u$diameter, classpon150a$diameter, alternative = "greater")

t.test(classpin300c$diameter, classfir150c$diameter, alternative = "greater") 
t.test(classpin300a$diameter, classfir150a$diameter, alternative = "greater") 
t.test(classpin300u$diameter, classfir150u$diameter, alternative = "greater") 
t.test(classpin300c$diameter, classfir150u$diameter, alternative = "greater") 
t.test(classpin300c$diameter, classfir150a$diameter, alternative = "greater") 
t.test(classpin300a$diameter, classfir150u$diameter, alternative = "greater")
t.test(classpin300u$diameter, classfir150c$diameter, alternative = "greater") 
t.test(classpin300a$diameter, classfir150c$diameter, alternative = "greater") 
t.test(classpin300u$diameter, classfir150a$diameter, alternative = "greater")

t.test(classpon300a$diameter, classfir150a$diameter, alternative = "greater") 
t.test(classpon300u$diameter, classfir150u$diameter, alternative = "greater")
t.test(classpon300a$diameter, classfir150u$diameter, alternative = "greater")
t.test(classpon300u$diameter, classfir150c$diameter, alternative = "greater") 
t.test(classpon300a$diameter, classfir150c$diameter, alternative = "greater") 
t.test(classpon300u$diameter, classfir150a$diameter, alternative = "greater")



t.test(classjun300c$diameter, classpin150300c$diameter, alternative = "greater") #p < 0.05 all
t.test(classjun300a$diameter, classpin150300a$diameter, alternative = "greater") 
t.test(classjun300u$diameter, classpin150300u$diameter, alternative = "greater")
t.test(classjun300c$diameter, classpin150300u$diameter, alternative = "greater") 
t.test(classjun300c$diameter, classpin150300a$diameter, alternative = "greater") 
t.test(classjun300a$diameter, classpin150300u$diameter, alternative = "greater")
t.test(classjun300u$diameter, classpin150300c$diameter, alternative = "greater") 
t.test(classjun300a$diameter, classpin150300c$diameter, alternative = "greater") 
t.test(classjun300u$diameter, classpin150300a$diameter, alternative = "greater")

t.test(classjun300a$diameter, classpon150300a$diameter, alternative = "greater") 
t.test(classjun300u$diameter, classpon150300u$diameter, alternative = "greater")
t.test(classjun300a$diameter, classpon150300u$diameter, alternative = "greater")
t.test(classjun300u$diameter, classpon150300c$diameter, alternative = "greater") 
t.test(classjun300a$diameter, classpon150300c$diameter, alternative = "greater") 
t.test(classjun300u$diameter, classpon150300a$diameter, alternative = "greater")

t.test(classjun300c$diameter, classfir150300c$diameter, alternative = "greater") 
t.test(classjun300a$diameter, classfir150300a$diameter, alternative = "greater") 
t.test(classjun300u$diameter, classfir150300u$diameter, alternative = "greater") 
t.test(classjun300c$diameter, classfir150300u$diameter, alternative = "greater") 
t.test(classjun300c$diameter, classfir150300a$diameter, alternative = "greater") 
t.test(classjun300a$diameter, classfir150300u$diameter, alternative = "greater")
t.test(classjun300u$diameter, classfir150300c$diameter, alternative = "greater") 
t.test(classjun300a$diameter, classfir150300c$diameter, alternative = "greater") 
t.test(classjun300u$diameter, classfir150300a$diameter, alternative = "greater")

t.test(classpin300a$diameter, classpon150300a$diameter, alternative = "greater") 
t.test(classpin300u$diameter, classpon150300u$diameter, alternative = "greater") 
t.test(classpin300c$diameter, classpon150300u$diameter, alternative = "greater") 
t.test(classpin300c$diameter, classpon150300a$diameter, alternative = "greater") 
t.test(classpin300a$diameter, classpon150300u$diameter, alternative = "greater") 
t.test(classpin300u$diameter, classpon150300a$diameter, alternative = "greater")

t.test(classpin300c$diameter, classfir150300c$diameter, alternative = "greater") 
t.test(classpin300a$diameter, classfir150300a$diameter, alternative = "greater") 
t.test(classpin300u$diameter, classfir150300u$diameter, alternative = "greater") 
t.test(classpin300c$diameter, classfir150300u$diameter, alternative = "greater") 
t.test(classpin300c$diameter, classfir150300a$diameter, alternative = "greater") 
t.test(classpin300a$diameter, classfir150300u$diameter, alternative = "greater")
t.test(classpin300u$diameter, classfir150300c$diameter, alternative = "greater") 
t.test(classpin300a$diameter, classfir150300c$diameter, alternative = "greater") 
t.test(classpin300u$diameter, classfir150300a$diameter, alternative = "greater")

t.test(classpon300a$diameter, classfir150300a$diameter, alternative = "greater") 
t.test(classpon300u$diameter, classfir150300u$diameter, alternative = "greater")
t.test(classpon300a$diameter, classfir150300u$diameter, alternative = "greater")
t.test(classpon300u$diameter, classfir150300c$diameter, alternative = "greater") 
t.test(classpon300a$diameter, classfir150300c$diameter, alternative = "greater") 
t.test(classpon300u$diameter, classfir150300a$diameter, alternative = "greater")









#figure 4 and s2##########################################################################################
excavations <- read.csv(file = "excavations.csv")
excavations$leafarea <- excavations$leafarea + (excavations$leafarea*0.31)
excavations$rsarea <- excavations$rootarea/(excavations$stemarea+excavations$leafarea)




#juniper#############################################################################
jleafarea <- subset(excavations, Species == "Juniper")
jleafarea$shootarea <- jleafarea$leafarea + jleafarea$stemarea
jleafarea$shootmass <- jleafarea$leafmass + jleafarea$stemmass
jleafarea$sla <- jleafarea$leafarea/jleafarea$leafmass

jleafarea150 <- subset(jleafarea, hclass == "<150")
jleafarea150300 <- subset(jleafarea, hclass == "150-300")
jleafarea300 <- subset(jleafarea, hclass == "300+")


#testing correlation between height and diameter and leaf and shoot area/weight
#la v height
#LA = 0.23*h, R2=0.84, p<0.0001 RMSE=10.12 cm2
y <- jleafarea150$leafarea
x <- jleafarea150$height
z <- lm(y ~ 0 + x)
summary(z)
sqrt(mean(z$residuals^2))


#LA = 1.892*h - 144.73, R2=0.26 p=0.11 RMSE=214.91 cm2
x <- jleafarea300$height
y <- jleafarea300$leafarea
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))


#LW v h
#LW = 0.0051*h, R2=0.87 p=1.83e-05 RMSE=0.202 g
y <- jleafarea150$leafmass
x <- jleafarea150$height
z <- lm(y ~ 0 + x)
summary(z)
sqrt(mean(z$residuals^2))


#LW = -4.603 + 0.054*h, R=0.42, p=0.049 RMSE=4.63 g
x <- jleafarea300$height
y <- jleafarea300$leafmass
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))






#la v diameter
#LA = 33.92*d + 2.04, R2=0.25, p=0.02 RMSE=197.93 cm2
y <- jleafarea150$leafarea
x <- jleafarea150$diameter
z <- lm(y ~ 0 + x)
summary(z)
sqrt(mean(z$residuals^2))



#LW v d
#LW = 0.22*d, R2=0.88, p=1.26e-05 RMSE=0.194 g
y <- jleafarea150$leafmass
x <- jleafarea150$diameter
z <- lm(y ~ 0 + x)
summary(z)
sqrt(mean(z$residuals^2))


#LW = 0.282 + 1.152*d, R2=0.45 p=0.04 RMSE=4.51 g
x <- jleafarea300$diameter
y <- jleafarea300$leafmass
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))









alltrees$leafarea <- alltrees$leafarea+(alltrees$leafarea*0.31)
juniper <- subset(alltrees, Species == "Juniper")
juniper$sla <- juniper$leafarea/juniper$leafmass

juniper$hclass <- gsub("<150", "< 150", juniper$hclass)
juniper$hclass <- factor(juniper$hclass, levels = c("< 150", "150-300", "300+"))

juniper150m <- subset(juniper, hclass == "< 150" & rootcat == "Measured")

juniper150300m <- subset(juniper, hclass == "150-300" & rootcat == "Measured")

juniper300m <- subset(juniper, hclass == "300+" & rootcat == "Measured")



xy <- data.frame(x=c(0.9,1.1), 
                 y=c(0.24,2.08))


fig4a <- ggplot() +
  geom_boxplot(juniper, mapping= aes(x=hclass, y=sla), fill="orangered2", outlier.shape = NA, width=0.6) +
  stat_boxplot(juniper, mapping= aes(x=hclass, y=sla), geom = "errorbar", width=0.35) +
  scale_y_continuous(expand = c(0,0), breaks = c(0,50,100,150,200,250), limits = c(0,250)) +
  #scale_y_log10(expand=c(0.0001,0.0001), breaks = c(0.1,1,10,100), limits = c(0.1,100))+
  ylab(bquote("Specific leaf area (" *cm^2~g^-1*")")) +
  labs(title = "Juniper") +
  theme_classic() +
  theme(legend.title = element_blank(), legend.position = c(0.34,0.77), legend.background = 
          element_rect(linetype = "solid", color = "black"), legend.text = element_text(size = 14),
        legend.key.width = unit(2.5, "line"), 
        axis.text.x = element_text(size = 26), axis.title.x = element_blank(),
        plot.title = element_text(size = 28), axis.text = element_text(size=26), axis.title = element_text(size=26)) +
  #geom_rect(data=xy, aes(xmin=0.5, xmax=1.5, ymin=0.24, ymax=2.08), fill="gray80", color="gray80", size=1) +
  #geom_rect(data=xy, aes(xmin=1.5, xmax=2.5, ymin=1.8, ymax=44.75), fill="gray80", color="gray80", size=1) +
  #geom_rect(data=xy, aes(xmin=2.5, xmax=3.5, ymin=1.8, ymax=44.75), fill="gray80", color="gray80", size=1) +
  #geom_boxplot(juniper, mapping= aes(x=hclass, y=sla), fill="orangered2", outlier.shape = NA) +
  #stat_boxplot(juniper, mapping= aes(x=hclass, y=sla), geom = "errorbar", width=0.4) +
  annotate(geom = "text", x=1, y=80, label="b", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=2, y=60, label="n=6", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=3, y=220, label="n=2", fontface='italic', color="black", size=10)

t.test(juniper150300m$sla, juniper150m$sla, alternative = "less") 





fig4aa <- ggplot() +
  geom_boxplot(juniper, mapping= aes(x=hclass, y=leafmass), fill="orangered2", outlier.shape = NA, width=0.6) +
  stat_boxplot(juniper, mapping= aes(x=hclass, y=leafmass), geom = "errorbar", width=0.35) +
  #scale_y_continuous(expand = c(0,0), breaks = c(0,50,100,150,200), limits = c(0,200)) +
  scale_y_log10(expand=c(0.0001,0.0001), breaks = c(0.1,1,10,100), limits = c(0.1,100)) +
  ylab("Leaf mass (g)") +
  labs(title = "Juniper") +
  theme_classic() +
  theme(legend.title = element_blank(), legend.position = c(0.34,0.77), legend.background = 
          element_rect(linetype = "solid", color = "black"), legend.text = element_text(size = 14),
        legend.key.width = unit(2.5, "line"),
        axis.text.x = element_text(size = 26), axis.title.x = element_blank(),
        plot.title = element_text(size = 28), axis.text = element_text(size=26), axis.title = element_text(size=26)) +
  geom_rect(data=xy, aes(xmin=0.5, xmax=1.5, ymin=0.24, ymax=2.08), fill="gray80", color="gray80", size=1) +
  geom_rect(data=xy, aes(xmin=1.5, xmax=2.5, ymin=1.8, ymax=44.75), fill="gray80", color="gray80", size=1) +
  geom_rect(data=xy, aes(xmin=2.5, xmax=3.5, ymin=1.8, ymax=44.75), fill="gray80", color="gray80", size=1) +
  geom_boxplot(juniper, mapping= aes(x=hclass, y=leafmass), fill="orangered2", outlier.shape = NA, width=0.6) +
  stat_boxplot(juniper, mapping= aes(x=hclass, y=leafmass), geom = "errorbar", width=0.35) +
  annotate(geom = "text", x=1, y=2.5, label="ab", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=2, y=40, label="n=6", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=3, y=50, label="n=2", fontface='italic', color="black", size=10)

#excavations only
t.test(juniper150300m$leafmass, juniper150m$leafmass, alternative = "greater") 
t.test(juniper300m$leafmass, juniper150300m$leafmass, alternative = "greater") 
t.test(juniper300m$leafmass, juniper150m$leafmass, alternative = "greater") 



fig4bb <- ggplot() +
  geom_boxplot(juniper, mapping= aes(x=hclass, y=leafarea), fill="orangered2", outlier.shape = NA, width=0.6) +
  stat_boxplot(juniper, mapping= aes(x=hclass, y=leafarea), geom = "errorbar", width=0.35) +
  scale_y_log10(expand = c(0.0001,0.0001), breaks = c(1,10,100,1000,10000), limits = c(1,10000)) +
  labs(y=expression(paste("Leaf area (cm"^{2},")")), title = "Juniper") +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_blank(), legend.position = "none",
        legend.background = element_blank(),
        axis.text.x = element_text(size = 26), axis.title.x = element_blank(), axis.title.y = element_text(size = 26),
        legend.spacing.y = unit(1, "mm"),
        plot.title = element_text(size = 28), axis.text = element_text(size=26)) +
  geom_rect(data=xy, aes(xmin=0.5, xmax=1.5, ymin=13.8, ymax=323.6), color="gray80", fill="gray80", size=1) +
  geom_rect(data=xy, aes(xmin=1.5, xmax=2.5, ymin=46.48, ymax=1154), color="gray80", fill="gray80", size=1) +
  geom_rect(data=xy, aes(xmin=2.5, xmax=3.5, ymin=46.48, ymax=1154), color="gray80", fill="gray80", size=1) +
  geom_boxplot(juniper, mapping= aes(x=hclass, y=leafarea), fill="orangered2", outlier.shape = NA, width=0.6) +
  stat_boxplot(juniper, mapping= aes(x=hclass, y=leafarea), geom = "errorbar", width=0.35) +
  annotate(geom = "text", x=1, y=200, label="a", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=2, y=1900, label="n=6", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=3, y=2300, label="n=2", fontface='italic', color="black", size=10)

# #excavations only
t.test(juniper150300m$leafarea, juniper150m$leafarea, alternative = "greater") #p=0.04
t.test(juniper300m$leafarea, juniper150300m$leafarea, alternative = "greater") #not sig
t.test(juniper300m$leafarea, juniper150m$leafarea, alternative = "greater") #not sig








#pinon###############################################################################
pinleafarea <- subset(excavations, Species == "Piñon")
pinleafarea$sla <- pinleafarea$leafarea/pinleafarea$leafmass

pinleafarea$shootarea <- pinleafarea$leafarea + pinleafarea$stemarea
pinleafarea$shootmass <- pinleafarea$leafmass + pinleafarea$stemmass

pinleafarea150 <- subset(pinleafarea, hclass == "<150")
pinleafarea150300 <- subset(pinleafarea, hclass == "150-300")
pinleafarea300 <- subset(pinleafarea, hclass == "300+")

#testing correlation between height and diameter and leaf area/weight
#sLA v h
#sLA = -0.225*h + 91.55, R2=0.39 p=2.5e-05 RMSE=8.85 cm2/g
y <- pinleafarea150$sla
x <- pinleafarea150$height
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))


#sLA = -0.065*h + 64.18, R2=0.34 p=0.005 RMSE=6.57 cm^2/g
y <- pinleafarea300$sla
x <- pinleafarea300$height
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))





#LA v h
#LA = 0.638*h, R2=0.73 p<0.0001 RMSE=39.46 cm2/g
y <- pinleafarea150$leafarea
x <- pinleafarea150$height
z <- lm(y ~ 0 + x)
summary(z)
sqrt(mean(z$residuals^2))


#LA = 1.746*h, R2=0.67 p<0.0001 RMSE=293.158 cm^2/g
y <- pinleafarea300$leafarea
x <- pinleafarea300$height
z <- lm(y ~ 0 + x)
summary(z)
sqrt(mean(z$residuals^2))





#LW v h
#LW = 0.0079*h, R2=0.71 p = 7.1e-06 RMSE=0.51 g
y <- pinleafarea150$leafmass
x <- pinleafarea150$height
z <- lm(y ~ 0 + x)
summary(z)
sqrt(mean(z$residuals^2))


#LW = 0.0336*h, R2=0.57 p=7.5e-05 RMSE=7.02 g
y <- pinleafarea300$leafmass
x <- pinleafarea300$height
z <- lm(y ~ 0 + x)
summary(z)
sqrt(mean(z$residuals^2))






#SLA v d
#SLA = -5.215*d + 89.18, R2=0.3 p=0.0003 RMSE=9.52 cm^2/g
y <- pinleafarea150$sla
x <- pinleafarea150$diameter
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))


#SLA = 5.468*d, R2=0.76 p=2.9e-07 RMSE=23.41 cm^2/g
y <- pinleafarea300$sla
x <- pinleafarea300$diameter
z <- lm(y ~ 0 + x)
summary(z)
sqrt(mean(z$residuals^2))




#LA v d
#LA = 17.15*d, R2=0.78 p<0.0001 RMSE=35.452 cm^2
y <- pinleafarea150$leafarea
x <- pinleafarea150$diameter
z <- lm(y ~ 0 + x)
summary(z)
sqrt(mean(z$residuals^2))


#LA = 56.503*d, R2=0.71 p<0.0001 RMSE=275.26 cm^2
y <- pinleafarea300$leafarea
x <- pinleafarea300$diameter
z <- lm(y ~ 0 + x)
summary(z)
sqrt(mean(z$residuals^2))






#LW v d
#LW = 0.21*d, R2=0.75 p=1.88e-12 RMSE=0.465 g
y <- pinleafarea150$leafmass
x <- pinleafarea150$diameter
z <- lm(y ~ 0 + x)
summary(z)
sqrt(mean(z$residuals^2))

#LW = 1.094*d, R2=0.61, p=2.9e-05 RMSE=6.68 g
y <- pinleafarea300$leafmass
x <- pinleafarea300$diameter
z <- lm(y ~ 0 + x)
summary(z)
sqrt(mean(z$residuals^2))









pinon <- subset(alltrees, Species == "Piñon")
pinon$sla <- pinon$leafarea/pinon$leafmass

pinon$hclass <- gsub("<150", "< 150", pinon$hclass)
pinon$hclass <- factor(pinon$hclass, levels = c("< 150", "150-300", "300+"))


pinon150m <- subset(pinon, hclass == "< 150" & rootcat == "Measured")

pinon150300m <- subset(pinon, hclass == "150-300" & rootcat == "Measured")

pinon300m <- subset(pinon, hclass == "300+" & rootcat == "Measured")







fig4b <- ggplot() +
  geom_boxplot(pinon, mapping= aes(x=hclass, y=sla),  fill="orange3", outlier.shape = NA, width=0.6) + 
  stat_boxplot(pinon, mapping= aes(x=hclass, y=sla), geom = "errorbar", width=0.35) +
  #scale_y_log10(expand=c(0.0001,0.0001), limits = c(0.1,100), breaks = c(0.01,0.1,1,10,100)) +
  scale_y_continuous(expand = c(0,0), breaks = c(0,50,100,150,200), limits = c(0,200)) +
  labs(title = "Piñon pine") +
  theme_classic() +
  theme(legend.title = element_blank(), legend.position = "none", legend.text = element_blank(),
        axis.text.x = element_text(size = 26), axis.title.x = element_blank(), 
        legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 28), axis.text = element_text(size=26), 
        axis.title = element_text(size=22), axis.title.y.left = element_blank()) +
  #geom_rect(data=xy, aes(xmin=0.5, xmax=1.5, ymin=13.93, ymax=87.72), color="gray80", fill="gray80", size=1) +
  #geom_rect(data=xy, aes(xmin=1.5, xmax=2.5, ymin=13.83, ymax=205.49), color="gray80", fill="gray80", size=1) +
  #geom_rect(data=xy, aes(xmin=2.5, xmax=3.5, ymin=13.83, ymax=205.49), color="gray80", fill="gray80", size=1) +
  #geom_boxplot(pinon, mapping= aes(x=hclass, y=sla),  fill="orange3", outlier.shape = NA) + 
  #stat_boxplot(pinon, mapping= aes(x=hclass, y=sla), geom = "errorbar", width=0.4) +
  annotate(geom = "text", x=1, y=155, label="a", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=2, y=107, label="df", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=3, y=85, label="n=4", fontface='italic', color="black", size=10)

t.test(pinon150300m$sla, pinon150m$sla, alternative = "less") #p<0.05
t.test(pinon300m$sla, pinon150300m$sla, alternative = "greater") 
t.test(pinon300m$sla, pinon150m$sla, alternative = "greater") 




fig4cc <- ggplot() +
  geom_boxplot(pinon, mapping= aes(x=hclass, y=leafmass),  fill="orange3", outlier.shape = NA, width=0.6) + 
  stat_boxplot(pinon, mapping= aes(x=hclass, y=leafmass), geom = "errorbar", width=0.35) +
  scale_y_log10(expand=c(0.0001,0.0001), limits = c(0.1,100), breaks = c(0.01,0.1,1,10,100)) +
  #scale_y_continuous(expand = c(0,0), breaks = c(0,50,100,150,200), limits = c(0,200)) +
  ylab("Leaf mass (g)") +
  labs(title = "Piñon pine") +
  theme_classic() +
  theme(legend.title = element_blank(), legend.position = "none", legend.text = element_blank(),
        axis.text.x = element_text(size = 26), axis.title.x = element_blank(), 
        legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 28), axis.text = element_text(size=26), 
        axis.title = element_text(size=26)) +
  geom_rect(data=xy, aes(xmin=0.5, xmax=1.5, ymin=0.13, ymax=3.03), color="gray80", fill="gray80", size=1) +
  geom_rect(data=xy, aes(xmin=1.5, xmax=2.5, ymin=2.77, ymax=41.11), color="gray80", fill="gray80", size=1) +
  geom_rect(data=xy, aes(xmin=2.5, xmax=3.5, ymin=2.77, ymax=41.11), color="gray80", fill="gray80", size=1) +
  geom_boxplot(pinon, mapping= aes(x=hclass, y=leafmass),  fill="orange3", outlier.shape = NA, width=0.6) + 
  stat_boxplot(pinon, mapping= aes(x=hclass, y=leafmass), geom = "errorbar", width=0.35) +
  annotate(geom = "text", x=1, y=5, label="ac", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=2, y=17, label="d", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=3, y=45, label="n=4", fontface='italic', color="black", size=10)

# t.test(pinon150300m$leafmass, pinon150m$leafmass, alternative = "greater")
# t.test(pinon300m$leafmass, pinon150300m$leafmass, alternative = "greater")
# t.test(pinon300m$leafmass, pinon150m$leafmass, alternative = "greater") 




fig4dd <- ggplot() +
  geom_boxplot(pinon, mapping= aes(x=hclass, y=leafarea), fill="orange3", outlier.shape = NA, width=0.6) +
  stat_boxplot(pinon, mapping= aes(x=hclass, y=leafarea), geom = "errorbar", width=0.35) +
  scale_y_log10(expand = c(0.0001,0.0001), breaks = c(1,10,100,1000,10000), limits = c(1,10000)) +
  labs(y=expression(paste("Leaf area (cm"^{2},")")), title = "Piñon pine") +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_blank(), legend.position = "none",
        legend.background = element_blank(), axis.text.x = element_text(size = 26), axis.title.x = element_blank(),
        axis.title.y = element_text(size = 26), legend.spacing.y = unit(1, "mm"),
        plot.title = element_text(size = 28), axis.text = element_text(size=26)) +
  geom_rect(data=xy, aes(xmin=0.5, xmax=1.5, ymin=10.85, ymax=247.47), color="gray80", fill="gray80", size=1) +
  geom_rect(data=xy, aes(xmin=1.5, xmax=2.5, ymin=142.95, ymax=2123.38), color="gray80", fill="gray80", size=1) +
  geom_rect(data=xy, aes(xmin=2.5, xmax=3.5, ymin=142.95, ymax=2123.38), color="gray80", fill="gray80", size=1) +
  geom_boxplot(pinon, mapping= aes(x=hclass, y=leafarea), fill="orange3", outlier.shape = NA, width=0.6) +
  stat_boxplot(pinon, mapping= aes(x=hclass, y=leafarea), geom = "errorbar", width=0.35) +
  annotate(geom = "text", x=1, y=600, label="b", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=2, y=1400, label="c", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=3, y=4600, label="n=4", fontface='italic', color="black", size=10)


# t.test(pinon150300m$leafarea, pinon150m$leafarea, alternative = "greater") #p=0.0005
# t.test(pinon300m$leafarea, pinon150300m$leafarea, alternative = "greater") #p=0.03
# t.test(pinon300m$leafarea, pinon150m$leafarea, alternative = "greater") #p=0.02






#ponderosa##############################################################################
ponleafarea <- subset(excavations, Species == "Ponderosa")
ponleafarea$shootarea <- ponleafarea$leafarea + ponleafarea$stemarea
ponleafarea$shootmass <- ponleafarea$leafmass + ponleafarea$stemmass
ponleafarea$sla <- ponleafarea$leafarea/ponleafarea$leafmass

ponleafarea150 <- subset(ponleafarea, hclass == "<150")
ponleafarea150300 <- subset(ponleafarea, hclass == "150-300")
ponleafarea300 <- subset(ponleafarea, hclass == "300+")

#testing correlation between height and diameter and leaf area/mass
#SLA v h
#SLA = -0.148*h + 73.6, R2=0.13 p=0.001 RMSE=11.28 cm^2/g
y <- ponleafarea150$sla
x <- ponleafarea150$height
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))


#SLA = -0.074*h + 68.54, R2=0.10 p=0.06 RMSE=6.84 cm^2/g
y <- ponleafarea150300$sla
x <- ponleafarea150300$height
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))


#SLA = -0.029*h + 60.2, R2=0.08 p=0.21 RMSE=7.63 cm^2/g
y <- ponleafarea300$sla
x <- ponleafarea300$height
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))




#LA v h
#LA = 0.49*h, R2=0.74 p<0.0001 RMSE=23.13 cm^2
y <- ponleafarea150$leafarea
x <- ponleafarea150$height
z <- lm(y ~ 0 + x)
summary(z)
sqrt(mean(z$residuals^2))


#LA = 1.48*h - 141.094, R2=0.35 p=0.001 RMSE=71.85 cm^2
y <- ponleafarea150300$leafarea
x <- ponleafarea150300$height
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))


#LA = 4*h - 1129.744, R2=0.56 p=0.008 RMSE=393.53 cm^2
y <- ponleafarea300$leafarea
x <- ponleafarea300$height
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))







#LW v h
#LW = 0.0068*h, R2=0.63 p<2.2e-16 RMSE=0.41 g
y <- ponleafarea150$leafmass
x <- ponleafarea150$height
z <- lm(y ~ 0 + x)
summary(z)
sqrt(mean(z$residuals^2))


#LW = -2.781 + 0.0256*h, R2=0.31 p=0.002 RMSE=1.36 g 
y <- ponleafarea150300$leafmass
x <- ponleafarea150300$height
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))


#LW = -22.403 + 0.076*h, R2=0.54 p=0.009 RMSE=7.85
y <- ponleafarea300$leafmass
x <- ponleafarea300$height
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))








#SLA v d
#SLA = -4.232*d + 75.66, R2=0.16 p=0.0002 RMSE=11.06 cm^2/g
y <- ponleafarea150$sla
x <- ponleafarea150$diameter
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))


#SLA = -2.807*d + 70.5, R2=0.27 p=0.005 RMSE=6.2 cm^2/g
y <- ponleafarea150300$sla
x <- ponleafarea150300$diameter
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))


#sLA = 3.013*d, R2=0.69 p=0.001 RMSE=25.16 cm^2
y <- ponleafarea300$sla
x <- ponleafarea300$diameter
z <- lm(y ~ 0 + x)
summary(z)
sqrt(mean(z$residuals^2))





#LA v d
#LA = 12.203*d, R2=0.78 p<0.0001 RMSE=21.09 cm^2
y <- ponleafarea150$leafarea
x <- ponleafarea150$diameter
z <- lm(y ~ 0 + x)
summary(z)
sqrt(mean(z$residuals^2))


#LA = 28.338*d, R2=0.87 p<0.0001 RMSE=67.53 cm^2
y <- ponleafarea150300$leafarea
x <- ponleafarea150300$diameter
z <- lm(y ~ 0 + x)
summary(z)
sqrt(mean(z$residuals^2))


#LA = 61.54*d, R2=0.78 p=0.0002 RMSE=409.18 cm^2
y <- ponleafarea300$leafarea
x <- ponleafarea300$diameter
z <- lm(y ~ 0 + x)
summary(z)
sqrt(mean(z$residuals^2))






#LW v diameter
#LW = 0.173*d, R2=0.69 p<2.2e-16 RMSE=0.373g
y <- ponleafarea150$leafmass
x <- ponleafarea150$diameter
z <- lm(y ~ 0 + x)
summary(z)
sqrt(mean(z$residuals^2))


#LW = 0.443*d, R2=0.83, p=7.6e-11 RMSE=1.264
y <- ponleafarea150300$leafmass
x <- ponleafarea150300$diameter
z <- lm(y ~ 0 + x)
summary(z)
sqrt(mean(z$residuals^2))


#LW = 1.1237*d, R2=0.76, p=0.0003 RMSE=8.058
y <- ponleafarea300$leafmass
x <- ponleafarea300$diameter
z <- lm(y ~ 0 + x)
summary(z)
sqrt(mean(z$residuals^2))






ponderosa <- subset(alltrees, Species == "Ponderosa")
ponderosa$sla <- ponderosa$leafarea/ponderosa$leafmass


ponderosa$hclass <- gsub("<150", "< 150", ponderosa$hclass)
ponderosa$hclass <- factor(ponderosa$hclass, levels = c("< 150", "150-300", "300+"))


ponderosa150m <- subset(ponderosa, hclass == "< 150" & rootcat == "Measured")

ponderosa150300m <- subset(ponderosa, hclass == "150-300" & rootcat == "Measured")

ponderosa300m <- subset(ponderosa, hclass == "300+" & rootcat == "Measured")



fig4c <- ggplot() +
  geom_boxplot(ponderosa, mapping= aes(x=hclass, y=sla), fill="orchid3", outlier.shape = NA, width=0.6) +
  stat_boxplot(ponderosa, mapping= aes(x=hclass, y=sla), geom = "errorbar", width=0.35) +
  scale_y_continuous(expand = c(0,0), breaks = c(0,50,100,150), limits = c(0,150)) +
  #scale_y_log10(limits = c(0.1,107), expand=c(0,0)) +
  ylab(bquote("Specific leaf area ("*cm^2~g^-1*")")) +
  labs(x="Height class (mm)", title = "Ponderosa pine") +
  theme_classic() +
  theme(legend.title = element_blank(), legend.position = "none", 
        axis.text.x = element_text(size = 26), 
        legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 28), axis.text = element_text(size=26), 
        axis.title = element_text(size=26)) +
  #geom_rect(data=xy, aes(xmin=0.5, xmax=1.5, ymin=37.15, ymax=71.51), color="gray80", fill="gray80", size=1) +
  #geom_rect(data=xy, aes(xmin=1.5, xmax=2.5, ymin=21.18, ymax=68.2), color="gray80", fill="gray80", size=1) +
  #geom_rect(data=xy, aes(xmin=2.5, xmax=3.5, ymin=19.95, ymax=210.31), color="gray80", fill="gray80", size=1) +
  #geom_boxplot(ponderosa, mapping= aes(x=hclass, y=sla), fill="orchid3", outlier.shape = NA) +
  #stat_boxplot(ponderosa, mapping= aes(x=hclass, y=sla), geom = "errorbar", width=0.4) +
  annotate(geom = "text", x=1, y=137, label="ce", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=2, y=115, label="d", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=3, y=115, label="f", fontface='italic', color="black", size=10)

# t.test(ponderosa150300m$sla, ponderosa150m$sla, alternative = "less") #p<0.05
# t.test(ponderosa300m$sla, ponderosa150300m$sla, alternative = "less") #p<0.05
# t.test(ponderosa300m$sla, ponderosa150m$sla, alternative = "less") #p<0.05


fig4ee <- ggplot() +
  geom_boxplot(ponderosa, mapping= aes(x=hclass, y=leafmass), fill="orchid3", outlier.shape = NA, width=0.6) +
  stat_boxplot(ponderosa, mapping= aes(x=hclass, y=leafmass), geom = "errorbar", width=0.35) +
  #scale_y_continuous(expand = c(0,0), breaks = c(0,50,100,150,200), limits = c(0,200)) +
  scale_y_log10(limits = c(0.1,115), expand=c(0,0)) +
  ylab("Leaf mass (g)") +
  labs(title = "Ponderosa pine") +
  theme_classic() +
  theme(legend.title = element_blank(), legend.position = "none", 
        axis.text.x = element_text(size = 26), 
        legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 28), axis.text = element_text(size=26), 
        axis.title = element_text(size=26), axis.title.x = element_blank()) +
  geom_rect(data=xy, aes(xmin=0.5, xmax=1.5, ymin=0.16, ymax=1.57), color="gray80", fill="gray80", size=1) +
  geom_rect(data=xy, aes(xmin=1.5, xmax=2.5, ymin=0.36, ymax=7.78), color="gray80", fill="gray80", size=1) +
  geom_rect(data=xy, aes(xmin=2.5, xmax=3.5, ymin=0.47, ymax=106.42), color="gray80", fill="gray80", size=1) +
  geom_boxplot(ponderosa, mapping= aes(x=hclass, y=leafmass), fill="orchid3", outlier.shape = NA, width=0.6) +
  stat_boxplot(ponderosa, mapping= aes(x=hclass, y=leafmass), geom = "errorbar", width=0.35) +
  annotate(geom = "text", x=1, y=8, label="ab", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=2, y=20, label="d", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=3, y=90, label="e", fontface='italic', color="black", size=10)

# t.test(ponderosa150300m$leafmass, ponderosa150m$leafmass, alternative = "greater") 
# t.test(ponderosa300m$leafmass, ponderosa150300m$leafmass, alternative = "greater") 
# t.test(ponderosa300m$leafmass, ponderosa150m$leafmass, alternative = "greater") 


fig4ff <- ggplot() +
  geom_boxplot(ponderosa, mapping= aes(x=hclass, y=leafarea), fill="orchid3", outlier.shape = NA, width=0.6) +
  stat_boxplot(ponderosa, mapping= aes(x=hclass, y=leafarea), geom = "errorbar", width=0.35) +
  scale_y_log10(expand = c(0.0001,0.0001), breaks = c(1,10,100,1000,10000), limits = c(1,12000)) +
  labs(y=expression(paste("Leaf area (cm"^{2},")")), title = "Ponderosa pine") +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_blank(), legend.position = "none",
        legend.background = element_blank(), axis.text.x = element_text(size = 26), axis.title.x = element_blank(),
        axis.title.y = element_text(size = 26), legend.spacing.y = unit(1, "mm"),
        plot.title = element_text(size = 28), axis.text = element_text(size=26)) +
  geom_rect(data=xy, aes(xmin=0.5, xmax=1.5, ymin=11.76, ymax=111.05), color="gray80", fill="gray80", size=1) +
  geom_rect(data=xy, aes(xmin=1.5, xmax=2.5, ymin=23.24, ymax=497.9), color="gray80", fill="gray80", size=1) +
  geom_rect(data=xy, aes(xmin=2.5, xmax=3.5, ymin=74.26, ymax=5650.26), color="gray80", fill="gray80", size=1) +
  geom_boxplot(ponderosa, mapping= aes(x=hclass, y=leafarea), fill="orchid3", outlier.shape = NA, width=0.6) +
  stat_boxplot(ponderosa, mapping= aes(x=hclass, y=leafarea), geom = "errorbar", width=0.35) +
  annotate(geom = "text", x=1, y=500, label="a", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=2, y=1400, label="c", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=3, y=6000, label="d", fontface='italic', color="black", size=10)

# t.test(ponderosa150300m$leafarea, ponderosa150m$leafarea, alternative = "greater") #p=5.3e-08
# t.test(ponderosa300m$leafarea, ponderosa150300m$leafarea, alternative = "greater") #p=0.02
# t.test(ponderosa300m$leafarea, ponderosa150m$leafarea, alternative = "greater") #p=0.007







#fir###################################################################################
fleafarea <- subset(excavations, Species == "Fir")
fleafarea$shootarea <- fleafarea$leafarea + fleafarea$stemarea
fleafarea$shootmass <- fleafarea$leafmass + fleafarea$stemmass
fleafarea$sla <- fleafarea$leafarea/fleafarea$leafmass

fleafarea150 <- subset(fleafarea, hclass == "<150")
fleafarea150300 <- subset(fleafarea, hclass != "<150")

#testing correlation between height and diameter and leaf area/mass
#SLA v h
#SLA = 0.007*h + 68.51, R2=0 p=0.93 RMSE=12.23 cm2
y <- fleafarea150$sla
x <- fleafarea150$height
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))


#SLA = -0.016*h + 70.55, R2=0 p=0.92 RMSE=16.4 cm2
y <- fleafarea150300$sla
x <- fleafarea150300$height
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))

 


#LA v h
#LA = 0.392*h, R2=0.78 p<0.0001 RMSE=17.4 cm2
y <- fleafarea150$leafarea
x <- fleafarea150$height
z <- lm(y ~ 0 + x)
summary(z)
sqrt(mean(z$residuals^2))


#LA = 3.637*h - 526.747, R2=0.63 p=0.004 RMSE=93.06 cm2
y <- fleafarea150300$leafarea
x <- fleafarea150300$height
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))






#LW v h
#LW = 0.0046*h, R2=0.72 p=9.7e-10 RMSE=0.24 g
y <- fleafarea150$leafmass
x <- fleafarea150$height
z <- lm(y ~ 0 + x)
summary(z)
sqrt(mean(z$residuals^2))


#LW = -8.21 + 0.055*h, R2=0.45 p=0.02 RMSE=1.91 g
y <- fleafarea150300$leafmass
x <- fleafarea150300$height
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))







#SLA v d
#SLA = -2.21*d + 75.3, R2=0 p=0.31 RMSE=12.01 cm2
y <- fleafarea150$sla
x <- fleafarea150$diameter
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))


#SLA = -6.02*d + 98.5, R2=0.35 p=0.04 RMSE=12.47 cm2
y <- fleafarea150300$sla
x <- fleafarea150300$diameter
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))




#LA v d
#LA = 11.086*d, R2=0.76 p<0.0001 RMSE=18.32 cm2
y <- fleafarea150$leafarea
x <- fleafarea150$diameter
z <- lm(y ~ 0 + x)
summary(z)
sqrt(mean(z$residuals^2))


#LA = 37.16*d, R2=0.75 p=0.0004 RMSE=109.9 cm2
y <- fleafarea150300$leafarea
x <- fleafarea150300$diameter
z <- lm(y ~ 0 + x)
summary(z)
sqrt(mean(z$residuals^2))






#LW v d
#LW = 0.131*d, R2=0.71 p<2.1e-09 RMSE=0.257 g
y <- fleafarea150$leafmass
x <- fleafarea150$diameter
z <- lm(y ~ 0 + x)
summary(z)
sqrt(mean(z$residuals^2))


#LW = 0.53*d, R2=0.65, p=0.002 RMSE=1.98 g
y <- fleafarea150300$leafmass
x <- fleafarea150300$diameter
z <- lm(y ~ 0 + x)
summary(z)
sqrt(mean(z$residuals^2))








fir <- subset(alltrees, Species == "Fir")
fir$sla <- fir$leafarea/fir$leafmass

fir$hclass <- gsub("<150", "< 150", fir$hclass)
fir$hclass <- factor(fir$hclass, levels = c("< 150", "150-300", "300+"))
fir <- fir[-c(355), ]


fir150m <- subset(fir, hclass == "< 150" & rootcat == "Measured")

fir150300m <- subset(fir, hclass == "150-300" & rootcat == "Measured")





fig4d <- ggplot() +
  geom_boxplot(fir, mapping= aes(x=hclass, y=sla), fill="palegreen4", outlier.shape = NA, width=0.6) +
  stat_boxplot(fir, mapping= aes(x=hclass, y=sla), geom = "errorbar", width=0.35) +
  scale_y_continuous(expand = c(0,0), breaks = c(0,50,100,150,200), limits = c(0,200)) +
  #scale_y_log10(expand = c(0.0001,0.0001), breaks = c(0.1,1,10,100), limits = c(0.07,100)) +
  labs(x="Height class (mm)", title = "Fir") +
  theme_classic() +
  theme(legend.title = element_blank(), legend.position = "none",
        axis.text.x = element_text(size = 26), 
        axis.title.x = element_text(size = 26), 
        legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 28), axis.text = element_text(size=26), 
        axis.title = element_text(size=26), axis.title.y.left = element_blank()) +
  #geom_rect(data=xy, aes(xmin=0.5, xmax=1.5, ymin=0.11, ymax=1.3), color="gray80", fill="gray80", size=1) +
  #geom_rect(data=xy, aes(xmin=1.5, xmax=2.5, ymin=12.9, ymax=89.65), color="gray80", fill="gray80", size=1) +
  #geom_boxplot(fir, mapping= aes(x=hclass, y=sla), fill="palegreen4", outlier.shape = NA) +
  #stat_boxplot(fir, mapping= aes(x=hclass, y=sla), geom = "errorbar", width=0.4) +
  annotate(geom = "text", x=1, y=150, label="a", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=2, y=150, label="ae", fontface='italic', color="black", size=10)

# t.test(fir150300m$sla, fir150m$sla, alternative = "less") #not sig



fig4gg <- ggplot() +
  geom_boxplot(fir, mapping= aes(x=hclass, y=leafmass), fill="palegreen4", outlier.shape = NA, width=0.6) +
  stat_boxplot(fir, mapping= aes(x=hclass, y=leafmass), geom = "errorbar", width=0.35) +
  #scale_y_continuous(expand = c(0,0), breaks = c(0,50,100,150,200), limits = c(0,200)) +
  scale_y_log10(expand = c(0.0001,0.0001), breaks = c(0.1,1,10,100), limits = c(0.07,100)) +
  ylab("Leaf mass (g)") +
  labs(x="Height class (mm)", title = "Fir") +
  theme_classic() +
  theme(legend.title = element_blank(), legend.position = "none",
        axis.text.x = element_text(size = 26), 
        axis.title.x = element_text(size = 26), 
        legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 28), axis.text = element_text(size=26), 
        axis.title = element_text(size=26)) +
  geom_rect(data=xy, aes(xmin=0.5, xmax=1.5, ymin=0.11, ymax=1.3), color="gray80", fill="gray80", size=1) +
  geom_rect(data=xy, aes(xmin=1.5, xmax=2.5, ymin=0.42, ymax=8.29), color="gray80", fill="gray80", size=1) +
  geom_boxplot(fir, mapping= aes(x=hclass, y=leafmass), fill="palegreen4", outlier.shape = NA, width=0.6) +
  stat_boxplot(fir, mapping= aes(x=hclass, y=leafmass), geom = "errorbar", width=0.35) +
  annotate(geom = "text", x=1, y=4, label="b", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=2, y=30, label="cd", fontface='italic', color="black", size=10)


# t.test(fir150300m$leafmass, fir150m$leafmass, alternative = "greater")



fig4hh <- ggplot() +
  geom_boxplot(fir, mapping= aes(x=hclass, y=leafarea),  fill="palegreen4", outlier.shape = NA, width=0.6) +
  stat_boxplot(fir, mapping= aes(x=hclass, y=leafarea), geom = "errorbar", width=0.35) +
  scale_y_log10(expand = c(0.0001,0.0001), breaks = c(1,10,100,1000,10000), limits = c(1,10000)) +
  labs(x="Height class (mm)", y=expression(paste("Leaf area (cm"^{2},")")), title = "Fir") +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_blank(), legend.position = "none",
        legend.background = element_blank(), axis.text.x = element_text(size = 26),
        axis.title.x = element_text(size = 26),
        axis.title.y = element_text(size = 26), legend.spacing.y = unit(1, "mm"),
        plot.title = element_text(size = 28), axis.text = element_text(size=26)) +
  geom_rect(data=xy, aes(xmin=0.5, xmax=1.5, ymin=9.8, ymax=110.42), color="gray80", fill="gray80", size=1) +
  geom_rect(data=xy, aes(xmin=1.5, xmax=2.5, ymin=44.26, ymax=564.35), color="gray80", fill="gray80", size=1) +
  geom_boxplot(fir, mapping= aes(x=hclass, y=leafarea),  fill="palegreen4", outlier.shape = NA, width=0.6) +
  stat_boxplot(fir, mapping= aes(x=hclass, y=leafarea), geom = "errorbar", width=0.35) +
  annotate(geom = "text", x=1, y=350, label="a", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=2, y=2000, label="c", fontface='italic', color="black", size=10)

# t.test(fir150300m$leafarea, fir150m$leafarea, alternative = "greater") #p=0.01




#figure 4
plot_grid(fig4a, fig4b, fig4c, fig4d, nrow = 2,
          labels = c("a", "b", "c", "d"),
          label_x = 0.92, label_y = 0.97, label_size = 40)

#figure s2
plot_grid(fig4aa, fig4bb, fig4cc, fig4dd, fig4ee, fig4ff, fig4gg, fig4hh, nrow = 4,
          labels = c("a", "b", "c", "d", "e", "f", "g", "h"),
          label_x = 0.92, label_y = 0.97, label_size = 40)





#comparison between species 
t.test(juniper150m$sla, pinyon150m$sla, alternative = "less") #p<0.05
t.test(juniper150m$sla, ponderosa150m$sla, alternative = "less") #p<0.05
t.test(juniper150m$sla, fir150m$sla, alternative = "less") #p<0.05
t.test(pinyon150m$sla, ponderosa150m$sla, alternative = "greater") #p<0.05
t.test(pinyon150m$sla, fir150m$sla, alternative = "greater") 
t.test(ponderosa150m$sla, fir150m$sla, alternative = "less") #p<0.05


t.test(pinyon150300m$sla, ponderosa150300m$sla, alternative = "less") 
t.test(pinyon150300m$sla, fir150300m$sla, alternative = "less") #p<0.05
t.test(ponderosa150300m$sla, fir150300m$sla, alternative = "less") #p<0.05


t.test(pinyon150300m$sla, juniper150m$sla, alternative = "greater") #p<0.0001
t.test(pinyon150300m$sla, ponderosa150m$sla, alternative = "less") #p<0.0001
t.test(pinyon150300m$sla, fir150m$sla, alternative = "less") #p<0.0001
t.test(ponderosa150300m$sla, juniper150m$sla, alternative = "greater") #p<0.0001
t.test(ponderosa150300m$sla, pinyon150m$sla, alternative = "less") #p<0.0001
t.test(ponderosa150300m$sla, fir150m$sla, alternative = "less") #p<0.0001
t.test(fir150300m$sla, juniper150m$sla, alternative = "greater") #p<0.0001
t.test(fir150300m$sla, pinyon150m$sla, alternative = "less")
t.test(fir150300m$sla, ponderosa150m$sla, alternative = "greater") 


t.test(ponderosa300m$sla, juniper150m$sla, alternative = "greater") #p<0.05
t.test(ponderosa300m$sla, pinyon150m$sla, alternative = "less") #p<0.05
t.test(ponderosa300m$sla, fir150m$sla, alternative = "less") #p<0.05


t.test(ponderosa300m$sla, pinyon150300m$sla, alternative = "less") 
t.test(ponderosa300m$sla, fir150300m$sla, alternative = "less") #p<0.05






t.test(juniper150m$leafarea, pinyon150m$leafarea, alternative = "less") #p=0.0002, pinon sig greater
t.test(juniper150m$leafarea, ponderosa150m$leafarea, alternative = "less") 
t.test(juniper150m$leafarea, fir150m$leafarea, alternative = "less") 
t.test(pinyon150m$leafarea, ponderosa150m$leafarea, alternative = "greater") #p=0.003
t.test(pinyon150m$leafarea, fir150m$leafarea, alternative = "greater") #p=0.002
t.test(ponderosa150m$leafarea, fir150m$leafarea, alternative = "greater") 


t.test(juniper150300m$leafarea, pinyon150300m$leafarea, alternative = "greater") 
t.test(juniper150300m$leafarea, ponderosa150300m$leafarea, alternative = "greater") 
t.test(juniper150300m$leafarea, fir150300m$leafarea, alternative = "greater") 
t.test(pinyon150300m$leafarea, ponderosa150300m$leafarea, alternative = "greater") 
t.test(pinyon150300m$leafarea, fir150300m$leafarea, alternative = "greater") 
t.test(ponderosa150300m$leafarea, fir150300m$leafarea, alternative = "greater") 


t.test(pinyon150300m$leafarea, juniper150m$leafarea, alternative = "greater") #p<0.05
t.test(pinyon150300m$leafarea, ponderosa150m$leafarea, alternative = "greater")#p<0.05 
t.test(pinyon150300m$leafarea, fir150m$leafarea, alternative = "greater") #p<0.05
t.test(ponderosa150300m$leafarea, juniper150m$leafarea, alternative = "greater") #p<0.05
t.test(ponderosa150300m$leafarea, pinyon150m$leafarea, alternative = "greater") #p<0.05
t.test(ponderosa150300m$leafarea, fir150m$leafarea, alternative = "greater") #p<0.05
t.test(fir150300m$leafarea, juniper150m$leafarea, alternative = "greater") #p<0.05
t.test(fir150300m$leafarea, pinyon150m$leafarea, alternative = "greater") #p<0.05
t.test(fir150300m$leafarea, ponderosa150m$leafarea, alternative = "greater")  #p<0.05


t.test(ponderosa300m$leafarea, juniper150m$leafarea, alternative = "greater") #p<0.05
t.test(ponderosa300m$leafarea, pinyon150m$leafarea, alternative = "greater") #p<0.05
t.test(ponderosa300m$leafarea, fir150m$leafarea, alternative = "greater") #p<0.05


t.test(ponderosa300m$leafarea, pinyon150300m$leafarea, alternative = "greater") #p<0.05
t.test(ponderosa300m$leafarea, fir150300m$leafarea, alternative = "greater") #p<0.05









#figure 5########################################################################################
#juniper#########################################################
jrootarea <- subset(excavations, Species == "Juniper")

jrootarea150 <- subset(jrootarea, hclass == "<150")
jrootarea300 <- subset(jrootarea, hclass != "<150")


#testing correlation between height and diameter and leaf area/mass
#RA v h
#RA = 0.0689*h, R2=0.91 p=3.4e-06 RMSE=2.21 cm2
y <- jrootarea150$rootarea
x <- jrootarea150$height
z <- lm(y ~ 0 + x)
summary(z)
sqrt(mean(z$residuals^2))


#RA = -18.864 + 0.45*h, R2=0.05 p=0.6 RMSE=164.264 cm2
y <- jrootarea300$rootarea
x <- jrootarea300$height
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))





#RW v h
#RW = 0.0074*h, R2=0.56 p=0.005 RMSE=0.65 g
y <- jrootarea150$rootmass
x <- jrootarea150$height
z <- lm(y ~ 0 + x)
summary(z)
sqrt(mean(z$residuals^2))


#RW = 0.0401*h, R2=0.48 p=0.02 RMSE = 9.64 g
y <- jrootarea300$rootmass
x <- jrootarea300$height
z <- lm(y ~ 0 + x)
summary(z)
sqrt(mean(z$residuals^2))







#RA v d
#RA = 1.215 + 2.478*d, R2=0.30 p=0.056 RMSE=2.37 cm2
y <- jrootarea150$rootarea
x <- jrootarea150$diameter
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))


#RA = 45.36 + 6.524*d, R2=0.02 p=0.7 RMSE = 166.64 cm2
y <- jrootarea300$rootarea
x <- jrootarea300$diameter
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))








#RW v d
#RW = 0.272*d, R2=0.38 p=0.02 RMSE=0.766 g
y <- jrootarea150$rootmass
x <- jrootarea150$diameter
z <- lm(y ~ 0 + x)
summary(z)
sqrt(mean(z$residuals^2))


#RW = 1.45*d, R2=0.71, p=0.003 RMSE=7.21 g
y <- jrootarea300$rootmass
x <- jrootarea300$diameter
z <- lm(y ~ 0 + x)
summary(z)
sqrt(mean(z$residuals^2))







options(scipen=10000) #so values don't appear in scientific notation

fig5a <- ggplot() +
  geom_boxplot(juniper, mapping= aes(x=hclass, y=rootmass), fill="orangered2", outlier.shape = NA, width=0.6) +
  stat_boxplot(juniper, mapping= aes(x=hclass, y=rootmass), geom = "errorbar", width=0.35) +
  scale_y_log10(limits = c(0.01,130), breaks = c(0.01,0.1,1,10,100), expand=c(0.0001,0.0001),
                labels = c("0.01","0.1","1.0","10","100")) +
  labs(y="Root mass (g)", title = "Juniper") +
  theme_classic() +
  theme(legend.title = element_blank(), legend.position = c(0.35,0.9), legend.background = 
          element_rect(linetype = "solid", color = "black"), legend.text = element_text(size = 26), 
        legend.key.width = unit(2.5, "line"),
        axis.text.x = element_text(size = 26), axis.title.x = element_blank(), legend.spacing.y = unit(1, "mm"), 
        plot.title = element_text(size = 28), axis.text = element_text(size=26), axis.title = element_text(size=26)) +
  geom_rect(data=xy, aes(xmin=0.5, xmax=1.5, ymin=0.3, ymax=2.58), color="gray80", fill="gray80", size=1) +
  geom_rect(data=xy, aes(xmin=1.5, xmax=2.5, ymin=1.9, ymax=49.24), color="gray80", fill="gray80", size=1) +
  geom_rect(data=xy, aes(xmin=2.5, xmax=3.5, ymin=1.9, ymax=49.24), color="gray80", fill="gray80", size=1) +
  geom_boxplot(juniper, mapping= aes(x=hclass, y=rootmass), fill="orangered2", outlier.shape = NA, width=0.6) +
  stat_boxplot(juniper, mapping= aes(x=hclass, y=rootmass), geom = "errorbar", width=0.35) +
  annotate(geom = "text", x=1, y=7, label="abc", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=2, y=14, label="n=6", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=3, y=75, label="n=2", fontface='italic', color="black", size=10)

t.test(juniper150300m$rootmass, juniper150m$rootmass, alternative = "greater") #p=0.01
t.test(juniper300m$rootmass, juniper150300m$rootmass, alternative = "greater") 
t.test(juniper300m$rootmass, juniper150m$rootmass, alternative = "greater") 



fig5b <- ggplot() +
  geom_boxplot(juniper, mapping= aes(x=hclass, y=rootarea), fill="orangered2", outlier.shape = NA, width=0.6) +
  stat_boxplot(juniper, mapping= aes(x=hclass, y=rootarea), geom = "errorbar", width=0.35) +
  scale_y_log10(expand = c(0.0001,0.0001), limits= c(1,300), breaks = c(1,3,10,30,100,300)) +
  labs(y=expression(paste("Root area (cm"^{2},")")), title = "Juniper") +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_blank(), legend.position = "none",
        legend.background = element_blank(), 
        axis.text.x = element_text(size = 26), axis.title.x = element_blank(), axis.title.y = element_text(size = 26),
        legend.spacing.y = unit(1, "mm"), 
        plot.title = element_text(size = 28), axis.text = element_text(size=26)) +
  geom_rect(data=xy, aes(xmin=0.5, xmax=1.5, ymin=3.96, ymax=24.71), color="gray80", fill="gray80", size=1) +
  geom_boxplot(juniper, mapping= aes(x=hclass, y=rootarea), fill="orangered2", outlier.shape = NA, width=0.6) +
  stat_boxplot(juniper, mapping= aes(x=hclass, y=rootarea), geom = "errorbar", width=0.35) +
  annotate(geom = "text", x=1, y=23, label="ab", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=2, y=55, label="n=6", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=3, y=130, label="n=2", fontface='italic', color="black", size=10)

t.test(juniper150300m$rootarea, juniper150m$rootarea, alternative = "greater") 
t.test(juniper150300m$rootarea, juniper300m$rootarea, alternative = "greater") 
t.test(juniper300m$rootarea, juniper150m$rootarea, alternative = "greater") 








#pinon###################################################################################
pinrootarea <- subset(excavations, Species == "Pinyon")

pinrootarea150 <- subset(pinrootarea, hclass == "<150")
pinrootarea150300 <- subset(pinrootarea, hclass == "150-300")
pinrootarea300 <- subset(pinrootarea, hclass != "<150")

#testing correlation between height and diameter and leaf area/mass
#RA v h
#RA = 0.975 + 0.0843*h, R2=0.34 p=0.0001 RMSE=3.67 cm2
y <- pinrootarea150$rootarea
x <- pinrootarea150$height
z<- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))


#RA = 2.193 + 0.0889*h, R2=0.02 p=0.28 RMSE=8.66 cm2
y <- pinrootarea150300$rootarea
x <- pinrootarea150300$height
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))


#RA = -15.25 + 0.174*h, R2=0.57 p=0.0001 RMSE=11.27 cm2
y <- pinrootarea300$rootarea
x <- pinrootarea300$height
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))







#RW v h
#RW = 0.00846*h, R2=0.79 p=1.2e-13 RMSE=0.438 g
y <- pinrootarea150$rootmass
x <- pinrootarea150$height
z <- lm(y ~ 0 + x)
summary(z)
sqrt(mean(z$residuals^2))


#RW = -1.581 + 0.02*h, R2=0.19 p=0.058 RMSE=1.054 g
y <- pinrootarea150300$rootmass
x <- pinrootarea150300$height
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))


#RW = -7.059 + 0.0472*h, R2=0.72 p=2.4e-06 RMSE=2.22 g
y <- pinrootarea300$rootmass
x <- pinrootarea300$height
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))







#RA v d
#RA = 2.495*d, R2=0.92 p<2.2e-16 RMSE=2.858 cm2
y <- pinrootarea150$rootarea
x <- pinrootarea150$diameter
z <- lm(y ~ 0 + x)
summary(z)
sqrt(mean(z$residuals^2))


#RA = 3.219*d, R2=0.90 p=1.4e-08 RMSE=6.76 cm2
y <- pinrootarea150300$rootarea
x <- pinrootarea150300$diameter
z <- lm(y ~ 0 + x)
summary(z)
sqrt(mean(z$residuals^2))


#RA = 3.762*d, R2=0.89 p=3.1e-10 RMSE=10.28 cm2
y <- pinrootarea300$rootarea
x <- pinrootarea300$diameter
z <- lm(y ~ 0 + x)
summary(z)
sqrt(mean(z$residuals^2))







#RW v d
#RW = 0.227*d, R2=0.84 p=6.4e-16 RMSE=0.378 g
y <- pinrootarea150$rootmass
x <- pinrootarea150$diameter
z <- lm(y ~ 0 + x)
summary(z)
sqrt(mean(z$residuals^2))


#RW = 0.402*d, R2=0.91, p=5.2e-09 RMSE=0.782 g
y <- pinrootarea150300$rootmass
x <- pinrootarea150300$diameter
z <- lm(y ~ 0 + x)
summary(z)
sqrt(mean(z$residuals^2))


#RW = 0.663*d, R2=0.75, p=4.7e-07 RMSE=2.93 g
y <- pinrootarea300$rootmass
x <- pinrootarea300$diameter
z <- lm(y ~ 0 + x)
summary(z)
sqrt(mean(z$residuals^2))





fig5c <- ggplot() +
  geom_boxplot(pinon, mapping= aes(x=hclass, y=rootmass), fill="orange3", outlier.shape = NA, width=0.6) + 
  stat_boxplot(pinon, mapping= aes(x=hclass, y=rootmass), geom = "errorbar", width=0.35) +
  scale_y_log10(limits = c(0.01,100), breaks = c(0.01, 0.1, 1, 10, 100), expand=c(0.0001,0.0001),
                labels = c("0.01","0.1","1.0","10","100")) +
  labs(y="Root mass (g)", title = "Piñon pine") +
  theme_classic() +
  theme(legend.title = element_blank(), legend.position = "none", 
        axis.text.x = element_text(size = 26), axis.title.x = element_blank(), 
        legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 28), axis.text = element_text(size=26), 
        axis.title = element_text(size=26)) +
  geom_rect(data=xy, aes(xmin=0.5, xmax=1.5, ymin=0.14, ymax=3.27), color="gray80", fill="gray80", size=1) +
  geom_rect(data=xy, aes(xmin=1.5, xmax=2.5, ymin=0.021, ymax=24.91), color="gray80", fill="gray80", size=1) +
  geom_rect(data=xy, aes(xmin=2.5, xmax=3.5, ymin=0.021, ymax=24.91), color="gray80", fill="gray80", size=1) +
  geom_boxplot(pinon, mapping= aes(x=hclass, y=rootmass), fill="orange3", outlier.shape = NA, width=0.6) + 
  stat_boxplot(pinon, mapping= aes(x=hclass, y=rootmass), geom = "errorbar", width=0.35) +
  annotate(geom = "text", x=1, y=5.5, label="a", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=2, y=15, label="e", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=3, y=40, label="n=4", fontface='italic', color="black", size=10)

t.test(pinyon150300m$rootmass, pinyon150m$rootmass, alternative = "greater") #p=6.8e-05
t.test(pinyon300m$rootmass, pinyon150300m$rootmass, alternative = "greater") #p=0.04
t.test(pinyon300m$rootmass, pinyon150m$rootmass, alternative = "greater") #p=0.02



fig5d <- ggplot() +
  geom_boxplot(pinon, mapping= aes(x=hclass, y=rootarea), fill="orange3", outlier.shape = NA, width=0.6) +
  stat_boxplot(pinon, mapping= aes(x=hclass, y=rootarea), geom = "errorbar", width=0.35) +
  scale_y_log10(expand = c(0.0001,0.0001), breaks = c(1,3,10,30,100,300), limits = c(1,300)) +
  labs(y=expression(paste("Root area (cm"^{2},")")), title = "Piñon pine") +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_blank(), legend.position = "none",
        legend.background = element_blank(), axis.text.x = element_text(size = 26), axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 26), legend.spacing.y = unit(1, "mm"), 
        plot.title = element_text(size = 28), axis.text = element_text(size=26)) +
  geom_rect(data=xy, aes(xmin=0.5, xmax=1.5, ymin=1.9, ymax=36), color="gray80", fill="gray80", size=1) +
  geom_rect(data=xy, aes(xmin=1.5, xmax=2.5, ymin=9.52, ymax=141.37), color="gray80", fill="gray80", size=1) +
  geom_rect(data=xy, aes(xmin=2.5, xmax=3.5, ymin=9.52, ymax=141.37), color="gray80", fill="gray80", size=1) +
  geom_boxplot(pinon, mapping= aes(x=hclass, y=rootarea), fill="orange3", outlier.shape = NA, width=0.6) +
  stat_boxplot(pinon, mapping= aes(x=hclass, y=rootarea), geom = "errorbar", width=0.35) +
  annotate(geom = "text", x=1, y=40, label="a", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=2, y=70, label="c", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=3, y=125, label="n=4", fontface='italic', color="black", size=10)

t.test(pinyon150300m$rootarea, pinyon150m$rootarea, alternative = "greater") #p=0.0002
t.test(pinyon300m$rootarea, pinyon150300m$rootarea, alternative = "greater") 
t.test(pinyon300m$rootarea, pinyon150m$rootarea, alternative = "greater") #p=0.02




#ponderosa###########################################################################
ponrootarea <- subset(excavations, Species == "Ponderosa")

ponrootarea150 <- subset(ponrootarea, hclass == "<150")
ponrootarea150300 <- subset(ponrootarea, hclass == "150-300")
ponrootarea300 <- subset(ponrootarea, hclass == ">300")

#testing correlation between height and diameter and leaf area/mass
#RA v h
#RA = 0.0807*h, R2=0.78 p<2.2e-16 RMSE=3.38 cm2
y <- ponrootarea150$rootarea
x <- ponrootarea150$height
z <- lm(y ~ 0 + x)
summary(z)
sqrt(mean(z$residuals^2))


#RA = 1.513 + 0.08*h, R2=0.18 p=0.02 RMSE=5.82 cm2
y <- ponrootarea150300$rootarea
x <- ponrootarea150300$height
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))


#RA = -27.903 + 0.153*h, R2=0.65 p=0.003 RMSE=12.7 cm2
y <- ponrootarea300$rootarea
x <- ponrootarea300$height
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))







#RW v h
#RW = 0.005*h, R2=0.59 p=6.1e-16 RMSE=0.327 g
y <- ponrootarea150$rootmass
x <- ponrootarea150$height
z <- lm(y ~ 0 + x)
summary(z)
sqrt(mean(z$residuals^2))


#RW = 0.95 + 0.0033*h, R2=0.009 p=0.64 RMSE=1.284 g
y <- ponrootarea150300$rootmass
x <- ponrootarea150300$height
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))


#RW = -16.172 + 0.0553*h, R2=0.36 p=0.04 RMSE=7.936 g
y <- ponrootarea300$rootmass
x <- ponrootarea300$height
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))






#RA v d
#RA = 2.014*d, R2=0.83 p<2.2e-16 RMSE=2.97 cm2
y <- ponrootarea150$rootarea
x <- ponrootarea150$diameter
z <- lm(y ~ 0 + x)
summary(z)
sqrt(mean(z$residuals^2))


#RA = 0.0546 + 2.937*d, R2=0.39 p=0.0005 RMSE=5.025 cm2
y <- ponrootarea150300$rootarea
x <- ponrootarea150300$diameter
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))


#RA = 3.45*d, R2=0.95 p=1.4e-07 RMSE=9.59 cm2
y <- ponrootarea300$rootarea
x <- ponrootarea300$diameter
z <- lm(y ~ 0 + x)
summary(z)
sqrt(mean(z$residuals^2))






#RW v d
#RW = 0.13*d, R2=0.67 p<2.2e-16 RMSE=0.292 g
y <- ponrootarea150$rootmass
x <- ponrootarea150$diameter
z <- lm(y ~ 0 + x)
summary(z)
sqrt(mean(z$residuals^2))


#RW = 0.275*d, R2=0.69, p=8.2e-08 RMSE=1.14 g
y <- ponrootarea150300$rootmass
x <- ponrootarea150300$diameter
z <- lm(y ~ 0 + x)
summary(z)
sqrt(mean(z$residuals^2))


#RW = 0.867*d, R2=0.69, p=0.0009 RMSE=7.23 g
y <- ponrootarea300$rootmass
x <- ponrootarea300$diameter
z <- lm(y ~ 0 + x)
summary(z)
sqrt(mean(z$residuals^2))






fig5e <- ggplot() +
  geom_boxplot(ponderosa, mapping= aes(x=hclass, y=rootmass), fill="orchid3", outlier.shape = NA, width=0.6) +
  stat_boxplot(ponderosa, mapping= aes(x=hclass, y=rootmass), geom = "errorbar", width=0.35) +
  scale_y_log10(limits = c(0.01,130), breaks = c(0.01, 0.1, 1, 10, 100), expand=c(0.0001,0.0001),
                labels = c("0.01","0.1","1.0","10","100")) +
  labs(y="Root mass (g)", title = "Ponderosa pine") +
  theme_classic() +
  theme(legend.title = element_blank(), legend.position = "none", 
        axis.text.x = element_text(size = 26), axis.title.x = element_blank(), 
        legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 28), axis.text = element_text(size=26), 
        axis.title = element_text(size=26)) +
  geom_rect(data=xy, aes(xmin=0.5, xmax=1.5, ymin=0.12, ymax=1.18), color="gray80", fill="gray80", size=1) +
  geom_rect(data=xy, aes(xmin=1.5, xmax=2.5, ymin=0.22, ymax=4.83), color="gray80", fill="gray80", size=1) +
  geom_rect(data=xy, aes(xmin=2.5, xmax=3.5, ymin=0.47, ymax=77.56), color="gray80", fill="gray80", size=1) +
  geom_boxplot(ponderosa, mapping= aes(x=hclass, y=rootmass), fill="orchid3", outlier.shape = NA, width=0.6) +
  stat_boxplot(ponderosa, mapping= aes(x=hclass, y=rootmass), geom = "errorbar", width=0.35) +
  annotate(geom = "text", x=1, y=6, label="b", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=2, y=10, label="d", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=3, y=70, label="f", fontface='italic', color="black", size=10)

t.test(ponderosa150300m$rootmass, ponderosa150m$rootmass, alternative = "greater") #p=2.3e-05
t.test(ponderosa300m$rootmass, ponderosa150300m$rootmass, alternative = "greater") #p=0.03
t.test(ponderosa300m$rootmass, ponderosa150m$rootmass, alternative = "greater") #p=0.02



fig5f <- ggplot() +
  geom_boxplot(ponderosa, mapping= aes(x=hclass, y=rootarea), fill="orchid3", outlier.shape = NA, width=0.6) +
  stat_boxplot(ponderosa, mapping= aes(x=hclass, y=rootarea), geom = "errorbar", width=0.35) +
  scale_y_log10(expand = c(0.0001,0.0001), breaks = c(1,3,10,30,100,300), limits = c(1,300)) +
  labs(y=expression(paste("Root area (cm"^{2},")")), title = "Ponderosa pine") +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_blank(), legend.position = "none",
        legend.background = element_blank(), axis.text.x = element_text(size = 26), axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 26), legend.spacing.y = unit(1, "mm"), 
        plot.title = element_text(size = 28), axis.text = element_text(size=26)) +
  geom_rect(data=xy, aes(xmin=0.5, xmax=1.5, ymin=1.94, ymax=18.33), color="gray80", fill="gray80", size=1) +
  geom_rect(data=xy, aes(xmin=1.5, xmax=2.5, ymin=2.46, ymax=51.65), color="gray80", fill="gray80", size=1) +
  geom_rect(data=xy, aes(xmin=2.5, xmax=3.5, ymin=18.15, ymax=240.81), color="gray80", fill="gray80", size=1) +
  geom_boxplot(ponderosa, mapping= aes(x=hclass, y=rootarea), fill="orchid3", outlier.shape = NA, width=0.6) +
  stat_boxplot(ponderosa, mapping= aes(x=hclass, y=rootarea), geom = "errorbar", width=0.35) +
  annotate(geom = "text", x=1, y=45, label="b", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=2, y=80, label="c", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=3, y=160, label="d", fontface='italic', color="black", size=10)




t.test(ponderosa150300m$rootarea, ponderosa150m$rootarea, alternative = "greater") #p=3e-10
t.test(ponderosa300m$rootarea, ponderosa150300m$rootarea, alternative = "greater") #p=0.007
t.test(ponderosa300m$rootarea, ponderosa150m$rootarea, alternative = "greater") #p=0.0005






#fir#############################################################################
frootarea <- subset(excavations, Species == "Fir")

frootarea150 <- subset(frootarea, hclass == "<150")
frootarea300 <- subset(frootarea, hclass != "<150")

#testing correlation between height and diameter and root area/mass
#RA v h
#RA = 1.641 + 0.054*h, R2=0.37 p=0.0002 RMSE=1.88 cm2
y <- frootarea150$rootarea
x <- frootarea150$height
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))


#RA = -29.863 + 0.247*h, R2=0.51 p=0.01 RMSE=7.74 cm2
y <- frootarea300$rootarea
x <- frootarea300$height
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))







#RW v h
#RW = 0.0069*h, R2=0.79 p=1.5e-11 RMSE=0.297 g
y <- frootarea150$rootmass
x <- frootarea150$height
z <- lm(y ~ 0 + x)
summary(z)
sqrt(mean(z$residuals^2))


#RW = -6.43 + 0.046*h, R2=0.51 p=0.012 RMSE=1.44 g
y <- frootarea300$rootmass
x <- frootarea300$height
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))






#RA v d
#RA = 0.604 + 1.948*d, R2=0.72 p=1.9e-09 RMSE=1.25 cm2
y <- frootarea150$rootarea
x <- frootarea150$diameter
z<- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))


#RA = 3.597*d, R2=0.87 p=1.6e-05 RMSE=7.045 cm2
y <- frootarea300$rootarea
x <- frootarea300$diameter
z <- lm(y ~ 0 + x)
summary(z)
sqrt(mean(z$residuals^2))






#RW v d
#RW = 0.208*d, R2=0.88 p=2e-15 RMSE=0.22 g
y <- frootarea150$rootmass
x <- frootarea150$diameter
z <- lm(y ~ 0 + x)
summary(z)
sqrt(mean(z$residuals^2))


#RW = 0.53*d, R2=0.76, p=0.0003 RMSE=1.51 g
y <- frootarea300$rootmass
x <- frootarea300$diameter
z <- lm(y ~ 0 + x)
summary(z)
sqrt(mean(z$residuals^2))








fig5g <- ggplot() +
  geom_boxplot(fir, mapping= aes(x=hclass, y=rootmass), fill="palegreen4", outlier.shape = NA, width=0.6) +
  stat_boxplot(fir, mapping= aes(x=hclass, y=rootmass), geom = "errorbar", width=0.35) +
  scale_y_log10(limits = c(0.01,100), breaks = c(0.01, 0.1, 1, 10, 100), expand=c(0.0001,0.0001),
                labels = c("0.01","0.1","1.0","10","100")) +
  labs(x="Height class (mm)", y="Root mass (g)", title = "Fir") +
  theme_classic() +
  theme(legend.title = element_blank(), legend.position = "none", axis.text.x = element_text(size = 26), 
        axis.title.x = element_text(size = 26), legend.text = element_text(size = 26),
        legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 28), axis.text = element_text(size=26), 
        axis.title = element_text(size=26)) +
  geom_rect(data=xy, aes(xmin=0.5, xmax=1.5, ymin=0.17, ymax=2.07), color="gray80", fill="gray80", size=1) +
  geom_rect(data=xy, aes(xmin=1.5, xmax=2.5, ymin=0.78, ymax=7.54), color="gray80", fill="gray80", size=1) +
  geom_boxplot(fir, mapping= aes(x=hclass, y=rootmass), fill="palegreen4", outlier.shape = NA, width=0.6) +
  stat_boxplot(fir, mapping= aes(x=hclass, y=rootmass), geom = "errorbar", width=0.35) +
  annotate(geom = "text", x=1, y=4, label="c", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=2, y=20, label="de", fontface='italic', color="black", size=10) 

t.test(fir150300m$rootmass, fir150m$rootmass, alternative = "greater") #p=0.02




fig5h <- ggplot() +
  geom_boxplot(fir, mapping= aes(x=hclass, y=rootarea), fill="palegreen4", outlier.shape = NA, width=0.6) +
  stat_boxplot(fir, mapping= aes(x=hclass, y=rootarea), geom = "errorbar", width=0.35) +
  scale_y_log10(expand = c(0.0001,0.0001), breaks = c(1,3,10,30,100), limits = c(1,100)) +
  labs(x="Height class (mm)", y=expression(paste("Root area (cm"^{2},")")), title = "Fir") +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_blank(), legend.position = "none",
        legend.background = element_blank(), axis.text.x = element_text(size = 26), 
        axis.title.x = element_text(size = 26), 
        axis.title.y = element_text(size = 26), legend.spacing.y = unit(1, "mm"), 
        plot.title = element_text(size = 28), axis.text = element_text(size=26)) +
  geom_rect(data=xy, aes(xmin=0.5, xmax=1.5, ymin=2.53, ymax=20.01), color="gray80", fill="gray80", size=1) +
  geom_rect(data=xy, aes(xmin=1.5, xmax=2.5, ymin=5.29, ymax=51.15), color="gray80", fill="gray80", size=1) +
  geom_boxplot(fir, mapping= aes(x=hclass, y=rootarea), fill="palegreen4", outlier.shape = NA, width=0.6) +
  stat_boxplot(fir, mapping= aes(x=hclass, y=rootarea), geom = "errorbar", width=0.35) +
  annotate(geom = "text", x=1, y=20, label="b", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=2, y=48, label="c", fontface='italic', color="black", size=10)

t.test(fir150300m$rootarea, fir150m$rootarea, alternative = "greater") #p=0.01




plot_grid(fig5a, fig5b, fig5c, fig5d, fig5e, fig5f, fig5g, fig5h, nrow = 4,
          labels = c("a", "b", "c", "d", "e", "f", "g", "h"),
          label_x = 0.93, label_y = 0.95, label_size = 40)






#comparison between species 
t.test(juniper150m$rootmass, pinyon150m$rootmass, alternative = "less") 
t.test(juniper150m$rootmass, ponderosa150m$rootmass, alternative = "greater") 
t.test(juniper150m$rootmass, fir150m$rootmass, alternative = "greater") 
t.test(pinyon150m$rootmass, ponderosa150m$rootmass, alternative = "greater") #p<0.05
t.test(pinyon150m$rootmass, fir150m$rootmass, alternative = "greater") #p<0.05
t.test(ponderosa150m$rootmass, fir150m$rootmass, alternative = "less") #p<0.05, fir sig greater


t.test(pinyon150300m$rootmass, ponderosa150300m$rootmass, alternative = "greater") #p<0.05
t.test(pinyon150300m$rootmass, fir150300m$rootmass, alternative = "greater") 
t.test(ponderosa150300m$rootmass, fir150300m$rootmass, alternative = "greater") 


t.test(pinyon150300m$rootmass, juniper150m$rootmass, alternative = "greater") #p<0.05
t.test(pinyon150300m$rootmass, ponderosa150m$rootmass, alternative = "greater") #p<0.05
t.test(pinyon150300m$rootmass, fir150m$rootmass, alternative = "greater") #p<0.05
t.test(ponderosa150300m$rootmass, juniper150m$rootmass, alternative = "greater") #p<0.05
t.test(ponderosa150300m$rootmass, pinyon150m$rootmass, alternative = "greater") #p<0.05
t.test(ponderosa150300m$rootmass, fir150m$rootmass, alternative = "greater") #p<0.05
t.test(fir150300m$rootmass, juniper150m$rootmass, alternative = "greater") #p<0.05
t.test(fir150300m$rootmass, pinyon150m$rootmass, alternative = "greater") #p<0.05
t.test(fir150300m$rootmass, ponderosa150m$rootmass, alternative = "greater") #p<0.05


t.test(ponderosa300m$rootmass, juniper150m$rootmass, alternative = "greater") #p<0.05
t.test(ponderosa300m$rootmass, pinyon150m$rootmass, alternative = "greater") #p<0.05
t.test(ponderosa300m$rootmass, fir150m$rootmass, alternative = "greater") #p<0.05


t.test(ponderosa300m$rootmass, pinyon150300m$rootmass, alternative = "greater") #p<0.05
t.test(ponderosa300m$rootmass, fir150300m$rootmass, alternative = "greater") #p<0.05






t.test(juniper150m$rootarea, pinyon150m$rootarea, alternative = "less") 
t.test(juniper150m$rootarea, ponderosa150m$rootarea, alternative = "greater") 
t.test(juniper150m$rootarea, fir150m$rootarea, alternative = "greater") 
t.test(pinyon150m$rootarea, ponderosa150m$rootarea, alternative = "greater") #p<0.05
t.test(pinyon150m$rootarea, fir150m$rootarea, alternative = "greater") #p<0.05
t.test(ponderosa150m$rootarea, fir150m$rootarea, alternative = "greater") 


t.test(juniper150300m$rootarea, pinyon150300m$rootarea, alternative = "greater") 
t.test(juniper150300m$rootarea, ponderosa150300m$rootarea, alternative = "greater") 
t.test(juniper150300m$rootarea, fir150300m$rootarea, alternative = "greater") 
t.test(pinyon150300m$rootarea, ponderosa150300m$rootarea, alternative = "greater") 
t.test(pinyon150300m$rootarea, fir150300m$rootarea, alternative = "greater") 
t.test(ponderosa150300m$rootarea, fir150300m$rootarea, alternative = "greater") 


t.test(pinyon150300m$rootarea, juniper150m$rootarea, alternative = "greater") #p<0.05
t.test(pinyon150300m$rootarea, ponderosa150m$rootarea, alternative = "greater") #p<0.05
t.test(pinyon150300m$rootarea, fir150m$rootarea, alternative = "greater") #p<0.05
t.test(ponderosa150300m$rootarea, juniper150m$rootarea, alternative = "greater") #p<0.05
t.test(ponderosa150300m$rootarea, pinyon150m$rootarea, alternative = "greater") #p<0.05
t.test(ponderosa150300m$rootarea, fir150m$rootarea, alternative = "greater") #p<0.05
t.test(fir150300m$rootarea, juniper150m$rootarea, alternative = "greater") #p<0.05
t.test(fir150300m$rootarea, pinyon150m$rootarea, alternative = "greater") #p<0.05
t.test(fir150300m$rootarea, ponderosa150m$rootarea, alternative = "greater") #p<0.05


t.test(ponderosa300m$rootarea, juniper150m$rootarea, alternative = "greater") #p<0.05
t.test(ponderosa300m$rootarea, pinyon150m$rootarea, alternative = "greater") #p<0.05
t.test(ponderosa300m$rootarea, fir150m$rootarea, alternative = "greater") #p<0.05


t.test(ponderosa300m$rootarea, pinyon150300m$rootarea, alternative = "greater") #p<0.05
t.test(ponderosa300m$rootarea, fir150300m$rootarea, alternative = "greater") #p<0.05






#figure s3######################################################################################################
#juniper####################################################################
absbreaks <- c(0,-100,-200,-300,-400, -500)
junroot$hclass <- gsub(">300", "300+", junroot$hclass)
junroot$hclass <- factor(junroot$hclass, levels = c("<150", "150-300", "300+"))


names(junroot150)[names(junroot150) == "hclass"] <- "hclass1"
names(junroot150300)[names(junroot150300) == "hclass"] <- "hclass1"
names(junroot300)[names(junroot300) == "hclass"] <- "hclass1"

junroot$hclass <- gsub("<150", "< 150", junroot$hclass)

figs3a <- ggplot() +
  geom_point(junroot, mapping = aes(x=rdiameter, y=depth, color=hclass), size=3) + 
  scale_color_manual(labels = c("< 150 mm","150-300 mm", "300+ mm"), values=c("black","snow3", "gray45")) +
  new_scale("color") +
  stat_smooth(junroot150, mapping = aes(x=rdiameter,y=depth, color=hclass1), method = "nls", 
              formula = y~a*log(x) + b, se=FALSE, 
              method.args =  list(start=list(a=203.67,b=-173.82)),size=3) +
  stat_smooth(junroot150300, mapping = aes(x=rdiameter,y=depth, color=hclass1), method = "nls", 
              formula = y~a*log(x) + b, se=FALSE, 
              method.args =  list(start=list(a=174.62,b=-347.181)),size=3) +
  stat_smooth(junroot300, mapping = aes(x=rdiameter,y=depth, color=hclass1), method = "nls", 
              formula = y~a*log(x) + b, se=FALSE, 
              method.args =  list(start=list(a=153.675,b=-388.353)),size=3) +
  scale_color_manual(labels = c("< 150 mm","150-300 mm", "300+ mm"), values = c("red3", "orangered2","violetred4")) +
  theme_classic() +
  labs(y="Rooting depth (mm)", title = "Juniper") +
  scale_x_continuous(expand = c(0,0), breaks = c(0,5,10,15,20,25), limits = c(0,26)) +
  scale_y_continuous(expand = c(0,0), breaks = absbreaks, labels = abs(absbreaks), limits = c(-500,5)) +
  theme(text = element_text(size = 26), axis.text = element_text(size=26),
        axis.title = element_text(size = 26), legend.position = c(0.85,0.3), legend.background = 
          element_rect(linetype = "solid", color = "black"),
        legend.title = element_blank(), legend.text = element_text(size = 26), plot.title = element_text(size = 28),
        axis.title.x = element_blank(), panel.spacing.x = unit(1.5, "lines"))



#polynomial component
#stat_smooth(junroot150, mapping = aes(x=rdiameter,y=depth), color="red3", method = "nls", 
 #           formula = y~a*(x^3) + b*(x^2) + c*(x) + d, se=FALSE, 
  #          method.args =  list(start=list(a=-39.1,b=41.84,c=262.18,d=-443.94)),size=2,
   #         xseq = seq(0,3, length=80)) +
  #stat_smooth(junroot150300, mapping = aes(x=rdiameter,y=depth), color="orangered2", method = "nls", 
   #           formula = y~a*(x^2) + b*(x) + c, se=FALSE, 
    #          method.args =  list(start=list(a=1.13,b=51.61,c=-329.167)),size=2,
     #         xseq = seq(0,3, length=80)) +
  #stat_smooth(junroot300, mapping = aes(x=rdiameter,y=depth), color="gray30", method = "nls", 
   #           formula = y~a*(x^2) + b*(x) + c, se=FALSE, 
    #          method.args =  list(start=list(a=-1.5,b=52.48,c=-369.88)),size=2,
     #         xseq = seq(0,3, length=80)) +
#annotate(geom = "text", label = (~R^2 ==~"0.86 ± 0.091, RMSE = 20.19 ± 10.02 mm"), x=15,y=-320, color="black", size=4.5) +
 # annotate(geom="text", label = "150-300: Depth = 1.13*d² + 51.61*d - 329.167", x=14.5, y=-360, color="black",
  #         size=4.5) +
  #annotate(geom = "text", label = (~R^2 ==~"0.93 ± 0.008, RMSE = 17.56 ± 4.07 mm"), x=15,y=-380, color="black", size=4.5) +
  #annotate(geom="text", label = ">300: Depth = -1.5*d² + 52.48*d - 369.88", x=15, y=-420, color="black",
   #        size=4.5) +
  #annotate(geom = "text", label = (~R^2 ==~"0.89 ± 0.019, RMSE = 30.65 ± 20.41 mm"), x=15,y=-440, color="black", size=4.5)


#pinon##############################################################################
pinroot$hclass <- gsub(">300", "300+", pinroot$hclass)
pinroot$hclass <- factor(pinroot$hclass, levels = c("<150", "150-300", "300+"))

names(pinroot150)[names(pinroot150) == "hclass"] <- "hclass1"
names(pinroot150300)[names(pinroot150300) == "hclass"] <- "hclass1"
names(pinroot300)[names(pinroot300) == "hclass"] <- "hclass1"

pinroot$hclass <- gsub("<150", "< 150", pinroot$hclass)

absbreakspin <- c(0,-100,-200,-300,-400)
figs3b <- ggplot() +
  geom_point(data=pinroot, aes(x=rdiameter, y=depth, color=hclass), show.legend=FALSE, size=3) +
  scale_color_manual(values = c("black","snow3", "gray45")) +
  new_scale("color") +
  stat_smooth(pinroot150, mapping = aes(x=rdiameter,y=depth, color=hclass1), method = "nls", 
              formula = y~a*log(x) + b, se=FALSE, 
              method.args =  list(start=list(a=150.354,b=-203.5)),size=3) +
  stat_smooth(pinroot150300, mapping = aes(x=rdiameter,y=depth, color=hclass1), method = "nls", 
              formula = y~a*log(x) + b, se=FALSE, 
              method.args =  list(start=list(a=163.832,b=-323.482)),size=3) +
  stat_smooth(pinroot300, mapping = aes(x=rdiameter,y=depth, color=hclass1), method = "nls", 
              formula = y~a*log(x) + b, se=FALSE, 
              method.args =  list(start=list(a=153.3,b=-411.867)),size=3) +
    scale_color_manual(labels = c("< 150","150-300", "300+"), values = c("orange3", "wheat2", "tan1")) +
    theme_classic() +
  theme_classic() +
  labs(x="Root diameter (mm)", y="Rooting depth (mm)", title = "Piñon pine") +
  scale_x_continuous(expand = c(0,0), breaks = c(0,5,10,15,20,25), limits = c(0,26)) +
  scale_y_continuous(expand = c(0,0), breaks = absbreakspin, labels = abs(absbreakspin), limits = c(-400,5)) +
  theme(text = element_text(size = 26), axis.text = element_text(size=26),
        axis.title = element_blank(),
        legend.title = element_blank(), legend.position = c(0.85,0.35), legend.background = 
          element_rect(linetype = "solid", color = "black"),
        legend.text = element_text(size = 26), plot.title = element_text(size = 28),
        panel.spacing.x = unit(1.5, "lines"))


#polynomial component
#stat_smooth(pinroot150, mapping = aes(x=rdiameter,y=depth), color="orange3", method = "nls", 
 #           formula = y~a*(x^3) + b*(x^2) + c*(x) + d, se=FALSE, 
  #          method.args =  list(start=list(a=-3.87,b=4.6,c=94.44,d=-275.72)),size=2,
   #         xseq = seq(0,3, length=80)) +
  #stat_smooth(pinroot150300, mapping = aes(x=rdiameter,y=depth), color="tan1", method = "nls", 
   #           formula = y~a*(x^2) + b*(x) + c, se=FALSE, 
    #          method.args =  list(start=list(a=-2.72,b=63.14,c=-323.35)),size=2,
     #         xseq = seq(0,3, length=80)) +
  #stat_smooth(pinroot300, mapping = aes(x=rdiameter,y=depth), color="gray35", method = "nls", 
   #           formula = y~a*(x^2) + b*(x) + c, se=FALSE, 
    #          method.args =  list(start=list(a=0.297,b=16.96,c=-251.504)),size=2,
     #         xseq = seq(0,3, length=80)) +







#ponderosa#############################################################################################
names(ponroot150)[names(ponroot150) == "hclass"] <- "hclass1"
names(ponroot150300)[names(ponroot150300) == "hclass"] <- "hclass1"
names(ponroot300)[names(ponroot300) == "hclass"] <- "hclass1"

ponroot$hclass <- gsub("<150", "< 150", ponroot$hclass)

figs3c <- ggplot() +
  geom_point(data=ponroot, aes(x=rdiameter, y=depth, color=hclass), show.legend=FALSE, size=3) +
  scale_color_manual(values = c("black","snow3", "gray45")) +
  new_scale("color") +
  stat_smooth(ponroot150, mapping = aes(x=rdiameter,y=depth, color=hclass1), method = "nls", 
              formula = y~a*log(x) + b, se=FALSE, 
              method.args =  list(start=list(a=86.841,b=-91.07)),size=3) +
  stat_smooth(ponroot150300, mapping = aes(x=rdiameter,y=depth, color=hclass1), method = "nls", 
              formula = y~a*log(x) + b, se=FALSE, 
              method.args =  list(start=list(a=128.91,b=-227.505)),size=3) +
  stat_smooth(ponroot300, mapping = aes(x=rdiameter,y=depth, color=hclass1), method = "nls", 
              formula = y~a*log(x) + b, se=FALSE, 
              method.args =  list(start=list(a=152.411,b=-363.804)),size=3) +
  scale_color_manual(labels = c("< 150","150-300", "300+"), values = c("darkviolet", "orchid3", "mediumpurple")) +
  theme_classic() +
  labs(x="Root diameter (mm)", y="Rooting depth (mm)", title = "Ponderosa pine") +
  scale_x_continuous(expand = c(0,0), breaks = c(0,5,10,15,20,25), limits = c(0,26)) +
  scale_y_continuous(expand = c(0,0), breaks = absbreaks, labels = abs(absbreaks), limits = c(-515,5)) +
  theme(text = element_text(size = 22), axis.text = element_text(size=26),
        axis.title = element_text(size = 26),
        legend.title = element_blank(), legend.position = c(0.85,0.4), legend.background = 
          element_rect(linetype = "solid", color = "black"),
        legend.text = element_text(size = 26), plot.title = element_text(size = 28),
        panel.spacing.x = unit(1.5, "lines"))


#polynomial component
#stat_smooth(ponroot150, mapping = aes(x=rdiameter,y=depth), color="darkviolet", method = "nls", 
 #           formula = y~a*(x^3) + b*(x^2) + c*(x) + d, se=FALSE, 
  #          method.args =  list(start=list(a=-2.8,b=-30.68,c=159.29,d=-217.26)),size=2,
   #         xseq = seq(0,3, length=80)) +
  #stat_smooth(ponroot150300, mapping = aes(x=rdiameter,y=depth), color="orchid3", method = "nls", 
   #           formula = y~a*(x^2) + b*(x) + c, se=FALSE, 
    #          method.args =  list(start=list(a=-12.29,b=122.375,c=-336.24)),size=2,
     #         xseq = seq(0,3, length=80)) +
  #stat_smooth(ponroot300, mapping = aes(x=rdiameter,y=depth), color="gray40", method = "nls", 
   #           formula = y~a*(x^2) + b*(x) + c, se=FALSE, 
    #          method.args =  list(start=list(a=-8.3,b=101.61,c=-461.63)),size=2,
     #         xseq = seq(0,3, length=80)) +




#fir######################################################################################
firroot$hclass <- factor(firroot$hclass, levels = c("<150", "150-300", ">300"))
names(firroot150)[names(firroot150) == "hclass"] <- "hclass1"
names(firroot150300)[names(firroot150300) == "hclass"] <- "hclass1"

firroot$hclass <- gsub("<150", "< 150", firroot$hclass)

figs3d <- ggplot() +
  geom_point(data=firroot, aes(x=rdiameter, y=depth, color=hclass), show.legend = FALSE, size=3) + 
  scale_color_manual(values = c("black","snow3")) +
  new_scale("color") +
  stat_smooth(firroot150, mapping = aes(x=rdiameter,y=depth, color=hclass1), method = "nls", 
              formula = y~a*log(x) - b, se=FALSE, 
              method.args =  list(start=list(a=160.658,b=147.582)),size=3) +
  stat_smooth(firroot150300, mapping = aes(x=rdiameter,y=depth, color=hclass1), method = "nls", 
              formula = y~a*log(x) - b, se=FALSE, 
              method.args =  list(start=list(a=172.085,b=320.819)),size=3) +
  scale_color_manual(labels = c("< 150","150-300"), values = c("palegreen4", "palegreen3")) +
  theme_classic() +
  labs(x="Root diameter (mm)", y="Rooting depth (mm)", title = "Fir") +
  scale_x_continuous(expand = c(0,0), breaks = c(0,5,10,15,20,25), limits = c(0,26)) +
  scale_y_continuous(expand = c(0,0), breaks = absbreakspin, labels = abs(absbreakspin), limits = c(-320,5)) +
  theme(text = element_text(size = 26), axis.text = element_text(size=26),
        axis.title = element_text(size = 26),
        legend.title = element_blank(), legend.position = c(0.85,0.4), legend.background = 
          element_rect(linetype = "solid", color = "black"), axis.title.y = element_blank(),
        legend.text = element_text(size = 26), plot.title = element_text(size = 28),
        panel.spacing.x = unit(1.5, "lines"))

#polynomial component
#stat_smooth(firroot150, mapping = aes(x=rdiameter,y=depth), color="palegreen4", method = "nls", 
 #           formula = y~a*(x^3) + b*(x^2) + c*(x) + d, se=FALSE, 
  #          method.args =  list(start=list(a=1.09,b=1.91,c=99.31,d=-265.12)), size=2,
   #         xseq = seq(0,3, length=80)) +
  #stat_smooth(firroot150300, mapping = aes(x=rdiameter,y=depth), color="palegreen3", method = "nls", 
   #           formula = y~a*(x^2) + b*(x) + c, se=FALSE, 
    #          method.args =  list(start=list(a=-3.1,b=87.88,c=-365.25)), size=2,
     #         xseq = seq(0,3, length=80)) +




plot_grid(figs3a, figs3b, figs3c, figs3d, nrow = 2,
          labels = c("a", "b", "c", "d"),
          label_x = 0.95, label_y = 0.97, label_size = 40)











#figure 6#########################################################################################
#juniper#########################################################################################
jmaxdep <- subset(excavations, Species == "Juniper")
jmaxdep$maxroot305 <- alljuncoef1$rootdep35
jmaxdep$maxroot405 <- alljuncoef1$rootdep45

jmaxdep150 <- subset(jmaxdep, hclass == "<150")
jmaxdep300 <- subset(jmaxdep, hclass != "<150")

#testing correlation between height and diameter and maximum rooting depth
#seedlings <150mm - max root depth v d
#RD = -139-70.56*(d), R2=0.02 p=0.3 RMSE=137.35 mm
x <- jmaxdep150$diameter
y <- jmaxdep150$maxroot405
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))

#seedlings <150mm - max root depth v h
#RD = -264.734-1.492*h, R2=0.0004 p=0.95 RMSE=147.282 mm
x <- jmaxdep150$height
y <- jmaxdep150$maxroot405
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))




#seedlings 300mm - max root depth v d
#RD = -204.364 - 14.548*d, R2=0.25 p=0.11 RMSE=79.97 mm
x <- jmaxdep300$diameter
y <- jmaxdep300$maxroot305
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))

#seedlings 300mm - max root depth v h
#RD = -278.513 - 0.139*h, R2=0.01 p=0.78 RMSE=99.38 mm
x <- jmaxdep300$height
y <- jmaxdep300$maxroot305
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))




ggscatter(jmaxdep, x="height", y="maxroot", color = "hclass")
ggscatter(jmaxdep, x="diameter", y="maxroot", color = "hclass")



absbreaks <- c(0,-250,-500,-750,-1000)


fig6a <- ggplot() +
  geom_boxplot(juniper, mapping= aes(x=hclass, y=maxroot), fill="orangered2", outlier.shape = NA, width=0.6) +
  stat_boxplot(juniper, mapping= aes(x=hclass, y=maxroot), geom = "errorbar", width=0.4) +
  scale_y_continuous(expand = c(0,0), breaks = absbreaks, labels = abs(absbreaks), limits = c(-750,0)) +
  labs(x="Height class (mm)", y="Maximum rooting depth (mm)", title = "Juniper") +
  theme_classic() +
  theme(legend.title = element_blank(), legend.position = c(0.2,0.2), 
        legend.text = element_text(size = 26), legend.key.width = unit(2.5, "line"),
        legend.background = element_rect(linetype = "solid", color = "black"), 
        axis.text.x = element_text(size = 26), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 26), legend.spacing.y = unit(1, "mm"), 
        plot.title = element_text(size = 28), axis.text = element_text(size=26)) +
  annotate(geom = "text", x=1, y=-130, label="bcd", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=2, y=-190, label="n=6", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=3, y=-220, label="n=2", fontface='italic', color="black", size=10)

t.test(juniper150300m$maxroot, juniper150m$maxroot, alternative = "less") 
t.test(juniper300m$maxroot, juniper150300m$maxroot, alternative = "less") 
t.test(juniper300m$maxroot, juniper150m$maxroot, alternative = "less") #p=0.0001








#pinon######################################################################################
pinmaxdep <- subset(excavations, Species == "Pinon")
pinmaxdep$maxroot305 <- allpincoef1$rootdep35
pinmaxdep$maxroot405 <- allpincoef1$rootdep45

pinmaxdep150 <- subset(pinmaxdep, hclass == "<150")
pinmaxdep150300 <- subset(pinmaxdep, hclass == "150-300")
pinmaxdep300 <- subset(pinmaxdep, hclass != "<150")

#testing correlation between height and diameter and maximum rooting depth
#seedlings <150mm - max root depth v d
#RD = -100.87-33.95*d, R2=0.11 p=0.02 RMSE=107.8 mm
x <- pinmaxdep150$diameter
y <- pinmaxdep150$maxroot405
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))

#seedlings <150mm - max root depth v h
#RD = -172.52-0.57*h, R2=0.03 p=0.34 RMSE=114.31 mm
x <- pinmaxdep150$height
y <- pinmaxdep150$maxroot405
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))




#seedlings 150-300mm - max root depth v d
#RD = -144.31-25.53*d, R2=0.12 p=0.11 RMSE=86.83 mm
x <- pinmaxdep150300$diameter
y <- pinmaxdep150300$maxroot305
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))

#seedlings 150-300mm - max root depth v h
#RD = -116.657-0.938*h, R2=0.02 p=0.3 RMSE=91.83 mm
x <- pinmaxdep150300$height
y <- pinmaxdep150300$maxroot305
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))




#seedlings >300mm - max root depth v d
#RD = -250.14-5.681*d, R2=0.02 p=0.6 RMSE=112.2 mm
x <- pinmaxdep300$diameter
y <- pinmaxdep300$maxroot305
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))

#seedlings >300mm - max root depth v h
#RD = -1.1*h, R2=0.76 p=2.8e-07 RMSE=148.15 mm
x <- pinmaxdep300$height
y <- pinmaxdep300$maxroot305
z <- lm(y ~ 0 + x)
summary(z)
sqrt(mean(z$residuals^2))



ggscatter(pinmaxdep, x="height", y="maxroot305", color = "hclass")
ggscatter(pinmaxdep, x="diameter", y="maxroot305", color = "hclass")


classpin1$hclass <- factor(classpin1$hclass, levels = c("<150", "150-300", ">300"))
pinmaxdep$hclass <- factor(pinmaxdep$hclass, levels = c("<150", "150-300", ">300"))




fig6b <- ggplot() +
  geom_boxplot(pinon, mapping= aes(x=hclass, y=maxroot), fill="orange3", outlier.shape = NA, width=0.6) +
  stat_boxplot(pinon, mapping= aes(x=hclass, y=maxroot), geom = "errorbar", width=0.4) +
  scale_y_continuous(expand = c(0,0), breaks = absbreaks, labels = abs(absbreaks), limits = c(-1000,70)) +
  labs(x="Height class (mm)", y="Maximum rooting depth (mm)", title = "Piñon pine") +
  theme_classic() +
  theme(legend.title = element_blank(), legend.position = "none", 
        axis.text.x = element_text(size = 26), axis.title = element_blank(), 
        plot.title = element_text(size = 28), axis.text = element_text(size=26)) +
  geom_rect(data=xy, aes(xmin=0.5, xmax=1.5, ymin=-590.77, ymax=-126.67), color="gray80", fill="gray80", size=1) +
  geom_rect(data=xy, aes(xmin=1.5, xmax=2.5, ymin=-957, ymax=-165), color="gray80", fill="gray80", size=1) +
  geom_rect(data=xy, aes(xmin=2.5, xmax=3.5, ymin=-957, ymax=-165), color="gray80", fill="gray80", size=1) +
  geom_boxplot(pinon, mapping= aes(x=hclass, y=maxroot), fill="orange3", outlier.shape = NA, width=0.6) +
  stat_boxplot(pinon, mapping= aes(x=hclass, y=maxroot), geom = "errorbar", width=0.4) +
  annotate(geom = "text", x=1, y=27, label="b", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=2, y=-130, label="cd", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=3, y=-20, label="n=4", fontface='italic', color="black", size=10)

t.test(pinyon150300m$maxroot, pinyon150m$maxroot, alternative = "less") #p=0.005
t.test(pinyon300m$maxroot, pinyon150300m$maxroot, alternative = "less") 
t.test(pinyon300m$maxroot, pinyon150m$maxroot, alternative = "less") #p=0.02







#ponderosa################################################################################
ponmaxdep <- subset(excavations, Species == "Ponderosa")
ponmaxdep$maxroot305 <- allponcoef1$rootdep35
ponmaxdep$maxroot405 <- allponcoef1$rootdep45

ponmaxdep150 <- subset(ponmaxdep, hclass == "<150")
ponmaxdep150300 <- subset(ponmaxdep, hclass == "150-300")
ponmaxdep300 <- subset(ponmaxdep, hclass == ">300")

#testing correlation between height and diameter and maximum rooting depth
#seedlings <150mm - max root depth v d
#RD = -36.43-35.78*d, R2=0.39 p=1.7e-09 RMSE=52.52 mm
x <- ponmaxdep150$diameter
y <- ponmaxdep150$maxroot405
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))

#seedlings <150mm - max root depth v h
#RD = -51.52-1.288*h, R2=0.34 p=4.3e-08 RMSE=54.92 mm
x <- ponmaxdep150$height
y <- ponmaxdep150$maxroot405
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))



#seedlings 150-300mm - max root depth v d
#RD = -242.94-5.65*d, R2=0.01 p=0.71 RMSE=103.91 mm
x <- ponmaxdep150300$diameter
y <- ponmaxdep150300$maxroot305
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))

#seedlings 150-300mm - max root depth v h
#RD = -264.7-0.064*h, R2=0.0005 p=0.9 RMSE=104.2 mm
x <- ponmaxdep150300$height
y <- ponmaxdep150300$maxroot305
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))



#seedlings >300mm - max root depth v d
#RD = -212.68-12.57*d, R2=0.74 p=0.0009 RMSE=38.06 mm
x <- ponmaxdep300$diameter
y <- ponmaxdep300$maxroot305
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))

#seedlings >300mm - max root depth v h
#RD = -168.85-0.439*h, R2=0.41 p=0.03 RMSE=56.96 mm
x <- ponmaxdep300$height
y <- ponmaxdep300$maxroot305
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))

 
ggscatter(ponmaxdep, x="height", y="maxroot", color = "hclass")
ggscatter(ponmaxdep, x="diameter", y="maxroot", color = "hclass")


classpon1$hclass <- factor(classpon1$hclass, levels = c("<150", "150-300", ">300"))



fig6c <- ggplot() +
  geom_boxplot(ponderosa, mapping= aes(x=hclass, y=maxroot), fill="orchid3", outlier.shape = NA, width=0.6) +
  stat_boxplot(ponderosa, mapping= aes(x=hclass, y=maxroot), geom = "errorbar", width=0.4) +
  scale_y_continuous(expand = c(0,0), breaks = absbreaks, labels = abs(absbreaks), limits = c(-1100,50)) +
  labs(x="Height class (mm)", y="Maximum rooting depth (mm)", title = "Ponderosa pine") +
  theme_classic() +
  theme(legend.title = element_blank(), legend.position = "none", 
        axis.text.x = element_text(size = 26), 
        axis.title.x = element_text(size = 26), axis.title.y = element_text(size = 26),
        plot.title = element_text(size = 28), axis.text = element_text(size=26)) +
  geom_rect(data=xy, aes(xmin=0.5, xmax=1.5, ymin=-362.03, ymax=-71.49), color="gray80", fill="gray80", size=1) +
  geom_rect(data=xy, aes(xmin=2.5, xmax=3.5, ymin=-1090.07, ymax=-295.89), color="gray80", fill="gray80", size=1) +
  geom_boxplot(ponderosa, mapping= aes(x=hclass, y=maxroot), fill="orchid3", outlier.shape = NA, width=0.6) +
  stat_boxplot(ponderosa, mapping= aes(x=hclass, y=maxroot), geom = "errorbar", width=0.4) +
  annotate(geom = "text", x=1, y=30, label="a", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=2, y=-60, label="c", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=3, y=-225, label="d", fontface='italic', color="black", size=10)

t.test(ponderosa150300m$maxroot, ponderosa150m$maxroot, alternative = "less") #p=1e-06
t.test(ponderosa300m$maxroot, ponderosa150300m$maxroot, alternative = "less") #p=0.007
t.test(ponderosa300m$maxroot, ponderosa150m$maxroot, alternative = "less") #p=3.4e-06






#fir###############################################################################################
fmaxdep <- subset(excavations, Species == "Fir")
fmaxdep$maxroot305 <- allfircoef1$rootdep35
fmaxdep$maxroot405 <- allfircoef1$rootdep45

fmaxdep150 <- subset(fmaxdep, hclass == "<150")
fmaxdep150300 <- subset(fmaxdep, hclass == "150-300")

#testing correlation between height and diameter and maximum rooting depth
#seedlings <150mm - max root depth v d
#RD = -117.28-34.84*d, R2=0.07 p=0.08 RMSE=109.92 mm
x <- fmaxdep150$diameter
y <- fmaxdep150$maxroot405
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))

#testing correlation between height and diameter and maximum rooting depth
#seedlings <150mm - max root depth v h
#RD = -102.74-1.38*h, R2=0.08 p=0.07 RMSE=109.48 mm
x <- fmaxdep150$height
y <- fmaxdep150$maxroot405
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))



#seedlings 150-300mm - max root depth v d
#RD = -63.19*d, R2=0.92 p=1.5e-06 RMSE=93.33 mm
x <- fmaxdep150300$diameter
y <- fmaxdep150300$maxroot305
z <- lm(y ~ 0 + x)
summary(z)
sqrt(mean(z$residuals^2))

#seedlings 150-300mm - max root depth v h
#RD = -62.23-1.365*h, R2=0.1 p=0.37 RMSE=144.17 mm
x <- fmaxdep150300$height
y <- fmaxdep150300$maxroot305
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))


ggscatter(fmaxdep, x="height", y="maxroot305", color = "hclass")
ggscatter(fmaxdep, x="diameter", y="maxroot405", color = "hclass")





fig6d <- ggplot() +
  geom_boxplot(fir, mapping= aes(x=hclass, y=maxroot), fill="palegreen4", outlier.shape = NA, width=0.6) +
  stat_boxplot(fir, mapping= aes(x=hclass, y=maxroot), geom = "errorbar", width=0.4) +
  scale_y_continuous(expand = c(0,0), breaks = absbreaks, labels = abs(absbreaks), limits = c(-500,0)) +
  labs(x="Height class (mm)", y="Maximum rooting depth (mm)", title = "Fir") +
  theme_classic() +
  theme(legend.title = element_blank(), legend.position = "none", 
        axis.text.x = element_text(size = 26), 
        axis.title.x = element_text(size = 26), axis.title.y = element_blank(),
        plot.title = element_text(size = 28), axis.text = element_text(size=26)) +
  geom_rect(data=xy, aes(xmin=1.5, xmax=2.5, ymin=-471.73, ymax=-276.53), color="gray80", fill="gray80", size=1) +
  geom_boxplot(fir, mapping= aes(x=hclass, y=maxroot), fill="palegreen4", outlier.shape = NA, width=0.6) +
  stat_boxplot(fir, mapping= aes(x=hclass, y=maxroot), geom = "errorbar", width=0.4) +
  annotate(geom = "text", x=1, y=-30, label="b", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=2, y=-140, label="bcd", fontface='italic', color="black", size=10)

t.test(fir150300m$maxroot, fir150m$maxroot, alternative = "less") 





plot_grid(fig6a, fig6b, fig6c, fig6d, nrow = 2,
          labels = c("a", "b", "c", "d"),
          label_x = 0.92, label_y = 0.95, label_size = 40)



#comparisons between species
t.test(juniper150m$maxroot, pinyon150m$maxroot, alternative = "less") 
t.test(juniper150m$maxroot, ponderosa150m$maxroot, alternative = "less") #p=0.005
t.test(juniper150m$maxroot, fir150m$maxroot, alternative = "less") 
t.test(pinyon150m$maxroot, ponderosa150m$maxroot, alternative = "less") #p=0.0001
t.test(pinyon150m$maxroot, fir150m$maxroot, alternative = "less") 
t.test(ponderosa150m$maxroot, fir150m$maxroot, alternative = "greater") #p=0.0003, fir greater root dep


t.test(pinyon150300m$maxroot, ponderosa150300m$maxroot, alternative = "less") 
t.test(pinyon150300m$maxroot, fir150300m$maxroot, alternative = "greater") 
t.test(ponderosa150300m$maxroot, fir150300m$maxroot, alternative = "greater") 


t.test(pinyon150300m$maxroot, juniper150m$maxroot, alternative = "greater")
t.test(pinyon150300m$maxroot, ponderosa150m$maxroot, alternative = "less") #p<0.0001
t.test(pinyon150300m$maxroot, fir150m$maxroot, alternative = "less") #p=0.01
t.test(ponderosa150300m$maxroot, juniper150m$maxroot, alternative = "greater")
t.test(ponderosa150300m$maxroot, pinyon150m$maxroot, alternative = "less") #p=0.04
t.test(ponderosa150300m$maxroot, fir150m$maxroot, alternative = "less") #p=0.04
t.test(fir150300m$maxroot, juniper150m$maxroot, alternative = "less")
t.test(fir150300m$maxroot, pinyon150m$maxroot, alternative = "less")
t.test(fir150300m$maxroot, ponderosa150m$maxroot, alternative = "less") #p=0.003


t.test(ponderosa300m$maxroot, juniper150m$maxroot, alternative = "less")
t.test(ponderosa300m$maxroot, pinyon150m$maxroot, alternative = "less") #p=0.0002
t.test(ponderosa300m$maxroot, fir150m$maxroot, alternative = "less") #p=0.0002


t.test(ponderosa300m$maxroot, pinyon150300m$maxroot, alternative = "less")
t.test(ponderosa300m$maxroot, fir150300m$maxroot, alternative = "less")










#figure 7###############################################################################################
juniper$rsarea <- (juniper$rootarea)/(juniper$leafarea+juniper$stemarea)
pinon$rsarea <- (pinon$rootarea)/(pinon$leafarea+pinon$stemarea)
ponderosa$rsarea <- (ponderosa$rootarea)/(ponderosa$leafarea+ponderosa$stemarea)
fir$rsarea <- (fir$rootarea)/(fir$leafarea+fir$stemarea)

excavations$rsarea <- (excavations$rootarea)/(excavations$leafarea+excavations$stemarea)

#juniper###############################################################################################
jrs <- subset(excavations, Species == "Juniper")
jrs150 <- subset(jrs, hclass == "<150")
jrs300 <- subset(jrs, hclass != "<150")

#rsarea v h
#RS = -0.0006*h +0.366 R2=0.006 p=0.83 RMSE=0.21
x <- jrs150$height
y <- jrs150$rsarea
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))

#rsarea v d
#RS = -0.0248*d+0.367 R2=0.008, p=0.8 RMSE=0.21
x <- jrs150$diameter
y <- jrs150$rsarea
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))


#rsweight v h
#RS = 0.132 + 0.0073*h R2=0.08 p=0.43
x <- jrs150$height
y <- jrs150$rsweight
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))

#rsweight v d
#RS = 1.44 - 0.234*d R2=0.05, p=0.51 RMSE=0.734
x <- jrs150$diameter
y <- jrs150$rsweight
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))





#rsarea v h
#RS = 386.658 + 1.581*h R2=-0.09 p=0.66
x <- jrs300$height
y <- jrs300$rsarea
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))

#rsarea v d
#sa = 148.9 + 168.4*d R2=0.09, p=0.2
x <- jrs300$diameter
y <- jrs300$rsarea
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))


#rsweight v h
#stem mass = 0.0328 + 0.0025*h R2=0.187 p=0.12
x <- jrs300$height
y <- jrs300$rsweight
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))

#rsweight v d
#sm = 0.033 + 0.108*d R2=0.26, p=0.07
x <- jrs300$diameter
y <- jrs300$rsweight
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))



ggscatter(jmaxdep, x="height", y="stemarea", color = "hclass")
ggscatter(jmaxdep, x="height", y="stemmass", color = "hclass")
ggscatter(jmaxdep, x="diameter", y="stemarea", color = "hclass")
ggscatter(jmaxdep, x="diameter", y="stemmass", color = "hclass")



classjun1$rsarea <- (classjun1$rootarea)/(classjun1$leafarea+classjun1$stemarea)
classjun1$rsweight <- (classjun1$rootmass)/(classjun1$leafmass+classjun1$stemmass)



fig7a <- ggplot() +
  geom_boxplot(juniper, mapping= aes(x=hclass, y=rsweight), fill="orangered2", outlier.shape = NA, width=0.6) +
  stat_boxplot(juniper, mapping= aes(x=hclass, y=rsweight), geom = "errorbar", width=0.35) +
  scale_y_continuous(expand = c(0,0), breaks = c(0,0.5,1,1.5,2), limits = c(0,2),
                     labels = c("0.0","","1.0","","2.0")) +
  labs(y="R:S by weight", title = "Juniper") +
  theme_classic() +
  theme(legend.title = element_blank(), legend.position = c(0.62,0.87), legend.background = 
          element_rect(linetype = "solid", color = "black"), legend.text = element_text(size = 26),
        axis.title.x = element_blank(), legend.spacing.y = unit(1, "mm"), legend.key.width = unit(2.5, "line"),
        plot.title = element_text(size = 28), axis.text = element_text(size=26), axis.title = element_text(size=26)) +
  geom_rect(data=xy, aes(xmin=0.5, xmax=1.5, ymin=0.795, ymax=0.937), color="gray80", fill="gray80", size=1) +
  geom_rect(data=xy, aes(xmin=1.5, xmax=2.5, ymin=0.51, ymax=0.545), color="gray80", fill="gray80", size=1) +
  geom_rect(data=xy, aes(xmin=2.5, xmax=3.5, ymin=0.51, ymax=0.545), color="gray80", fill="gray80", size=1) +
  geom_boxplot(juniper, mapping= aes(x=hclass, y=rsweight), fill="orangered2", outlier.shape = NA, width=0.6) +
  stat_boxplot(juniper, mapping= aes(x=hclass, y=rsweight), geom = "errorbar", width=0.35) +
  annotate(geom = "text", x=1, y=0.97, label="abce", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=2, y=0.7, label="n=6", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=3, y=0.9, label="n=2", fontface='italic', color="black", size=10) 

t.test(juniper150m$rsweight, juniper150300m$rsweight, alternative = "greater") #p=0.04
t.test(juniper300m$rsweight, juniper150300m$rsweight, alternative = "greater") 
t.test(juniper150m$rsweight, juniper300m$rsweight, alternative = "greater") 


juniper$rsarea <- juniper$rootarea/(juniper$leafarea+juniper$stemarea)


fig7b <- ggplot() +
  geom_boxplot(juniper, mapping= aes(x=hclass, y=rsarea), fill="orangered2", outlier.shape = NA, width=0.6) +
  stat_boxplot(juniper, mapping= aes(x=hclass, y=rsarea), geom = "errorbar", width=0.35) +
  scale_y_continuous(expand = c(0,0), breaks = c(0,0.25,0.5,0.75,1), limits = c(0,1),
                     labels = c("0.0","","0.5","","1.0")) +
  labs(y="R:S by area", title="Juniper") +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_blank(), legend.position = "none",
        legend.background = element_blank(), plot.title = element_text(size = 28),
        axis.text.x = element_text(size = 26), axis.title.x = element_blank(), legend.spacing.y = unit(1, "mm"),
        axis.text = element_text(size=26), axis.title = element_text(size=26)) +
  geom_rect(data=xy, aes(xmin=0.5, xmax=1.5, ymin=0.222, ymax=0.252), color="gray80", fill="gray80", size=1) +
  geom_boxplot(juniper, mapping= aes(x=hclass, y=rsarea), fill="orangered2", outlier.shape = NA, width=0.6) +
  stat_boxplot(juniper, mapping= aes(x=hclass, y=rsarea), geom = "errorbar", width=0.35) +
  annotate(geom = "text", x=1, y=0.7, label="b", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=2, y=0.5, label="n=6", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=3, y=0.25, label="n=2", fontface='italic', color="black", size=10)

t.test(juniper150m$rsarea, juniper150300m$rsarea, alternative = "greater") 
t.test(juniper150300m$rsarea, juniper300m$rsarea, alternative = "greater") 
t.test(juniper150m$rsarea, juniper300m$rsarea, alternative = "greater") #p=0.02






#pinon###############################################################################################
pinrs <- subset(excavations, Species == "Pinon")
pinrs150 <- subset(pinrs, hclass == "<150")
pinrs150300 <- subset(pinrs, hclass == "150-300")
pinrs300 <- subset(pinrs, hclass != "<150")

#rsarea v h
#stem area = 386.658 + 1.581*h R2=-0.09 p=0.66
x <- pinrs150$height
y <- pinrs150$rsarea
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))

#rsarea v d
#sa = 148.9 + 168.4*d R2=0.09, p=0.2
x <- pinrs150$diameter
y <- pinrs150$rsarea
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))


#rsweight v h
#stem mass = 0.0328 + 0.0025*h R2=0.187 p=0.12
x <- pinrs150$height
y <- pinrs150$rsweight
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))

#rsweight v d
#sm = 0.033 + 0.108*d R2=0.26, p=0.07
x <- pinrs150$diameter
y <- pinrs150$rsweight
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))





#rsarea v h
#stem area = 386.658 + 1.581*h R2=-0.09 p=0.66
x <- pinrs150300$height
y <- pinrs150300$rsarea
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))

#rsarea v d
#sa = 148.9 + 168.4*d R2=0.09, p=0.2
x <- pinrs150300$diameter
y <- pinrs150300$rsarea
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))


#rsweight v h
#stem mass = 0.0328 + 0.0025*h R2=0.187 p=0.12
x <- pinrs150300$height
y <- pinrs150300$rsweight
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))

#rsweight v d
#sm = 0.033 + 0.108*d R2=0.26, p=0.07
x <- pinrs150300$diameter
y <- pinrs150300$rsweight
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))





#rsarea v h
#stem area = 386.658 + 1.581*h R2=-0.09 p=0.66
x <- pinrs300$height
y <- pinrs300$rsarea
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))

#rsarea v d
#sa = 148.9 + 168.4*d R2=0.09, p=0.2
x <- pinrs300$diameter
y <- pinrs300$rsarea
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))


#rsweight v h
#stem mass = 0.0328 + 0.0025*h R2=0.187 p=0.12
x <- pinrs300$height
y <- pinrs300$rsweight
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))

#rsweight v d
#sm = 0.033 + 0.108*d R2=0.26, p=0.07
x <- pinrs300$diameter
y <- pinrs300$rsweight
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))



ggscatter(pinmaxdep, x="height", y="stemarea", color = "hclass")
ggscatter(pinmaxdep, x="height", y="stemmass", color = "hclass")
ggscatter(pinmaxdep, x="diameter", y="stemarea", color = "hclass")
ggscatter(pinmaxdep, x="diameter", y="stemmass", color = "hclass")




classpin1$hclass <- factor(classpin1$hclass, levels = c("<150", "150-300", ">300"))
classpin1$rsarea <- (classpin1$rootarea)/(classpin1$leafarea+classpin1$stemarea)
classpin1$rsweight <- (classpin1$rootmass)/(classpin1$leafmass+classpin1$stemmass)



fig7c <- ggplot() +
  geom_boxplot(pinon, mapping= aes(x=hclass, y=rsweight), fill="orange3", outlier.shape = NA, width=0.6) + 
  stat_boxplot(pinon, mapping= aes(x=hclass, y=rsweight), geom = "errorbar", width=0.35) +
  scale_y_continuous(expand = c(0,0), breaks = c(0,0.5,1,1.5,2), limits = c(0,2),
                     labels = c("0.0","","1.0","","2.0")) +
  labs(y="R:S by weight", title = "Piñon pine") +
  theme_classic() +
  theme(legend.title = element_blank(), legend.position = "none", 
        axis.text.x = element_text(size = 26), axis.title.x = element_blank(), legend.spacing.y = unit(1, "mm"), 
        plot.title = element_text(size = 28), axis.text = element_text(size=26), axis.title = element_text(size=26)) +
  geom_rect(data=xy, aes(xmin=0.5, xmax=1.5, ymin=0.615, ymax=0.627), color="gray80", fill="gray80", size=1) +
  geom_rect(data=xy, aes(xmin=1.5, xmax=2.5, ymin=0.002, ymax=0.525), color="gray80", fill="gray80", size=1) +
  geom_rect(data=xy, aes(xmin=2.5, xmax=3.5, ymin=0.002, ymax=0.525), color="gray80", fill="gray80", size=1) +
  geom_boxplot(pinon, mapping= aes(x=hclass, y=rsweight), fill="orange3", outlier.shape = NA, width=0.6) + 
  stat_boxplot(pinon, mapping= aes(x=hclass, y=rsweight), geom = "errorbar", width=0.35) +
  annotate(geom = "text", x=1, y=1.55, label="be", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=2, y=1.2, label="cd", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=3, y=0.7, label="n=4", fontface='italic', color="black", size=10)

t.test(pinyon150m$rsweight, pinyon150300m$rsweight, alternative = "greater") #p=8e-05
t.test(pinyon150300m$rsweight, pinyon300m$rsweight, alternative = "greater") #p=0.02
t.test(pinyon150m$rsweight, pinyon300m$rsweight, alternative = "greater") #p=3.1e-05








pinon$rsarea <- pinon$rootarea/(pinon$leafarea+pinon$stemarea)


fig7d <- ggplot() +
  geom_boxplot(pinon, mapping= aes(x=hclass, y=rsarea), fill="orange3", outlier.shape = NA, width=0.6) + 
  stat_boxplot(pinon, mapping= aes(x=hclass, y=rsarea), geom = "errorbar", width=0.35) +
  scale_y_continuous(expand = c(0,0), breaks = c(0,0.25,0.50), limits = c(0,0.5),
                     labels = c("0.0","","0.5")) +
  labs(y="R:S by area", title = "Piñon pine") +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_blank(), legend.position = "none",
        legend.background = element_blank(), plot.title = element_text(size = 28),
        axis.text.x = element_text(size = 26), axis.title.x = element_blank(), legend.spacing.y = unit(1, "mm"),
        axis.text = element_text(size=26), axis.title = element_text(size=26)) +
  geom_rect(data=xy, aes(xmin=0.5, xmax=1.5, ymin=0.128, ymax=0.2), color="gray80", fill="gray80", size=1) +
  geom_rect(data=xy, aes(xmin=1.5, xmax=2.5, ymin=0.039, ymax=0.0737), color="gray80", fill="gray80", size=1) +
  geom_rect(data=xy, aes(xmin=2.5, xmax=3.5, ymin=0.039, ymax=0.0737), color="gray80", fill="gray80", size=1) +
  geom_boxplot(pinon, mapping= aes(x=hclass, y=rsarea), fill="orange3", outlier.shape = NA, width=0.6) + 
  stat_boxplot(pinon, mapping= aes(x=hclass, y=rsarea), geom = "errorbar", width=0.35) +
  annotate(geom = "text", x=1, y=0.4, label="ab", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=2, y=0.26, label="c", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=3, y=0.15, label="n=4", fontface='italic', color="black", size=10)

t.test(pinyon150m$rsarea, pinyon150300m$rsarea, alternative = "greater") #p=0.0005
t.test(pinyon150300m$rsarea, pinyon300m$rsarea, alternative = "greater") #p=0.003
t.test(pinyon150m$rsarea, pinyon300m$rsarea, alternative = "greater") #p=1.6e-06







#ponderosa###############################################################################################
ponrs <- subset(excavations, Species == "Ponderosa")
ponrs150 <- subset(ponrs, hclass == "<150")
ponrs150300 <- subset(ponrs, hclass == "150-300")
ponrs300 <- subset(ponrs, hclass == ">300")

#rsarea v h
#stem area = 386.658 + 1.581*h R2=-0.09 p=0.66
x <- ponrs150$height
y <- ponrs150$rsarea
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))

#rsarea v d
#sa = 148.9 + 168.4*d R2=0.09, p=0.2
x <- ponrs150$diameter
y <- ponrs150$rsarea
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))


#rsweight v h
#stem mass = 0.0328 + 0.0025*h R2=0.187 p=0.12
x <- ponrs150$height
y <- ponrs150$rsweight
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))

#rsweight v d
#sm = 0.033 + 0.108*d R2=0.26, p=0.07
x <- ponrs150$diameter
y <- ponrs150$rsweight
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))





#rsarea v h
#stem area = 386.658 + 1.581*h R2=-0.09 p=0.66
x <- ponrs150300$height
y <- ponrs150300$rsarea
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))

#rsarea v d
#sa = 148.9 + 168.4*d R2=0.09, p=0.2
x <- ponrs150300$diameter
y <- ponrs150300$rsarea
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))


#rsweight v h
#stem mass = 0.0328 + 0.0025*h R2=0.187 p=0.12
x <- ponrs150300$height
y <- ponrs150300$rsweight
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))

#rsweight v d
#sm = 0.033 + 0.108*d R2=0.26, p=0.07
x <- ponrs150300$diameter
y <- ponrs150300$rsweight
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))





#rsarea v h
#stem area = 386.658 + 1.581*h R2=-0.09 p=0.66
x <- ponrs300$height
y <- ponrs300$rsarea
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))

#rsarea v d
#sa = 148.9 + 168.4*d R2=0.09, p=0.2
x <- ponrs300$diameter
y <- ponrs300$rsarea
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))


#rsweight v h
#stem mass = 0.0328 + 0.0025*h R2=0.187 p=0.12
x <- ponrs300$height
y <- ponrs300$rsweight
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))

#rsweight v d
#sm = 0.033 + 0.108*d R2=0.26, p=0.07
x <- ponrs300$diameter
y <- ponrs300$rsweight
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))



ggscatter(ponmaxdep, x="height", y="stemarea", color = "hclass")
ggscatter(ponmaxdep, x="height", y="stemmass", color = "hclass")
ggscatter(ponmaxdep, x="diameter", y="stemarea", color = "hclass")
ggscatter(ponmaxdep, x="diameter", y="stemmass", color = "hclass")




fig7e <- ggplot() +
  geom_boxplot(ponderosa, mapping= aes(x=hclass, y=rsweight), fill="orchid3", outlier.shape = NA, width=0.6) +
  stat_boxplot(ponderosa, mapping= aes(x=hclass, y=rsweight), geom = "errorbar", width=0.35) +
  scale_y_continuous(expand = c(0,0), breaks = c(0,0.5,1,1.5,2), limits = c(0,2),
                     labels = c("0.0","","1.0","","2.0")) +
  labs(y="R:S by weight", title = "Ponderosa pine") +
  theme_classic() +
  theme(legend.title = element_blank(), legend.position = "none",
        axis.text.x = element_text(size = 26), axis.title.x = element_blank(), legend.spacing.y = unit(1, "mm"), 
        plot.title = element_text(size = 28), axis.text = element_text(size=26), axis.title = element_text(size=26)) +
  geom_rect(data=xy, aes(xmin=0.5, xmax=1.5, ymin=0.455, ymax=0.462), color="gray80", fill="gray80", size=1) +
  geom_rect(data=xy, aes(xmin=1.5, xmax=2.5, ymin=0.357, ymax=0.358), color="gray80", fill="gray80", size=1) +
  geom_rect(data=xy, aes(xmin=2.5, xmax=3.5, ymin=0.332, ymax=0.9), color="gray80", fill="gray80", size=1) +
  geom_boxplot(ponderosa, mapping= aes(x=hclass, y=rsweight), fill="orchid3", outlier.shape = NA, width=0.6) +
  stat_boxplot(ponderosa, mapping= aes(x=hclass, y=rsweight), geom = "errorbar", width=0.35) +
  annotate(geom = "text", x=1, y=1.25, label="cf", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=2, y=0.75, label="d", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=3, y=0.67, label="d", fontface='italic', color="black", size=10)

t.test(ponderosa150m$rsweight, ponderosa150300m$rsweight, alternative = "greater") #p=0.001
t.test(ponderosa150300m$rsweight, ponderosa300m$rsweight, alternative = "greater") 
t.test(ponderosa150m$rsweight, ponderosa300m$rsweight, alternative = "greater") #p=0.002




ponderosa$rsarea <- ponderosa$rootarea/(ponderosa$leafarea+ponderosa$stemarea)




fig7f <- ggplot() +
  geom_boxplot(ponderosa, mapping= aes(x=hclass, y=rsarea), fill="orchid3", outlier.shape = NA, width=0.6) +
  stat_boxplot(ponderosa, mapping= aes(x=hclass, y=rsarea), geom = "errorbar", width=0.35) +
  scale_y_continuous(expand = c(0,0), breaks = c(0,0.25,0.5,0.75,1), limits = c(0,1),
                     labels = c("0.0","","0.5","","1.0")) +
  labs(y="R:S by area", title = "Ponderosa pine") +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_blank(), legend.position = "none",
        legend.background = element_blank(), plot.title = element_text(size = 28),
        axis.text.x = element_text(size = 26), axis.title.x = element_blank(), legend.spacing.y = unit(1, "mm"),
        axis.text = element_text(size=26), axis.title = element_text(size=26)) +
  geom_rect(data=xy, aes(xmin=0.5, xmax=1.5, ymin=0.149, ymax=0.15), color="gray80", fill="gray80", size=1) +
  geom_rect(data=xy, aes(xmin=1.5, xmax=2.5, ymin=0.076, ymax=0.151), color="gray80", fill="gray80", size=1) +
  geom_rect(data=xy, aes(xmin=2.5, xmax=3.5, ymin=0.0389, ymax=0.192), color="gray80", fill="gray80", size=1) +
  geom_boxplot(ponderosa, mapping= aes(x=hclass, y=rsarea), fill="orchid3", outlier.shape = NA, width=0.6) +
  stat_boxplot(ponderosa, mapping= aes(x=hclass, y=rsarea), geom = "errorbar", width=0.35) +
  annotate(geom = "text", x=1, y=0.55, label="ab", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=2, y=0.3, label="c", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=3, y=0.22, label="c", fontface='italic', color="black", size=10)

t.test(ponderosa150m$rsarea, ponderosa150300m$rsarea, alternative = "greater") #p=0.002
t.test(ponderosa150300m$rsarea, ponderosa300m$rsarea, alternative = "greater") 
t.test(ponderosa150m$rsarea, ponderosa300m$rsarea, alternative = "greater") #p=5.1e-05






#fir###############################################################################################
frs <- subset(excavations, Species == "Fir")
frs150 <- subset(frs, hclass == "<150")
frs150300 <- subset(frs, hclass == "150-300")
frs300 <- subset(frs, hclass != "<150")

#rsarea v h
#stem area = 386.658 + 1.581*h R2=-0.09 p=0.66
x <- frs150$height
y <- frs150$rsarea
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))

#rsarea v d
#sa = 148.9 + 168.4*d R2=0.09, p=0.2
x <- frs150$diameter
y <- frs150$rsarea
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))


#rsweight v h
#stem mass = 0.0328 + 0.0025*h R2=0.187 p=0.12
x <- frs150$height
y <- frs150$rsweight
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))

#rsweight v d
#sm = 0.033 + 0.108*d R2=0.26, p=0.07
x <- frs150$diameter
y <- frs150$rsweight
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))





#rsarea v h
#stem area = 386.658 + 1.581*h R2=-0.09 p=0.66
x <- frs300$height
y <- frs300$rsarea
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))

#rsarea v d
#sa = 148.9 + 168.4*d R2=0.09, p=0.2
x <- frs300$diameter
y <- frs300$rsarea
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))


#rsweight v h
#stem mass = 0.0328 + 0.0025*h R2=0.187 p=0.12
x <- frs300$height
y <- frs300$rsweight
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))

#rsweight v d
#sm = 0.033 + 0.108*d R2=0.26, p=0.07
x <- frs300$diameter
y <- frs300$rsweight
z <- lm(y ~ x)
summary(z)
sqrt(mean(z$residuals^2))




ggscatter(fmaxdep, x="height", y="stemarea", color = "hclass")
ggscatter(fmaxdep, x="height", y="stemmass", color = "hclass")
ggscatter(fmaxdep, x="diameter", y="stemarea", color = "hclass")
ggscatter(fmaxdep, x="diameter", y="stemmass", color = "hclass")



classfir1$hclass <- factor(classfir1$hclass, levels = c("<150", "150-300", ">300"))
classfir1$rsarea <- (classfir1$rootarea)/(classfir1$leafarea+classfir1$stemarea)
classfir1$rsweight <- (classfir1$rootmass)/(classfir1$leafmass+classfir1$stemmass)




fig7g <- ggplot() +
  geom_boxplot(fir, mapping= aes(x=hclass, y=rsweight), fill="palegreen4", outlier.shape = NA, width=0.6) +
  stat_boxplot(fir, mapping= aes(x=hclass, y=rsweight), geom = "errorbar", width=0.35) +
  scale_y_continuous(expand = c(0,0), breaks = c(0,0.5,1,1.5,2), limits = c(0,2),
                     labels = c("0.0","","1.0","","2.0")) +
  labs(x="Height class (mm)", y="R:S by weight", title = "Fir") +
  theme_classic() +
  theme(legend.title = element_blank(), legend.position = "none",
        axis.text.x = element_text(size = 26), axis.title.x = element_text(size = 26), 
        legend.spacing.y = unit(1, "mm"), 
        plot.title = element_text(size = 28), axis.text = element_text(size=26), axis.title = element_text(size=26)) +
  geom_rect(data=xy, aes(xmin=0.5, xmax=1.5, ymin=0.852, ymax=0.878), color="gray80", fill="gray80", size=1) +
  geom_rect(data=xy, aes(xmin=1.5, xmax=2.5, ymin=0.19, ymax=0.927), color="gray80", fill="gray80", size=1) +
  geom_boxplot(fir, mapping= aes(x=hclass, y=rsweight), fill="palegreen4", outlier.shape = NA, width=0.6) +
  stat_boxplot(fir, mapping= aes(x=hclass, y=rsweight), geom = "errorbar", width=0.35) +
  annotate(geom = "text", x=1, y=1.92, label="a", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=2, y=1.6, label="ef", fontface='italic', color="black", size=10) 

t.test(fir150m$rsweight, fir150300m$rsweight, alternative = "greater") #p=0.03






fir$rsarea <- fir$rootarea/(fir$leafarea+fir$stemarea)


fig7h <- ggplot() +
  geom_boxplot(fir, mapping= aes(x=hclass, y=rsarea), fill="palegreen4", outlier.shape = NA, width=0.6) +
  stat_boxplot(fir, mapping= aes(x=hclass, y=rsarea), geom = "errorbar", width=0.35) +
  scale_y_continuous(expand = c(0,0), breaks = c(0,0.25,0.5,0.75,1), limits = c(0,1),
                     labels = c("0.0","","0.5","","1.0")) +
  labs(x="Height class (mm)", y="R:S by area", title = "Fir") +
  theme_classic() +
  theme(legend.title = element_blank(), legend.position="none", plot.title = element_text(size=28),
        axis.text.x = element_text(size = 26), axis.title.x = element_text(size = 26), 
        axis.text = element_text(size=26), axis.title = element_text(size=26)) +
  geom_rect(data=xy, aes(xmin=0.5, xmax=1.5, ymin=0.147, ymax=0.27), color="gray80", fill="gray80", size=1) +
  geom_rect(data=xy, aes(xmin=1.5, xmax=2.5, ymin=0.0729, ymax=0.138), color="gray80", fill="gray80", size=1) +
  geom_boxplot(fir, mapping= aes(x=hclass, y=rsarea), fill="palegreen4", outlier.shape = NA, width=0.6) +
  stat_boxplot(fir, mapping= aes(x=hclass, y=rsarea), geom = "errorbar", width=0.35) +
  annotate(geom = "text", x=1, y=0.62, label="ab", fontface='italic', color="black", size=10) +
  annotate(geom = "text", x=2, y=0.3, label="ac", fontface='italic', color="black", size=10)

t.test(fir150m$rsarea, fir150300m$rsarea, alternative = "greater")




plot_grid(fig7a, fig7b, fig7c, fig7d, fig7e, fig7f, fig7g, fig7h, nrow = 4,
          labels = c("a", "b", "c", "d", "e", "f", "g", "h"),
          label_x = 0.94, label_y = 0.95, label_size = 40)



t.test(juniper150m$rsweight, pinyon150m$rsweight, alternative = "greater") 
t.test(juniper150m$rsweight, ponderosa150m$rsweight, alternative = "greater") 
t.test(juniper150m$rsweight, fir150m$rsweight, alternative = "less") 
t.test(pinyon150m$rsweight, ponderosa150m$rsweight, alternative = "greater") #p<0.05
t.test(pinyon150m$rsweight, fir150m$rsweight, alternative = "less") #p<0.05, fir greater
t.test(ponderosa150m$rsweight, fir150m$rsweight, alternative = "less") #p<0.05, fir greater


t.test(pinyon150300m$rsweight, ponderosa150300m$rsweight, alternative = "greater") 
t.test(pinyon150300m$rsweight, fir150300m$rsweight, alternative = "less") #p<0.05
t.test(ponderosa150300m$rsweight, fir150300m$rsweight, alternative = "less") #p<0.05


t.test(pinyon150300m$rsweight, juniper150m$rsweight, alternative = "less") 
t.test(pinyon150300m$rsweight, ponderosa150m$rsweight, alternative = "less") 
t.test(pinyon150300m$rsweight, fir150m$rsweight, alternative = "less") #p<0.05
t.test(ponderosa150300m$rsweight, juniper150m$rsweight, alternative = "less") #p<0.05
t.test(ponderosa150300m$rsweight, pinyon150m$rsweight, alternative = "less") #p<0.05
t.test(ponderosa150300m$rsweight, fir150m$rsweight, alternative = "less") #p<0.05
t.test(fir150300m$rsweight, juniper150m$rsweight, alternative = "less") 
t.test(fir150300m$rsweight, pinyon150m$rsweight, alternative = "less") 
t.test(fir150300m$rsweight, ponderosa150m$rsweight, alternative = "greater") 


t.test(ponderosa300m$rsweight, juniper150m$rsweight, alternative = "less") #p<0.05
t.test(ponderosa300m$rsweight, pinyon150m$rsweight, alternative = "less") #p<0.05
t.test(ponderosa300m$rsweight, fir150m$rsweight, alternative = "less") #p<0.05


t.test(ponderosa300m$rsweight, pinyon150300m$rsweight, alternative = "less") 
t.test(ponderosa300m$rsweight, fir150300m$rsweight, alternative = "less") #p<0.05






t.test(juniper150m$rsarea, pinyon150m$rsarea, alternative = "greater") 
t.test(juniper150m$rsarea, ponderosa150m$rsarea, alternative = "greater") 
t.test(juniper150m$rsarea, fir150m$rsarea, alternative = "greater") 
t.test(pinyon150m$rsarea, ponderosa150m$rsarea, alternative = "less") 
t.test(pinyon150m$rsarea, fir150m$rsarea, alternative = "less") 
t.test(ponderosa150m$rsarea, fir150m$rsarea, alternative = "less") 


t.test(pinyon150300m$rsarea, ponderosa150300m$rsarea, alternative = "less") 
t.test(pinyon150300m$rsarea, fir150300m$rsarea, alternative = "less") 
t.test(ponderosa150300m$rsarea, fir150300m$rsarea, alternative = "less") 


t.test(pinyon150300m$rsarea, juniper150m$rsarea, alternative = "less") #p<0.05
t.test(pinyon150300m$rsarea, ponderosa150m$rsarea, alternative = "less") #p<0.05
t.test(pinyon150300m$rsarea, fir150m$rsarea, alternative = "less") #p<0.05
t.test(ponderosa150300m$rsarea, juniper150m$rsarea, alternative = "less") #p<0.05
t.test(ponderosa150300m$rsarea, pinyon150m$rsarea, alternative = "less") #p<0.05
t.test(ponderosa150300m$rsarea, fir150m$rsarea, alternative = "less") #p<0.05
t.test(fir150300m$rsarea, juniper150m$rsarea, alternative = "less") #p<0.05
t.test(fir150300m$rsarea, pinyon150m$rsarea, alternative = "less") 
t.test(fir150300m$rsarea, ponderosa150m$rsarea, alternative = "less") 


t.test(ponderosa300m$rsarea, juniper150m$rsarea, alternative = "less") #p<0.05
t.test(ponderosa300m$rsarea, pinyon150m$rsarea, alternative = "less") #p<0.05
t.test(ponderosa300m$rsarea, fir150m$rsarea, alternative = "less") #p<0.05


t.test(ponderosa300m$rsarea, pinyon150300m$rsarea, alternative = "less") 
t.test(ponderosa300m$rsarea, fir150300m$rsarea, alternative = "less") 







#calculation of height ratio####################################################################
classjun150 <- subset(classjun, hclass == "<150")
classjun150$meanh <- mean(classjun150$height)
classjun150$meand <- mean(classjun150$diameter)
classjun150$hratio <- classjun150$height/classjun150$meanh
classjun150$hratio <- classjun150$hratio
classjun150$dratio <- classjun150$diameter/classjun150$meand

classjun150300 <- subset(classjun, hclass == "150-300")
classjun150300$meanh <- mean(classjun150300$height)
classjun150300$meand <- mean(classjun150300$diameter)
classjun150300$hratio <- classjun150300$height/classjun150300$meanh
classjun150300$hratio <- classjun150300$hratio
classjun150300$dratio <- classjun150300$diameter/classjun150300$meand

classjun300 <- subset(classjun, hclass == "300+")
classjun300$meanh <- mean(classjun300$height)
classjun300$meand <- mean(classjun300$diameter)
classjun300$hratio <- classjun300$height/classjun300$meanh
classjun300$hratio <- classjun300$hratio
classjun300$dratio <- classjun300$diameter/classjun300$meand

classjun <- rbind(classjun150,classjun150300,classjun300)
classjun <- classjun[order(classjun[,'core']), ]
write.csv(classjun, file = "juniper/classjun.csv")
classjun <- read.csv(file = "juniper/classjun.csv")



classpin150 <- subset(classpin, hclass == "<150")
classpin150$meanh <- mean(classpin150$height)
classpin150$meand <- mean(classpin150$diameter)
classpin150$hratio <- classpin150$height/classpin150$meanh
classpin150$hratio <- classpin150$hratio
classpin150$dratio <- classpin150$diameter/classpin150$meand

classpin150300 <- subset(classpin, hclass == "150-300")
classpin150300$meanh <- mean(classpin150300$height)
classpin150300$meand <- mean(classpin150300$diameter)
classpin150300$hratio <- classpin150300$height/classpin150300$meanh
classpin150300$hratio <- classpin150300$hratio
classpin150300$dratio <- classpin150300$diameter/classpin150300$meand

classpin300 <- subset(classpin, hclass == "300+")
classpin300$meanh <- mean(classpin300$height)
classpin300$meand <- mean(classpin300$diameter)
classpin300$hratio <- classpin300$height/classpin300$meanh
classpin300$hratio <- log(classpin300$hratio)
classpin300$dratio <- classpin300$diameter/classpin300$meand

classpin <- rbind(classpin150,classpin150300,classpin300)
classpin <- classpin[order(classpin[,'core']), ]
write.csv(classpin, file = "pinon/classpin.csv")
classpin <- read.csv(file = "pinon/classpin.csv")




classpon150 <- subset(classpon, hclass == "<150")
classpon150$meanh <- mean(classpon150$height)
classpon150$meand <- mean(classpon150$diameter)
classpon150$hratio <- classpon150$height/classpon150$meanh
classpon150$hratio <- log(classpon150$hratio)
classpon150$dratio <- classpon150$diameter/classpon150$meand

classpon150300 <- subset(classpon, hclass == "150-300")
classpon150300$meanh <- mean(classpon150300$height)
classpon150300$meand <- mean(classpon150300$diameter)
classpon150300$hratio <- classpon150300$height/classpon150300$meanh
classpon150300$hratio <- log(classpon150300$hratio)
classpon150300$dratio <- classpon150300$diameter/classpon150300$meand

classpon300 <- subset(classpon, hclass == "300+")
classpon300$meanh <- mean(classpon300$height)
classpon300$meand <- mean(classpon300$diameter)
classpon300$hratio <- classpon300$height/classpon300$meanh
classpon300$hratio <- log(classpon300$hratio)
classpon300$dratio <- classpon300$diameter/classpon300$meand

classpon <- rbind(classpon150,classpon150300,classpon300)
classpon <- classpon[order(classpon[,'core']), ]
write.csv(classpon, file = "ponderosa/classpon.csv")
classpon <- read.csv(file = "ponderosa/classpon.csv")




classfir150 <- subset(classfir, hclass == "<150")
classfir150$meanh <- mean(classfir150$height)
classfir150$meand <- mean(classfir150$diameter)
classfir150$hratio <- classfir150$height/classfir150$meanh
classfir150$hratio <- log(classfir150$hratio)
classfir150$dratio <- classfir150$diameter/classfir150$meand

classfir150300 <- subset(classfir, hclass == "150-300")
classfir150300$meanh <- mean(classfir150300$height)
classfir150300$meand <- mean(classfir150300$diameter)
classfir150300$hratio <- classfir150300$height/classfir150300$meanh
classfir150300$hratio <- log(classfir150300$hratio)
classfir150300$dratio <- classfir150300$diameter/classfir150300$meand

classfir300 <- subset(classfir, hclass == "300+")
classfir300$meanh <- mean(classfir300$height)
classfir300$meand <- mean(classfir300$diameter)
classfir300$hratio <- classfir300$height/classfir300$meanh
classfir300$hratio <- log(classfir300$hratio)
classfir300$dratio <- classfir300$diameter/classfir300$meand

classfir <- rbind(classfir150,classfir150300,classfir300)
classfir <- classfir[order(classfir[,'core']), ]
write.csv(classfir, file = "fir/classfir.csv")
classfir <- read.csv(file = "fir/classfir.csv")


allspec <- rbind(classjun, classpin, classpon, classfir)
allspec <- allspec[with(allspec, order(core, plot)), ]
allspec$core <- factor(allspec$core, levels = c("AZFG", "COMT", "NMSV", 'NMVC', "NVSH", "NVSP"))
write.csv(allspec, file = "allspec.csv")
allspec <- read.csv(file = "allspec.csv")



allspeccomt <- subset(allspec, core == "COMT")
allspecnmsv <- subset(allspec, core == "NMSV")
allspecnvsp <- subset(allspec, core == "NVSP")


anovah <- aov(hratio ~ plot, data = allspec)
summary(anovah)





meanhratio <- with(allspec, by(hratio,plot,mean))
meandratio <- with(allspec, by(dratio,plot,mean))
elevation <- with(allspec, by(elevation,plot,mean))

df <- data.frame(cbind(elevation, meanhratio, meandratio))
write.csv(df, file = "df.csv")
df <- read.csv(file = "df.csv")
df$core <- c("AZFG", "AZFG", "AZFG", "COMT", "COMT", "COMT", "COMT", "COMT", "NMSV", "NMSV", "NMSV", "NMSV", "NMVC",
             "NMVC", "NVSH", "NVSH", "NVSH", "NVSP", "NVSP", "NVSP", "NVSP", "NVSP", "NVSP", "NVSP")
df$sig <- factor(df$sig, levels = c("Low", "Average", "High"))
dfh <- df[-c(3,4,7,9,10,12,13,14,15,16,17,18,20,22,23,24), ] 
dfd <- df[-c(1,3,7,10,11,14,15,16,17,18,20,22), ]


#scatterplots of abiotic and biotic site characteristics to hratio#######################################
#abiotic site characteristics####################################################################
#elevation - no correlations to plots or sites
#run stats on all plots to see if sig relationships
#highlight the ones that are sig different from each other by using diff colors on one end and the other
x <- df$elevation
y <- df$meanhratio
z <- lm(y ~ x)
summary(z)

x <- df$elevation
y <- df$meandratio
z <- lm(y ~ x)
summary(z)

#color groups differently
ggplot() +
  geom_point(dfh, mapping= aes(x=elevation, y=meanhratio, color=core), size=3) +
  geom_hline(yintercept = 1.05) +
  geom_segment(aes(x = 2800, xend=2800, y = 1.06, yend = 1.16), size=1.2,
               arrow = arrow(length = unit(0.5, "cm"))) +
  geom_segment(aes(x = 2800, xend=2800, y = 1.04, yend = 0.94), size=1.2,
               arrow = arrow(length = unit(0.5, "cm"))) +
  annotate(geom = "text", x=2700, y=1.1, label="a, ab", fontface='italic', color="black", size=5.5) +
  annotate(geom = "text", x=2700, y=1, label="cde, de, e", fontface='italic', color="black", size=5.5) +
  annotate(geom = "text", x=2350, y=1.12, label="abc", fontface='italic', color="black", size=5.5) +
  guides(colour = guide_legend(nrow = 1)) +
  labs(x="Elevation (m)", y="Height ratio") +
  scale_x_continuous(expand = c(0,0), breaks = c(2000,2200,2400,2600,2800), limits = c(1940,2850)) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size=18), legend.position = "bottom",
        legend.background = element_blank(), legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 18),
        axis.text = element_text(size=16), axis.title = element_text(size = 16))

ggplot() +
  geom_point(dfd, mapping= aes(x=elevation, y=meandratio, color=core), size=3) +
  guides(colour = guide_legend(nrow = 1)) +
  labs(x="Elevation (m)", y="Diameter ratio") +
  scale_x_continuous(expand = c(0,0), breaks = c(2000,2200,2400,2600,2800), limits = c(1940,2850)) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size=18), legend.position = "bottom",
        legend.background = element_blank(), legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 18),
        axis.text = element_text(size=16), axis.title = element_text(size = 16))




#aspect
x <- df$aspect
y <- df$meanhratio
z <- lm(y ~ x)
summary(z)

x <- df$aspect
y <- df$meandratio
z <- lm(y ~ x)
summary(z)


ggplot() +
  geom_point(dfh, mapping= aes(x=aspect, y=meanhratio, color=core), size=3) +
  guides(color = guide_legend(nrow = 1)) +
  labs(x="Aspect (deg)", y="Height ratio") +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size=18), legend.position = "bottom",
        legend.background = element_blank(), legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 18),
        axis.text = element_text(size=16), axis.title = element_text(size = 16), axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 16))

ggplot() +
  geom_point(dfd, mapping= aes(x=aspect, y=meandratio, color=core), size=3) +
  guides(color = guide_legend(nrow = 1)) +
  labs(x="Aspect (deg)", y="Diameter ratio") +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size=18), legend.position = "bottom",
        legend.background = element_blank(), legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 18),
        axis.text = element_text(size=16), axis.title = element_text(size = 16), axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 16))


#slope
x <- df$slope
y <- df$meanhratio
z <- lm(y ~ x)
summary(z)

x <- df$slope
y <- df$meandratio
z <- lm(y ~ x)
summary(z)


ggplot() +
  geom_point(dfh, mapping= aes(x=slope, y=meanhratio, color=core), size=3) +
  guides(color = guide_legend(nrow = 1)) +
  labs(x="Slope (deg)", y="Height ratio") +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size=18), legend.position = "bottom",
        legend.background = element_blank(), legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 18),
        axis.text = element_text(size=16), axis.title = element_text(size = 16), axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 16))

ggplot() +
  geom_point(dfd, mapping= aes(x=slope, y=meandratio, color=core), size=3) +
  guides(color = guide_legend(nrow = 1)) +
  labs(x="Slope (deg)", y="Diameter ratio") +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size=18), legend.position = "bottom",
        legend.background = element_blank(), legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 18),
        axis.text = element_text(size=16), axis.title = element_text(size = 16), axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 16))



#percent sand
x <- df$sand
y <- df$meanhratio
z <- lm(y ~ x)
summary(z)

x <- df$sand
y <- df$meandratio
z <- lm(y ~ x)
summary(z)


ggplot() +
  geom_point(dfh, mapping= aes(x=sand, y=meanhratio, color=core), size=3) +
  guides(color = guide_legend(nrow = 1)) +
  scale_x_continuous(breaks = c(20,40,60,80), limits = c(20,80)) +
  labs(x="Percent sand", y="Height ratio") +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size=18), legend.position = "bottom",
        legend.background = element_blank(), legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 18),
        axis.text = element_text(size=16), axis.title = element_text(size = 16), axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 16))

ggplot() +
  geom_point(dfd, mapping= aes(x=sand, y=meandratio, color=core), size=3) +
  guides(color = guide_legend(nrow = 1)) +
  labs(x="Percent sand", y="Diameter ratio") +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size=18), legend.position = "bottom",
        legend.background = element_blank(), legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 18),
        axis.text = element_text(size=16), axis.title = element_text(size = 16), axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 16))



#percent clay
x <- df$clay
y <- df$meanhratio
z <- lm(y ~ x)
summary(z)

x <- df$clay
y <- df$meandratio
z <- lm(y ~ x)
summary(z)


ggplot() +
  geom_point(dfh, mapping= aes(x=clay, y=meanhratio, color=core), size=3) +
  guides(color= guide_legend(nrow = 1)) +
  labs(x="Percent clay", y="Height ratio") +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size=18), legend.position = "bottom",
        legend.background = element_blank(), legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 18),
        axis.text = element_text(size=16), axis.title = element_text(size = 16), axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 16))

ggplot() +
  geom_point(dfd, mapping= aes(x=clay, y=meandratio, color=core), size=3) +
  guides(color = guide_legend(nrow = 1)) +
  scale_x_continuous(breaks = c(0,10,20,30), limits = c(0,30)) +
  labs(x="Percent clay", y="Diameter ratio") +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size=18), legend.position = "bottom",
        legend.background = element_blank(), legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 18),
        axis.text = element_text(size=16), axis.title = element_text(size = 16), axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 16))



#percent om
x <- df$om
y <- df$meanhratio
z <- lm(y ~ x)
summary(z)

x <- df$om
y <- df$meandratio
z <- lm(y ~ x)
summary(z)


ggplot() +
  geom_point(dfh, mapping= aes(x=om, y=meanhratio, color=core), size=3) +
  guides(color = guide_legend(nrow = 1)) +
  labs(x="Percent OM", y="Height ratio") +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size=18), legend.position = "bottom",
        legend.background = element_blank(), legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 18),
        axis.text = element_text(size=16), axis.title = element_text(size = 16), axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 16))

ggplot() +
  geom_point(dfd, mapping= aes(x=om, y=meandratio, color=core), size=3) +
  guides(color = guide_legend(nrow = 1)) +
  labs(x="Percent OM", y="Diameter ratio") +
  scale_x_continuous(breaks = c(2,4,6,8,10), limits = c(2,10)) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size=18), legend.position = "bottom",
        legend.background = element_blank(), legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 18),
        axis.text = element_text(size=16), axis.title = element_text(size = 16), axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 16))



#climate########################################################################
#MAP - no correlation
x <- df$map
y <- df$meanhratio
z <- lm(y ~ x)
summary(z)

x <- df$map
y <- df$meandratio
z <- lm(y ~ x)
summary(z)


ggplot() +
  geom_point(dfh, mapping=aes(x=map, y=meanhratio, color=core), size=3) +
  guides(color = guide_legend(nrow = 1)) +
  labs(x="MAP (mm)", y="Height ratio") +
  scale_x_continuous(breaks = c(30,40,50,60), limits = c(30,60)) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size=18), legend.position = "bottom",
        legend.background = element_blank(), legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 18),
        axis.text = element_text(size=16), axis.title = element_text(size = 16), axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 16))

ggplot() +
  geom_point(dfd, mapping=aes(x=map, y=meandratio, color=core), size=3) +
  guides(color = guide_legend(nrow = 1)) +
  labs(x="MAP (mm)", y="Diameter ratio") +
  scale_x_continuous(breaks = c(30,40,50,60), limits = c(30,60)) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size=18), legend.position = "bottom",
        legend.background = element_blank(), legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 18),
        axis.text = element_text(size=16), axis.title = element_text(size = 16), axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 16))


#MAT - no correlation
x <- df$mat
y <- df$meanhratio
z <- lm(y ~ x)
summary(z)

x <- df$mat
y <- df$meandratio
z <- lm(y ~ x)
summary(z)


ggplot() +
  geom_point(dfh, mapping=aes(x=mat, y=meanhratio, color=core), size=3) +
  guides(color = guide_legend(nrow = 1)) +
  labs(x="MAT (deg C)", y="Height ratio") +
  scale_x_continuous(breaks = c(5,7.5,10,12.5), limits = c(5,12.5)) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size=18), legend.position = "bottom",
        legend.background = element_blank(), legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 18),
        axis.text = element_text(size=16), axis.title = element_text(size = 16), axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 16))

ggplot() +
  geom_point(dfd, mapping=aes(x=mat, y=meandratio, color=core), size=3) +
  guides(color = guide_legend(nrow = 1)) +
  labs(x="MAT (deg C)", y="Diameter ratio") +
  scale_x_continuous(breaks = c(5,7.5,10,12.5), limits = c(5,12.5)) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size=18), legend.position = "bottom",
        legend.background = element_blank(), legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 18),
        axis.text = element_text(size=16), axis.title = element_text(size = 16), axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 16))





#biotic site characteristics###################################################################
#basal area
x <- df$ba
y <- df$meanhratio
z <- lm(y ~ x)
summary(z)

x <- df$ba
y <- df$meandratio
z <- lm(y ~ x)
summary(z)


ggplot() +
  geom_point(dfh, mapping=aes(x=basal, y=meanhratio, color=core), size=3) +
  guides(color = guide_legend(nrow = 1)) +
  labs(x="Basal area (m2 ha-1)", y="Height ratio") +
  scale_x_continuous(breaks = c(0,10,20,30,40), limits = c(0,40)) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size=18), legend.position = "bottom",
        legend.background = element_blank(), legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 18),
        axis.text = element_text(size=16), axis.title = element_text(size = 16), axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 16))

ggplot() +
  geom_point(dfd, mapping=aes(x=basal, y=meandratio, color=core), size=3) +
  geom_smooth(data=dfd, mapping = aes(x=basal, y=meandratio), color="black", size=1.2, 
              method = "lm", 
              formula = y~(x), se=FALSE) +
  guides(color = guide_legend(nrow = 1)) +
  labs(x="Basal area (m2 ha-1)", y="Diameter ratio") +
  scale_x_continuous(breaks = c(10,20,30,40), limits = c(10,40)) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size=18), legend.position = "bottom",
        legend.background = element_blank(), legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 18),
        axis.text = element_text(size=16), axis.title = element_text(size = 16), axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 16))



#canopy cover
x <- df$canopy
y <- df$meanhratio
z <- lm(y ~ x)
summary(z)

x <- df$canopy
y <- df$meandratio
z <- lm(y ~ x)
summary(z)


ggplot() +
  geom_point(dfh, mapping=aes(x=canopy, y=meanhratio, color=core), size=3) +
  guides(color = guide_legend(nrow = 1)) +
  labs(x="Canopy cover (%)", y="Height ratio") +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size=18), legend.position = "bottom",
        legend.background = element_blank(), legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 18),
        axis.text = element_text(size=16), axis.title = element_text(size = 16), axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 16))

ggplot() +
  geom_point(dfd, mapping=aes(x=canopy, y=meandratio, color=core), size=3) +
  guides(color = guide_legend(nrow = 1)) +
  labs(x="Canopy cover (%)", y="Diameter ratio") +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size=18), legend.position = "bottom",
        legend.background = element_blank(), legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 18),
        axis.text = element_text(size=16), axis.title = element_text(size = 16), axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 16))



#number juveniles - no correlation
x <- df$juveniles
y <- df$meanhratio
z <- lm(y ~ x)
summary(z)

x <- df$juveniles
y <- df$meandratio
z <- lm(y ~ x)
summary(z)


ggplot() +
  geom_point(dfh, mapping=aes(x=juveniles, y=meanhratio, color=core), size=3) +
  guides(color = guide_legend(nrow = 1)) +
  labs(x="Number of juveniles", y="Height ratio") +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size=18), legend.position = "bottom",
        legend.background = element_blank(), legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 18),
        axis.text = element_text(size=16), axis.title = element_text(size = 16), axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 16))

ggplot() +
  geom_point(dfd, mapping=aes(x=juveniles, y=meandratio, color=core), size=3) +
  guides(color = guide_legend(nrow = 1)) +
  labs(x="Number of juveniles", y="Diameter ratio") +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size=18), legend.position = "bottom",
        legend.background = element_blank(), legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 18),
        axis.text = element_text(size=16), axis.title = element_text(size = 16), axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 16))



#juvenile age based on excavated samples
x <- df$age
y <- df$meanhratio
z <- lm(y ~ x)
summary(z)

x <- df$age
y <- df$meandratio
z <- lm(y ~ x)
summary(z)


ggplot() +
  geom_point(dfh, mapping=aes(x=age, y=meanhratio, color=core), size=3) +
  
  guides(color = guide_legend(nrow = 1)) +
  labs(x="Average age (yr)", y="Height ratio") +
  scale_x_continuous(breaks = c(2.5,5,7.5,10,12.5), limits = c(2.5,13)) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size=18), legend.position = "bottom",
        legend.background = element_blank(), legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 18),
        axis.text = element_text(size=16), axis.title = element_text(size = 16), axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 16))

ggplot() +
  geom_point(dfd, mapping=aes(x=age, y=meandratio, color=core), size=3) +
  guides(color = guide_legend(nrow = 1)) +
  labs(x="Average age (yr)", y="Diameter ratio") +
  scale_x_continuous(breaks = c(2.5,5,7.5,10,12.5), limits = c(2.5,12.5)) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size=18), legend.position = "bottom",
        legend.background = element_blank(), legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 18),
        axis.text = element_text(size=16), axis.title = element_text(size = 16), axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 16))



#number adults
x <- df$adults
y <- df$meanhratio
z <- lm(y ~ x)
summary(z)

x <- df$adults
y <- df$meandratio
z <- lm(y ~ x)
summary(z)


ggplot() +
  geom_point(dfh, mapping=aes(x=adults, y=meanhratio, color=core), size=3) +
  guides(color = guide_legend(nrow = 1)) +
  labs(x="Number of adults", y="Height ratio") +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size=18), legend.position = "bottom",
        legend.background = element_blank(), legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 18),
        axis.text = element_text(size=16), axis.title = element_text(size = 16), axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 16))

ggplot() +
  geom_point(dfd, mapping=aes(x=adults, y=meandratio, color=core), size=3) +
  guides(color = guide_legend(nrow = 1)) +
  labs(x="Number of adults", y="Diameter ratio") +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size=18), legend.position = "bottom",
        legend.background = element_blank(), legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 18),
        axis.text = element_text(size=16), axis.title = element_text(size = 16), axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 16))



#herb cover
x <- df$herb
y <- df$meandratio
z <- lm(y ~ x)
summary(z)


ggplot() +
  geom_point(dfd, mapping=aes(x=herb, y=meandratio, color=core), size=3) +
  guides(color = guide_legend(nrow = 1)) +
  labs(x="Herbaceous cover (%)", y="Diameter ratio") +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size=18), legend.position = "bottom",
        legend.background = element_blank(), legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 18),
        axis.text = element_text(size=16), axis.title = element_text(size = 16), axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 16))


#litter cover
x <- df$littercov
y <- df$meanhratio
z <- lm(y ~ x)
summary(z)

x <- df$littercov
y <- df$meandratio
z <- lm(y ~ x)
summary(z)


ggplot() +
  geom_point(dfh, mapping=aes(x=littercov, y=meanhratio, color=core), size=3) +
  guides(color = guide_legend(nrow = 1)) +
  labs(x="Litter cover (%)", y="Height ratio") +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size=18), legend.position = "bottom",
        legend.background = element_blank(), legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 18),
        axis.text = element_text(size=16), axis.title = element_text(size = 16), axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 16))

ggplot() +
  geom_point(dfd, mapping=aes(x=littercov, y=meandratio, color=core), size=3) +
  guides(color = guide_legend(nrow = 1)) +
  labs(x="Litter cover (%)", y="Diameter ratio") +
  scale_x_continuous(breaks = c(0,25,50,75,100), limits = c(0,100)) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size=18), legend.position = "bottom",
        legend.background = element_blank(), legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 18),
        axis.text = element_text(size=16), axis.title = element_text(size = 16), axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 16))



#litter depth
x <- df$litterdep
y <- df$meanhratio
z <- lm(y ~ x)
summary(z)

x <- df$litterdep
y <- df$meandratio
z <- lm(y ~ x)
summary(z)


ggplot() +
  geom_point(dfh, mapping=aes(x=litterdep, y=meanhratio, color=core), size=3) +
  guides(color = guide_legend(nrow = 1)) +
  labs(x="Litter depth (cm)", y="Height ratio") +
  scale_x_continuous(breaks = c(0,1,2,3,4,5), limits = c(0,5)) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size=18), legend.position = "bottom",
        legend.background = element_blank(), legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 18),
        axis.text = element_text(size=16), axis.title = element_text(size = 16), axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 16))

ggplot() +
  geom_point(dfd, mapping=aes(x=litterdep, y=meandratio, color=core), size=3) +
  geom_smooth(data=dfd, mapping = aes(x=litterdep, y=meandratio), color="black", size=1.2, 
              method = "lm", 
              formula = y~(x), se=FALSE) +
  guides(color = guide_legend(nrow = 1)) +
  labs(x="Litter depth (cm)", y="Diameter ratio") +
  scale_x_continuous(breaks = c(0,1,2,3,4,5), limits = c(0,5)) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size=18), legend.position = "bottom",
        legend.background = element_blank(), legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 18),
        axis.text = element_text(size=16), axis.title = element_text(size = 16), axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 16))








#figure 8########################################################################################

fig8a <- ggplot() +
  geom_hline(yintercept = 1, linetype=2) +
  geom_boxplot(allspec, mapping= aes(x=num, y=hratio, group=plot, fill=core), 
               width=1.5, outlier.shape = NA) +
  stat_boxplot(allspec, mapping= aes(x=num, y=hratio, group=plot), geom = "errorbar", width=0.75) +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x="Study sites", y="Height ratio") +
  scale_y_continuous(expand = c(0,0), breaks = c(0.0,1.0,2.0,3.0), limits = c(0.0,3.0)) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size=34), legend.position = "top",
        legend.margin = margin(0,0,0,0),
        legend.background = element_blank(), legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 28),
        axis.text = element_text(size=28), axis.title = element_text(size = 28), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), plot.margin = unit(c(0.5,0.5,0.25,0.95),"cm")) +
  annotate(geom = "text", x=1, y=2.2, label="ab", fontface='italic', color="tomato3", size=9) + #AZFG-CFPP-1
  annotate(geom = "text", x=4, y=2.3, label="a", fontface='italic', color="tomato3", size=9) + #AZFG-CFPP-2
  annotate(geom = "text", x=7, y=1.9, label="abcde", fontface='italic', color="darkgray", size=9) + #AZFG-CFPP-3
  annotate(geom = "text", x=10, y=1.65, label="abcde", fontface='italic', color="darkgray", size=9) + #COMT-MTPP-1
  annotate(geom = "text", x=13, y=1.45, label="de", fontface='italic', color="steelblue3", size=9) + #COMT-MTPP-2
  annotate(geom = "text", x=16, y=1.2, label="e", fontface='italic', color="steelblue3", size=9) + #COMT-MTPP-3
  annotate(geom = "text", x=19, y=1.7, label="abcde", fontface='italic', color="darkgray", size=9) + #COMT-MTPP-5
  annotate(geom = "text", x=22, y=2, label="abc", fontface='italic', color="darkgray", size=9) + #COMT-MTPP-6
  annotate(geom = "text", x=25, y=1.85, label="abc", fontface='italic', color="darkgray", size=9) + #NMSV-MAPJ-1
  annotate(geom = "text", x=28, y=1.6, label="bcde", fontface='italic', color="darkgray", size=9) + #NMSV-MAPJ-2
  annotate(geom = "text", x=31, y=2, label="ab", fontface='italic', color="tomato3", size=9) + #NMSV-MAPJ-3
  annotate(geom = "text", x=34, y=1.35, label="abc", fontface='italic', color="darkgray", size=9) + #NMSV-SVPJ-1
  annotate(geom = "text", x=37, y=1.55, label="bcde", fontface='italic', color="darkgray", size=9) + #NMVC-VCPP-1
  annotate(geom = "text", x=40, y=1.85, label="abc", fontface='italic', color="darkgray", size=9) + #NMVC-VCPP-2
  annotate(geom = "text", x=43, y=1.66, label="abcde", fontface='italic', color="darkgray", size=9) + #NVSH-DNPJ-1
  annotate(geom = "text", x=46, y=1.85, label="abc", fontface='italic', color="darkgray", size=9) + #NVSH-DNPJ-2
  annotate(geom = "text", x=49, y=1.4, label="bcde", fontface='italic', color="darkgray", size=9) + #NVSH-DNPJ-3
  annotate(geom = "text", x=52, y=1.75, label="abc", fontface='italic', color="darkgray", size=9) + #NVSP-BCPP-1
  annotate(geom = "text", x=55, y=1.6, label="cde", fontface='italic', color="steelblue3", size=9) + #NVSP-BCPP-4
  annotate(geom = "text", x=58, y=1.4, label="abcde", fontface='italic', color="darkgray", size=9) + #NVSP-BCPP-5
  annotate(geom = "text", x=61, y=1.65, label="cde", fontface='italic', color="steelblue3", size=9) + #NVSP-BTPJ-1
  annotate(geom = "text", x=64, y=1.8, label="abc", fontface='italic', color="darkgray", size=9) + #NVSP-BTPJ-2
  annotate(geom = "text", x=67, y=1.65, label="abcde", fontface='italic', color="darkgray", size=9) + #NVSP-DLPP-1
  annotate(geom = "text", x=70, y=1.9, label="abcde", fontface='italic', color="darkgray", size=9) #NVSP-WWPP-1





#only significant relationship across 25 plots - herbaceous cover and h ratio
x <- df$herb
y <- df$meanhratio
z <- lm(y ~ x)
summary(z)

fig8b <- ggplot() +
  geom_point(df, mapping=aes(x=herb, y=meanhratio, color=sig), size=5) +
  geom_smooth(data=df, mapping = aes(x=herb, y=meanhratio), color="black", size=2, 
              method = "lm", 
              formula = y~(x), se=FALSE) +
  scale_colour_manual(values = c("steelblue3", "darkgray", "tomato3")) +
  labs(x="Herbaceous cover (%)", y="Height ratio") +
  scale_x_continuous(expand = c(0,0), breaks = c(0,10,20,30,40,50,60), limits = c(0,60)) +
  scale_y_continuous(expand = c(0,0), breaks = c(0.7,0.8,0.9,1,1.1,1.2), limits = c(0.7,1.21)) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size=28), legend.position = c(0.92,0.2),
        legend.background = element_rect(linetype = "solid", colour = "black"), legend.spacing.y = unit(1, "mm"), 
        plot.title = element_text(size = 28),
        axis.text = element_text(size=28), axis.title = element_text(size = 28), axis.title.x = element_text(size = 28),
        axis.text.x = element_text(size = 28), plot.margin = unit(c(0.25,0.7,0.25,0.45),"cm")) +
  annotate(geom = "text", x=13, y=0.75, label="HR = 0.004*HC + 0.92", color="black", size=11) + 
  annotate(geom = "text", x=32, y=0.75, label=~R^2==~"0.25, p = 0.007", color="black", size=11)





#COMT diffs between plots - h ratio
#ba
panel1 <- ggplot() +
  geom_hline(yintercept = 1, linetype=2) +
  geom_boxplot(allspeccomt, mapping= aes(x=ba, y=hratio, group=plot), width=1.8, fill="gold3", outlier.shape = NA) +
  stat_boxplot(allspeccomt, mapping= aes(x=ba, y=hratio, group=plot), width=1.2,  geom = "errorbar") +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x = expression(paste("Basal area ", (m^2 * ha^-1))), y="Height ratio", title="COMT") +
  scale_x_continuous(expand = c(0,0), breaks = c(30,34,38), limits = c(28,40)) +
  scale_y_continuous(expand = c(0,0), breaks = c(0,1,2,3), limits = c(0,3.1)) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size=26), legend.position = "bottom",
        legend.background = element_blank(), legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 28),
        axis.text = element_text(size=26), axis.title = element_text(size = 26), axis.text.x = element_text(size = 26), 
        plot.margin = unit(c(0.5,0,0,0),"cm")) +
  annotate(geom = "text", x=30, y=2.2, label="abc", fontface='italic', color="darkgray", size=10) +
  annotate(geom = "text", x=34, y=1.85, label="de", fontface='italic', color="steelblue3", size=10) +
  annotate(geom = "text", x=38, y=1.6, label="e", fontface='italic', color="steelblue3", size=10) #+
  #annotate(geom = "text", x=39.5, y=2.8, label="c", fontface='bold', color="black", size=13)
  


#litter depth
panel2 <- ggplot() +
  geom_hline(yintercept = 1, linetype=2) +
  geom_boxplot(allspeccomt, mapping= aes(x=litterdep, y=hratio, group=plot), fill="gold3", outlier.shape = NA,
               width=0.5) +
  stat_boxplot(allspeccomt, mapping= aes(x=litterdep, y=hratio, group=plot), width=0.35,  geom = "errorbar") +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x = "Litter depth (cm)", y="Height ratio", title="COMT") +
  scale_x_continuous(expand = c(0,0), breaks = c(2,3,4,5), limits = c(2,5.5)) +
  scale_y_continuous(expand = c(0,0), breaks = c(0,1,2,3), limits = c(0,3.2)) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size=26), legend.position = "bottom",
        legend.background = element_blank(), legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 28),
        axis.text = element_text(size=26), axis.title = element_text(size = 26), axis.text.x = element_text(size = 26),
        axis.title.y = element_blank()) +
  annotate(geom = "text", x=2.6, y=2, label="abc", fontface='italic', color="darkgray", size=10) +
  annotate(geom = "text", x=3.8, y=1.65, label="de", fontface='italic', color="steelblue3", size=10) +
  annotate(geom = "text", x=4.8, y=1.4, label="e", fontface='italic', color="steelblue3", size=10) #+
  #annotate(geom = "text", x=5.4, y=2.75, label="d", fontface='bold', color="black", size=13)


#organic matter
panel3 <- ggplot() +
  geom_hline(yintercept = 1, linetype=2) +
  geom_boxplot(allspeccomt, mapping= aes(x=om, y=hratio, group=plot), fill="gold3", width=1.2, outlier.shape = NA) +
  stat_boxplot(allspeccomt, mapping= aes(x=om, y=hratio, group=plot), width=0.8,  geom = "errorbar") +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x = "Soil organic carbon (%)", y="Height ratio", title="COMT") +
  scale_x_continuous(expand = c(0,0), breaks = c(4,6,8,10), limits = c(3,10.5)) +
  scale_y_continuous(expand = c(0,0), breaks = c(0,1,2,3), limits = c(0,3.1)) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size=26), legend.position = "bottom",
        legend.background = element_blank(), legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 28),
        axis.text = element_text(size=26), axis.title = element_text(size = 26), axis.text.x = element_text(size = 26)) +
  annotate(geom = "text", x=8.78, y=2, label="abc", fontface='italic', color="darkgray", size=10) +
  annotate(geom = "text", x=6.4, y=1.7, label="de", fontface='italic', color="steelblue3", size=10) +
  annotate(geom = "text", x=4.16, y=1.4, label="e", fontface='italic', color="steelblue3", size=10) #+
  #annotate(geom = "text", x=10.17, y=2.8, label="e", fontface='bold', color="black", size=13)


#% clay
panel4 <- ggplot() +
  geom_hline(yintercept = 1, linetype=2) +
  geom_boxplot(allspeccomt, mapping= aes(x=clay, y=hratio, group=plot), fill="gold3", outlier.shape = NA, width=2.1) +
  stat_boxplot(allspeccomt, mapping= aes(x=clay, y=hratio, group=plot), width=1.3,  geom = "errorbar") +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x = "Clay content (%)", y="Height ratio", title="COMT") +
  scale_x_continuous(expand = c(0,0), breaks = c(6,10,14,18), limits = c(4,20)) +
  scale_y_continuous(expand = c(0,0), breaks = c(0,1,2,3), limits = c(0,3.2)) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size=20), legend.position = "bottom",
        legend.background = element_blank(), legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 28),
        axis.text = element_text(size=26), axis.title = element_text(size = 26), axis.text.x = element_text(size = 26),
        axis.title.y = element_blank()) +
  annotate(geom = "text", x=18, y=2, label="abc", fontface='italic', color="darkgray", size=10) +
  annotate(geom = "text", x=10, y=1.7, label="de", fontface='italic', color="steelblue3", size=10) +
  annotate(geom = "text", x=6, y=1.4, label="e", fontface='italic', color="steelblue3", size=10) #+
  #annotate(geom = "text", x=19.5, y=2.7, label="f", fontface='bold', color="black", size=13)




fig8cc <- plot_grid(panel1,panel2,panel3,panel4, nrow = 2, align = "vh")
fig8c <- fig8cc + theme(plot.margin = unit(c(0,0.5,0,0.7),"cm"))

plot_grid(fig8a,fig8b,fig8c, nrow = 3, rel_heights = c(0.95,0.77,1.1),
          labels = c("a", "b", ""),
          label_x = 0.92, label_y = 1, label_size = 40)








#figure 9#########################################################################################################
#test for sig diffs between plots - d ratio
options(repr.plot.width = 10, repr.plot.height = 8)
fig9a <- ggplot() +
  geom_hline(yintercept = 1, linetype=2) +
  geom_boxplot(allspec, mapping= aes(x=num, y=dratio, group=plot, fill=core), width=1.5, outlier.shape = NA) +
  stat_boxplot(allspec, mapping= aes(x=num, y=dratio, group=plot), geom = "errorbar", width=0.75) +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x="Study sites", y="Diameter ratio") +
  scale_y_continuous(expand = c(0,0), breaks = c(0,1,2,3), limits = c(0,3)) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size=34), legend.position = "top",
        legend.margin = margin(0,0,0,0),
        legend.background = element_blank(), legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 28),
        axis.text = element_text(size=28), axis.title = element_text(size = 28), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), plot.margin = unit(c(0.5,0.5,0.25,0.5),"cm")) +
  annotate(geom = "text", x=1, y=1.9, label="abcd", fontface='italic', color="darkgray", size=9) + #AZFG-CFPP-1
  annotate(geom = "text", x=4, y=2.1, label="ab", fontface='italic', color="tomato3", size=9) + #AZFG-CFPP-2
  annotate(geom = "text", x=7, y=1.65, label="bcde", fontface='italic', color="darkgray", size=9) + #AZFG-CFPP-3
  annotate(geom = "text", x=10, y=2.2, label="a", fontface='italic', color="tomato3", size=9) + #COMT-MTPP-1
  annotate(geom = "text", x=13, y=1.3, label="cde", fontface='italic', color="steelblue3", size=9) + #COMT-MTPP-2
  annotate(geom = "text", x=16, y=1.15, label="e", fontface='italic', color="steelblue3", size=9) + #COMT-MTPP-3
  annotate(geom = "text", x=19, y=2.3, label="abcde", fontface='italic', color="darkgray", size=9) + #COMT-MTPP-5
  annotate(geom = "text", x=22, y=1.6, label="bcde", fontface='italic', color="darkgray", size=9) + #COMT-MTPP-6
  annotate(geom = "text", x=25, y=1.8, label="bcde", fontface='italic', color="darkgray", size=9) + #NMSV-MAPJ-1
  annotate(geom = "text", x=28, y=2, label="abcde", fontface='italic', color="darkgray", size=9) + #NMSV-MAPJ-2
  annotate(geom = "text", x=31, y=2.2, label="abcde", fontface='italic', color="darkgray", size=9) + #NMSV-MAPJ-3
  annotate(geom = "text", x=34, y=2.4, label="a", fontface='italic', color="tomato3", size=9) + #NMSV-SVPJ-1
  annotate(geom = "text", x=37, y=1.5, label="de", fontface='italic', color="steelblue3", size=9) + #NMVC-VCPP-1
  annotate(geom = "text", x=40, y=2, label="bcde", fontface='italic', color="darkgray", size=9) + #NMVC-VCPP-2
  annotate(geom = "text", x=43, y=1.5, label="abcde", fontface='italic', color="darkgray", size=9) + #NVSH-DNPJ-1
  annotate(geom = "text", x=46, y=2.15, label="abcde", fontface='italic', color="darkgray", size=9) + #NVSH-DNPJ-2
  annotate(geom = "text", x=49, y=1.35, label="bcde", fontface='italic', color="darkgray", size=9) + #NVSH-DNPJ-3
  annotate(geom = "text", x=52, y=2.05, label="abcde", fontface='italic', color="darkgray", size=9) + #NVSP-BCPP-1
  annotate(geom = "text", x=55, y=1.85, label="bcde", fontface='italic', color="darkgray", size=9) + #NVSP-BCPP-4
  annotate(geom = "text", x=58, y=1.6, label="abcde", fontface='italic', color="darkgray", size=9) + #NVSP-BCPP-5
  annotate(geom = "text", x=61, y=1.8, label="bcde", fontface='italic', color="darkgray", size=9) + #NVSP-BTPJ-1
  annotate(geom = "text", x=64, y=2.45, label="abc", fontface='italic', color="darkgray", size=9) + #NVSP-BTPJ-2
  annotate(geom = "text", x=67, y=2.75, label="a", fontface='italic', color="tomato3", size=9) + #NVSP-DLPP-1
  annotate(geom = "text", x=70, y=1.95, label="bcde", fontface='italic', color="darkgray", size=9) #NVSP-WWPP-1

options(repr.plot.width = 10, repr.plot.height = 8)


plot_grid(fig8a, fig8b, fig9a, nrow = 3,
          labels = c("a", "b", "c"), align = "v",
          label_x = 0.94, label_y = 1, label_size = 46)


#COMT diffs between plots - d ratio
#canopy
panel1d <- ggplot() +
  geom_hline(yintercept = 1, linetype=2) +
  geom_boxplot(allspeccomt, mapping= aes(x=canopy, y=dratio, group=plot), width=6, fill="gold3", outlier.shape = NA) +
  stat_boxplot(allspeccomt, mapping= aes(x=canopy, y=dratio, group=plot), width=4,  geom = "errorbar") +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x = "Canopy cover (%)", y="Diameter ratio", title="COMT") +
  scale_x_continuous(expand = c(0,0), breaks = c(50,60,70,80,90), limits = c(49.5,91)) +
  scale_y_continuous(expand = c(0,0), breaks = c(0,1,2,3), limits = c(0,3)) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size=20), legend.position = "bottom",
        legend.background = element_blank(), legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 28),
        axis.text = element_text(size=26), axis.title = element_text(size = 26), axis.text.x = element_text(size = 26)) +
  annotate(geom = "text", x=54.1, y=2.7, label="a", fontface='italic', color="tomato3", size=10) +
  annotate(geom = "text", x=69.7, y=1.7, label="cde", fontface='italic', color="steelblue3", size=10) +
  annotate(geom = "text", x=82.1, y=1.6, label="e", fontface='italic', color="steelblue3", size=10) #+
  #annotate(geom = "text", x=89.7, y=2.65, label="b", fontface='bold', color="black", size=13)



#NVSP
#om
panel2d <- ggplot() +
  geom_hline(yintercept = 1, linetype=2) +
  geom_boxplot(allspecnvsp, mapping= aes(x=om, y=dratio, group=plot), fill="orchid2", outlier.shape = NA) +
  stat_boxplot(allspecnvsp, mapping= aes(x=om, y=dratio, group=plot),  geom = "errorbar") +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x = "Soil organic carbon (%)", y="Diameter ratio", title="NVSP") +
  scale_x_continuous(expand = c(0,0), breaks = c(5,6,7,8,9), limits = c(5,9.2)) +
  scale_y_continuous(expand = c(0,0), breaks = c(0,1,2,3), limits = c(0,3.3)) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size=20), legend.position = "bottom",
        legend.background = element_blank(), legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 28),
        axis.text = element_text(size=26), axis.title = element_text(size = 26), axis.text.x = element_text(size = 26),
        axis.title.y = element_blank()) +
  annotate(geom = "text", x=8.75, y=3.1, label="a", fontface='italic', color="tomato3", size=10) +
  annotate(geom = "text", x=6.79, y=2.2, label="bcde", fontface='italic', color="darkgray", size=10) +
  annotate(geom = "text", x=6.18, y=2.9, label="bcde", fontface='italic', color="darkgray", size=10) +
  annotate(geom = "text", x=5.78, y=2.25, label="bcde", fontface='italic', color="darkgray", size=10) #+
  #annotate(geom = "text", x=9.1, y=2.8, label="c", fontface='bold', color="black", size=13)



#NMSV
#canopy
panel3d <- ggplot() +
  geom_hline(yintercept = 1, linetype=2) +
  geom_boxplot(allspecnmsv, mapping= aes(x=canopy, y=dratio, group=plot), width=6, fill="green3", outlier.shape = NA) +
  stat_boxplot(allspecnmsv, mapping= aes(x=canopy, y=dratio, group=plot), width=4, geom = "errorbar") +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x = "Canopy cover (%)", y="Diameter ratio", title="NMSV") +
  scale_x_continuous(expand = c(0,0), breaks = c(15,25,35,45,55), limits = c(14.75,56)) +
  scale_y_continuous(expand = c(0,0), breaks = c(0,1,2,3), limits = c(0,3)) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size=20), legend.position = "bottom",
        legend.background = element_blank(), legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 28),
        axis.text = element_text(size=26), axis.title = element_text(size = 26), axis.text.x = element_text(size = 26)) +
  annotate(geom = "text", x=45.4, y=2.7, label="a", fontface='italic', color="tomato3", size=10) +
  annotate(geom = "text", x=19, y=2.3, label="bcde", fontface='italic', color="darkgray", size=10) #+
  #annotate(geom = "text", x=55, y=2.6, label="d", fontface='bold', color="black", size=13)


#number of adults
panel4d <- ggplot() +
  geom_hline(yintercept = 1, linetype=2) +
  geom_boxplot(allspecnmsv, mapping= aes(x=adults, y=dratio, group=plot), width=6, fill="green3", outlier.shape = NA) +
  stat_boxplot(allspecnmsv, mapping= aes(x=adults, y=dratio, group=plot), width=4, geom = "errorbar") +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x = "# adult trees", y="Diameter ratio", title="NMSV") +
  scale_x_continuous(expand = c(0,0), breaks = c(50,60,70,80,90), limits = c(49,91)) +
  scale_y_continuous(expand = c(0,0), breaks = c(0,1,2,3), limits = c(0,3)) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size=20), legend.position = "bottom",
        legend.background = element_blank(), legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 28),
        axis.text = element_text(size=26), axis.title = element_text(size = 26), axis.text.x = element_text(size = 26),
        axis.title.y = element_blank()) +
  annotate(geom = "text", x=53, y=2.7, label="a", fontface='italic', color="tomato3", size=10) +
  annotate(geom = "text", x=85, y=2.3, label="bcde", fontface='italic', color="darkgray", size=10) #+
  #annotate(geom = "text", x=90, y=2.8, label="e", fontface='bold', color="black", size=13)


#herb cover
panel5d <- ggplot() +
  geom_hline(yintercept = 1, linetype=2) +
  geom_boxplot(allspecnmsv, mapping= aes(x=herb, y=dratio, group=plot), width=5, fill="green3", outlier.shape = NA) +
  stat_boxplot(allspecnmsv, mapping= aes(x=herb, y=dratio, group=plot), width=3, geom = "errorbar") +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x = "Herbaceous cover (%)", y="Diameter ratio", title="NMSV") +
  scale_x_continuous(expand = c(0,0), breaks = c(10,20,30,40), limits = c(9.75,41)) +
  scale_y_continuous(expand = c(0,0), breaks = c(0,1,2,3), limits = c(0,3)) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size=20), legend.position = "bottom",
        legend.background = element_blank(), legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 28),
        axis.text = element_text(size=26), axis.title = element_text(size = 26), axis.text.x = element_text(size = 26)) +
  annotate(geom = "text", x=32.62, y=2.7, label="a", fontface='italic', color="tomato3", size=10) +
  annotate(geom = "text", x=16.75, y=2.3, label="bcde", fontface='italic', color="darkgray", size=10) #+
  #annotate(geom = "text", x=40, y=2.55, label="f", fontface='bold', color="black", size=13)


#litter cov
panel6d <- ggplot() +
  geom_hline(yintercept = 1, linetype=2) +
  geom_boxplot(allspecnmsv, mapping= aes(x=littercov, y=dratio, group=plot), width=9, fill="green3", outlier.shape = NA) +
  stat_boxplot(allspecnmsv, mapping= aes(x=littercov, y=dratio, group=plot), width=7, geom = "errorbar") +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x = "Litter cover (%)", y="Diameter ratio", title="NMSV") +
  scale_x_continuous(expand = c(0,0), breaks = c(0,20,40,60), limits = c(0,63)) +
  scale_y_continuous(expand = c(0,0), breaks = c(0,1,2,3), limits = c(0,3)) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size=20), legend.position = "bottom",
        legend.background = element_blank(), legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 28),
        axis.text = element_text(size=26), axis.title = element_text(size = 26), axis.text.x = element_text(size = 26),
        axis.title.y = element_blank()) +
  annotate(geom = "text", x=11.87, y=2.7, label="a", fontface='italic', color="tomato3", size=10) +
  annotate(geom = "text", x=54.37, y=2.3, label="bcde", fontface='italic', color="darkgray", size=10) #+
  #annotate(geom = "text", x=61, y=2.8, label="g", fontface='bold', color="black", size=13)


fig9bb <- plot_grid(panel1d,panel2d,panel3d,panel4d,panel5d,panel6d, nrow = 3)
fig9b <- fig9bb + theme(plot.margin = unit(c(0,0.5,0.25,0.5),"cm"))

plot_grid(fig9a, fig9b, nrow = 2, rel_heights = c(1/2,1),
          labels = c("a", ""),
          label_x = 0.95, label_y = 1, label_size = 40)


plot_grid(panel1, panel2, panel3, panel4, panel1d, panel2d, panel3d, panel4d, panel5d, panel6d,
          nrow = 5,
          labels = c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j"),
          label_x = 0.92, label_y = 1, label_size = 40)








#correlation matrices figure 10 and 11#############################################################################
install.packages("corrplot")
library(corrplot)


#height############################################################################################

#juniper######################################
jcor150 <- jleafarea150[, c(6,14,11,16,13,10)]
colnames(jcor150) <- c("height", "leaf mass", "leaf area", "root mass", "root area", "rooting depth")
jcor300 <- jleafarea300[, c(6,14,11,16,13,10)]
colnames(jcor300) <- c("height", "leaf mass", "leaf area", "root mass", "root area", "rooting depth")


juncor150 <- cor(jcor150[-1], jcor150$height)
juncor150[1:5,1] <- c(0.93,0.92,0.75,0.95,0)

juncor300 <- cor(jcor300[-1], jcor300$height)
juncor300[1:5,1] <- c(0.65,0,0.69,0,0)

juncor<- cbind(juncor150, juncor300)
colnames(juncor) <- c("< 150", "150-300+")


corrplot(juncor, method = "ellipse", sig.level = 0.05, insig = "blank", tl.col = "black", 
         tl.srt = 0, tl.cex = 2,  tl.offset = 0.5, tl.pos = "lt",
         cl.cex = 2, cl.ratio = 0.5, mar = c(0.5,0.5,4,1))
text(1.58, 6, labels = "Juniper", cex = 2.5)
par(xpd=TRUE)






#pinon#####################################################
picor150 <- pinleafarea150[, c(6,14,11,16,13,10)]
colnames(picor150) <- c("height", "leaf mass", "leaf area", "root mass", "root area", "rooting depth")
picor300 <- pinleafarea300[, c(6,14,11,16,13,10)]
colnames(picor300) <- c("height", "leaf mass", "leaf area", "root mass", "root area", "rooting depth")


pincor150 <- cor(picor150[-1], picor150$height)
pincor150[1:5,1] <- c(0.84,0.85,0.89,0.58,0)

pincor300 <- cor(picor300[-1], picor300$height)
pincor300[1:5,1] <- c(0.75,0.82,0.85,0.75,-0.87)

pincor<- cbind(pincor150, pincor300)
colnames(pincor) <- c("< 150", "150-300+")


corrplot(pincor, method = "ellipse", sig.level = 0.05, insig = "blank", tl.col = "black", 
         tl.srt = 0, tl.cex = 2,  tl.offset = 0.5, tl.pos = "lt",
         cl.cex = 2, cl.ratio = 0.5, mar = c(0.5,0.5,4,1))
text(1.58, 6, labels = "Piñon pine", cex = 2.5)





#ponderosa################################################
pocor150 <- ponleafarea150[, c(6,14,11,16,13,10)]
colnames(pocor150) <- c("height", "leaf mass", "leaf area", "root mass", "root area", "rooting depth")
pocor150300 <- ponleafarea150300[, c(6,14,11,16,13,10)]
colnames(pocor150300) <- c("height", "leaf mass", "leaf area", "root mass", "root area", "rooting depth")
pocor300 <- ponleafarea300[, c(6,14,11,16,13,10)]
colnames(pocor300) <- c("height", "leaf mass", "leaf area", "root mass", "root area", "rooting depth")


poncor150 <- cor(pocor150[-1], pocor150$height)
poncor150[1:5,1] <- c(0.79,0.85,0.77,0.88,-0.58)

poncor150300 <- cor(pocor150300[-1], pocor150300$height)
poncor150300[1:5,1] <- c(0.56,0.59,0,0.42,0)

poncor300 <- cor(pocor300[-1], pocor300$height)
poncor300[1:5,1] <- c(0.73,0.75,0.6,0.81,-0.64)

poncor<- cbind(poncor150, poncor150300, poncor300)
colnames(poncor) <- c("< 150", "150-300", "300+")


corrplot(poncor, method = "ellipse", sig.level = 0.05, insig = "blank", tl.col = "black", 
         tl.srt = 0, tl.cex = 1.75,  tl.offset = 0.5, tl.pos = "lt",
         cl.cex = 1.75, cl.ratio = 0.5, mar = c(0.5,0.5,1,1))
text(2, 6, labels = "Ponderosa", cex = 2.5)





#fir##############################################
fcor150 <- fleafarea150[, c(6,14,11,16,13,10)]
colnames(fcor150) <- c("height", "leaf mass", "leaf area", "root mass", "root area", "rooting depth")
fcor150300 <- fleafarea150300[, c(6,14,11,16,13,10)]
colnames(fcor150300) <- c("height", "leaf mass", "leaf area", "root mass", "root area", "rooting depth")


fircor150 <- cor(fcor150[-1], fcor150$height)
fircor150[1:5,1] <- c(0.85,0.88,0.89,0.61,0)

fircor150300 <- cor(fcor150300[-1], fcor150300$height)
fircor150300[1:5,1] <- c(0.67,0.79,0.71,0.71,0)

fircor<- cbind(fircor150, fircor150300)
colnames(fircor) <- c("< 150", "150-300")


corrplot(fircor, method = "ellipse", sig.level = 0.05, insig = "blank", tl.col = "black", 
         tl.srt = 0, tl.cex = 2,  tl.offset = 0.5, tl.pos = "lt",
         cl.cex = 2, cl.ratio = 0.5, mar = c(0.5,0.5,4,1))
text(1.58, 6, labels = "Fir", cex = 2.5)



cor <- cbind(juncor, pincor, poncor, fircor)

corrplot(cor, method = "ellipse", tl.col = "black", col = col1(25),
         tl.srt = 0, tl.cex = 2,  tl.offset = 0.25, tl.pos = "lt",
         cl.cex = 2, cl.ratio = 0.1, mar = c(0.5,0.5,4,2))
text(1.58, 6.1, labels = "Juniper", cex = 3)
text(3.5, 6.1, labels = "Piñon pine", cex = 3)
text(6, 6.1, labels = "Ponderosa pine", cex = 3)
text(8.5, 6.1, labels = "Fir", cex = 3)
rect(0.55, 0.55, 1.45, 1.45, border = FALSE, col = "white")
rect(1.55, 1.55, 2.45, 2.45, border = FALSE, col = "white")
rect(1.55, 0.55, 2.45, 1.45, border = FALSE, col = "white")
rect(1.55, 3.55, 2.45, 4.45, border = FALSE, col = "white")
rect(2.55, 0.55, 3.45, 1.45, border = FALSE, col = "white")
rect(5.55, 0.55, 6.45, 1.45, border = FALSE, col = "white")
rect(5.55, 2.55, 6.45, 3.45, border = FALSE, col = "white")
rect(7.55, 0.55, 8.45, 1.45, border = FALSE, col = "white")
rect(8.55, 0.55, 9.45, 1.45, border = FALSE, col = "white")
rect(0.5, 0.5, 2.5, 5.5, border = TRUE, lwd=3)
rect(2.5, 0.5, 4.5, 5.5, border = TRUE, lwd=3)
rect(4.5, 0.5, 7.5, 5.5, border = TRUE, lwd=3)
rect(7.5, 0.5, 9.5, 5.5, border = TRUE, lwd=3)
text(5.5, 0.06, labels = "Height (mm)", cex = 3)
text(-1.1, 3.2, srt=90, labels = "Variable", cex = 3, xpd=TRUE)







#diameter############################################################################################

#juniper###########################################
jcord150 <- jleafarea150[, c(8,14,11,16,13,10)]
colnames(jcord150) <- c("diameter", "leaf mass", "leaf area", "root mass", "root area", "rooting depth")
jcord300 <- jleafarea300[, c(8,14,11,16,13,10)]
colnames(jcord300) <- c("diameter", "leaf mass", "leaf area", "root mass", "root area", "rooting depth")


juncord150 <- cor(jcord150[-1], jcord150$diameter)
juncord150[1:5,1] <- c(0.94,0.5,0.62,0.55,0)

juncord300 <- cor(jcord300[-1], jcord300$diameter)
juncord300[1:5,1] <- c(0.67,0.5,0.84,0,0)

juncord <- cbind(juncord150, juncord300)
colnames(juncord) <- c("< 150", "150-300+")


corrplot(juncord, method = "ellipse", sig.level = 0.05, insig = "blank", tl.col = "black", 
         tl.srt = 0, tl.cex = 1.75,  tl.offset = 0.5, tl.pos = "lt",
         cl.cex = 1.75, cl.ratio = 0.5, mar = c(0.5,0.5,4,1))
text(1.58, 6, labels = "Juniper", cex = 2.5)
par(xpd=TRUE)





#pinon################################################
picord150 <- pinleafarea150[, c(8,14,11,16,13,10)]
colnames(picord150) <- c("diameter", "leaf mass", "leaf area", "root mass", "root area", "rooting depth")
picord300 <- pinleafarea300[, c(8,14,11,16,13,10)]
colnames(picord300) <- c("diameter", "leaf mass", "leaf area", "root mass", "root area", "rooting depth")


pincord150 <- cor(picord150[-1], picord150$diameter)
pincord150[1:5,1] <- c(0.87,0.88,0.92,0.96,-0.33)

pincord300 <- cor(picord300[-1], picord300$diameter)
pincord300[1:5,1] <- c(0.78,0.84,0.87,0.94,0)

pincord<- cbind(pincord150, pincord300)
colnames(pincord) <- c("< 150", "150-300+")


corrplot(pincord, method = "ellipse", sig.level = 0.05, insig = "blank", tl.col = "black", 
         tl.srt = 0, tl.cex = 1.75,  tl.offset = 0.5, tl.pos = "lt",
         cl.cex = 1.75, cl.ratio = 0.5, mar = c(0.5,0.5,1,1))
text(1.58, 6, labels = "Piñon pine", cex = 2.5)





#ponderosa#############################################
pocord150 <- ponleafarea150[, c(8,14,11,16,13,10)]
colnames(pocord150) <- c("diameter", "leaf mass", "leaf area", "root mass", "root area", "rooting depth")
pocord150300 <- ponleafarea150300[, c(8,14,11,16,13,10)]
colnames(pocord150300) <- c("diameter", "leaf mass", "leaf area", "root mass", "root area", "rooting depth")
pocord300 <- ponleafarea300[, c(8,14,11,16,13,10)]
colnames(pocord300) <- c("diameter", "leaf mass", "leaf area", "root mass", "root area", "rooting depth")


poncord150 <- cor(pocord150[-1], pocord150$diameter)
poncord150[1:5,1] <- c(0.83,0.88,0.82,0.91,-0.62)

poncord150300 <- cor(pocord150300[-1], pocord150300$diameter)
poncord150300[1:5,1] <- c(0.91,0.93,0.83,0.62,0)

poncord300 <- cor(pocord300[-1], pocord300$diameter)
poncord300[1:5,1] <- c(0.87,0.88,0.83,0.97,-0.86)

poncord<- cbind(poncord150, poncord150300, poncord300)
colnames(poncord) <- c("< 150", "150-300", "300+")


corrplot(poncord, method = "ellipse", sig.level = 0.05, insig = "blank", tl.col = "black", 
         tl.srt = 0, tl.cex = 1.75,  tl.offset = 0.5, tl.pos = "lt",
         cl.cex = 1.75, cl.ratio = 0.5, mar = c(0.5,0.5,1,1))
text(2, 6, labels = "Ponderosa", cex = 2.5)





#fir##################################################
fcord150 <- fleafarea150[, c(8,14,11,16,13,10)]
colnames(fcord150) <- c("diameter", "leaf mass", "leaf area", "root mass", "root area", "rooting depth")
fcord150300 <- fleafarea150300[, c(8,14,11,16,13,10)]
colnames(fcord150300) <- c("diameter", "leaf mass", "leaf area", "root mass", "root area", "rooting depth")


fircord150 <- cor(fcord150[-1], fcord150$diameter)
fircord150[1:5,1] <- c(0.84,0.87,0.94,0.85,0)

fircord150300 <- cor(fcord150300[-1], fcord150300$diameter)
fircord150300[1:5,1] <- c(0.81,0.87,0.87,0.93,-0.96)

fircord<- cbind(fircord150, fircord150300)
colnames(fircord) <- c("< 150", "150-300")


corrplot(fircord, method = "ellipse", sig.level = 0.05, insig = "blank", tl.col = "black", 
         tl.srt = 0, tl.cex = 1.75,  tl.offset = 0.5, tl.pos = "lt",
         cl.cex = 1.75, cl.ratio = 0.5, mar = c(0.5,0.5,1,1))
text(1.58, 6, labels = "Fir", cex = 2.5)



cord <- cbind(juncord, pincord, poncord, fircord)
col1 <- colorRampPalette(c("tomato", "gold", "dodgerblue"))

corrplot(cord, method = "ellipse", tl.col = "black", col = col1(25),
         tl.srt = 0, tl.cex = 2,  tl.offset = 0.25, tl.pos = "lt",
         cl.cex = 2, cl.ratio = 0.1, mar = c(0.5,0.5,4,2))
text(1.58, 6.1, labels = "Juniper", cex = 3)
text(3.5, 6.1, labels = "Piñon pine", cex = 3)
text(6, 6.1, labels = "Ponderosa pine", cex = 3)
text(8.5, 6.1, labels = "Fir", cex = 3)
rect(0.55, 0.55, 1.45, 1.45, border = FALSE, col = "white")
rect(1.55, 1.55, 2.45, 2.45, border = FALSE, col = "white")
rect(1.55, 0.55, 2.45, 1.45, border = FALSE, col = "white")
rect(3.55, 0.55, 4.45, 1.45, border = FALSE, col = "white")
rect(5.55, 0.55, 6.45, 1.45, border = FALSE, col = "white")
rect(7.55, 0.55, 8.45, 1.45, border = FALSE, col = "white")
rect(0.5, 0.5, 2.5, 5.5, border = TRUE, lwd=3)
rect(2.5, 0.5, 4.5, 5.5, border = TRUE, lwd=3)
rect(4.5, 0.5, 7.5, 5.5, border = TRUE, lwd=3)
rect(7.5, 0.5, 9.5, 5.5, border = TRUE, lwd=3)
text(5.5, 0.06, labels = "Diameter (mm)", cex = 3)
text(-1.1, 3.2, srt=90, labels = "Variable", cex = 3, xpd=TRUE)
















#ponderosa pine only for paper 2#############################################################
ponratio <- subset(allspec, Species == "Ponderosa") #& plot != "NMVC-VCPP-2" & core != "NVSH" &
#site != "DLPP") #remove sites with only 1 individual
ponratio$num <- gsub(37, 25, ponratio$num)
ponratio$num <- gsub(52, 28, ponratio$num)
ponratio$num <- gsub(55, 31, ponratio$num)
ponratio$num <- gsub(58, 34, ponratio$num)
ponratio$num <- as.integer(ponratio$num)
write.csv(ponratio, file = "ponratio.csv")
ponratio <- read.csv(file = "ponratio.csv")

# ponratioazfg <- subset(ponratio, core == "AZFG")
# ponrationmvc <- subset(ponratio, core == "NMVC")
# ponrationvsp <- subset(ponratio, core == "NVSP")

meanhratio <- with(ponratio, by(hratio,plot,mean))
meandratio <- with(ponratio, by(dratio,plot,mean))
dfpon <- data.frame(cbind(meanhratio, meandratio))
dfpon <- read.csv(file = "dfpon.csv")
dfpon$sigh <- factor(dfpon$sigh, levels = c("Low", "Average", "High"))
dfpon$sigd <- factor(dfpon$sigd, levels = c("Low", "Average", "High"))


#hratio#############################################################################################################
x <- dfpon$meanhratio
y <- dfpon$juvden
z <- lm(y ~ x)
summary(z)

x <- dfpon$meandratio
y <- dfpon$juvden
z <- lm(y ~ x)
summary(z)

#sig relationships across plots
x <- dfpon$herb
y <- dfpon$meanhratio
z <- lm(y ~ x)
summary(z)

x <- dfpon$mat
y <- dfpon$meanhratio
z <- lm(y ~ x)
summary(z)




ponanovah <- aov(hratio~plot, data = ponratio) #ponderosa only
tukeyh <- TukeyHSD(ponanovah)

generate_label_df <- function(tukeyh, variable){
  
  # Extract labels and factor levels from Tukey post-hoc 
  tukey.levels <- tukeyh[[variable]][,4]
  tukey.labels <- data.frame(multcompLetters(tukey.levels)['Letters'])
  
  #I need to put the labels in the same order as in the boxplot :
  tukey.labels$plot=rownames(tukey.labels)
  tukey.labels=tukey.labels[order(tukey.labels$plot) , ]
  return(tukey.labels)
}

# Apply the function on my dataset
labelsh=generate_label_df(tukeyh , "plot")



ggplot() +
  geom_hline(yintercept = 1, linetype=2) +
  geom_boxplot(ponratio, mapping= aes(x=num, y=hratio, group=plot, fill=core), 
               width=1.5, outlier.shape = NA) +
  stat_boxplot(ponratio, mapping= aes(x=num, y=hratio, group=plot), geom = "errorbar", width=0.75) +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x="", y="Height ratio") +
  scale_x_continuous(limits = c(0,35)) +
  scale_y_continuous(expand = c(0,0), breaks = c(0.0,1.0,2.0,3.0), limits = c(0.0,3.0)) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size=26), legend.position = "top",
        legend.margin = margin(0,0,0,0),
        legend.background = element_blank(), legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 28),
        axis.text = element_text(size=26), axis.title = element_text(size = 26), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), plot.margin = unit(c(0.5,1,0.25,0.95),"cm")) +
  annotate(geom = "text", x=1, y=2.2, label="ad", fontface='italic', color="tomato3", size=8.5) + #AZFG-CFPP-1
  annotate(geom = "text", x=4, y=2.3, label="a", fontface='italic', color="tomato3", size=8.5) + #AZFG-CFPP-2
  annotate(geom = "text", x=7, y=1.9, label="abcd", fontface='italic', color="darkgray", size=8.5) + #AZFG-CFPP-3
  annotate(geom = "text", x=10, y=1.65, label="abcd", fontface='italic', color="darkgray", size=8.5) + #COMT-MTPP-1
  annotate(geom = "text", x=13, y=1.45, label="c", fontface='italic', color="steelblue3", size=8.5) + #COMT-MTPP-2
  annotate(geom = "text", x=16, y=1.2, label="c", fontface='italic', color="steelblue3", size=8.5) + #COMT-MTPP-3
  annotate(geom = "text", x=19, y=1.7, label="abcd", fontface='italic', color="darkgray", size=8.5) + #COMT-MTPP-5
  annotate(geom = "text", x=22, y=2, label="abd", fontface='italic', color="tomato3", size=8.5) + #COMT-MTPP-6
  annotate(geom = "text", x=25, y=1.85, label="bc", fontface='italic', color="darkgray", size=8.5) + #NMVC-VCPP-1
  annotate(geom = "text", x=28, y=1.6, label="abcd", fontface='italic', color="darkgray", size=8.5) + #NVSP-BCPP-1
  annotate(geom = "text", x=31, y=2, label="bcd", fontface='italic', color="darkgray", size=8.5) + #NVSP-BCPP-4
  annotate(geom = "text", x=34, y=1.35, label="abcd", fontface='italic', color="darkgray", size=8.5) #NVSP-BCPP-5

sitemeanhratio <- ponratio %>% group_by(plot) %>% summarize(meanhratio = mean(hratio)) %>% as.data.frame
quantile(sitemeanhratio$meanhratio)
#1st quartile: 0.8475 cutoff - COMT-PP-2,3; NVSP-MC-2
#4th quartile: 1.0724 cuttoff - AZFG-PP-1,2; COMT-PP-6



#sig relationships across plots######################################################################################
#herb cover
fig1a <- ggplot() +
  geom_point(df, mapping=aes(x=herb, y=meanhratio, color=sigh), size=3) +
  geom_smooth(data=df, mapping = aes(x=herb, y=meanhratio), color="black", size=1.2, 
              method = "lm", 
              formula = y~(x), se=FALSE) +
  scale_colour_manual(values = c("steelblue3", "darkgray", "tomato3")) +
  labs(x="Herbaceous cover (%)", y="Height ratio") +
  scale_x_continuous(expand = c(0,0), breaks = c(0,10,20,30,40,50), limits = c(0,50)) +
  scale_y_continuous(expand = c(0,0), breaks = c(0.7,0.8,0.9,1,1.1,1.2), limits = c(0.7,1.2)) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size=26), legend.position = c(0.2,0.85),
        legend.background = element_rect(linetype = "solid", colour = "black"), legend.spacing.y = unit(1, "mm"), 
        plot.title = element_text(size = 28),
        axis.text = element_text(size=26), axis.title = element_text(size = 26), axis.title.x = element_text(size = 26),
        axis.text.x = element_text(size = 26), plot.margin = unit(c(0.75,0.8,0.25,0.45),"cm")) +
  annotate(geom = "text", x=28, y=0.75, label="HR = 0.006*HC + 0.84", color="black", size=10) + 
  annotate(geom = "text", x=42, y=0.75, label=~R^2==~"0.49, p = 0.007", color="black", size=10)


#mat
fig1b <- ggplot() +
  geom_point(df, mapping=aes(x=mat, y=meanhratio, color=sigh), size=3) +
  geom_smooth(data=df, mapping = aes(x=mat, y=meanhratio), color="black", size=1.2, 
              method = "lm", 
              formula = y~(x), se=FALSE) +
  scale_colour_manual(values = c("steelblue3", "darkgray", "tomato3")) +
  labs(x="Mean annual temperature (° C)", y="Height ratio") +
  scale_x_continuous(expand = c(0,0), breaks = c(5.5,6,6.5,7,7.5), limits = c(5.3,7.5)) +
  scale_y_continuous(expand = c(0,0), breaks = c(0.7,0.8,0.9,1,1.1,1.2), limits = c(0.7,1.2)) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_blank(), legend.position = "none",
        legend.background = element_blank(), legend.spacing.y = unit(1, "mm"), 
        plot.title = element_text(size = 28),
        axis.text = element_text(size=26), axis.title = element_text(size = 26), axis.title.x = element_text(size = 26),
        axis.text.x = element_text(size = 26), plot.margin = unit(c(0.75,0.75,0.25,0.45),"cm")) +
  annotate(geom = "text", x=6.5, y=0.75, label="HR = 0.12*MAT + 0.25", color="black", size=10) + 
  annotate(geom = "text", x=7.1, y=0.75, label=~R^2==~"0.29, p = 0.04", color="black", size=10)




#sig diffs within COMT-2,3,6############################################################################################
#figure 2
ponratiocomt <- subset(ponratio, core == "COMT")
ponratiocomt2 <- subset(ponratiocomt, plotnum == 2)
ponratiocomt3 <- subset(ponratiocomt, plotnum == 3)
ponratiocomt6 <- subset(ponratiocomt, plotnum == 6)
ponratiocomt <- as.data.frame(rbind(ponratiocomt2,ponratiocomt3,ponratiocomt6))

#ba
fig2a <- ggplot() +
  geom_hline(yintercept = 1, linetype=2) +
  geom_boxplot(ponratiocomt, mapping= aes(x=ba, y=hratio, group=plot), width=1.8, fill="gold3", outlier.shape = NA) +
  stat_boxplot(ponratiocomt, mapping= aes(x=ba, y=hratio, group=plot), width=1.2,  geom = "errorbar") +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x = expression(paste("Basal area ", (m^2 * ha^-1))), y="Height ratio", title="COMT") +
  scale_x_continuous(expand = c(0,0), breaks = c(30,34,38), limits = c(28,40)) +
  scale_y_continuous(expand = c(0,0), breaks = c(0,1,2,3), limits = c(0,3)) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size=30), legend.position = "bottom",
        legend.background = element_blank(), legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 36),
        axis.text = element_text(size=30), axis.title = element_text(size = 30), 
        axis.text.x = element_text(size = 30, vjust = -0.1), 
        plot.margin = unit(c(0.5,0,0,0),"cm")) +
  annotate(geom = "text", x=30, y=2, label="abd", fontface='italic', color="darkgray", size=12) +
  annotate(geom = "text", x=34, y=1.65, label="c", fontface='italic', color="steelblue3", size=12) +
  annotate(geom = "text", x=38, y=1.4, label="c", fontface='italic', color="steelblue3", size=12) 
#annotate(geom = "text", x=39.5, y=2.8, label="c", fontface='bold', color="black", size=13)


#clay
fig2b <- ggplot() +
  geom_hline(yintercept = 1, linetype=2) +
  geom_boxplot(ponratiocomt, mapping= aes(x=clay, y=hratio, group=plot), width=1.8, fill="gold3", outlier.shape = NA) +
  stat_boxplot(ponratiocomt, mapping= aes(x=clay, y=hratio, group=plot), width=1.2,  geom = "errorbar") +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x = "Clay (%)", y="Height ratio", title="COMT") +
  scale_x_continuous(expand = c(0,0), breaks = c(6,10,14,18), limits = c(5,19)) +
  scale_y_continuous(expand = c(0,0), breaks = c(0,1,2,3), limits = c(0,3)) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size=30), legend.position = "bottom",
        axis.title.y = element_blank(),
        legend.background = element_blank(), legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 36),
        axis.text = element_text(size=30), axis.title = element_text(size = 30), 
        axis.text.x = element_text(size = 30, vjust = -0.1), 
        plot.margin = unit(c(0.5,0,0,0),"cm")) +
  annotate(geom = "text", x=18, y=2, label="abd", fontface='italic', color="darkgray", size=12) +
  annotate(geom = "text", x=10, y=1.65, label="c", fontface='italic', color="steelblue3", size=12) +
  annotate(geom = "text", x=6, y=1.4, label="c", fontface='italic', color="steelblue3", size=12) 
#annotate(geom = "text", x=39.5, y=2.8, label="c", fontface='bold', color="black", size=13)


#litter depth
fig2c <- ggplot() +
  geom_hline(yintercept = 1, linetype=2) +
  geom_boxplot(ponratiocomt, mapping= aes(x=litterdep, y=hratio, group=plot), fill="gold3", outlier.shape = NA,
               width=0.5) +
  stat_boxplot(ponratiocomt, mapping= aes(x=litterdep, y=hratio, group=plot), width=0.35,  geom = "errorbar") +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x = "Litter depth (cm)", y="Height ratio", title="COMT") +
  scale_x_continuous(expand = c(0,0), breaks = c(2,3,4,5), limits = c(2,5.3)) +
  scale_y_continuous(expand = c(0,0), breaks = c(0,1,2,3), limits = c(0,3)) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size=30), legend.position = "bottom",
        legend.background = element_blank(), legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 36),
        axis.text = element_text(size=30), axis.title = element_text(size = 30), 
        axis.text.x = element_text(size = 30, vjust = -0.1),
        axis.title.y = element_text(size=30)) +
  annotate(geom = "text", x=2.6, y=2, label="abd", fontface='italic', color="darkgray", size=12) +
  annotate(geom = "text", x=3.8, y=1.65, label="c", fontface='italic', color="steelblue3", size=12) +
  annotate(geom = "text", x=4.8, y=1.4, label="c", fontface='italic', color="steelblue3", size=12) 
#annotate(geom = "text", x=5.4, y=2.75, label="d", fontface='bold', color="black", size=13)


#soc
fig2d <- ggplot() +
  geom_hline(yintercept = 1, linetype=2) +
  geom_boxplot(ponratiocomt, mapping= aes(x=om, y=hratio, group=plot), width=1, fill="gold3", outlier.shape = NA) +
  stat_boxplot(ponratiocomt, mapping= aes(x=om, y=hratio, group=plot), width=0.6,  geom = "errorbar") +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x = "Soil organic carbon (%)", y="Height ratio", title="COMT") +
  scale_x_continuous(expand = c(0,0), breaks = c(4,6,8,10), limits = c(3,10.5)) +
  scale_y_continuous(expand = c(0,0), breaks = c(0,1,2,3), limits = c(0,3)) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size=30), legend.position = "bottom",
        legend.background = element_blank(), legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 36),
        axis.title.y = element_blank(),
        axis.text = element_text(size=30), axis.title = element_text(size = 30), 
        axis.text.x = element_text(size = 30, vjust = -0.1)) +
  annotate(geom = "text", x=8.78, y=2, label="abd", fontface='italic', color="darkgray", size=12) +
  annotate(geom = "text", x=6.4, y=1.7, label="c", fontface='italic', color="steelblue3", size=12) +
  annotate(geom = "text", x=4.16, y=1.4, label="c", fontface='italic', color="steelblue3", size=12) 
#annotate(geom = "text", x=10.17, y=2.8, label="e", fontface='bold', color="black", size=13)


plot_grid(fig2a,fig2c,fig2b,fig2d, nrow = 2, #rel_heights = c(0.95,0.77,1.1),
          labels = c("a", "b", "c", "d"), align = "vh",
          label_x = 0.94, label_y = 1, label_size = 46)





#dratio#################################################################################################
ponanovad <- aov(dratio~plot, data = ponratio) #ponderosa only
tukeyd <- TukeyHSD(ponanovad)
install.packages("multcompView")
library(multcompView)
ponratio$plot <- gsub("-", "_", ponratio$plot)

generate_label_df <- function(tukeyd, variable){
  
  # Extract labels and factor levels from Tukey post-hoc 
  tukey.levels <- tukeyd[[variable]][,4]
  tukey.labels <- data.frame(multcompLetters(tukey.levels)['Letters'])
  
  #I need to put the labels in the same order as in the boxplot :
  tukey.labels$plot=rownames(tukey.labels)
  tukey.labels=tukey.labels[order(tukey.labels$plot) , ]
  return(tukey.labels)
}

# Apply the function on my dataset
labelsd=generate_label_df(tukeyd , "plot")

ggplot() +
  geom_hline(yintercept = 1, linetype=2) +
  geom_boxplot(ponratio, mapping= aes(x=num, y=dratio, group=plot, fill=core), 
               width=1.5, outlier.shape = NA) +
  stat_boxplot(ponratio, mapping= aes(x=num, y=dratio, group=plot), geom = "errorbar", width=0.75) +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x="", y="Diameter ratio") +
  scale_x_discrete(breaks = c(1,4,7,10,13,16,19,22,25,28,31,34)) +
  scale_y_continuous(expand = c(0,0), breaks = c(0.0,1.0,2.0,3.0), limits = c(0.0,3.0)) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size=26), legend.position = "top",
        legend.margin = margin(0,0,0,0),
        legend.background = element_blank(), legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 28),
        axis.text = element_text(size=26), axis.title = element_text(size = 26), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), plot.margin = unit(c(0.5,0.5,0.25,0.95),"cm")) +
  annotate(geom = "text", x=1, y=2, label="abe", fontface='italic', color="tomato3", size=8.5) + #AZFG-CFPP-1
  annotate(geom = "text", x=4, y=2.15, label="ab", fontface='italic', color="tomato3", size=8.5) + #AZFG-CFPP-2
  annotate(geom = "text", x=7, y=1.7, label="cde", fontface='italic', color="darkgray", size=8.5) + #AZFG-CFPP-3
  annotate(geom = "text", x=10, y=2.3, label="a", fontface='italic', color="tomato3", size=8.5) + #COMT-MTPP-1
  annotate(geom = "text", x=13, y=1.4, label="cd", fontface='italic', color="steelblue3", size=8.5) + #COMT-MTPP-2
  annotate(geom = "text", x=16, y=1.2, label="d", fontface='italic', color="steelblue3", size=8.5) + #COMT-MTPP-3
  annotate(geom = "text", x=19, y=2.4, label="abce", fontface='italic', color="darkgray", size=8.5) + #COMT-MTPP-5
  annotate(geom = "text", x=22, y=1.8, label="bcde", fontface='italic', color="darkgray", size=8.5) + #COMT-MTPP-6
  annotate(geom = "text", x=25, y=1.55, label="cd", fontface='italic', color="steelblue3", size=8.5) + #NMVC-VCPP-1
  annotate(geom = "text", x=28, y=1.8, label="abe", fontface='italic', color="tomato3", size=8.5) + #NVSP-BCPP-1
  annotate(geom = "text", x=31, y=2, label="bcde", fontface='italic', color="darkgray", size=8.5) + #NVSP-BCPP-4
  annotate(geom = "text", x=34, y=2.7, label="ab", fontface='italic', color="tomato3", size=8.5) #NVSP-BCPP-5

sitemeandratio <- ponratio %>% group_by(plot) %>% summarize(meandratio = mean(dratio)) %>% as.data.frame
quantile(sitemeandratio$meandratio)
#1st quartile: 0.8251 cutoff - COMT-PP-2,3; NMVC-PP-1
#4th quartile: 1.2031 cuttoff - COMT-PP-1; NVSP-MC-1,2

#characteristics across plots########################################################

corhratio <- df[, c(32,13,12,18,23)]
colnames(corhratio) <- c("height ratio", "herbaceous cover", "juvenile density", "MAT", "Max winter T")
cordratio <- df[, c(33,7,14)]
colnames(cordratio) <- c("diameter ratio", "canopy cover", "litter cover")


corhratio <- cor(corhratio[-1], corhratio$`height ratio`)
corhratio[1:4,1] <- c(0.7,0.52,0.54,0.53)

cordratio <- cor(cordratio[-1], cordratio$`diameter ratio`)
cordratio[1:2,1] <- c(-0.56,-0.62)

cor <- rbind(corhratio, cordratio)
colnames(cor) <- c("")


corrplot(cor, method = "ellipse", tl.col = "black", col = col1(25),
         tl.srt = 0, tl.cex = 2,  tl.offset = 0.25, tl.pos = "lt",
         cl.cex = 2, cl.ratio = 0.6, mar = c(1,0.5,1,2))
rect(0.5, 0.5, 1.5, 2.5, border = TRUE, lwd=4)
rect(0.5, 2.5, 1.5, 6.5, border = TRUE, lwd=4)
text(-1.5, 4.5, srt=90, labels = "Height ratio", cex = 3, xpd=TRUE)
text(-1.5, 1.5, srt=90, labels = "Diameter ratio", cex = 3, xpd=TRUE)


#canopy cover
x <- df$canopy
y <- df$meandratio
z <- lm(y ~ x)
summary(z)

fig1c <- ggplot() +
  geom_point(df, mapping=aes(x=canopy, y=meandratio, color=sigd), size=3) +
  geom_smooth(data=df, mapping = aes(x=canopy, y=meandratio), color="black", size=1.2, 
              method = "lm", 
              formula = y~(x), se=FALSE) +
  scale_colour_manual(values = c("steelblue3", "darkgray", "tomato3")) +
  labs(x="Canopy cover (%)", y="Diameter ratio") +
  scale_x_continuous(expand = c(0,0), breaks = c(40,50,60,70,80), limits = c(35,83)) +
  scale_y_continuous(expand = c(0,0), breaks = c(0.6,0.8,1,1.2,1.4,1.6), limits = c(0.6,1.6)) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_blank(), 
        legend.position = "none",
        legend.background = element_blank(), 
        legend.spacing.y = unit(1, "mm"), 
        plot.title = element_text(size = 28),
        axis.text = element_text(size=26), axis.title = element_text(size = 26), 
        axis.title.x = element_text(size = 26),
        axis.text.x = element_text(size = 26), 
        plot.margin = unit(c(0.5,0.5,0.25,0.45),"cm")) +
  annotate(geom = "text", x=43, y=0.7, label="DR = -0.011*CC + 1.7", color="black", size=10) + 
  annotate(geom = "text", x=55, y=0.7, label=~R^2==~"0.31, p = 0.03", color="black", size=10)


#litter cover
x <- df$littercov
y <- df$meandratio
z <- lm(y ~ x)
summary(z)

fig1d <- ggplot() +
  geom_point(df, mapping=aes(x=littercov, y=meandratio, color=sigd), size=3) +
  geom_smooth(data=df, mapping = aes(x=littercov, y=meandratio), color="black", size=1.2, 
              method = "lm", 
              formula = y~(x), se=FALSE) +
  scale_colour_manual(values = c("steelblue3", "darkgray", "tomato3")) +
  labs(x="Litter cover (%)", y="Diameter ratio") +
  scale_x_continuous(expand = c(0,0), breaks = c(70,80,90,100), limits = c(65,100)) +
  scale_y_continuous(expand = c(0,0), breaks = c(0.6,0.8,1,1.2,1.4,1.6), limits = c(0.6,1.6)) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_blank(), 
        legend.position = "none",
        legend.background = element_blank(), 
        legend.spacing.y = unit(1, "mm"), 
        plot.title = element_text(size = 28),
        axis.text = element_text(size=26), axis.title = element_text(size = 26), axis.title.x = element_text(size = 26),
        axis.text.x = element_text(size = 26), 
        plot.margin = unit(c(0.5,0.75,0.25,0.45),"cm")) +
  annotate(geom = "text", x=71, y=0.7, label="DR = -0.016*LC + 2.46", color="black", size=10) + 
  annotate(geom = "text", x=80, y=0.7, label=~R^2==~"0.38, p = 0.02", color="black", size=10)


plot_grid(fig1a,fig1b,fig1c,fig1d, nrow = 4, #rel_heights = c(0.95,0.77,1.1),
          labels = c("a", "b", "c", "d"),
          label_x = 0.95, label_y = 0.95, label_size = 40)


#sig diff within plots###############################################
#AZFG 2 vs 3
ponratioazfg <- subset(ponratio, core == "AZFG")
ponratioazfg2 <- subset(ponratioazfg, plotnum == 2)
ponratioazfg3 <- subset(ponratioazfg, plotnum == 3)
ponratioazfg <- as.data.frame(rbind(ponratioazfg2,ponratioazfg3))

#clay, sand, ba, litter cov possibly, canopy cover, herb
#clay
fig3a <- ggplot() +
  geom_hline(yintercept = 1, linetype=2) +
  geom_boxplot(ponratioazfg, mapping= aes(x=clay, y=dratio, group=plot), width=1.8, 
               fill="salmon", outlier.shape = NA) +
  stat_boxplot(ponratioazfg, mapping= aes(x=clay, y=dratio, group=plot), width=0.9,  
               geom = "errorbar") +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x = "Clay (%)", y="Diameter ratio", title="AZFG") +
  scale_x_continuous(expand = c(0,0), breaks = c(12,16,20,24), limits = c(12,24)) +
  scale_y_continuous(expand = c(0,0), breaks = c(0,1,2,3), limits = c(0,3)) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size = 30), legend.position = "bottom",
        legend.background = element_blank(), legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 36),
        axis.text = element_text(size = 30), axis.title = element_text(size = 30), 
        axis.text.x = element_text(size = 30, vjust=-0.1), 
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm")) +
  annotate(geom = "text", x=13, y=2.4, label="ab", fontface='italic', color="tomato3", size=12) +
  annotate(geom = "text", x=21, y=1.95, label="cde", fontface='italic', color="darkgray", size=12) 
#annotate(geom = "text", x=6, y=1.6, label="c", fontface='italic', color="steelblue3", size=10) 
#annotate(geom = "text", x=39.5, y=2.8, label="c", fontface='bold', color="black", size=13)

#sand
fig3b <- ggplot() +
  geom_hline(yintercept = 1, linetype=2) +
  geom_boxplot(ponratioazfg, mapping= aes(x=sand, y=dratio, group=plot), width=1.4, 
               fill="salmon", outlier.shape = NA) +
  stat_boxplot(ponratioazfg, mapping= aes(x=sand, y=dratio, group=plot), width=0.7,  
               geom = "errorbar") +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x = "Sand (%)", y="Diameter ratio", title="AZFG") +
  scale_x_continuous(expand = c(0,0), breaks = c(24,26,28,30,32), limits = c(24,32)) +
  scale_y_continuous(expand = c(0,0), breaks = c(0,1,2,3), limits = c(0,3)) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size = 30), legend.position = "bottom",
        legend.background = element_blank(), legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 36),
        axis.text = element_text(size = 30), axis.title = element_text(size = 30), axis.title.y = element_blank(),
        axis.text.x = element_text(size = 30, vjust=-0.1), 
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm")) +
  annotate(geom = "text", x=31, y=2.4, label="ab", fontface='italic', color="tomato3", size=12) +
  annotate(geom = "text", x=25, y=1.95, label="cde", fontface='italic', color="darkgray", size=12) 
#annotate(geom = "text", x=6, y=1.6, label="c", fontface='italic', color="steelblue3", size=10) 
#annotate(geom = "text", x=39.5, y=2.8, label="c", fontface='bold', color="black", size=13)


#ba
fig3c <- ggplot() +
  geom_hline(yintercept = 1, linetype=2) +
  geom_boxplot(ponratioazfg, mapping= aes(x=ba, y=dratio, group=plot), width=2.6, 
               fill="salmon", outlier.shape = NA) +
  stat_boxplot(ponratioazfg, mapping= aes(x=ba, y=dratio, group=plot), width=1.3,  
               geom = "errorbar") +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x = expression(paste("Basal area ", (m^2 * ha^-1))), y="Diameter ratio", 
       title="AZFG") +
  scale_x_continuous(expand = c(0,0), breaks = c(15,20,25,30), limits = c(15,30)) +
  scale_y_continuous(expand = c(0,0), breaks = c(0,1,2,3), limits = c(0,3)) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size = 30), legend.position = "bottom",
        legend.background = element_blank(), legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 36),
        axis.text = element_text(size = 30), axis.title = element_text(size = 30), 
        axis.text.x = element_text(size = 30, vjust=-0.1), 
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm")) +
  annotate(geom = "text", x=26, y=2.4, label="ab", fontface='italic', color="tomato3", size=12) +
  annotate(geom = "text", x=18, y=1.95, label="cde", fontface='italic', color="darkgray", size=12) 
#annotate(geom = "text", x=6, y=1.6, label="c", fontface='italic', color="steelblue3", size=10) 
#annotate(geom = "text", x=39.5, y=2.8, label="c", fontface='bold', color="black", size=13)


#canopy cover
fig3d <- ggplot() +
  geom_hline(yintercept = 1, linetype=2) +
  geom_boxplot(ponratioazfg, mapping= aes(x=canopy, y=dratio, group=plot), width=5, 
               fill="salmon", outlier.shape = NA) +
  stat_boxplot(ponratioazfg, mapping= aes(x=canopy, y=dratio, group=plot), width=2.5,  
               geom = "errorbar") +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x = "Canopy cover (%)", y="Diameter ratio", 
       title="AZFG") +
  scale_x_continuous(expand = c(0,0), breaks = c(40,50,60,70), limits = c(40,70)) +
  scale_y_continuous(expand = c(0,0), breaks = c(0,1,2,3), limits = c(0,3)) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size = 30), legend.position = "bottom",
        legend.background = element_blank(), legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 36),
        axis.text = element_text(size = 30), axis.title = element_text(size = 30), axis.title.y = element_blank(),
        axis.text.x = element_text(size = 30, vjust=-0.1), 
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm")) +
  annotate(geom = "text", x=45.4, y=2.4, label="ab", fontface='italic', color="tomato3", size=12) +
  annotate(geom = "text", x=66.4, y=1.95, label="cde", fontface='italic', color="darkgray", size=12) 
#annotate(geom = "text", x=6, y=1.6, label="c", fontface='italic', color="steelblue3", size=10) 
#annotate(geom = "text", x=39.5, y=2.8, label="c", fontface='bold', color="black", size=13)


#litter cover
fig3e <- ggplot() +
  geom_hline(yintercept = 1, linetype=2) +
  geom_boxplot(ponratioazfg, mapping= aes(x=littercov, y=dratio, group=plot), width=5, 
               fill="salmon", outlier.shape = NA) +
  stat_boxplot(ponratioazfg, mapping= aes(x=littercov, y=dratio, group=plot), width=2.5,  
               geom = "errorbar") +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x = "Litter cover (%)", y="Diameter ratio", 
       title="AZFG") +
  scale_x_continuous(expand = c(0,0), breaks = c(70,80,90,100), limits = c(70,100)) +
  scale_y_continuous(expand = c(0,0), breaks = c(0,1,2,3), limits = c(0,3)) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size = 30), legend.position = "bottom",
        legend.background = element_blank(), legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 36),
        axis.text = element_text(size = 30), axis.title = element_text(size = 30), 
        axis.text.x = element_text(size = 30, vjust=-0.1), 
        plot.margin = unit(c(0.5,0.85,0.5,0.5),"cm")) +
  annotate(geom = "text", x=77.75, y=2.4, label="ab", fontface='italic', color="tomato3", size=12) +
  annotate(geom = "text", x=92, y=1.95, label="cde", fontface='italic', color="darkgray", size=12) 
#annotate(geom = "text", x=6, y=1.6, label="c", fontface='italic', color="steelblue3", size=10) 
#annotate(geom = "text", x=39.5, y=2.8, label="c", fontface='bold', color="black", size=13)


#herbaceous cover
fig3f <- ggplot() +
  geom_hline(yintercept = 1, linetype=2) +
  geom_boxplot(ponratioazfg, mapping= aes(x=herb, y=dratio, group=plot), width=9, 
               fill="salmon", outlier.shape = NA) +
  stat_boxplot(ponratioazfg, mapping= aes(x=herb, y=dratio, group=plot), width=4.5,  
               geom = "errorbar") +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x = "Herbaceous cover (%)", y="Diameter ratio", 
       title="AZFG") +
  scale_x_continuous(expand = c(0,0), breaks = c(0,20,40,60), limits = c(0,60)) +
  scale_y_continuous(expand = c(0,0), breaks = c(0,1,2,3), limits = c(0,3)) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size = 30), legend.position = "bottom",
        legend.background = element_blank(), legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 36),
        axis.text = element_text(size = 30), axis.title = element_text(size = 30), axis.title.y = element_blank(),
        axis.text.x = element_text(size = 30, vjust=-0.1),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm")) +
  annotate(geom = "text", x=50.62, y=2.4, label="ab", fontface='italic', color="tomato3", size=12) +
  annotate(geom = "text", x=5.5, y=1.95, label="cde", fontface='italic', color="darkgray", size=12) 
#annotate(geom = "text", x=6, y=1.6, label="c", fontface='italic', color="steelblue3", size=10) 
#annotate(geom = "text", x=39.5, y=2.8, label="c", fontface='bold', color="black", size=13)




#COMT 1 vs 2, 3
ponratiocomt <- subset(ponratio, core == "COMT")
ponratiocomt1 <- subset(ponratiocomt, plotnum == 1)
ponratiocomt2 <- subset(ponratiocomt, plotnum == 2)
ponratiocomt3 <- subset(ponratiocomt, plotnum == 3)
ponratiocomt <- as.data.frame(rbind(ponratiocomt1,ponratiocomt2,ponratiocomt3))

#canopy
#ba
fig3g <- ggplot() +
  geom_hline(yintercept = 1, linetype=2) +
  geom_boxplot(ponratiocomt, mapping= aes(x=ba, y=dratio, group=plot), width=3.5, 
               fill="gold3", outlier.shape = NA) +
  stat_boxplot(ponratiocomt, mapping= aes(x=ba, y=dratio, group=plot), width=2,  
               geom = "errorbar") +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x = expression(paste("Basal area ", (m^2 * ha^-1))), y="Diameter ratio", title="COMT") +
  scale_x_continuous(expand = c(0,0), breaks = c(15,20,25,30,35,40), limits = c(15,40)) +
  scale_y_continuous(expand = c(0,0), breaks = c(0,1,2,3), limits = c(0,3)) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size=20), legend.position = "bottom",
        legend.background = element_blank(), legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 36),
        axis.text = element_text(size = 30), axis.title = element_text(size = 30), 
        axis.text.x = element_text(size = 30, vjust=-0.1), 
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm")) +
  annotate(geom = "text", x=18, y=2.7, label="a", fontface='italic', color="tomato3", size=12) +
  annotate(geom = "text", x=34, y=1.7, label="cd", fontface='italic', color="steelblue3", size=12) +
  annotate(geom = "text", x=38, y=1.6, label="d", fontface='italic', color="steelblue3", size=12) 
#annotate(geom = "text", x=89.7, y=2.65, label="b", fontface='bold', color="black", size=13)



fig3h <- ggplot() +
  geom_hline(yintercept = 1, linetype=2) +
  geom_boxplot(ponratiocomt, mapping= aes(x=canopy, y=dratio, group=plot), width=6, 
               fill="gold3", outlier.shape = NA) +
  stat_boxplot(ponratiocomt, mapping= aes(x=canopy, y=dratio, group=plot), width=3,  geom = "errorbar") +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x = "Canopy cover (%)", y="Diameter ratio", title="COMT") +
  scale_x_continuous(expand = c(0,0), breaks = c(50,60,70,80,90), limits = c(50,90)) +
  scale_y_continuous(expand = c(0,0), breaks = c(0,1,2,3), limits = c(0,3)) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size=20), legend.position = "bottom",
        legend.background = element_blank(), legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 36),
        axis.text = element_text(size = 30), axis.title = element_text(size = 30), 
        axis.text.x = element_text(size = 30, vjust=-0.1), axis.title.y = element_blank(),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm")) +
  annotate(geom = "text", x=54.1, y=2.7, label="a", fontface='italic', color="tomato3", size=12) +
  annotate(geom = "text", x=69.7, y=1.7, label="cd", fontface='italic', color="steelblue3", size=12) +
  annotate(geom = "text", x=82.1, y=1.6, label="d", fontface='italic', color="steelblue3", size=12) 
#annotate(geom = "text", x=89.7, y=2.65, label="b", fontface='bold', color="black", size=13)




#litter cover
fig3i <- ggplot() +
  geom_hline(yintercept = 1, linetype=2) +
  geom_boxplot(ponratiocomt, mapping= aes(x=littercov, y=dratio, group=plot), width=3.5, 
               fill="gold3", outlier.shape = NA) +
  stat_boxplot(ponratiocomt, mapping= aes(x=littercov, y=dratio, group=plot), width=1.5,  
               geom = "errorbar") +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x = "Litter cover (%)", y="Diameter ratio", title="COMT") +
  scale_x_continuous(expand = c(0,0), breaks = c(80,90,100), limits = c(80,102)) +
  scale_y_continuous(expand = c(0,0), breaks = c(0,1,2,3), limits = c(0,3)) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size=20), legend.position = "bottom",
        legend.background = element_blank(), legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 36),
        axis.text = element_text(size = 30), axis.title = element_text(size = 30), 
        axis.text.x = element_text(size = 30, vjust=-0.1), axis.title.y = element_text(size = 30),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm")) +
  annotate(geom = "text", x=84, y=2.7, label="a", fontface='italic', color="tomato3", size=12) +
  annotate(geom = "text", x=93.5, y=1.7, label="cd", fontface='italic', color="steelblue3", size=12) +
  annotate(geom = "text", x=99.87, y=1.6, label="d", fontface='italic', color="steelblue3", size=12) 
#annotate(geom = "text", x=89.7, y=2.65, label="b", fontface='bold', color="black", size=13)


#sand
fig3j <- ggplot() +
  geom_hline(yintercept = 1, linetype=2) +
  geom_boxplot(ponratiocomt, mapping= aes(x=sand, y=dratio, group=plot), width=5, 
               fill="gold3", outlier.shape = NA) +
  stat_boxplot(ponratiocomt, mapping= aes(x=sand, y=dratio, group=plot), width=2.5,  
               geom = "errorbar") +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x = "Sand (%)", y="Diameter ratio", title="COMT") +
  scale_x_continuous(expand = c(0,0), breaks = c(60,70,80,90), limits = c(60,90)) +
  scale_y_continuous(expand = c(0,0), breaks = c(0,1,2,3), limits = c(0,3)) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size=20), legend.position = "bottom",
        legend.background = element_blank(), legend.spacing.y = unit(1, "mm"), plot.title = element_text(size = 36),
        axis.text = element_text(size = 30), axis.title = element_text(size = 30), axis.title.y = element_blank(),
        axis.text.x = element_text(size = 30, vjust=-0.1), 
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm")) +
  annotate(geom = "text", x=81, y=2.7, label="a", fontface='italic', color="tomato3", size=12) +
  annotate(geom = "text", x=75, y=1.7, label="cd", fontface='italic', color="steelblue3", size=12) +
  annotate(geom = "text", x=65, y=1.6, label="d", fontface='italic', color="steelblue3", size=12) 
#annotate(geom = "text", x=89.7, y=2.65, label="b", fontface='bold', color="black", size=13)



plot_grid(fig3a,fig3d,fig3e,fig3b,fig3c,fig3f,fig3g,fig3h,fig3i,fig3j, nrow = 5, #rel_heights = c(0.95,0.77,1.1),
          labels = c("a","b","c","d","e","f","g","h","i","j","k"), align = "vh",
          label_x = 0.95, label_y = 0.94, label_size = 46)


########################################################################################################

