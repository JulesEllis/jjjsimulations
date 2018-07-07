Produces tables and figures of paper.

How to generate tables and figures of the paper:

1. PRODUCE FWER TABLE

sizes1<-c(2,3,4,5,6, 7,8,9,10,15,20,25,50,75, 100, 1000, 10000, 100000,1000000)
sizes2<-c(2,3,4,5,6, 7,8,9,10,15,20,25,50,75, 100, 1000)

set.seed(2019)
P1<-Study_FWER("P1", sizes2, 100, 10000)
set.seed(2020)
P2<-Study_FWER("P2", sizes1, 100, 10000)
set.seed(2021)
P3<-Study_FWER("P3", sizes1, 100, 10000)
set.seed(2022)
N1<-Study_FWER("N1", sizes2, 100, 10000)
set.seed(2023)
N2<-Study_FWER("N2", sizes1, 100, 10000)
set.seed(2024)
N3<-Study_FWER("N3", sizes1, 100, 10000)

2. PRODUCE POWER TABLE

Study_RelativePower()

3. PRODUCE POWER FIGURES

set.seed(2025)
nsim <- 10000
a<- Study_Power_Fig2 (nsim)
a<-t(a)
a<-data.frame(a)
colnames(a)<-c("ncp0","ncp1", "sigma0","sigma1","n0", "n1", "Bonferroni","FGS","Fisher","Tippet","Iplus","LR","Conditional Bonferroni","Conditional FGS","Conditional Fisher","Conditional Tippet","Conditional Iplus","Conditional LR")
Fig2 <-a
write.table(Fig2,file="Fig2.txt")
save(Fig2,file="Fig2.Rda")

set.seed(2027)
nsim = 10000
a<-Study_FDR(nsim)
dim(a) <- c(4, 20)
rownames(a) <- c("fdru", "tpru", "fdrc", "tprc")
a<-t(a) 
a<- data.frame(a)
s <-1:20
FigFdr<-cbind(s,a)

a<- cbind(Fig2, FigFdr)


plot(a$n0, a$Bonferroni, type = "b", xlab = "Number of true hypotheses", ylab = "Power or TPR", lty=1, pch=1, ylim = c(0,1))
title(main ="Unconditionalized")
points(a$n0, a$FGS, type = "b", pch=2, lty=1)
points(a$n0, a$Fisher, type = "b", pch=3, lty=1)
points(a$n0, a$Tippet, type = "b", pch=4, lty=1)
points(a$n0, a$Iplus, type = "b", pch=5, lty=1)
points(a$n0, a$LR, type = "b", pch=6, lty=1)
points(a$n0, a$tpru, type = "b",pch=19, bg=1, lty=1)

plot(a$n0, a$`Conditional Bonferroni`, type = "b", xlab = "Number of true hypotheses", ylab = "Power or TPR", lty=1, pch=1, ylim = c(0,1))
title(main ="Conditionalized")
points(a$n0, a$`Conditional FGS`, type = "b",pch=2, lty=1)
points(a$n0, a$`Conditional Fisher`, type = "b",pch=3, lty=1)
points(a$n0, a$`Conditional Tippet`, type = "b",pch=4, lty=1)
points(a$n0, a$`Conditional Iplus`, type = "b",pch=5, lty=1)
points(a$n0, a$`Conditional LR`, type = "b",pch=6, lty=1)
points(a$n0, a$tprc, type = "b",pch=19, lty=1)

legend("bottomleft", inset=.02, legend = c("Bonferroni", "FGS", "Fisher", "Tippet", "Iplus", "LR", "BH"),lty = 1,pch=c(1:6,19), box.lty=1, cex=0.8, ncol=2)



set.seed(2026)
nsim <- 10000
a<- Study_Power_Fig3 (nsim)
a<-t(a)
a<-data.frame(a)
colnames(a)<-c("ncp0","ncp1", "sigma0","sigma1","n0", "n1", "Bonferroni","FGS","Fisher","Tippet","Iplus","LR","Conditional Bonferroni","Conditional FGS","Conditional Fisher","Conditional Tippet","Conditional Iplus","Conditional LR")
Fig3 <-a
write.table(Fig3,file="Fig3.txt")
save(Fig3,file="Fig3.Rda")

plot(a$n0, a$Bonferroni, type = "b", xlab = "Percentage of true hypotheses", ylab = "Power", lty=1, ylim = c(0,1))
points(a$n0, a$FGS, type = "b", pch=2, lty=1)
points(a$n0, a$`Conditional Bonferroni`, type = "b", pch=3, lty=1)
points(a$n0, a$`Conditional FGS`, type = "b",pch=4, lty=1)
legend("bottomleft", inset=.02, legend = c("Bonferroni", "FGS", "Conditional Bonferroni", "Conditional FGS"),lty = 1,pch=1:4, box.lty=1, cex=0.8)

4. COMPUTE MINIMAX VALUE OF LAMBDA

StudyMiniMax(100, 0.5, 1.0, 100, 10, 100, 10, 500, 0.5)

5. PRODUCE FWER PLOTS (RUN SECTION 1 FIRST)
 
plot(log(P1$m,2),P1$rej/nsim2, xlab = "log(size,2)", ylab= "FWER", ylim = c(0,0.1))
title(main ="Model P1")
abline(h = 0.05, lty = 2)

plot(log(P2$m,2),P2$rej/nsim2, xlab = "log(size,2)", ylab= "FWER", ylim = c(0,0.1))
title(main ="Model P2")
abline(h = 0.05, lty = 2)

plot(log(P3$m,2),P3$rej/nsim2, xlab = "log(size,2)", ylab= "FWER", ylim = c(0,0.1))
title(main ="Model P3")
abline(h = 0.05, lty = 2)

plot(log(N1$m,2),N1$rej/nsim2, xlab = "log(size,2)", ylab= "FWER", ylim = c(0,0.1))
title(main ="Model N1")
abline(h = 0.05, lty = 2)

plot(log(N2$m,2),N2$rej/nsim2, xlab = "log(size,2)", ylab= "FWER", ylim = c(0,0.1))
title(main ="Model N2")
abline(h = 0.05, lty = 2)

plot(log(N3$m,2),N3$rej/nsim2, xlab = "log(size,2)", ylab= "FWER", ylim = c(0,0.1))
title(main ="Model N3")
abline(h = 0.05, lty = 2)
plot(P1$minr,P1$rej/nsim2, xlab = "minimum correlation", ylab= "FWER", ylim = c(0,0.1))
title(main ="Model P1")
abline(h = 0.05, lty = 2)
plot(P2$a/(P2$a + P2$b),P2$rej/nsim2, xlab = "mean autocorrelation", ylab= "FWER", ylim = c(0,0.1))
title(main ="Model P2")
abline(h = 0.05, lty = 2)

plot(P3$minr,P3$rej/nsim2, xlab = "minimum correlation", ylab= "FWER", ylim = c(0,0.1))
title(main ="Model P3")
abline(h = 0.05, lty = 2)

plot(N1$minr,N1$rej/nsim2, xlab = "minimum correlation", ylab= "FWER", ylim = c(0,0.1))
title(main ="Model N1")
abline(h = 0.05, lty = 2)
plot(N2$a/(N2$a + N2$b),N2$rej/nsim2, xlab = "mean absolute autocorrelation", ylab= "FWER", ylim = c(0,0.1))
title(main ="Model N2")
abline(h = 0.05, lty = 2)


plot(N3$minr,N3$rej/nsim2, xlab = "correlation", ylab= "FWER", ylim = c(0,0.1))
title(main ="Model N3")
abline(h = 0.05, lty = 2)

