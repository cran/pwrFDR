### R code from vignette source 'pwrFDR-vignette.Rnw'

###################################################
### code chunk number 1: header
###################################################
options(keep.source = TRUE, width = 90, digits=4)
pwrFDRpd <- packageDescription("pwrFDR")
pct <- function(x, dig=2) format(round(10^dig*100*x)/10^dig) %,% "\\\\%"
one.in.k <-
function(x, out="c")
{
  nmch <- c("two","three","four","five","six","seven","eight","nine","ten")
  kk <- 2:10
  y <- kk*x
  err <- abs(y-round(y))
  idx <- which(err==min(err))
  rslt <- kk[idx]
  if(out=="c") rslt <- nmch[idx]
  rslt
}


###################################################
### code chunk number 2: seqtstdat
###################################################
SeqTstDat10 <- 
    as.data.frame(t(matrix(
        c(1,3.93066049143225,0.000166075536538912,0.005,1,
          1,2.74360567522164,0.00733385249064922,0.01,1,
          1,2.41434752771541,0.0177875327883759,0.015,0,
          1,2.38482062693501,0.0191857939241782,0.02,1,
          1,2.23237296410553,0.028071740070398,0.025,0,
          0,-1.90361989911666,0.0601552703570536,0.03,0,
          0,1.12448792534619,0.263796263653014,0.035,0,
          0,-1.00660629048567,0.316822712495861,0.04,0,
          0,0.932733736292402,0.353453014154629,0.045,0,
          0,0.90140806501516,0.369777327548002,0.05,0), 5, 10)))
names(SeqTstDat10) <- c("xi", "X", "P", "BHFDR.crit", "Marked")


###################################################
### code chunk number 3: load-libs
###################################################
library(pwrFDR)
library(ggplot2)
library(TableMonster)


###################################################
### code chunk number 4: show-SeqTstDat10
###################################################
basic.tmPrint(SeqTstDat10, lbl="tbl:SeqTstDat10", cptn="Sorted p-values, their test statistics, " %,%
              "population indicators, and BH-FDR threshold in 10 simultaneous tests")


###################################################
### code chunk number 5: vinfo
###################################################
rhome <- Sys.getenv("R_HOME")
sep<-c("\\\\", "/")[1+(substring(rhome, 1, 1)=="/")]
rscript.path <- rhome %,% sep %,% "site-library" %,% sep %,% "pwrFDR" %,% sep %,% "doc" %,% sep %,% "pwrFDR-vignette.R"


###################################################
### code chunk number 6: avgpwr_minf
###################################################
ss.fdr.r05 <- pwrFDR(effect.size=0.79, alpha=0.15, r.1=0.05, average.power=0.80)


###################################################
### code chunk number 7: solve-for
###################################################
find.avp <- pwrFDR(effect.size=0.79, alpha=0.15, n.sample=ss.fdr.r05$n.sample, r.1=0.05)
find.ss <- pwrFDR(effect.size=0.79, alpha=0.15, r.1=0.05, average.power=0.80)
find.es <- pwrFDR(alpha=0.15, n.sample=ss.fdr.r05$n.sample, r.1=0.05, average.power=0.80)
find.r1 <- pwrFDR(effect.size=0.79, alpha=0.15, n.sample=ss.fdr.r05$n.sample, average.power=0.80)


###################################################
### code chunk number 8: avgpwr_r10_minf
###################################################
ss.fdr.r10 <- update(ss.fdr.r05, r.1=0.10)


###################################################
### code chunk number 9: avgpwr-tbl
###################################################
print(ss.fdr.r05, label="tbl:minf", result="tex", cptn="$m=\\infty$")


###################################################
### code chunk number 10: avgpwr-tbl-join
###################################################
print(join.tbl(ss.fdr.r05, ss.fdr.r10), label="tbl:minf-r05-r10",
               result="tex", cptn="$m=\\infty, r_1=0.05, 0.10$")


###################################################
### code chunk number 11: avgpwr-sim
###################################################
avgpwr.fdr.sim.r10.m1e5 <- pwrFDR(effect.size=0.79, alpha=0.15, r.1=0.10,
                                  n.sample=ss.fdr.r10$n.sample, N.tests=10000,
                                  meth="sim")

avgpwr.fdr.sim.r10.m1e3 <- update(avgpwr.fdr.sim.r10.m1e5, N.tests=1000)

avgpwr.fdr.sim.r10.m100 <- update(avgpwr.fdr.sim.r10.m1e5, N.tests=100)


###################################################
### code chunk number 12: genplotVoRr10m100
###################################################
print(join.tbl(avgpwr.fdr.sim.r10.m1e5, avgpwr.fdr.sim.r10.m1e3, avgpwr.fdr.sim.r10.m100), label="tbl:sim_cmpst",
      result="tex", cptn="Results of simulation calls with varying `m'.")


###################################################
### code chunk number 13: cmpt-ss
###################################################
ss <- pwrFDR(effect.size = 0.79, average.power=0.80, r.1 = 0.10, alpha = 0.15)
avgp <- update(ss, average.power=NULL, n.sample=ss$n.sample)


###################################################
### code chunk number 14: gen-FDP-violin-plot
###################################################
sim <- list()
FDX20.tbl <- TPX70.tbl <- NULL
m.lst <- list(m10000=10000, m2000=2000, m1000=1000, m500=500, m250=250, m100=100)
n.mlst <- length(m.lst)
for(k in 1:n.mlst)
{
  sim[[names(m.lst)[k]]] <- update(avgp, method = "sim", N.tests=m.lst[[k]])
  FDX20.tbl <- c(FDX20.tbl, with(detail(sim[[names(m.lst)[k]]])$reps, mean((R-T)%over%R>=0.20)))
  TPX70.tbl <- c(TPX70.tbl, with(detail(sim[[names(m.lst)[k]]])$reps, mean( T %over% M1 < 0.70)))
}
FDX20.tbl <- as.data.frame(rbind(FDX20.tbl))
names(FDX20.tbl) <- m.lst
FDX20.tbl <- format(FDX20.tbl)

TPX70.tbl <- as.data.frame(rbind(TPX70.tbl))
names(TPX70.tbl) <- m.lst
TPX70.tbl <- format(TPX70.tbl)

FDX20.TPX70.tbl <- rbind(FDX20.tbl, TPX70.tbl)
FDX20.TPX70.tbl <- cbind(` `=c("$P(\\mrm{FDP} \\geq 0.20)$", "$P(\\mrm{TPP} < 0.70)$"), FDX20.TPX70.tbl)
    
n.sim <- nrow(detail(sim[["m10000"]])$reps)
FDP.dat <- TPP.dat <- NULL
for(k in 1:length(m.lst))
{
    FDP.dat <- c(FDP.dat, with(detail(sim[[names(m.lst)[k]]])$reps, (R-T)%over%R))
    TPP.dat <- c(TPP.dat, with(detail(sim[[names(m.lst)[k]]])$reps, T %over% M1))
}

FDP.dat <- data.frame(FDP=FDP.dat)
FDP.dat$group <- rep(unlist(m.lst), each=n.sim)
FDP.dat$group <- ordered(FDP.dat$group, levels=FDP.dat$group[1000*(0:(n.mlst-1))+1])
FDP.plot <- ggplot()+geom_violin(position="dodge", alpha=0.5, aes(y=FDP, x=group), data=FDP.dat) + ylim(c(0, 0.4))

TPP.dat <- data.frame(TPP=TPP.dat)
TPP.dat$group <- rep(unlist(m.lst), each=n.sim)
TPP.dat$group <- ordered(TPP.dat$group, levels=TPP.dat$group[1000*(0:(n.mlst-1))+1])
TPP.plot <- ggplot()+geom_violin(position="dodge", alpha=0.5, aes(y=TPP, x=group), data=TPP.dat) + ylim(c(0.6, 1))


###################################################
### code chunk number 15: prnt.excdnc
###################################################
basic.tmPrint(FDX20.TPX70.tbl, cptn="Simulation estimates of indicated probabilities of FDX and TPX for indicated " %,%
                                    "values of $m$ when $\\alpha=0.15, r_1=0.20$, effect size 0.79 with average " %,%
                                    "power 80\\%", lbl="tbl:FDX20")


###################################################
### code chunk number 16: plot-FDP-violin-plot
###################################################
FDP.plot


###################################################
### code chunk number 17: plot-TPP-violin-plot
###################################################
TPP.plot


###################################################
### code chunk number 18: cmpt-ss-Rom-Inf
###################################################
ss.Rom <- pwrFDR(effect.size = 0.79, average.power=0.80, r.1 = 0.20, alpha = 0.15,
                 FDP.control.method="Romano")


###################################################
### code chunk number 19: cmpt-avgp-Rom-Inf
###################################################
avgp.Rom <- pwrFDR(effect.size = 0.79, n.sample=ss.Rom$n.sample, r.1 = 0.20, alpha = 0.15,
                   FDP.control.method="Romano")
avgp.Rom.500.sim <- update(avgp.Rom, N.tests=500, method="sim")
avgp.Rom.250.sim <- update(avgp.Rom, N.tests=250, method="sim")
avgp.Rom.100.sim <- update(avgp.Rom, N.tests=100, method="sim")


###################################################
### code chunk number 20: prt-Rom
###################################################
print(join.tbl(avgp.Rom.500.sim, avgp.Rom.250.sim, avgp.Rom.100.sim), result="tex", cptn=" ", label=" ")


###################################################
### code chunk number 21: cmpt-BHFDX
###################################################
ss.BHFDX.500 <- pwrFDR(effect.size = 0.79, average.power=0.80, r.1 = 0.20, alpha = 0.15,
                       FDP.control.method="BHFDX", N.tests=500)
ss.BHFDX.250 <- update(ss.BHFDX.500, N.tests=250)
ss.BHFDX.100 <- update(ss.BHFDX.500, N.tests=100)

avgp.BHFDX.500.sim <- update(ss.BHFDX.500, n.sample=ss.BHFDX.500$n.sample, average.power=NULL,
                             method="sim")
avgp.BHFDX.250.sim <- update(ss.BHFDX.250, n.sample=ss.BHFDX.250$n.sample, average.power=NULL,
                             method="sim")
avgp.BHFDX.100.sim <- update(ss.BHFDX.100, n.sample=ss.BHFDX.100$n.sample, average.power=NULL,
                             method="sim")


###################################################
### code chunk number 22: prt-BHFDX
###################################################
print(join.tbl(avgp.BHFDX.500.sim, avgp.BHFDX.250.sim, avgp.BHFDX.100.sim), result="tex", cptn=" ", label=" ")


