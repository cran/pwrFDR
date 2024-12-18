### R code from vignette source 'pwrFDR-vignette.Rnw'

###################################################
### code chunk number 1: pwrFDRpd
###################################################
options(keep.source = TRUE, width = 90, digits=4)
pwrFDRpd <- packageDescription("pwrFDR")
pct <- function(x, dig=2) format(ceiling(10^dig*100*x)/10^dig) %,% "\\\\%"
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
### code chunk number 2: load-libs
###################################################
library(pwrFDR)
library(ggplot2)
library(TableMonster)


###################################################
### code chunk number 3: vinfo
###################################################
rhome <- Sys.getenv("R_HOME")
sep<-c("\\\\", "/")[1+(substring(rhome, 1, 1)=="/")]
rscript.path <- rhome %,% sep %,% "site-library" %,% sep %,% "pwrFDR" %,% sep %,% "doc" %,% sep %,% "pwrFDR-vignette.R"


###################################################
### code chunk number 4: avgpwr_minf
###################################################
avgpwr.fdr.r05 <- pwrFDR(effect.size=0.79, alpha=0.15, r.1=0.05, average.power=0.80)


###################################################
### code chunk number 5: avgpwr_r10_minf
###################################################
avgpwr.fdr.r10 <- update(avgpwr.fdr.r05, r.1=0.10)


###################################################
### code chunk number 6: avgpwr-tbl
###################################################
print(avgpwr.fdr.r05, label="tbl:minf", result="tex", cptn="$m=\\infty$")


###################################################
### code chunk number 7: avgpwr-tbl-join
###################################################
print(join.tbl(avgpwr.fdr.r05, avgpwr.fdr.r10), label="tbl:minf-r05-r10",
               result="tex", cptn="$m=\\infty, r_1=0.05, 0.10$")


###################################################
### code chunk number 8: avgpwr-sim
###################################################
avgpwr.fdr.sim.r05.m1e5 <- pwrFDR(effect.size=0.79, alpha=0.15, r.1=0.05,
                                n.sample=avgpwr.fdr.r05$n.sample, N.tests=10000,
                                  meth="sim")

avgpwr.fdr.sim.r05.m1e3 <- update(avgpwr.fdr.sim.r05.m1e5, N.tests=1000)

avgpwr.fdr.sim.r05.m100 <- update(avgpwr.fdr.sim.r05.m1e5, N.tests=100)

avgpwr.fdr.r10.sim.m100 <- update(avgpwr.fdr.r10, n.sample=avgpwr.fdr.r10$n.sample,
                                  average.power=NULL, method="sim", N.tests=100)


###################################################
### code chunk number 9: genplotVoRr10m100
###################################################
print(join.tbl(avgpwr.fdr.sim.r05.m1e5, avgpwr.fdr.sim.r05.m1e3, avgpwr.fdr.sim.r05.m100, avgpwr.fdr.r10.sim.m100), 
      label="tbl:sim_cmpst", result="tex", cptn="Results of simulation calls with varying `m' and `r.1'.")


###################################################
### code chunk number 10: gen-FDP-violin-plot
###################################################
sim <- list()
m.lst <- list(m10000=10000, m2000=2000, m1000=1000, m500=500, m100=100)
for(k in 1:length(m.lst))
sim[[names(m.lst)[k]]] <- pwrFDR(effect.size = 0.79, n.sample = 47, r.1 = 0.20, alpha = 0.15, method = "sim", N.tests=m.lst[[k]])

n.sim <- nrow(detail(sim[["m10000"]])$reps)
FDP.dat <- TPP.dat <- NULL
for(k in 1:length(m.lst))
{
    FDP.dat <- c(FDP.dat, with(detail(sim[[names(m.lst)[k]]])$reps, (R-T)%over%R))
    TPP.dat <- c(TPP.dat, with(detail(sim[[names(m.lst)[k]]])$reps, T %over% M1))
}

FDP.dat <- data.frame(FDP=FDP.dat)
FDP.dat$group <- rep(unlist(m.lst), each=n.sim)
FDP.dat$group <- ordered(FDP.dat$group, levels=FDP.dat$group[1000*(0:4)+1])
FDP.plot <- ggplot()+geom_violin(position="dodge", alpha=0.5, aes(y=FDP, x=group), data=FDP.dat) + ylim(c(0, 0.4))

TPP.dat <- data.frame(TPP=TPP.dat)
TPP.dat$group <- rep(unlist(m.lst), each=n.sim)
TPP.dat$group <- ordered(TPP.dat$group, levels=TPP.dat$group[1000*(0:4)+1])
TPP.plot <- ggplot()+geom_violin(position="dodge", alpha=0.5, aes(y=TPP, x=group), data=TPP.dat) + ylim(c(0.8, 1))


###################################################
### code chunk number 11: plot-FDP-violin-plot
###################################################
FDP.plot


###################################################
### code chunk number 12: plot-TPP-violin-plot
###################################################
TPP.plot


