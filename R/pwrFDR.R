"pwrFDR" <- function(effect.size, n.sample, r.1, alpha, delta=NULL, groups=2, N.tests,
         average.power, TPX.power, lambda, type=c("paired","balanced","unbalanced"),
         grpj.per.grp1=NULL, corr.struct=list(type=c("CS-Blocks", "Toeplitz-Blocks"), block.size=NA, rho=NA),
         FDP.control.method=c("BHFDR","BHFDX","Romano","Auto","both","Holm","Hochberg","Bonferroni"), distopt=1,
         method=c("Theoretical", "simulation"), n.sim=1000, temp.file,
         control=list(tol=1e-8, max.iter=c(1000,20), sim.level=2,
             low.power.stop=TRUE, FDP.meth.thresh=FDP.cntl.mth.thrsh.def,
                      ast.le.a=TRUE, verb=FALSE, show.footer=FALSE))
{
    
    .call. <- m <- match.call()
    pfx <- as.character(m[[1]])
    
    frmls <- names(formals(pwrFDR))
    sppld <- names(m)[-1]

    m <- evald.call <- args.def.err.chk(m, sppld, frmls)
    err.msg <- attr(m, "err.msg")
    
    if(!is.null(err.msg)) stop(err.msg)
    
    sfx <- pwrfdrmthd.sfx[m$method]
    m[[1]] <- as.name(pfx %,% "." %,% sfx)

    m$method <- m$temp.file <- NULL 

    args.full <- evald.call
    args.full[[1]] <- as.name("list")
    args.full <- eval(args.full)
    args.full$call <- evald.call

    verb <- args.full$control$verb
    if(verb>=3) browser()

    out <- eval(m)

    if(verb>=2) browser()

    fdpcmthd.spcd <- attr(out,"passback")$FDP.control.method
    fdpcmthd.spcd <- ifelse(is.null(fdpcmthd.spcd), "BHFDR", fdpcmthd.spcd)
    FDP.control.method <- ifelse(sfx=="th", fdpcmthd.spcd, args.full$FDP.control.method)

    dtl <- attr(out, "detail")
    
    if(verb>=1) browser()
    
    msng <- as.list(args.full$control$is.msng)
    
    ## add N.tests to the position 5 of the output list if available
    N.tests <- ifelse(msng$N.tests, Inf, args.full$N.tests)
    out <- append(out, list(N.tests=N.tests), after=4)
    
    ## when calculating sample size, effect.size, r.1 or alpha, then we must add
    ## n.sample, effect.size, r.1 and alpha manually to the head of the output list.
    if(msng$average.power && msng$TPX.power)
        out <- append(out, with(args.full, list(n.sample=n.sample, effect.size=effect.size, r.1=r.1, alpha=alpha)), 0)

    args.out <- args.full
    is.Auto <- !is.null(out$Auto)
 
    is.BHFDX <- FDP.control.method == "BHFDX"
    is.Romano <- FDP.control.method == "Romano"
    is.Bonferroni <- FDP.control.method=="Bonferroni"
    is.Holm <- FDP.control.method=="Holm"
    is.Hochberg <- FDP.control.method=="Hochberg"
    is.pTPX <- !is.null(args.out$lambda)
    is.sim <-  substring(args.out$call$method,1,3)=="sim"
    is.n <- !is.null(args.out$n.sample)
    
    pwr.type <- c(" average power ", " pTPX power ")[is.pTPX+1]
    pwr <- round(100*c(out$average.power, out$TPX.power)[is.pTPX+1],1)%,%"%"
    n.smp <- ifelse(!is.n, out$n.sample, args.out$n.sample)

    cntl.type <- " BH-FDR control"
    if(is.BHFDX)
        cntl.type <-" pFDX control, alpha=delta="%,%(100*round(args.out$delta,4))%,%"%"

     if(is.Romano)
        cntl.type <- " Romano FDX control"

    if(is.Bonferroni)
        cntl.type <-" Bonferroni FWER control"

    if(is.Holm)
        cntl.type <-" Holm FWER control"

    if(is.Hochberg)
        cntl.type <-" Hochberg FWER control"

    footer <- c("","Sample size of "%,% n.smp %,%" required for "%,%pwr.type%,% pwr%,% 
                   " under" %,% cntl.type %,% ".")

    if(is.BHFDX)
      footer <- c(footer, 
              "NOTE:" %,% cntl.type %,% " is guaranteed if BH-FDR is applied at alpha=" %,%
              round(out$alpha.star,4))
    
    attr(out, "arg.vals") <- args.out
    show.footer <- control$show.footer
    if(is.null(show.footer))show.footer <- FALSE
    if(show.footer) attr(out, "footer") <- sprintf("%s",footer)
    class(out) <- "pwrFDR"
    out$call <- .call.
    
    out.sv <- out
    args <- arg.vals(out)
    msng <- as.list(args$control$is.msng)
    is.Ntst <- !msng$N.tests
    is.FDPcnt <- !msng$FDP.control.method
    
    vals.tpx <- NULL
    if(is.pTPX)
    {
        val.tpx <- ifelse(msng$TPX.power, out$TPX.power, args$call$TPX.power)
        val.lmbd <- arg.vals(out)$lambda
        vals.tpx <- list(val.tpx, val.lmbd)
        names(vals.tpx) <- c("TPX.power", "lambda")
    }
    nms.sigma <- list(c("sigma.rtm.Rom","sigma.rtm.VoR", "sigma.rtm.ToM"),
                      c("se.Rom","se.VoR", "se.ToM"))[[1+is.Ntst]]
    
    nms.reord.1 <- c("call", "N.tests","r.1","n.sample","effect.size","alpha")
    
    if(is.sim) nms.reord.3 <- c("emp.FDR", "average.power", "gamma", nms.sigma)
    if(!is.sim) nms.reord.3 <- c("average.power", "gamma", nms.sigma)
    
    ind <- which(c(nms.reord.1,nms.reord.3)==c("alpha","emp.FDR")[1+is.sim])

    out <- out[c(nms.reord.1, nms.reord.3)]
    
    if(is.pTPX)
        out <- append(out, vals.tpx, ind)

    cntl.n.vals <- c(`BHFDR`=1, `BHFDX`=2, `Romano`=3, `Hochberg`=4, `Holm`=5, `Bonferroni`=6)
    
    cntl <- FDP.control.method
    if(cntl=="both") cntl <- "BHFDX"
    cntl.val <- cntl.n.vals[cntl]
    
    ind <- which(names(out)==c("alpha","emp.FDR")[1+is.sim])
    out <- switch(cntl,
                  BHFDR = append(out, list(FDP.cnt=cntl.val), ind),
                  BHFDX = append(out, list(FDP.cnt=cntl.val, delta=args$delta, alpha.star=out.sv$alpha.st), ind),
                  Romano = append(out, list(FDP.cnt=cntl.val), ind),
                  Bonferroni = append(out, list(FDP.cnt=cntl.val), ind),
                  Holm = append(out, list(FDP.cnt=cntl.val), ind),
                  Hochberg = append(out, list(FDP.cnt=cntl.val), ind))

    attr(out,"arg.vals") <- args

    ## Deterimines the content for the main output for the print method
    ## There's only one possibility for all choices and all methods except
    ## the simulation method for which FDP.control.method="both" is an option
    ## which computes both the BHFDX and the Romano methods. Everything is
    ## in the reps component of the detail attribute, but the components of
    ## the result in this case (FDP.control.method="both") will 
    is.BHFDX <- is.BHFDX || (args$FDP.control.method=="both")
    if(is.sim && is.BHFDX)
    {
        out$average.power <- out.sv$average.power.st
        out$emp.FDR <- out.sv$emp.FDR.st
        out$gamma <- out.sv$gamma.st
        out$se.Rom <- out.sv$se.Rom.st
        out$se.VoR <- out.sv$se.VoR.st
        out$se.ToM <- out.sv$se.ToM.st
        out$TPX.power <- out.sv$TPX.power.st
    }
    if(is.sim && (is.Romano || is.Holm || is.Hochberg))
    {
        out$average.power <- out.sv$average.power.R
        out$emp.FDR <- out.sv$emp.FDR.R
        out$gamma <- out.sv$gamma.R
        out$se.Rom <- out.sv$se.Rom.R
        out$se.VoR <- out.sv$se.VoR.R
        out$se.ToM <- out.sv$se.ToM.R
        out$TPX.power <- out.sv$TPX.power.R
    }
    class(out) <- class(out.sv)
    attr(out, "detail") <- dtl
    out
}

"pwrFDR.grid" <-
function(effect.size, n.sample, r.1, alpha, delta, groups, N.tests, average.power,
             TPX.power, lambda, type, grpj.per.grp1, corr.struct, FDP.control.method, distopt, control)
{
    .call. <- m <- fd <- match.call()
    pfx <- as.character(m[[1]])
    frmls <- names(formals(pwrFDR.grid))
    sppld <- names(m)[-1]
    n.sppld <- length(sppld)
    is.msng <- !(frmls %in% sppld)
    names(is.msng) <- frmls
    msng.nms <- names(is.msng)[is.msng]
    eval.env <- sys.parent
    for (k in 1:n.sppld)
    {
      m[[sppld[k]]] <- eval(m[[sppld[k]]], eval.env())
      fd[[sppld[k]]] <- length(m[[sppld[k]]])
    }
    is.pwr <- (!is.msng["average.power"]) || (!is.msng["TPX.power"])
    nera <- c("n.sample", "effect.size", "r.1", "alpha")
   
    # check that either power is unspecified and nera are all specified
    # or that power is unspecified and exactly one nera is missing.

    # are we doing the power specified or power unspecified case.
    pwr.spcfd <- !is.msng["average.power"]||!is.msng["TPX.power"]

    # which power type?
    use.L.pwr <- !is.msng["TPX.power"]
    pwr.nm <- c("average.power", "TPX.power")[1+use.L.pwr]
    msng.nera <- nera[which(is.msng[nera])]

    fd[[1]] <- as.name("factorial.design")

    fd <- as.data.frame(eval(fd))
    n.conds <- nrow(fd)
 
    pwrFDR.call <- as.call(expression(pwrFDR))

    names(fd) <- sppld

    rslts <- list()
    cnds <- NULL
    for(j in 1:n.conds)
    {
      cnds.j <- NULL
      for(k in 1:n.sppld)
      {
          cnds.j <- c(cnds.j, m[[sppld[k]]][fd[j,k]])
          pwrFDR.call[[sppld[k]]] <- m[[sppld[k]]][fd[j,k]]
      }
      cnds <- rbind(cnds, cnds.j)
      rslts[[j]] <- try(eval(pwrFDR.call, eval.env()), silent=TRUE)
    }
    dimnames(cnds)[[1]] <- NULL
    cnds <- as.data.frame(cnds)
    names(cnds) <- sppld
    
    list(conditions=cnds, results=rslts)
}

"pwrFDR.th" <-
function(effect.size, n.sample, r.1, alpha, delta, groups, N.tests,
         average.power, TPX.power, lambda, type, grpj.per.grp1,
         corr.struct, FDP.control.method, distopt, method, n.sim, temp.file, control)
{
    ## pwr, n.sample, effect.size, r.1, alpha
    .call. <- m.sv <- m <- match.call()
    verb <- control$verb
    ee <- function(x)exp(exp(x))-1
    ll <- function(x)log(log(1 + x))
    ijumps <- FALSE
    Auto <- ""
    iter.max <- control$max.iter[1]
    
    is.msng <- control$is.msng

    pfx <- "pwrFDR." %,% ifelse(substring(FDP.control.method, 1, 4) != "Bonf", "FP", "BIN")

    # pfx <- as.character(m[[1]])
    
    is.N <- !is.msng["N.tests"]
    pwr.spcfd <- !is.msng["average.power"]||!is.msng["TPX.power"]

    if(!pwr.spcfd)
    {
      m[[1]] <- as.name(pfx %,% ".1X")
      #  m[[1]] <- as.name(as.character(m[[1]]) %,% ".1X")
      m$average.power <- m$TPX.power <- NULL
      out <- eval(m)
    }
    if(pwr.spcfd)
    {
      use.L.pwr <- !is.msng["TPX.power"]
      pwr.nm <- c("average.power", "TPX.power")[1+use.L.pwr]
      msng.nera <- nera[which(is.msng[nera])]
      ## 'tau' is a list of functions with one component for each of the actual parameters
      ## n.sample, effect.size, r.1 and alpha, mapping each of these actaul parameter to the real line
      ## 'tau.inv' is a list of functions, each component being the corresponding inverse of that
      ## component the list of functions, 'tau'.
      tau.inv <- list(n.sample=ee, effect.size=ee, r.1=logitInv, alpha=logitInv)[[msng.nera]]
      tau <- list(n.sample=ll, effect.size=ll, r.1=logit, alpha=logit)[[msng.nera]]
      if(control$verb>=3) browser()
      OBJ <-
      function(x, m, misc)
      {
        rslt.call <- m
        rslt.call[[1]] <- as.name(pfx %,% ".1X")
        rslt.call[[misc$msng.nera]] <- misc$tau.inv(x)

        rslt.call$average.power <- rslt.call$TPX.power <- NULL
        rslt <- suppressWarnings(eval(rslt.call))

        if(!use.L.pwr)
        {
          pwr <- rslt$pwr <- rslt$average.power
          pwr. <- max(min(pwr, 1-1e-16),1e-16)
          out <- out.n <- logit(pwr.) - logit(average.power)
        }
        if(use.L.pwr)
        {
          pwr <- rslt$pwr <- rslt$TPX.power
          pwr. <- min(max(pwr, 1e-16), 1-1e-16)
          out <- out.n <- logit(pwr.) - logit(TPX.power)
        }        
        attr(out, "detail") <- rslt
        if(control$verb>0)
          cat(sprintf("(obj, %s, %s)=(%g, %g, %g)\n", misc$msng.nera, misc$pwr.nm, out, misc$tau.inv(x), pwr.))
        if(control$verb>0 && is.na(out)) browser()            
        out
      }
      idistopt <- distopt
      is.pos <- (dists$minv[[1+idistopt]]==0)

      if(use.L.pwr) pwr <- TPX.power
      if(!use.L.pwr) pwr <- average.power

      l.nera <- c(n.sample=3,      effect.size=0.01, r.1=0.00001, alpha=0.00001)
      u.nera <- c(n.sample=5000,   effect.size=3,    r.1=0.99999, alpha=0.99999)
      l.x <- tau(l.nera[msng.nera])
      u.x <- tau(u.nera[msng.nera])
      if(msng.nera=="n.sample")
      {
          grps <- groups
          r.0 <- 1 - r.1
          alpha.0 <- r.0*alpha
          if(use.L.pwr) pwr <- TPX.power
          if(!use.L.pwr) pwr <- average.power
          gma <- r.1*pwr/(1-alpha.0)

          done <- FALSE
          while(!done)
          {
              obj.ux <-OBJ(u.x,m=m.sv, misc=list(msng.nera=msng.nera, pwr.nm=pwr.nm, tau.inv=tau.inv))
              done <- !is.na(obj.ux)
              if(!done) u.x <- tau(tau.inv(u.x) -1)
          }
          lu.k <- c(l.x, u.x)
          err <- 1
          x <- u.x
          iter <- 0
          brk.1to2 <- FALSE
          while(err > 1e-3 && iter <= iter.max)
          {
            while(err > 1e-3 && iter < 20)
            {
              brk.1to2 <- FALSE  
              obj <- OBJ(x, m=m.sv, misc=list(msng.nera=msng.nera, pwr.nm=pwr.nm, tau.inv=tau.inv))
              pos <- obj > 0
              err <- abs(obj)
              if(control$verb > 0) cat(sprintf("n: %d, x: %g, lu.k[1]: %g, lu.k[2]: %g, obj: %g, err: %g\n",
                                               ceiling(tau.inv(x)), x, lu.k[1], lu.k[2], obj, err))
              lu.k <- pos*c(lu.k[1], x) + (1-pos)*c(x, lu.k[2])
              x <- mean(lu.k)
              iter <- iter + 1
            }
            if(err > 1e-3 && abs(lu.k[1]-lu.k[2])<1e-5)
            {
              lu.k[2] <- 10*lu.k[2]
              iter <- 0
              brk.1to2 <- TRUE
            }
            if(err > 1e-3 && !brk.1to2)
            {
              ## we got here because the "Auto" option can't make up its mind. In this case just decde based
              ## upon N.tests
              FDP.meth.thresh <- control$FDP.meth.thresh  
              is.cs <- !missing(corr.struct)
              N.eff <- N.tests
              if(is.cs) N.eff <- with(corr.struct, N.tests/(block.size*rho))  
              lu.k <- lu.k <- c(l.x, u.x)
              fdp.mth <- 1*(N.eff > FDP.meth.thresh[2]) + 2*(N.eff <= FDP.meth.thresh[2])
              m.sv$FDP.control.method <- Auto <- c("BHFDX", "Romano")[fdp.mth]
              ijumps <- TRUE
              iter <- 0
            }
          }
          if(err > 1e-3)
            stop("No Solution, " %,% msng.nera %,% " too close to " %,%
                 l.nera[msng.nera] %,% " or " %,% u.nera[msng.nera] %,% " and " %,%
                 pwr.nm %,% " is " %,% signif(attr(obj,"detail")$pwr,3))
          
          ans <- list(f.root=obj)
          y <- ceiling(tau.inv(x))
      }
      if(msng.nera!="n.sample")
      {
          ans <- uniroot(f=OBJ, lower=l.x, upper=u.x, m=m.sv, tol=.Machine$double.eps^0.35,
                          misc=list(msng.nera=msng.nera, pwr.nm=pwr.nm, tau.inv=tau.inv))
          obj <- ans$f.root
          if(abs(obj)>1e-3) stop("No Solution, " %,% msng.nera %,% " too close to " %,%
                                 l.nera[msng.nera] %,% " or " %,% u.nera[msng.nera] %,% " and " %,%
                                 pwr.nm %,% " is " %,% signif(attr(obj,"detail")$pwr,3))
          y <- tau.inv(ans$root)
      }

      ans <- append(ans, attr(obj, "detail"))
      
      assign(msng.nera, y)

      nii.sample <- n.sample
      if(groups>=2 && type==3) nii.sample <- n.sample*grpj.per.grp1

      .DF. <- eval(DF[type])
      .NCP. <- eval(NCP[type])
      
      pars0 <- eval(dists$pars0[[1+idistopt]])
      pars1 <- eval(dists$pars1[[1+idistopt]])

      r.0 <- 1 - r.1
      alpha.0 <- r.0*alpha
      gma <- ifelse(FDP.control.method=="Bonferroni", ans$gamma, ans$r.1*ans$average.power/(1-alpha.0))
      c.g <- dists$qdist[[1+idistopt]](1 - gma*alpha/2^(!is.pos), pars0)
      eIII <- dists$pdist[[1+idistopt]](-c.g, pars1)
      dtl <- attr(ans$f.root, "detail")
      dtl$pwr <- NULL
      out <- numeric(0)
      out[msng.nera] <- as.list(get(msng.nera))
      out[nera[nera!=msng.nera]] <- sapply(nera[nera!=msng.nera], FUN=\(x)get(x))
      if(ijumps) dtl[["Auto"]] <- Auto
      out <- c(out, dtl)
    }
    if(verb>=2) browser()
    out
}

"pwrFDR.FP.1X" <-
function(effect.size, n.sample, r.1, alpha, delta, groups, N.tests,
         average.power, TPX.power, lambda, type, grpj.per.grp1,
         corr.struct, FDP.control.method, distopt, method, n.sim, temp.file, control)
{
    .call. <- m <- m.sv <- match.call()

    use.Hlm.Hch <- FDP.control.method %in% c("Holm","Hochberg")
    if(use.Hlm.Hch)
    {
        FDP.control.method <- "Romano"
        delta <- m$delta <- 1/N.tests*(1 - 1e-3)
    }
        
    is.msng <- control$is.msng
    ast.le.a <- control$ast.le.a
    if(is.msng["alpha"]) .call.$delta <- m$delta <- m.sv$delta <- m$alpha

    is.N <- !missing(N.tests)
    is.cs <- !missing(corr.struct)
    
    idistopt <- distopt
    
    nii.sample <- n.sample
    if(groups>=2 && type==3) nii.sample <- n.sample*grpj.per.grp1

    ## if(groups==1) type is 1 (paired)
    ## if(groups==2) type is 1 (paired), 2 (balanced) or 3 (unbalanced)
    ## if(groups>=3) type is 2 (balanced) or 3 (unbalanced)
    .DF. <- eval(DF[type])
    .NCP. <- eval(NCP[type])
    
    pars0 <- eval(dists$pars0[[1+idistopt]])
    pars1 <- eval(dists$pars1[[1+idistopt]])
    
    is.pos <- (dists$minv[[1+idistopt]]==0)
    
    Auto <- se.by.a <- gma.Ntsts <- NULL
    
    n <- n.sample
    r.0 <- 1 - r.1
    alpha.0 <- r.0*alpha
    alpha.st <- NULL
    do.TPX.power <- !missing(lambda)
    
    FDP.cnt.mthd <- ifelse(FDP.control.method %in% c("Holm","Hochberg"), "Romano", FDP.control.method)

    if(FDP.cnt.mthd == "Auto")
    {        
      ## determine whether to use BHFDR, BHFDX or Romano
      ## compute se.VoR under BHFDR
      m.fxd.pt <- m  
      m.fxd.pt[[1]] <- as.name("CDF.Pval.au.eq.u")
      m.fxd.pt$N.tests <- m.fxd.pt$lambda <- m.fxd.pt$corr.struct <- m.fxd.pt$FDP.control.method <- m.fxd.pt$delta <- NULL

      gma <- eval(m.fxd.pt, sys.parent())$gamma

      c.g <- qdist(1 - gma*m$alpha/2^(!is.pos), pars0)
 
      ## Now call variance functions.
      var.VoR.call <- as.call(expression(var.VoR))
      var.VoR.call$r.1 <- m$r.1
      var.VoR.call$alpha <- m$alpha
      var.VoR.call$delta <- m$delta
      var.VoR.call$gamma <- gma
      var.VoR.call$n.sample <- n.sample
      var.VoR.call$nii.sample <- nii.sample
      var.VoR.call$groups <- groups 
      var.VoR.call$effect.size <- effect.size
      var.VoR.call$type <- type
      var.VoR.call$distopt <- distopt
      var.VoR.call$FDP.control.method <- FDP.control.method

      N.eff <- N.tests
      if(is.cs)
      {
        var.VoR.call$corr.struct <- corr.struct
        N.eff <- with(corr.struct, N.tests/(block.size*rho))
      }
      se.VoR <- (eval(var.VoR.call)/N.tests)^0.5
        
      # based upon se[V/R] / alpha   and N.tests
      FDP.meth.thresh <- control$FDP.meth.thresh
      se.by.a <- se.VoR/alpha

      do.bhfdr <- (se.by.a <= FDP.meth.thresh[1] && N.eff >= FDP.meth.thresh[2])
      do.bhclt <- (se.by.a > FDP.meth.thresh[1] && N.eff >= FDP.meth.thresh[2])
      do.romano <- !do.bhfdr && !do.bhclt
        
      mthd.chc <- 1*do.romano + 2*do.bhclt + 3*do.bhfdr
      FDP.control.method <- FDP.cnt.mthd <- c("Romano", "BHFDX", "BHFDR")[mthd.chc]
      Auto <- FDP.control.method
    }

    if(FDP.cnt.mthd == "BHFDR")
    {
      m.fxd.pt <- m  
      m.fxd.pt[[1]] <- as.name("CDF.Pval.au.eq.u")
      m.fxd.pt$N.tests <- m.fxd.pt$lambda <- m.fxd.pt$corr.struct <- m.fxd.pt$FDP.control.method <- m.fxd.pt$delta <- NULL
      gma <- gma.st <- eval(m.fxd.pt, sys.parent())$gamma
      c.g <- qdist(1 - gma*alpha/2^(!is.pos), pars0)
      average.power <- pi.1 <- 1 - pdist(c.g, pars1) + (!is.pos) * pdist(-c.g, pars1)
 
      ## Prepare calls to variance functions.
      var.VoR.call <- as.call(expression(var.VoR))
      var.VoR.call$r.1 <- m$r.1
      var.VoR.call$alpha <- alpha
      var.VoR.call$delta <- eval(m$delta)
      var.VoR.call$gamma <- gma
      var.VoR.call$n.sample <- n.sample
      var.VoR.call$nii.sample <- nii.sample
      var.VoR.call$groups <- groups 
      var.VoR.call$effect.size <- effect.size
      var.VoR.call$type <- type
      var.VoR.call$distopt <- distopt
      var.VoR.call$FDP.control.method <- FDP.control.method

      if(is.cs) var.VoR.call$corr.struct <- corr.struct
      
      var.Rom.call <- var.ToM.call <- var.VoR.call
      var.Rom.call[[1]] <- as.name("var.Rom")
      var.ToM.call[[1]] <- as.name("var.ToM")

      ## Now call variance functions.
      sigma.rtm.VoR <- eval(var.VoR.call)^0.5
      tau2 <- eval(var.Rom.call)
      sigma.rtm.ToM <- eval(var.ToM.call)^0.5
    }
    if(FDP.cnt.mthd == "BHFDX")
    {
      m.cntlfdp <- m  
      m.cntlfdp[[1]] <- as.name("controlFDP")
      
      m.cntlfdp$lambda <- m.cntlfdp$FDP.control.method <- NULL
      alpha.st <- eval(m.cntlfdp, sys.parent())$alpha.star
      if(ast.le.a) alpha.st <- min(alpha.st, alpha)
      if(is.na(alpha.st)) FDP.control.method <- FDP.cnt.mthd <- Auto <- "Romano"
      if(!is.na(alpha.st))
      {
        m.fxd.pt <- m          
        m.fxd.pt[[1]] <- as.name("CDF.Pval.au.eq.u")
        m.fxd.pt$N.tests <- m.fxd.pt$lambda <- m.fxd.pt$corr.struct <- m.fxd.pt$FDP.control.method <- m.fxd.pt$delta <- NULL
        m.fxd.pt$alpha <- alpha.st
        gma <- gma.st <- eval(m.fxd.pt, sys.parent())$gamma
        c.g <- qdist(1 - gma.st*alpha.st/2^(!is.pos), pars0)
        average.power <- pi.1 <- 1 - pdist(c.g, pars1) + (!is.pos) * pdist(-c.g, pars1)

      ## Prepare calls to variance functions.
        var.VoR.call <- as.call(expression(var.VoR))
        var.VoR.call$r.1 <- m$r.1
        var.VoR.call$alpha <- alpha.st
        var.VoR.call$delta <- eval(m$delta)
        var.VoR.call$gamma <- gma.st
        var.VoR.call$n.sample <- n.sample
        var.VoR.call$nii.sample <- nii.sample
        var.VoR.call$groups <- groups 
        var.VoR.call$effect.size <- effect.size
        var.VoR.call$type <- type
        var.VoR.call$distopt <- distopt
        var.VoR.call$FDP.control.method <- FDP.control.method

        if(is.cs)
          var.VoR.call$corr.struct <- corr.struct
      
        var.Rom.call <- var.ToM.call <- var.VoR.call
        var.Rom.call[[1]] <- as.name("var.Rom")
        var.ToM.call[[1]] <- as.name("var.ToM")
 
        ## Now call variance functions.
        sigma.rtm.VoR <- eval(var.VoR.call)^0.5
        tau2 <- eval(var.Rom.call)
        sigma.rtm.ToM <- eval(var.ToM.call)^0.5
      }
    }
    if(FDP.cnt.mthd == "Romano")
    {
      m.fxd.pt <- m
      m.fxd.pt[[1]] <- as.name("CDF.Pval.apsi.eq.u")
      m.fxd.pt$N.tests <- m.fxd.pt$lambda <- m.fxd.pt$corr.struct <- m.fxd.pt$FDP.control.method <- NULL
      gma <- gma.R <- eval(m.fxd.pt, sys.parent())$gamma
    
      psi. <- gma.R*delta/(1-(1-delta)*gma.R)
      psi.pr <- delta/(1-(1-delta)*gma.R)^2
      u <- psi.*alpha
      c.g <- qdist(1 - u/2^(!is.pos), pars0)
      average.power <- pi.1 <- 1 - pdist(c.g, pars1) + (!is.pos)*pdist(-c.g, pars1)
      
      ## Prepare calls to variance functions.

      var.VoR.call <- as.call(expression(var.VoR))
      var.VoR.call$r.1 <- m$r.1
      var.VoR.call$alpha <- alpha
      var.VoR.call$delta <- eval(m$delta)
      var.VoR.call$gamma <- gma.R
      var.VoR.call$n.sample <- n.sample
      var.VoR.call$nii.sample <- nii.sample
      var.VoR.call$groups <- groups 
      var.VoR.call$effect.size <- effect.size
      var.VoR.call$type <- type
      var.VoR.call$distopt <- distopt
      var.VoR.call$FDP.control.method <- FDP.control.method

      if(is.cs)
        var.VoR.call$corr.struct <- corr.struct
      
      var.Rom.call <- var.ToM.call <- var.VoR.call
      var.Rom.call[[1]] <- as.name("var.Rom")
      var.ToM.call[[1]] <- as.name("var.ToM")

      ## Now call variance functions.
      sigma.rtm.VoR <- eval(var.VoR.call)^0.5
      tau2 <- eval(var.Rom.call)
      sigma.rtm.ToM <- eval(var.ToM.call)^0.5
    }

    if(control$verb>=3) browser()
    
    eIII <- dists$pdist[[1+idistopt]](-c.g, pars1)
    out <- list(average.power = average.power, c.g = c.g, gamma = gma, err.III=eIII)
    if(!is.N)
    {
      out$sigma.rtm.Rom <- max(tau2^0.5, 1e-10)
      out$sigma.rtm.VoR <- max(sigma.rtm.VoR, 1e-10)
      out$sigma.rtm.ToM <- max(sigma.rtm.ToM, 1e-10)
    }
    if(is.N)
    {
      out$se.Rom <- max((tau2/N.tests)^0.5, 1e-10)
      out$se.VoR <- max(sigma.rtm.VoR/N.tests^0.5, 1e-10)
      out$se.ToM <- max(sigma.rtm.ToM/N.tests^0.5, 1e-10)
    }
    
    if(!is.null(Auto)) out$Auto <- Auto
    if(!is.null(alpha.st)&&!is.na(alpha.st)) out$alpha.star <- alpha.st
    if(do.TPX.power)
    {
      se.ToM <- sigma.rtm.ToM/N.tests^0.5
      TPX.power <- pnorm((average.power - lambda)/se.ToM)
      L.eq <- average.power - se.ToM * qnorm(average.power)
      out$TPX.power <- TPX.power
      out$L.eq <- L.eq
    }
    attr(out, "passback") <- list(pars=list(pars0=pars0,pars1=pars1), FDP.control.method=FDP.control.method)
    
    out
}

"pwrFDR.BIN.1X" <-
function(effect.size, n.sample, r.1, alpha, delta=NULL, groups=2, N.tests,
         average.power, TPX.power, lambda, type, grpj.per.grp1, corr.struct,
         FDP.control.method, distopt=1, method, n.sim, temp.file, control)
{
    ## calculates TPX power
    ## specify lambda. Default, lambda=1, total power.

    ## default to Bonferroni (gamma=1/m) when gamma is missing
    ## otherwise you can try to approximate the distribution of a sequential method
    ## via plugin estimate for the fixed point, gamma=R/m 
    
    bad <- missing(alpha) || missing(effect.size) || missing(n.sample) || missing(N.tests)
    if(bad) stop("Arguments 'alpha', 'effect.size', 'n.sample', 'N.tests' are required")

    ## to do: allow arbitrary fixed threshold
    # if(missing(gamma)) gamma <- 1/N.tests 
    # if(missing(lambda)) lambda <- 1
    if(missing(groups)) groups <- 2

    do.tp.pwr <- !missing(lambda)
    do.avg.pwr <- missing(lambda)

    if(do.tp.pwr && missing(r.1)) stop("Argument 'r.1' is required for when calculating 'TPX.power'")

    ## wait, think about this
    ## Bonferroni is applied under the global null hypothesis but is still conservative under
    ## mixed null/alternative. Is this why var(T/M) doesn't aggree with simulation? 
    ## if(!do.tp.pwr && !missing(r.1)) stop("Don't specify argument 'r.1' unless calculating 'TPX.power'")
   
    m <- N.tests
    n <- n.sample
    d <- groups
    
    ## wait think about this:
    ## if(missing(r.1)) r.1 <- 0
    
    nii.sample <- n.sample
    if(groups>=2 && type==3) nii.sample <- n.sample*grpj.per.grp1

    .DF. <- eval(DF[type])
    .NCP. <- eval(NCP[type])
    
    idistopt <- distopt
    pars0 <- eval(dists$pars0[[1+idistopt]])
    pars1 <- eval(dists$pars1[[1+idistopt]])
    
    is.pos <- (dists$minv[[1+idistopt]]==0)
    c.g <- qdist(1 - alpha/(m*2^(!is.pos)), pars0)
    pwr.per.test <- 1-pdist(c.g, pars1) + (!is.pos) * pdist(-c.g, pars1)
    G.1 <- r.1 * pwr.per.test

    ## P( T/M >= lambda )
    ## = sum_{k=1}^m P(M=k) P(T/M >= lambda | M=k)
    ## = sum_{k=1}^m P(M=k) P(T >= k*lambda | M=k)
    ## = sum_{k=1}^m P(M=k) (1 - P(T < k*lambda | M=k))
    ## = sum_{k=1}^m P(M=k) (1 - P(T %in% c(0,1,..., max(ceiling(k*lambda)-1,0) | M=k))
    ## = sum_{k=1}^m P(M=k) (1 - pbinom( max( ceiling(k*lambda)-1, 0), size=k, prob=p))
    ## = sum_{k=1}^m dbinom(k, size=m, prob=r.1) * pbinom( max( ceiling(k*lambda)-1, 0), size=k, prob=p, lower.tail=FALSE)
    if(do.tp.pwr)
    {
      s <- 0
      for(k in 1:m)
        s <- s + dbinom(k, size=m, prob=r.1) * pbinom(max(ceiling(k*lambda)-1,0), size=k, prob=pwr.per.test, lower.tail=FALSE)
      out <- list(average.power=pwr.per.test, TPX.power=s)
    }
    if(do.avg.pwr) out <- list(average.power=pwr.per.test)

    ## For Bonferroni, two different gamma values:
    ##   Inside of G:  threshold alpha*gamma   gamma=1/m
    ##   and Outside of G:  R/m -> gamma, gamma = G 
    r.0 <- 1-r.1
    G.0 <- r.0*alpha/m  ## inside
    G <- G.0 + G.1
    pi.1 <- pwr.per.test

    out$gamma <- gamma <- G  ## outside

    ## R/m -> gamma
    ## V/m -> G.0
    ## T/m -> G.1
    ## M/m -> r.1
    ## 
    ## var(R/m) = G*(1-G)/m
    ## var(V/m) = G.0*(1-G.0)/m
    ## var(T/m) = G.1*(1-G.1)/m
    ##
    ## cov(V/m,R/m) = var(R/m) = G*(1-G)/m
    ## cov(T/m,M/m) = var(T/m) = G.1*(1-G.1)/m
    ##
    ## the other two, var(V/R) and var(T/M), are ratio estimates
    ## var(X/Y)
    ##  = var(1/y*(X-x) -x/y^2*(Y-y))
    ##  = 1/y^2*var(X) - 2*x/y^3*cov(X,Y) + x^2/y^4*var(Y)
    ##
    out$se.Rom <- (G*(1-G)/m)^0.5
    out$se.VoR <- (G.0*(1-G.0)/(m*gamma^2) + G.0^2*G*(1-G)/(m*gamma^4) - 2*G.0*G.0*(1-G.0)/(m*gamma^3))^0.5
    out$se.ToM <- (G.1*(1-G.1)/(m*r.1^2) + G.1^2*r.1*(1-r.1)/(m*r.1^4) - 2*G.1*G.1*(1-G.1)/(m*r.1^3))^0.5
    if(control$verb>=3) browser()
    out
}

"pwrFDR.sim" <-
function(effect.size, n.sample, r.1, alpha, delta, groups, N.tests, lambda, type,
         grpj.per.grp1, corr.struct, FDP.control.method, distopt, n.sim=1000, control)
{
    .call. <- m <- match.call()

    is.cs <- !missing(corr.struct)
   
    m.FP.1X <- m
    m.FP.1X$n.sim <- NULL
    m.FP.1X$do.Romano <- NULL
    m.FP.1X[[1]] <- as.name("pwrFDR.FP.1X")
    m.FP.1X$FDP.control.method <- "BHFDR"
    rslt.FP.1X <- eval(m.FP.1X, sys.parent())
    low.power.stop <- control$low.power.stop
    if(is.null(control$low.power.stop)) low.power.stop <- TRUE
    if(rslt.FP.1X$average.power < 0.50 && low.power.stop)
        stop("You don't want to run a simulation on a set of inputs with average power < 0.50")
    verb <- control$verb

    ast.le.a <- control$ast.le.a
    Auto <- NULL
    
    nii.sample <- n.sample
    if(groups>=2 && type==3) nii.sample <- n.sample*grpj.per.grp1
    
    ## if(groups==1) type is 1 (paired)
    ## if(groups==2) type is 1 (paired), 2 (balanced) or 3 (unbalanced)
    ## if(groups>=3) type is 2 (balanced) or 3 (unbalanced)
    is.pos <- (dists$minv[[1+distopt]]==0)

    nsim <- n.sim
    do.TPX.power <- !missing(lambda)
    alpha.st <- alpha
    BHFDX.lvl <- control$sim.level*(FDP.control.method=="BHFDX" || FDP.control.method=="both")
    Rmno <- (FDP.control.method %in% c("Romano","Holm","Hochberg") || FDP.control.method=="both")
    Bonf <- substring(FDP.control.method, 1, 4)=="Bonf"
    if(Bonf)BHFDX.lvl <- 0
    step.down <- Rmno && FDP.control.method != "Hochberg"
    
    if(control$verb>=4) browser()

    if(BHFDX.lvl>=1)
    {
        m.cntlFDP <- m
        m.cntlFDP$n.sim <- m.cntlFDP$lambda <- m.cntlFDP$FDP.control.method <- m.cntlFDP$do.Romano <- NULL
        m.cntlFDP[[1]] <- as.name("controlFDP")        
        alpha.st <- eval(m.cntlFDP, sys.parent())$alpha.star
        alpha.st <- ifelse(!is.na(alpha.st), alpha.st, -1)
        if(ast.le.a) alpha.st <- min(alpha.st, alpha)
        BHFDX.lvl <- ifelse(alpha.st > 0, BHFDX.lvl, 0)
        Rmno <- Rmno || (alpha.st < 0)
        if(Rmno) Auto <- "Romano"
        step.down <- Rmno && FDP.control.method != "Hochberg"
    }
    if(!is.cs)
    {
      if(verb>0)
      {
        cat(sprintf("nsim=%d, alpha=%g, N.tests=%d, r.1=%g, n.sample=%d, effect.size=%g\n",
                    nsim, alpha, N.tests, r.1, n.sample, effect.size))
      }
      if(!Bonf)
        rslt <- pwrFDRsim.fn(n.sim,alpha,delta,BHFDX.lvl,Rmno,alpha.st,N.tests,r.1,n.sample,nii.sample,effect.size,
                             groups,type,distopt,step.down,control)
      if(Bonf)
        rslt <- pwrFDRsimBNF.fn(n.sim,alpha,N.tests,r.1,n.sample,nii.sample,effect.size,
                                groups,type,distopt,control)     
    }
    if(is.cs)
    {
      k.bs <- corr.struct$block.size
      cty <- corr.struct$type
      n.cty <- list(`CS-Blocks`=1, `Toeplitz-Blocks`=2)[[cty]]
      if(cty=="CS-Blocks") RHO.5 <- corr.struct$rho^0.5
      if(cty=="Toeplitz-Blocks") RHO.5 <- chol(toeplitz(c(1,corr.struct$rho)))

      if(verb>0)
      {
        cat(sprintf("nsim=%d, alpha=%g, N.tests=%d, r.1=%g, n.sample=%d, effect.size=%g, blk.sz=%d, type=%s, rho=",
                     nsim, alpha, N.tests, r.1, n.sample, effect.size, k.bs, type))
        cat(sprintf("%g\n", corr.struct$rho))
      }
      rslt <- pwrFDRsimCS.fn(n.sim,alpha,delta,BHFDX.lvl,Rmno,alpha.st,N.tests,r.1,n.sample,nii.sample,effect.size,
                             groups,type,distopt,RHO.5,k.bs,n.cty,step.down,control)
    }
    var.0 <- function(x)ifelse(length(x)>1, var(x), 0)
    reps <- data.frame(M1=rslt$M1, R=rslt$R, T=rslt$T)
    average.power <- with(reps, mean(T %over% M1))
    emp.FDR <- with(reps, mean((R-T) %over% R))
    gamma <- with(reps, mean(R)/N.tests)
    se.Rom <- with(reps, var.0(R/N.tests)^0.5)
    se.VoR <- with(reps, var.0((R - T) %over% R)^0.5)
    se.ToM <- with(reps, var.0(T %over% M1)^0.5)
    out <- list(average.power=average.power, emp.FDR=emp.FDR, gamma=gamma,
                se.Rom=se.Rom, se.VoR=se.VoR, se.ToM=se.ToM)
    if(do.TPX.power)
    {
      TPX.power <- with(reps, mean(T %over% M1 >= lambda))
      out$TPX.power <- TPX.power
      ## reorder elements of out so that TPX.power and lambda are just
      ## ahead of average.power
      nms.out <- names(out)
      out.apnd <- list(TPX.power=out$TPX.power, lambda=.call.$lambda)
      ind.avgpwr <- which(nms.out=="average.power")
      out <- out[-which(nms.out=="TPX.power")]
      out <- append(out, out.apnd, after=ind.avgpwr-1)
    }
    
    se.VoR.st <- se.ToM.st <- se.VoR.st.ht <- se.ToM.st.ht <- se.VoR.R <-
        se.ToM.R <- NULL
    if(BHFDX.lvl >= 1)
    {
      reps <- cbind(reps, R.st=rslt$R.st, T.st=rslt$T.st)
      average.power.st <- with(reps, mean(T.st %over% M1))
      emp.FDR.st <- with(reps, mean((R.st-T.st) %over% R.st))
      P.st <- with(reps, mean((R.st - T.st) %over% R.st > alpha))
      se.Rom.st <-  with(reps, var.0(R.st/N.tests)^0.5)
      se.VoR.st <- with(reps, var.0((R.st - T.st) %over% R.st)^0.5)
      se.ToM.st <- with(reps, var.0(T.st %over% M1)^0.5)
      gamma.st <- with(reps, mean(R.st)/N.tests)
      out <- c(out, average.power.st=average.power.st, emp.FDR.st=emp.FDR.st,
               gamma.st=gamma.st, se.Rom.st=se.Rom.st, se.VoR.st=se.VoR.st,
               se.ToM.st=se.ToM.st, alpha.star=alpha.st)
      if(do.TPX.power)
      {
        TPX.power.st <- with(reps, mean(T.st %over% M1 >= lambda))
        out$TPX.power.st <- TPX.power.st
        ## reorder elements of out so that TPX.power and lambda are just
        ## ahead of average.power
        nms.out <- names(out)
        out.apnd <- list(TPX.power=out$TPX.power, lambda=.call.$lambda)
        ind.avgpwr <- which(nms.out=="average.power")
        out <- out[-which(nms.out=="TPX.power")]
        out <- append(out, out.apnd, after=ind.avgpwr-1)
      }
    }
    ##if(BHFDX.lvl == 2)
    ##{
    ##  reps <- cbind(reps, p0.ht=pmin(rslt$p0.ht,1))
    ##  out <- c(out, p0.ht.avg=mean(rslt$p0.ht), sd.p0.ht=var(rslt$p0.ht)^0.5)
    ##}
    if(Rmno)
    {
      reps <- cbind(reps, R.R=rslt$R.R, T.R=rslt$T.R)
      average.power.R <- with(reps, mean(T.R %over% M1))
      emp.FDR.R <- with(reps, mean((R.R-T.R) %over% R.R))
      P.R <- with(reps, mean((R.R - T.R) %over% R.R > alpha))
      se.Rom.R <-  with(reps, var.0(R.R/N.tests)^0.5)
      se.VoR.R <- with(reps, var.0((R.R - T.R) %over% R.R)^0.5)
      se.ToM.R <- with(reps, var.0(T.R %over% M1)^0.5)
      gamma.R <- with(reps, mean(R.R)/N.tests)
      out <- c(out, average.power.R=average.power.R, emp.FDR.R=emp.FDR.R,
               gamma.R=gamma.R,se.Rom.R=se.Rom.R, se.VoR.R=se.VoR.R, se.ToM.R=se.ToM.R)
      if(do.TPX.power)
      {
        TPX.power.R <- with(reps, mean(T.R %over% M1 >= lambda))
        out$TPX.power.R <- TPX.power.R
      }
    }
    
    dtl <- list(reps=reps, dat1x=rslt$dat1x)
    attr(out, "detail") <- dtl

    .DF. <- eval(DF[type])
    .NCP. <- eval(NCP[type])

    pars0 <- eval(dists$pars0[[1+distopt]])
    pars1 <- eval(dists$pars1[[1+distopt]])
    
    if(!is.null(Auto)) out$Auto <- Auto
    attr(out, "passback") <- list(pars=list(pars0=pars0,pars1=pars1))
    out$call <- .call.
    class(out) <- "pwrFDR"
    out
}

"pwrFDRsim.fn" <-
function(n.sim,alpha,delta,BHFDX.lvl,Rmno,alpha.st,N.tests,r.1,n.sample,nii.sample,effect.size,groups,
         type,distopt,step.down,control)
{
  BHFDX <-BHFDX.lvl
  do.Rmno <- Rmno
  alpha.star <- alpha.st
  Ng <- N.tests
  r1 <- r.1
  ##  n <- xn <- n.sample
  ##  theta <- effect.size

  .DF. <- eval(DF[type])
  .NCP. <- eval(NCP[type])

  pars0 <- eval(dists$pars0[[1+distopt]])
  pars1 <- eval(dists$pars1[[1+distopt]])
  rdist <- eval(dists$rdist[[1+distopt]])
  pdist <- eval(dists$pdist[[1+distopt]])
  
  is.pos <- (dists$minv[[1+distopt]]==0)
  verb <- control$verb

  if(verb>0) browser()
  rslt <-
    t(replicate(n.sim, {
      M1 <- rbinom(1, size=N.tests, prob=r.1)
      X1 <- rdist(M1, pars1)
      X0 <- rdist(N.tests-M1, pars0)
      xi <- c(rep(1, M1), rep(0, N.tests-M1))
      X <- c(X1, X0)
      P <- 2*(1-pdist(abs(X), pars0))
      dat <- cbind(xi=xi, X=X, P=P)
      dimnames(dat)[[2]] <- c("xi", "X", "P")
      dat.o <- dat[order(dat[,"P"]),]
      R. <- suppressWarnings(if.absInf.x(max(which(dat.o[,"P"] <= alpha*(1:N.tests)/N.tests)), 0))
      T. <- suppressWarnings(if.na.x(sum(xi*(P<=alpha*R./N.tests)),0))
      R.st <- T.st <- R.R <- T.R <- NULL
      if(BHFDX.lvl>=1)
      {
        R.st <- suppressWarnings(if.absInf.x(max(which(dat.o[,"P"] <= alpha.st*(1:N.tests)/N.tests)), 0))
        T.st <- suppressWarnings(if.na.x(sum(xi*(P<=alpha.st*R.st/N.tests)),0))
      }
      if(do.Rmno && step.down)
      {
        R.R <- suppressWarnings(if.absInf.x(min(which(dat.o[,"P"] > alpha*psi.m(1:N.tests, delta, N.tests))) - 1,0))
        T.R <- suppressWarnings(if.absInf.x(sum(dat.o[min(1, R.R):R.R,"xi"]), 0))
      }
      if(do.Rmno && !step.down)
      {
        R.R <- suppressWarnings(if.absInf.x(max(which(dat.o[,"P"] <= alpha*psi.m(1:N.tests, delta, N.tests))), 0))
        T.R <- suppressWarnings(if.na.x(sum(xi*(P <= alpha*psi.m(R.R, delta, N.tests))),0))
      }
      c(M1=M1, R=R., T=T., R.st=R.st, T.st=T.st, R.R=R.R, T.R=T.R)
  }))
  ## do 1X and save sorted p-values for output
  
  M1 <- rbinom(1, size=N.tests, prob=r.1)
  X1 <- rdist(M1, pars1)
  X0 <- rdist(N.tests-M1, pars0)
  xi <- c(rep(1, M1), rep(0, N.tests-M1))
  X <- c(X1, X0)
  P <- 2*(1-pdist(abs(X), pars0))
  dat <- cbind(xi=xi, X=X, P=P)
  dimnames(dat)[[2]] <- c("xi", "X", "P")
  dat.o <- dat[order(dat[,"P"]),]
  
  out <- list()
  out$n.sim <- n.sim
  out$alpha <- alpha
  out$BHFDX <- BHFDX.lvl
  out$do.Rmno <- Rmno
  out$alpha.star <- alpha.st
  out$Ng <- N.tests
  out$r.1 <- r.1
  out$n <- n.sample
  out$nii.sample <- nii.sample
  out$theta <- effect.size
  out$distopt <- distopt
  out$groups <- groups
  out$delta <- delta
  out$pverb <- verb
  out$M1 <- rslt[,1]
  out$R <- rslt[,2]
  out$T <- rslt[,3]
  cols.bhclt <- 0
  if(BHFDX.lvl>=1)
  {
    out$R.st <- rslt[,"R.st"]
    out$T.st <- rslt[,"T.st"]
    cols.bhclt <- 2
  }
  if(Rmno)
  {
    out$R.R <- rslt[,"R.R"]
    out$T.R <- rslt[,"T.R"]
  }
  out$dat1x <- as.data.frame(dat.o)
  out
}

"pwrFDRsimCS.fn" <-
function(n.sim,alpha,delta,BHFDX.lvl,Rmno,alpha.st,N.tests,r.1,n.sample,nii.sample,effect.size,groups,
         type,distopt,RHO.5,k.bs,n.cty,step.down,control)
{
  BHFDX <-BHFDX.lvl
  do.Rmno <- Rmno
  alpha.star <- alpha.st
  Ng <- m <- N.tests
  r1 <- r.1
  n.ty <- n.cty
  pverb <- control$verb
  
  .DF. <- eval(DF[type])
  .NCP. <- eval(NCP[type])

  pars0 <- eval(dists$pars0[[1+distopt]])
  pars1 <- eval(dists$pars1[[1+distopt]])
  rdist <- eval(dists$rdist[[1+distopt]])
  
  is.pos <- (dists$minv[[1+distopt]]==0)
  verb <- control$verb

  ncp <- pars1[2]
  ## n.cty <- list(`CS-Blocks`=1, `Toeplitz-Blocks`=2)[[cty]]
  frc <- N.tests/k.bs - floor(N.tests/k.bs)
  Ntst.rnd <- k.bs*floor(N.tests/k.bs)
  if(verb>0) browser()
  if(n.cty==1)
  {
    rslt <-
        t(replicate(n.sim, {
          if(verb==-768) browser()
          tau <- (1 - RHO.5^2)^0.5
          err <- rnorm(N.tests, sd=tau)
          Z <- rep(rnorm(Ntst.rnd/k.bs, sd=RHO.5), each=k.bs)
          if(frc>0) Z <- c(Z, rep(rnorm(1, sd=RHO.5), round(frc*k.bs)))
          xi <- 1*(runif(N.tests) <= r.1)
          M1 <- sum(xi)
          MU <- xi*ncp
          X <- MU + Z + err
          P <- 2*(1-pnorm(abs(X)))
          dat <- cbind(xi=xi, Z, X, P=P)
          dimnames(dat)[[2]] <- c("xi", "Z", "X", "P")
          dat.o <- dat[order(dat[,"P"]),]
          R. <- suppressWarnings(if.absInf.x(max(which(dat.o[,"P"] <= alpha*(1:N.tests)/N.tests)),0))
          T. <- suppressWarnings(if.na.x(sum(xi*(P<=alpha*R./N.tests)),0))
          R.st <- T.st <- R.R <- T.R <- NULL
          if(BHFDX.lvl>=1)
          {
            R.st <- suppressWarnings(if.absInf.x(max(which(dat.o[,"P"] <= alpha.st*(1:N.tests)/N.tests)),0))
            T.st <- suppressWarnings(if.na.x(sum(xi*(P<=alpha.st*R.st/N.tests)),0))
          }
          if(do.Rmno && step.down)
          {
            R.R <- suppressWarnings(if.absInf.x(min(which(dat.o[,"P"] > alpha*psi.m(1:N.tests, delta, N.tests))) - 1,0))
            T.R <- suppressWarnings(if.na.x(sum(xi*(P<=alpha*psi.m(R.R, delta, N.tests))),0))
          }
          if(do.Rmno && !step.down)
          {
            R.R <- suppressWarnings(if.absInf.x(max(which(dat.o[,"P"] <= alpha*psi.m(1:N.tests, delta, N.tests))), 0))
            T.R <- suppressWarnings(if.na.x(sum(xi*(P <= alpha*psi.m(R.R, delta, N.tests))),0))
          }
          c(M1=M1, R=R., T=T., R.st=R.st, T.st=T.st, R.R=R.R, T.R=T.R)
        }))
      tau <- (1 - RHO.5^2)^0.5
      err <- rnorm(N.tests, sd=tau)
      Z <- rep(rnorm(Ntst.rnd/k.bs, sd=RHO.5), each=k.bs)
      if(frc>0) Z <- c(Z, rep(rnorm(1, sd=RHO.5), round(frc*k.bs)))
      xi <- 1*(runif(N.tests) <= r.1)
      M1 <- sum(xi)
      MU <- xi*ncp
      X <- MU + Z + err
      P <- 2*(1-pnorm(abs(X)))
      dat <- cbind(xi=xi, Z, X, P=P)
      dimnames(dat)[[2]] <- c("xi", "Z", "X", "P")
      dat.o <- dat[order(dat[,"P"]),]
   }
  
  if(n.cty==2)
  {
    rslt <-
      t(replicate(n.sim, {
          Z0 <- matrix(rnorm(Ntst.rnd), Ntst.rnd/k.bs, k.bs)
          Z <- c(Z0%*%RHO.5)
          if(frc>0) Z <- c(Z, rnorm(round(frc*k.bs))%*%RHO.5[1:round(frc*k.bs), 1:round(frc*k.bs)])
          xi <- 1*(runif(N.tests) <= r.1)
          M1 <- sum(xi)
          MU <- xi*ncp
          X <- MU + Z
          P <- 2*(1-pnorm(abs(X)))
          dat <- cbind(xi=xi, Z, X, P=P)
          dimnames(dat)[[2]] <- c("xi", "Z", "X", "P")
          dat.o <- dat[order(dat[,"P"]),]
          R. <- suppressWarnings(if.absInf.x(max(which(dat.o[,"P"] <= alpha*(1:N.tests)/N.tests)),0))
          T. <- suppressWarnings(if.na.x(sum(xi*(P<=alpha*R./N.tests)),0))
          R.st <- T.st <- R.R <- T.R <- NULL
          if(BHFDX.lvl>=1)
          {
            R.st <- suppressWarnings(if.absInf.x(max(which(dat.o[,"P"] <= alpha.st*(1:N.tests)/N.tests)),0))
            T.st <- suppressWarnings(if.na.x(sum(xi*(P<=alpha.st*R.st/N.tests)),0))
          }
          if(do.Rmno && step.down)
          {
            R.R <- suppressWarnings(if.absInf.x(min(which(dat.o[,"P"] > alpha*psi.m(1:N.tests, delta, N.tests))) - 1,0))
            T.R <- suppressWarnings(if.absInf.x(sum(dat.o[min(1, R.R):R.R,"xi"]), 0))
          }
          if(do.Rmno && !step.down)
          {
            R.R <- suppressWarnings(if.absInf.x(max(which(dat.o[,"P"] <= alpha*psi.m(1:N.tests, delta, N.tests))), 0))
            T.R <- suppressWarnings(if.na.x(sum(xi*(P <= alpha*psi.m(R.R, delta, N.tests))),0))              
          }
          c(M1=M1, R=R., T=T., R.st=R.st, T.st=T.st, R.R=R.R, T.R=T.R)
      }))
      Z0 <- matrix(rnorm(Ntst.rnd), Ntst.rnd/k.bs, k.bs)
      Z <- c(Z0%*%RHO.5)
      if(frc>0) Z <- c(Z, rnorm(frc*k.bs)%*%RHO.5[1:round(frc*k.bs), 1:round(frc*k.bs)])
      xi <- 1*(runif(N.tests) <= r.1)
      M1 <- sum(xi)
      MU <- xi*ncp
      X <- MU + Z
      P <- 2*(1-pnorm(abs(X)))
      dat <- cbind(xi=xi, Z, X, P=P)
      dimnames(dat)[[2]] <- c("xi", "Z", "X", "P")
      dat.o <- dat[order(dat[,"P"]),]
   }

  out <- list()
  out$n.sim <- n.sim
  out$alpha <- alpha
  out$BHFDX <- BHFDX.lvl
  out$do.Rmno <- Rmno
  out$step.down <- step.down
  out$alpha.star <- alpha.st
  out$Ng <- N.tests
  out$r.1 <- r.1
  out$n <- n.sample
  out$nii.sample <- nii.sample
  out$theta <- effect.size
  out$rho.5 <- RHO.5
  out$k.bs <- k.bs
  out$n.ty <- n.cty
  out$delta <- delta
  out$pverb <- verb
  out$M1 <- rslt[,1]
  out$R <- rslt[,2]
  out$T <- rslt[,3]
  cols.bhclt <- 0
  if(BHFDX.lvl>=1)
  {
    out$R.st <- rslt[,"R.st"]
    out$T.st <- rslt[,"T.st"]
    cols.bhclt <- 2
  }
  if(Rmno)
  {
    out$R.R <- rslt[,"R.R"]
    out$T.R <- rslt[,"T.R"]
  }
  out$dat1x <- as.data.frame(dat.o)
  out
}

pwrFDRsimBNF.fn <- function(n.sim,alpha,N.tests,r.1,n.sample,nii.sample,effect.size, groups,type,distopt,control)     
{
  .DF. <- eval(DF[type])
  .NCP. <- eval(NCP[type])

  pars0 <- eval(dists$pars0[[1+distopt]])
  pars1 <- eval(dists$pars1[[1+distopt]])
  rdist <- eval(dists$rdist[[1+distopt]])
  pdist <- eval(dists$pdist[[1+distopt]])
  
  is.pos <- (dists$minv[[1+distopt]]==0)
  verb <- control$verb

  if(verb>0) browser()
  rslt <-
    t(replicate(n.sim, {
      M1 <- rbinom(1, size=N.tests, prob=r.1)
      X1 <- rdist(M1, pars1)
      X0 <- rdist(N.tests-M1, pars0)
      xi <- c(rep(1, M1), rep(0, N.tests-M1))
      X <- c(X1, X0)
      P <- 2*(1-pdist(abs(X), pars0))
      dat <- cbind(xi=xi, X=X, P=P)
      dimnames(dat)[[2]] <- c("xi", "X", "P")
      dat.o <- dat[order(dat[,"P"]),]
      R. <- suppressWarnings(if.absInf.x(max(which(dat.o[,"P"] <= alpha/N.tests)), 0))
      T. <- suppressWarnings(if.na.x(sum(xi*(P<=alpha/N.tests)),0))
      c(M1=M1, R=R., T=T.)
  }))
  ## do 1X and save sorted p-values for output
  
  M1 <- rbinom(1, size=N.tests, prob=r.1)
  X1 <- rdist(M1, pars1)
  X0 <- rdist(N.tests-M1, pars0)
  xi <- c(rep(1, M1), rep(0, N.tests-M1))
  X <- c(X1, X0)
  P <- 2*(1-pdist(abs(X), pars0))
  dat <- cbind(xi=xi, X=X, P=P)
  dimnames(dat)[[2]] <- c("xi", "X", "P")
  dat.o <- dat[order(dat[,"P"]),]
  
  out <- list()
  out$n.sim <- n.sim
  out$alpha <- alpha
  out$Ng <- N.tests
  out$r.1 <- r.1
  out$n <- n.sample
  out$nii.sample <- nii.sample
  out$theta <- effect.size
  out$distopt <- distopt
  out$groups <- groups
  out$pverb <- verb
  out$M1 <- rslt[,1]
  out$R <- rslt[,2]
  out$T <- rslt[,3]
  out$dat1x <- as.data.frame(dat.o)
  out
}

seFDP.o.alpha <-
function(effect.size, n.sample, r.1, alpha, groups, N.tests, type, grpj.per.grp1,
         distopt, rho, k.bs)
{
    m <- match.call()
    m[[1]] <- as.name("pwrFDR")
    m$rho <- m$k.bs <- NULL
    is.distopt <- !missing(distopt)
    is.rho <- !missing(rho)
    if(is.distopt) m$distopt <- distopt
    if(is.rho) m$corr.struct <- list(type="CS-Blocks", block.size=k.bs, rho=rho)
    rslt <- eval(m, sys.parent())
    rslt$se.VoR/alpha
}        

backsolve.seFDPoalpha <-
function(seFDPoalpha, effect.size, n.sample, r.1, alpha, groups=2, N.tests, type="balanced",
         grpj.per.grp1=1, distopt=1, rho, k.bs)
{
    if(missing(groups)) groups <- 2
    if(missing(type)) type <- "balanced"
    if(missing(grpj.per.grp1)) grpj.per.grp1 <- 1
    OBJ <-
    function(x, seFDPoalpha, effect.size, n.sample, r.1, alpha, groups, N.tests, type,
             grpj.per.grp1, distopt, rho, k.bs, msng)
    {
        asgn <- as.call(expression("f"))
        asgn[[1]] <- as.name("<-")
        asgn[[2]] <- msng
        if(!(msng%in%c("r.1","N.tests"))) asgn[[3]] <- exp(x)
        if(msng=="N.tests") asgn[[3]] <- ceiling(exp(x))
        if(msng=="r.1") asgn[[3]] <- logitInv(x)
        eval(asgn)
        seFDP.o.alpha.call <- as.call(expression(seFDP.o.alpha))
        seFDP.o.alpha.call$effect.size <- effect.size
        seFDP.o.alpha.call$n.sample <- n.sample
        seFDP.o.alpha.call$r.1 <- r.1
        seFDP.o.alpha.call$alpha <- alpha
        seFDP.o.alpha.call$N.tests <- N.tests
        is.distopt <- !missing(distopt)
        if(is.distopt) seFDP.o.alpha.call$distopt <- distopt
        is.rho <- !missing(rho)
        if(is.rho)
        {
          seFDP.o.alpha.call$rho <- rho
          seFDP.o.alpha.call$k.bs <- k.bs
        }
        seFDPoalpha - suppressWarnings(eval(seFDP.o.alpha.call))
    }
    m <- match.call()
    all.args <- names(formals(seFDP.o.alpha))
    all.args <- all.args[-which(all.args %in% c("rho","k.bs"))]
    sup.args <- names(m)
    sup.args <- sup.args[-which(sup.args %in% c("seFDPoalpha", "rho", "k.bs"))]
    extr.args <- c("groups", "type", "grpj.per.grp1", "distopt")
    msng <- setdiff(all.args, c(sup.args, extr.args))
    bad <- FALSE
    if(length(msng)!=1)bad <- TRUE
    if(!bad && msng=="alpha") bad <- TRUE
    if(bad)stop("You must supply all but one argument, to include 'seFDPoalpha' and 'alpha'")
    intvl <- list(N.tests=log(c(10,1e9)), effect.size=log(c(1e-5,2)), n.sample=log(c(5,1e9)), r.1=logit(c(1e-6,1-1e-6)))
    soln <- uniroot(OBJ, interval=intvl[[msng]],seFDPoalpha,effect.size,n.sample,r.1,alpha,groups,N.tests,type,
                    grpj.per.grp1,distopt,rho,k.bs,msng)
    if(!(msng%in%c("r.1","N.tests"))) out <- exp(soln$root)
    if(msng=="N.tests") out <- ceiling(exp(soln$root))
    if(msng=="r.1") out <- logitInv(soln$root)
    out.un <- out
    names(out) <- msng
    m[[msng]] <- out.un
    m[[1]] <- as.name("pwrFDR")
    is.distopt <- !missing(distopt)
    if(is.distopt) m$distopt <- distopt
    is.rho <- !missing(rho)
    if(is.rho) m$corr.struct <- list(type="CS-Blocks", block.size=k.bs, rho=rho)
    m$seFDPoalpha <- m$rho <- m$k.bs <- NULL
    pwr <- eval(m)
    out <- c(out, average.power=pwr$average.power, se.VoR=pwr$se.VoR, value=soln$f.root)
    out
}

seTPP.o.avgpwr <-
function(effect.size, n.sample, r.1, alpha, groups, N.tests, type, grpj.per.grp1,
         distopt, rho, k.bs)
{
    m <- match.call()
    m[[1]] <- as.name("pwrFDR")
    m$rho <- m$k.bs <- NULL
    is.distopt <- !missing(distopt)
    is.rho <- !missing(rho)
    if(is.distopt) m$distopt <- distopt
    if(is.rho) m$corr.struct <- list(type="CS-Blocks", block.size=k.bs, rho=rho)
    rslt <- eval(m, sys.parent())
    rslt$se.ToM/rslt$average.power
}        

backsolve.seTPPoavgpwr <-
function(seTPPoavgpwr, effect.size, n.sample, r.1, alpha, groups=2, N.tests, type="balanced",
         grpj.per.grp1=1, distopt=1, rho, k.bs)
{
    if(missing(groups)) groups <- 2
    if(missing(type)) type <- "balanced"
    if(missing(grpj.per.grp1)) grpj.per.grp1 <- 1
    OBJ <-
    function(x, seTPPoavgpwr, effect.size, n.sample, r.1, alpha, groups, N.tests, type,
             grpj.per.grp1, distopt, rho, k.bs, msng)
    {
        asgn <- as.call(expression("f"))
        asgn[[1]] <- as.name("<-")
        asgn[[2]] <- msng
        if(!(msng%in%c("r.1","N.tests"))) asgn[[3]] <- exp(x)
        if(msng=="N.tests") asgn[[3]] <- ceiling(exp(x))
        if(msng=="r.1") asgn[[3]] <- logitInv(x)
        eval(asgn)

        seTPP.o.avgpwr.call <- as.call(expression(seTPP.o.avgpwr))
        seTPP.o.avgpwr.call$effect.size <- effect.size
        seTPP.o.avgpwr.call$n.sample <- n.sample
        seTPP.o.avgpwr.call$r.1 <- r.1
        seTPP.o.avgpwr.call$alpha <- alpha
        seTPP.o.avgpwr.call$N.tests <- N.tests
        is.distopt <- !missing(distopt)
        if(is.distopt) seTPP.o.avgpwr.call$distopt <- distopt
        is.rho <- !missing(rho)
        if(is.rho)
        {
          seTPP.o.avgpwr.call$rho <- rho
          seTPP.o.avgpwr.call$k.bs <- k.bs
        }
        seTPPoavgpwr - suppressWarnings(eval(seTPP.o.avgpwr.call))
    }
    m <- match.call()
    all.args <- names(formals(seTPP.o.avgpwr))
    all.args <- all.args[-which(all.args %in% c("rho","k.bs"))]
    sup.args <- names(m)
    sup.args <- sup.args[-which(sup.args %in% c("seTPPoavgpwr", "rho", "k.bs"))]
    extr.args <- c("groups", "type", "grpj.per.grp1", "distopt")
    msng <- setdiff(all.args, c(sup.args, extr.args))
    bad <- FALSE
    if(length(msng)!=1)bad <- TRUE
    if(!bad && msng=="alpha") bad <- TRUE
    if(bad)stop("You must supply all but one argument, to include 'seTPPoavgpwr' and 'alpha'")
    intvl <- list(N.tests=log(c(1,1e9)), effect.size=log(c(1e-5,2)), n.sample=log(c(5,1e9)), r.1=logit(c(1e-6,1-1e-6)))
    soln <- try(uniroot(OBJ,interval=intvl[[msng]],seTPPoavgpwr,effect.size,n.sample,r.1,alpha,groups,N.tests,type,
                        grpj.per.grp1,distopt,rho,k.bs,msng), silent=TRUE)
    is.err <- class(soln)=="try-error"
    if(!is.err)
    {
      if(!(msng%in%c("r.1","N.tests"))) out <- exp(soln$root)
      if(msng=="N.tests") out <- ceiling(exp(soln$root))
      if(msng=="r.1") out <- logitInv(soln$root)
      out.un <- out
      names(out) <- msng
      m[[msng]] <- out.un
      m[[1]] <- as.name("pwrFDR")
      is.distopt <- !missing(distopt)
      if(is.distopt) m$distopt <- distopt
      is.rho <- !missing(rho)
      if(is.rho) m$corr.struct <- list(type="CS-Blocks", block.size=k.bs, rho=rho)
      m$seTPPoavgpwr <- m$rho <- m$k.bs <- NULL
      pwr <- eval(m)
      out <- c(out, average.power=pwr$average.power, se.ToM=pwr$se.ToM, value=soln$f.root)
    }else out <-paste(str_split(soln, '\n')[[1]], collapse="")
    out
}

"criterion" <-
function(alpha, delta, N.tests, FDP.control.method=c("BHFDR","Romano"))
{
    .call. <- m <- match.call()

    frmls <- names(formals(criterion))
    sppld <- names(m)[-1]

    is.msng <- !(frmls %in% sppld)
    names(is.msng) <- frmls

    err.msng <- any(c("alpha", "N.tests", "FDP.control.method") %in% names(is.msng[is.msng]))

    if(err.msng) stop("Arguments 'alpha', 'N.tests' and 'FDP.control.method' are required")
    
    if(is.msng["delta"]) delta <- alpha
    
    errs <- list(alpha=(alpha <= 0 || alpha >= 1),
                 delta=(delta <= 0 || delta >= 1),
                 N.tests=(!is.int(N.tests) || N.tests <= 0),
                 FDP.control.method=!(FDP.control.method %in% c("BHFDR", "Romano")))

    msgs <- list(alpha="Argument, 'alpha' must be between 0 and 1\n",
                 delta="Argument, 'delta' must be between 0 and 1\n",
                 N.tests="Argument, 'N.tests' must be integral and positive\n",
                 FDP.control.method="Argument, 'FDP.control.method' must be either 'BHFDR' or 'Romano'\n")
    
    err <- any(unlist(errs))
    msg <- paste(msgs[names(which(unlist(errs)))], collapse="")
    if(err) stop(msg)

    ii <- 1:m$N.tests
    
    if(FDP.control.method=="BHFDR") crit <- m$alpha*ii/m$N.tests
    if(FDP.control.method=="Romano") crit <- (floor(delta*ii) + 1)*alpha/(N.tests + floor(delta*ii) + 1 - ii)

    crit
}

"es.ROC" <-
function(FPR0, FPR1=NULL, TPR0, TPR1=NULL, b=NULL)
{
    m <- match.call()
    fn.nm <- as.character(m[[1]])
    is.FPR1 <- !missing(FPR1)
    is.TPR1 <- !missing(TPR1)
    if(!is.FPR1&!is.TPR1) stop("You must specify exactly one of the arguments 'FPR1' and 'TPR1'")
    if(is.TPR1) fn.nm <- fn.nm %,% ".fxFPR"
    if(is.FPR1) fn.nm <- fn.nm %,% ".fxTPR"
    m[[1]] <- as.name(fn.nm)
    eval(m)
}

"cc.ROC" <-
function(FPR0, FPR1=NULL, TPR0, TPR1=NULL, b=NULL)
{
    m <- match.call()
    fn.nm <- as.character(m[[1]])
    is.FPR1 <- !missing(FPR1)
    is.TPR1 <- !missing(TPR1)
    if(!is.FPR1&!is.TPR1) stop("You must specify exactly one of the arguments 'FPR1' and 'TPR1'")
    if(is.TPR1) fn.nm <- fn.nm %,% ".fxFPR"
    if(is.FPR1) fn.nm <- fn.nm %,% ".fxTPR"
    m[[1]] <- as.name(fn.nm)
    eval(m)
}

"es.ROC.fxTPR" <-
function(FPR0, FPR1, TPR0, b=NULL)
{
  if(missing(b)) b <- 1
  v.tpr0 <- TPR0*(1-TPR0)
  v.fpr1 <- FPR1*(1-FPR1)
  r <- b*dnorm(qnorm(1-TPR0))/dnorm(qnorm(1-FPR1))
  k <- (v.fpr1/v.tpr0)^0.5/r
#  k <- ceiling(k)
  nu <- (FPR1-FPR0)/(v.fpr1 + r^2/k*v.tpr0)^0.5
  abs(nu)
}

"cc.ROC.fxTPR" <-
function(FPR0, FPR1, TPR0, b=NULL)
{
  if(missing(b)) b <- 1
  v.tpr0 <- TPR0*(1-TPR0)
  v.fpr1 <- FPR1*(1-FPR1)
  r <- b*dnorm(qnorm(1-TPR0))/dnorm(qnorm(1-FPR1))
  k <- (v.fpr1/v.tpr0)^0.5/r
#  k <- ceiling(k)
  k
}

"es.ROC.fxFPR" <-
function(FPR0, TPR0, TPR1, b=NULL)
{
  if(missing(b)) b <- 1
  v.fpr0 <- FPR0*(1-FPR0)
  v.tpr1 <- TPR1*(1-TPR1)
  r <- b*dnorm(qnorm(TPR1))/dnorm(qnorm(FPR0))
  k <- (v.tpr1/v.fpr0)^0.5/r
#  k <- ceiling(k)
  nu <- (TPR1-TPR0)/(v.tpr1 + k*r^2*v.fpr0)^0.5
  abs(nu)
}

"cc.ROC.fxFPR" <-
function(FPR0, TPR0, TPR1, b=NULL)
{
  if(missing(b)) b <- 1
  v.fpr0 <- FPR0*(1-FPR0)
  v.tpr1 <- TPR1*(1-TPR1)
  r <- b*dnorm(qnorm(TPR1))/dnorm(qnorm(FPR0))
  k <- (v.tpr1/v.fpr0)^0.5/r
#  k <- ceiling(k)
  k
}

"CDF.Pval" <-
function(u, effect.size, n.sample, r.1, groups=2, type="balanced", grpj.per.grp1=1, distopt=1, control)
{
  .call. <- m <- match.call()
  pfx <- as.character(m[[1]])
    
  frmls <- names(formals(CDF.Pval))
  sppld <- names(m)[-1]

  m <- evald.call <- args.def.err.chk(m, sppld, frmls, other.rules=FALSE)
  err.msg <- attr(m, "err.msg")
    
  if(!is.null(err.msg)) stop(err.msg)

  u <- m$u
  groups <- m$groups
  r.1 <- m$r.1
  effect.size <- m$effect.size
  n.sample <- m$n.sample
  type <- m$type
  grpj.per.grp1 <- m$grpj.per.grp1
  control <- m$control
  
  nii.sample <- n.sample
  if(groups>=2 && type==3) nii.sample <- n.sample*grpj.per.grp1

  ## if(groups==1) type is 1 (paired)
  ## if(groups==2) type is 1 (paired), 2 (balanced) or 3 (unbalanced)
  ## if(groups>=3) type is 2 (balanced) or 3 (unbalanced)
  idistopt <- distopt
  .DF. <- eval(DF[type])
  .NCP. <- eval(NCP[type])
  
  pars0 <- eval(dists$pars0[[1+idistopt]])
  pars1 <- eval(dists$pars1[[1+idistopt]])
  
  u <- abs(u)
  is.pos <- (dists$minv[[1+idistopt]]==0)
  c.g <- dists$qdist[[1+idistopt]](1-u/2^(!is.pos), pars0)
  ans <- (1-r.1)*u + r.1*(1-dists$pdist[[1+idistopt]](c.g, pars1) + (!is.pos)*dists$pdist[[1+idistopt]](-c.g, pars1))
  df <- as.data.frame(list(u=u, CDF.Pval=ans))
  out <- list(CDF.Pval=df, call=.call.)
  attr(out, "passback") <- list(pars=list(pars0=pars0,pars1=pars1))
  
  class(out) <- "cdf"
  out
}

"CDF.Pval.HA" <-
function(u, effect.size, n.sample, r.1, groups=2, type="balanced", grpj.per.grp1=1, distopt=1, control)
{
  .call. <- m <- match.call()
  pfx <- as.character(m[[1]])
    
  frmls <- names(formals(CDF.Pval.HA))
  sppld <- names(m)[-1]

  m <- evald.call <- args.def.err.chk(m, sppld, frmls, other.rules=FALSE)
  err.msg <- attr(m, "err.msg")
    
  if(!is.null(err.msg)) stop(err.msg)

  u <- m$u
  groups <- m$groups
  r.1 <- m$r.1
  effect.size <- m$effect.size
  n.sample <- m$n.sample
  type <- m$type
  grpj.per.grp1 <- m$grpj.per.grp1
  control <- m$control

  nii.sample <- n.sample
  if(groups>=2 && type==3) nii.sample <- n.sample*grpj.per.grp1

  ## if(groups==1) type is 1 (paired)
  ## if(groups==2) type is 1 (paired), 2 (balanced) or 3 (unbalanced)
  ## if(groups>=3) type is 2 (balanced) or 3 (unbalanced)
  idistopt <- distopt  
  .DF. <- eval(DF[type])
  .NCP. <- eval(NCP[type])
  
  pars0 <- eval(dists$pars0[[1+idistopt]])
  pars1 <- eval(dists$pars1[[1+idistopt]])

  u <- abs(u)

  is.pos <- (dists$minv[[1+idistopt]]==0)
  c.g <- dists$qdist[[1+idistopt]](1-u/2^(!is.pos), pars0)
  ans <- 1-dists$pdist[[1+idistopt]](c.g, pars1) + (!is.pos)*dists$pdist[[1+idistopt]](-c.g, pars1)
  df <- as.data.frame(list(u=u, CDF.Pval.HA=ans))
  out <- list(CDF.Pval.HA=df, call=.call.)
  attr(out, "passback") <- list(pars=list(pars0=pars0,pars1=pars1))
  
  class(out) <- "cdf"
  out  
}

"CDF.Pval.au.eq.u" <-
function(effect.size, n.sample, r.1, alpha, groups=2, type="balanced", grpj.per.grp1=1, distopt=1, control)
{
    calling.env <- sys.parent()

    eval.env <- ifelse(calling.env==0, topenv, sys.parent)    
    m <- .call. <- match.call()

    n.args <- length(m)-1
    arg.nms <- names(m)[-1]

    is.rqd <- all(nera %in% arg.nms)
    if(!is.rqd) stop(list.a(nera, "Arguments ", " are required by CDF.Pval.au.eq.u"))

    if(missing(groups)) m$groups <- 2
    if(missing(type)) m$type <- 2
    if(missing(grpj.per.grp1)) m$grpj.per.grp1 <- 1
    
    if(missing(control)) m$control <- list(distopt=1, tol=1e-8)
    
    alpha <- eval(m$alpha, eval.env())
    groups <- eval(m$groups, eval.env())
    r.1 <- eval(m$r.1, eval.env())
    effect.size <- eval(m$effect.size, eval.env())
    n.sample <- eval(m$n.sample, eval.env())
    type <- eval(m$type, eval.env())
    grpj.per.grp1 <- eval(m$grpj.per.grp1, eval.env())
    control <- eval(m$control, eval.env())
        
    idistopt <- distopt

    nii.sample <- n.sample
    if(groups>=2 && type==3) nii.sample <- n.sample*grpj.per.grp1

    ## if(groups==1) type is 1 (paired)
    ## if(groups==2) type is 1 (paired), 2 (balanced) or 3 (unbalanced)
    ## if(groups>=3) type is 2 (balanced) or 3 (unbalanced)
    .DF. <- eval(DF[type])
    .NCP. <- eval(NCP[type])

    pars0 <- eval(dists$pars0[[1+idistopt]])
    pars1 <- eval(dists$pars1[[1+idistopt]])

    is.pos <- (dists$minv[[1+idistopt]]==0)
    
    conv <- is.na(alpha)
    gma.old <- r.1
    gma.new <- 1

    while(!conv)
    {
       u <- gma.old*alpha
       c.g <- dists$qdist[[1+idistopt]](1-u/2^(!is.pos), pars0)
       gma.new <- (1-r.1)*u + r.1*(1-dists$pdist[[1+idistopt]](c.g, pars1) +
                                       (!is.pos)*dists$pdist[[1+idistopt]](-c.g, pars1))
       obj <- abs(gma.new-gma.old)
       conv <- (obj < control$tol)
       gma.old <- gma.new
    }
    gamma <- gma.new
    out <- list(gamma=gamma, call=.call.)
    attr(out, "passback") <- list(pars=list(pars0=pars0,pars1=pars1))
    
    class(out) <- "cdf"
    out
}

"CDF.Pval.apsi.eq.u" <-
function(effect.size, n.sample, r.1, alpha, delta, groups=2, type="balanced", grpj.per.grp1=1, distopt=1, control)
{
    calling.env <- sys.parent()

    eval.env <- ifelse(calling.env==0, topenv, sys.parent)    
    m <- .call. <- match.call()

    n.args <- length(m)-1
    arg.nms <- names(m)[-1]

    is.rqd <- all(nera %in% arg.nms)
    if(!is.rqd) stop(list.a(nera, "Arguments ", " are required by CDF.Pval.apsi.eq.u"))
    
    if(missing(groups)) m$groups <- 2
    if(missing(delta)) m$delta <- alpha
    if(missing(type)) m$type <- 2
    if(missing(grpj.per.grp1)) m$grpj.per.grp1 <- 1
    
    if(missing(control)) m$control <- list(distopt=1, tol=1e-8)
    
    alpha <- eval(m$alpha, eval.env())
    delta <- eval(m$delta, eval.env())
    groups <- eval(m$groups, eval.env())
    r.1 <- eval(m$r.1, eval.env())
    effect.size <- eval(m$effect.size, eval.env())
    n.sample <- eval(m$n.sample, eval.env())
    type <- eval(m$type, eval.env())
    grpj.per.grp1 <- eval(m$grpj.per.grp1, eval.env())
    control <- eval(m$control, eval.env())
    
    idistopt <- distopt
    
    nii.sample <- n.sample
    if(groups>=2 && type==3) nii.sample <- n.sample*grpj.per.grp1

    ## if(groups==1) type is 1 (paired)
    ## if(groups==2) type is 1 (paired), 2 (balanced) or 3 (unbalanced)
    ## if(groups>=3) type is 2 (balanced) or 3 (unbalanced)
    .DF. <- eval(DF[type])
    .NCP. <- eval(NCP[type])
    
    pars0 <- eval(dists$pars0[[1+idistopt]])
    pars1 <- eval(dists$pars1[[1+idistopt]])
  
    is.pos <- (dists$minv[[1+idistopt]]==0)
    
    conv <- is.na(alpha)
    gma.old <- r.1
    gma.new <- 1

    while(!conv)
    {
        u <- psi(gma.old, delta)*alpha
        c.g <- dists$qdist[[1+idistopt]](1-u/2^(!is.pos), pars0)
        gma.new <- (1-r.1)*u + r.1*(1-dists$pdist[[1+idistopt]](c.g, pars1) +
                                                  (!is.pos)*dists$pdist[[1+idistopt]](-c.g, pars1))
        obj <- abs(gma.new-gma.old)
        conv <- (obj < control$tol)
        gma.old <- gma.new
    }
    gamma <- gma.new
    out <- list(gamma=gamma, call=.call.)
    attr(out, "passback") <- list(pars=list(pars0=pars0,pars1=pars1))
    
    class(out) <- "cdf"
    out
}

"sd.rtm.Rom" <-
function(object) 
{
  fn.nm <- as.character(object$call[[1]])
  frml.vals <- formals(get(fn.nm))
  frmls <- names(frml.vals)
  sppld <- names(object$call)[-1]
  is.msng <- !(frmls %in% sppld)
  names(is.msng) <- frmls
  if(is.msng["method"]) object$call$method <- "Theoretical"
  method <- match.arg(object$call$method, eval(frml.vals$method))
  sfx <- NULL
  if(method=="simulation")
  {
    fdpsfx <- c(`BHFDR`="",`BHFDX`=".st", `Romano`=".R")
    if(is.msng["FDP.control.method"]) object$call$FDP.control.method <- "BHFDR"
    sfx <- fdpsfx[object$call$FDP.control.method]
  }
  if(!is.msng["N.tests"]) ans <- object[["se.Rom" %,% sfx]]*eval(object$call$N.tests)^0.5
  if(is.msng["N.tests"]) ans <- object[["sigma.rtm.Rom" %,% sfx]]
  ans
}

"sd.rtm.VoR" <-
function(object)
{
  fn.nm <- as.character(object$call[[1]])
  frml.vals <- formals(get(fn.nm))
  frmls <- names(frml.vals)
  sppld <- names(object$call)[-1]
  is.msng <- !(frmls %in% sppld)
  names(is.msng) <- frmls
  if(is.msng["method"]) object$call$method <- "Theoretical"
  method <- match.arg(object$call$method, eval(frml.vals$method))
  sfx <- NULL
  if(method=="simulation")
  {
    fdpsfx <- c(`BHFDR`="",`BHFDX`=".st", `Romano`=".R")
    if(is.msng["FDP.control.method"]) object$call$FDP.control.method <- "BHFDR"
    sfx <- fdpsfx[object$call$FDP.control.method]
  }
  if(!is.msng["N.tests"]) ans <- object[["se.VoR" %,% sfx]]*eval(object$call$N.tests)^0.5
  if(is.msng["N.tests"]) ans <- object[["sigma.rtm.VoR" %,% sfx]]
  ans
}

"sd.rtm.ToM" <-
function(object)
{
  fn.nm <- as.character(object$call[[1]])
  frml.vals <- formals(get(fn.nm))
  frmls <- names(frml.vals)
  sppld <- names(object$call)[-1]
  is.msng <- !(frmls %in% sppld)
  names(is.msng) <- frmls
  if(is.msng["method"]) object$call$method <- "Theoretical"
  method <- match.arg(object$call$method, eval(frml.vals$method))
  sfx <- NULL
  if(method=="simulation")
  {
    fdpsfx <- c(`BHFDR`="",`BHFDX`=".st", `Romano`=".R")
    if(is.msng["FDP.control.method"]) object$call$FDP.control.method <- "BHFDR"
    sfx <- fdpsfx[object$call$FDP.control.method]
  }
  if(!is.msng["N.tests"]) ans <- object[["se.ToM" %,% sfx]]*eval(object$call$N.tests)^0.5
  if(is.msng["N.tests"]) ans <- object[["sigma.rtm.ToM" %,% sfx]]
  ans
}

# idea.
#    P( V_m / R_m  > \lambda )
#  = P( m^0.5 ( V_m / R_m - (1-r)alpha )/v^0.5  > m^0.5(\lambda - (1-r)alpha)/(v/m)^0.5
#  = 1 - pnorm(m^0.5(\lambda - (1-r)alpha)/v^0.5)
#
# if you want to control the tail probability instead of the expected value, set the above equal to (1-r) alpha
# and invert for \lambda, obtaining

# (1-r)alpha =  1 - pnorm(m^0.5(\lambda - (1-r)alpha)/v^0.5)
# pnorm(m^0.5(\lambda - (1-r)alpha)/v^0.5) = 1-(1-r)alpha
# m^0.5(\lambda - (1-r)alpha)/v^0.5 = qnorm(1-(1-r)alpha)
# \lambda_{r,alpha} = (1-r)alpha + (v/m)^0.5 qnorm(1-(1-r)alpha)

# so BH-FDR procedure at FDR=alpha bounds the probability that V_m / R_m  is in excess of \lambda_{r, alpha}
# by (1-r) alpha. This is good to have the probability bounded instead of the expected value, but what
# if \lambda_{r,alpha} is unacceptably large. For alpha < 1/2, \lambda_{r,alpha} > alpha 

# how about finding alpha* so that for v=v(alpha*)

# (1-r)alpha = \lambda_{r,alpha*} = (1-r)alpha* + (v(alpha*)/m)^0.5 qnorm(1-(1-r)alpha*)


# algorithm
# 1. consider r, theta=effect.size, n.sample, known
# 2. stipulate alpha
# 3. find alpha* such that:  (1-r)alpha = (1-r)alpha* + (v(alpha*)/m)^0.5 qnorm(1-(1-r)alpha*)


"controlFDP" <-
function(effect.size, n.sample, r.1, alpha, delta, groups=2, N.tests, 
         type=2, grpj.per.grp1=1, corr.struct=NULL, control,
         formula, data, distopt=1)
{
    calling.env <- sys.parent()
    eval.env <- ifelse(calling.env==0, topenv, sys.parent)
    m <- .call. <- match.call()
    
    m.sv <- m
    n.args <- length(m)-1
    arg.nms <- names(m)[-1]

    if(missing(type)) m$type <- type <- 2
    if(missing(distopt)) m$distopt <- distopt <- 1
    if(missing(groups)) m$groups <- groups <- 2
    if(missing(grpj.per.grp1)) m$grpj.per.grp1 <- 1
        
    nii.sample <- n.sample
    if(groups>=2 && type==3) nii.sample <- n.sample*grpj.per.grp1

    is.form <- !missing(formula)
    if(is.form)
    {
      m[[1]] <- as.name("controlFDP.formula")  
      dfa <- c("data","formula","alpha")
      is.rqd <- all(dfa %in% arg.nms)
      if(!is.rqd) stop(list.a(dfa, "Aguments ", " are required by the formula " %,%
                                                "'controlFDP' method "))
      if(missing(delta)) m$delta <- alpha
      if(missing(control)) m$control <- list(verb=0, tol=1e-8)
    }
    if(!is.form)
    {
      m[[1]] <- as.name("controlFDP.theoret")
      is.rqd <- all(Nnera %in% arg.nms)
      if(!is.rqd) stop(list.a(Nnera, "Arguments ", " are required by the asymptotic " %,%
                                                   "'controlFDP' method"))
      if(missing(delta)) m$delta <- alpha
      if(missing(control)) m$control <- list(verb=0, tol=1e-8)
    }

    out <- eval(m)
    out$call <- .call.
    class(out) <- "vvv"
    out   
}

"controlFDP.formula" <-
function(effect.size, n.sample, r.1, alpha, delta, groups=2, N.tests, 
         type, grpj.per.grp1, corr.struct, control, formula, data, distopt)
{
    .call. <- m <- match.call()

    calling.env <- sys.parent()
    eval.env <- ifelse(calling.env==0, topenv, sys.parent)

    nii.sample <- n.sample
    if(groups>=2 && type==3) nii.sample <- n.sample*grpj.per.grp1

    mf <- as.call(expression(model.frame))
    mf$groups <- mf$effect.size <- mf$n.sample <- mf$r.1 <- mf$N.tests <- mf$control <- mf$delta <- NULL
    mf$formula <- formula
    mf$data <- data
    mf <- eval(mf, eval.env())
    pval <- model.extract(mf, "response")
    N.tests <- length(pval)
    pvalgthlf <- 1*(pval > 0.5)
    r.0 <- 2*mean(pvalgthlf)
    r.1 <- 1-r.0

    conv <- neg <- FALSE
    msg <- NULL
    ast.old <- alpha
    a.star <- gma.st <- obj <- sqrt.v <- P.star <- NA
    
    is.cs <- !missing(corr.struct)
    var.VoR.call <- as.call(expression(var.VoR))
    var.VoR.call$r.1 <- m$r.1
    var.VoR.call$delta <- m$delta
    var.VoR.call$n.sample <- n.sample
    var.VoR.call$nii.sample <- nii.sample
    var.VoR.call$groups <- groups 
    var.VoR.call$effect.size <- effect.size
    var.VoR.call$type <- type
    var.VoR.call$distopt <- distopt
    var.VoR.call$FDP.control.method <- "BHFDX"

    if(is.cs) var.VoR.call$corr.struct <- corr.struct  

    while(!(conv||neg))
    {
        gma.st <- max(which(pval <= (1:N.tests)/N.tests*ast.old))/N.tests
        var.VoR.call$gamma <- gma.st
        var.VoR.call$alpha <- ast.old
        v <- eval(var.VoR.call)/N.tests
        ast.new <- (delta - v^0.5*qnorm(1-alpha))/r.0
        neg <- ast.new < 0 
        if(control$verb>1) cat(sprintf("ast.new=%f\n", ast.new))
        obj <- abs(ast.new-ast.old)
        conv <- (obj < 1e-4*alpha)
        ast.old <- ast.new
    }
    if(neg)
    {
        a.star <- gma.st <- obj <- sqrt.v <- P.star <- NA
        msg <- "Solution not attainable."
    }
    if(!neg)
    {
        a.star <- ast.new
        gma.st <- max(which(pval <= (1:N.tests)/N.tests*ast.old))/N.tests
        var.VoR.call$gamma <- gma.st
        var.VoR.call$alpha <- a.star

        v <- eval(var.VoR.call)/N.tests
        sqrt.v <- v^0.5
        P.star <- 1-pnorm((delta - r.0*a.star)/sqrt.v)
    }
    out <- list(alpha.star = a.star, gamma=gma.st, obj = obj, se.FDF=sqrt.v, P.star = P.star)
    if(!is.null(msg)) out$note <- msg
    out
}

"controlFDP.theoret" <-
function(effect.size, n.sample, r.1, alpha, delta, groups=2, N.tests, 
         type, grpj.per.grp1, corr.struct, control, formula, data, distopt)
{
    .call. <- m <- match.call()
    
    calling.env <- sys.parent()
    eval.env <- ifelse(calling.env==0, topenv, sys.parent)

    nii.sample <- n.sample
    if(groups>=2 && type==3) nii.sample <- n.sample*grpj.per.grp1

    ast.old <- alpha
    conv <- neg <- FALSE
    msg <- NULL
    a.star <- gma.st <- obj <- sqrt.v <- P.star <- NA

    gma.call <- match.call()
    gma.call[[1]] <- as.name("CDF.Pval.au.eq.u")
    
    gma.call$formula <- gma.call$data <- gma.call$corr.struct <- gma.call$N.tests <- gma.call$delta <- NULL
    gma0 <- eval(gma.call)
    r.0 <- 1-r.1

    is.cs <- !missing(corr.struct)
    var.VoR.call <- as.call(expression(var.VoR))
    var.VoR.call$r.1 <- m$r.1
    var.VoR.call$delta <- m$delta
    var.VoR.call$n.sample <- n.sample
    var.VoR.call$nii.sample <- nii.sample
    var.VoR.call$groups <- groups 
    var.VoR.call$effect.size <- effect.size
    var.VoR.call$type <- type
    var.VoR.call$distopt <- distopt
    var.VoR.call$FDP.control.method <- "BHFDX"
    
    if(is.cs) var.VoR.call$corr.struct <- corr.struct  

    while(!(conv || neg))
    {
        gma.st <- update(gma0, alpha=ast.old)$gamma
        var.VoR.call$alpha <- ast.old
        var.VoR.call$gamma <- gma.st                
        v <- eval(var.VoR.call)/N.tests
        ast.new <- (delta - v^0.5*qnorm(1 - alpha))/r.0
        neg <- ast.new < 0 
        if(control$verb>1) cat(sprintf("ast.new=%f\n", ast.new))
        obj <- abs(ast.new-ast.old)
        conv <- (obj < 1e-4*alpha)
        ast.old <- ast.new
    }
    if(neg)
    {
        a.star <- gma.st <- obj <- sqrt.v <- P.star <- NA
        msg <- "Solution not attainable."
    }
    if(!neg)
    {
        a.star <- ast.new
        gma.st <- update(gma0, alpha=a.star)$gamma
        var.VoR.call$alpha <- a.star
        var.VoR.call$gamma <- gma.st                
        v <- eval(var.VoR.call)/N.tests
        sqrt.v <- v^0.5
        P.star <- 1-pnorm((delta - r.0*a.star)/sqrt.v)
    }
    out <- list(alpha.star = a.star, gamma=gma.st, obj = obj, se.FDF=sqrt.v, P.star = P.star)
    if(!is.null(msg)) out$note <- msg
    out
}

RR <-
function(ss, tt, u, v, r.1, corr.struct, n.sample, nii.sample, groups, effect.size, type,
         distopt, FDP.control.method)
{
    ## if(groups==1) type is 1 (paired)
    ## if(groups==2) type is 1 (paired), 2 (balanced) or 3 (unbalanced)
    ## if(groups>=3) type is 2 (balanced) or 3 (unbalanced)
    .DF. <- eval(DF[type])
    .NCP. <- eval(NCP[type])
    idistopt <- distopt
    is.pos <- (dists$minv[[1+idistopt]]==0)
    pars0 <- eval(dists$pars0[[1+idistopt]])
    pars1 <- eval(dists$pars1[[1+idistopt]])
    pars.u <- (1-u)*pars0 + u*pars1
    pars.v <- (1-v)*pars0 + v*pars1

    r.0 <- 1-r.1
    r.u <- (1-u)*r.0 + u*r.1
    r.v <- (1-v)*r.0 + v*r.1
    
    #### from var.VoR
    is.cs <- !missing(corr.struct)
    z.s <- dists$qdist[[1+idistopt]](1-ss/2^(!is.pos), pars0)
    z.t <- dists$qdist[[1+idistopt]](1-tt/2^(!is.pos), pars0)
    z.minst <- dists$qdist[[1+idistopt]](1-pmin(ss,tt)/2^(!is.pos), pars0)
      
    cPhi.zs.u <- (1-dists$pdist[[1+idistopt]](z.s, pars.u)) + (!is.pos)*dists$pdist[[1+idistopt]](-z.s, pars.u)
    cPhi.zt.v <- (1-dists$pdist[[1+idistopt]](z.t, pars.v)) + (!is.pos)*dists$pdist[[1+idistopt]](-z.t, pars.v)
    cPhi.zminst.u <- (1-dists$pdist[[1+idistopt]](z.minst, pars.u)) + (!is.pos)*dists$pdist[[1+idistopt]](-z.minst, pars.u)
    
    sumGMA <- 0
    if(is.cs)
    {
      cty <- corr.struct$type
      k.bs <- corr.struct$block.size
      rhoX <- corr.struct$rho
        
      if(cty=="CS-Blocks")
      {
        RHO <- matrix(c(1, rhoX, rhoX, 1), 2, 2)
          
        ## Outer 4 corners of tic tac toe board. Clockwise, starting from upper right.
        ## Written as double integrals against the bivariate density, we have 
        ## UR:  int_{z.s - u MU}^{Inf} int_{z.t - v MU}{Inf}
        ## LR:  int_{z.s - u MU}^{Inf} int_{-Inf}{-z.t - v MU}
        ## LL:  int_{-Inf}^{-z.s - u MU} int_{-Inf}{-z.t - v MU}
        ## UL:  int_{-Inf}^{-z.s - u MU} int_{z.t - v MU}{Inf}
        ## The same comments apply below in the "Toeplitz-Blocks" section
        
        tails.mvn <- pmvnorm(lower=c(z.s-u*.NCP., z.t-v*.NCP.), corr=RHO) + ## UR
           (!is.pos)*pmvnorm(lower=c(z.s-u*.NCP., -Inf), upper=c(Inf, -z.t-v*.NCP.), corr=RHO) + ## LR
           (!is.pos)*pmvnorm(upper=c(-z.s-u*.NCP., -z.t-v*.NCP.), corr=RHO) + ## LL
           (!is.pos)*pmvnorm(lower=c(-Inf, z.t-v*.NCP.), upper=c(-z.s-u*.NCP., Inf), corr=RHO) ## UL

        GMA <- tails.mvn - cPhi.zs.u*cPhi.zt.v
        sumGMA <- (k.bs-1)*GMA
      }
      if(cty=="Toeplitz-Blocks")
      {
        j <- 1:(k.bs-1)
        GMA <- NULL
        for(i in 1:(k.bs-1))
        {
          RHO <- matrix(c(1,rhoX[i],rhoX[i],1),2,2)
            
          tails.mvn <- pmvnorm(lower=c(z.s-u*.NCP., z.t-v*.NCP.), corr=RHO) + ## UR
             (!is.pos)*pmvnorm(lower=c(z.s-u*.NCP., -Inf), upper=c(Inf, -z.t-v*.NCP.), corr=RHO) + ## LR
             (!is.pos)*pmvnorm(upper=c(-z.s-u*.NCP., -z.t-v*.NCP.), corr=RHO) + ## LL
             (!is.pos)*pmvnorm(lower=c(-Inf, z.t-v*.NCP.), upper=c(-z.s-u*.NCP., Inf), corr=RHO) ## UL

          GMA <- c(GMA, tails.mvn - cPhi.zs.u*cPhi.zt.v)
        }
        sumGMA <- 2/k.bs*sum((k.bs-j)*GMA)
      }
    } 
    ##
    I(u==v)*r.u*cPhi.zminst.u - r.u*r.v*cPhi.zs.u*cPhi.zt.v + r.u*r.v*sumGMA
}

var.VoR <-
function(gamma, alpha, delta, r.1, corr.struct, n.sample, nii.sample, groups, effect.size, type,
         distopt, FDP.control.method)
{
   ## if(groups==1) type is 1 (paired)
   ## if(groups==2) type is 1 (paired), 2 (balanced) or 3 (unbalanced)
   ## if(groups>=3) type is 2 (balanced) or 3 (unbalanced)
   .DF. <- eval(DF[type])
   .NCP. <- eval(NCP[type])

   r.0 <- 1-r.1

   flg.psi <- eval(flg.psi.expr)
   gma <- gamma 
   psi. <- gamma
   psi.pr <- 1
   if(flg.psi)
   {
       psi. <- delta*gamma/(1-(1-delta)*gamma)
       psi.pr <- delta/(1-(1-delta)*gamma)^2
   }
    
   idistopt <- distopt
   is.pos <- (dists$minv[[1+idistopt]]==0)

   pars0 <- eval(dists$pars0[[1+idistopt]])
   pars1 <- eval(dists$pars1[[1+idistopt]])

   c.g <- qdist(1 - alpha*psi./2^(!is.pos), pars0)
   G1.pr <- r.1*(ddist(c.g, pars1) + (!is.pos)*ddist(-c.g, pars1))/(ddist(c.g, pars0) + (!is.pos)*ddist(-c.g, pars0))
   G.pr <- r.0 + G1.pr

   denom <- 1-alpha*psi.pr*G.pr
   aa <- r.0*alpha*(psi.pr - psi./gma)/denom
   bb <- 1+aa

   R.00 <- RR(alpha*psi., alpha*psi., 0, 0, r.1, corr.struct, n.sample, nii.sample, groups, effect.size, type,
         distopt, FDP.control.method)
   R.01 <- RR(alpha*psi., alpha*psi., 0, 1, r.1, corr.struct, n.sample, nii.sample, groups, effect.size, type,
         distopt, FDP.control.method)
   R.11 <- RR(alpha*psi., alpha*psi., 1, 1, r.1, corr.struct, n.sample, nii.sample, groups, effect.size, type,
         distopt, FDP.control.method)
    
   v.VoR <- (bb^2*R.00 + 2*bb*aa*R.01 + aa^2*R.11)/gamma^2
   v.VoR
}

var.Rom <-
function(gamma, alpha, delta, r.1, corr.struct, n.sample, nii.sample, groups, effect.size, type,
         distopt, FDP.control.method)
{
   ## if(groups==1) type is 1 (paired)
   ## if(groups==2) type is 1 (paired), 2 (balanced) or 3 (unbalanced)
   ## if(groups>=3) type is 2 (balanced) or 3 (unbalanced)
   .DF. <- eval(DF[type])
   .NCP. <- eval(NCP[type])

   r.0 <- 1-r.1

   flg.psi <- eval(flg.psi.expr)
   gma <- gamma 
   psi. <- gamma
   psi.pr <- 1
   if(flg.psi)
   {
       psi. <- delta*gamma/(1-(1-delta)*gamma)
       psi.pr <- delta/(1-(1-delta)*gamma)^2
   }

   idistopt <- distopt
   is.pos <- (dists$minv[[1+idistopt]]==0)

   pars0 <- eval(dists$pars0[[1+idistopt]])
   pars1 <- eval(dists$pars1[[1+idistopt]])

   c.g <- qdist(1 - alpha*psi./2^(!is.pos), pars0)
   G1.pr <- r.1*(ddist(c.g, pars1) + (!is.pos)*ddist(-c.g, pars1))/(ddist(c.g, pars0) + (!is.pos)*ddist(-c.g, pars0))
   G.pr <- r.0 + G1.pr

   R.00 <- RR(alpha*psi., alpha*psi., 0, 0, r.1, corr.struct, n.sample, nii.sample, groups, effect.size, type,
         distopt, FDP.control.method)
   R.01 <- RR(alpha*psi., alpha*psi., 0, 1, r.1, corr.struct, n.sample, nii.sample, groups, effect.size, type,
         distopt, FDP.control.method)
   R.11 <- RR(alpha*psi., alpha*psi., 1, 1, r.1, corr.struct, n.sample, nii.sample, groups, effect.size, type,
         distopt, FDP.control.method)
    
   denom <- 1-alpha*psi.pr*G.pr
   v.Rom <- (R.00 + 2*R.01 + R.11)/denom^2
   v.Rom
}

var.ToM <-
function(gamma, alpha, delta, r.1, corr.struct, n.sample, nii.sample, groups, effect.size, type,
         distopt, FDP.control.method)
{
   ## if(groups==1) type is 1 (paired)
   ## if(groups==2) type is 1 (paired), 2 (balanced) or 3 (unbalanced)
   ## if(groups>=3) type is 2 (balanced) or 3 (unbalanced)
   .DF. <- eval(DF[type])
   .NCP. <- eval(NCP[type])

   r.0 <- 1-r.1

   flg.psi <- eval(flg.psi.expr)
   gma <- gamma 
   psi. <- gamma
   psi.pr <- 1
   if(flg.psi)
   {
       psi. <- delta*gamma/(1-(1-delta)*gamma)
       psi.pr <- delta/(1-(1-delta)*gamma)^2
   }

   idistopt <- distopt
   is.pos <- (dists$minv[[1+idistopt]]==0)

   pars0 <- eval(dists$pars0[[1+idistopt]])
   pars1 <- eval(dists$pars1[[1+idistopt]])

   c.g <- qdist(1 - alpha*psi./2^(!is.pos), pars0)
   pi.1 <- 1 - pdist(c.g, pars1) + (!is.pos) * pdist(-c.g, pars1)

   G1.pr <- r.1*(ddist(c.g, pars1) + (!is.pos)*ddist(-c.g, pars1))/(ddist(c.g, pars0) + (!is.pos)*ddist(-c.g, pars0))
   G.pr <- r.0 + G1.pr

   R.00 <- RR(alpha*psi., alpha*psi., 0, 0, r.1, corr.struct, n.sample, nii.sample, groups, effect.size, type,
         distopt, FDP.control.method)
   R.01 <- RR(alpha*psi., alpha*psi., 0, 1, r.1, corr.struct, n.sample, nii.sample, groups, effect.size, type,
         distopt, FDP.control.method)
   R.11 <- RR(alpha*psi., alpha*psi., 1, 1, r.1, corr.struct, n.sample, nii.sample, groups, effect.size, type,
         distopt, FDP.control.method)

   R.11.11 <-  RR(1, 1, 1, 1, r.1, corr.struct, n.sample, nii.sample, groups, effect.size, type,
         distopt, FDP.control.method)
   R.01.az.1 <-  RR(alpha*psi., 1, 0, 1, r.1, corr.struct, n.sample, nii.sample, groups, effect.size, type,
         distopt, FDP.control.method)
   R.11.az.1 <-  RR(alpha*psi., 1, 1, 1, r.1, corr.struct, n.sample, nii.sample, groups, effect.size, type,
         distopt, FDP.control.method)
   
   denom <- 1-alpha*psi.pr*G.pr
   tau2 <- (R.00 + 2*R.01 + R.11)/denom^2

   var.Y1 <- R.11 + alpha^2*psi.pr^2*G1.pr^2*tau2 + 2*alpha*psi.pr*G1.pr*(R.11 + R.01)/denom
   var.W1.1 <- R.11.11
   cov.Y1.W1.1 <- R.11.az.1 + alpha*psi.pr*G1.pr*(R.01.az.1+R.11.az.1)/denom
   
   v.ToM <- 1/r.1^2*(var.Y1 + pi.1^2*var.W1.1 - 2*pi.1*cov.Y1.W1.1)
   v.ToM
}

"cCDF.Rom" <-
function(u, effect.size, n.sample, r.1, alpha, delta, groups=2, N.tests,
         type=c("paired","balanced","unbalanced"), grpj.per.grp1=NULL, FDP.control.method="BHFDR", 
         distopt, control=list(tol=1e-8, max.iter=c(1000,20), sim.level=2,
                      low.power.stop=TRUE, FDP.meth.thresh=FDP.cntl.mth.thrsh.def, verb=FALSE))
{
  m <- .call. <- match.call()
  m[[1]] <- as.name("pwrFDR")
  m$u <- NULL  
  avgpwr <- eval(m, sys.parent())
  gma <- avgpwr$gamma
  sdrtmRom <- sd.rtm.Rom(avgpwr)
  ans <- 1-pnorm(N.tests^0.5*(u - gma)/sdrtmRom)
  df <- as.data.frame(list(u=u, cCDF.Rom=ans))
  out <- list(cCDF.Rom=df, call=.call.)
  class(out) <- "cdf"
  out
}

"cCDF.ToM" <-
function(u, effect.size, n.sample, r.1, alpha, delta, groups=2, N.tests,
         type=c("paired","balanced","unbalanced"), grpj.per.grp1=NULL, FDP.control.method="BHFDR", 
         distopt, control=list(tol=1e-8, max.iter=c(1000,20), sim.level=2,
                      low.power.stop=TRUE, FDP.meth.thresh=FDP.cntl.mth.thrsh.def, verb=FALSE))
{
  m <- .call. <- match.call()
  m[[1]] <- as.name("pwrFDR")
  m$u <- NULL  
  avgpwr <- eval(m, sys.parent())
  pi <- avgpwr$average.power
  sdrtmToM <- sd.rtm.ToM(avgpwr)
  ans <- 1-pnorm(N.tests^0.5*(u - pi)/sdrtmToM)
  df <- as.data.frame(list(lambda=u, cCDF.ToM=ans))
  out <- list(cCDF.ToM=df, call=.call.)
  class(out) <- "cdf"
  out
}

"cCDF.VoR" <-
function(u, effect.size, n.sample, r.1, alpha, delta, groups=2, N.tests, 
         type=c("paired","balanced","unbalanced"), grpj.per.grp1=NULL, FDP.control.method="BHFDR", 
         distopt, control=list(tol=1e-8, max.iter=c(1000,20), sim.level=2,
                      low.power.stop=TRUE, FDP.meth.thresh=FDP.cntl.mth.thrsh.def, verb=FALSE))
{
  m <- .call. <- match.call()
  m[[1]] <- as.name("pwrFDR")
  m$u <- NULL  
  avgpwr <- eval(m, sys.parent())
  sdrtmVoR <- sd.rtm.VoR(avgpwr)
  r.0 <- 1-r.1
  do.bhclt <- do.R <- FALSE
  verb <- control$verb
  if(FDP.control.method=="BHFDX") do.bhclt <- TRUE
  if(FDP.control.method=="Romano") do.R <- TRUE
  if(FDP.control.method=="Auto")
  {
      if(names(avgpwr$FDP.cnt)=="BHFDX") do.bhclt <- TRUE
      if(names(avgpwr$FDP.cnt)=="Romano") do.R <- TRUE
  }
  alpha. <- alpha
  psi.o.gma <- 1
  if((do.bhclt || do.R) && missing(delta)) delta <- attr(avgpwr, "arg.vals")$delta
  if(do.bhclt) alpha. <- avgpwr$alpha.star
  if(do.R)
  {
      psi <- function(gma,a)gma*a/(1-(1-a)*gma)
      gma <- avgpwr$gamma
      .psi. <- psi(gma, delta)
      psi.o.gma <- .psi./gma
  }
  ans <- 1-pnorm(N.tests^0.5*(u - r.0*psi.o.gma*alpha.)/sdrtmVoR)
  df <- as.data.frame(list(delta=u, cCDF.VoR=ans))
  out <- list(cCDF.VoR=df, call=.call.)
  class(out) <- "cdf"
  out
}

"detail" <-
function(obj)
{
  attr(obj, "detail")
}

"%,%" <- function(x,y)paste(x,y,sep="")
"DX" <- function(x)c(x[1],diff(x))
  
"ddist" <- 
function(x, pars)
{
  idistopt <- pars[1]
  dists$ddist[[1+idistopt]](x, pars)
}

"pdist" <-
function(x, pars)
{
  idistopt <- pars[1]
  dists$pdist[[1+idistopt]](x, pars)
}

"qdist" <-
function(x, pars)
{
  idistopt <- pars[1]
  dists$qdist[[1+idistopt]](x, pars)
}

"gentempfilenm" <-
function(prfx="temp", sfx=".txt")
{
    alphanum <- c(letters, toupper(letters), 0:9)
    n <- length(alphanum)
    prfx %,% paste(alphanum[sample(n, 5)], collapse="") %,% sfx
}

"print.pwrFDR" <-
function(x, result=c("raw","tex","html"), cptn=NULL, label=NULL, ...)
{
  if(arg.vals(x)$control$verb >=3) browser()
  ## original first line
  if(missing(result)) result <- "raw"
  fn.nm <- "printpwrFDR" %,% result  
  m <- as.call(expression(prnt, x, cptn=cptn, label=label, ...))
  m[[1]] <- as.name(fn.nm)
  eval(m)
}

"print.join.pwrFDR" <-
function(x, result=c("raw","tex","html"), cptn=NULL, label=NULL, ...)
{
  .call. <- match.call()
  if(missing(result)) result <- "raw"
  fn.nm <- "printjoinpwrFDR" %,% result
  m <- as.call(expression(prnt, x))
  if(!missing(cptn)) m$cptn <- cptn
  if(!missing(label)) m$label <- label
  m[[1]] <- as.name(fn.nm)
  eval(m)
}

"printpwrFDRraw" <-
function(x, ...)
{
  show.footer <- attr(x,"arg.vals")$control$show.footer
  obj <- x
  cat("Call:\n")
  print(x$call)
  x.df <- join.tbl(x)
  m <- as.call(expression(print, x.df, row.names=FALSE))
  eval(m)
  if(show.footer) writeLines(attr(x, "footer"))
  invisible(obj)
}

"printpwrFDRtex" <-
function(x, cptn=NULL, label=NULL, ...)
{
  obj <- x
  m <- match.call(expand.dots=TRUE)
  m[[1]] <- as.name("tm.pwrFDR")
  m$x <- as.call(expression(join.tbl, x))
  m$label <- label
  eval(m)
  invisible(obj)
}

"printpwrFDRhtml" <-
function(x, cptn=NULL, ...)
{
  obj <- x
  m <- as.call(expression(flextable))
  m$data <- as.call(expression(join.tbl, x))
  ft <- eval(m)
  if(!missing(cptn)) ft <- set_caption(ft, cptn)
  print(ft)
  invisible(obj)
}

"join.tbl" <-
function(...)
{
    m <- match.call(expand.dots=TRUE)
    n.args <- length(m) - 1    ## 4
    tbl.lst <- list()
    cls <- NULL
    nc <- NULL
    for(k in 1:n.args)  ## 1:4
    {
      tbl.lst[[k]] <- eval(m[[k+1]], sys.parent())
      cls <- c(cls, class(tbl.lst[[k]]))
      nc <- c(nc, length(tbl.lst[[k]]))
    }
    bad <- any(cls!="pwrFDR") || (length(unique(nc))>1)
    if(bad)stop("All arguments must be of class 'pwrFDR' and of the same length")
    TBL <- AsDataFrame_pwrFDR(tbl.lst[[1]])

    if(n.args>1)
    for(k in 2:n.args)
      TBL <- cbind(TBL, AsDataFrame_pwrFDR(tbl.lst[[k]])[,2])
    
    names(TBL)[2:(n.args+1)] <- "result " %,% letters[1:n.args]
    class(TBL) <- "join.pwrFDR"
    TBL
}

printjoinpwrFDRraw <-
function(x)
{
  class(x) <- "data.frame"
  print(x, row.names=FALSE)
  invisible(x)
}

printjoinpwrFDRtex <-
function(x, cptn=NULL, label=NULL, ...)
{
  class(x) <- "data.frame"
  m <- as.call(expression(basic.tmPrint, x))
  if(!missing(cptn)) m$cptn <- cptn
  if(!missing(label)) m$lbl <- label
  eval(m)
  invisible(x)
}
printjoinpwrFDRhtml <-
function(x, cptn=NULL, ...)
{
  class(x) <- "data.frame"
  m <- as.call(expression(flextable))
  m$data <- x
  ft <- eval(m)
  if(!missing(cptn)) ft <- set_caption(ft, cptn)
  print(ft)
  invisible(x)
}

AsDataFrame_pwrFDR <-
function(x)
{
    verb <- arg.vals(x)$control$verb
    if(verb>=3) browser()
    .call. <- x$call
    ind.call <- which(names(x)=="call")
    x <- x[-ind.call]
    TBL <- as.data.frame(list(names(x), unlist(cbind(x))))
    names(TBL) <- c("Parameter", "Value")
    vals <- vals.sv <- TBL$Value
    nms.vals <- TBL$Parameter

    is.inf <- abs(vals)==Inf
    idx.inf <- which(is.inf)
    any.inf <- any(is.inf)
    if(any.inf) nms.vals <- nms.vals[-idx.inf]
    
    vals <- vals[as.numeric(sapply(nms.vals, FUN=\(x, y)which(x==y), y=TBL$Parameter))]
    vals <- sapply(vals, 
                        function(x){
                          is.char <- 1*is.character(x)
                          idx <- is.int <- 0
                          if(!is.char)
                          {
                            is.int <- 2*(abs(x)-floor(x) < 1e-12)
                            idx <- (!is.int)*(3*(x-floor(x) <= 1e-3) + 4*(x-floor(x) > 1e-3))
                          }
                          idx <- is.char + is.int + idx
                          lbl <- c("char", "int","sci","dec")[idx]
                          int.expr <- as.call(expression(format, floor(x)))
                          sci.expr <- as.call(expression(format, x, scientific=TRUE, digits=4))
                          dec.expr <- as.call(expression(format, round(x, 4)))
                          switch(lbl,
                                 char=x,
                                 int=eval(int.expr),
                                 sci=eval(sci.expr),
                                 dec=eval(dec.expr))
                        })
    if(any.inf) vals <- append(vals, vals.sv[idx.inf], idx.inf-1)
    TBL$Value <- unlist(vals)

    ## translate the FDP.control.method from numeric back to character
    cntl.nms <- c(`1`="BHFDR", `2`="BHFDX", `3`="Romano", `4`="Hochberg", `5`="Holm", `6`="Bonferroni")
    TBL$Value[which(TBL$Parameter=="FDP.cnt")] <- cntl.nms[TBL$Value[which(TBL$Parameter=="FDP.cnt")]]
    TBL
}

tm.pwrFDR <-
function(x, cptn=NULL, label=NULL, ...)
{
    .call. <- match.call(expand.dots=TRUE)
    nms.call <- names(.call.)
    x.call. <- x$call

    vals <- x$Value
    is.inf <- vals=="Inf"
    idx.inf <- which(is.inf)
    any.inf <- any(is.inf)
    if(any.inf) x$Value[idx.inf] <- "$\\infty$"
    hdr <- as.list(rep("", 2))
    names(hdr) <- names(x)
    tmHeadings(x) <- hdr
    tmCtypes(x) <- c("c", "n")
    tmDigits(x) <- c(0, 4)
    if(!missing(cptn)) tmCaption(x) <- cptn
    class(x) <- "TableMonster"
    prnt <- as.call(expression(print, x))
    if(!missing(label)) prnt$label <- label
    prnt$sanitize.text.function <- I
    eval(prnt)
}

"print.vvv" <-
 function (x, ...) 
{
    y <- x
    cat("Call:\n")
    print(x$call)
    class(x) <- NULL
    x$call <- NULL
    x <- as.data.frame(x)
    dimnames(x)[[1]] <- " "
    print(x)
    invisible(x)
}

"print.cdf" <-
function(x, ...)
{
    y <- x
    cat("Call:\n")
    print(x$call)
    class(x) <- NULL
    x$call <- NULL
    print(x[[1]])
    invisible(x)    
}

"is.int" <-
function(x)
{
  abs(x - floor(x)) < 1e-10
}

`+.pwrFDR` <-
  function(x,y)
  {
    xpwr <- x
    ypwr <- y
    if(is(x,"pwrFDR")) xpwr <- ifelse(!is.null(x$call$lambda), x$TPX.power, x$average.power)
    if(is(y,"pwrFDR")) ypwr <- ifelse(!is.null(y$call$lambda), y$TPX.power, y$average.power)
    xpwr + ypwr
  }

`-.pwrFDR` <-
  function(x,y)
  {
    xpwr <- x
    ypwr <- y
    if(is(x,"pwrFDR")) xpwr <- ifelse(!is.null(x$call$lambda), x$TPX.power, x$average.power)
    if(is(y,"pwrFDR")) ypwr <- ifelse(!is.null(y$call$lambda), y$TPX.power, y$average.power)
    xpwr - ypwr
  }

`*.pwrFDR` <-
  function(x,y)
  {
    xpwr <- x
    ypwr <- y
    if(is(x,"pwrFDR")) xpwr <- ifelse(!is.null(x$call$lambda), x$TPX.power, x$average.power)
    if(is(y,"pwrFDR")) ypwr <- ifelse(!is.null(y$call$lambda), y$TPX.power, y$average.power)
    xpwr * ypwr
  }

`/.pwrFDR` <-
  function(x,y)
  {
    xpwr <- x
    ypwr <- y
    if(is(x,"pwrFDR")) xpwr <- ifelse(!is.null(x$call$lambda), x$TPX.power, x$average.power)
    if(is(y,"pwrFDR")) ypwr <- ifelse(!is.null(y$call$lambda), y$TPX.power, y$average.power)
    xpwr / ypwr
  }

`^.pwrFDR` <-
  function(x,y)
  {
    xpwr <- x
    ypwr <- y
    if(is(x,"pwrFDR")) xpwr <- ifelse(!is.null(x$call$lambda), x$TPX.power, x$average.power)
    if(is(y,"pwrFDR")) ypwr <- ifelse(!is.null(y$call$lambda), y$TPX.power, y$average.power)
    xpwr ^ ypwr
  }

logit <- function(mu)UseMethod("logit")
logit.pwrFDR <- function(mu)logit.default(ifelse(!is.null(mu$call$lambda), mu$TPX.power, mu$average.power))
logit.default <- function(mu)binomial("logit")$linkfun(mu)

logitInv <- function(eta)UseMethod("logitInv")
logitInv.default <- function(eta)binomial("logit")$linkinv(eta)
logitInv.pwrFDR <- function(eta)logitInv.default(ifelse(!is.null(eta$call$lambda), eta$TPX.power, eta$average.power))

exp.pwrFDR <- function(x) exp(ifelse(!is.null(x$call$lambda), x$TPX.power, x$average.power))

log.pwrFDR <- function(x, base)log(ifelse(!is.null(x$call$lambda), x$TPX.power, x$average.power), base)

arg.vals <- function(object)attr(object, "arg.vals")
"%over%" <-
function(x,y)
{
  dx <- dim(x)
  dy <- dim(y)
  bad <- any(dx!=dy) && dx[1]!=length(y)
  if(bad) stop("Incompatable dimensions")
  ans <- 0*x
  bad <- c(which(y==0),which(is.na(y)))
  if(length(bad)>0) ans[-bad] <- x[-bad]/y[-bad]
  if(length(bad)==0) ans <- x/y
  structure(ans, dim=dx)
}

##  Sample size and power for omnibus F-test, 'groups' groups
##
##  H_0  mu_j - mu_avg = 0
##  H_A  mu_j = mu_avg = -theta/2, 0, ..., 0, theta/2
##
##  NOTE: this will be the case (H_A) when
##  mu = (mu_0, mu_0 + theta/2, ..., mu_0 + theta/2, mu_0 + theta)
##  In this case, the ncp = n*(mu - mu_avg)^2/2 = n*theta^2/2
##  which is the same as the square of the two sample t-test.

"f.power" <-
function(n.sample, groups, effect.size, e.I, power)
{
  no.n <- missing(n.sample)
  no.g <- missing(groups)
  no.e <- missing(effect.size)
  no.p <- missing(power)
  if(no.g || no.e) stop("Arguements 'groups' and 'effect.size' are required")
  bad <- (no.n && no.p) || (!no.n && !no.p)
  if(bad) stop("You must specify _either_ 'n.sample' or 'power'")
  do.p <- no.p
  do.ss <- no.n
  if(do.ss)
  {
    pwr <- power    
    OBJ <-
    function(x, groups, effect.size, e.I, power)
    {
      n.sample <- exp(x)
      ncp <- n.sample*effect.size^2/2
      df1 <- groups - 1
      df2 <- groups*(n.sample-1)
      Fcrit <- qf(1-e.I, df1=df1, df2=df2)
      ((pwr - (1-pf(Fcrit, ncp=ncp, df1=df1, df2=df2)))^2)^(1/1.25)
    }
    norm.ss <- 2*(qnorm(power) + qnorm(1-e.I))^2/effect.size^2
    l <- log(1/10*norm.ss)
    u <- log(10*norm.ss)
    opt <- optimize(f=OBJ, lower=l, upper=u, groups=groups, effect.size=effect.size, e.I=e.I, power=power)
    n.sample <- round(exp(opt$minimum))
    obj <- opt$objective
  }
  if(do.p)
  {
    obj <- 0
    ncp <- n.sample*effect.size^2/2
    df1 <- groups - 1
    df2 <- groups*(n.sample-1)
    Fcrit <- qf(1-e.I, df1=df1, df2=df2)
    pwr <- 1-pf(Fcrit, ncp=ncp, df1=df1, df2=df2)
  }
  ncp <- n.sample*effect.size^2/2
  df1 <- groups - 1
  df2 <- groups*(n.sample-1)
  Fcrit <- qf(1-e.I, df1=df1, df2=df2)
  out <- list(power=pwr, n.sample=n.sample, Fcrit=Fcrit, ncp=ncp, df1=df1, df2=df2, objective=obj)
  as.data.frame(out)
}

"factorial.design" <-
function(...)
{
    m <- match.call()
    cc <- m
    cc[[1]] <- as.name("c")
    arg.vals <- eval(cc, sys.parent())
    n.args <- length(arg.vals)
    
    n.l <- 0
    n.r <- n.args-1
    fwd.cpd <- c(1, cumprod(arg.vals))
    rev.cpd <- c(cumprod(arg.vals[n.args:2])[(n.args-1):1],1)
    "%,%" <- paste0
    main.call <- as.call(expression(cbind))
    for(k in 1:n.args)
    {
        arg.k <- eval(m[[1+k]], sys.parent())
        rep.call <- as.call(expression(rep))
        rep.call$x <- 1:arg.k
        if(n.l>0) rep.call$times <- fwd.cpd[k]
        if(n.r>0) rep.call$each <- rev.cpd[k]
        main.call[["x" %,% k]] <- rep.call
        n.l <- n.l + 1
        n.r <- n.r - 1
    }
    eval(main.call)
}

"nna" <-
function(x)
{
  par.env.nms <- ls.str(, envir=parent.frame())
  m <- match.call()
  vnm <- as.character(m[[2]])
  exts <- vnm %in% par.env.nms
  n <- length(get(par.env.nms[1], envir=parent.frame()))
  ans <- rep(NA, n)
  if(exts) ans <- get(vnm, envir=parent.frame())
  ans
}

is.formula <- function(x)class(x)=="formula"

"if.na.x" <- function(x, x0=FALSE)
{
  d.x <- dim(x)
  x. <- c(as.matrix(x))
  y <- rep(x0, length(x.))
  y[!is.na(x.)] <- x.[!is.na(x.)]
  structure(y, dim=d.x)
}

if.y.z <-
function(x, y=0,z=1)
{
  ans <- x
  ans[x==y] <- z
  ans
}

if.absInf.x <-
function(x, z)
{
  ans <- x
  ans[abs(x)==Inf] <- z
  ans
}

if.0.rm <- function(x)x[x!=0]


#######################################################################
## argument checking/default argument block                          ##
#######################################################################
## helper fns
"is.int" <- function(x)(floor(x)==x)
"cma.and" <- function(s){n.s <- length(s); c(rep(", ", n.s-2), " and ", "")}
"cma.or" <- function(s){n.s <- length(s); c(rep(", ", n.s-2), " or ", "")}
"list.a" <-
function(lst, pfx=NULL, sfx=NULL)
  pfx %,% paste(c(t(cbind("\'" %,% lst %,% "\'", cma.and(lst)))), collapse="") %,% sfx
"list.o" <-
function(lst, pfx=NULL, sfx=NULL)
    pfx %,% paste(c(t(cbind("\'" %,% lst %,% "\'", cma.or(lst)))), collapse="") %,% sfx

## arg checking defaults
pwrfdrmthd.sfx <- c(`Theoretical`="th", `simulation`="sim")
pwrfdrmthd.vals <- names(pwrfdrmthd.sfx)
FDP.cntl.mth.thrsh.def <- c(0.20, 50)
FDP.mthd.vals <- c("BHFDR","BHFDX","Romano","Auto","both","Holm","Hochberg","Bonferroni")
distopt.vals <- 0:2
type.vals <- c("paired", "balanced", "unbalanced")
ctype.vals <- c("CS-Blocks", "Toeplitz-Blocks")

nera <- c("n.sample", "effect.size", "r.1", "alpha")
Nnera <- c(nera, "N.tests")
dNnera <- c("data", "formula", "alpha")

## default values for missing arguments
default <- list(u=function(...)seq(from=0,to=1,len=100),
                alpha=function(...)NA,
                groups=function(...)2,
                control=function(groups=groups, alpha=alpha, corr.struct=corr.struct)
                        list(tol=1e-08, max.iter=c(1000,20), sim.level=2, FDP.meth.thresh=FDP.cntl.mth.thrsh.def,
                             verb=FALSE, low.power.stop=TRUE, ast.le.a=TRUE, show.footer=FALSE),
                distopt=function(groups=groups, alpha=alpha, corr.struct=corr.struct)
                (0*(groups<=2 && (length(corr.struct$rho)>0)) + 1*(groups<=2 && (length(corr.struct)==0)) + 2*(groups>2)),
                type=function(groups=groups,alpha=alpha, corr.struct=corr.struct)
                    c("paired","balanced")[min(max(floor(eval(groups)),1),2)],
                grpj.per.grp1=function(...)1,
                method=function(...)"Theoretical",
                delta=function(groups=groups, alpha=alpha, corr.struct=corr.struct)alpha)

valid.arg.dict <- list(u=function(x)all(x>=0&x<=1),
                       groups=function(x)is.int(x)&(x>0),
                       effect.size=function(x)(x>0),
                       n.sample=function(x)is.int(x)&(x>0),
                       r.1=function(x)(x>0&x<1),
                       alpha=function(x){if(is.na(x)) x <- -1; (x>0 & x < 1)},
                       N.tests=function(x)is.int(x)&(x>0),
                       average.power=function(x)(x>0.5&x<1),
                       TPX.power=function(x)(x>0.5&x<1),
                       lambda=function(x)(x>0&x<1),
                       corr.struct=function(x, cty.vls=ctype.vals)
                       {
                           good <- all(sort(names(x))[1:3] == c("block.size","rho","type"))
                           if(good)
                           {
                               cty <- try(match.arg(x$type, cty.vls))
                               good <- class(cty)!="try-error"
                               cm <- tp <- fct <- TRUE
                               rho <- x$rho
                               if(x$type=="CS-Blocks") cm <- (length(x$rho)==1)
                               if(x$type=="Toeplitz-Blocks") tp <- (length(x$rho)==(x$block.size-1))
                               bs <- is.int(x$block.size)&(x$block.size>0)
                               rh <- all(rho > 0 & rho < 1)
                               good <- good & bs & rh & cm & tp & fct
                           }
                           good
                       },
                       FDP.control.method=function(x, vals=FDP.mthd.vals)
                       {
                           fdp.mthd <- try(match.arg(x, vals))
                           class(fdp.mthd)!="try-error"
                       },
                       distopt=function(x, vals=distopt.vals)
                       {
                           x %in% vals
                       },
                       method=function(x, vals=pwrfdrmthd.vals)
                       {
                           mthd <- try(match.arg(x, vals))
                           class(mthd)!="try-error"
                       },
                       type=function(x, vals=type.vals)
                       {
                           ty <- try(match.arg(x, vals))
                           class(ty)!="try-error"
                       },
                       grpj.per.grp1=function(x, vals=1)TRUE, # Checking this argument in 'other rules'
                       n.sim=function(x)is.int(x)&(x>0),
                       delta=function(x)(x>0&x<1),
                       tempfile=function(x)is.character(x)&length(x)==1)

"valid.arg.msg" <-
function(x, vrbl)
{
  msg <- 
    switch(vrbl,
         u=ifelse(!x, "All components of 'u' must be between 0 and 1 inclusive\n", NA),
         groups=ifelse(!x, "Argument 'groups' must be integral and > 0\n", NA),
         effect.size=ifelse(!x, "Argument 'effect.size' must be > 0\n", NA),
         n.sample=ifelse(!x, "Argument 'n.sample' must be integral and > 0\n", NA),
         r.1=ifelse(!x, "Argument 'r.1' must be between 0 and 1\n", NA),
         alpha=ifelse(!x, "Argument 'alpha' must be between 0 and 1\n", NA),
         N.tests=ifelse(!x, "Argument 'N.tests' must be integral and > 0\n", NA),
         average.power=ifelse(!x,"Argument 'average.power' must be between 1/2 and 1\n",NA),
         TPX.power=ifelse(!x, "Argument 'TPX.power' must be between 1/2 and 1\n", NA),
         lambda=ifelse(!x, "Argument 'lambda' must be between 1/2 and 1\n", NA),
         corr.struct=ifelse(!x, "Argument 'corr.struct' is mis-specified--consult documentation\n", NA),
         FDP.control.method=ifelse(!x, list.o(FDP.mthd.vals,
                                              "Argument 'FDP.control.method' must be " %,%
                                              "set to one of the values \n"), NA),
         distopt=ifelse(!x, "distopt must be 0, 1, or 2", NA),
         method=ifelse(!x, list.o(pwrfdrmthd.vals, 
                                  "Argument 'method' must be set to one " %,%
                                  "of the values \n"), NA),
         n.sim=ifelse(!x, "Argument 'n.sim' must be integral and > 0\n", NA),
         type=ifelse(!x, "Argument 'type' must be 'paired', 'balanced' or 'unbalanced'\n", NA),
         grpj.per.grp1=ifelse(!x, "this gets checked in other rules", NA),
         delta=ifelse(!x, "Argument 'delta' must be between 0 and 1\n", NA),
         tempfile=ifelse(!x, "Argument 'tempfile' must be a character string\n", NA))
  if(is.na(msg)) msg <- NULL
  msg
}
                       
other.rules <-
function(sppld, frmls, m)
{
  is.msng <- !(frmls %in% sppld)
  names(is.msng) <- frmls
  no.n.sample <- is.msng["n.sample"]
  is.pwr <- (!is.msng["average.power"]) || (!is.msng["TPX.power"])

  pwr.nera.chk <- (!is.pwr && all(!is.msng[nera])) || (is.pwr && (sum(is.msng[nera])==1))
  
  fdp.needs.N <- c("Auto", "BHFDX","Romano")
  
  msg.1 <- msg.2 <- msg.3 <- msg.4 <- msg.5 <- msg.6 <- msg.7 <- msg.8 <- msg.9 <- msg.10 <-
      msg.11 <- msg.12 <- msg.13 <- NULL
  
  ##  pwr.nera.chk failed
  if(!pwr.nera.chk)
  {
      msg.1 <- list.o(nera,"When neither 'average.power' or 'TPX.power' is specified, " %,%
                           "you must specify each of the arguemtns,", "\n")
      msg.1 <- msg.1 %,% list.o(nera, "If 'average.power' or 'TPX.power' is specified, " %,%
                                      "then all but one of arguemtns.", "must be specified")
  }
  ## Specification of 'TPX.power' requires specification of 'lambda'
  if(!is.msng["TPX.power"] && is.msng["lambda"])
    msg.2 <-"Specification of argument 'TPX.power' requires specification of argument 'lambda' "%,%
            "TPX.power = P(TPP > lambda)\n"

  ## 'simulation' method not for sample size calculation
  if(is.pwr && m$method=="simulation")
    msg.3 <- "Sample size calculation not implemented in 'simulation' method.\n"
  
  if(m$method=="simulation" && !is.msng["grpj.per.grp1"])
    msg.4 <- "You specified 'grpj.per.grp1'. The statistic 'type' is unbalanced not " %,%
             "implemented in 'simulation' method.\n"
  
  ## 'simulation' method, or 'FP' with 'Auto', 'BHFDX', or 'Romano' FDP.control.method's need N.tests
  if((m$FDP.control.method %in% fdp.needs.N || m$method=="simulation" || !is.msng["lambda"]
                                            || !is.msng["TPX.power"])&&is.msng["N.tests"])
    msg.5 <- "Argument, 'N.tests' is required for the simulation method, for computing " %,%
             "TPX.power, and for the " %,%
             list.a(fdp.needs.N, "", " FDP.control.method's\n")
  
  ## 'both' FDP.control.method is for simulation only
  if(m$FDP.control.method=="both" && m$method=="Theoretical")
    msg.6 <- "FDP.control.method 'both' is for the simulation method only\n"
  
  ## 'Auto' FDP.control.method is for FixedPoint method only
  if(m$FDP.control.method=="Auto" && m$method=="simulation")
      msg.7 <- "FDP.control.method 'Auto' is for the FixedPoint method only\n"

  ## If 'type' is 'balanced' then argument 'grpj.per.grp1' must be unspecified
  if(m$type=="balanced" && !is.msng["grpj.per.grp1"])
      msg.8 <- "If 'type' is 'balanced' then argument 'grpj.per.grp1' must be unspecified"

  ## If 'type' is 'unbalanced' then argument 'grpj.per.grp1' must be supplied and > 0
  if(m$type=="unbalanced" && (is.msng["grpj.per.grp1"] || any(m$grpj.per.grp1 <= 0)))
      msg.9 <- "If 'type' is 'unbalanced' then argument 'grpj.per.grp1' must be supplied and > 0" 

  ## if 'groups'==1 then 'type' must be unspecified
  if(m$groups==1 && !is.msng["type"])
      msg.10 <- "When argument, 'groups' is 1, then argument 'type' must be unspecified\n"

  ## if 'groups'>=2 and 'type'=='unbalanced' then 'grpj.per.grp1' must be of length groups - 1
  if(m$groups>=2 && m$type=="unbalanced" && length(m$grpj.per.grp1)!=(m$groups-1))
      msg.11 <- "If 'groups'>=2 and 'type'=='unbalanced' then 'grpj.per.grp1' must be of length groups - 1"
     
  ## if 'groups>=3' then 'type' must be balanced or unbalanced
  if(m$groups>=3 && !(m$type %in% c("balanced", "unbalanced")))
      msg.12 <- "If 'groups>=3' then 'type' must be balanced or unbalanced"

  ## if 'FDP.control.method=="Bonferroni" then 'corr.struct' not supported
  if(substring(m$FDP.control.method,1,4)=="Bonf" && !is.null(m$corr.struct)) 
      msg.13 <- "If 'FDP.control.method==\"Bonferroni\" then 'corr.struct' not supported"

  msg <- c(msg.1, msg.2, msg.3, msg.4, msg.5, msg.6, msg.7, msg.8, msg.9, msg.10, msg.11, msg.12, msg.13)
  msg
}

"args.def.err.chk" <-
function(m, sppld, frmls, other.rules=TRUE, eval.env)
{
    if(missing(eval.env))
    {
      calling.env <- sys.parent()
      eval.env <- ifelse(calling.env==0, topenv, sys.parent)
    }
    
    n.sppld <- length(sppld)

    is.msng <- !(frmls %in% sppld)
    names(is.msng) <- frmls
    msng.nms <- names(is.msng)[is.msng]
    
    fill.in.default.arg.nms <- match.arg(msng.nms, names(default), several.ok=TRUE)
    for(v in fill.in.default.arg.nms) m[[v]] <- default[[v]](m$groups, m$alpha, m$corr.struct)

    for(k in 1:n.sppld) m[[sppld[k]]] <- eval(m[[sppld[k]]], eval.env())

    if(!("control" %in% msng.nms) && ("control" %in% frmls))
    {
        pfdr.cntl.call <- as.call(expression(pwrFDR.control))
        pfdr.cntl.call$groups <- m$groups
        pfdr.cntl.call$alpha <- m$alpha
        nms.cntl <- names(m$control)
        n.c <- length(m$control)
        for(k in 1:n.c) pfdr.cntl.call[[nms.cntl[k]]] <- m$control[[nms.cntl[k]]]
        m[["control"]] <- eval(pfdr.cntl.call)
    }
    
    err.msg <- NULL
    chk.valid.arg.nms <- match.arg(sppld, names(valid.arg.dict), several.ok=TRUE)
    for(v in chk.valid.arg.nms)
        err.msg <- c(err.msg, valid.arg.msg(valid.arg.dict[[v]](m[[v]]), v))

    m[["method"]] <- match.arg(m$method, pwrfdrmthd.vals)
    m[["FDP.control.method"]] <- match.arg(m$FDP.control.method, FDP.mthd.vals)
    if(m$method!="simulation") m$n.sim <- NULL

    if(other.rules) err.msg <- c(err.msg, other.rules(sppld, frmls, m))
    
    m$type <- list(`paired`=1, `balanced`=2, `unbalanced`=3)[[m$type]]
    m[["control"]]$is.msng <- is.msng

    attr(m, "err.msg") <- err.msg
    m
}

"pwrFDR.control" <-
function(tol, max.iter, sim.level, FDP.meth.thresh, verb, low.power.stop, ast.le.a, show.footer, groups, alpha)
{
    calling.env <- sys.parent()
    eval.env <- ifelse(calling.env==0, topenv, sys.parent)
    m <- match.call()
    m$groups <- m$alpha <- NULL
    frmls <- names(formals(pwrFDR.control))
    frmls <- frmls[-which(frmls=="groups")]
    frmls <- frmls[-which(frmls=="alpha")]
    sppld <- names(m)[-1]
    n.sppld <- length(sppld)
    is.msng <- !(frmls %in% sppld)
    names(is.msng) <- frmls
    msng.nms <- names(is.msng)[is.msng]
    default.pwrFDR.cntl <- default$control(groups=groups, alpha=alpha)
    fill.in.default.arg.nms <- match.arg(msng.nms, names(default.pwrFDR.cntl), several.ok=TRUE)
    for(v in fill.in.default.arg.nms) m[[v]] <- default.pwrFDR.cntl[[v]]
    for(k in 1:n.sppld) m[[sppld[k]]] <- eval(m[[sppld[k]]], eval.env())
    m[[1]] <- as.name("list")
    eval(m, eval.env())
}
#######################################################################
## end argument checking/default argument block                      ##
#######################################################################



#######################################################################
## test statistic distribution and psi block                         ##
#######################################################################
## the 3rd component reduces to d*(n-1) when nii = n = const
## otherwise it gives effective degrees when n, nii describe unbalanced
## groups.
DF <- c(expression(n.sample-1),
        expression(groups*(n.sample-1)),
        expression((1/n.sample + sum(1/nii.sample))^2/(1/(n.sample^2*(n.sample-1)) + sum(1/(nii.sample^2*(nii.sample-1))))))

NCP <- c(expression(n.sample^0.5*effect.size),
         expression((n.sample/groups)^0.5*effect.size),
         expression(((n.sample-1)/(1 + sum((n.sample-1)/(nii.sample-1))))^0.5*effect.size))

"dists" <-
  as.data.frame(rbind(c(### Normal with 2 groups
  pars0=as.call(expression(c, 0, ncp=0,   sd=1)),
  pars1=as.call(expression(c, 0, ncp=.NCP., sd=1)),
  minv=-Inf, 
  ddist=function(x, par) dnorm(x, mean=par[2], sd=par[3]),
  pdist=function(x, par) pnorm(x, mean=par[2], sd=par[3]),
  qdist=function(x, par) qnorm(x, mean=par[2], sd=par[3]),
  rdist=function(n, par) rnorm(n, mean=par[2], sd=par[3])),
c(### t with 2 groups
  pars0=as.call(expression(c,1,ncp=0,     ndf=.DF.)),
  pars1=as.call(expression(c,1,ncp=.NCP., ndf=.DF.)),
  minv=-Inf, 
  ddist=function(x, par) dt(x, ncp=par[2], df=par[3]),
  pdist=function(x, par) pt(x, ncp=par[2], df=par[3]),
  qdist=function(x, par) qt(x, ncp=par[2], df=par[3]),
  rdist=function(n, par) rt(n, ncp=par[2], df=par[3])),
c(### F with 'groups' groups, effect.size=theta*c(0, 0.5, 0.5, ..., 0.5, 1)
  pars0=as.call(expression(c,2,ncp=0,     ndf1=groups-1, ndf2=.DF.)),
  pars1=as.call(expression(c,2,ncp=.NCP.^2, ndf1=groups-1, ndf2=.DF.)),
  minv=0,
  ddist=function(x, par) df(x, ncp=par[2], df1=par[3], df2=par[4]),
  pdist=function(x, par) pf(x, ncp=par[2], df1=par[3], df2=par[4]),
  qdist=function(x, par) qf(x, ncp=par[2], df1=par[3], df2=par[4]),
  rdist=function(n, par) rf(n, ncp=par[2], df1=par[3], df2=par[4]))))

psi <- function(u, delta) delta*u/(1-(1-delta)*u)
psi.m <- function(j, delta, m) (floor(delta*j) + 1)/(m + floor(delta*j) + 1 - j)

flg.psi.expr <- as.call(expression(`==`, FDP.control.method, "Romano"))
#######################################################################
## END test statistic distribution and psi block                     ##
#######################################################################

.onAttach <- function(libname, pkgname)
{
    options(stringsAsFactors=FALSE)
    ver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname), fields="Version")
    msg <- paste(pkgname, ver) %,% "\n\n" %,% "Type ?pwrFDR"
    msg <- msg %,% "\n\n" %,% "Please cite this package in your work use citation(\"pwrFDR\") \n" %,% 
                              "or toBibtex(citation(\"pwrFDR\"))"
    packageStartupMessage(msg)
}

