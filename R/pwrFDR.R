"dists" <-
  as.data.frame(rbind(c( ### Normal with 2 groups
                        pars0=as.call(expression(c,0,ncp=0,                            sd=1)),
                        pars1=as.call(expression(c,0,ncp=(n.sample/2)^0.5*effect.size, sd=1)),
                        minv=-Inf, 
                        ddist=function(x, par) dnorm(x, mean=par[2], sd=par[3]),
                        pdist=function(x, par) pnorm(x, mean=par[2], sd=par[3]),
                        qdist=function(x, par) qnorm(x, mean=par[2], sd=par[3])), 
                      c( ### t with 2 groups
                        pars0=as.call(expression(c,1,ncp=0,                            ndf=2*n.sample - 2)),
                        pars1=as.call(expression(c,1,ncp=(n.sample/2)^0.5*effect.size, ndf=2*n.sample - 2)),
                        minv=-Inf, 
                        ddist=function(x, par) dt(x, ncp=par[2], df=par[3]),
                        pdist=function(x, par) pt(x, ncp=par[2], df=par[3]),
                        qdist=function(x, par) qt(x, ncp=par[2], df=par[3])), 
                      c(### F with 'groups' groups, effect.size=theta*c(0, 0.5, 0.5, ..., 0.5, 1)
                        pars0=as.call(expression(c,2,ncp=0,                        ndf1=groups-1, ndf2=groups*(n.sample-1))),
                        pars1=as.call(expression(c,2,ncp=n.sample*effect.size^2/2, ndf1=groups-1, ndf2=groups*(n.sample-1))),
                        minv=0,
                        ddist=function(x, par) df(x, ncp=par[2], df1=par[3], df2=par[4]),
                        pdist=function(x, par) pf(x, ncp=par[2], df1=par[3], df2=par[4]),
                        qdist=function(x, par) qf(x, ncp=par[2], df1=par[3], df2=par[4]))))

"pwrFDR" <-
function(groups=2, effect.size, n.sample, r.1, FDR, N.tests, average.power, L.power, lambda,
         method=c("approximate","simulation", "JL", "Iz"), control=list(version=0, tol=1e-8, max.iter=1000, distopt=1, CS=list(NULL)),
         n.sim=1000, temp.file)
{
    .call. <- m <- match.call()
    if(missing(groups)) groups <- m$groups <- 2
    if(!missing(groups))
    {
      if(groups > 2) control$distopt <- 2
    }
    if(!missing(effect.size))
    {
      effect.size <- m$effect.size <- eval(m$effect.size, sys.parent())
      if(effect.size <= 0) stop("Argument 'effect.size' must be between positive")
    }
    if(!missing(n.sample))
    {
      n.sample <- m$n.sample <- eval(m$n.sample, sys.parent())
      if(n.sample <= 0) stop("Argument 'n.sample' must be between positive")
    }
    if(!missing(r.1))
    {
      r.1 <- m$r.1 <- eval(m$r.1, sys.parent())
      if(r.1 >=1 || r.1 <= 0) stop("Argument 'r.1' must be between 0 and 1")
    }
    if(!missing(FDR))
    {
      FDR <- m$FDR <- eval(m$FDR, sys.parent())
      if(FDR >=1 || FDR <= 0) stop("Argument 'FDR' must be between 0 and 1")
    }
    if(!missing(N.tests))
    {
      N.tests <- m$N.tests <- eval(m$N.tests, sys.parent())
      if(N.tests <= 0) stop("Argument 'N.tests' must be between positive")
    }
    if(!missing(average.power))
    {
      average.power <- m$average.power <- eval(m$average.power, sys.parent())
      if(average.power >= 1 || average.power <= 0.5) stop("Argument 'average.power' must be between 0.50 and 1")
    }
    if(!missing(L.power))
    {
      L.power <- m$L.power <- eval(m$L.power, sys.parent())
      if(L.power >= 1 || L.power <= 0.5) stop("Argument 'L.power' must be between 0.50 and 1")
    }
    if(!missing(lambda)) lambda <- m$lambda <- eval(m$lambda, sys.parent())
    if(missing(method)) method <-"Iz"
    if(missing(control)) m$control <- list(version=0, tol=1e-8, max.iter=1000, distopt=1, CS=list(NULL))
    if(!missing(control))
    {
      if(!is.list(control)) stop("Argument 'control' must be a list")
      if(is.null(control$version)) control$version <- 0
      if(is.null(control$tol)) control$tol <- 1e-8
      if(is.null(control$max.iter)) control$max.iter <- 1000
      if(is.null(control$distopt)) control$distopt <- 1
      if(is.null(control$CS)) control$CS <- list(NULL)
      m$control <- control
    }
    nch.mthd <- nchar(method)
    apprx <- substring("approximate", 1, nch.mthd)
    smltn <- substring("simulation", 1, nch.mthd)
    if(method!=apprx && method!=smltn && method!="JL" && method!="Iz")
        stop("Argument 'method' must be set to \"approximate\", \"simulation\", " %,%
             "any substing thereof or to \"JL\", or \"Iz\".")

    long.method <- c("approximate","simulation", "JL", "Iz")[1*(method==apprx)+2*(method==smltn)+3*(method=="JL") + 4*(method=="Iz")]
    sfx <- c("apx","sim","JL", "Iz")[1*(method==apprx)+2*(method==smltn)+3*(method=="JL") + 4*(method=="Iz")]

    pwrFDR.fnc.nm <- "pwrFDR." %,% sfx
    no.n.sample <- missing(n.sample)
    no.power <- missing(average.power) && missing(L.power)
    if((no.n.sample && no.power)||(!no.n.sample && !no.power))
        stop("You must specify exactly one of the arguments 'n.sample', 'average.power' or 'L.power'.")
    if(!no.power && long.method=="simulation") stop("Sample size calculation not implemented in 'simulation' method.")

    if(missing(r.1)||missing(FDR)||missing(effect.size)) stop("Arguments 'r.1', 'FDR' and 'effect.size' are required for all methods")
    if(missing(N.tests)&&(long.method=="approximate"||long.method=="simulation"))
      stop("Argument 'N.tests' is required for the 'approximate' and 'simulation' methods")

    m$method <- m$temp.file <- NULL 
    if(long.method!="simulation") m$n.sim <- NULL
#    if(long.method %in% c("JL", "Iz")) m$N.tests <- NULL
    m[[1]] <- as.name(pwrFDR.fnc.nm)
    ans <- suppressWarnings(eval(m))

    ans$call <- .call.
    class(ans) <- "pwr"
    ans
}

"pwrFDR.Iz" <-
function (groups, effect.size, n.sample, average.power, L.power, lambda, r.1, FDR, N.tests, control)
{
    .call. <- m <- match.call()
    do.ss <- !missing(average.power)||!missing(L.power)
    if(!do.ss)
    {
      m[[1]] <- as.name(as.character(m[[1]]) %,% ".1X")
      m$average.power <- m$l.power <- NULL
      out <- eval(m)
    }
    if(do.ss)
    {
      use.L.pwr <- !missing(L.power)
      OBJ <-
      function (x, groups, effect.size, average.power, L.power, use.L.pwr, r.1, FDR, control)
      {
        n <- exp(x)
        if(!use.L.pwr)
        {
          rslt <- pwrFDR.Iz.1X(groups=groups, effect.size=effect.size, n.sample=n, r.1=r.1, FDR=FDR, control=control)
          pwr <- rslt$average.power
          out <- ((average.power - pwr)^2)^(1/1.25)
        }
        
        if(use.L.pwr)
        {
          rslt <- pwrFDR.Iz.1X(groups=groups, effect.size=effect.size, n.sample=n, r.1=r.1, FDR=FDR, N.tests=N.tests, lambda=lambda,
                               control=control)
          pwr <- rslt$L.power
          out <- ((L.power - pwr)^2)^(1/1.25)
        }
        
        attr(out, "detail") <- rslt
        out
      }
      
      idistopt <- control$distopt
      is.pos <- (dists$minv[[1+idistopt]]==0)

      grps <- groups
      r.0 <- 1 - r.1
      f <- FDR
      f.0 <- r.0*f
      if(use.L.pwr) pwr <- L.power
      if(!use.L.pwr) pwr <- average.power
      gma <- r.1*pwr/(1-f.0)
      ftest.ss <- f.power(power=pwr, groups=grps, effect.size=effect.size, e.I=gma*f)
      l <- log(max(0.75*ftest.ss$n.sample, 5))
      u <- log(min(1.5*ftest.ss$n.sample, 10000000))
      ans <- optimize(f=OBJ, lower=l, upper=u, tol=control$tol, groups=groups, r.1=r.1, FDR=FDR, effect.size=effect.size,
                      average.power=average.power, L.power=L.power, use.L.pwr=use.L.pwr, control=control)
      obj <- ans$objective

      n.sample <- ceiling(exp(ans$minimum))
      
      pars0 <- eval(dists$pars0[[1+idistopt]])
      pars1 <- eval(dists$pars1[[1+idistopt]])
      
      r.0 <- 1 - r.1
      f <- FDR
      f.0 <- r.0*f
      gma <- r.1*ans$average.power/(1-f.0)
      c.g <- qdist(1 - gma*f/2^(!is.pos), pars0)
      eIII <- pdist(-c.g, pars1)

      out <- c(n.sample=n.sample, attr(ans$objective, "detail"))
    }
    out
}

"pwrFDR.Iz.1X" <-
function (groups, effect.size, n.sample, r.1, FDR, N.tests, lambda, control) 
{
    .call. <- match.call()
    n <- n.sample
    r.0 <- 1 - r.1
    f <- FDR
    f.0 <- r.0*f
    do.L.power <- !missing(lambda)
    OBJ <-
    function(x, FDR, r.1, pars0, pars1, control)
    {
        r.0 <- 1 - r.1
        f <- FDR
        f.0 <- r.0*f
        avg.pwr <- logitInv(x)
        
        idistopt <- control$distopt
        is.pos <- (dists$minv[[1+idistopt]]==0)
        
        c.g <- qdist(1-r.1*avg.pwr*f/(2^(!is.pos)*(1-f.0)), pars0)
        F1c <- 1 - pdist(c.g, pars1) + (!is.pos)*pdist(-c.g, pars1)
        ((avg.pwr - F1c)^2)^(1/1.25)
    }

    idistopt <- control$distopt
    pars0 <- eval(dists$pars0[[1+idistopt]])
    pars1 <- eval(dists$pars1[[1+idistopt]])
    is.pos <- (dists$minv[[1+idistopt]]==0)

    res <- optimize(f = OBJ, lower = -5, upper = 10, tol=control$tol, FDR = FDR, r.1 = r.1,
                    pars0=pars0, pars1=pars1, control=control)

    average.power <- logitInv(res$minimum)
    obj <- res$objective
    gma <- r.1*average.power/(1-f.0)
    
    c.g <- qdist(1 - gma*f/2^(!is.pos), pars0)
    if(average.power < 0.25)
    {
        m <- .call.
        m[[1]] <- as.name("pwrFDR.apx.1X")
        m$N.tests <- 500
        m$lambda <- m$L.power <- NULL
        avg.pwr <- eval(m)$average.power
        if(average.power < avg.pwr)
        {
            m <- .call.
            m[[1]] <- as.name("pwrFDR.JL.1X")
            m$control <- eval(m$control)
            m$control$version <- 1
            m$lambda <- m$L.power <- NULL
            
            res <- eval(m)
            average.power <- res$average.power
            if(average.power < avg.pwr) average.power <- avg.pwr
            obj <- res$objective
            c.g <- res$c.g
        }
    }
    gma <- r.1*average.power/(1-f.0)
    c.g <- qdist(1 - gma*f/2^(!is.pos), pars0)
    eIII <- pdist(-c.g, pars1)

    ### var.rtm.SoM <-
    ### function(effect.size,n.sample,groups,r.1,FDR,N.tests,control=list(version=0,tol=1e-8,max.iter=1000,
    
    out <- list(average.power = average.power, c.g = c.g, gamma = gma, objective = obj, err.III=eIII)
    out$call <- .call.
    class(out) <- "pwr"
    v <- var.rtm.SoM(out)$var
    out$sigma.rtm.SoM <- v^0.5
    out <- unclass(out)
    out$call <- NULL
    if(do.L.power)
    {
      v.SoM <- v/N.tests
      average.power <- out$average.power
      L.power <- pnorm((average.power - lambda)/v.SoM^0.5)
      L.eq <- average.power - v.SoM^0.5 * qnorm(average.power)
      out$L.power <- L.power
      out$L.eq <- L.eq
    }
    out
}

"pwrFDR.apx" <-
function (groups, effect.size, n.sample, r.1, FDR, N.tests, average.power, control)
{
    .call. <- m <- match.call()
    m$average.power <- NULL
    do.ss <- !missing(average.power)
    
    idistopt <- control$distopt
    is.pos <- (dists$minv[[1+idistopt]]==0)   
    if(!do.ss)
    {
      m[[1]] <- as.name(as.character(m[[1]]) %,% ".1X")
      out <- eval(m)
    }
    if(do.ss)
    {
      OBJ <-
      function (x, groups, r.1, FDR, N.tests, effect.size, control)
      {
        n <- exp(x)
        (average.power - pwrFDR.apx.1X(groups=groups, effect.size=effect.size, n.sample=n, r.1=r.1, FDR=FDR, N.tests=N.tests,
                                       control=control)$average.power)^2
      }
      
      pars0 <- eval(dists$pars0[[1+idistopt]])
      pars1 <- eval(dists$pars1[[1+idistopt]])

      l <- log(3)
      u <- log(300000)
      ans <- optimize(f=OBJ, lower=l, upper=u, tol=control$tol, r.1=r.1, FDR=FDR, N.tests=N.tests, effect.size=effect.size, control=control)
      out <- list(n.sample=ceiling(exp(ans$minimum)), average.power=average.power, objective=ans$objective)
    }
    out
}

"pwrFDR.apx.1X" <-
function(groups, effect.size, n.sample, r.1, FDR, N.tests, control)
{
    # The event { S > s } is equal to the event { P_{1,(s) <= J f/N.g }
    # The approximation replaces J with s,  { P_{1,(s) <= s f/N.g }
    # which makes calculation very simple. The complement of the second
    # event within the first event is  { P_{1,(s) <= J f/N.g } \  { P_{1,(s) <= s f/N.g }
    # The error probability is  P{ P_{1,(s) <= J f/N.g } -  P{ P_{1,(s) <= s f/N.g }
    # Once again the presence of J makes calcuation of this intractible.
    # This is resolved by replacing J with  m.1/(1-f) which stochastically (if
    # not almost surely) bounds it. Thus a tractible error bound is given by
    # P{ P_{1,(s) <= m1*f/(N.g*(1-f)) } -  P{ P_{1,(s) <= s f/N.g }

    .call. <- m <- match.call()

    f <- FDR
    r.0 <- 1-r.1
    f.0 <- r.0*f

    p.m <- dbinom(1:N.tests, N.tests, r.1)
    m.rng <- which(p.m > 1e-10)
    p.m.rng <- dbinom(m.rng, N.tests, r.1)
    m.mean <- round(r.1*N.tests)
    m.min <- min(m.rng)
    s <- 1:m.min

    idistopt <- control$distopt    
    pars0 <- eval(dists$pars0[[1+idistopt]])
    pars1 <- eval(dists$pars1[[1+idistopt]])

    is.pos <- (dists$minv[[1+idistopt]]==0)
    
    CCDF.s <- rep(0, m.mean)
    CCDF.s.g.m   <- t(1-pbeta(pdist(qdist(1-s*f/(2^(!is.pos)*N.tests), pars0), pars1), outer(-s+1, m.rng, FUN="+"), s))
    CCDF.s[s] <- c(p.m.rng%*%CCDF.s.g.m)
    average.power <- sum(p.m.rng%*%(CCDF.s.g.m/m.rng))
    
    CCDF.s.g.m.U <- t(1-pbeta(pdist(qdist(1 - r.1*f/(2^(!is.pos)*(1-f.0)), pars0), pars1), outer(-s+1, m.rng, FUN="+"), s))
    average.power.U <- sum(p.m.rng%*%(CCDF.s.g.m.U/m.rng))

    out <- list(average.power=average.power, err.bdd=average.power.U-average.power, call=.call.)
    class(out) <- "pwr"
    out
}

"pwrFDR.JL" <-
function (groups, effect.size, n.sample, r.1, FDR, N.tests, average.power, lambda, control)
{
    .call. <- m <- match.call()
    m$average.power <- NULL
    do.ss <- !missing(average.power)
    
    idistopt <- control$distopt
    is.pos <- (dists$minv[[1+idistopt]]==0)   
    do.L.power <- !missing(lambda)
    if(do.L.power) stop("L.power not implemented in the 'JL' method")
    m$lambda <- NULL
    if(!do.ss)
    {
      m[[1]] <- as.name(as.character(m[[1]]) %,% ".1X")
      out <- eval(m)
    }
    if(do.ss)
    {
      OBJ <-
      function (x, groups, r.1, FDR, effect.size, average.power, control)
      {
        n <- exp(x)
        out <- pwrFDR.JL.1X(groups=groups, r.1=r.1, FDR=FDR, n.sample=n, effect.size=effect.size, control=control)
        ((average.power - out$average.power)^2)^(1/1.25)
      }
      l <- log(3)
      u <- log(300000)
      ans <- optimize(f=OBJ, lower=l, upper=u, tol=control$tol, groups=groups, r.1=r.1, FDR=FDR, effect.size=effect.size,
                      average.power=average.power, control=control)

      n.sample <- ceiling(exp(ans$minimum))
      
      pars0 <- eval(dists$pars0[[1+idistopt]])
      pars1 <- eval(dists$pars1[[1+idistopt]])

      out <- list(n.sample=n.sample, average.power=average.power, c.g=qdist(average.power, pars1), objective=ans$objective)
    }
    out
}

"pwrFDR.JL.1X" <-
function (groups, effect.size, n.sample, r.1, FDR, N.tests, control) 
{
    idistopt <- control$distopt    
    pars0 <- eval(dists$pars0[[1+idistopt]])
    pars1 <- eval(dists$pars1[[1+idistopt]])
    is.pos <- (dists$minv[[1+idistopt]]==0)
    OBJ <-
    function(x, groups, FDR, r.1, pars0, pars1, control)
    {
      r.0 <- 1 - r.1
      f.0 <- f <- FDR
      if(control$version==1) f.0 <- r.0*f
      
      c.g <- exp(x)
      
      idistopt <- control$distopt    
      is.pos <- (dists$minv[[1+idistopt]]==0)
      F0c <- 2^(!is.pos)*(1-pdist(c.g, pars0))
      .F1c. <- tryCatch(1-pdist(c.g, pars1) + (!is.pos)*pdist(-c.g, pars1), warning=function(w)w, error=function(e)e)
      if(class(.F1c.)[1] != "numeric") .F1c. <- 1 - pnorm(abs(c.g-pars1[2]))
      F1c <- .F1c.
      ans <- ((f.0 - r.0*F0c/(r.0*F0c + r.1*F1c))^2)^(1/1.25)
      ans
    }

    res <- optimize(f = OBJ, lower = -1, upper = 5, tol=control$tol, FDR = FDR, r.1 = r.1, pars0=pars0, pars1=pars1, control=control)

    c.g <- exp(res$minimum)
    eIII <- pdist(-c.g, pars1)
    obj <- res$objective
    average.power <- F1c <- 1-pdist(c.g, pars1) + pdist(-c.g, pars1)
    out <- list(average.power = average.power, c.g = c.g, objective = obj,err.III=eIII)
    out
}

"pwrFDR.sim" <-
function (groups, effect.size, n.sample, r.1, FDR, N.tests, control, lambda, n.sim=1000)
{
    .call. <- match.call()
    nsim <- n.sim
    do.L.power <- !missing(lambda)

    is.CS <- !is.null(control$CS[[1]])
    if(!is.CS)
    {
      rslt <- .C("pwrFDRsim",
                 nsim    = as.integer(nsim),
                 FDR     = as.double(FDR),
                 Ng      = as.integer(N.tests),
                 r1      = as.double(r.1),
                 n       = as.integer(n.sample),
                 theta   = as.double(effect.size),
                 distopt = as.integer(control$distopt),
                 groups  = as.double(groups),
                 M1      = integer(nsim),
                 J       = integer(nsim),
                 S       = integer(nsim),
                 X       = double(N.tests),
                 PACKAGE = "pwrFDR")
    }
    if(is.CS)
    {
      CS <- control$CS    
      bad <- is.null(CS$rho)||is.null(CS$n.WC)
      if(!bad)
      {
        rho <- control$CS$rho
        n.WC <- control$CS$n.WC
        bad <- (rho <= -1 || rho >= 1 || n.WC <= 0)
        if(!bad) bad <- (N.tests %% n.WC) != 0
        if(bad) stop("list argument 'control' must contain a component, 'CS', of type list, " %,%
                     "with components 'rho' and 'n.WC', which must satisfy -1 < rho < 1 " %,%
                     "and 'n.WC' divides 'N.tests' evenly")
        # cat(sprintf("is.CS=%d, rho=%g, n.WC=%d\n",is.CS, rho,n.WC))
      }
      rslt <- .C("pwrFDRsimCS",
                 nsim    = as.integer(nsim),
                 FDR     = as.double(FDR),
                 Ng      = as.integer(N.tests),
                 r1      = as.double(r.1),
                 n       = as.integer(n.sample),
                 theta   = as.double(effect.size),
                 rho     = as.double(rho),
                 n.WC    = as.integer(n.WC),
                 M1      = integer(nsim),
                 J       = integer(nsim),
                 S       = integer(nsim),
                 X       = double(N.tests),
                 PACKAGE = "pwrFDR")     
    }

    ## void pwrFDRsim(int *pnsim, double *pFDR, int *pNg, double *pr1, int *pn, double *ptheta,
    ##	              int *pdistopt, int *pgroups, int *pM, int *pJ, int *pS, double *pX)

    J <- rslt$J
    S <- rslt$S
    M1 <- rslt$M1

    u <- sort(unique(S/M1))
    u.fill.left <- (1:(floor(u[1]*10)-1))/10
    u.fill.right <- (ceiling(max(u)*10):10)/10
    u <- c(u.fill.left, u, u.fill.right)
    n.u <- length(u)
    DX <- function(x)c(x[1], diff(x))
    CCDF.SoM <- NULL
    for(i in 1:n.u)
      CCDF.SoM <- c(CCDF.SoM, mean(S/M1 >= u[i]))

    average.power <- mean(S %over% M1)
    v.SoM.emp <- var(S %over% M1)
    
    out <- list(average.power=average.power, v.SoM.emp = v.SoM.emp, call=.call.)

    if(do.L.power)
    {
      L.power <- mean(S/M1 >= lambda)
      out$L.power <- L.power
    }
       
    dtl <- list(reps=data.frame(M1=M1, J=J, S=S), CCDF=data.frame(u=u, CCDF=CCDF.SoM))
    dtl$X <- rslt$X
    attr(out, "detail") <- dtl
    class(out) <- "pwr"
    out
}

CDF.Pval <-
function(u, groups=2, r.1, effect.size, n.sample, control)
{
  m <- .call. <- match.call()
  u <- eval(m$u, sys.parent())
  r.1 <- eval(m$r.1, sys.parent())
  effect.size <- eval(m$effect.size, sys.parent())
  n.sample <- eval(m$n.sample, sys.parent())
  if(missing(control)) idistopt <- 1
  if(!missing(control)) idistopt <- control$distopt
  pars0 <- eval(dists$pars0[[1+idistopt]])
  pars1 <- eval(dists$pars1[[1+idistopt]])
  
  u <- abs(u)

  is.pos <- (dists$minv[[1+idistopt]]==0)
  c.g <- qdist(1-u/2^(!is.pos), pars0)
  ans <- (1-r.1)*u + r.1*(1-pdist(c.g, pars1) + (!is.pos)*pdist(-c.g, pars1))
  out <- list(u=u, CDF.Pval=ans, call=.call.)
  class(out) <- "vvv"
  out
}

CDF.Pval.HA <-
function (u, groups=2, r.1, effect.size, n.sample, control)
{
  m <- .call. <- match.call()
  u <- eval(m$u, sys.parent())
  r.1 <- eval(m$r.1, sys.parent())
  effect.size <- eval(m$effect.size, sys.parent())
  n.sample <- eval(m$n.sample, sys.parent())
  if(!missing(groups)) groups <- eval(m$groups, sys.parent())
  if(missing(control)) idistopt <- 1
  if(!missing(control)) idistopt <- control$distopt  
  pars0 <- eval(dists$pars0[[1+idistopt]])
  pars1 <- eval(dists$pars1[[1+idistopt]])

  u <- abs(u)

  is.pos <- (dists$minv[[1+idistopt]]==0)
  c.g <- qdist(1-u/2^(!is.pos), pars0)
  ans <- 1-pdist(c.g, pars1) + (!is.pos)*pdist(-c.g, pars1)
  out <- list(u=u, CDF.Pval.HA=ans, call=.call.)
  class(out) <- "vvv"
  out  
}

var.J.o.rtm <-
function(x=NULL, groups=2, effect.size, n.sample, r.1, FDR, N.tests, control)
{
    .call. <- pwrcall <- vcall <- match.call()
    is.x <- !missing(x)
    x.is.pwr <- is.pi <- FALSE
    if(is.x) if(class(x)=="pwr") x.is.pwr <- TRUE
    if(x.is.pwr)
    {
      pwrcall <- xcall <- x$call
      pi.1 <- x$average.power
      c.g <- x$c.g
      is.pi <- TRUE

      if(is.null(xcall$groups))
      {
        groups <- 2
      }
      if(!is.null(xcall$groups))
      {
        groups <- eval(xcall$groups, sys.parent())
        if(groups < 2) stop("Argument 'groups' must be >=2")    
      }
      if(!is.null(xcall$effect.size))
      {
        effect.size <- eval(xcall$effect.size, sys.parent())
        if(effect.size <= 0) stop("Argument 'effect.size' must be between positive")
      }
      if(!is.null(xcall$n.sample))
      {
        n.sample <- eval(xcall$n.sample, sys.parent())
        if(n.sample <= 0) stop("Argument 'n.sample' must be between positive")
      }
      if(!is.null(xcall$r.1))
      {
        r.1 <- eval(xcall$r.1, sys.parent())
        if(r.1 >=1 || r.1 <= 0) stop("Argument 'r.1' must be between 0 and 1")
      }
      if(!is.null(xcall$FDR))
      {
        FDR <- eval(xcall$FDR, sys.parent())
        if(FDR >=1 || FDR <= 0) stop("Argument 'FDR' must be between 0 and 1")
      }
      if(!is.null(xcall$N.tests))
      {
         N.tests <- eval(xcall$N.tests, sys.parent())
         if(N.tests <= 0) stop("Argument 'N.tests' must be between positive")
      }
      if(is.null(xcall$control)) control <- list(version=0,tol=1e-8,max.iter=1000,distopt=1*(groups<=2) + 2*(groups>2),
                                                    CS=list(rho=NULL,n.WC=NULL))
      if(!is.null(xcall$control)) control <- eval(xcall$control, sys.parent())
      if(!is.null(xcall$control))
      {
        if(is.null(control$version)) control$version <- 0
        if(is.null(control$tol)) control$tol <- 1e-8
        if(is.null(control$max.iter)) control$max.iter <- 1000
        if(is.null(control$distopt)) control$distopt <- 1*(groups<=2) + 2*(groups > 2)
        if(is.null(control$CS)) control$CS <- list(rho=NULL, n.WC=NULL)
      }
    }
    else{
     if(missing(groups))
     {
       groups <- 2
     }
     if(!is.null(vcall$groups))
     {
       groups <- eval(vcall$groups, sys.parent())
       if(groups < 2) stop("Argument 'groups' must be >=2")    
     }
     if(!missing(effect.size))
      {
        effect.size <- eval(vcall$effect.size, sys.parent())
        if(effect.size <= 0) stop("Argument 'effect.size' must be between positive")
      }
      if(!missing(n.sample))
      {
        n.sample <- eval(vcall$n.sample, sys.parent())
        if(n.sample <= 0) stop("Argument 'n.sample' must be between positive")
      }
      if(!missing(r.1))
      {
        r.1 <- eval(vcall$r.1, sys.parent())
        if(r.1 >=1 || r.1 <= 0) stop("Argument 'r.1' must be between 0 and 1")
      }
      if(!missing(FDR))
      {
        FDR <- eval(vcall$FDR, sys.parent())
        if(FDR >=1 || FDR <= 0) stop("Argument 'FDR' must be between 0 and 1")
      }
      if(!missing(N.tests))
      {
        N.tests <- eval(vcall$N.tests, sys.parent())
        if(N.tests <= 0) stop("Argument 'N.tests' must be between positive")
      }
      if(missing(control)) control <- list(version=0, tol=1e-8, max.iter=1000, distopt=1, CS=list(rho=NULL,n.WC=NULL))
      if(!missing(control))
      {
        if(is.null(control$version)) control$version <- 0
        if(is.null(control$tol)) control$tol <- 1e-8
        if(is.null(control$max.iter)) control$max.iter <- 1000
        if(is.null(control$distopt)) control$distopt <- 1*(groups<=2) + 2*(groups>2) 
        if(is.null(control$CS)) control$CS <- list(NULL)
      }
    }
    n <- n.sample
    r.0 <- 1 - r.1
    f <- FDR
    f.0 <- r.0*f

    idistopt <- control$distopt
    pars0 <- eval(dists$pars0[[1+idistopt]])
    pars1 <- eval(dists$pars1[[1+idistopt]])
    is.pos <- (dists$minv[[1 + idistopt]] == 0)

    
    pwrcall[[1]] <- as.name("pwrFDR")
    pwrcall$x <- NULL
    if(!is.pi)
    {
      pwr <- eval(pwrcall, sys.parent())   
      pi.1 <- pwr$average.power
      c.g <- pwr$c.g
    }
      
    gma <- r.1*pi.1/(1-f.0)

    G.pr <- (1-r.1) + r.1*(ddist(c.g, pars1) + (!is.pos)*ddist(-c.g, pars1))/(ddist(c.g, pars0) + (!is.pos)*ddist(-c.g, pars0))
    tau2 <- gma*(1-gma)/(1-G.pr*f)^2
    
    ans <- list(var.J.o.rtm = tau2, average.power = pi.1, gamma = gma, c.g = c.g, call = .call.)

    class(ans) <- "vvv"
    ans
}

var.rtm.SoM <-
function(x=NULL, groups=2, effect.size, n.sample, r.1, FDR, N.tests, control)
{
    .call. <- pwrcall <- vcall <- match.call()
    is.x <- !missing(x)
    x.is.pwr <- is.pi <- FALSE
    if(is.x) if(class(x)=="pwr") x.is.pwr <- TRUE
    if(x.is.pwr)
    {
      pwrcall <- xcall <- x$call
      pi.1 <- x$average.power
      c.g <- x$c.g
      is.pi <- TRUE

      if(is.null(xcall$groups))
      {
        groups <- 2
      }
      if(!is.null(xcall$groups))
      {
        groups <- eval(xcall$groups, sys.parent())
        if(groups < 2) stop("Argument 'groups' must be >=2")    
      }
      if(!is.null(xcall$effect.size))
      {
        effect.size <- eval(xcall$effect.size, sys.parent())
        if(effect.size <= 0) stop("Argument 'effect.size' must be between positive")
      }
      if(!is.null(xcall$n.sample))
      {
        n.sample <- eval(xcall$n.sample, sys.parent())
        if(n.sample <= 0) stop("Argument 'n.sample' must be between positive")
      }
      if(!is.null(xcall$r.1))
      {
        r.1 <- eval(xcall$r.1, sys.parent())
        if(r.1 >=1 || r.1 <= 0) stop("Argument 'r.1' must be between 0 and 1")
      }
      if(!is.null(xcall$FDR))
      {
        FDR <- eval(xcall$FDR, sys.parent())
        if(FDR >=1 || FDR <= 0) stop("Argument 'FDR' must be between 0 and 1")
      }
      if(!is.null(xcall$N.tests))
      {
         N.tests <- eval(xcall$N.tests, sys.parent())
         if(N.tests <= 0) stop("Argument 'N.tests' must be between positive")
      }
      if(!is.null(xcall$groups)) groups <- xcall$groups <- eval(pwrcall$groups, sys.parent())
      if(!is.null(xcall$control)) control <- eval(xcall$control, sys.parent())
      if(is.null(xcall$control)) control <- list(version=0,tol=1e-8,max.iter=1000, distopt=1*(groups<=2) + 2*(groups > 2),
                                                     CS=list(rho=NULL,n.WC=NULL))
      if(!is.null(xcall$control))
      {
        if(is.null(control$version)) control$version <- 0
        if(is.null(control$tol)) control$tol <- 1e-8
        if(is.null(control$max.iter)) control$max.iter <- 1000
        if(is.null(control$distopt)) control$distopt <- 1*(groups<=2) + 2*(groups > 2)
        if(is.null(control$CS)) control$CS <- list(NULL)
      }
    }
    else{
      if(missing(groups))
      {
        groups <- 2
      }
      if(!is.null(vcall$groups))
      {
       groups <- eval(vcall$groups)
       if(groups < 2) stop("Argument 'groups' must be >=2")    
      }
      if(!missing(effect.size))
      {
        effect.size <- eval(vcall$effect.size, sys.parent())
        if(effect.size <= 0) stop("Argument 'effect.size' must be between positive")
      }
      if(!missing(n.sample))
      {
        n.sample <- eval(vcall$n.sample, sys.parent())
        if(n.sample <= 0) stop("Argument 'n.sample' must be between positive")
      }
      if(!missing(r.1))
      {
        r.1 <- eval(vcall$r.1, sys.parent())
        if(r.1 >=1 || r.1 <= 0) stop("Argument 'r.1' must be between 0 and 1")
      }
      if(!missing(FDR))
      {
        FDR <- eval(vcall$FDR, sys.parent())
        if(FDR >=1 || FDR <= 0) stop("Argument 'FDR' must be between 0 and 1")
      }
      if(!missing(N.tests))
      {
        N.tests <- eval(vcall$N.tests, sys.parent())
        if(N.tests <= 0) stop("Argument 'N.tests' must be between positive")
      }      
      if(missing(control)) control <- list(version=0,tol=1e-8,max.iter=1000, distopt=1*(groups<=2) + 2*(groups > 2),
                                           CS=list(rho=NULL,n.WC=NULL))
      if(!missing(control))
      {
        if(is.null(control$version)) control$version <- 0
        if(is.null(control$tol)) control$tol <- 1e-8
        if(is.null(control$max.iter)) control$max.iter <- 1000
        if(is.null(control$distopt)) control$distopt <- 1
        if(is.null(control$CS)) control$CS <- list(NULL)
      }
    }

    r.0 <- 1 - r.1
    f <- FDR
    f.0 <- r.0*f

    idistopt <- control$distopt    
    pars0 <- eval(dists$pars0[[1+idistopt]])
    pars1 <- eval(dists$pars1[[1+idistopt]])
    is.pos <- (dists$minv[[1+idistopt]]==0)
    
    pwrcall[[1]] <- as.name("pwrFDR")
    pwrcall$x <- NULL
    if(!is.pi)
    {
      pwr <- eval(pwrcall, sys.parent())   
      pi.1 <- pwr$average.power
      c.g <- pwr$c.g
    }
      
    gma <- r.1*pi.1/(1-f.0)

    G1.pr <- r.1*(ddist(c.g, pars1) + (!is.pos)*ddist(-c.g, pars1))/(ddist(c.g, pars0) + (!is.pos)*ddist(-c.g, pars0))
    G.pr <- (1-r.1) + G1.pr

    tau2 <- gma*(1-gma)/(1-G.pr*f)^2
    v.W1 <- r.1*pi.1 - r.1^2*pi.1^2
    c.W0.W1 <- - r.0*r.1*gma*f*pi.1
    
    v.top <- v.W1 + G1.pr^2*f^2*tau2 + 2*G1.pr*f*(v.W1 + c.W0.W1)/(1-G.pr*f)
    v.bot <- r.0*r.1
    c.top.bot <- r.0*r.1*(pi.1 + (G1.pr*f*gma*f + G1.pr*f*pi.1)/(1-G.pr*f))
    # T/B - t/b = 1/b (T-t) - t/b^2 (B - b)
    # var(T/B) = 1/b^2 var(T) - 2*t/b^3 cov(T, B) + t^2/b^4 var(B)
    #          = 1/b^2 * (var(T) - t/b * cov(T, B) + (t/b)^2*var(B)

    ans <- 1/r.1^2*(v.top - 2*pi.1*c.top.bot + pi.1^2*v.bot)
    ans <- list(var.rtm.SoM = ans, average.power = pi.1, gamma = gma, c.g = c.g, call = .call.)
    
    class(ans) <- "vvv"
    ans
}

var.rtm.ToJ <-
function(x=NULL, groups=2, effect.size, n.sample, r.1, FDR, N.tests, control)
{
    .call. <- pwrcall <- vcall <- match.call()
    is.x <- !missing(x)
    x.is.pwr <- is.pi <- FALSE
    if(is.x) if(class(x)=="pwr") x.is.pwr <- TRUE
    if(x.is.pwr)
    {
      pwrcall <- xcall <- x$call
      pi.1 <- x$average.power
      c.g <- x$c.g
      is.pi <- TRUE

      if(!is.null(xcall$effect.size))
      {
        effect.size <- eval(xcall$effect.size, sys.parent())
        if(effect.size <= 0) stop("Argument 'effect.size' must be between positive")
      }
      if(!is.null(xcall$n.sample))
      {
        n.sample <- eval(xcall$n.sample, sys.parent())
        if(n.sample <= 0) stop("Argument 'n.sample' must be between positive")
      }
      if(!is.null(xcall$r.1))
      {
        r.1 <- eval(xcall$r.1, sys.parent())
        if(r.1 >=1 || r.1 <= 0) stop("Argument 'r.1' must be between 0 and 1")
      }
      if(!is.null(xcall$FDR))
      {
        FDR <- eval(xcall$FDR, sys.parent())
        if(FDR >=1 || FDR <= 0) stop("Argument 'FDR' must be between 0 and 1")
      }
      if(!is.null(xcall$N.tests))
      {
         N.tests <- eval(xcall$N.tests, sys.parent())
         if(N.tests <= 0) stop("Argument 'N.tests' must be between positive")
      }
      if(!is.null(xcall$control)) control <- eval(xcall$control, sys.parent())
      if(is.null(xcall$control)) control <- list(version=0,tol=1e-8,max.iter=1000,distopt=1*(groups<=2) + 2*(groups > 2),
                                                    CS=list(rho=NULL,n.WC=NULL))
      if(!is.null(xcall$control))
      {
        if(is.null(control$version)) control$version <- 0
        if(is.null(control$tol)) control$tol <- 1e-8
        if(is.null(control$max.iter)) control$max.iter <- 1000
        if(is.null(control$distopt)) control$distopt <- 1*(groups<=2) + 2*(groups > 2)
        if(is.null(control$CS)) control$CS <- list(NULL)
      }
    }
    else{
      if(missing(groups))
      {
        groups <- 2
      }
      if(!is.null(vcall$groups))
      {
       groups <- eval(vcall$groups)
       if(groups < 2) stop("Argument 'groups' must be >=2")    
      }
      if(!missing(effect.size))
      {
        effect.size <- eval(vcall$effect.size, sys.parent())
        if(effect.size <= 0) stop("Argument 'effect.size' must be between positive")
      }
      if(!missing(n.sample))
      {
        n.sample <- eval(vcall$n.sample, sys.parent())
        if(n.sample <= 0) stop("Argument 'n.sample' must be between positive")
      }
      if(!missing(r.1))
      {
        r.1 <- eval(vcall$r.1, sys.parent())
        if(r.1 >=1 || r.1 <= 0) stop("Argument 'r.1' must be between 0 and 1")
      }
      if(!missing(FDR))
      {
        FDR <- eval(vcall$FDR, sys.parent())
        if(FDR >=1 || FDR <= 0) stop("Argument 'FDR' must be between 0 and 1")
      }
      if(!missing(N.tests))
      {
        N.tests <- eval(vcall$N.tests, sys.parent())
        if(N.tests <= 0) stop("Argument 'N.tests' must be between positive")
      }      
      if(missing(control)) control <- list(version=0, tol=1e-8, max.iter=1000, distopt=1, CS=list(rho=NULL,n.WC=NULL))
      if(!missing(control))
      {
        if(is.null(control$version)) control$version <- 0
        if(is.null(control$tol)) control$tol <- 1e-8
        if(is.null(control$max.iter)) control$max.iter <- 1000
        if(is.null(control$distopt)) control$distopt <- 1*(groups<=2) + 2*(groups > 2)
        if(is.null(control$CS)) control$CS <- list(NULL)
      }
    }

    r.0 <- 1 - r.1
    f <- FDR
    f.0 <- r.0*f

    idistopt <- control$distopt    
    pars0 <- eval(dists$pars0[[1+idistopt]])
    pars1 <- eval(dists$pars1[[1+idistopt]])
    is.pos <- (dists$minv[[1+idistopt]]==0)

    pwrcall[[1]] <- as.name("pwrFDR")
    pwrcall$x <- NULL
    if(!is.pi)
    {
      pwr <- eval(pwrcall, sys.parent())
      pi.1 <- pwr$average.power
      c.g <- pwr$c.g
    }
      
    gma <- r.1*pi.1/(1-f.0)

    ans <- (1-r.1)*f*(1-(1-r.1)*f*gma)/gma

    ans <- list(var.rtm.ToJ = ans, average.power = pi.1, gamma = gma, c.g = c.g, call = .call.)

    class(ans) <- "vvv"
    ans
}

# idea.
#    P( T_m / J_m  > \lambda )
#  = P( m^0.5 ( T_m / J_m - (1-r)f )/v^0.5  > m^0.5(\lambda - (1-r)f)/(v/m)^0.5
#  = 1 - pnorm(m^0.5(\lambda - (1-r)f)/v^0.5)
#
# if you want to control the tail probability instead of the expected value, set the above equal to (1-r) f
# and invert for \lambda, obtaining


# (1-r)f =  1 - pnorm(m^0.5(\lambda - (1-r)f)/v^0.5)
# pnorm(m^0.5(\lambda - (1-r)f)/v^0.5) = 1-(1-r)f
# m^0.5(\lambda - (1-r)f)/v^0.5 = qnorm(1-(1-r)f)
# \lambda_{r,f} = (1-r)f + (v/m)^0.5 qnorm(1-(1-r)f)

# so BH-FDR procedure at FDR=f bounds the probability that T_m / J_m  is in excess of \lambda_{r, f}
# by (1-r) f. This is good to have the probability bounded instead of the expected value, but what
# if \lambda_{r,f} is unacceptably large. For f < 1/2, \lambda_{r,f} > f 

# how about finding f* so that for v=v(f*)

# (1-r)f = \lambda_{r,f*} = (1-r)f* + (v(f*)/m)^0.5 qnorm(1-(1-r)f*)


# algorithm
# 1. consider r, theta=effect.size, n.sample, known
# 2. stipulate f
# 3. find f* such that:  (1-r)f = (1-r)f* + (v(f*)/m)^0.5 qnorm(1-(1-r)f*)

"find.f.star" <-
function(groups=2, FDR, r.1, N.tests, effect.size, n.sample)
{
  .call. <- match.call()
  OBJ <-
  function(x, groups, FDR, r.1, N.tests, effect.size, n.sample)
  {
    f.star <- logitInv(x)
    pwr <- pwrFDR(groups=groups, effect.size=effect.size, n.sample=n.sample, r.1=r.1, FDR=f.star)
    v <- var.rtm.ToJ(pwr)$var
    ( ( (1-r.1)*FDR - ((1-r.1)*f.star + (v/N.tests)^0.5 * qnorm(1 -(1-r.1)*f.star))  )^2 )/1.25
  }
  rslt <- 
  optimize(f=OBJ, lower=-8, upper=8, groups=groups, FDR=FDR, r.1=r.1, N.tests=N.tests, effect.size=effect.size, n.sample=n.sample)
  f.star <- logitInv(rslt$min)
  pwr <- pwrFDR(groups=groups, effect.size=effect.size, n.sample=n.sample, r.1=r.1, FDR=f.star)
  v <- var.rtm.ToJ(pwr)$var
  lambda.star <- ((1-r.1)*f.star + (v/N.tests)^0.5 * qnorm(1 -(1-r.1)*f.star))
  prob.star <- 1-pnorm(N.tests^0.5*(lambda.star - (1-r.1)*f.star)/v^0.5)
  obj <- rslt$objective
  out <- as.data.frame(list(f.star=f.star, obj=obj, L.star=lambda.star, P.star=prob.star))
  out <- c(out, pwr)
  out$call <- .call.
  class(out) <- "vvv"
  out
}

"cCDF.SoM" <-
function(lambda, x=NULL, groups, effect.size, n.sample, r.1, FDR, N.tests, control)
{
  .call. <- pwrcall <- vcall <- match.call()
  is.x <- !missing(x)
  x.is.pwr <- is.pi <- FALSE
  if(is.x) if(class(x)=="pwr") x.is.pwr <- TRUE
  if(x.is.pwr)
  {
    pwrcall <- xcall <- x$call
    pwrrslt <- x
    pi.1 <- pwrrslt$average.power
    sgma <- pwrrslt$sigma.rtm.SoM
    c.g <- pwrrslt$c.g
    is.pi <- TRUE
    if(is.null(xcall$groups))
    {
      groups <- 2
    }
    if(!is.null(xcall$groups))
    {
      groups <- eval(xcall$groups, sys.parent())
      if(groups < 2) stop("Argument 'groups' must be >=2")    
    }
    if(!is.null(xcall$effect.size))
    {
      effect.size <- eval(xcall$effect.size, sys.parent())
      if(effect.size <= 0) stop("Argument 'effect.size' must be between positive")
    }
    if(!is.null(xcall$n.sample))
    {
      n.sample <- eval(xcall$n.sample, sys.parent())
      if(n.sample <= 0) stop("Argument 'n.sample' must be between positive")
    }
    if(!is.null(xcall$r.1))
    {
      r.1 <- eval(xcall$r.1, sys.parent())
      if(r.1 >=1 || r.1 <= 0) stop("Argument 'r.1' must be between 0 and 1")
    }
    if(!is.null(xcall$FDR))
    {
      FDR <- eval(xcall$FDR, sys.parent())
      if(FDR >=1 || FDR <= 0) stop("Argument 'FDR' must be between 0 and 1")
    }
    if(!is.null(xcall$N.tests))
    {
       N.tests <- eval(xcall$N.tests, sys.parent())
       if(N.tests <= 0) stop("Argument 'N.tests' must be between positive")
    }
    if(is.null(xcall$N.tests) && missing(N.tests)) 
      stop("Argument 'N.tests' must be specified in either this function call " %,%
           "or in the call producing the supplied argument 'x'")
    if(!is.null(xcall$control)) control <- eval(xcall$control, sys.parent())
    if(is.null(xcall$control)) control <- list(version=0,tol=1e-8,max.iter=1000,distopt=1*(groups<=2) + 2*(groups > 2),
                                               CS=list(rho=NULL,n.WC=NULL))
    if(!is.null(xcall$control))
    {
      if(is.null(control$version)) control$version <- 0
      if(is.null(control$tol)) control$tol <- 1e-8
      if(is.null(control$max.iter)) control$max.iter <- 1000
      if(is.null(control$distopt)) control$distopt <- 1*(groups<=2) + 2*(groups > 2)
      if(is.null(control$CS)) control$CS <- list(NULL)
    }
  }
  else{
   if(missing(groups))
   {
     groups <- 2
   }
   if(!is.null(vcall$groups))
   {
     groups <- eval(vcall$groups, sys.parent())
     if(groups < 2) stop("Argument 'groups' must be >=2")    
   }
   if(!missing(effect.size))
    {
      effect.size <- eval(vcall$effect.size, sys.parent())
      if(effect.size <= 0) stop("Argument 'effect.size' must be between positive")
    }
    if(!missing(n.sample))
    {
      n.sample <- eval(vcall$n.sample, sys.parent())
      if(n.sample <= 0) stop("Argument 'n.sample' must be between positive")
    }
    if(!missing(r.1))
    {
      r.1 <- eval(vcall$r.1, sys.parent())
      if(r.1 >=1 || r.1 <= 0) stop("Argument 'r.1' must be between 0 and 1")
    }
    if(!missing(FDR))
    {
      FDR <- eval(vcall$FDR, sys.parent())
      if(FDR >=1 || FDR <= 0) stop("Argument 'FDR' must be between 0 and 1")
    }
    if(!missing(N.tests))
    {
      N.tests <- eval(vcall$N.tests, sys.parent())
      if(N.tests <= 0) stop("Argument 'N.tests' must be between positive")
    }
    if(missing(control)) control <- list(version=0, tol=1e-8, max.iter=1000, distopt=1, CS=list(rho=NULL,n.WC=NULL))
    if(!missing(control))
    {
      if(is.null(control$version)) control$version <- 0
      if(is.null(control$tol)) control$tol <- 1e-8
      if(is.null(control$max.iter)) control$max.iter <- 1000
      if(is.null(control$distopt)) control$distopt <- 1*(groups<=2) + 2*(groups > 2)
      if(is.null(control$CS)) control$CS <- list(NULL)
    }
  }
  n <- n.sample
  r.0 <- 1 - r.1
  f <- FDR
  f.0 <- r.0*f

  idistopt <- control$distopt
  pars0 <- eval(dists$pars0[[1+idistopt]])
  pars1 <- eval(dists$pars1[[1+idistopt]])
  is.pos <- (dists$minv[[1 + idistopt]] == 0)
  
  pwrcall[[1]] <- as.name("pwrFDR")
  pwrcall$x <- NULL
  if(!is.pi)
  {
    pwrrslt <- eval(pwrcall, sys.parent())   
    pi.1 <- pwrrslt$average.power
    c.g <- pwrrslt$c.g
    sgma <- pwrrslt$sigma.rtm.SoM
  }
  m <- N.tests
  ans <- 1-pnorm(m^0.5*(lambda - pi.1)/sgma)
  pwrrslt$call <- NULL 
  out <- c(cCDF.SoM=ans, unclass(pwrrslt))
  out$call <- .call.

  class(out) <- "vvv"
  out
}

"cCDF.ToJ" <-
function(lambda, x=NULL, groups, effect.size, n.sample, r.1, FDR, N.tests, control)
{
  .call. <- pwrcall <- vcall <- match.call()
  is.x <- !missing(x)
  x.is.pwr <- is.pi <- FALSE
  if(is.x) if(class(x)=="pwr") x.is.pwr <- TRUE
  if(x.is.pwr)
  {
    pwrcall <- xcall <- x$call
    pwrrslt <- x
    pi.1 <- pwrrslt$average.power
    sgma <- pwrrslt$sigma.rtm.SoM
    c.g <- pwrrslt$c.g
    is.pi <- TRUE
    if(is.null(xcall$groups))
    {
      groups <- 2
    }
    if(!is.null(xcall$groups))
    {
      groups <- eval(xcall$groups, sys.parent())
      if(groups < 2) stop("Argument 'groups' must be >=2")    
    }
    if(!is.null(xcall$effect.size))
    {
      effect.size <- eval(xcall$effect.size, sys.parent())
      if(effect.size <= 0) stop("Argument 'effect.size' must be between positive")
    }
    if(!is.null(xcall$n.sample))
    {
      n.sample <- eval(xcall$n.sample, sys.parent())
      if(n.sample <= 0) stop("Argument 'n.sample' must be between positive")
    }
    if(!is.null(xcall$r.1))
    {
      r.1 <- eval(xcall$r.1, sys.parent())
      if(r.1 >=1 || r.1 <= 0) stop("Argument 'r.1' must be between 0 and 1")
    }
    if(!is.null(xcall$FDR))
    {
      FDR <- eval(xcall$FDR, sys.parent())
      if(FDR >=1 || FDR <= 0) stop("Argument 'FDR' must be between 0 and 1")
    }
    if(!is.null(xcall$N.tests))
    {
       N.tests <- eval(xcall$N.tests, sys.parent())
       if(N.tests <= 0) stop("Argument 'N.tests' must be between positive")
    }
    if(is.null(xcall$N.tests) && missing(N.tests)) 
      stop("Argument 'N.tests' must be specified in either this function call " %,%
           "or in the call producing the supplied argument 'x'")
    if(!is.null(xcall$control)) control <- eval(xcall$control, sys.parent())
    if(is.null(xcall$control)) control <- list(version=0,tol=1e-8,max.iter=1000,distopt=1*(groups<=2) + 2*(groups > 2),
                                               CS=list(rho=NULL,n.WC=NULL))
    if(!is.null(xcall$control))
    {
      if(is.null(control$version)) control$version <- 0
      if(is.null(control$tol)) control$tol <- 1e-8
      if(is.null(control$max.iter)) control$max.iter <- 1000
      if(is.null(control$distopt)) control$distopt <- 1*(groups<=2) + 2*(groups > 2)
      if(is.null(control$CS)) control$CS <- list(NULL)
    }
  }
  else{
   if(missing(groups))
   {
     groups <- 2
   }
   if(!is.null(vcall$groups))
   {
     groups <- eval(vcall$groups, sys.parent())
     if(groups < 2) stop("Argument 'groups' must be >=2")    
   }
   if(!missing(effect.size))
    {
      effect.size <- eval(vcall$effect.size, sys.parent())
      if(effect.size <= 0) stop("Argument 'effect.size' must be between positive")
    }
    if(!missing(n.sample))
    {
      n.sample <- eval(vcall$n.sample, sys.parent())
      if(n.sample <= 0) stop("Argument 'n.sample' must be between positive")
    }
    if(!missing(r.1))
    {
      r.1 <- eval(vcall$r.1, sys.parent())
      if(r.1 >=1 || r.1 <= 0) stop("Argument 'r.1' must be between 0 and 1")
    }
    if(!missing(FDR))
    {
      FDR <- eval(vcall$FDR, sys.parent())
      if(FDR >=1 || FDR <= 0) stop("Argument 'FDR' must be between 0 and 1")
    }
    if(!missing(N.tests))
    {
      N.tests <- eval(vcall$N.tests, sys.parent())
      if(N.tests <= 0) stop("Argument 'N.tests' must be between positive")
    }
    if(missing(control)) control <- list(version=0, tol=1e-8, max.iter=1000, distopt=1, CS=list(rho=NULL,n.WC=NULL))
    if(!missing(control))
    {
      if(is.null(control$version)) control$version <- 0
      if(is.null(control$tol)) control$tol <- 1e-8
      if(is.null(control$max.iter)) control$max.iter <- 1000
      if(is.null(control$distopt)) control$distopt <- 1*(groups<=2) + 2*(groups > 2)
      if(is.null(control$CS)) control$CS <- list(NULL)
    }
  }
  n <- n.sample
  r.0 <- 1 - r.1
  f <- FDR
  f.0 <- r.0*f

  idistopt <- control$distopt
  pars0 <- eval(dists$pars0[[1+idistopt]])
  pars1 <- eval(dists$pars1[[1+idistopt]])
  is.pos <- (dists$minv[[1 + idistopt]] == 0)
  
  pwrcall[[1]] <- as.name("pwrFDR")
  pwrcall$x <- NULL
  if(!is.pi)
  {
    pwrrslt <- eval(pwrcall, sys.parent())   
    pi.1 <- pwrrslt$average.power
    c.g <- pwrrslt$c.g
    sgma <- var.J.o.rtm(pwrrslt, N.tests=N.tests)$var^0.5
  }
  m <- N.tests

  ans <- 1-pnorm(m^0.5*(lambda - (1-r.1)*FDR)/sgma)

  pwrrslt$call <- NULL 
  pwrrslt <- unclass(pwrrslt)
  pwrrslt$sigma.rtm.SoM <- pwrrslt$L.power <- pwrrslt$L.eq <- NULL
  out <- c(cCDF.ToJ=ans, unclass(pwrrslt))
  out$sigma.rtm.ToJ <- sgma
  outcall <- .call.

  class(out) <- "vvv"
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

"print.pwr" <- 
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

"is.int" <-
function(x)
{
  abs(x - floor(x)) < 1e-10
}

`+.pwr` <-
  function(x,y)
  {
    xpwr <- x
    ypwr <- y
    if(is(x,"pwr")) xpwr <- ifelse(!is.null(x$call$lambda), x$L.power, x$average.power)
    if(is(y,"pwr")) ypwr <- ifelse(!is.null(y$call$lambda), y$L.power, y$average.power)
    xpwr + ypwr
  }

`-.pwr` <-
  function(x,y)
  {
    xpwr <- x
    ypwr <- y
    if(is(x,"pwr")) xpwr <- ifelse(!is.null(x$call$lambda), x$L.power, x$average.power)
    if(is(y,"pwr")) ypwr <- ifelse(!is.null(y$call$lambda), y$L.power, y$average.power)
    xpwr - ypwr
  }

`*.pwr` <-
  function(x,y)
  {
    xpwr <- x
    ypwr <- y
    if(is(x,"pwr")) xpwr <- ifelse(!is.null(x$call$lambda), x$L.power, x$average.power)
    if(is(y,"pwr")) ypwr <- ifelse(!is.null(y$call$lambda), y$L.power, y$average.power)
    xpwr * ypwr
  }

`/.pwr` <-
  function(x,y)
  {
    xpwr <- x
    ypwr <- y
    if(is(x,"pwr")) xpwr <- ifelse(!is.null(x$call$lambda), x$L.power, x$average.power)
    if(is(y,"pwr")) ypwr <- ifelse(!is.null(y$call$lambda), y$L.power, y$average.power)
    xpwr / ypwr
  }

`^.pwr` <-
  function(x,y)
  {
    xpwr <- x
    ypwr <- y
    if(is(x,"pwr")) xpwr <- ifelse(!is.null(x$call$lambda), x$L.power, x$average.power)
    if(is(y,"pwr")) ypwr <- ifelse(!is.null(y$call$lambda), y$L.power, y$average.power)
    xpwr ^ ypwr
  }

exp.pwr <- function(x) exp(ifelse(!is.null(x$call$lambda), x$L.power, x$average.power))
log.pwr <- function(x, base)log(ifelse(!is.null(x$call$lambda), x$L.power, x$average.power), base)
logitInv.default <- binomial("logit")$linkinv
logit.default <- binomial("logit")$linkfun
logit <- function(mu)UseMethod("logit")
logitInv <- function(eta)UseMethod("logitInv")
logit.pwr <- function(mu)logit.default(ifelse(!is.null(mu$call$lambda), mu$L.power, mu$average.power))
logitInv.pwr <- function(eta)logitInv.default(ifelse(!is.null(eta$call$lambda), eta$L.power, eta$average.power))
"%over%" <-
function(x,y){
  ans <- 0*x
  ans[y!=0] <- x[y!=0]/y[y!=0]
  ans
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

.onAttach <- function(libname, pkgname)
{
    options(stringsAsFactors=FALSE)
    ver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                    fields="Version")
    msg <- paste(pkgname, ver) %,% "\n\n" %,%
           "Type ?pwrFDR"
    packageStartupMessage(msg)
}

