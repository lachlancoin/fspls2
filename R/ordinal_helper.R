## FOR ORDINAL FIT WITH LAMBDA Aradana ortaz

#F(η) = 1 − (λ·e^η + 1)^(−1/λ)
#Solving for η gives:
  # η = log( ((1 − p)^(−λ) − 1) / λ )




pao_ <- function(eta, lambda) {
  if (lambda < 1e-8) return(-expm1(-exp(eta)))          # cloglog limit
  -expm1(-log1p(lambda * exp(eta)) / lambda)
}

# Inverse link: probability -> latent (reuses pao_ from before)
qao_ <- function(p, lambda) {
  if (lambda < 1e-8) return(log(-log1p(-p)))             # cloglog limit
  log(expm1(-lambda * log1p(-p)) / lambda)
}

# Density of the Aranda-Ordaz distribution (needed for the conditional mean)
dao_ <- function(eps, lambda) {
  if (lambda < 1e-8) return(exp(eps - exp(eps)))
  exp(eps - (1 / lambda + 1) * log1p(lambda * exp(eps)))
}

# Closed-form expected information for (theta_1..theta_{K-1}, beta_1..beta_p)
ao_info <- function(theta, beta, X, lambda) {
  X <- as.matrix(X); n <- nrow(X); K <- length(theta) + 1; p <- ncol(X)
  npar <- K - 1 + p
  eta  <- drop(X %*% beta)
  thf  <- c(-Inf, theta, Inf)
  Fd   <- matrix(0, n, K + 1)                       # density at each cutpoint (0 at +-Inf)
  for (j in seq_len(K - 1)) Fd[, j + 1] <- dao_(theta[j] - eta, lambda)
  I <- matrix(0, npar, npar)
  for (j in seq_len(K)) {
    fu <- Fd[, j + 1]; fl <- Fd[, j]
    pj <- pmax(pao_(thf[j + 1] - eta, lambda) - pao_(thf[j] - eta, lambda), 1e-300)
    G  <- matrix(0, n, npar)
    if (j <= K - 1) G[, j]     <- fu
    if (j >= 2)     G[, j - 1] <- -fl
    G[, (K - 1) + seq_len(p)]  <- -(fu - fl) * X
    I <- I + crossprod(G / sqrt(pj))
  }
  I
}

clm_ao <- function(formula, data, lambda,lower=NULL, firth = FALSE, steps = 6, maxit = 1000) {
  stopifnot(is.numeric(lambda), length(lambda) == 1, lambda > 0)
  formula <- as.formula(formula)
  
  # Rows with NAs are dropped here, as clm does
  mf <- model.frame(formula, data)
  resp <- mf[[1]]
  if (!is.factor(resp)) stop("response must be a factor (ideally ordered)")
  resp  <- droplevels(resp)
  mf[[1]] <- resp
  y <- as.integer(resp); K <- nlevels(resp); n <- length(y)
  if (K < 2) stop("response needs at least 2 levels")
  
  # Design matrix without intercept; drop rank-deficient columns like clm does
  full <- model.matrix(formula, mf)
  q    <- qr(full)
  keep <- sort(q$pivot[seq_len(q$rank)])
  if (q$rank < ncol(full))
    warning("dropped rank-deficient column(s): ",
            paste(colnames(full)[-keep], collapse = ", "))
  X <- full[, setdiff(keep, 1L), drop = FALSE]
  p <- ncol(X); npar <- K - 1 + p
  
  lb <- setNames(rep(-Inf, p), colnames(X))
  if (!is.null(lower)) {
    if (length(lower) == 1 && is.null(names(lower))) {
      lb[] <- lower                                   # scalar: applies to every beta
    } else {
      if (is.null(names(lower))) stop("`lower` must be a scalar or a named vector")
      bad <- setdiff(names(lower), colnames(X))
      if (length(bad)) stop("no such coefficient(s): ", paste(bad, collapse = ", "))
      lb[names(lower)] <- lower
    }
  }
  lower_all <- c(rep(-Inf, K - 1), lb)   # thresholds are never bounded
  
  
  # Parameters: first threshold, log-increments (keeps thresholds ordered), betas
  th_of <- function(par) cumsum(c(par[1], exp(par[seq_len(K - 2) + 1])))
  
  nll0_fn <- function(lam) function(par) {          # unpenalised, as before
    th  <- c(-Inf, th_of(par), Inf)
    eta <- drop(X %*% par[(K - 1) + seq_len(p)])
    pr  <- pao_(th[y + 1] - eta, lam) - pao_(th[y] - eta, lam)
    -sum(log(pmax(pr, 1e-300)))
  }
  
  info_fn <- function(par, lam) ao_info(th_of(par), par[(K - 1) + seq_len(p)], X, lam)
  
  # Expected Fisher information on the (theta, beta) scale
  info_fn_old <- function(par, lam) {
    th   <- th_of(par)
    beta <- par[(K - 1) + seq_len(p)]
    eta  <- drop(X %*% beta)
    thf  <- c(-Inf, th, Inf)
    Fd   <- matrix(0, n, K + 1)                      # density at each cutpoint (0 at +-Inf)
    for (j in seq_len(K - 1)) Fd[, j + 1] <- dao_(th[j] - eta, lam)
    I <- matrix(0, npar, npar)
    for (j in seq_len(K)) {                          # category j: (th_{j-1}, th_j]
      fu <- Fd[, j + 1]; fl <- Fd[, j]
      pj <- pmax(pao_(thf[j + 1] - eta, lam) - pao_(thf[j] - eta, lam), 1e-300)
      G  <- matrix(0, n, npar)                       # d p_ij / d(theta, beta)
      if (j <= K - 1) G[, j]     <- fu
      if (j >= 2)     G[, j - 1] <- -fl
      G[, (K - 1) + seq_len(p)]  <- -(fu - fl) * X
      I <- I + crossprod(G / sqrt(pj))
    }
    I
  }
 
  nll_fn <- function(lam) {
    f0 <- nll0_fn(lam)
    if (!firth) return(f0)
    function(par) {
      ld <- determinant(info_fn(par, lam), logarithm = TRUE)
      if (ld$sign <= 0 || !is.finite(ld$modulus)) return(1e20)
      f0(par) - 0.5 * as.numeric(ld$modulus)         # penalised negative log-likelihood
    }
  }
  
  # signature: add ridge = 0
 # nll_pen <- function(lam) {
#    f0 <- nll0_fn(lam)
#    function(par) f0(par) + 0.5 * ridge * sum(par[(K - 1) + seq_len(p)]^2)
#  }
  # then use nll_pen(l) in place of nll_fn(l) inside the nlminb loop
  
  # Start from marginal logit cutpoints and zero betas
  cp <- qlogis(cumsum(tabulate(y, K))[-K] / n)
  #st <- c(cp[1], log(diff(cp)), rep(0, p))
  st <- c(cp[1], log(diff(cp)), pmax(rep(0, p), lb))
  
  # Continuation from logit (lambda = 1) down/up to the target, with warm starts
  path <- if (lambda == 1) 1 else exp(seq(0, log(lambda), length.out = steps + 1))
  path[length(path)] <- lambda
  for (l in path) {
    #  fit <- nlminb(st, nll_fn(l), control = list(iter.max = maxit, eval.max = 2 * maxit))
    fit <- nlminb(st, nll_fn(l), lower = lower_all,
                  control = list(iter.max = maxit, eval.max = 2 * maxit))
    st  <- fit$par
  }
  if (fit$convergence != 0) warning("nlminb did not converge: ", fit$message)
  
  # Covariance on the (theta, beta) scale via the delta method
  at_bound <- is.finite(lower_all) & fit$par <= lower_all + 1e-6
  free <- which(!at_bound)
  H <- optimHess(fit$par, nll_fn(lambda))
  J <- diag(npar)
  J[seq_len(K - 1), 1] <- 1
  if (K > 2) for (j in 2:(K - 1)) for (k in 2:j) J[j, k] <- exp(fit$par[k])
  V <- matrix(NA_real_, npar, npar)
 
  if (firth) {
    Ith <- info_fn(fit$par, lambda)
    V[free, free] <- tryCatch(solve(Ith[free, free, drop = FALSE]),
                              error = function(e) matrix(NA_real_, length(free), length(free)))
  } else {
    # ...your existing optimHess / Jacobian code, unchanged
 
      
      
      Jf <- J[free, free, drop = FALSE]
      V[free, free] <- tryCatch(Jf %*% solve(H[free, free, drop = FALSE]) %*% t(Jf),
                                error = function(e) matrix(NA_real_, length(free), length(free)))
  }
  #H <- optimHess(fit$par, nll_fn(lambda))
  #J <- diag(npar)
  #J[seq_len(K - 1), 1] <- 1
  #if (K > 2) for (j in 2:(K - 1)) for (k in 2:j) J[j, k] <- exp(fit$par[k])
  #V <- tryCatch(J %*% solve(H) %*% t(J),
  #              error = function(e) matrix(NA_real_, npar, npar))
  
  theta <- th_of(fit$par)
  beta  <- fit$par[(K - 1) + seq_len(p)]
  names(theta) <- paste(levels(resp)[-K], levels(resp)[-1], sep = "|")
  names(beta)  <- colnames(X)
  se <- suppressWarnings(sqrt(diag(V)))
  names(se) <- c(names(theta), names(beta))
  
  eta <- drop(X %*% beta)
  cum <- sapply(theta, function(t) pao_(t - eta, lambda))
  if (!is.matrix(cum)) cum <- matrix(cum, nrow = n)
  fitted.prob <- cbind(cum, 1) - cbind(0, cum)
  colnames(fitted.prob) <- levels(resp)
  
  beta1 = rep(0, ncol(data)-1); names(beta1) = colnames(data)[-1]
  beta1[match( names(beta),names(beta1))] = beta
  logLik     = -nll0_fn(lambda)(fit$par)
  attr(logLik,"df")=npar
  structure(list(lambda = lambda, beta = beta1, alpha = theta,
                 se_alpha = se[1:length(theta)], se_beta = se[-(1:length(theta))], vcov = V,
                logLik=logLik, pen.logLik = -fit$objective,
                 npar = npar,firth = firth,
              #   removed = colnames(full)[-keep],
                 nobs = n, convergence = fit$convergence,
                 fitted.prob = fitted.prob, call = match.call()),
            class = "clm_ao")
}

logLik.clm_ao <- function(object, ...)
  structure(object$logLik, df = object$npar, nobs = object$nobs, class = "logLik")

print.clm_ao <- function(x, ...) {
  cat("Aranda-Ordaz ordinal regression, lambda =", x$lambda, "(fixed)\n")
  cat("n =", x$nobs, "  logLik =", format(x$logLik, digits = 8),
      "  AIC =", format(AIC(x), digits = 8), "\n\n")
  est <- c(x$alpha, x$coefficients)
  z   <- est / x$se
  printCoefmat(cbind(Estimate = est, `Std. Error` = x$se, `z value` = z,
                     `Pr(>|z|)` = 2 * pnorm(-abs(z))))
  invisible(x)
}


ao_probs <- function(X, alpha, beta, lambda, type = c("prob", "cum")) {
  type <- match.arg(type)
  X   <- as.matrix(X)
  #print("HH")
  #print(dim(X)); print(dim(beta))
  eta <- drop(X %*% beta)                                   # latent location per row
  cum <- outer(eta, alpha, function(e, t) pao_(t - e, lambda))  # P(Y <= j), n x (K-1)
  if (type == "cum") return(cum)
  pr <- cbind(cum, 1) - cbind(0, cum)                       # P(Y = j), n x K
  colnames(pr) <- c(sub("\\|.*", "", names(alpha)),
                    sub(".*\\|", "", names(alpha)[length(alpha)]))
  pr
}






latent_y <- function(y, lambda, type = c("mean", "mid"), levels = NULL) {
  type <- match.arg(type)
  stopifnot(lambda > 0)
  if (is.null(levels))
    levels <- if (is.factor(y)) levels(y) else sort(unique(y[!is.na(y)]))
  yf <- factor(y, levels = levels)
  K  <- length(levels)
  
  n    <- tabulate(as.integer(yf), K)
  prop <- n / sum(n)
  cum  <- cumsum(prop)
  
  Fl <- c(0, cum[-K])                      # cumulative probability at each category's lower edge
  Fu <- c(cum[-K], 1)                      # ... and upper edge
  lo <- qao_(Fl, lambda)                   # latent cutpoints (-Inf for the first category)
  hi <- qao_(Fu, lambda)                   # (+Inf for the last)
  
  mid <- qao_((Fl + Fu) / 2, lambda)       # probability midpoint of each category
  
  if (type == "mid") {
    latent <- mid
  } else {
    latent <- vapply(seq_len(K), function(j) {
      den <- Fu[j] - Fl[j]
      if (den < 1e-12) return(mid[j])      # empty category: fall back
      num <- tryCatch(integrate(function(e) e * dao_(e, lambda), lo[j], hi[j])$value,
                      error = function(e) NA_real_)
      num / den                            # E[Z | lo < Z <= hi]
    }, numeric(1))
  }
  
  if (any(n == 0)) warning("empty categories: ", paste(levels[n == 0], collapse = ", "))
  z = unname(latent[as.integer(yf)])
  return(z)
#  tab <- data.frame(level = levels, n = n, prop = prop, lower = lo, upper = hi, latent = latent)
#  list(table = tab,z=z, lambda = lambda, type = type)
}

clm_ao_grid <- function(formula, data, lambdas, cores = 1, ...) {
  one <- function(l) {
    tryCatch(clm_ao(formula, data, lambda = l, ...),
             error = function(e) structure(list(lambda = l, error = conditionMessage(e)),
                                           class = "clm_ao_fail"))
  }
  fits <- if (cores > 1) parallel::mclapply(lambdas, one, mc.cores = cores)
  else lapply(lambdas, one)
  names(fits) <- format(lambdas)
  
  ok  <- !vapply(fits, inherits, logical(1), "clm_ao_fail")
  tab <- data.frame(
    lambda = lambdas,
    logLik = NA_real_, pen.logLik = NA_real_, AIC = NA_real_, conv = NA_integer_
  )
  tab$logLik[ok]     <- vapply(fits[ok], function(f) f$logLik, numeric(1))
  tab$pen.logLik[ok] <- vapply(fits[ok], function(f) if (is.null(f$pen.logLik)) NA_real_ else f$pen.logLik, numeric(1))
  tab$AIC[ok]        <- -2 * tab$logLik[ok] + 2 * vapply(fits[ok], function(f) f$npar, numeric(1))
  tab$conv[ok]       <- vapply(fits[ok], function(f) f$convergence, numeric(1))
  list(fits = fits, table = tab)
}


ao_link <- function(lambda, reflected = FALSE) {
  stopifnot(lambda > 0)
  eps <- .Machine$double.eps
  F_  <- function(e) -expm1(-log1p(lambda * exp(e)) / lambda)          # CDF
  Q_  <- function(p) log(expm1(-lambda * log1p(-p)) / lambda)          # inverse
  f_  <- function(e) exp(e - (1 / lambda + 1) * log1p(lambda * exp(e)))# density
  
  if (!reflected) {
    linkfun <- function(mu) Q_(mu)
    linkinv <- function(eta) pmin(pmax(F_(eta), eps), 1 - eps)
    mu.eta  <- function(eta) pmax(f_(eta), eps)
  } else {
    linkfun <- function(mu) -Q_(1 - mu)
    linkinv <- function(eta) pmin(pmax(1 - F_(-eta), eps), 1 - eps)
    mu.eta  <- function(eta) pmax(f_(-eta), eps)
  }
  structure(list(linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta,
                 valideta = function(eta) TRUE,
                 name = paste0("AO(", lambda, if (reflected) ",reflected" else "", ")")),
            class = "link-glm")
}





##note linkinv is for predicting, linkinv(X*Beta);  this aligns with powfunc
pow_link <- function(pow, norm = 1, offset = 1e-4) {
  stopifnot(pow > 0, norm > 0)
  structure(list(
    linkfun  = function(mu)  sign(mu) * abs(mu)^(1 / pow) * norm - offset,
    linkinv  = function(eta) { x1 <- (eta + offset) / norm; sign(x1) * abs(x1)^pow },
    mu.eta   = function(eta) {
      x1 <- (eta + offset) / norm
      pmin(pmax(pow / norm * abs(x1)^(pow - 1), 1e-10), 1e10)
    },
    valideta = function(eta) TRUE,
    name = paste0("pow2(", pow, ")")),
    class = "link-glm")
}

rand_link <- function(seed) {
  stopifnot(seed %% 1 == 0)
  structure(list(
    linkfun  = function(mu)  { set.seed(seed); inds <- sample.int(length(mu));  mu[order(inds)]} ,
    linkinv  = function(eta) {   set.seed(seed); inds <- sample.int(length(eta)); eta[inds] },
    mu.eta   = function(eta) return(NA),
    valideta = function(eta) TRUE,
    name = paste0("rand(", seed, ")")),
    class = "link-glm")
}


perm_link <- function(seed) {
  stopifnot(seed %% 1 == 0)
  structure(list(
    linkfun  = function(mu)  { rotate(mu, -seed)} ,
    linkinv  = function(eta) {   rotate(eta, seed) },
    mu.eta   = function(eta) return(NA),
    valideta = function(eta) TRUE,
    name = paste0("perm(", seed, ")")),
    class = "link-glm")
}




ordinal_validation<-function(){
  
  set.seed(42)
  lambda <- 0.5
  K <- 4; p <- 2; n <- 5000
  theta0 <- c(-1, 0.3, 1.5); beta0 <- c(0.8, -0.5)
  X <- matrix(rnorm(n * p), n, p, dimnames = list(NULL, c("x1", "x2")))
  eta0 <- drop(X %*% beta0)
  
  sim_y <- function() {                       # Y = j  <=>  theta_{j-1} < eta + eps <= theta_j
    z <- eta0 + qao_(runif(n), lambda)
    findInterval(z, theta0) + 1L
  }
  y <- sim_y()
  par0 <- c(theta0, beta0); npar <- length(par0)
  
  prob_fn <- function(par) {                  # n x K category probabilities
    theta <- par[seq_len(K - 1)]; beta <- par[(K - 1) + seq_len(p)]
    eta <- drop(X %*% beta)
    cum <- outer(eta, theta, function(e, t) pao_(t - e, lambda))
    cbind(cum, 1) - cbind(0, cum)
  }
  
  I_closed <- ao_info(theta0, beta0, X, lambda)
  
  ## ---- Test 1 (exact): closed form vs finite-difference derivatives of the probabilities
  P  <- prob_fn(par0)
  dP <- array(0, c(n, K, npar))
  for (k in seq_len(npar)) {
    h <- 1e-5; up <- dn <- par0; up[k] <- up[k] + h; dn[k] <- dn[k] - h
    dP[, , k] <- (prob_fn(up) - prob_fn(dn)) / (2 * h)
  }
  I_fd <- matrix(0, npar, npar)
  for (j in seq_len(K)) {
    Gj <- matrix(dP[, j, ], n, npar)
    I_fd <- I_fd + crossprod(Gj / sqrt(P[, j]))
  }
  cat("Test 1 (closed form vs finite differences): max rel diff =",
      max(abs(I_closed - I_fd)) / max(abs(I_fd)), "  [expect < 1e-6]\n")
  
  ## ---- Test 2 (statistical): E[score score'] = I, and mean score = 0
  idx <- cbind(seq_len(n), y)
  S <- t(sapply(seq_len(n), function(i) dP[i, y[i], ] / P[i, y[i]]))
  cat("Test 2a (outer product of scores vs info): max rel diff =",
      max(abs(crossprod(S) - I_closed)) / max(abs(I_closed)), "  [expect ~ 0.05 or less]\n")
  cat("Test 2b (mean score, should be ~0):", round(colMeans(S), 4), "\n")
  
  ## ---- Test 3 (statistical): observed Hessian at the truth vs expected information
  nll <- function(par) -sum(log(pmax(prob_fn(par)[idx], 1e-300)))
  H <- optimHess(par0, nll)
  cat("Test 3 (observed Hessian vs info): max rel diff =",
      max(abs(H - I_closed)) / max(abs(I_closed)), "  [expect ~ 0.05 or less]\n")
  
  ## ---- Test 4 (statistical): sampling variance of the unpenalised MLE vs diag(solve(I))
  R <- 200
  pred_var <- diag(solve(I_closed))[(K - 1) + seq_len(p)]
  est <- matrix(NA_real_, R, p); se2 <- matrix(NA_real_, R, p)
  bad <- character(0)
  for (r in seq_len(R)) {
    print(r)
    d <- data.frame(y = factor(sim_y(), levels = 1:K, ordered = TRUE), X)
    
    #str(d)                                        # expect columns y, x1, x2
    ##colnames(model.matrix(y ~ x1 + x2, model.frame(y ~ x1 + x2, d)))   # expect "(Intercept)" "x1" "x2"
    #f <- clm_ao(y ~ x1 + x2, d, lambda = 0.5, firth = FALSE, steps = 1)
    #str(f$coefficients); f$alpha
    
    
    f <- tryCatch(clm_ao(y ~ x1 + x2, d, lambda = lambda, firth = FALSE, steps = 1),
                  error = function(e) e)
    if (inherits(f, "error")) { print("bad"); bad <- c(bad, paste(r, "error:", conditionMessage(f))); next }
    if (length(f$beta) != p || length(f$alpha) != K - 1) {
      bad <- c(bad, paste(r, "layout: nbeta =", length(f$beta), "nalpha =", length(f$alpha)))
      next
    }
    est[r, ] <- f$beta
    se2[r, ] <- f$se_beta^2
    
   
  }
  length(bad); head(bad)
  
  cat("Test 4: empirical var / predicted var =",
      round(apply(est, 2, var, na.rm = TRUE) / pred_var, 3), "  [expect 1 +- ~0.2]\n")
  cat("        mean reported SE^2 / predicted =",
      round(colMeans(se2, na.rm = TRUE) / pred_var, 3), "  [expect ~1]\n")
  ##NOW TEST FIRTH
  n_small <- 20; R <- 200
  Xs <- X[1:n_small, ]; es <- drop(Xs %*% beta0)
  one <- function(firth) {
    z <- es + qao_(runif(n_small), lambda)
    d <- data.frame(y = factor(findInterval(z, theta0) + 1L, levels = 1:K, ordered = TRUE), Xs)
    tryCatch(clm_ao(y ~ x1 + x2, d, lambda = lambda, firth = firth, steps = 1)$beta,
             error = function(e) rep(NA_real_, p))
  }
  set.seed(1); ml <- t(replicate(R, one(FALSE)))
  set.seed(1); fi <- t(replicate(R, one(TRUE)))        # same data for both
  rbind(truth = beta0,
        ml_mean = colMeans(ml, na.rm = TRUE), firth_mean = colMeans(fi, na.rm = TRUE),
        ml_sd   = apply(ml, 2, sd, na.rm = TRUE), firth_sd = apply(fi, 2, sd, na.rm = TRUE))
  
}


TEST=FALSE;
if(TEST){
  
  ordinal_validation()
  
  fam <- pow_link(0.1)
  fam$linkinv(2)
  
  fam <- rand_link(10)
  fam$linkinv(1:10)
  
  
  fam <- perm_link(1)
  fam$linkinv(1:10)
  
  fam<-ao_link(0.5, reflected = TRUE)
  plot(-10:10,fam$linkinv(-10:10))

  
  df2 = df1
  df2$y = ifelse(df1$y >median(as.numeric(as.character(df1$y))), 1,0)
  df2$y = factor(df2$y, ordered=TRUE)
  m <- glm(func , data = df2, family = binomial(link = ao_link(0.5, reflected=TRUE)))
  m1<-clm_ao(func, df2,lambda = 0.5, firth = FALSE);
  cbind(m$coefficients, c(m1$alpha,m1$beta))
  c(logLik(m), m1$logLik)
  
m_ao_0 <- clm_ao(func, df1,lower=NULL, lambda = 0.25, firth = FALSE);
m_ao_1 <- clm_ao(func, df1,lower=0, lambda = 0.25, firth = TRUE)

m_ao_0$beta
m_ao_1$beta

m_ao
aop=ao_probs(X, m_ao$alpha, m_ao$coefficients, 0.25, type="prob")

#reverse
z <- latent_y(df1$y, lambda = 0.25, type = "mean")

lambdas = seq(0,1,by =0.01); names(lambdas)=lambdas
p <- seq(0,1,by=0.001); names(p) = p
toplot = .merge1_new(lapply(lambdas, function(lambda){
  df = data.frame(cbind(p,qao_(p, lambda)))
  names(df) = c("p","q")
  df
}), addName="lambda")
toplot$lambda = as.numeric(toplot$lambda)
#toplot$lambda = factor(toplot$lambda, ordered=TRUE)
unique(toplot$lambda)
ggplot(toplot, mapping = aes(x=p, y=q, color=lambda))+geom_point()
for (lam in c(0.1, 0.5, 1, 2)) print(all.equal(pao_(qao_(p, lam), lam), p))

# Covariate-free version: latent cutpoints implied by the marginal proportions
#qao_(cumsum(prop.table(table(df1$y)))[-nlevels(df1$y)], lambda = 0.25)

}

