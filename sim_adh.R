library(dplyr)
library(ggplot2)
library(parallel)
library(survival)
library(tidyr)

options(mc.cores = detectCores())

RNGkind("L'Ecuyer-CMRG")
set.seed(666)

# Pareto cumulative distribution function
ppareto <- function(q, k, m) ifelse(q < m, 0, 1 - (m / q)^k)

# Number of subjects
n <- 10^3

# Follow-up time
m <- 80

# Scenarios

S <- unlist(recursive = FALSE, mclapply(c(.6, .7, .8), function(l) {
  dis <- function(x) x[-1] - x[-length(x)]
  parT <- c(4, 1)
  pT <- dis(plnorm(0:m, parT[1], parT[2]))
  fn <- function(x) {
    cumsum((1 - cumsum(pT)) * dis(plnorm(0:m, 2.5, x)))[80] - l
  }
  parC <- c(2.5, uniroot(fn, c(0, 10^3))$root)
  pC <- dis(plnorm(0:m, parC[1], parC[2]))
  parCstr <- sprintf("%#.2f", parC)
  parCstr <- paste0("C~logN(", parCstr[1], ",", parCstr[2], ")")
  rowlab <- paste("Censoring =", l)
  mclapply(c(.6, .9), function(p) {
    pI <- rep(p, m)
    collab <- paste("Implementation =", p)
    list(pT = pT, pC = pC, pI = pI, parT = parT, parC = parC,
         parCstr = parCstr, distC = "lnorm", rowlab = rowlab,
         collab = collab)
  })
}))
S <- lapply(1:length(S), function(k) append(S[[k]], list(scenario = k)))

# Figure of the distributions and bias
prob <- do.call(rbind, mclapply(S, function(s) {
  tibble(
    scenario = s$scenario,
    rowlab = s$rowlab,
    collab = s$collab,
    t = 1:m,
    pT = s$pT,
    pC = s$pC,
    sT = 1 - cumsum(pT),
    sC = 1 - cumsum(pC),
    pI = s$pI,
    pA = pI * sT,
    pO = sC + cumsum(pC * cumsum(pT)),
    pM = 1 - pO,
    bias = -pI * sT * cumsum(pC * cumsum(pT)) / pO
  )
}))
fig_dist <- prob %>%
  select(rowlab, collab, t, sT, sC, pM) %>%
  pivot_longer(c(sT, sC, pM), names_to = "x", values_to = "p") %>%
  mutate(
    x = factor(x, c("sT", "sC", "pM"),
               c("P(T>t)", "P(C>t)", "P(C<T,C≤t)"))
  ) %>%
  ggplot(aes(x = t, y = p, colour = x)) +
  geom_line() +
  facet_grid(row = vars(rowlab), col = vars(collab)) +
  labs(x = "Time (t)", y ="Probability") +
  theme_classic() +
  theme(legend.title = element_blank())
fig_bias <- ggplot(prob, aes(x = t, y = bias)) +
  geom_line() +
  facet_grid(row = vars(rowlab), col = vars(collab)) +
  labs(x = "Time (t)", y ="Bias") +
  theme_classic()
if (FALSE) {
  x11(); fig_dist
  x11(); fig_bias
}

# Function that calculates adherence
adherence <- function(time, comp, id, cens) {
  # ------------------------------------------------------------------------
  # Inputs:
  #   time  individual (and time dependent) times
  #   comp  individual and time dependent compliances (0/1)
  #   id    ids in the longitudinal dataset
  #   cens   censoring indicator in the longitudinal dataset
  #
  # Outputs: estimated compliance and adherence as a list of two dataframes
  #   ad: a dataframe containing three estimates of adherence
  #     time  all observed time points
  #     imp   empirical implementation
  #     pers  estimatated persistence
  #     w0    weights for "non-compliant" observation
  #     w1    weights for "compliant" observation
  #     adh1  empirical estimate
  #     adh2  indirect estimate
  #     adh3  weighted estimate
  #   DatLong: a dataframe containing individual and time dependent compliance
  #            data for adherence
  #     id    patients ids
  #     time  observation time
  #     comp  observed compliance
  #     w     associated weight (to correct adherence estimations)
  # ------------------------------------------------------------------------
  # Compute the observed implementation with the implementation function
  imp <- aggregate(comp, list(time), mean, na.rm = TRUE)
  names(imp) <- c("time", "imp")
  imp <- imp[order(imp$time), ]
  # Estimate the persistence with the Kaplan-Meier estimator
  ## 1) Generate the survival data
  last <- aggregate(time, list(id), max)
  names(last) <- c("id", "surv")
  last$surv <- last$surv + 1 # add one unit to get survival time
  ### check that censoring indicator is uniquely defined for each id
  if(any(aggregate(cens, list(id), function(z) length(unique(z)) > 1)[[2]])) {
    stop("censoring indicator is not uniquely defined for each id")
  }
  ### get censoring indicator for each id
  ind <- aggregate(cens, list(id), unique)
  names(ind) <- c("id", "cens")
  ### merge survival time and censoring indicator
  sdat <- merge(last, ind, by = "id")
  ## 2) Estimate the persistence with the Kaplan-Meier estimator
  m <- survfit(Surv(sdat$surv, sdat$cens) ~ 1)
  pers <- data.frame(time = m$time, pers = m$surv)
  ## 3) Complete persistence data
  pers <- do.call(rbind, lapply(imp$time, function(u) {
    if (any(pers$time <= u)) {
      r <- data.frame(
        time = u,
        pers = pers[pers$time == max(pers$time[pers$time <= u]), "pers"]
      )
    } else {
      r <- data.frame(time = u, pers = 1)
    }
  }))
  # Missing values and discontinuations
  # In order to estimate the adherence, we 
  #   i)  add missing values;
  #   ii) add "compliance=0" in the data when a discontinuation is observed.
  ## 1) Put compliances in a dataframe with times and ids
  cdat <- data.frame(id = id, time = time, comp = comp)
  ## 2) Add missing compliances (merge cdat with all ids and observed times
  ##    combinations using the option all = TRUE)
  cdat <- merge(cdat, expand.grid(id = unique(id), time = unique(time)),
                by = c("id", "time"), all = TRUE)
  ## 3) Replace missing values of comp by 0 in case of discontinuation
  for (i in sdat[sdat$cens == 1 & sdat$surv <= max(time), "id"]) {
    cdat[cdat$id == i & cdat$time >= sdat[sdat$id == i, "surv"], "comp"] <- 0
  }
  # Probability to simultaneously be non-/observed and compliant
  # ... will be used to compute the weights
  c0obs <- aggregate(comp ~ time, cdat, function(x) mean(!is.na(x) & x == 0))
  names(c0obs) <- c("time", "c0obs")
  c1obs <- aggregate(comp ~ time, cdat, function(x) mean(!is.na(x) & x == 1))
  names(c1obs) <- c("time", "c1obs")
  # Empirical (biased) adherence
  emp_adh <- aggregate(comp ~ time, cdat, mean, na.rm = TRUE)
  names(emp_adh) <- c("time", "adh1")
  # Put all together
  adh <- merge(imp, pers, by = "time")
  adh <- merge(adh, c0obs, by = "time")
  adh <- merge(adh, c1obs, by = "time")
  adh <- merge(adh, emp_adh, by = "time")
  # Compute the weights
  adh$w0 <- with(adh, c1obs / (imp * pers))
  adh$w1 <- with(adh, ifelse(imp * pers < 1, c0obs / (1 - imp * pers), 1))
  # Compute the corrected adherence
  ## Indirect adherence
  adh$adh2 <- with(adh, imp * pers)
  ## Weighted adherence
  adh$adh3 <- with(adh, w1 * adh1 / (w1 * adh1 + w0 * (1 - adh1)))
  # Format adh
  adh <- adh[order(adh$time),
             c("time", "imp", "pers", "w0", "w1", "adh1", "adh2", "adh3")]
  # Add weights to cdat
  cdat <- merge(cdat, adh[c("time", "w0", "w1")], by = "time", all.x = TRUE)
  cdat$w <- with(cdat, (1 - comp) * w0 + comp * w1)
  cdat <- cdat[order(cdat$id, cdat$time), c("id", "time", "comp", "w")]
  # Return the results
  return(list(ad = adh, DatLong = cdat))
}

# Simulation of the adherence
simAdh <- function(nsim = 10^3, pT, pC, pI, rhoI = .2) {
  m <- min(length(pT), length(pC), length(pI))
  # Conditional distribution of the implementation to generate an autoregessive
  # process.
  # If X and Y are two binary variables, then
  # P(X = i, Y = j) = P(X = i) * P (Y = j) + (-1)^{i+j} * Cov(X, Y), i,j=0,1
  # pI1[t] = P(I_{t+1} = 1 | I_t = 1) =
  #          P(I_{t+1} = 1) + Cov(I_t, I_{t+1}) / P(I_t = 1), t = 1,...,m-1
  # pI0[t] = P(I_{t+1} = 1 | I_t = 0) =
  #          P(I_{t+1} = 1) - Cov(I_t, I_{t+1}) / P(I_t = 0), t = 1,...,m-1
  # Cov[t] = Cov(I_t, I_{t+1}), t = 1,...,m-1
  Cov <- rhoI * sqrt(pI[-m] * (1 - pI[-m]) * pI[-1] * (1 - pI[-1]))
  pI1 <- pI[-1] + Cov / pI[-m]
  pI0 <- pI[-1] - Cov / (1 - pI[-m])
  # Simulations of the adherence
  do.call(rbind, lapply(1:nsim, function(i) {
    # T: Discontinuation time
    # A: Adherence, n*m-matrix (rows=subjects, cols=times)
    T <- ifelse(runif(n) <= pT[1], 1, Inf)
    A <- matrix((T > 1) * (runif(n) <= pI[1]), ncol = 1)
    for (k in 1:(m - 1)) {
      T <- ifelse((T > k) & runif(n) <= pT[k + 1] / (1 - cumsum(pT)[k]),
                  k + 1, T)
      A <- cbind(A, (T > k + 1) * (A[, k] * (runif(n) < pI1[k]) +
                                     (1 - A[, k]) * (runif(n) < pI0[k])))
    }
    # True adherence
    adh <- data.frame(sim = i, time = 1:m, adh.true = apply(A, 2, mean))
    # Censorship times
    C <- sapply(1:n, function(i) 1 + sum(cumsum(pC) <= runif(1)))
    # Data is censored and formatted to analyze it with the package `adherence`
    # Subjects with C=1 are ignored
    dtaC <- do.call(rbind, lapply((1:n)[C>1], function(i) {
      time = 1:(min(T[i], C[i]) - 1)
      data.frame(id = i, time = time, comp = A[i, time],
                 dis = 1 * (T[i] <= C[i]))
    }))
    # Censored adherence - Naive and weighted estimates
    a <- adherence(dtaC$time, dtaC$comp, dtaC$id, dtaC$dis)
    adhC <- with(dtaC, {
      adherence(time, comp, id, dis)$ad[c("time", "adh1", "adh3")]
    })
    names(adhC)[names(adhC) == "adh1"] <- "adh.naive"
    names(adhC)[names(adhC) == "adh3"] <- "adh.weighted"
    adh <- merge(adh, adhC, by = "time")
  }))
}
adh <- do.call(rbind, mclapply(S, function(s) {
  adh <- simAdh(nsim = 10^3, pT = s$pT, pC = s$pC, pI = s$pI)
  cbind(as.data.frame(s[c("scenario", "rowlab", "collab")]), adh)
}))

# Figures of the bias
tmp <- adh %>%
  left_join(select(prob, scenario, time = t, pA),
            by = c("scenario", "time")) %>%
  group_by(rowlab, collab, time) %>%
  summarise(
    bias.mean__naive  = mean(adh.naive - pA),
    bias.sd__naive = sd(adh.naive - pA),
    bias.mean__weighted = mean(adh.weighted - pA),
    bias.sd__weighted = sd(adh.weighted - pA),
    .groups = "drop"
  ) %>%
  pivot_longer(
    starts_with("bias"),
    names_to = c(".value", "method"),
    names_pattern = "(.+)__(.+)"
  )
fig_sim <- ggplot(tmp, aes(x = time, y = bias.mean, colour = method)) +
  geom_point() +
  geom_errorbar(aes(ymin = bias.mean - bias.sd, ymax = bias.mean + bias.sd)) +
  facet_grid(row = vars(rowlab), col = vars(collab)) +
  labs(x = "Time", y ="Bias") +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "bottom")
rm(tmp)
if (FALSE) fig_sim

# Save results
save.image("sim_adh.rda")







