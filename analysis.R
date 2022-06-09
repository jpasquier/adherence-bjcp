# R packages
library(survival)
library(geeM)


# Read data
rm(list=ls())
data <- read.table("blurred_OncoData.txt", header = TRUE)


# ----------------------------- IMPLEMENTATION ------------------------------ #

# Empirical
t <- unique(data$drel)
m <- length(t)
imp <- impF <- impM <- rep(NA,m) 
for (i in 1:m) {
  imp[i] <- mean(data$z[data$drel == t[i]], na.rm = TRUE)
  impF[i] <- mean(data$z[data$drel == t[i] & data$SEX == 0], na.rm = TRUE)
  impM[i] <- mean(data$z[data$drel == t[i] & data$SEX == 1], na.rm = TRUE)
}
Imp <- as.data.frame(cbind(t, imp, impF, impM))

# Predicted (unit=one month)
data$drel1 <- data$drel / 30
out1.gee <- geem(z ~ SEX + drel1 + SEX:drel1, id = id, family = binomial,
                 corstr = "ar1", data = data)
beta <- as.matrix(out1.gee$beta)  
se <- summary(out1.gee)$se.model
t1 <- t / 30
X0 <- cbind(1, 0, t1, 0 * t1)
X1 <- cbind(1, 1, t1, 1 * t1)
imp.hatF <- as.vector(plogis(X0 %*% beta))
imp.hatH <- as.vector(plogis(X1 %*% beta))

colF <- "#00BFC4"
colM <- "#F8766D"

# GRAPH IMPLEMENTATION
plot(impF ~ t, Imp, type = "l", ylim = c(0, 1), xlab = "Time (months)",
     ylab = "", col = colF, lwd = 1, lty = 1, axes = FALSE) 
lines(impM ~ t, Imp, col = colM, lwd = 1, lty = 1) 
lines(imp.hatF ~ t, col = colF, lwd = 2, lty = 1)
lines(imp.hatH ~ t, col = colM, lwd = 2, lty = 1)
axis(1, c(0, 365/4, 365/2, 365/4 * 3, 365), c(0, 3, 6, 9, 12))
axis(2, c(0, 0.2, 0.4, 0.6, 0.8, 1),
     c("0%", "20%", "40%", "60%", "80%", "100%"),
     las = 1)
legend(20, 0.4,
       c("Women empirical implementation", "Women estimated implementation",
         "Men empirical implementation", "Men estimated implementation"),
       col = c(colF, colF, colM, colM), lwd = c(1, 2, 1, 2),
       lty = c(1, 1, 1, 1), bty = "n", y.intersp = 1.5, cex = 0.8)

# ------------------------------- PERSISTENCE ------------------------------- #

id <- unique(data$id)
l <- length(id)
surv <- cens <- sex <- rep(NA,l)
for (i in 1:l) {
  surv[i] <- 1 + length(data$id[data$id == id[i]])
  cens[i] <- data$cens[data$id == id[i]][1]
  sex[i] <- data$SEX[data$id == id[i]][1]
}
dataS <- as.data.frame(cbind(id, surv, cens, sex))
names(dataS) <- c("id", "surv", "cens", "sex")

# -------------------------------- ADHERENCE -------------------------------- #

# Indirect: a = i*p
t <- unique(data$drel)
end <- max(t)

# Women
m1F <- survfit(Surv(surv, cens) ~ 1, dataS[dataS$sex == 0, ])
ti <- c(0, summary(m1F)$time, end)
S <- summary(m1F)$surv
su <- c(1, S, S[length(S)])
persF <- rep(NA,end)
for (i in 1:end) { 
  for (j in 1:(length(ti) - 1)) {
    if (ti[j] <= t[i] & t[i] < ti[j + 1]) {
      persF[i] <- su[j]
     }
  }
}
persF[end] <- S[length(S)]
adhF <- persF * impF

# Men
m1M <- survfit(Surv(surv, cens) ~ 1,dataS[dataS$sex == 1, ])
ti <- c(0, summary(m1M)$time, end)
S <- summary(m1M)$surv
su <- c(1, S, S[length(S)])
persM <- rep(NA, end)
for (i in 1:end) {
  for (j in 1:(length(ti) - 1)) {
    if (ti[j] <= t[i] & t[i] < ti[j + 1]) {
      persM[i] <- su[j]
     }
   }
}
persM[end] <- S[length(S)]
adhM <- persM * impM

# --------------------------------- Weights --------------------------------- #

# Dataset long
id.a <- sort(rep(1:l, end))
drel.a <- rep(1:end, l) 
za <- list()
for (i in 1:l) {
  if (length(data$z[data$id == unique(data$id)[i]]) == end) {
    za[[i]] <- data$z[data$id == id[i]]
  } else {
    za[[i]] <- c(
      data$z[data$id == id[i]],
      rep(100 * (1 - cens[i]), (end - length(data$z[data$id == id[i]])))
    )
  }
} 
z.a <- unlist(za)
z.a[z.a == 100] <- NA
dataA <- as.data.frame(cbind(id.a, drel.a, z.a))
names(dataA) <- c("id", "drel", "z")
dataA$z1 <- dataA$z
dataA$z1[is.na(dataA$z)] <- 2
dataA$SEX <- NA
for (i in 1:nrow(dataA)) {
  dataA$SEX[i] <- data$SEX[data$id==dataA$id[i]][1]
}

# Women
pr0F <- pr1F <- rep(NA, end)
for (i in 1:end) {
  pr0F[i] <- mean(dataA$z1[dataA$drel == i & dataA$SEX == 0] == 0)
  pr1F[i] <- mean(dataA$z1[dataA$drel == i & dataA$SEX == 0] == 1)
}
PrF <-  adhF
w1F <- pr0F / (1 - PrF) 
w1F[is.na(w1F)] <- 1 
w0F <- pr1F / PrF

# Men
pr0H <- pr1H <- rep(NA, end)
for (i in 1:end) {
  pr0H[i] <- mean(dataA$z1[dataA$drel == i & dataA$SEX == 1] == 0)
  pr1H[i] <- mean(dataA$z1[dataA$drel == i & dataA$SEX == 1] == 1)
}
PrH <- adhM
w1H <- pr0H / (1 - PrH) 
w1H[is.na(w1H)] <- 1
w0H <- pr1H / PrH

# merge
ww <- rep(NA, nrow(dataA))
for (i in 1:end) {
  ww[dataA$z == 1 & dataA$drel == i & dataA$SEX == 0] <- w1F[i]
  ww[dataA$z == 1 & dataA$drel == i & dataA$SEX == 1] <- w1H[i] 
  ww[dataA$z == 0 & dataA$drel == i & dataA$SEX == 0] <- w0F[i]
  ww[dataA$z == 0 & dataA$drel == i & dataA$SEX == 1] <- w0H[i]
}
dataA$w <- ww


# weighted adherence = indirect adherence
adhM1 <- adhF1 <- rep(100,end)
for (i in 1:end) {
  adhM1[i] <- weighted.mean(dataA$z[dataA$drel == t[i] & dataA$SEX == 1],
  dataA$w[dataA$drel == t[i] & dataA$SEX == 1], na.rm = TRUE)
  adhF1[i] <- weighted.mean(dataA$z[dataA$drel == t[i] & dataA$SEX == 0],
  dataA$w[dataA$drel == t[i] & dataA$SEX == 0], na.rm = TRUE)
}

# predicted
data1 <- dataA[dataA$z1 < 2,]
data1$drel1 <- data1$drel / 30

out.gee <- geem(z ~ SEX + drel1 + I(drel1^2) + I(drel1^3) + SEX:drel1 + 
                  SEX:I(drel1^2) + SEX:I(drel1^3), id = id, family = binomial,
                corstr = "ar1", weights = w, data = data1) 
beta <- as.matrix(out.gee$beta)
se <- summary(out.gee)$se.model
t1 <- t / 30
X0 <- cbind(1, 0, t1, t1^2, t1^3, 0 * t1, 0 * t1^2, 0 * t1^3)
X1 <- cbind(1, 1, t1, t1^2, t1^3, 1 * t1, 1 * t1^2, 1 * t1^3)
PrF <- as.vector(plogis(X0 %*% beta))
PrH <- as.vector(plogis(X1 %*% beta))


# biased adherence
outb <- geem(z ~ SEX + drel1 + I(drel1^2) + I(drel1^3) + SEX:drel1 +
               SEX:I(drel1^2) + SEX:I(drel1^3), id = id, family = binomial,
             corstr = "ar1", data = data1) 
betab <- as.matrix(outb$beta)
PrFb <- as.vector(plogis(X0%*%betab))
PrHb <- as.vector(plogis(X1%*%betab))

# ----------------------------- graph adherence ----------------------------- #

x11()
plot(survfit(Surv(surv, cens) ~ sex, dataS), lwd = 2, conf.int = FALSE,
     xlab = "Time (months)", ylab = "", mark.time = TRUE, mark = 19, 
     col = c(colF, colM), lty = 2, axes = FALSE)
axis(1, c(0, 365/4, 365/2, 365/4 * 3, 365), c(0, 3, 6, 9, 12))
axis(2, c(0, 0.2, 0.4, 0.6, 0.8, 1),
     c("0%", "20%", "40%", "60%", "80%", "100%"),
     las = 1)
lines(t, adhM, col = colM, lwd = 1, lty = 1)
lines(t, adhF, col = colF, lwd = 1, lty = 1)
lines(PrF ~ t, col = colF, lwd = 2, lty = 1)
lines(PrH ~ t, col = colM, lwd = 2, lty = 1)
lines(PrFb ~ t, col = "#D39200", lwd = 2, lty = 4)
lines(PrHb ~ t, col = "#D39200", lwd = 2, lty = 4)
legend(20, 0.4, 
       c("Women persistence",
         "Women empirical adherence",
         "Women estimated adherence", 
         "Men persistence",
         "Men empirical adherence",
         "Men estimated adherence"),
        col = c(colF, colF, colF, colM, colM, colM),
        lwd = c(2, 1, 2, 2, 1, 2), lty = c(2, 1, 1, 2, 1, 1), 
        bty = "n", y.intersp = 1.5, cex = 0.8)
legend("bottomright", c("biased adherence"), col = "#D39200", lwd = 2,
       lty = 4, bty = "n")
