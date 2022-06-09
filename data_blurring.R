# Set seed
set.seed(readChar("secret_seed.txt", file.info("secret_seed.txt")$size - 1))

# Read original data
data <- read.table("OncoData.txt", header = TRUE)

# Calculate compliances
data$z <- 1 * (data$observed >= data$expected)

# Removes useless columns
data$observed <- NULL
data$expected <- NULL

# For each subject, inverse censorship indicator with a probability of 0.05
cens <- unique(data[c("id", "cens")])
cens$cens <- abs(cens$cens - rbinom(1:nrow(cens), 1, .05))
data <- merge(data[names(data) != "cens"], cens, by = "id")

# For each subject, change the sex with a probability of 0.05
sex <- unique(data[c("id", "SEX")])
sex$SEX <- abs(sex$SEX - rbinom(1:nrow(sex), 1, .05))
data <- merge(data[names(data) != "SEX"], sex, by = "id")

# Implementation of each subject
imp <- aggregate(z ~ id, data, mean, na.rm = TRUE)

# Fill missing compliances
for (i in unique(data$id)) {
  b <- data$id == i & is.na(data$z)
  data[b, "z"] <- rbinom(sum(b), 1, imp[imp$id == i, "z"])
}

# For each subject, the length of the follow-up (before discontinuation or
# censorship) is multiplied by a random value between 0.95 and 1.05.
fu <- aggregate(z ~ id, data, length)
fu$n <- round(fu$z - fu$z * runif(nrow(fu), 0.95, 1.05))
for (i in fu$id) {
  n <- fu[fu$id == i, "n"]
  max_drel <- max(data[data$id == i, "drel"])
  if (n < 0) {
    data <- data[data$id != i | data$drel <= max_drel + n, ]
  } else if (n > 0) {
    add <- data.frame(
      id = i,
      drel = max_drel + 1:n,
      cens = data[data$id == i, "cens"][1],
      z = rbinom(n, 1, imp[imp$id == i, "z"]),
      SEX = data[data$id == i, "cens"][1]
    )
    data <- rbind(data, add)
  }
}
data <- data[order(data$id, data$drel), ]

# Inverse compliance with a probability of 0.025
data$z <- abs(data$z - rbinom(1:nrow(data), 1, .025))

# Save the data
write.table(data, "blurred_OncoData.txt")
