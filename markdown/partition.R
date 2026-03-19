library(fields)
s <- cbind(c(-1, 1), c(0, 0))

reps <- 1000000
prob <- 0
plot(
  s,
  main = "Nearest-knot overlap check (r = 1, 1 knot per quadrant)",
  xlab = "x",
  ylab = "y"
)
abline(h = 0, v = 0)
legend(
  "topright",
  legend = c("Sites s", "Knots when same nearest knot"),
  pch = c(1, 19),
  col = c("black", "dodgerblue"),
  bty = "n"
)

for (i in 1:reps) {
  # get a point in every quadrant of the circle
  p1 <- p2 <- p3 <- p4 <- rep(0, 2)
  
  p1[1] <- runif(1, 0, 1)
  p1[2] <- runif(1, 0, sqrt(1 - p1[1]^2))
  p2[1] <- runif(1, -1, 0)
  p2[2] <- runif(1, 0, sqrt(1 - p2[1]^2))
  p3[1] <- runif(1, -1, 0)
  p3[2] <- runif(1, -sqrt(1 - p3[1]^2), 0)
  p4[1] <- runif(1, 0, 1)
  p4[2] <- runif(1, -sqrt(1 - p4[1]^2), 0)
  
  knots <- rbind(p1, p2, p3, p4)
  
  d <- rdist(s, knots)
  idx1 <- which(d[1, ] == min(d[1, ]))
  idx2 <- which(d[2, ] == min(d[2, ]))
  
  if (idx1 == idx2) {
    points(knots, col = "dodgerblue", pch = 19)
  }
  prob <- prob + (idx1 == idx2) / reps
}



r <- 2
s <- cbind(c(-r, r), c(0, 0))
area <- pi * r^2

reps <- 500000
prob <- 0
plot(
  s,
  main = "Nearest-knot overlap check (r = 2, 2 knots per quadrant)",
  xlab = "x",
  ylab = "y"
)
abline(h = 0, v = 0)
legend(
  "topright",
  legend = c("Sites s", "Knots when same nearest knot"),
  pch = c(1, 19),
  col = c("black", "dodgerblue"),
  bty = "n"
)
set.seed(1)
for (i in 1:reps) {
  this.knot <- 0
  # nknots <- rpois(4, area / 4)
  # get 2 points in every quadrant of the circle
  nknots <- rep(2, 4)
  knots <- matrix(0, nrow = sum(nknots), ncol = 2)
  
  if (nknots[1] > 0) {
    for (k in 1:nknots[1]) {  # quadrant 1
      this.knot <- this.knot + 1
      knot.x <- runif(1, 0, r)
      knot.y <- runif(1, 0, sqrt(r^2 - knot.x^2))
      knots[this.knot, ] <- c(knot.x, knot.y)
    }
  }
  
  if (nknots[2] > 0) {
    for (k in 1:nknots[2]) {  # quadrant 2
      this.knot <- this.knot + 1
      knot.x <- runif(1, -r, 0)
      knot.y <- runif(1, 0, sqrt(r^2 - knot.x^2))
      knots[this.knot, ] <- c(knot.x, knot.y)
    }
  }
  
  if (nknots[3] > 0) {
    for (k in 1:nknots[3]) {  # quadrant 3
      this.knot <- this.knot + 1
      knot.x <- runif(1, -r, 0)
      knot.y <- runif(1, -sqrt(r^2 - knot.x^2), 0)
      knots[this.knot, ] <- c(knot.x, knot.y)
    }
  }
  
  if (nknots[4] > 0) {
    for (k in 1:nknots[4]) {  # quadrant 4
      this.knot <- this.knot + 1
      knot.x <- runif(1, 0, r)
      knot.y <- runif(1, -sqrt(r^2 - knot.x^2), 0)
      knots[this.knot, ] <- c(knot.x, knot.y)
    }
  }
  
  d <- rdist(s, knots)
  idx1 <- which(d[1, ] == min(d[1, ]))
  idx2 <- which(d[2, ] == min(d[2, ]))
  
  if (idx1 == idx2) {
    points(knots, col = "dodgerblue", pch = 19)
  }
  prob <- prob + (idx1 == idx2) / reps
}


r <- 2
s <- cbind(c(-r, r), c(0, 0))
area <- pi * r^2

reps <- 500000
prob <- 0
plot(
  s,
  ylim = c(-r, r),
  xlim = c(-r, r),
  main = "Nearest-knot overlap check (r = 2, 8 knots over full disk)",
  xlab = "x",
  ylab = "y"
)
abline(h = 0, v = 0)
legend(
  "topright",
  legend = c("Sites s", "First found overlap case"),
  pch = c(1, 19),
  col = c("black", "dodgerblue"),
  bty = "n"
)
set.seed(1)
for (i in 1:reps) {

  nknots <- 8
  knots <- matrix(0, nrow = nknots, ncol = 2)
  
  for (k in 1:nknots) {  # quadrant 1
    this.knot <- this.knot + 1
    knot.x <- runif(1, -r, r)
    knot.y <- runif(1, -sqrt(r^2 - knot.x^2), sqrt(r^2 - knot.x^2))
    knots[k, ] <- c(knot.x, knot.y)
  }
  
  d <- rdist(s, knots)
  idx1 <- which(d[1, ] == min(d[1, ]))
  idx2 <- which(d[2, ] == min(d[2, ]))
  hit <- as.integer(idx1 == idx2)
  prob <- prob + hit / reps
  
  if (hit == 1) {
    points(knots, col = "dodgerblue", pch = 19)
    message("found one case")
    break
  }
}
cat("Estimated overlap probability up to break:", prob, "\n")


r <- 2
s <- cbind(c(-r, r), c(0, 0))
area <- pi * r^2

reps <- 10000000
prob <- 0
plot(
  s,
  ylim = c(-r, r),
  xlim = c(-r, r),
  main = "Nearest-knot overlap check (r = 2, 2 knots per quadrant, long run)",
  xlab = "x",
  ylab = "y"
)
abline(h = 0, v = 0)
legend(
  "topright",
  legend = c("Sites s", "Knots when same nearest knot"),
  pch = c(1, 19),
  col = c("black", "dodgerblue"),
  bty = "n"
)
set.seed(1)
for (i in 1:reps) {
  this.knot <- 0
  # nknots <- rpois(4, area / 4)
  # get 2 points in every quadrant of the circle
  nknots <- c(2, 2, 2, 2)
  knots <- matrix(0, nrow = sum(nknots), ncol = 2)
  
  if (nknots[1] > 0) {
    for (k in 1:nknots[1]) {  # quadrant 1
      this.knot <- this.knot + 1
      knot.x <- runif(1, 0, r)
      knot.y <- runif(1, 0, sqrt(r^2 - knot.x^2))
      knots[this.knot, ] <- c(knot.x, knot.y)
    }
  }
  
  if (nknots[2] > 0) {
    for (k in 1:nknots[2]) {  # quadrant 2
      this.knot <- this.knot + 1
      knot.x <- runif(1, -r, 0)
      knot.y <- runif(1, 0, sqrt(r^2 - knot.x^2))
      knots[this.knot, ] <- c(knot.x, knot.y)
    }
  }
  
  if (nknots[3] > 0) {
    for (k in 1:nknots[3]) {  # quadrant 3
      this.knot <- this.knot + 1
      knot.x <- runif(1, -r, 0)
      knot.y <- runif(1, -sqrt(r^2 - knot.x^2), 0)
      knots[this.knot, ] <- c(knot.x, knot.y)
    }
  }
  
  if (nknots[4] > 0) {
    for (k in 1:nknots[4]) {  # quadrant 4
      this.knot <- this.knot + 1
      knot.x <- runif(1, 0, r)
      knot.y <- runif(1, -sqrt(r^2 - knot.x^2), 0)
      knots[this.knot, ] <- c(knot.x, knot.y)
    }
  }
  
  d <- rdist(s, knots)
  idx1 <- which(d[1, ] == min(d[1, ]))
  idx2 <- which(d[2, ] == min(d[2, ]))
  
  if (idx1 == idx2) {
    points(knots, col = "dodgerblue", pch = 19)
  }
  prob <- prob + (idx1 == idx2) / reps
}