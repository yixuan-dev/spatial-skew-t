updateTauTS_AR2 <- function(phi1, phi2, tau, taug, g, res, nparts.tau, prec,
                            z, tau.alpha, tau.beta, skew, att, acc, mh,
                            att.phi1, acc.phi1, mh.phi1,
                            att.phi2, acc.phi2, mh.phi2) {
    ns <- nrow(res)
    nt <- ncol(res)
    nknots <- nrow(tau)

    # 檢查 AR(2) 平穩性條件
    if (!check_ar2_stability(phi1, phi2)) {
        warning(
            "AR(2) parameters (phi1=", phi1, ", phi2=", phi2,
            ") do not satisfy stationarity conditions. Using default values."
        )
        phi1 <- 0
        phi2 <- 0
    }

    tau.star <- gamma.cop(tau, tau.alpha, tau.beta)

    if (nknots == 1) {
        for (t in 1:nt) {
            att[1, t] <- att[1, t] + 1
            cur.lly <- 0.5 * ns * log(tau[1, t]) -
                0.5 * tau[1, t] * quad.form(prec, res[, t])

            can.tau.star <- rnorm(1, tau.star[1, t], mh[1, t])

            # transform back to R+
            can.tau <- gamma.invcop(can.tau.star, tau.alpha, tau.beta)
            if (can.tau < 1e-6) {
                can.tau <- 1e-6
            }

            can.lly <- 0.5 * ns * log(can.tau) -
                0.5 * can.tau * quad.form(prec, res[, t])

            if (skew) {
                cur.llz <- 0.5 * log(tau[1, t]) - 0.5 * tau[1, t] * z[1, t]^2
                can.llz <- 0.5 * log(can.tau) - 0.5 * can.tau * z[1, t]^2
            } else {
                cur.llz <- can.llz <- 0
            }

            R <- can.lly - cur.lly + can.llz - cur.llz

            # AR(2) 時間序列先驗分佈
            # 計算 AR(2) 的創新標準差 (innovation variance)
            gamma0 <- 1 # 標準化邊際方差
            gamma1 <- phi1 * gamma0 / (1 - phi2)
            gamma2 <- phi1 * gamma1 + phi2 * gamma0
            innov_var <- gamma0 - phi1 * gamma1 - phi2 * gamma2

            # 確保創新方差為正
            if (innov_var <= 1e-6) {
                innov_var <- 1e-6
            }
            innov_sd <- sqrt(innov_var)

            if (t > 2) {
                mean <- phi1 * tau.star[1, (t - 1)] + phi2 * tau.star[1, (t - 2)]
                sd <- innov_sd
            } else if (t == 2) {
                mean <- phi1 * tau.star[1, 1]
                sd <- innov_sd # AR(2) 的創新標準差
            } else { # t == 1
                mean <- 0
                sd <- sqrt(gamma0) # 邊際標準差
            }

            # evaluate the prior
            R <- R + dnorm(can.tau.star, mean, sd, log = TRUE) -
                dnorm(tau.star[1, t], mean, sd, log = TRUE)

            # 考慮下一時刻的影響
            if (t < nt) {
                tau.star.next <- tau.star[1, t + 1]

                if (t == 1) {
                    next_mean_can <- phi1 * can.tau.star
                    next_mean_cur <- phi1 * tau.star[1, t]
                } else {
                    next_mean_can <- phi1 * can.tau.star + phi2 * tau.star[1, t - 1]
                    next_mean_cur <- phi1 * tau.star[1, t] + phi2 * tau.star[1, t - 1]
                }

                R <- R + dnorm(tau.star.next, next_mean_can, 1, log = TRUE) -
                    dnorm(tau.star.next, next_mean_cur, 1, log = TRUE)
            }

            if (!is.na(R)) {
                if (log(runif(1)) < R) {
                    tau.star[1, t] <- can.tau.star
                    tau[1, t] <- can.tau
                    taug[, t] <- rep(can.tau, ns)
                    acc[1, t] <- acc[1, t] + 1
                }
            }
        } # end t
    } else { # nknots > 1
        for (t in 1:nt) {
            cur.lly <- 0.5 * sum(log(taug[, t])) -
                0.5 * quad.form(prec, sqrt(taug[, t]) * res[, t])

            for (k in 1:nknots) {
                att[k, t] <- att[k, t] + 1
                these <- which(g[, t] == k)
                nparts <- length(these)
                nparts.tau[k, t] <- nparts

                can.tau.star <- tau.star[, t]
                can.tau.star[k] <- rnorm(1, tau.star[k, t], mh[k, t])

                # transform back to R+
                can.tau <- gamma.invcop(can.tau.star, tau.alpha, tau.beta)
                if (can.tau[k] < 1e-6) {
                    can.tau[k] <- 1e-6
                }
                can.taug <- can.tau[g[, t]]

                can.lly <- 0.5 * sum(log(can.taug)) -
                    0.5 * quad.form(prec, sqrt(can.taug) * res[, t])

                if (skew) {
                    cur.llz <- 0.5 * log(tau[k, t]) - 0.5 * tau[k, t] * z[k, t]^2
                    can.llz <- 0.5 * log(can.tau[k]) - 0.5 * can.tau[k] * z[k, t]^2
                } else {
                    cur.llz <- can.llz <- 0
                }

                R <- can.lly - cur.lly + can.llz - cur.llz

                # AR(2) 時間序列先驗分佈 (修正: 使用 [k, t] 索引)
                # 計算 AR(2) 的創新標準差 (innovation variance)
                gamma0 <- 1 # 標準化邊際方差
                gamma1 <- phi1 * gamma0 / (1 - phi2)
                gamma2 <- phi1 * gamma1 + phi2 * gamma0
                innov_var <- gamma0 - phi1 * gamma1 - phi2 * gamma2

                # 確保創新方差為正
                if (innov_var <= 1e-6) {
                    innov_var <- 1e-6
                }
                innov_sd <- sqrt(innov_var)

                if (t > 2) {
                    mean <- phi1 * tau.star[k, (t - 1)] + phi2 * tau.star[k, (t - 2)]
                    sd <- innov_sd
                } else if (t == 2) {
                    mean <- phi1 * tau.star[k, 1]
                    sd <- innov_sd # AR(2) 的創新標準差
                } else { # t == 1
                    mean <- 0
                    sd <- sqrt(gamma0) # 邊際標準差
                }

                R <- R + dnorm(can.tau.star[k], mean, sd, log = TRUE) -
                    dnorm(tau.star[k, t], mean, sd, log = TRUE)

                # 考慮下一時刻的影響 (修正: 使用 [k, ...] 索引)
                if (t < nt) {
                    tau.star.next <- tau.star[k, t + 1]

                    if (t == 1) {
                        next_mean_can <- phi1 * can.tau.star[k]
                        next_mean_cur <- phi1 * tau.star[k, t]
                    } else {
                        next_mean_can <- phi1 * can.tau.star[k] + phi2 * tau.star[k, t - 1]
                        next_mean_cur <- phi1 * tau.star[k, t] + phi2 * tau.star[k, t - 1]
                    }

                    R <- R + dnorm(tau.star.next, next_mean_can, innov_sd, log = TRUE) -
                        dnorm(tau.star.next, next_mean_cur, innov_sd, log = TRUE)
                }

                if (!is.na(R)) {
                    if (log(runif(1)) < R) {
                        acc[k, t] <- acc[k, t] + 1
                        tau.star[k, t] <- can.tau.star[k]
                        tau[k, t] <- can.tau[k]
                        taug[these, t] <- can.tau[k]
                        cur.lly <- can.lly
                    }
                }
            } # end k
        } # end t
    } # fi nknots > 1

    # 更新 phi1 和 phi2
    phi.update <- updatePhiAR2TS(
        data = tau.star, phi1 = phi1, phi2 = phi2, day.mar = 2,
        att.phi1 = att.phi1, acc.phi1 = acc.phi1, mh.phi1 = mh.phi1,
        att.phi2 = att.phi2, acc.phi2 = acc.phi2, mh.phi2 = mh.phi2
    )
    phi1 <- phi.update$phi1
    phi2 <- phi.update$phi2
    acc.phi1 <- phi.update$acc.phi1
    att.phi1 <- phi.update$att.phi1
    acc.phi2 <- phi.update$acc.phi2
    att.phi2 <- phi.update$att.phi2

    results <- list(
        tau = tau, taug = taug,
        phi1 = phi1, phi2 = phi2,
        acc = acc, att = att,
        acc.phi1 = acc.phi1, att.phi1 = att.phi1,
        acc.phi2 = acc.phi2, att.phi2 = att.phi2
    )

    return(results)
}

# 更新 AR(2) 的 phi 參數
updatePhiAR2TS <- function(data, phi1, phi2, day.mar,
                           att.phi1, acc.phi1, mh.phi1,
                           att.phi2, acc.phi2, mh.phi2,
                           prior.mean = 0, prior.sd = 0.5) {
    att.phi1 <- att.phi1 + 1
    att.phi2 <- att.phi2 + 1

    nt <- dim(data)[day.mar]

    # 如果時間點數量不足,直接返回
    if (nt <= 2) {
        return(list(
            phi1 = phi1, phi2 = phi2,
            att.phi1 = att.phi1, acc.phi1 = acc.phi1,
            att.phi2 = att.phi2, acc.phi2 = acc.phi2
        ))
    }

    # 提取不同落後階數的數據
    if (day.mar == 2) {
        data.current <- data[, 3:nt, drop = FALSE]
        data.lag1 <- data[, 2:(nt - 1), drop = FALSE]
        data.lag2 <- data[, 1:(nt - 2), drop = FALSE]
    } else if (day.mar == 3) {
        data.current <- data[, , 3:nt, drop = FALSE]
        data.lag1 <- data[, , 2:(nt - 1), drop = FALSE]
        data.lag2 <- data[, , 1:(nt - 2), drop = FALSE]
    }

    # Determine prior means for phi1 and phi2
    p.m.1 <- if (length(prior.mean) >= 1) prior.mean[1] else 0
    p.m.2 <- if (length(prior.mean) >= 2) prior.mean[2] else p.m.1
    p.s.1 <- if (length(prior.sd) >= 1) prior.sd[1] else 0.5
    p.s.2 <- if (length(prior.sd) >= 2) prior.sd[2] else p.s.1

    # 當前模型的條件均值
    cur.mean <- phi1 * data.lag1 + phi2 * data.lag2

    # 更新 phi1
    can.phi1 <- rnorm(1, phi1, mh.phi1)

    # 檢查穩定性
    if (check_ar2_stability(can.phi1, phi2)) {
        can.mean <- can.phi1 * data.lag1 + phi2 * data.lag2

        R <- sum(dnorm(data.current, can.mean, 1, log = TRUE)) -
            sum(dnorm(data.current, cur.mean, 1, log = TRUE)) +
            dnorm(can.phi1, p.m.1, p.s.1, log = TRUE) -
            dnorm(phi1, p.m.1, p.s.1, log = TRUE)

        if (!is.na(R) && log(runif(1)) < R) {
            phi1 <- can.phi1
            acc.phi1 <- acc.phi1 + 1
            cur.mean <- can.mean # 更新當前均值供 phi2 使用
        }
    }

    # 更新 phi2 (使用更新後的 phi1)
    can.phi2 <- rnorm(1, phi2, mh.phi2)

    # 檢查穩定性
    if (check_ar2_stability(phi1, can.phi2)) {
        can.mean <- phi1 * data.lag1 + can.phi2 * data.lag2

        R <- sum(dnorm(data.current, can.mean, 1, log = TRUE)) -
            sum(dnorm(data.current, cur.mean, 1, log = TRUE)) +
            dnorm(can.phi2, p.m.2, p.s.2, log = TRUE) -
            dnorm(phi2, p.m.2, p.s.2, log = TRUE)

        if (!is.na(R) && log(runif(1)) < R) {
            phi2 <- can.phi2
            acc.phi2 <- acc.phi2 + 1
        }
    }

    results <- list(
        phi1 = phi1, phi2 = phi2,
        att.phi1 = att.phi1, acc.phi1 = acc.phi1,
        att.phi2 = att.phi2, acc.phi2 = acc.phi2
    )

    return(results)
}

updateZTS_AR2 <- function(z, zg, y, lambda, x.beta,
                          phi1, phi2, tau, taug, g, prec,
                          acc, att, mh,
                          acc.phi1, att.phi1, mh.phi1,
                          acc.phi2, att.phi2, mh.phi2) {
    nt <- ncol(z)
    nknots <- nrow(z)

    # transform via copula to normal
    sig <- 1 / sqrt(tau)
    z.star <- hn.cop(x = z, sig = sig)

    # 檢查 AR(2) 平穩性條件
    if (!check_ar2_stability(phi1, phi2)) {
        warning(
            "AR(2) parameters (phi1=", phi1, ", phi2=", phi2,
            ") do not satisfy stationarity conditions. Using default values."
        )
        phi1 <- 0
        phi2 <- 0
    }

    # 計算 AR(2) 的創新標準差 (innovation variance)
    gamma0 <- 1 # 標準化邊際方差
    gamma1 <- phi1 * gamma0 / (1 - phi2)
    gamma2 <- phi1 * gamma1 + phi2 * gamma0
    innov_var <- gamma0 - phi1 * gamma1 - phi2 * gamma2

    # 確保創新方差為正
    if (innov_var <= 1e-6) {
        innov_var <- 1e-6
    }
    innov_sd <- sqrt(innov_var)

    for (t in 1:nt) {
        taug.t <- sqrt(taug[, t])
        mu.t <- x.beta[, t] + lambda * zg[, t]
        cur.res <- y[, t] - mu.t
        cur.lly <- -0.5 * quad.form(prec, taug.t * cur.res)

        for (k in 1:nknots) {
            att[k, t] <- att[k, t] + 1
            these <- which(g[, t] == k)
            can.z.star <- z.star[, t]
            can.z.star[k] <- rnorm(1, z.star[k, t], mh[k, t])

            # transform back to R+
            can.z <- hn.invcop(x = can.z.star, sig = sig[, t])
            if (can.z[k] < 1e-6) { # numerical stability
                can.z[k] <- 1e-6
            }
            can.zg <- can.z[g[, t]] # ns long

            can.mu.t <- x.beta[, t] + lambda * can.zg
            can.res <- y[, t] - can.mu.t
            can.lly <- -0.5 * quad.form(prec, taug.t * can.res)

            # prior
            if (t > 2) {
                mean <- phi1 * z.star[k, (t - 1)] + phi2 * z.star[k, (t - 2)]
                sd <- innov_sd
            } else if (t == 2) {
                mean <- phi1 * z.star[k, 1]
                sd <- innov_sd
            } else {
                mean <- 0
                sd <- sqrt(gamma0)
            }

            R <- can.lly - cur.lly +
                dnorm(can.z.star[k], mean, sd, log = TRUE) -
                dnorm(z.star[k, t], mean, sd, log = TRUE)

            if (t < nt) {
                z.star.next <- z.star[k, t + 1]

                if (t == 1) {
                    next_mean_can <- phi1 * can.z.star[k]
                    next_mean_cur <- phi1 * z.star[k, t]
                } else {
                    next_mean_can <- phi1 * can.z.star[k] + phi2 * z.star[k, t - 1]
                    next_mean_cur <- phi1 * z.star[k, t] + phi2 * z.star[k, t - 1]
                }

                R <- R + dnorm(z.star.next, next_mean_can, innov_sd, log = TRUE) -
                    dnorm(z.star.next, next_mean_cur, innov_sd, log = TRUE)
            }

            if (!is.na(R)) {
                if (log(runif(1)) < R) {
                    acc[k, t] <- acc[k, t] + 1
                    z[k, t] <- can.z[k]
                    zg[these, t] <- can.z[k]
                    z.star[k, t] <- can.z.star[k]
                    cur.lly <- can.lly
                }
            }
        }
    }

    # 更新 phi1 和 phi2
    phi.update <- updatePhiAR2TS(
        data = z.star, phi1 = phi1, phi2 = phi2, day.mar = 2,
        att.phi1 = att.phi1, acc.phi1 = acc.phi1, mh.phi1 = mh.phi1,
        att.phi2 = att.phi2, acc.phi2 = acc.phi2, mh.phi2 = mh.phi2
    )
    phi1 <- phi.update$phi1
    phi2 <- phi.update$phi2
    acc.phi1 <- phi.update$acc.phi1
    att.phi1 <- phi.update$att.phi1
    acc.phi2 <- phi.update$acc.phi2
    att.phi2 <- phi.update$att.phi2

    results <- list(
        z = z, zg = zg,
        phi1 = phi1, phi2 = phi2,
        acc = acc, att = att,
        acc.phi1 = acc.phi1, att.phi1 = att.phi1,
        acc.phi2 = acc.phi2, att.phi2 = att.phi2
    )

    return(results)
}

updateKnotsTS_AR2 <- function(phi1, phi2, knots, g, ts, tau, z, s, min.s,
                              max.s, x.beta, lambda, y, prec, att,
                              acc, mh, update.prop = 1,
                              att.phi1, acc.phi1, mh.phi1,
                              att.phi2, acc.phi2, mh.phi2) {
    ns <- nrow(y)
    nt <- ncol(y)
    nknots <- dim(knots)[1]
    avgparts <- rep(0, nt)
    taug <- zg <- matrix(NA, ns, nt)

    # recalculate knots.star at the beginning
    knots.star <- array(NA, dim = c(nknots, 2, nt))
    knots.star[, 1, ] <- transform$probit(knots[, 1, ],
        lower = min.s[1],
        upper = max.s[1]
    )
    knots.star[, 2, ] <- transform$probit(knots[, 2, ],
        lower = min.s[2],
        upper = max.s[2]
    )

    if (!ts) { # will be returning these with the function results
        phi1 <- phi2 <- 0
        att.phi1 <- acc.phi1 <- mh.phi1 <- 0
        att.phi2 <- acc.phi2 <- mh.phi2 <- 0
    } else {
        # 檢查 AR(2) 平穩性條件
        if (!check_ar2_stability(phi1, phi2)) {
            warning(
                "AR(2) parameters (phi1=", phi1, ", phi2=", phi2,
                ") do not satisfy stationarity conditions. Using default values."
            )
            phi1 <- 0
            phi2 <- 0
        }
    }

    # 計算 AR(2) 的創新標準差 (innovation variance)
    gamma0 <- 1 # 標準化邊際方差
    gamma1 <- phi1 * gamma0 / (1 - phi2)
    gamma2 <- phi1 * gamma1 + phi2 * gamma0
    innov_var <- gamma0 - phi1 * gamma1 - phi2 * gamma2

    # 確保創新方差為正
    if (innov_var <= 1e-6) {
        innov_var <- 1e-6
    }
    innov_sd <- sqrt(innov_var)

    for (t in 1:nt) {
        att[, t] <- att[, t] + 1
        taug.t <- tau[g[, t], t]
        y.t <- y[, t]
        x.beta.t <- x.beta[, t]
        zg.t <- z[g[, t], t]
        res.t <- y.t - x.beta.t - lambda * zg.t
        cur.lly <- 0.5 * sum(log(taug.t)) -
            0.5 * quad.form(prec, sqrt(taug.t) * res.t)

        can.knots.star <- cur.knots.star <- knots.star[, , t]
        for (k in 1:nknots) {
            can.knots.star[k, ] <- cur.knots.star[k, ] + mh[k, t] * rnorm(2)
        }
        can.knots <- matrix(NA, nknots, 2)
        can.knots[, 1] <- transform$inv.probit(can.knots.star[, 1],
            lower = min.s[1],
            upper = max.s[1]
        )
        can.knots[, 2] <- transform$inv.probit(can.knots.star[, 2],
            lower = min.s[2],
            upper = max.s[2]
        )

        # recreate the partition
        can.g <- mem(s, can.knots)
        can.taug <- tau[can.g, t]
        can.zg <- z[can.g, t]
        can.res <- y.t - x.beta.t - lambda * can.zg
        can.lly <- 0.5 * sum(log(can.taug)) -
            0.5 * quad.form(prec, sqrt(can.taug) * can.res)

        # remember, when not a TS, phi = 0
        if (ts & (t > 2)) {
            mean <- phi1 * knots.star[, , (t - 1)] + phi2 * knots.star[, , (t - 2)]
            sd <- innov_sd
        } else if (ts & (t == 2)) {
            mean <- phi1 * knots.star[, , 1]
            sd <- innov_sd
        } else {
            mean <- 0
            sd <- sqrt(gamma0)
        }

        R <- can.lly - cur.lly +
            sum(dnorm(can.knots.star, mean, sd, log = TRUE)) -
            sum(dnorm(cur.knots.star, mean, sd, log = TRUE))

        # time series also needs to adjust R to account for next day
        if (ts & (t < nt)) {
            knots.star.next <- knots.star[, , t + 1]

            if (t == 1) {
                next_mean_can <- phi1 * can.knots.star
                next_mean_cur <- phi1 * cur.knots.star
            } else {
                next_mean_can <- phi1 * can.knots.star + phi2 * knots.star[, , t - 1]
                next_mean_cur <- phi1 * cur.knots.star + phi2 * knots.star[, , t - 1]
            }

            R <- R + sum(dnorm(knots.star.next, next_mean_can, innov_sd, log = TRUE)) -
                sum(dnorm(knots.star.next, next_mean_cur, innov_sd, log = TRUE))
        }

        if (!is.na(R)) {
            if (log(runif(1)) < R) {
                acc[, t] <- acc[, t] + 1
                knots.star[, , t] <- can.knots.star
                knots[, , t] <- can.knots
                g[, t] <- can.g
            }
        }

        zg[, t] <- z[g[, t], t]
        taug[, t] <- tau[g[, t], t]
    }

    if (ts) {
        phi.update <- updatePhiAR2TS(
            data = knots.star, phi1 = phi1, phi2 = phi2, day.mar = 3,
            att.phi1 = att.phi1, acc.phi1 = acc.phi1, mh.phi1 = mh.phi1,
            att.phi2 = att.phi2, acc.phi2 = acc.phi2, mh.phi2 = mh.phi2
        )
        phi1 <- phi.update$phi1
        phi2 <- phi.update$phi2
        acc.phi1 <- phi.update$acc.phi1
        att.phi1 <- phi.update$att.phi1
        acc.phi2 <- phi.update$acc.phi2
        att.phi2 <- phi.update$att.phi2
    }

    results <- list(
        knots.star = knots.star, knots = knots, g = g, taug = taug,
        zg = zg, acc = acc, att = att,
        phi1 = phi1, phi2 = phi2,
        acc.phi1 = acc.phi1, att.phi1 = att.phi1,
        acc.phi2 = acc.phi2, att.phi2 = att.phi2
    )

    return(results)
}
