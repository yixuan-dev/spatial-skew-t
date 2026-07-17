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
    ar2.moments <- get_ar2_standard_moments(phi1, phi2)

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
            lag1 <- if (t >= 2) tau.star[1, t - 1] else NULL
            lag2 <- if (t >= 3) tau.star[1, t - 2] else NULL
            R <- R +
                ar2_conditional_log_density(
                    can.tau.star, t, lag1, lag2,
                    phi1 = phi1, phi2 = phi2, moments = ar2.moments
                ) -
                ar2_conditional_log_density(
                    tau.star[1, t], t, lag1, lag2,
                    phi1 = phi1, phi2 = phi2, moments = ar2.moments
                )

            # 考慮下一時刻的影響 (tau.star[1, t] 在 t+1 條件密度中為 lag-1)
            if (t < nt) {
                tau.star.next <- tau.star[1, t + 1]
                next.lag2 <- if (t >= 2) tau.star[1, t - 1] else NULL

                R <- R +
                    ar2_conditional_log_density(
                        tau.star.next, t + 1, can.tau.star, next.lag2,
                        phi1 = phi1, phi2 = phi2, moments = ar2.moments
                    ) -
                    ar2_conditional_log_density(
                        tau.star.next, t + 1, tau.star[1, t], next.lag2,
                        phi1 = phi1, phi2 = phi2, moments = ar2.moments
                    )
            }

            # Fix 2: tau.star[1, t] 在 t+2 條件密度中為 lag-2，缺此項會破壞
            # AR(2) 先驗的細緻平衡 (lag-2 ACF 會偏離 gamma_2)。
            # The companion t+1 block above evaluates p(y_2 | y_1) via the
            # t == 2 branch of get_ar2_conditional_params, which now returns
            # the stationary form N(gamma_1*y_1, 1 - gamma_1^2) -- see Fix 1.
            if ((t + 2) <= nt) {
                tau.star.next2 <- tau.star[1, t + 2]
                R <- R +
                    ar2_conditional_log_density(
                        tau.star.next2, t + 2,
                        tau.star[1, t + 1], can.tau.star,
                        phi1 = phi1, phi2 = phi2, moments = ar2.moments
                    ) -
                    ar2_conditional_log_density(
                        tau.star.next2, t + 2,
                        tau.star[1, t + 1], tau.star[1, t],
                        phi1 = phi1, phi2 = phi2, moments = ar2.moments
                    )
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
                lag1 <- if (t >= 2) tau.star[k, t - 1] else NULL
                lag2 <- if (t >= 3) tau.star[k, t - 2] else NULL
                R <- R +
                    ar2_conditional_log_density(
                        can.tau.star[k], t, lag1, lag2,
                        phi1 = phi1, phi2 = phi2, moments = ar2.moments
                    ) -
                    ar2_conditional_log_density(
                        tau.star[k, t], t, lag1, lag2,
                        phi1 = phi1, phi2 = phi2, moments = ar2.moments
                    )

                # 考慮下一時刻的影響 (tau.star[k, t] 在 t+1 為 lag-1)
                if (t < nt) {
                    tau.star.next <- tau.star[k, t + 1]
                    next.lag2 <- if (t >= 2) tau.star[k, t - 1] else NULL

                    R <- R +
                        ar2_conditional_log_density(
                            tau.star.next, t + 1, can.tau.star[k], next.lag2,
                            phi1 = phi1, phi2 = phi2, moments = ar2.moments
                        ) -
                        ar2_conditional_log_density(
                            tau.star.next, t + 1, tau.star[k, t], next.lag2,
                            phi1 = phi1, phi2 = phi2, moments = ar2.moments
                        )
                }

                # Fix 2: tau.star[k, t] 在 t+2 為 lag-2
                # The t == 2 branch of get_ar2_conditional_params now
                # returns the stationary N(gamma_1*y_1, 1 - gamma_1^2);
                # see Fix 1.
                if ((t + 2) <= nt) {
                    tau.star.next2 <- tau.star[k, t + 2]
                    R <- R +
                        ar2_conditional_log_density(
                            tau.star.next2, t + 2,
                            tau.star[k, t + 1], can.tau.star[k],
                            phi1 = phi1, phi2 = phi2, moments = ar2.moments
                        ) -
                        ar2_conditional_log_density(
                            tau.star.next2, t + 2,
                            tau.star[k, t + 1], tau.star[k, t],
                            phi1 = phi1, phi2 = phi2, moments = ar2.moments
                        )
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

    # 提取不同落後階數的數據 (含 t = 1, 2 初始切片)
    if (day.mar == 2) {
        data.current <- data[, 3:nt, drop = FALSE]
        data.lag1 <- data[, 2:(nt - 1), drop = FALSE]
        data.lag2 <- data[, 1:(nt - 2), drop = FALSE]
        data.t1 <- data[, 1]
        data.t2 <- data[, 2]
    } else if (day.mar == 3) {
        data.current <- data[, , 3:nt, drop = FALSE]
        data.lag1 <- data[, , 2:(nt - 1), drop = FALSE]
        data.lag2 <- data[, , 1:(nt - 2), drop = FALSE]
        data.t1 <- data[, , 1]
        data.t2 <- data[, , 2]
    }

    # Determine prior means for phi1 and phi2
    p.m.1 <- if (length(prior.mean) >= 1) prior.mean[1] else 0
    p.m.2 <- if (length(prior.mean) >= 2) prior.mean[2] else p.m.1
    p.s.1 <- if (length(prior.sd) >= 1) prior.sd[1] else 0.5
    p.s.2 <- if (length(prior.sd) >= 2) prior.sd[2] else p.s.1

    # 精確平穩似然: t >= 3 的轉移密度 + 依賴 phi 的初始塊 p(Y_2 | Y_1, phi)。
    # 只用轉移密度 (條件似然) 在 AR(1) 是恆等式、在 AR(2) 會使 phi-block 與
    # latent-block 目標不同的分佈 (Gibbs blocks 不相容); 詳見
    # tex/ar2_phi_exact_likelihood/. Y_1 的邊際 N(0,1) 不含 phi, 於接受比中消去。
    cur.moments <- get_ar2_standard_moments(phi1, phi2)
    cur.ll <- ar2_stationary_loglik(
        data.current, data.lag1, data.lag2, data.t1, data.t2,
        phi1 = phi1, phi2 = phi2, moments = cur.moments
    )

    # 更新 phi1 (silent-rejection scheme)
    # Propose phi1.can ~ N(phi1, mh.phi1^2) untruncated. If the candidate
    # falls outside the AR(2) stationarity region we silently leave the
    # state unchanged. Per Prop. 3.3 of tex/where_we_are/where_we_are.tex,
    # this kernel has the posterior pi restricted to S as its invariant
    # distribution -- no Hastings adjustment is needed for a symmetric
    # untruncated proposal combined with silent rejection.
    can.phi1 <- rnorm(1, phi1, mh.phi1)

    if (check_ar2_stability(can.phi1, phi2)) {
        can.moments <- tryCatch(
            get_ar2_standard_moments(can.phi1, phi2),
            error = function(e) NULL
        )

        if (!is.null(can.moments)) {
            can.ll <- ar2_stationary_loglik(
                data.current, data.lag1, data.lag2, data.t1, data.t2,
                phi1 = can.phi1, phi2 = phi2, moments = can.moments
            )

            R <- can.ll - cur.ll +
                dnorm(can.phi1, p.m.1, p.s.1, log = TRUE) -
                dnorm(phi1, p.m.1, p.s.1, log = TRUE)

            if (!is.na(R) && log(runif(1)) < R) {
                phi1 <- can.phi1
                acc.phi1 <- acc.phi1 + 1
                cur.ll <- can.ll
            }
        }
    }

    # 更新 phi2 (使用更新後的 phi1; same silent-rejection scheme)
    can.phi2 <- rnorm(1, phi2, mh.phi2)

    if (check_ar2_stability(phi1, can.phi2)) {
        can.moments <- tryCatch(
            get_ar2_standard_moments(phi1, can.phi2),
            error = function(e) NULL
        )

        if (!is.null(can.moments)) {
            can.ll <- ar2_stationary_loglik(
                data.current, data.lag1, data.lag2, data.t1, data.t2,
                phi1 = phi1, phi2 = can.phi2, moments = can.moments
            )

            R <- can.ll - cur.ll +
                dnorm(can.phi2, p.m.2, p.s.2, log = TRUE) -
                dnorm(phi2, p.m.2, p.s.2, log = TRUE)

            if (!is.na(R) && log(runif(1)) < R) {
                phi2 <- can.phi2
                acc.phi2 <- acc.phi2 + 1
            }
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
    ar2.moments <- get_ar2_standard_moments(phi1, phi2)

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

            lag1 <- if (t >= 2) z.star[k, t - 1] else NULL
            lag2 <- if (t >= 3) z.star[k, t - 2] else NULL
            R <- can.lly - cur.lly +
                ar2_conditional_log_density(
                    can.z.star[k], t, lag1, lag2,
                    phi1 = phi1, phi2 = phi2, moments = ar2.moments
                ) -
                ar2_conditional_log_density(
                    z.star[k, t], t, lag1, lag2,
                    phi1 = phi1, phi2 = phi2, moments = ar2.moments
                )

            # z.star[k, t] 在 t+1 條件密度中為 lag-1
            if (t < nt) {
                z.star.next <- z.star[k, t + 1]
                next.lag2 <- if (t >= 2) z.star[k, t - 1] else NULL

                R <- R +
                    ar2_conditional_log_density(
                        z.star.next, t + 1, can.z.star[k], next.lag2,
                        phi1 = phi1, phi2 = phi2, moments = ar2.moments
                    ) -
                    ar2_conditional_log_density(
                        z.star.next, t + 1, z.star[k, t], next.lag2,
                        phi1 = phi1, phi2 = phi2, moments = ar2.moments
                    )
            }

            # Fix 2: z.star[k, t] 在 t+2 條件密度中為 lag-2
            # The t == 2 branch of get_ar2_conditional_params now returns
            # the stationary N(gamma_1*y_1, 1 - gamma_1^2); see Fix 1.
            if ((t + 2) <= nt) {
                z.star.next2 <- z.star[k, t + 2]
                R <- R +
                    ar2_conditional_log_density(
                        z.star.next2, t + 2,
                        z.star[k, t + 1], can.z.star[k],
                        phi1 = phi1, phi2 = phi2, moments = ar2.moments
                    ) -
                    ar2_conditional_log_density(
                        z.star.next2, t + 2,
                        z.star[k, t + 1], z.star[k, t],
                        phi1 = phi1, phi2 = phi2, moments = ar2.moments
                    )
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
    ar2.moments <- get_ar2_standard_moments(phi1, phi2)

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

        if (ts) {
            lag1 <- if (t >= 2) knots.star[, , t - 1] else NULL
            lag2 <- if (t >= 3) knots.star[, , t - 2] else NULL
            prior.can <- ar2_conditional_log_density(
                can.knots.star, t, lag1, lag2,
                phi1 = phi1, phi2 = phi2, moments = ar2.moments
            )
            prior.cur <- ar2_conditional_log_density(
                cur.knots.star, t, lag1, lag2,
                phi1 = phi1, phi2 = phi2, moments = ar2.moments
            )
        } else {
            prior.can <- dnorm(can.knots.star, 0, sqrt(ar2.moments$gamma0), log = TRUE)
            prior.cur <- dnorm(cur.knots.star, 0, sqrt(ar2.moments$gamma0), log = TRUE)
        }

        R <- can.lly - cur.lly +
            sum(prior.can) - sum(prior.cur)

        # time series also needs to adjust R to account for next day
        # (knots.star[, , t] 在 t+1 為 lag-1)
        if (ts & (t < nt)) {
            knots.star.next <- knots.star[, , t + 1]
            next.lag2 <- if (t >= 2) knots.star[, , t - 1] else NULL

            R <- R +
                sum(ar2_conditional_log_density(
                    knots.star.next, t + 1, can.knots.star, next.lag2,
                    phi1 = phi1, phi2 = phi2, moments = ar2.moments
                )) -
                sum(ar2_conditional_log_density(
                    knots.star.next, t + 1, cur.knots.star, next.lag2,
                    phi1 = phi1, phi2 = phi2, moments = ar2.moments
                ))
        }

        # Fix 2: knots.star[, , t] 在 t+2 為 lag-2
        # The t == 2 branch of get_ar2_conditional_params now returns the
        # stationary N(gamma_1*y_1, 1 - gamma_1^2); see Fix 1.
        if (ts & ((t + 2) <= nt)) {
            knots.star.next2 <- knots.star[, , t + 2]
            R <- R +
                sum(ar2_conditional_log_density(
                    knots.star.next2, t + 2,
                    knots.star[, , t + 1], can.knots.star,
                    phi1 = phi1, phi2 = phi2, moments = ar2.moments
                )) -
                sum(ar2_conditional_log_density(
                    knots.star.next2, t + 2,
                    knots.star[, , t + 1], cur.knots.star,
                    phi1 = phi1, phi2 = phi2, moments = ar2.moments
                ))
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
