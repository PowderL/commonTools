#' @export
cal_RMSE <- function(obs, pre) {
	# n = length(obs)
	sel <- complete.cases(cbind(obs, pre))
	obs <- obs[sel]
	pre <- pre[sel]
	n = length(obs)
	sqrt(sum((pre - obs)^2) / n)
}

# cal_RMSE <- function(obs, pre, pop) {
# 	n = length(obs)
# 	sqrt(sum(pop*(pre - obs)^2) / sum(pop))
# }
#' @export
cal_cvRMSE <- function(obs, pre) {
	# n = length(obs)
	sel <- complete.cases(cbind(obs, pre))
	obs <- obs[sel]
	pre <- pre[sel]
	n = length(obs)
	sqrt(sum((pre - obs)^2) / n) / mean(obs)
}
#' @export
cal_R2 <- function(obs, pre) {
	sel <- complete.cases(cbind(obs, pre))
	obs <- obs[sel]
	pre <- pre[sel]
	(cor(obs, pre))^2
}

# population weighted
# cal_R2 <- function(obs, pre) {
# 	(cor(obs, pre))^2
# }
#' @export
cal_MFB <- function(obs, pre) { # mean fractional bias
	sel <- complete.cases(cbind(obs, pre))
	obs <- obs[sel]
	pre <- pre[sel]

	mean((pre-obs)/(obs+pre)*2)
}
#' @export
cal_MFE <- function(obs, pre) { # mean fractional error
	sel <- complete.cases(cbind(obs, pre))
	obs <- obs[sel]
	pre <- pre[sel]

	mean(abs(pre-obs)/(obs+pre)*2)
}
#' @export
cal_MNB <- function(obs, pre) { # mean normalized bias
	sel <- complete.cases(cbind(obs, pre))
	obs <- obs[sel]
	pre <- pre[sel]

	mean((pre-obs)/obs)
}
#' @export
cal_MNE <- function(obs, pre) { # mean normalized error
	sel <- complete.cases(cbind(obs, pre))
	obs <- obs[sel]
	pre <- pre[sel]

	mean(abs(pre-obs)/obs)
}
#' @export
cal_NMB <- function(obs, pre) {
	sel <- complete.cases(cbind(obs, pre))
	obs <- obs[sel]
	pre <- pre[sel]

	sum(pre-obs)/sum(obs)
}
#' @export
cal_NME <- function(obs, pre) {
	sel <- complete.cases(cbind(obs, pre))
	obs <- obs[sel]
	pre <- pre[sel]

	sum(abs(pre-obs))/sum(obs)
}
#' @export
cal_slope <- function(obs, pre) {
	sel <- complete.cases(cbind(obs, pre))
	obs <- obs[sel]
	pre <- pre[sel]

	lm.sol <- summary(lm(pre ~ obs))
	lm.sol$coefficients[2,1]
}