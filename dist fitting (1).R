setwd("C:\\Users\\HYEWON\\Desktop\\Lab\\감염병 모델링")
rm=ls()

library(stats4); library(dplyr); library(ggplot2); library(MASS); library(fitdistrplus)

load(file="distribution.rda")
dist<-as.data.frame(dist)

expanded_periods <- rep(dist$incubation.period, dist$frequency)
plot(expanded_periods)
plot(dist)
# test 1
short <- expanded_periods[expanded_periods <= 150]
#plot(short)

long <- expanded_periods[expanded_periods > 150]
#plot(long)

gamma_fit <- fitdist(short, distr = "gamma", method = "mle")
gamma_params <- gamma_fit$estimate
print(gamma_params) #shape=0.046356146 , rate=0.003334965 
plot(gamma_fit)

#  Log-normal 분포 MLE
lognorm_fit <- fitdist(short, "lnorm", method = "mle")
lognorm_params <- lognorm_fit$estimate
print(lognorm_params) #meanlog=-11.63254, sdlog=12.28972 
plot(lognorm_fit)

#  Weibull 분포 MLE
weibull_fit <- fitdist(short, "weibull", method = "mle")
weibull_params <- weibull_fit$estimate
print(weibull_params) #shape=0.092768993, scale=0.004444189
plot(weibull_fit)


normal_fit <- fitdist(long, "norm", method='mle')
normal_params <- normal_fit$estimate

##########################################################################
# 여기부터!!!!!
# 저장할 결과 변수 초기화
best_cutoff <- NULL
min_total_aic <- Inf  # 최소 AIC 값을 저장할 변수

for (cutoff in seq(150, 210, by=5)) {
  
  # short과 long 데이터 분리
  short <- expanded_periods[expanded_periods <= cutoff]
  long <- expanded_periods[expanded_periods > cutoff]
  
  # short 부분의 분포 적합 및 AIC 계산
  gamma_fit <- fitdist(short, distr = "gamma", method = "mle")
  lognorm_fit <- fitdist(short, "lnorm", method = "mle")
  weibull_fit <- fitdist(short, "weibull", method = "mle")
  
  gamma_aic <- AIC(gamma_fit)
  lognorm_aic <- AIC(lognorm_fit)
  weibull_aic <- AIC(weibull_fit)
  
  # 가장 낮은 AIC 값을 가진 short의 분포 선택
  short_aic_values <- c(gamma = gamma_aic, lognorm = lognorm_aic, weibull = weibull_aic)
  best_short_dist <- names(which.min(short_aic_values))  # 최적 분포 이름
  short_best_aic <- min(short_aic_values)
  
  # long 부분의 분포 적합 및 AIC 계산
  normal_fit <- fitdist(long, "norm", method='mle')
  gamma_long_fit<-fitdist(long, "gamma", method='mle')
  lognorm_long_fit<-fitdist(long, "lnorm", method='mle')
  
  normal_long_aic <- AIC(normal_fit)
  gamma_long_aic <- AIC(gamma_long_fit)
  lognorm_long_aic <- AIC(lognorm_long_fit)
  
  long_aic_values <- c(normal = normal_long_aic, gamma = gamma_long_aic, lognorm = lognorm_long_aic)
  best_long_dist <- names(which.min(long_aic_values))  # 최적 분포 이름
  long_best_aic <- min(long_aic_values)
  
  # 전체 AIC 합산
  total_aic <- short_best_aic + long_best_aic
  
  # 최소 AIC 갱신
  if (total_aic < min_total_aic) {
    min_total_aic <- total_aic
    best_cutoff <- cutoff
    best_short_fit <- switch(best_short_dist,
                             gamma = gamma_fit,
                             lognorm = lognorm_fit,
                             weibull = weibull_fit)
    best_long_fit <- switch(best_long_dist,
                            normal = normal_fit,
                            gamma = gamma_long_fit,
                            lognorm = lognorm_long_fit)
    best_short_dist_name <- paste(best_short_dist, best_long_dist, sep = "_")
  }
}

# 최적 cutoff 출력
cat(sprintf("Best Cutoff: %d with Minimum Total AIC: %.2f\n", best_cutoff, min_total_aic))

# 최적 short 분포의 파라미터 출력
cat(sprintf("Best Short Distribution: %s\n", best_short_fit))
print(best_short_fit$estimate)

# 최적 long 분포 (정규 분포) 파라미터 출력
cat(sprintf("Best long Distribution: %s\n", best_long_fit))
print(best_long_fit$estimate)


short <- expanded_periods[expanded_periods <= best_cutoff]
long <- expanded_periods[expanded_periods > best_cutoff]

# x축 범위 설정
#x_range <- seq(min(expanded_periods), max(expanded_periods), length.out = 430)
x_range <- seq(0, 430, by = 1)

# short 부분의 PDF 계산
best_short_dist_extracted <- sub("_.*", "", best_short_dist_name)
if (best_short_dist_extracted == "gamma") {
  short_pdf <- dgamma(x_range, shape = best_short_fit$estimate["shape"], rate = best_short_fit$estimate["rate"])
} else if (best_short_dist_extracted == "lognorm") {
  short_pdf <- dlnorm(x_range, meanlog = best_short_fit$estimate["meanlog"], sdlog = best_short_fit$estimate["sdlog"])
} else if (best_short_dist_extracted == "weibull") {
  short_pdf <- dweibull(x_range, shape = best_short_fit$estimate["shape"], scale = best_short_fit$estimate["scale"])
}

# long 부분의 정규 분포 PDF 계산
#long_pdf <- dnorm(x_range, mean = best_long_fit$estimate["mean"], sd = best_long_fit$estimate["sd"])
best_long_dist_extracted <- sub(".*_", "", best_short_dist_name)
if (best_long_dist_extracted == "normal") {
  long_pdf <- dnorm(x_range, mean = best_long_fit$estimate["mean"], sd = best_long_fit$estimate["sd"])
} else if (best_long_dist_extracted == "gamma") {
  long_pdf <- dgamma(x_range, shape = best_long_fit$estimate["shape"], rate = best_long_fit$estimate["rate"])
} else if (best_long_dist_extracted == "lognorm") {
  long_pdf <- dlnorm(x_range, meanlog = best_long_fit$estimate["meanlog"], sdlog = best_long_fit$estimate["sdlog"])
}


# 데이터 프레임 생성
plot_data <- data.frame(
  x = rep(x_range, 2),
  y = c(short_pdf, long_pdf),
  type = rep(c("Short", "Long"), each = length(x_range))
)

# 히스토그램과 PDF 그래프 그리기
ggplot() +
  geom_histogram(aes(x = expanded_periods, y = ..density..), bins = 30, fill = "grey", alpha = 0.5) +
  geom_line(data = plot_data, aes(x = x, y = y, color = type), size = 1.2) +
  geom_vline(xintercept = best_cutoff, linetype = "dashed", color = "black", size = 1) +
  scale_color_manual(values = c("Short" = "blue", "Long" = "red")) +
  labs(title = "Best Fit Distributions for Short and Long Periods",
       x = "Incubation Period",
       y = "Density") +
  theme_minimal()

#############################################################
# uniform dist
xgrid<-0:9
uniform_pdf <- ifelse(xgrid >= 4 & xgrid <= 7, 1 / (7 - 4 + 1), 0)

conv_result <- convolve(short_pdf, uniform_pdf, type = "open")
plot(short_pdf, type='l', col='blue', lwd=2)
lines(conv_result, type='l', col='red',lwd=2)
legend("topright", legend=c("Short PDF", "Convolution Result"), col=c("blue", "red"), lwd=2)

conv_result2 <- convolve(long_pdf, uniform_pdf, type = "open")
plot(long_pdf, type='l', col='blue', lwd=2)
lines(conv_result2, type='l', col='red',lwd=2)
legend("topleft", legend=c("Long PDF", "Convolution Result"), col=c("blue", "red"), lwd=2)
