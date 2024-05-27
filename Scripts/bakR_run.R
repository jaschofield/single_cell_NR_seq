### Half life analysis with bakR
library(bakR)
library(tidyverse)

## input processed cB file
cB_time <- read.csv('//path_to_cB/cB.csv')
cB_time <- cB_time %>% dplyr::select(-XF) %>% mutate(XF = GF)

time_metadf <- data.frame(c(0, 1, 1, 1),
                          c(1, 1, 1, 1))
colnames(time_metadf) <- c("tl", "Exp_ID")
rownames(time_metadf) <- c("JS2211_0_GEX_S1_R2", "JS2211_10_GEX_S3_R2", "JS2211_15_GEX_S4_R2", "JS2211_5_GEX_S2_R2")

bakRData_time <- bakRData(cB_time, time_metadf)

time_fit_out <- bakRFit(bakRData_time, StanRateEst = FALSE)

Fn_estimates <- time_fit_out$Fast_Fit$Fn_Estimates
ggplot(Fn_estimates, aes(x = logit_fn, y = logit_fn_se)) + geom_point()

## Calculate fnew from logit_fn
Fn_estimates <- Fn_estimates %>%
  mutate(fnew = 1/(1+exp(-logit_fn)))

## Perform basic filtering and reformat data
half_life_DF <- Fn_estimates %>% filter(logit_fn_se < 0.145) %>% dplyr::select(XF, fnew, sample) %>% pivot_wider(names_from = sample, values_from = fnew)

## Calculate half life from fnew
half_life_DF <- half_life_DF %>% mutate(HL_5min = 5*log(2)/-log(1-JS2211_5_GEX_S2_R2),
                                        HL_10min = 10*log(2)/-log(1-JS2211_10_GEX_S3_R2),
                                        HL_15min = 15*log(2)/-log(1-JS2211_15_GEX_S4_R2))



ggplot(half_life_DF, aes(x = HL_10min)) + geom_histogram(bins = 100) + theme_bw() + coord_cartesian(xlim = c(0, 100))
