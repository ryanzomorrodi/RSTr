#
#
#
#
#
#
#
#
#
#
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
library(RSTr)
```
#
#
#
data_u <- lapply(miheart, \(x) x[, "75-84", "1988", drop = FALSE])
mod_ucar <- ucar("my_test_model", data_u, miadj, tempdir(), seed = 1234)
#
#
#
#
#
#
#
estimates <- get_estimates(mod_ucar)
estimates_supp <- estimates[estimates$relative_precision > 1 & estimates$events < 10, ]
plot(estimates$events, estimates$relative_precision, xlab = "Events", ylab = "Relative Precision")
points(estimates_supp$events, estimates_supp$relative_precision, col = "red")
abline(h = 1, col = "blue")
abline(v = 10, col = "blue")
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
mod_eucar <- eucar("my_test_model", data_u, miadj, tempdir(), seed = 1234, A = 6)
#
#
#
#
#
estimates_eucar <- get_estimates(mod_eucar)
plot(estimates_eucar$events, estimates_eucar$relative_precision, xlab = "Events", ylab = "Relative Precision", col = "purple")
points(estimates$events, estimates$relative_precision)
abline(h = 1, col = "blue")
abline(v = 10, col = "blue")

#
#
#
#
#
library(ggplot2)
ggplot(mishp) +
  geom_sf(aes(fill = estimates$medians)) +
  labs(
    title = "Spatially Smoothed Estimates, Unrestricted UCAR Model",
    fill = "Deaths per 100,000"
  ) +
  scale_fill_viridis_c() +
  theme_void()
ggplot(mishp) +
  geom_sf(aes(fill = estimates_eucar$medians)) +
  labs(
    title = "Spatially Smoothed Estimates, EUCAR Model",
    fill = "Deaths per 100,000"
  ) +
  scale_fill_viridis_c() +
  theme_void()
#
#
#
#
#
#
#
#
#
mod_eucar <- suppress_estimates(mod_eucar)
#
#
#
#
#
est_crude <- data_u$Y / data_u$n * 1e5
est_crude[data_u$Y < 10] = NA # Suppression criteria for CDC WONDER 
est_eucar <- get_estimates(mod_eucar)
ggplot(mishp) +
  geom_sf(aes(fill = est_crude)) +
  labs(
    title = "Crude Estimates",
    fill = "Deaths per 100,000"
  ) +
  scale_fill_viridis_c() +
  theme_void()
ggplot(mishp) +
  geom_sf(aes(fill = est_eucar$medians_suppressed)) +
  labs(
    title = "Spatially Smoothed Estimates, EUCAR Model",
    fill = "Deaths per 100,000"
  ) +
  scale_fill_viridis_c() +
  theme_void()
#
#
#
#
#
#
#
#
#
#
#
data_u <- lapply(miheart, \(x) x[, c("65-74", "75-84", "85+"), "1988", drop = FALSE])
A <- 6 * colSums(data_u$Y) / sum(data_u$Y)
mod_eucar <- eucar("test_eucar", data_u, miadj, tempdir(), seed = 1234, A = A)
#
#
#
#
#
std_pop <- c(68775, 34116, 9888)
mod_eucar <- age_standardize(mod_eucar, std_pop, "65up")
mod_eucar <- suppress_estimates(mod_eucar)
est_eucar <- get_estimates(mod_eucar)
ggplot(mishp) +
  geom_sf(aes(fill = est_eucar$medians_suppressed)) +
  labs(
    title = "Age-Standardized Spatially Smoothed Estimates, EUCAR Model",
    fill = "Deaths per 100,000"
  ) +
  scale_fill_viridis_c() +
  theme_void()
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
