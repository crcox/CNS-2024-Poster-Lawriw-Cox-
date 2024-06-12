library(dplyr)
library(ggplot2)
library(purrr)
library(readr)
library(rlang)
library(tidyr)

load_many_as_list <- function(x, labels = NULL) {
    e <- rlang::env()
    walk(x, load, envir = e)
    lst <- as.list(e)
    if (!is.null(labels)) {
        names(lst) <- labels
    } else if (!is.null(names(x))) {
        names(lst) <- names(x)
    }
    return(lst)
}

d <- read_csv("data/Fitted_values_psycholing.csv") %>%
    rename(
        aoa_raw = aoa,
        wf_raw = wf,
        wl_raw = wl,
        cd_raw = cd
    ) %>%
    pivot_longer(
        all_of(c("aoa_raw", "wf_raw", "wl_raw", "cd_raw", "aoa_fitted", "wf_fitted", "wl_fitted", "cd_fitted")),
        names_to = c("metric", ".value"),
        names_pattern = "([a-z]+)_([a-z]+)"
    )


standard_error <- function(x, na.rm = FALSE) {
    if (na.rm) {
        x <- na.omit(x)
    }
    sd(x) / sqrt(length(x))
}

itemwise <- d %>%
    group_by(metric, R, G, C, Cue) %>%
    summarize(x = mean(raw, na.rm = TRUE)) %>%
    mutate(R = as.factor(R + 2), G = factor(G, levels = c(-0.5, 0.5), labels = c("ASD", "non-ASD")), C = factor(C, levels = c(-.5, .5), labels = c("adult", "child")))

avg <- itemwise %>%
    group_by(metric, R, G, C) %>%
    summarize(
        m = mean(x, na.rm = TRUE),
        s = sd(x, na.rm = TRUE),
        se = standard_error(x, na.rm = TRUE)
    )

ggplot(avg, aes(x = R, y = m, color = interaction(G,C))) +
    geom_pointrange(aes(ymin = m - se, ymax = m + se), position = position_dodge(width = .2)) +
    facet_wrap(~metric, scales = "free_y")



itemwise_diff <- itemwise %>%
    pivot_wider(id_cols = c("metric", "R", "G", "Cue"), names_from = "C", values_from = "x") %>%
    mutate(diff = adult - child) %>%
    group_by(metric, R, G) %>%
    summarize(
        m = mean(diff, na.rm = TRUE),
        s = sd(diff, na.rm = TRUE),
        se = standard_error(diff, na.rm = TRUE)
    )

ggplot(itemwise_diff, aes(x = R, y = m, color = G)) +
    geom_pointrange(aes(ymin = m - se, ymax = m + se)) +
    geom_hline(yintercept = 0, color = NA) +
    facet_wrap(~metric, scales = "free_y")


load('data/null-repsim-cor-ASD-cond.Rdata')
load('data/null-repsim-cor-TD-cond.Rdata')
load('data/null-repsim-cor-cond.Rdata')
load('data/null-repsim-cor-group.Rdata')
load('data/overall_cor_ASD_cond.Rdata')
load('data/overall_cor_TD_cond.Rdata')
load('data/overall_cor_stan_child.Rdata')

## Add in effects of groupXcond


## Rerun analysis on data collected from TD group on full set of CDI words but reduced to the words that we include here. Look at semantic similarity among categories.
x <- list(ASD_cond = overall_cor_ASD_cond, TD_cond = overall_cor_TD_cond, group = overall_cor_ASD_TD, cond = overall_cor_stan_child)
list_null <- list(ASD_cond = null_repsim_cor_ASD_cond, TD_cond = null_repsim_cor_TD_cond, group = null_repsim_cor_GROUP, cond = null_repsim_cor_COND )
z_score <- function(correlation, null){
  return((correlation - mean(null))/sd(null))
}


stats_func <- function(x,y){
  fisher_x <- atanh(x)
  fisher_y <- atanh(y)
  z <- z_score(fisher_x,fisher_y)
  p <- (sum(x>y)+1)/length(y)
  mu <- mean(fisher_y)
  sigma <- sd(fisher_y)
  shapiro <- shapiro.test(fisher_y)
  return(list(p = p, z_score = z, mu = mu, sigma = sigma, shapiro_test = shapiro))
}

o <- mapply(stats_func,x, list_null, USE.NAMES = TRUE)
names(o) <- rep(c("p", "z_score", "mu", "sigma", "shapiro_test"),4)
ASD_cond <- o[1:5]
TD_cond <- o[6:10]
group <- o[11:15]
cond <- o[16:20]


library(ggplot2)

hist(null_repsim_cor_ASD_cond,xlim = c(0.50,0.8),xlab = "Pearson's r")

