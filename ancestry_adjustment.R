library(dplyr)
library(stringr)

prep_scores <- function(scores) {
    scores2 <- scores %>%
        select(IID=`IID`, ends_with("_SUM")) %>%
        select(IID, starts_with("PGS")) %>%
        mutate(IID = as.character(IID))
    names(scores2) <- str_replace(names(scores2), "_SUM", "")
    return(scores2)
}


fit_prs <- function(scores, pcs) {
    all_scores <- names(scores)[-1]
    mean_coef <- list()
    var_coef <- list()
    for (s in all_scores) {
        dat <- scores %>%
            select(IID, !!s) %>%
            inner_join(pcs, by="IID")
        pccols <- names(pcs)[str_detect(names(pcs), "^PC")]
        model_string <- paste(s, "~", paste(pccols, collapse="+"))
        pcmod <- lm(model_string, data=dat)
        mean_coef[[s]] <- pcmod$coefficients

        dat$residsq <- (pcmod$residuals)^2
        model_string <- paste("residsq ~", paste(pccols, collapse="+"))
        # use Gamma regression with log link to avoid negative values
        pcmod2 <- tryCatch({
                glm(model_string, family=Gamma(link = "log"), data=dat, control=list(maxit=1000))
            }, error = function(e) {
                #coefs <- rep(as.integer(NA), length(pcmod$coefficients))
                coefs <- c(log(var(pcmod$residuals)), rep(0, length(pccols)))
                return(list(coefficients=setNames(coefs, names(pcmod$coefficients))))
            }
        )
        var_coef[[s]] <- pcmod2$coefficients
    }
    prep_output <- function(x) {
        x <- bind_rows(x)
        names(x)[1] <- "Intercept"
        x <- bind_cols(tibble(score=all_scores), x)
    }
    out <- list(
        mean_coef = prep_output(mean_coef),
        var_coef = prep_output(var_coef)
    )
    return(out)
}


adjust_prs <- function(scores, pcs, mean_coef, var_coef) {
    samples <- intersect(scores$IID, pcs$IID)
    scores <- scores[match(samples, scores$IID),]
    pcs <- pcs[match(samples, pcs$IID),]
    pcmat <- pcs %>%
        select(starts_with("PC")) %>%
        as.matrix()

    all_scores <- names(scores)[-1]
    score_adj <- list()
    for (s in all_scores) {
        score <- scores %>%
            select(!!s) %>%
            unlist()
        mean_intercept <- mean_coef %>%
            filter(score == s) %>%
            select(Intercept) %>%
            unlist()
        mean_pcs <- mean_coef %>%
            filter(score == s) %>%
            select(starts_with("PC")) %>%
            unlist()
        var_intercept <- var_coef %>%
            filter(score == s) %>%
            select(Intercept) %>%
            unlist()
        var_pcs <- var_coef %>%
            filter(score == s) %>%
            select(starts_with("PC")) %>%
            unlist()
        score_adj[[s]] <- (score - (mean_intercept + as.vector(pcmat %*% mean_pcs))) / sqrt(exp(var_intercept + as.vector(pcmat %*% var_pcs)))
    }
    score_adj <- bind_cols(score_adj)
    score_adj <- bind_cols(tibble(IID=samples, score_adj))
    return(score_adj)
}


prep_groups <- function(groups) {
    groups2 <- groups[,1:2]
    names(groups2) <- c("IID", "group")
    groups2$IID <- as.character(groups2$IID)
    groups2$group <- as.character(groups2$group)
    return(groups2)
}


adjust_prs_empirical <- function(scores, groups) {
    dat <- inner_join(groups, scores)
    samples <- dat$IID

    all_scores <- names(scores)[-1]
    grp_list <- list()
    for (g in unique(dat$group)) {
        z_list <- list()
        tmp <- filter(dat, group == g)
        for (s in all_scores) {
            score <- tmp %>%
                select(!!s) %>%
                unlist()
            z_list[[s]] = (score - mean(score, na.rm=TRUE)) / sd(score, na.rm=TRUE)
        }
        z_list <- bind_cols(z_list)
        grp_list[[g]] <- bind_cols(tibble(IID=tmp$IID, group=g, z_list))
    }
    return(bind_rows(grp_list))
}
