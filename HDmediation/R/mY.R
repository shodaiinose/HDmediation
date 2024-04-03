mY <- function(data, npsem, family, folds, learners, ...) {
    q <- matrix(nrow = nrow(data), ncol = 2)
    colnames(q) <- c("Q(0,W)", "Q(1,W)")
    for (v in seq_along(folds)) {
        train <- origami::training(data, folds[[v]])
        try(train <- train[train[[npsem$S]] == 1, ], silent = TRUE)
        valid_1 <- valid_0 <- origami::validation(data, folds[[v]])
        valid_1[[npsem$A]] <- 1
        valid_0[[npsem$A]] <- 0
        
        if (!is.null(npsem$cens)) {
            obs <- train[, npsem$cens] == 1
        } else {
            obs <- rep(TRUE, nrow(train))
        }
        
        preds <- crossfit(train[obs, c(npsem$Y, npsem$W, npsem$A)],
                          list(valid_0, valid_1),
                          npsem$Y,
                          family,
                          learners = learners)
        
        q[folds[[v]]$validation_set, 1] <- preds[[1]]
        q[folds[[v]]$validation_set, 2] <- preds[[2]]
    }
    q
}
