# data is a numeric matrix, rows are samples, columns are features
# label is a 'cancer' or 'healthy' string vector. Its length must be equal to nrow of data
# return a numberic probility vector, 1 ~ cancer, 0 ~ healthy
#        a confusion matrix
#        the model itself


# need randomForest
# https://bit.ly/3fIJsxd
wrapper.randomForest <- function(data, label, data.test = NULL, mtry = NULL, ntree = 500) {
   result = list()

   if (nrow(data) != length(label)) {
      message(sprintf('The row number of %s, %s', nrow(data), length(label)))
      return(NULL)
   }

   if (any(is.na(label))) {
      message(sprintf('The row number of %s, %s', nrow(data), length(label)))
      return(NULL)
   } else {
      label.fa = factor(label, levels = c('healthy', 'cancer'))
   }

   if (any(is.na(data))) {
      data.ma = randomForest::rfImpute(x = data, y = label.fa, iter = 6)[, -1]
   } else {
      data.ma = data
   }


   # train the model
   if (!is.null(mtry)) {
      model = randomForest::randomForest(x = data.ma, y = label.fa, mtry = mtry, proximity = T, ntree = ntree)
   } else {
      model = randomForest::randomForest(x = data.ma, y = label.fa, proximity = T, ntree = ntree)
   }
   print(model)
   result$model = model

   err.ma = model$err.rate
   error.rate.df = data.frame(trees = rep(1:nrow(err.ma), ncol(err.ma)),
                              type = rep(colnames(err.ma), each = nrow(err.ma)),
                              error = as.vector(err.ma)
                              )

   plot = ggplot(data = error.rate.df) +
      geom_line(mapping = aes(x = trees, y = error, color = type)) +
      theme_classic() +
      theme(legend.position = 'bottom') +
      ylab('Error rate') +
      xlab('Trees number')
   print(plot)


   # predict the train data or test data if not null
   if (is.null(data.test)) {
      p <- predict(model, data.ma, type = 'vote',  norm.votes = T)
      predict.c = apply(p, MARGIN = 1, FUN = function(v) {names(v)[match(max(v), v)]})
      predict.fa = factor(predict.c, levels = c('healthy', 'cancer'))
      result$confusion = caret::confusionMatrix(predict.fa, reference = label.fa)
   } else {
      p <- predict(model, data.test, type = 'vote',  norm.votes = T)
      result$confusion = model$confusion
   }

   result$prob = p

   return(result)

}

# need gbm
# https://datascienceplus.com/gradient-boosting-in-r/
wrapper.gbm <- function(data, label, data.test = NULL) {

   result = list()

   label.c = ifelse(label == 'cancer', yes = 1, no = 0)

   model = gbm::gbm(label ~ .,
                    data = data.frame(data, label = label.c),
                    distribution = 'bernoulli',
                    n.trees = 10000,
                    shrinkage = 0.01,
                    interaction.depth = 4,
                    cv.folds = 4)
   result$model = model

   if (is.null(data.test)) {
      predict.c = predict(model, data.frame(data), type = 'response', n.trees = 10000)
      predict.fa = factor(ifelse(predict.c < 0.5, 'healthy', 'cancer'), levels = c('healthy', 'cancer'))
      result$confuse = caret::confusionMatrix(predict.fa, factor(label, levels = c('healthy', 'cancer')))
   } else {
      predict.c = predict(model, data.frame(data.test), type = 'response', n.trees = 10000)
   }

   result$prob = predict.c

   return(result)
}


# R in action third edition 13.2
wrapper.glm <- function(data, label, data.test = NULL) {
   result = list()

   label.c = ifelse(label == 'cancer', yes = 1, no = 0)

   family = binomial()
   model <- glm(label ~ ., data = data.frame(data, label = label.c), family = binomial())
   model.od <- glm(label ~ ., data = data.frame(data, label = label.c), family = quasibinomial())

   # test if overdispension, method 1
   p = pchisq(summary(model.od)$dispersion * model$df.residual, model$df.residual, lower = F)
   if (p <= 0.05) {
      model = model.od
      family = quasibinomial()
   } else { # test if overdispension, method 2
      p = model$deviance / model$df.residual
      if (p >= 2) {
         model = model.od
         family = quasibinomial()
         }
   }

   # extract significant coefficients
   coef.ma = summary(model)$coefficients[-1,]
   is.sign = coef.ma[, 4] <= 0.1
   if (any(is.sign)) {
      model.reduce <- glm(label ~ ., data = data.frame(data[, is.sign], label = label.c), family = family)
   } else {
      message() # no significant variables
   }
   print(summary(model.reduce))
   result$model = model.reduce

   plot(predict(model, type = 'response'), residuals(model, type = 'deviance'), xlab = 'Fitted values', ylab = 'residuals', main = 'Residuals vs. Fitted values')

   # predict
   if (!is.null(data.test)) {
      predict.c = predict(model.reduce, data.frame(data.test), type = "response")
   } else {
      predict.c = model$fitted.values
   }
   result$prob = predict.c
   predict.c = factor(ifelse(predict.c >= 0.5, yes = 'cancer', no = 'healthy'), levels = c('healthy', 'cancer'))
   result$confuse = caret::confusionMatrix(predict.c, factor(label, levels = c('healthy', 'cancer')))

   return(result)
}


# https://www.hackerearth.com/practice/machine-learning/machine-learning-algorithms/beginners-tutorial-on-xgboost-parameter-tuning-r/tutorial/
wrapper.xgboost <- function(data, label, data.test = NULL) {

   label.c = ifelse(label == 'cancer', yes = 1, no = 0)
   data.xgm = xgb.DMatrix(data = data, label = label.c)


   params <- list(booster = "gbtree",
                  objective = "binary:logistic",
                  eta = 0.3,
                  gamma = 0,
                  max_depth = 6,
                  min_child_weight = 1,
                  subsample = 1,
                  colsample_bytree = 1)

   model <- xgb.cv(params = params,
                   data = data.xgm,
                   nrounds = 100,
                   nfold = 4,
                   showsd = T,
                   stratified = T,
                   print_every_n = 10,
                   early_stop_round = 20,
                   maximize = F)

   min(model$test.error.mean)

}


# need h2o
# data is a numeric matrix, rows are samples, columns are features with row and column names
# label is a 'cancer' or 'healthy' string vector. Its length must be equal to nrow of data
# return a numberic probility vector, 1 ~ cancer, 0 ~ healthy
wrapper.stacked.ensemble <- function(data, label, data.test = NULL) {
   result = list()

   data = affairs.ma
   label = affairs.label

   if (nrow(data) != length(label)) {
      message(sprintf('The row number of %s, %s', nrow(data), length(label)))
      return(NULL)
   }

   if (any(is.na(label))) {
      message(sprintf('The row number of %s, %s', nrow(data), length(label)))
      return(NULL)
   } else {
      label.fa = factor(label, levels = c('healthy', 'cancer'))
   }

   if (any(is.na(data))) {
      message('Some items in data is NA')
      data.ma = randomForest::rfImpute(x = data, y = label.fa, iter = 6)
      data.df = data.frame(data.ma)
      data.df$y = factor(data.df$y)
   } else {
      data.df = data.frame(y = label.fa, data)
   }


   train = as.h2o(data.df)
   y <- "y"
   x <- setdiff(colnames(data.df), y)
   nfolds <- 4

   my_rf <- h2o.randomForest(x = x,
                             y = y,
                             training_frame = train,
                             ntrees = 50,
                             nfolds = nfolds,
                             keep_cross_validation_predictions = TRUE,
                             )

   my_xgb <- h2o.xgboost(x = x,
                         y = y,
                         training_frame = train,
                         ntrees = 50,
                         nfolds = nfolds,
                         keep_cross_validation_predictions = TRUE,
                         booster = "dart",
                         normalize_type = "tree"
                         )

}

# need h2o
wrapper.automl <- function(data, label, data.test = NULL) {
   h2o.init()

   result = list()

   #data = affairs.ma # for testing
   #label = affairs.label

   if (nrow(data) != length(label)) {
      message(sprintf('The row number of %s, %s', nrow(data), length(label)))
      return(NULL)
   }

   if (any(is.na(label))) {
      message(sprintf('The row number of %s, %s', nrow(data), length(label)))
      return(NULL)
   } else {
      label.fa = factor(label, levels = c('healthy', 'cancer'))
   }

   if (any(is.na(data))) {
      message('Some items in data is NA')
      data.ma = randomForest::rfImpute(x = data, y = label.fa, iter = 6)
      data.df = data.frame(data.ma)
      data.df$y = factor(data.df$y)
   } else {
      data.df = data.frame(y = label.fa, data)
   }

   train = as.h2o(data.df)
   y <- "y"
   x <- setdiff(colnames(data.df), y)

   aml <- h2o.automl(x = x,
                     y = y,
                     training_frame = train,
                     nfolds = 4,
                     #max_runtime_secs = 30,
                     keep_cross_validation_predictions = T)

   m <- h2o.get_best_model(aml)


   result$model = m
   result$confuse = h2o.confusionMatrix(m)
   result$prob = h2o.predict(m, newdata = train)

   return(result)

}


