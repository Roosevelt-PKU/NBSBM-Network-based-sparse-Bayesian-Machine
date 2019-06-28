#' Classification with SMVs and prior knowledge
#' 
#' pathClass is a collection of classification methods that
#' use information about how features are connected in the underlying
#' biological network as an additional source of information. This
#' additional knowledge is incorporated into the classification a
#' priori. Several authors have shown that this approach significantly
#' increases the classification performance.
#' 
#' @name pathClass-package
#' @aliases pathClass
#' @title Classification with SMVs and prior knowledge
#' @author Marc Johannes \email{JohannesMarc@@gmail.com}
#' @keywords package
#' @docType package
NULL


#' Performs cross-validation with a specified algorithm
#'
#' Performs a cross-validation using the specified algorithms.
#' If package parallel is loaded the cross-validation will be
#' performed in parallel. If the parallel package is loaded but a parallel
#' cross-validation is not wanted \code{parallel} can be set to \code{FALSE}.
#' If parallel cross-validation is desired the number of cores can be choosen by
#' using the \code{cores} parameter.
#'
#' @param x a p x n matrix of expression measurements with p samples and n genes.
#' @param y a factor of length p comprising the class labels.
#' @param theta.fit the method to learn a decision boundary. Currently available are \code{\link{fit.rrfe}}, \code{\link{fit.rfe}}, \code{\link{fit.graph.svm}}, \code{\link{fit.networkBasedSVM}}
#' @param folds number of folds to perform
#' @param repeats number of how often to repeat the x-fold cross-validation
#' @param parallel should the cross-validation be performed in parallel
#' i.e. on several cpu-cores. (see also \code{Details} section)
#' @param cores specify the number of cores that should be used for parallel cross-validation.
#' @param DEBUG should debugging information be plotted.
#' Defaults to n - 1 cores.
#' @param ... additional parameters to theta fit.
#' @return a list with the results of the cross-validation. See details for more information.
#' @export
#' @callGraphPrimitives
#' @note Parallel cross-validation can only be performed if the parallel-package
#' was loaded prior to calling this function.
#' @seealso \code{\link{fit.rrfe}}, \code{\link{fit.rfe}}, \code{\link{fit.graph.svm}}, \code{\link{fit.networkBasedSVM}}
#' @author Marc Johannes \email{JohannesMarc@@gmail.com}
#' @examples
#' set.seed(4321)
#' data(example_data)
#' res.rfe  <- crossval(x, y, DEBUG=TRUE, theta.fit=fit.rfe, folds=2, repeats=1, parallel=TRUE, Cs=10^(-3:3))
#' res.rrfe <- crossval(x, y, DEBUG=TRUE, theta.fit=fit.rrfe, folds=3, repeats=1, parallel=TRUE, Cs=10^(-3:3), mapping=mapping, Gsub=adjacency.matrix, d=1/2)
crossval.nbsbm <- function(x, y, folds=10, repeats=1, parallel = TRUE, cores = NULL, DEBUG=FALSE, ...){

  pp <- ("package:parallel" %in% search())
  
  if(pp == TRUE && parallel == TRUE){
    
    if(is.null(cores)) cores <- parallel:::detectCores()
    options(cores = cores - 1)
    
    cat("Detected ", cores," cores. Will use ", getOption("cores"), " of them.\n")
    parallel <- TRUE
  }
  else{

    if(parallel == TRUE) cat('Package \'parallel\' not loaded. Please, load it manually prior to calling this function if you want to run classification in parallel.\n',sep='')
    cat('Will continue with sequential crossvalidation.\n', sep='')
    parallel <- FALSE
  }

  if(!is.factor(y)) stop("y must be factor!\n")
  if(length(levels(y)) != 2) stop('y must be factor with 2 levels.\n')
  if(length(y) != nrow(x)) stop('y must have same length as nrow(x).\n')
  
  n     <- length(y)
  folds <- trunc(folds)
  print(folds)
  if (folds < 2) stop("folds should be greater than or equal to 2.\n")
  if (folds > n) stop("folds should be less than or equal to the number of observations.\n")

  cuts  <- cv.repeats <- list()

  print(cuts)
  
  for(r in 1:repeats){

    perm <- sample(1:n)
    repeat.models <- NULL
    
    for(k in 1:folds){
      tst <- perm[seq(k, n, by=folds)]
      trn <- setdiff(1:n, tst)      
      cuts[[k]] <- list(trn=trn, tst=tst)
    }

    pb <- txtProgressBar(min = 0, max = folds, style = 3)

    #print(pb)
    
    if(DEBUG) cat('Starting classification of repeat:',r,'\n')
    print(cuts)
    
    if(parallel)  repeat.models <- mclapply(1:folds, classify, cuts=cuts, pb=pb, x=x, y=y, theta.fit= theta.fit, cv.repeat=r, DEBUG=DEBUG, ...)
    else          repeat.models <-   lapply(1:folds, classify, cuts=cuts, pb=pb, x=x, y=y, theta.fit= theta.fit, cv.repeat=r, DEBUG=DEBUG, ...)

    close(pb)
    
    if(length(repeat.models) != folds){
      geterrmessage()
      stop("One or more processes did not return. May be due to lack of memory.\n")
    }
    if(DEBUG) cat('All models of repeat:',r,'have been trained.\n')
    cv.repeats[[r]] <- repeat.models
  }

  ## summarize results
  cv <- sapply(cv.repeats, function(cv.repeat) rowSums(sapply(cv.repeat, function(model) model$cv)))
  colnames(cv) <- paste("Repeat",1:repeats,sep="")

  auc <- sapply(cv.repeats, function(cv.repeat) sapply(cv.repeat, function(model) model$auc))
  colnames(auc) <- paste("Repeat",1:repeats,sep="")
  rownames(auc) <- paste("Fold",1:folds,sep="")

  fits <- lapply(cv.repeats,function(cv.repeat) lapply(cv.repeat, function(model) model$model))
  names(fits) <- paste("Repeat",1:repeats,sep="")
  fits <- lapply(fits, function(x)  {names(x) = paste("Fold", 1:folds, sep = ""); x })

  res <- list(cv=cv, auc=auc, fits=fits, labels=y)
  class(res) <- 'pathClassResult'
  return(res)
  
  ## return(list(cv=cv, auc=auc, fits=fits, labels=y))
}

classify <- function(fold, cuts, pb, x, y, theta.fit, cv.repeat, DEBUG=FALSE, ...){
  gc()
  if(DEBUG) cat('starting Fold:',fold,'\n')
  
  ## get training and test indices
  trn <- cuts[[fold]]$trn
  tst <- cuts[[fold]]$tst

  #print(trn)
  #print(tst)

  #print(fold)
  #print(pb)
  #print(theta.fit)
  #print(cv.repeat)
  #print("trainres")
  #print(x[trn,,drop=FALSE])
  #print(y[trn])
  print(DEBUG)
  #print(trained)
  ## train and test the model

  #gridBeta <- exp(seq(-10,-9,length.out = 10))

  ##gridBeta <- exp(seq(-20,2,length.out = 40))

  gridBeta <- exp(seq(-10,2,length.out = len))

  print(gridBeta)
  ret <- list()
  error <- matrix(0, length(gridBeta))
  
  trainData <- x[trn,,drop=FALSE]
  testData <-  x[tst,,drop=FALSE]
  m <- apply(trainData, 2, mean)
  s <- apply(trainData, 2, sd)
 
  scaledTrain <- (trainData - matrix(m, nrow = nrow(trainData), ncol = ncol(trainData), byrow = T))/matrix(s, nrow = nrow(trainData), ncol = ncol(trainData), byrow = T)
  scaledTest <- (testData - matrix(m, nrow = nrow(testData), ncol = ncol(testData), byrow = T))/matrix(s, nrow = nrow(testData), ncol = ncol(testData), byrow = T)
  scaledTrain <- cbind(scaledTrain, rep(1, nrow(scaledTrain)))
  scaledTest <- cbind(scaledTest, rep(1, nrow(scaledTest)))
  #print("scaledtrain")
  #print(scaledTrain)

  for (j in 1 : length(gridBeta)) {

        #resulttry <- try(ret[[ j ]] <- nbsbc_fast(x[trn,,drop=FALSE], y[trn], network, alpha = 0, beta = gridBeta[ j ]))
        resulttry <- try(ret[[ j ]] <- nbsbc_fast(scaledTrain, y[trn], network, alpha = 0, beta = gridBeta[ j ], b0 = b))
        if (class(resulttry) == "try-error")

	         ret[[ j ]] <- list(error = 1)
		
		   error[ j ]  <- ret[[ j ]]$error
		
		   cat(".", gridBeta[ j ], ret[[ j ]]$error, "\n")
	}

	cat(gridBeta[ which.min(error) ], ret[[ which.min(error) ]]$error, "\n")
	ret <- ret[[ which.min(error) ]]
      w <- ret$m; 
      #pred <- sign(x[trn,,drop=FALSE] %*% w); 

      # pred <- x[tst,,drop=FALSE] %*% w ;
      pred <- scaledTest %*% w;
  
      #trained <- theta.fit(x[trn,,drop=FALSE], y[trn], DEBUG=DEBUG, ...)
  
      #test <- predict(object=trained, newdata=x[tst,,drop=FALSE])
  
      #print("tested_restults")
      #print(test)

  ## save the test indices
  #trained[["tst.indices"]] <- tst
  
  ## calculate the AUC
  #label <- y[tst] 

  label <- sign(as.numeric(y[tst]) - 1.5) # because y is a factor
  print("test_label")
  print(label)
  print("expr_label")

  test <- pred

  print(pred)

  auc   <- calc.auc(test, label)

  print("auc")
  print(auc)
  #if(DEBUG) 
  cat(" Test AUC =", auc, "\n")

  cv      <- double(length=length(y))
  cv[tst] <- test

  if(DEBUG) {cat('Finished fold:',fold,'\n\n')}
  else {setTxtProgressBar(pb, getTxtProgressBar(pb) + 1)}

  gc()
  #list(fold=fold, model=trained, auc=auc, cv=cv)
  list(fold=fold, model=ret,auc=auc,cv=cv)
}


calc.auc <- function(prob,labels)
{
  ## this corrects a bug in ROCR:
  ## if all labels are from one group and there
  ## is no missclassification, ROCR is not able
  ## to calculate the auc
  ## patch => add a artificial prediction with prob = 0
  if(length(unique(labels)) == 1)
  {
    if(sign(labels[1]) == -1)
      labels <- c(labels,1)
    else
      labels <- c(labels,-1)
    prob <- c(prob,0)
  }
  print("prob")
  print(prob)
  print("labels")
  print(labels)
  #pred <- prediction(prob, labels)
  pred <- prediction(prob, labels, label.ordering=c(-1,1))
  print(pred)
  unlist(performance(pred, "auc")@y.values)
}

#' Prints the result of one or more cross-validation run(s)
#'
#' This function creates boxplots of the distribution of AUC for each reapeat of the cross-validation.
#' In a second plot the ROC curve of the AUCs is shown. If your result contains more than one cross-validation
#' result these are plotted one after the other.
#'
#' @param x A result of \code{crossval}.
#' @param label the main label of the plots.
#' @param toFile Should the results plotted into PDF file(s). If your result contains more than one cross-validation
#'               one PDF file is created for each result.
#' @param fname the name of the file to save the results in.
#' @param switchLabels If your AUC is below 0.5 you can switch the labels to get an AUC above 0.5.
#' @param avg the method for averaging the AUCs of several repeats. See \code{'\linkS4class{performance}'} for more information.
#' @param spread.estimate method to show the variation around the average of the ROC curve. See \code{'\linkS4class{performance}'} for more information.
#' @param ... currently ignored.
#' @export
#' @author Marc Johannes \email{JohannesMarc@@gmail.com}
#' @examples
#' \dontrun{
#' library(Biobase)
#' data(sample.ExpressionSet)
#' x <- t(exprs(sample.ExpressionSet))
#' y <- factor(pData(sample.ExpressionSet)$sex)
#' res.rfe <- crossval(x,y,DEBUG=TRUE,theta.fit=fit.rfe,folds=2,repeats=1,parallel=TRUE,Cs=10^(-3:3))
#' plot(res.rfe, toFile=FALSE)
#' }
plot.pathClassResult <- function(x, label='', toFile=TRUE, fname='Result', switchLabels=FALSE, avg="vertical", spread.estimate="boxplot",...){

  run <- label

  if(toFile){
    cat("Creating file ",fname,".pdf ...",sep="")
    pdf(file = paste(fname,".pdf",sep=""))
  }

  ## boxplot of test auc
  main <- paste(run)
  boxplot(x$auc,ylim=c(0,1),main=main, outline=FALSE)
  stripchart(as.data.frame(x$auc), method="jitter",jitter=0.05,add=T,pch=20,vertical=TRUE)

  ## open a second plotting device
  if(!toFile) dev.new()

  ## ROC Curve
  repeats <- ncol(x$cv)
  y.num <- sign(as.numeric(x$labels) - 1.5)
  if(switchLabels == T)
    y.num <- y.num * -1
  
  labels <- matrix(rep(y.num,repeats),ncol=repeats)


  labelswith =F

  ###auc.value.test=unlist(performance(prediction(x$cv,labels,label.ordering=c(-1,1)), "auc")@y.values)
  ###auc.value=unlist(performance(prediction(x$cv,labels,label.ordering=c(1,-1)), "auc")@y.values)

  auc.value.test=unlist(performance(prediction(x$cv,labels), "auc")@y.values)

  ###print(auc.values)
  ###print(auc.value)
  labels_reverse=NULL
  for(i in 1:repeats){
  ##if(length(which( auc.value.test[i]<0.5))>=1){
  if(auc.value.test[i]<0.5){
  labels_reverse <- cbind(labels_reverse, i)
  labelswith=T
  }else{
  #labels <- cbind(labels, y.num)
  }
  } 
  for (indexlable in labels_reverse){
  labels[,indexlable]=-labels[,indexlable]
  }

  ##labels <- matrix(rep(y.num,repeats),ncol=repeats)
  print(labelswith)
  
####if(false){##if(labelswith==F){

####  pred <- prediction(x$cv,labels,label.ordering=c(-1,1))

####  perf <- performance(pred,measure="tpr",x.measure="fpr")
  
  #auc <- mean(unlist(performance(prediction(x$cv,labels,label.ordering=c(-1,1)), "auc")@y.values))

####  auc <- mean(unlist(performance(prediction(x$cv,labels), "auc")@y.values))

  #print(auc.test)
  #print(labels)
  #print(auc)
  #print(unlist(performance(prediction(x$cv,labels), "auc")@y.values))
  #print(unlist(performance(prediction(x$cv,labels,label.ordering=c(1,-1)), "auc")@y.values))
####  }else{

####  pred <- prediction(x$cv,labels,label.ordering=c(1,-1))

####  perf <- performance(pred,measure="tpr",x.measure="fpr")
  
####  auc <- mean(unlist(performance(prediction(x$cv,labels,label.ordering=c(1,-1)), "auc")@y.values))

  ##auc <- mean(unlist(performance(prediction(x$cv,labels), "auc")@y.values))

####}
####}

  pred <- prediction(x$cv,labels)
  perf <- performance(pred,measure="tpr",x.measure="fpr")
  auc <- mean(unlist(performance(prediction(x$cv,labels), "auc")@y.values))

  main <- paste(run,"\nAUC = ",auc,sep="" )

  if(repeats > 1){
    plot(perf, avg=avg, spread.estimate=spread.estimate, main=main)
  }else{
    plot(perf,main=main)
  }

  if(toFile){
    dev.off()
  }

  cat("done\n",sep="")

}

#' Extracts features which have been choosen by the classifier(s).
#'
#' This function extracts the features which have been selected by the classifiers
#' during the cross-validation along with the number of times they have been choosen.
#' When, for example, performing a 5 times repeated 10-fold cross-validation the maximum
#' number a feature can be choosen is 50.
#'
#' @param res A result of \code{crossval}.
#' @param toFile Should the results be printed into a CSV-file.
#' @param fName the name of the file to save the results in.
#' @return a \code{data.frame} indicating the number of times a feature has been choosen.
#' @export
#' @author Marc Johannes \email{JohannesMarc@@gmail.com}
#' @examples
#' \dontrun{
#' library(Biobase)
#' data(sample.ExpressionSet)
#' x <- t(exprs(sample.ExpressionSet))
#' y <- factor(pData(sample.ExpressionSet)$sex)
#' res.rfe <- crossval(x,y,DEBUG=TRUE,theta.fit=fit.rfe,folds=2,repeats=1,parallel=TRUE,Cs=10^(-3:3))
#' extractFeatures(res.rfe, toFile=FALSE)
#' }
extractFeatures <- function(res, toFile=FALSE, fName='ClassificationFeatures.csv'){
  if(class(res) != 'pathClassResult') stop('\'res\' must be of class \'pathClassResult\'')
  
  ttt <- sort(table(unlist(lapply(res$fits, function(cv.repeat) lapply(cv.repeat, function(fold) fold$features)))), decreasing=TRUE)
  ttt <- data.frame(time.choosen=ttt)
  if(toFile){
    write.csv(ttt, file=fName)
    cat('Created file: ',getwd(),fName,'\n', sep='')
    invisible(ttt)
  }
  else
    return(ttt)
}


## INPUT
## CV       = a list of results of crossval.parallel
## res1     = result of first algorithm (corresponds to a index of CV)
## res2     = result of second algorithm (corresponds to a index of CV)
## alt      = alternative for wilcoxon test
## toFile   = should the result printed to a file
## filename = optional filename if not provided filename will be
##            res1vsres2
compare.auc <- function(res1, res2, alt="less", toFile=T, filename="", showPlots=T, switchLabels = F, name1 = "", name2 = ""){
  library(ROCR)
  
  repeats1 <- ncol(res1$cv)
  y.num1 <- sign(as.numeric(res1$labels) - 1.5)
  if(switchLabels == T){y.num1 = y.num1 * -1}
  labels1 <- matrix(rep(y.num1, repeats1), ncol=repeats1)
  pred1 <- prediction(res1$cv, labels1)
  AUCs1 <- unlist(performance(pred1, "auc")@y.values)
  perf1 <- performance(pred1, measure="tpr",x.measure="fpr")
 
  ## shifting auc values
  print("original auc list")
  print(unlist(performance(prediction(res1$cv,labels1), "auc")@y.values  )) 
  need_reverse_index_list = c()
  index = 1
  for(auc in unlist(performance(prediction(res1$cv,labels1), "auc")@y.values)){
    if (auc < 0.5){
      need_reverse_index_list = append(need_reverse_index_list, index)
    }
    index = index + 1
  }
  # print(need_reverse_index_list)
  for(index in need_reverse_index_list){
    labels1[,index] = -labels1[,index]
  }
  # print(labels)
  pred1 <- prediction(res1$cv,labels1)
  AUCs1 <- unlist(performance(pred1, "auc")@y.values)
  perf1 <- performance(pred1,measure="tpr",x.measure="fpr")
  print("final auc list")
  print(unlist(performance(prediction(res1$cv,labels1), "auc")@y.values  ))
  
  
  repeats2 <- ncol(res2$cv)
  y.num2 <- sign(as.numeric(res2$labels) - 1.5)
  if(switchLabels == T){y.num2 = y.num2 * -1}
  labels2 <- matrix(rep(y.num2, repeats2), ncol=repeats2)
  pred2 <- prediction(res2$cv, labels2)
  AUCs2 <- unlist(performance(pred2, "auc")@y.values)
  perf2 <- performance(pred2, measure="tpr",x.measure="fpr")
  
  pval <- round(wilcox.test(AUCs1, AUCs2, alternative=alt)$p.value, 5)
  print("pval")
  print(pval)
  
  ## turn alternative around if
  ## pVal is to big
  if(pval > 0.5){
    if(alt == "less") alt <- "greater"
    else alt <- "less"
  }
  pval <- round(wilcox.test(AUCs1, AUCs2, alternative=alt)$p.value, 5)

  main <- ""
  if(alt == "less")
    main <- paste("Wilcoxon p-Value:\n", name1, "(black line) < ", name2 ,"(blue line)","\n= ", pval, sep="")
  else
    main <- paste("Wilcoxon p-Value:\n" , name1, "(black line) > ", name2 ,"(blue line)","\n= ", pval, sep="")
  
  if(filename == "")
    filename <- paste("Comparison_", "nbsvm VS rfe", ".pdf", sep="")
  
  if(toFile){
    pdf(file=filename)
    cat(paste("Creating file ", filename, " ...", sep=""))
  }
  else if(showPlots){
    par(mfrow=c(2,1))
  }     
  
  if(showPlots | toFile){
    plot(perf1, avg="vertical", spread.estimate="boxplot", main=main)
    plot(perf2, avg="vertical", spread.estimate="boxplot", main=main, col="blue", add=T)
    boxplot(AUCs1,AUCs2, main=main)
  }
  if(toFile){
    dev.off()
    cat(" done\n")
  }
}
