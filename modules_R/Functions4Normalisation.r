###############################################################################
# Library of functions used in Normalisation 
# Author: Vadim Farztdinov
# Tested on R version 3.2.3 (2015-12-10), nickname Wooden Christmas-Tree
# ...
# Tested on R version 4.3.3 (2024-02-29), nickname Angel Food Cake  
###############################################################################
# Most of functions use object Data/DataObj, which is a list with three slots:
# 	Data$Matrix
# 	Data$MetaData
# 	Data$PepDesc
# DataDesc is a list with descriptions of factors in MetaData
# WorkDesc is a list with descriptions of the steps to be applied
###############################################################################
createDataObject<- function(DataMatrix,DataDesc){
  # create DataObject from DataMatrix in MeV format
  # Author: Vadim Farztdinov
  # Version: 1.0.0, 09-Sep-2013
  # Version: 1.0.1, 18-Feb-2014. Updated for case with N2=2
  
  # row.names(DataMatrix)<-DataMatrix[,DataDesc$col$peptideID];
  # Data$Matrix     <- as.data.frame(Data$Matrix)
  
  Data<-vector("list");
  if (DataDesc$fmt[2]>1) {
    Data$MetaData   <- DataMatrix[ (1:(DataDesc$fmt[1]-1)),-(1:(DataDesc$fmt[2]-1))];
    Data$Matrix     <- DataMatrix[-(1:(DataDesc$fmt[1]-1)),-(1:(DataDesc$fmt[2]-1))];
    Data$PepDesc    <- DataMatrix[-(1:(DataDesc$fmt[1]-1)), (1:(DataDesc$fmt[2]-1))];
    if (DataDesc$fmt[2]==2) {
      Data$PepDesc<-as.matrix(Data$PepDesc,nrow=length(Data$PepDesc));
      Data$PepDesc<-cbind(rownames(Data$Matrix),Data$PepDesc)
      colnames(Data$PepDesc)<-c("ID",colnames(DataMatrix)[1]);
    }
  } else {
    Data$MetaData   <- DataMatrix[ (1:(DataDesc$fmt[1]-1)),];
    Data$Matrix     <- DataMatrix[-(1:(DataDesc$fmt[1]-1)),];
    Data$PepDesc    <- rownames(DataMatrix)[-(1:(DataDesc$fmt[1]-1))];
    Data$PepDesc    <- as.matrix(Data$PepDesc, nrow=length(Data$PepDesc))
    rownames(Data$PepDesc)    <- Data$PepDesc;
  }
  NR=nrow(Data$Matrix); DM<- matrix(as.numeric(Data$Matrix), nrow=NR);
  colnames(DM)<- colnames(Data$MetaData); rownames(DM)<- rownames(Data$PepDesc)
  # check if the datamatrix should be extended by metadata rows fN
  if (!is.null(DataDesc$raw$fN)) {
    AddM<-matrix(as.numeric(Data$MetaData[DataDesc$raw$fN,]), ncol=ncol(Data$Matrix));
    DM<-rbind(AddM,DM)    
  }
  Data$Matrix<-DM
  return(Data)
}

###############################################################################
normalizeDataByGroups<- function(Data,DataDesc,normalize){
  # function for group based normalization of data matrix 
  # Author: Vadim Farztdinov
  # Version: 1.0.0, 30-Jan-2021
  
  if (isTRUE(normalize$groupBased)) {
    fq<-table(Data$MetaData[normalize$group,]);
    NG<-length(fq);GrN<-dimnames(fq);
    print("Group based normalisation of samples.  Sample groups:")
    print(fq)
    
    if (all(fq>1)) {
      DOii <- Data;
      for (ii in 1:NG) {
        cond <- Data$MetaData[normalize$group,]==GrN[[1]][ii]
        DOii$MetaData   <- Data$MetaData[,cond];
        DOii$Matrix     <- Data$Matrix[,cond];
        DOii<-normalizeData(DOii,DataDesc,normalize);
        Data$Matrix[,cond] <- DOii$Matrix;
      }
    } else { 
      print("Number of samples  at least in one of the sample groups")
      print("is not sufficient. Normalization skipped")
    }
  } else {
    Data<-normalizeData(Data,DataDesc,normalize);
  }
  return(Data)
}
###############################################################################
normalizeData<- function(Data,DataDesc,normalize){
  # normalize DataObject 
  # Author: Vadim Farztdinov
  # Version: 1.0.0, 09-Sep-2013
  # Version: 1.0.1, 12-Jan-2015. Keep quantile model
  # Version: 1.0.2, 19-Feb-2016. Updated scaling methods, introduced normalization on log2 scale
  # Version: 1.0.3, 13-Oct-2016. Updated scaling methods, introduced normalization on external factor
  # Version: 1.0.4, 07-Nov-2017. Updated scaling methods, introduced "weighted_mean" normalization
  # Version: 1.0.5, 11-Nov-2017. Updated affine and quantile normalization
  # Version: 1.0.6, 13-Nov-2017. Affine normalization set for log scale only
  # Version: 1.1.0, 30-Apr-2019. Added cyclicloess normalization with weights
  # Version: 1.2.0, 29-Jan-2021. Added vsn normalization
  # Version: 1.3.0, 01-Feb-2021. Update of the input and output scale trasformations. Added more diagnostics
  # Version: 1.3.1, 24-Jan-2023. Added quantile scaling (with normalize$qvalue giving quantile). Cyclicloess normaization set for log scale only.
  #
  # function for data normalization 
  # DataDesc$scale         <- "log2"; # scale of data in the input matrix "log2"|"intensity"|"ratio"; default -"log2"
  # normalize[["run"]]     <- TRUE;
  # normalize[["verbose"]] <- FALSE;
  # normalize[["onLogScale"]]<-FALSE=NULL(default)|TRUE;  # normalization is supposed to be carried out on intensity scale
  # normalize[["method"]]  <- "quantile";  # "quantile", "affine", "median", "mean", "trimmed", "sum", "logRMedian", "logRMean",  
  #                                        "scaleon", "scaleByFactor", "weighted_mean", "cyclicloess", "vsn", "q-scale"
  # normalize[["centre"]]    <- "mean";    # "mean"|"median"; "mean" by default scales with SD; "median" by default scales by mad; 
  # normalize[["weightBy"]]<- "Features_Weight";        # factor from PepDesc for "weighted_mean" method
  # normalize[["applyQNorm"]]<-FALSE=NULL(default)|TRUE;
  # normalize[["keepQNorm"]]<-FALSE=NULL(default)|TRUE;
  # normalize[["QNModel"]]<-NULL(default)|column name in Data$PepDesc[,normalize$QNModel]
  # normalize[["target"]]  <- NULL; # used in "affine" as target for center and width of distributions, in scaling methods as target  
  # normalize[["baseln"]] <-NULL| column number to be used as target for scaling
  # normalize[["scaleFactor"]]<- NULL; # Row in MetaData providing a vector Data$MetaData[normalize$scaleFactor,] to be used for scaling 
  # normalize[["plotBefore"]] <-FALSE
  # normalize[["plotAfter"]]  <-TRUE
  
  # require(aroma.light)
  if (is.null(normalize$verbose)) normalize$verbose<-FALSE;
  
  if (!is.double(Data$Matrix)) {
    NR=nrow(Data$Matrix); DM<- matrix(as.numeric(as.matrix(Data$Matrix)), nrow=NR)  
  } else {
    DM <- Data$Matrix
  }
  DM[!is.finite(DM)]<-NA;
  if (normalize$verbose) {
    print("Properties of input data matrix median distribution. Summary: ")
    print(summary(apply(DM,1,median, na.rm=TRUE)))
  }
  
  # Some normalization methods can be applied on both scales, some require specific scale. Fix the scale for them
  if (is.null(normalize$onLogScale)) {normalize[["onLogScale"]]<-FALSE}     # normalization is supposed to be carried out on intensity scale
  if (is.null(normalize$method))   normalize$method <- "median";    # median scaling can be applied by default
  if(normalize$method %in% c("median", "mean", "trimmed", "sum", "vsn", "q-scale") ) {
    normalize$onLogScale<-FALSE}                                            # some of normalization take place on intensity scale
  if(normalize$method %in% c("affine", "cyclicloess") ) {normalize$onLogScale<-TRUE}  # affine normalization takes place on log scale
  
  # Check the scale of input data - by default normalization is intended to be carried out on intensity/ratio scale
  # if ((min(DM, na.rm=TRUE)<0) || (!is.null(DataDesc$scale) && DataDesc$scale=="log2")) {  # used for testing (log)ratio data
  
  InputDataScale="log2";     # by default we assume that the input data are on the log2 scale
  if (!is.null(DataDesc$scale)) {
    if (DataDesc$scale %in% c("intensity", "ratio")) InputDataScale="intensity";     
  } else if (is.null(DataDesc$scale)) {
    cat("Warning! Input data scale for normalization was not specified (DataDesc$scale is null). Trying to define it automatically\n", file = "", sep = "")
    if (median(DM, na.rm=TRUE)>32) {
      InputDataScale="intensity"; 
    } else if (min(DM, na.rm=TRUE)<0) {
      InputDataScale="log2";     
    } else if (min(DM, na.rm=TRUE)>=0 && max(DM, na.rm = TRUE)<=1) {
      InputDataScale="ratio";
    } 
    cat("Input data scale was automatically defined as ", InputDataScale, "\n", file = "", sep = "")
    cat("Please specify it in the DataDesc$scale if it was defined incorrectly\n", file = "", sep = "")
  }

  backTrans2Int=FALSE; backTrans2Log=FALSE;
  if (isTRUE(normalize$onLogScale) && InputDataScale %in% c("intensity", "ratio")) {
    DM <- log2(DM); backTrans2Int=TRUE;  
    cat("Data internally transformed to log2 scale\n", file = "", sep = "")
    if (normalize$verbose) {
      print("Properties of log2 transformed input data matrix median distribution. Summary: ")
      print(summary(apply(DM,1,median, na.rm=TRUE)))
    }
  } else  if (!isTRUE(normalize$onLogScale) && InputDataScale == "log2") {
    DM <- 2^DM; backTrans2Log=TRUE;  
    cat("Data internally power transformed to intensity scale\n", file = "", sep = "")
    if (normalize$verbose) {
      print("Properties of power transformed input data matrix median distribution. Summary: ")
      print(summary(apply(DM,1,median, na.rm=TRUE)))
    }
  }

  # Apply normalization
  if (normalize$method == "quantile") {
    require(aroma.light)
    if ((is.null(normalize[["applyQNorm"]])) || (!normalize[["applyQNorm"]])) {
      DM <- normalizeQuantileRank(DM); # normalizeQuantileRank
      QNModel<- averageQuantile(as.list(as.data.frame(DM)))
      if (!is.null(normalize$keepQNorm) && (normalize$keepQNorm) && (!is.null(Data$PepDesc))) {
        Data$PepDesc<-  cbind(Data$PepDesc,QNModel); 
        colnames(Data$PepDesc)[ncol(Data$PepDesc)]<-"QNModel";
      }
      if (normalize$verbose) {
        print("Properties of data matrix quantile distribution. Summary: ")
        print(summary(QNModel))
        if(!isTRUE(normalize$onLogScale)) {
          print("Properties of log2 tranformed quantile distribution. Summary: ")
          print(summary(log2(QNModel)))
          mainText<-sprintf("log2(QNModel) mean = %g and sd = %g", mean(log2(QNModel)),sd(log2(QNModel))) 
          print(mainText)
        }
      }
    } else if (isTRUE(normalize[["applyQNorm"]])) {
      target<- sort(Data$PepDesc[,normalize$QNModel], decreasing = FALSE, na.last = NA)
      DM <- normalizeQuantileRank(DM, xTarget=target)
    }
    
  } else if (normalize$method == "cyclicloess") {
    # Normalize the columns of a matrix, cyclicly applying loess normalization to normalize each pair of columns to each other.
    require(limma)
    if (!is.null(normalize$weightBy)) {Weight<-as.numeric(Data$PepDesc[,normalize$weightBy])} else {Weight<-NULL}
    if (!is.null(normalize$option)) {Option <- normalize$option} else {Option<-"fast"}
    if (!is.null(normalize$iter)) {iter <- normalize$iter} else {iter<-3}
    DM <- normalizeCyclicLoess(DM, weights = Weight, span=0.7, iterations = iter, method = Option)
    # if (normalize$verbose) {
    #   print("Properties of data matrix cyclicloess reference distribution. Summary: ")
    #   print(summary(apply(DM,1,mean, na.rm=TRUE)))
    # }
    
  } else if (normalize$method == "vsn") {
    # Apply variance stabilizing normalization (vsn) to data objects. The input DM should contain raw intensities. 
    if (min(DM, na.rm=TRUE)<0) { DM <- 2^DM } # normalization is carried out on intensity/ratio scale
    require(limma)
    DM <- normalizeVSN(DM);  
    # switch off back transformations. We assume that the input/output data are on log scale
    backTrans2Log=FALSE; backTrans2Int=FALSE;
    if (InputDataScale %in% c("intensity","ratio")) backTrans2Int <- TRUE  # do transform data back to input intensity scale 
    # if (normalize$verbose) {
    #   print("Properties of data matrix vsn normalized mean distribution. Summary: ")
    #   print(summary(apply(DM,1,mean, na.rm=TRUE)))
    # }
  } else {
    require(aroma.light)
    
    if (normalize$method == "affine") {
      # robust centering and scaling of columns of a numeric matrix, worth to apply on log scale only
      # DM<- scale(DM, center =TRUE, scale = TRUE)
      centre <- ifelse(is.null(normalize$centre), "mean", normalize$centre)
      
      WD <- vector("list"); WD$standardize <- vector("list");
      WD$standardize$run   <-TRUE                # TRUE|FALSE
      WD$standardize$method<-centre;            # "mean"|"median"; "mean" by default scales with SD; "median" by default scales by mad; 
      WD$standardize$center<-TRUE;                # logical value 
      WD$standardize$scale <-TRUE;                # logical value 
      DM<- t(standardizeData(t(DM),DataDesc,WD))
      
      if (normalize$verbose) {
        print("Properties of data matrix sample distributions. Centres: "); print(attr(DM,"scaled:center"))
        print("Scale: "); print(attr(DM,"scaled:scale") )
      }
      mainText<-sprintf("Mean(Centre) = %g and mean(Scale) = %g", mean(attr(DM,"scaled:center")),mean(attr(DM,"scaled:scale"))) 
      print(mainText)
      
      if (!is.null(normalize$target)) {
        if (is.numeric(normalize$target)) {
          if (length(normalize$target)!=2) {
            normalize$target = c(0,1)
            cat("Warning!!! (pre)normalize$target is introduced incorrectly. Redefined to default =c(0,1) \n", file = "", sep = "")
          }
          DM <- DM*normalize$target[2]+normalize$target[1];  
          cat("Data matrix affine transformed to matrix with samples/column mean = ",normalize$target[1], " and sd = ", normalize$target[2], "\n", file = "", sep = "")
        } else if (normalize$target=="mean") {
          DM <- DM*mean(attr(DM,"scaled:scale")) + mean(attr(DM,"scaled:center"))
          cat("Data matrix affine transformed to matrix with samples/column mean = ",mean(attr(DM,"scaled:center")), " and sd = ", mean(attr(DM,"scaled:scale")), "\n", file = "", sep = "")
        } else {
          cat("Warning!!! (pre)normalize$target is introduced incorrectly. Redefined to default = c(0,1)\n", file = "", sep = "")
        }
      } else {
        print("Data matrix affine transformed to matrix with samples/column mean = 0 and sd = 1")
      }
      ############################  
    } else if (normalize$method %in% c("median", "mean", "trimmed", "sum", "weighted_mean", "q-scale")) {      
      # require(robustbase)
      require(matrixStats)
      if (normalize$method=="median") {
        meds<-colMedians(DM, na.rm=TRUE);
      } else if (normalize$method=="mean") {
        meds<-colMeans(DM, na.rm=TRUE);
      } else if (normalize$method=="weighted_mean") {
        # require(stats); 
        weight<-as.numeric(Data$PepDesc[,normalize$weightBy]);   # VF added as.numeric 12May2019
        meds<-apply(DM, 2, weighted.mean, w=weight, na.rm=TRUE);
      } else if (normalize$method=="trimmed") {
        meds<-apply(DM, 2, mean, trim=0.2, na.rm=TRUE);
      } else if (normalize$method=="sum") {
        meds<-colSums(DM, na.rm=TRUE);
      } else if (normalize$method=="q-scale") {
        if (is.null(normalize$qvalue)) normalize$qvalue <- 0.75
        meds<-apply(DM, 2, quantile, probs=normalize$qvalue, na.rm=TRUE);
      };     
      if (!is.null(normalize$target) && is.numeric(normalize$target)) {
        denom<-normalize$target;
      } else if (!is.null(normalize$baseln) && is.numeric(normalize$baseln)){
        denom<-meds[normalize$baseln]
      } else {
        denom<-median(meds);
      }
      scaleV<-meds/denom;
      DM<-scale(DM, center=FALSE, scale=scaleV)
      if (normalize$verbose) {
        if (!is.null(normalize$onLogScale) && normalize$onLogScale) {addText="Data scaled on log2 scale by vector\n"
        } else {addText="Data scaled by vector\n"}
        cat(addText, paste(scaleV, collapse=", "), "\nto common", normalize$method, "value",denom, "\n" , file = "", sep = " ")
        Data$MetaData<-rbind(Data$MetaData,scaleV)
      }      
      #     ############################  
      #     } else if (normalize$method %in% c("median", "mean", "sum")) {
      #       if (!is.null(normalize$target) && is.numeric(normalize$target)) {
      #         DM <- normalizeAverage(DM, avg=normalize$method, targetAvg=normalize$target, na.rm=TRUE)
      #       } else {
      #         DM <- normalizeDifferencesToAverage(as.list(as.data.frame(DM)),baseline=normalize$baseln, avg=normalize$method, na.rm=TRUE)  #  
      #         DM <- as.matrix(as.data.frame(DM))
      #       }
      #     ############################  
      ############################  
    } else if (normalize$method == "logRMedian") {
      DM <- median.scale.to.zero(DM);
    } else if (normalize$method == "logRMean") {
      DM <- mean.scale.to.zero(DM);
    } else if (normalize$method == "scaleByFactor") {
      scaleV<-as.numeric(Data$MetaData[normalize$scaleFactor,]);
      if (length(scaleV)!=ncol(DM)) {
        stop("number of factors does not coincide with number of columns in datamatrix, please adjust\n")
      }
      DM<-scale(DM, center=FALSE, scale=scaleV)
    } else if (normalize$method == "scaleon") {
      rownames(DM)<- rownames(Data$PepDesc)
      DM <- scaleon(DM,normalize);
    }
    
    detach("package:aroma.light")
  }
  
  # Normalization finished, get column and row names
  colnames(DM)<- colnames(Data$MetaData)
  rownames(DM)<- rownames(Data$PepDesc)
  # now apply back transformations
  if (backTrans2Int && backTrans2Log) {
    cat("Warning!!! Something went wrong during the normalization - cannot apply tranformation to log2 and intensity scales simultaneouosly\n", file = "", sep = "")
  } else if (backTrans2Log && !backTrans2Int) {   # Were input data on log scale? Then back to log scale
    DM <- log2(DM)
    if (normalize$verbose) {
      cat("Properties of back log2 transformed ",normalize$method," normalized data matrix median distribution. Summary: \n", file = "", sep = "")
      print(summary(apply(DM,1,median, na.rm=TRUE)))
    }
  } else if (backTrans2Int && !backTrans2Log) { 
    DM <- 2^(DM)
    if (normalize$verbose) {
      cat("Properties of back power transformed ",normalize$method," normalized data matrix median distribution. Summary: \n", file = "", sep = "")
      print(summary(apply(DM,1,median, na.rm=TRUE)))
    } 
  } else {
    cat("Properties of ",normalize$method," normalized data matrix median distribution. Summary: \n", file = "", sep = "")
    print(summary(apply(DM,1,median, na.rm=TRUE)))
  }
  
  
  Data$Matrix<-DM
  return(Data)
}

###############################################################################
normalizeFastLoess <- function(x, xTarget = NULL, weights = NULL, span=0.7, iterations = 3)
{
  #	Cyclic loess normalization of columns of matrix x to its rowMeans or to xTarget
  #	incorporating probes weights.
  #	Yunshun (Andy) Chen and Gordon Smyth
  #	14 April 2010.  Last modified 24 Feb 2012.
  #	05 July  2023.  Vadim Farztdinov added xTarget for normalisation to target distribution
  x <- as.matrix(x)
  n <- ncol(x)
  if (!is.null(xTarget)) {iterations  <- 1; a <- xTarget}
  
  for (k in 1:iterations) {
    if (is.null(xTarget)) {
      a <- rowMeans(x,na.rm=TRUE)
    }
    for (i in 1:n){
      m <- x[,i] - a
      f <- limma::loessFit(m, a, weights=weights, span=span)$fitted
      x[,i] <- x[,i] - f
    }
  }
  return(x)
}

###############################################################################
standardizeVector<- function(Vvalues,DataDesc=NULL,WorkDesc){
  # To standardize a vector 
  # Author: Vadim Farztdinov
  # Version: 1.0.0, 19-Sep-2020
  #
  # WorkDesc$standardize$method="mean"|"median"|anything else; "mean" by default
  #    scales with SD; "median" by default scales by mad; specify your own values
  #    for centre and scale otherwise 
  # WorkDesc$standardize$center=either a logical value or a numerical value
  # WorkDesc$standardize$scale =either a logical value or a numerical value
  
  # check definitions
  if (is.null(WorkDesc$standardize$method)) WorkDesc$standardize$method<-"mean"
  if (is.null(WorkDesc$standardize$center)) WorkDesc$standardize$center<-TRUE
  if (is.null(WorkDesc$standardize$scale))  WorkDesc$standardize$scale <-TRUE
  
  # standardize
  if (WorkDesc$standardize$method=="mean") {
    Vvalues<- scale(Vvalues, center = WorkDesc$standardize$center, scale = WorkDesc$standardize$scale)
  } else if (WorkDesc$standardize$method=="median") {
    if (WorkDesc$standardize$center==TRUE) {Vmedians=median(Vvalues,na.rm = TRUE)
    } else {
      Vmedians= WorkDesc$standardize$center
    }
    if (WorkDesc$standardize$scale==TRUE) {
      Vmads=mad(Vvalues,na.rm = TRUE); 
      Vvalues<- scale(Vvalues, center = Vmedians, scale = Vmads)
    } else 
      Vvalues<- scale(Vvalues, center = Vmedians, scale = WorkDesc$standardize$scale)
  }
  return(Vvalues)
}

###############################################################################
standardizeData<- function(Data,DataDesc=NULL,WorkDesc){
  # called function scale centres and/or scales the columns of a numeric matrix
  # To standardize over the samples (along rows) data matrix is transposed
  # Author: Vadim Farztdinov
  # Version: 1.0.0, 09-Sep-2013
  # Version: 1.0.1, 11-Jan-2015. Added return of centre and scale and saving them in Annotation slot
  # Version: 1.0.2, 19-Aug-2020. Updated for more control over parameters
  #
  # WorkDesc$standardize$method="mean"|"median"|anything else; "mean" by default
  #    scales with SD; "median" by default scales by mad; specify your own vectors
  #    for centre and scale otherwise 
  # WorkDesc$standardize$center=either a logical value or a numeric vector of
  #    length equal to the number of rows
  # WorkDesc$standardize$scale=either a logical value or a numeric vector of
  #    length equal to the number of rows
  
  # check definitions
  if (is.null(WorkDesc$standardize$method)) WorkDesc$standardize$method<-"mean"
  if (is.null(WorkDesc$standardize$center)) WorkDesc$standardize$center<-TRUE
  if (is.null(WorkDesc$standardize$scale))  WorkDesc$standardize$scale <-TRUE
  
  # standardize  
  if (is.list(Data)) DM<-t(Data$Matrix) else DM<-t(Data);
  if (WorkDesc$standardize$method=="mean") {
    DM<- scale(DM, center = WorkDesc$standardize$center, scale = WorkDesc$standardize$scale)
  } else if (WorkDesc$standardize$method=="median") {
    if (WorkDesc$standardize$center==TRUE) {DMmedians=apply(DM,2,median,na.rm = TRUE)
    } else {
      DMmedians= WorkDesc$standardize$center
    }
    if (WorkDesc$standardize$scale==TRUE) {
      DMmads=apply(DM,2,mad,na.rm = TRUE); 
      DM<- scale(DM, center = DMmedians, scale = DMmads)
    } else 
      DM<- scale(DM, center = DMmedians, scale = WorkDesc$standardize$scale)
  } else {
    DM<- scale(DM, center = WorkDesc$standardize$center, scale = WorkDesc$standardize$scale)
  }
  DMCenter<-attr(DM, "scaled:center")
  DMScale<-attr(DM, "scaled:scale")
  
  if (is.list(Data)) {
    Data$Matrix <-t(DM); 
    SA<-cbind(DMCenter, DMScale); colnames(SA)<-c("DMCenter","DMScale")
    Data$PepDesc<- cbind(Data$PepDesc,SA) 
  } else {Data<-t(DM)};
  return(Data)
}

###############################################################################
transformData<- function(Data,DataDesc,WorkDesc){
  # transform values of DataObject 
  # Author: Vadim Farztdinov
  # Version: 1.0.0, 09-Sep-2013
  # Version: 1.0.1, 28-Aug-2015. Added tanh and atanh
  # Version: 1.0.2, 18-Oct-2017. Added shifting by 1 before transformation and replacement infinite values by NAs
  
  if (!is.null(WorkDesc$transform$withShift) && WorkDesc$transform$withShift) {
    Data$Matrix <- Data$Matrix + 1
  }
  if (WorkDesc$transform$method == "log2") {
    Data$Matrix     <- log2(Data$Matrix);
  } else if (WorkDesc$transform$method == "log10") {
    Data$Matrix     <- log10(Data$Matrix);
  } else if (WorkDesc$transform$method == "log") {
    Data$Matrix     <- log(Data$Matrix);
  } else if (WorkDesc$transform$method == "exp") {
    Data$Matrix     <- exp(Data$Matrix);
  } else if (WorkDesc$transform$method == "powerOf2") {
    Data$Matrix     <- 2^(Data$Matrix);
  } else if (WorkDesc$transform$method == "powerOf10") {
    Data$Matrix     <- 10^(Data$Matrix);
  } else if (WorkDesc$transform$method == "tanh") {
    Data$Matrix     <- tanh(Data$Matrix);
  } else if (WorkDesc$transform$method == "atanh") {
    Data$Matrix     <- atanh(Data$Matrix);
  }
  Data$Matrix[!is.finite(Data$Matrix)]<-NA;
  return(Data)
}

###############################################################################