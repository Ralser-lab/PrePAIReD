###############################################################################
# Library of functions used in Preprocessing Workflows 
# Author: Vadim Farztdinov
# Tested on R version 3.2.3 (2015-12-10), nickname Wooden Christmas-Tree
# ...
# Tested on R version 4.3.3 (2024-02-29), nickname Angel Food Cake  
###############################################################################
# Version: 1.0.0, 09-Sep-2013. Vadim Farztdinov
# Version: 4.4.6, 03-Jan-2025. Updated reduceData: introduced option for sample names from Metadata
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
makeMeVMatr<- function(Data){
  # create DataMatrix from DataObject
  # Author: Vadim Farztdinov
  # Version: 1.0.0, 02-Oct-2013
  # Version: 1.0.1, 22-Jan-2014. Added as.matrix(data$PepDesc)
  # Version: 1.0.2, 02-Aug-2016, modified by Vadim Farztdinov, added row names for AddMatr
  nrowAdd<-nrow(Data$MetaData); ncolAdd<-ncol(Data$PepDesc); 
  AddMatr<-matrix(data=rownames(Data$MetaData),nrow=nrowAdd,ncol=ncolAdd,byrow=FALSE);
  colnames(AddMatr)<-colnames(Data$PepDesc); rownames(AddMatr)<-rownames(Data$MetaData);
  DM<-cbind( rbind(AddMatr,as.matrix(Data$PepDesc)), rbind(Data$MetaData,Data$Matrix) );
  return(DM)
}

###############################################################################
reduceData<- function(Data,DataDesc,WorkDesc){
  # reduce DataObject 
  # Author: Vadim Farztdinov
  # Version: 1.0.2, 09-Jan-2014
  # find subset using group information in row named as WorkDesc$subset
  # or remove samples specified in the list
  # Version: 1.0.3, 22-Jul-2015. Added control for WorkDesc$removeSamples$run being logical
  # Version: 1.1.0, 03-Jan-2025. Added reduction of samples by names in metadata
  cond<-TRUE;
  if (is.null(WorkDesc$removeSamples$byNames)) {sampleNames <- colnames(Data$Matrix)
  } else {sampleNames <- (Data$MetaData[WorkDesc$removeSamples$byNames,])}
  
  if (is.logical(WorkDesc$removeSamples$run) && (WorkDesc$removeSamples$run)) {
    cond<- cond & !(sampleNames %in% WorkDesc$removeSamples[["list"]])
  }
  if (!is.null(WorkDesc$subset$run)) {
    if (WorkDesc$subset$run) {
      icon<-which(names(WorkDesc$subset) %in% names(DataDesc$row))
      if (length(icon)>=1){
        for (ii in icon) {
          cond<-cond & (Data$MetaData[DataDesc$row[[ names(WorkDesc$subset)[ii] ]],] %in% WorkDesc$subset[[ii]])
        }
      }
    } 
  }
  Data$MetaData   <- Data$MetaData[,cond];
  Data$Matrix     <- Data$Matrix[,cond];
  return(Data)
}

###############################################################################
removeNAs<- function(Data,DataDesc,WorkDesc){
  # remove NAs from DataObject 
  # Author: Vadim Farztdinov
  # Version: 1.0.0, 09-Sep-2013
  # Version: 1.1.0, 09-Nov-2013. Group based removal of NAs. Works only with %AllowedNA
  # Version: 1.1.1, 18-Apr-2014. Work only with finite values, all others are replaced by NA
  # Version: 1.1.2, 25-Jun-2015. Introduced verbose option
  # Version: 1.1.3, 28-Jan-2021. Added possibility to clean groups with one sample
  # Version: 2.0.0, 29-Mar-2024. Cleaning by presence in number of groups for groupBased NA removal
  # verbose<-TRUE; if (!is.null(WorkDesc[["verbose"]]) &&  !WorkDesc[["verbose"]]) verbose<- FALSE;
  if (is.null(WorkDesc$removeNA$PresenceInNGroups)) {
    PresenceInNGroups=1
  } else {
    PresenceInNGroups <- WorkDesc$removeNA$PresenceInNGroups
  }
  DM<-Data$Matrix;
  DM[!is.finite(DM)]<-NA;
  
  if (!WorkDesc$removeNA$groupBased)  {
    # remove features by absense in total set
    if (!is.null(WorkDesc$removeNAs$allowNA)) {
      # Obsolete option - remove by number of NAs in feature
      NAllowed<-WorkDesc$removeNAs$allowNA
    } else {
      # default option - remove by % of NAs in feature
      NAllowed <- floor(WorkDesc$removeNAs$allowPrcNA*ncol(DM)/100);
    }
    if (is.data.frame(DM)) {
      cond<- (rowSums(is.na.data.frame(DM)) <= NAllowed);  
    } else {
      cond<- (rowSums(is.na(DM)) <= NAllowed); 
    }
    # if (verbose) {cat(sprintf("Number of NA allowed per feature = %i\n", NAllowed) , file = "", sep = " ")}
    # end of remove features by absense in total set 
  } else {
    # remove features by absense in a group - group based removal
    # if (verbose) {cat(sprintf("Group = %s based NA removal. Allowed NAs per feature = %3.1f percent\n", WorkDesc$removeNAs$group, WorkDesc$removeNAs$allowPrcNA) , file = "", sep = " ")}
    # get the number of samples in each group
    fq<-table(Data$MetaData[WorkDesc$removeNAs$group,]);
    NG<-length(fq);GrN<-dimnames(fq);
    presInGr <- rep_len(0, nrow(DM));
    # find features that satisfy the requirement in every group 
    for (ii in 1:NG) {
      DMii<-DM[,Data$MetaData[WorkDesc$removeNAs$group,]==GrN[[1]][ii] ];
      if (fq[ii]>1) {
        if (is.data.frame(DM)) {
          presInGr<- presInGr + as.numeric(rowSums(is.na.data.frame(DMii)) <= floor(WorkDesc$removeNAs$allowPrcNA*fq[ii]/100) );  
        } else {
          presInGr<- presInGr + as.numeric(rowSums(is.na(DMii)) <= floor(WorkDesc$removeNAs$allowPrcNA*fq[ii]/100) ); 
        }
      } else {  # if only one sample in a group
        presInGr<- (presInGr & as.numeric(!is.na(DMii))) 
      }
    }
    cond <- presInGr>=PresenceInNGroups
  } # end of remove features by absense in a group
  Data$PepDesc   <- Data$PepDesc[cond,];
  Data$Matrix     <- Data$Matrix[cond,];
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
summarizeData<- function(Data,DataDesc,WorkDesc){
  # summarize DataObject from peptide level to protein level
  # Author: Vadim Farztdinov
  # Version: 1.0.0, 09-Jan-2014
  # Version: 1.1.0, 13-Sep-2017. Added annotation columns to keep after summarization
  # Version: 1.1.1, 08-May-2018. Added replacement of generated NaNs by NAs
  # Version: 1.1.2, 20-Aug-2018. Added control on the number of peptide names (<=300) written into the column "PeptideNames"
  # Version: 1.2.0, 22-Feb-2019. Added flexibility for N_PSMs column name, added flexibility for Anno2keep; redesigned PA matrix
  # Version: 1.2.1, 23-Mar-2019. Added control on the number of item names (<=50) written into char column
  # Version: 1.2.2, 30-Apr-2019. Updated Anno2sum - performed before Anno2keep
  # Version: 1.2.3, 15-Aug-2019. Updated sumAnno
  # Version: 1.3.0, 23-Oct-2019. Introduced new methods medianPolish, PLM and PLMr
  # Version: 1.3.1, 24-Oct-2019. Added WPLM method
  # Version: 1.3.2, 28-Oct-2019. Updated WPLM method, introduced sample weights
  # Version: 1.3.3, 17-Aug-2020. Added summarization by Max
  # Version: 1.3.4, 21-Aug-2020. bebugged summarization by Max
  # Version: 1.3.5, 23-Sep-2020. Added maxLFQ summarization to function summarizeData
  # Version: 1.3.6, 27-Jan-2021. Updated maxLFQ summarization to allow summarization of only one numerical column in AM
  # Version: 1.3.7, 07-Apr-2021. Added as.matrix control of Data$PepDesc 
  # Version: 1.3.8, 07-Feb-2022. Added reduction of Data$PepDesc to required columns only, added gc at the protein cycle end
  # Version: 1.3.9, 13-Nov-2023. Added power for weight
  
  addAnno<-FALSE; addSumOfPSMs<-FALSE; # VF 13Sep2017
  sumAnno<-FALSE;
  Data$PepDesc <- as.matrix(Data$PepDesc)   # VF 07Apt2021
  
  if (!is.double(Data$Matrix)) {
    NR=nrow(Data$Matrix); DM<- matrix(as.numeric(as.matrix(Data$Matrix)), nrow=NR)  
  } else {
    DM <- Data$Matrix
  }
  # define default columns to keep in summarized annotation matrix
  PAcolnames<-c(WorkDesc$summarize$Proteins,DataDesc$col$ProteinDesc, "PeptideNN","PeptideNames","SumOfUniqueProteins"); 
  ncolPA_basic<-length(PAcolnames);
  
  # Take over the number of PSMs if provided
  N_PSMs_names<- c("N_PSMs", "NN_PSMs", "NN_PSM")
  N_PSMs_ind<-which(N_PSMs_names %in% colnames(Data$PepDesc))
  if (length(N_PSMs_ind)>0) { # VF 21Feb2019
    ncolPA<-ncolPA_basic+1; addSumOfPSMs<-TRUE; PSM_Col<-N_PSMs_names[N_PSMs_ind];
    PAcolnames<-c(PAcolnames,PSM_Col);
    } else {ncolPA<-ncolPA_basic}   
  
  # Check the additional annotation columns to sum and/or keep
  if (!is.null(DataDesc$col$Anno2sum)) { # VF 21Feb2019
    Anno_ind<-which(DataDesc$col$Anno2sum %in% colnames(Data$PepDesc))
    if (length(Anno_ind)==0) stop("Anno2sum slot is empty. Please specify or switch off")
    Anno2sum<-DataDesc$col$Anno2sum[Anno_ind]  # use only present in PepDesc columns
    Anno_i2rem<-which(Anno2sum %in% c(PAcolnames,N_PSMs_names)) 
    if (length(Anno_i2rem)>0) Anno2sum<-Anno2sum[-Anno_i2rem]; # remove overlaps if they exist
    AM<- matrix(as.numeric(as.matrix(Data$PepDesc[,Anno2sum])), ncol=length(Anno2sum))  
    ncolPA<-ncolPA+length(Anno2sum); sumAnno<-TRUE
    PAcolnames<-c(PAcolnames, Anno2sum);
  } 
  if (!is.null(DataDesc$col$Anno2keep)) { # updated VF 21Feb2019
    Anno_ind<-which(DataDesc$col$Anno2keep %in% colnames(Data$PepDesc))
    Anno2keep<-DataDesc$col$Anno2keep[Anno_ind]  # use only present in PepDesc columns
    Anno_i2rem<-which(Anno2keep %in% c(PAcolnames,N_PSMs_names)) 
    if (length(Anno_i2rem)>0) Anno2keep<-Anno2keep[-Anno_i2rem]; # remove overlaps if they exist
    ncolPA<-ncolPA+length(Anno2keep); addAnno<-TRUE;
    PAcolnames<-c(PAcolnames, Anno2keep);
    } 
  
  UniqueProt<-unique(Data$PepDesc[,WorkDesc$summarize$Proteins]);
  PM <- matrix(data = NA, nrow = length(UniqueProt), ncol = ncol(DM));
  colnames(PM)<- colnames(Data$MetaData); rownames(PM)<- UniqueProt;
  
  PA <- matrix(data = NA, nrow = length(UniqueProt), ncol = ncolPA);
  colnames(PA)<- PAcolnames; rownames(PA)<- UniqueProt;
  
  # reduce size of PepDesc to speed up for big data
  keepInd <- which(colnames(Data$PepDesc) %in% c(PAcolnames, WorkDesc$summarize$weightBy))
  Data$PepDesc <- Data$PepDesc[, unique(colnames(Data$PepDesc)[keepInd])]
  gc()
  
  # internal function
  # updated by Vadim Farztdinov 23Mar2019   # collapseUnique<- function(x) paste(unique(x), collapse = "+")
  collapseUnique<- function(x) ifelse(length(unique(x))<=50, paste(unique(x), collapse = "+"), paste(c(unique(x)[1:50], "and_other"), collapse = "+")) 

  if (WorkDesc$summarize$method %in% c("max","median","mean","trimmedmean","weighted_mean","maxLFQ","medianPolish", "PLM", "PLMr", "WPLM")) {
    for (ii in 1: length(UniqueProt)) {
      cond<-which(Data$PepDesc[,WorkDesc$summarize$Proteins] %in% UniqueProt[ii]);
      if (length(cond)<=300) {
        PA[ii,1:ncolPA_basic]<- c(Data$PepDesc[cond[1],c(WorkDesc$summarize$Proteins,DataDesc$col$ProteinDesc)], 
                       length(cond),paste(rownames(Data$PepDesc)[cond],collapse="+"),
                       collapseUnique(Data$PepDesc[cond,DataDesc$col$Protein]) );
      } else {
        PA[ii,1:ncolPA_basic]<- c(Data$PepDesc[cond[1],c(WorkDesc$summarize$Proteins,DataDesc$col$ProteinDesc)], 
                       length(cond),
                       paste(c(rownames(Data$PepDesc)[cond[1:300]],"other_Peptides"), collapse="+"),
                       collapseUnique(Data$PepDesc[cond,DataDesc$col$Protein]) );
      }
      if (addSumOfPSMs) { PA[ii,PSM_Col]<- sum(as.numeric(Data$PepDesc[cond, PSM_Col])) }    # VF 13Sep2017 & 22Feb2019
      if (length(cond)>1) {
        
        if (addAnno) { PA[ii,Anno2keep]<- apply(Data$PepDesc[cond, Anno2keep],2,collapseUnique) } # VF 13Sep2017 & 22Feb2019
        if (WorkDesc$summarize$method == "median") {
          PM[ii,]<- apply(DM[cond,],2,median, na.rm = TRUE); # , na.rm = TRUE
          if (sumAnno) { PA[ii,Anno2sum]<- apply(AM[cond, ],2, median, na.rm = TRUE) } # VF 22Feb2019
        } else if (WorkDesc$summarize$method == "max") {
          PM[ii,]<- apply(DM[cond,],2,max, na.rm = TRUE); 
          if (sumAnno) { PA[ii,Anno2sum]<- apply(AM[cond, ],2, max, na.rm = TRUE) } # VF 19Aug2020
        } else if (WorkDesc$summarize$method == "mean") {
          PM[ii,]<- apply(DM[cond,],2,mean, na.rm = TRUE); 
          if (sumAnno) { PA[ii,Anno2sum]<- apply(AM[cond, ],2, mean, na.rm = TRUE) } # VF 22Feb2019
        } else if (WorkDesc$summarize$method == "trimmedmean") {
          PM[ii,]<- apply(DM[cond,],2,mean, trim = WorkDesc$summarize$trim, na.rm = TRUE); 
          if (sumAnno) { PA[ii,Anno2sum]<- apply(AM[cond, ],2,mean, trim = WorkDesc$summarize$trim, na.rm = TRUE) } # VF 22Feb2019
        } else if (WorkDesc$summarize$method=="weighted_mean") {
          # require(stats); 
          NP <- ifelse(!is.null(WorkDesc$summarize$wgtPower), WorkDesc$summarize$wgtPower,1)
          FWeight<-(as.numeric(Data$PepDesc[cond,WorkDesc$summarize$weightBy]))^NP;
          FWeight[is.na(FWeight)] <- 0;
          PM[ii,]<- apply(DM[cond,],2, weighted.mean, w=FWeight, na.rm=TRUE);
          if (sumAnno) { PA[ii,Anno2sum]<- apply(AM[cond, ],2,weighted.mean,w=FWeight,na.rm=TRUE) } # VF 22Feb2019
        } else if (WorkDesc$summarize$method == "maxLFQ") {  # VF 23Sep2020
          PM[ii,]<- iq::maxLFQ(DM[cond,])$estimate
          # maxLFQ requires ion intensities in log2 space and cannot be applied PepDesc data with arbitrary scale. 
          # Median summarization is applied instead. VF 10.042021  
          if (sumAnno) { PA[ii,Anno2sum]<- apply(AM[cond, ],2, median, na.rm = TRUE) } 
        } else if (WorkDesc$summarize$method == "medianPolish") {  # VF 24Oct2019 
          require(preprocessCore)
          rcm <- rcModelMedianPolish(DM[cond,])
          PM[ii,]<- rcm$Estimates[1:ncol(DM)];  # ErrM[ii,]<- rcm$StdErrors[1:ncol(DM)]
          if (sumAnno) { PA[ii,Anno2sum]<- apply(AM[cond, ],2, median, na.rm = TRUE) } 
        } else if (WorkDesc$summarize$method == "PLM") {  # VF 24Oct2019 
          require(preprocessCore)
          rcm <- rcModelPLM(DM[cond,])
          PM[ii,]<- rcm$Estimates[1:ncol(DM)];  # ErrM[ii,]<- rcm$StdErrors[1:ncol(DM)]
          if (sumAnno) { PA[ii,Anno2sum]<- apply(AM[cond, ],2,mean, trim = 0.2, na.rm = TRUE) } 
        } else if (WorkDesc$summarize$method == "PLMr") {  # VF 24Oct2019 
          require(preprocessCore)
          rcm <- rcModelPLMr(DM[cond,])
          PM[ii,]<- rcm$Estimates[1:ncol(DM)];  # ErrM[ii,]<- rcm$StdErrors[1:ncol(DM)]
          if (sumAnno) { PA[ii,Anno2sum]<- apply(AM[cond, ],2,mean, trim = 0.2, na.rm = TRUE) } 
        } else if (WorkDesc$summarize$method == "WPLM") {  # VF 24Oct2019 
          require(preprocessCore)
          NP <- ifelse(!is.null(WorkDesc$summarize$wgtPower), WorkDesc$summarize$wgtPower,1)
          FWeight<-(as.numeric(Data$PepDesc[cond,WorkDesc$summarize$weightBy]))^NP;
          FWeight[is.na(FWeight)] <- 0;
          if (!is.null(WorkDesc$summarize$SamWeight)) {
            SWeight<-as.numeric(Data$MetaData[WorkDesc$summarize$SamWeight,])
          } else {SWeight <- rep(1,times=ncol(DM))}
          Weights<-outer(FWeight, SWeight,"*")
          rcm <- rcModelWPLM(DM[cond,], w=Weights)
          PM[ii,]<- rcm$Estimates[1:ncol(DM)];  # ErrM[ii,]<- rcm$StdErrors[1:ncol(DM)]
          if (sumAnno) { PA[ii,Anno2sum]<- apply(AM[cond, ],2,weighted.mean,w=FWeight,na.rm=TRUE) } 
        }
        
      } else {
        PM[ii,]<- DM[cond,]; 
        PA[ii,Anno2keep]<-Data$PepDesc[cond, Anno2keep]
        if (sumAnno) PA[ii,Anno2sum]<- AM[cond, ]    # added by Vadim Farztdinov 15-Aug-2019
      }   # end if (length(cond)>1)
      rm(cond); 
      gc()
    }   # end for
  }  # end if (WorkDesc$summarize$method...)
  
  # replace generated NaNs by NAs
  PM[is.nan(PM)] <- NA;
  PM[!is.finite(PM)] <- NA;
  
  Data$Matrix<-PM
  Data$PepDesc<-PA
  return(Data)
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
reduceFeatures<- function(Data,Selection){
  # reduce features in DataObj
  # Author: Vadim Farztdinov
  # Version: 1.0.0, 12-Sep-2013
  ###                   
  DOS<-Data;
  DOS$Matrix <- Data$Matrix[Selection,]
  DOS$PepDesc<- Data$PepDesc[Selection,]
  return(DOS)
}

###############################################################################
as.numeric.matrix <- function(Matrix){
  # function convert a character matrix to a numerical matrix
  # param Matrix is a matrix with character like "1"
  # Author: from internet?
  # Version: 1.0.0, 12-Sep-2013

  if(ncol(Matrix)<=nrow(Matrix)){
    nMatrix <- apply(Matrix, 2, as.numeric)
  }else{
    nMatrix <- t(apply(Matrix, 1, as.numeric))
  }
  if((dim(nMatrix)[1] != dim(Matrix)[1]) && (dim(nMatrix)[2] != dim(Matrix)[2])){
    stop("function as.numeric.matrix: Error by converting!")
  }else{
    dimnames(nMatrix) <- dimnames(Matrix)
  }
  return(nMatrix)
}

###############################################################################
batchEffectRemoval <- function(Data,batchER){
  # remove batch effect from data
  # Author: Vadim Farztdinov
  # Version: 1.0.0, 09-Sep-2015.
  # Version: 1.0.1, 23-Aug-2016. batchF3 defined as numeric 
  # Version: 1.1.0, 30-Dec-2023. added ComBat
  # Arguments: 
  #     Data - Data object with three slots, MetaData, Matrix and PepDesc
  #     batchER - settings list with slots:
  #     batchER[["method"]]<- "LIMMA";           # currently only one option "LIMMA", reserved for future versions
  #     batchER[["batch"]]<-  "MS_run_Cat";      # for affine norm, c(centre,sd) (e.g. c(0,1)) or "mean" to centre over mean col centre and mean col sd 
  #     batchER[["batch2"]]<-  NULL;             # factor indicating batches (row name from input datamatrix), default=NULL
  #     batchER[["covariates"]]<-NULL;           # list of covariates to be adjusted for, default=NULL 
  #     batchER[["Class"]]<-NULL;                # optional treatment conditions factor to be preserved, default=NULL 
  # Output:
  #     Data - Data object with three slots, MetaData, Matrix and PepDesc with batch effect removed from Matrix
  batchF1=NULL; batchF2=NULL; batchF3=NULL; Class=NULL;
  
  DM<-Data$Matrix
  if (batchER[["method"]]=="LIMMA") {
    require(limma)
    # check input variables
    if (!is.null(batchER$batch))      batchF1<-factor(Data$MetaData[batchER$batch,])
    if (!is.null(batchER$batch2))     batchF2<-factor(Data$MetaData[batchER$batch2,])
    if (!is.null(batchER$covariates)) batchF3<-as.numeric(Data$MetaData[batchER$covariates,])

    # call removeBatchEffect function from LIMMA package
    if (is.null(batchER$Class)) {
      DM     <-removeBatchEffect(DM, batch=batchF1, batch2=batchF2, covariates=batchF3);
    } else {
      Class  <-factor(Data$MetaData[batchER$Class,])
      designM<-model.matrix(~0+Class)
      DM     <-removeBatchEffect(DM, batch=batchF1, batch2=batchF2, covariates=batchF3, design=designM);
    }
  } else if (batchER[["method"]]=="ComBat") {
    DM = sva::ComBat(dat=DM, batch=Data$MetaData[batchER$batch,], mod=NULL, par.prior=TRUE, prior.plots=FALSE)
  }
  # make sure that row and col names are preserved
  colnames(DM)<- colnames(Data$Matrix); rownames(DM)<- rownames(Data$Matrix)
  Data$Matrix<-DM
  return(Data)  # return corrected for batch effect Data object
} 

###############################################################################
makeExpressionSet<- function(DataObj) {
  # make Bioconductor ExpressionSet from Vadim Far DataObj 
  # Author: Vadim Farztdinov
  # Version: 1.0.0, 29-Jul-2017
  require(Biobase)
  Metadata<-createDataFrame(t(DataObj$MetaData))
  Featuredata<-createDataFrame(DataObj$PepDesc)
  expSet <- ExpressionSet(assayData = DataObj$Matrix,
                          phenoData = new("AnnotatedDataFrame", data = Metadata), 
                          featureData=new("AnnotatedDataFrame", data = Featuredata))
  return(expSet)
}
###############################################################################
makeMSnSet<- function(DataObj) {
  # make MSnbase MSnSet from Vadim Far DataObj 
  # Author: Vadim Farztdinov
  # Version: 1.0.0, 14-Dec-2020
  require(MSnbase)
  Metadata<-createDataFrame(t(DataObj$MetaData))
  Featuredata<-createDataFrame(DataObj$PepDesc)
  expSet <- MSnSet(exprs = DataObj$Matrix,
                          pData = new("AnnotatedDataFrame", data = Metadata), 
                          fData=new("AnnotatedDataFrame", data = Featuredata))
  return(expSet)
}
###############################################################################
annFilter<- function(Data,annFilterDesc){
  # removes unwanted features based on information form PepDesc slot
  # Author: Vadim Farztdinov
  # Version: 1.0.0, 07-May-2018
  # Version: 1.1.0, 24-Jan-2022, added calculation of absolute value
  # Version: 1.2.0, 22-Jan-2024, added filtering by pattern
  #
  # Arguments
  # Data     - structure with slots MetaData, Matrix, PepDesc
  # annFilterDesc - list with annFilter steps description
  # Value
  # Data     - reduced structure with slots MetaData, Matrix, PepDesc
  #
  # annFilterDesc[["var"]]     <- c("pMedE_biological","pPres_biological");   # c("Quan_Info","NN_Protein_Groups") factor for selestion of PSMs, e.g "Quan Usage"; Check character replacement
  # annFilterDesc[["expres"]]  <- c(">= 6.6", ">= 34")      # condition(s) to pass filter, e.g. c("<= 50","> 0")
  # annFilterDesc[["factor"]]  <- NULL;          # "Proteotypic", "Quan_Info","NN_Protein_Groups" factor for selestion of PSMs, e.g "Quan Usage"; Check character replacement
  # annFilterDesc[["value"]]   <- c("1")         #   , "1"  "Unique"
  # annFilterDesc$column         <- "Precursor.Id"    # vector of column names to filter by pattern
  # annFilterDesc$pattern        <- "(SILAC-.-H)"     # character vector of patterns of the same length as annFilter$column vector
  # annFilterDesc$negateP        <- TRUE         # logical vector, TRUE if pattern to exclude, FALSE to keep 
  
  if (!is.null(annFilterDesc$var)) {
    cat("Filtering by var expression started. Variables: \n" , file = "", sep = " ", append=TRUE)  
    cat(paste(annFilterDesc[["var"]], collapse = ", "),"\n" , file = "", sep = " ", append=TRUE)  
    cat("Expressions: ",paste(annFilterDesc[["expres"]], collapse = ", "),"\n" , file = "", sep = " ", append=TRUE)  
    for (ii in 1:length(annFilterDesc$var)) {
      annSelection<- as.numeric(Data$PepDesc[,annFilterDesc$var[ii]]) ; 
      if (isTRUE(annFilterDesc$abs[ii])) annSelection <- abs(annSelection);
      annSelection<- eval(parse(text= paste("annSelection", annFilterDesc$expres[ii], collapse = ""))); 
      annSelection[is.na(annSelection)]<- FALSE; 
      print(table(annSelection));    
      Data <- reduceFeatures(Data,annSelection)
    }
  } 
  if (!is.null(annFilterDesc$factor)) {
    cat("Filtering by factor value started. Factor(s): \n" , file = "", sep = " ", append=TRUE)  
    cat(paste(annFilterDesc$factor, collapse = ", "),"\n" , file = "", sep = " ", append=TRUE)  
    cat("Kept values: ",paste(unlist(annFilterDesc$value), collapse = ", "),"\n" , file = "", sep = " ", append=TRUE)  
    for (ii in 1:length(annFilterDesc$factor)) {
      annSelection<-trimws(Data$PepDesc[,annFilterDesc$factor[ii]],which="both") %in% annFilterDesc$value[[ii]]; 
      print(table(annSelection));    
      Data <- reduceFeatures(Data,annSelection)
    }
  }
  if (!is.null(annFilterDesc$column)) {
    cat("Filtering by pattern started. Column(s): \n" , file = "", sep = " ", append=TRUE)  
    cat(paste(annFilterDesc$column, collapse = ", "),"\n" , file = "", sep = " ", append=TRUE)  
    cat("Patterns: ",paste(unlist(annFilterDesc$pattern), collapse = ", "),"\n" , file = "", sep = " ", append=TRUE)  
    for (ii in 1:length(annFilterDesc$column)) {
      annSelection<-stringr::str_detect(Data$PepDesc[,annFilterDesc$column[ii]], annFilterDesc$pattern[ii], negate = annFilterDesc$negateP[ii])
      print(table(annSelection));    
      Data <- reduceFeatures(Data,annSelection)
    }
  }
  return(Data)  # return reduced Data object
}

###############################################################################
getDataObj <- function(DataDesc, verbose=FALSE) {
  # get DataObj from file or RData 
  # Author: Vadim Farztdinov
  # Version: 1.0.0, 12-04-2022
  # DataDesc - Data Description list object  
  # Output - DataObj with three slots: MetaData, Matrix and PepDesc

  if (is.null(DataDesc[["fileFormat"]])) DataDesc[["fileFormat"]]<-"txt";
  if (DataDesc[["fileFormat"]]=="txt") {
    # read the Data file  
    DataMatrix <- read.DataMatrix(fileName=file.path(DataDesc$dir,DataDesc$file), SEP=DataDesc$fsep, checkNames = FALSE)
    cat("Reading input DataMatrix finished. Size =",dim(DataMatrix), ", format = ",DataDesc$fmt ,"\n" , file = "", sep = " ")
    # prepare Data object. 
    cat("Creation of Data object started\n" , file = "", sep = " ")
    DataObj <- createDataObject(DataMatrix,DataDesc)
    n1<- dim(DataObj$Matrix)[1]; n2<- dim(DataObj$Matrix)[2];
    cat("Creation of Data object finished. DataObj$Matrix Size =",c(n1,n2), ". Top left corner is \n" , file = "", sep = " ")
    print(DataObj$Matrix[1:min(3,n1),1:min(3,n2)])
    rm(DataMatrix)
    
  } else if (DataDesc[["fileFormat"]]=="RData") {
    ObjName<-ifelse (!is.null(DataDesc[["Obj"]]), DataDesc[["Obj"]], "DataObj")
    tmp.env <- new.env() # create a temporary environment
    load(file = file.path(DataDesc$dir,DataDesc$file), envir=tmp.env)
    DataObj <- get(ObjName, pos=tmp.env) 
    rm(tmp.env)
    n1<- dim(DataObj$Matrix)[1]; n2<- dim(DataObj$Matrix)[2];
    cat("Reading of Data object finished. DataObj$Matrix Size = (",n1,",",n2,"). Top left corner is \n" , file = "", sep = "")
    print(DataObj$Matrix[1:min(3,n1),1:min(3,n2)])
    
  } else if (DataDesc[["fileFormat"]]=="WoSp") {
    ObjName<-ifelse (!is.null(DataDesc[["Obj"]]), DataDesc[["Obj"]], "LIMMAOut")
    tmp.env <- new.env() # create a temporary environment
    load(file = file.path(DataDesc$dir,DataDesc$file), envir=tmp.env)
    DataObj <- get(ObjName, pos=tmp.env) 
    rm(tmp.env)
    if (verbose) {
      cat("Reading of Object from Working Space object finished. Its structure is \n" , file = "", sep = "")
      print(str((DataObj)))
    }

  } else if (DataDesc[["fileFormat"]]=="xlsx") {
    library(readxl)
    DataMatrix <- read_excel(path = file.path(DataDesc$dir,DataDesc$file), 
                             sheet = DataDesc[["sheet"]])
    cat("Reading input DataMatrix finished. Size =",dim(DataMatrix), ", format = ",DataDesc$fmt ,"\n" , file = "", sep = " ")
    # prepare Data object. 
    cat("Creation of Data object started\n" , file = "", sep = " ")
    DataObj <- createDataObject(DataMatrix,DataDesc)
    n1<- dim(DataObj$Matrix)[1]; n2<- dim(DataObj$Matrix)[2];
    cat("Creation of Data object finished. DataObj$Matrix Size =",c(n1,n2), ". Top left corner is \n" , file = "", sep = " ")
    print(DataObj$Matrix[1:min(3,n1),1:min(3,n2)])
    rm(DataMatrix)
  }
  return(DataObj)
}
###############################################################################
imputationTool<- function(Data,WorkDesc,Output=NULL){
  # impute NAs in Data. Data is either list with slot Matrix or numerical matrix  
  # Author: Vadim Farztdinov
  # Version: 1.0.0, 09-Sep-2013
  # Version: 1.1.0, 20-Apr-2015
  # Version: 1.1.1, 25-Jun-2015. Replaced removeNAs$group by impute$group
  # Version: 1.2.1, 21-Aug-2016. New imputation method added - iterative PCA 
  # Version: 1.2.2, 15-Oct-2016. Added warnings when imputation is skipped 
  # Version: 1.3.0, 14-Dec-2020. Added group based iPCA imputation. Added methods  
  # Version: 1.4.0, 14-Dec-2020. Added MSnbase methods "MLE","bpca", "QRILC", "MinDet", "MinProb", "min" 
  # Version: 1.4.1, 19-Jan-2021. Added impute$zeros 
  # Version: 2.0.0, 14-Nov-2022. New design of ImputeData function. Introduced iterative 
  #                              nested group and sliding window based imputation
  # Version: 2.1.0, 01-Dec-2022. Added imputation of OnOff features
  # Version: 2.1.1, 11-Jan-2023. Update treating PepDesc for case it is not a df
  # Version: 2.1.2, 19-Jul-2023. Added printing the number of OnOff features left for imputation
  # Version: 2.1.3, 21-Jul-2023. Updated impute[["zeros"]] - replaced by impute nonpositive values
  # Version: 2.1.4, 03-Nov-2023. Introduced group dependent imputation method 
  # Arguments
  # Data - data object with slot Data$Matrix containing numerical matrix
  # DataDesc - description of the data. Currently not in use, reserved for future versions
  # WorkDesc - list with slots:
  # WorkDesc[["impute"]]          <- vector("list");  # settings for imputation
  # WorkDesc$impute[["run"]]      <- TRUE;            # TRUE|FALSE
  # WorkDesc$impute[["groupBased"]]<- TRUE;           # TRUE|FALSE, imputation over group members only or over all samples
  # WorkDesc$impute$afterBatchER  <- FALSE;           # TRUE|FALSE = NULL - default if TRUE imputation is carried out after final subset selection
  # WorkDesc$impute[["zeros"]]    <- TRUE;            # TRUE|FALSE- default=NULL - replace zeros by NA for imputation
  # WorkDesc$impute[["Th_pre"]] <- 34;                 # % of NA when feature is considered present
  # WorkDesc$impute[["Th_abs"]] <- 66;                 # % of NA when feature is considered absent
  
  # WorkDesc$impute$level[[1]]        <- vector("list");  # settings for imputation
  # WorkDesc$impute$level[[1]][["group"]]  <- DataDesc$row$SubProject  # DataDesc$row$SubClass;  # default is DataDesc$row$Class
  # WorkDesc$impute$level[[1]][["method"]] <- "iPCA";          # "iPCA" - recommended, "knn", "mean", "median","MLE","bpca", "QRILC", "MinDet", "MinProb", "min"
  # WorkDesc$impute$level[[1]][["PctNA"]]  <-34;        # % of allowed NA
  # WorkDesc$impute$level[[1]][["ncpMax"]] <- 5;               # if (method=="iPCA") default = 5
  # WorkDesc$impute$level[[1]][["nn"]]     <- 3;               # if method=="knn"
  
  # WorkDesc$impute$level[[2]]        <- vector("list");  # settings for imputation
  # WorkDesc$impute$level[[2]][["group"]]  <- DataDesc$row$SampleGr;  # default is DataDesc$row$Class
  # WorkDesc$impute$level[[2]][["method"]] <- "iPCA";          # "iPCA" - recommended, "knn", "mean", "median","MLE","bpca", "QRILC", "MinDet", "MinProb", "min"
  # WorkDesc$impute$level[[2]][["PctNA"]]  <-34;        # % of allowed NA
  # WorkDesc$impute$level[[2]][["ncpMax"]] <- 5;               # if (method=="iPCA") default = 5
  # WorkDesc$impute$level[[2]][["nn"]]     <- 3;               # if method=="knn"
  
  # WorkDesc$impute[["fullSet"]]  <- vector("list");  # settings for imputation
  # WorkDesc$impute$fullSet[["method"]] <- "iPCA";          # "iPCA" - recommended, "knn", "mean", "median","MLE","bpca", "QRILC", "MinDet", "MinProb", "min"
  # WorkDesc$impute$fullSet[["PctNA"]]  <-34;               # % of allowed NA
  # WorkDesc$impute$fullSet[["ncpMax"]] <- 5;               # if (method=="iPCA") default = 5
  # WorkDesc$impute$fullSet[["nn"]]     <- 3;               # if method=="knn"
  
  # WorkDesc$impute[["OnOff"]]    <- TRUE;            # TRUE|FALSE- default=NULL - impute OnOff features
  # WorkDesc$impute[["plotAfter"]]<- TRUE             # TRUE- default|FALSE
  #
  # Output Value
  # Data - data object with slot Data$Matrix containing imputed numerical matrix; 
  #        the observed values are kept for the non-missing entries and the 
  #        missing values are replaced by the predicted ones. 
  
  if (is.list(Data)) {
    if (!is.double(Data$Matrix)) {
      NR=nrow(Data$Matrix); DM<- matrix(as.numeric(as.matrix(Data$Matrix)), nrow=NR)  
    } else {
      DM <- Data$Matrix
    }
  } else DM<-Data;
  
  DM[!is.finite(DM)]<-NA;
  NFeatures<-nrow(DM); NSamples<- ncol(DM);
  
  if(isTRUE(WorkDesc$impute[["zeros"]])) {
    cat(sprintf("%i features were found in DataMatrix with zero/negative values. Replaced by NAs for imputation\n",length(which(DM<=0))) , file = "", sep = " ")
    DM[DM<=0] <- NA
  }
  cat(sprintf("Total fraction of NAs found in DataMatrix = %.1f %s \n",100*sum(is.na(DM))/(NFeatures*NSamples), "%" ) , file = "", sep = " ")

  if (!WorkDesc$impute$groupBased) {
    impSettings <- WorkDesc$impute$fullSet
    cat(sprintf("Total set based missing values imputation using %s method\n",impSettings$method) , file = "", sep = " ")
    # get pct missing
    NApctTot<- 100*rowSums(is.na(DM))/ncol(DM)
    imputeFeatures <- NApctTot <= impSettings[["PctNA"]]
    fDM <- DM[imputeFeatures,]
    fDM <- imputeDM(fDM,impSettings)
    
    if(isTRUE(WorkDesc$impute[["zeros"]])) {
      cat(sprintf("Checking DataMatrix after imputation for zero/negative values resulted in %i such features\n",length(which(fDM<=0))) , file = "", sep = " ")
      cat(sprintf("They were replaced by min(fDM[fDM>0]) = %f \n", min(fDM[fDM>0])) , file = "", sep = " ")
      fDM[fDM<=0] <- min(fDM[fDM>0])
    }
    
    # check if input is DataObject
    if (is.list(Data)) {
      Data$PepDesc   <- Data$PepDesc[imputeFeatures, ];
      Data$Matrix    <- fDM;
    } else {Data<-fDM}
    
  } else {   # group based imputation; only for input DataObjects
    for (kk in 1:length(WorkDesc$impute$level)) { 
      impSettings <- WorkDesc$impute$level[[kk]]
      
      sampleGroup<-Data$MetaData[impSettings$group, ]
      fq<-table(sampleGroup);
      NG<-length(fq);GrN<-dimnames(fq);
      cat(sprintf("Group = %s based missing values imputation using %s method(s)\n   Group names = %s\n",
                  impSettings$group, paste0(impSettings$method,collapse=","), paste0(GrN[[1]],collapse=",")), file = "", sep = " ")
      if (all(impSettings$method %in% c("knn","mean","median","iPCA", "MLE","bpca", "QRILC", "MinDet", "MinProb", "min"))) {
        
        for (ii in 1:NG) {
          DMii<-DM[ , sampleGroup==GrN[[1]][ii] ];
          
          # get pct missing
          NApct<- 100*rowSums(is.na(DMii))/ncol(DMii)
          imputeFeatures <- NApct <= impSettings[["PctNA"]]
          fDMii <- DMii[imputeFeatures,]
          
          # get imputation method
          impSettings_ii <- impSettings
          impSettings_ii$method <- ifelse(length(impSettings$method)==NG, impSettings$method[ii],impSettings$method[1])
          if (length(impSettings_ii$method) > 1) stop("length of imputation method should be 1")
          #
          fDMii <- imputeDM(fDMii,impSettings_ii)
          #
          DMii[imputeFeatures,] <- fDMii
          DM[,sampleGroup==GrN[[1]][ii] ]<-DMii;  # backward assignment
          gc()
        }
      } else {
        cat(sprintf("%s method is not used for group based NA imputation. Imputation skipped\n",impSettings$method) , file = "", sep = " ")
      }
    }
    
    impSettings <- WorkDesc$impute$fullSet
    cat(sprintf("Final step of group based imputation.\n") , file = "", sep = " ")
    
    # WorkDesc$impute[["OnOff"]] 
    if (isTRUE(WorkDesc$impute$OnOff)) {
      cat(sprintf("Imputation of OnOff features in groups %s \nas left-censored missing data using %s\n", WorkDesc$analyseNAs$group,"random draws from Gaussian distribution") , file = "", sep = " ")
      if (is.null(WorkDesc$impute$BgDist)) {
        mnarDist <- getMissingDistPar(DM, WorkDesc$impute)
      } else{
        mnarDist <- WorkDesc$impute$BgDist
      }
      # cat(sprintf("Centre of mnar distribution = %4.2f and sigma = %4.2f\n", mnarDist$medBg,mnarDist$madBg) , file = "", sep = " ")
      cat("Number of OnOff features left for imputation =", sum(Data$PepDesc[,"is_OnOff"]),"\n", file = "", sep = " ")
      
      imputeFeatures <- as.logical(Data$PepDesc[,"is_OnOff"])
      fDMii <- DM[imputeFeatures,]
      fDMii <- imputeMNAR(fDMii, center = mnarDist$medBg, sigma = mnarDist$madBg) 
      #
      DM[imputeFeatures,] <- fDMii
    }
    # check missingness after all steps
    cat(sprintf("Total set based missing values imputation using %s method\n",impSettings$method) , file = "", sep = " ")
    # get pct missing
    NApctTot<- 100*rowSums(is.na(DM))/ncol(DM)
    imputeFeatures <- NApctTot <= impSettings[["PctNA"]]
    fDM <- DM[imputeFeatures,]
    fDM <- imputeDM(fDM,impSettings)
    
    if(isTRUE(WorkDesc$impute[["zeros"]])) {
      cat(sprintf("Checking DataMatrix after imputation for zero/negative values resulted in %i such features\n",length(which(fDM<=0))) , file = "", sep = " ")
      cat(sprintf("They were replaced by min(fDM[fDM>0]) = %f \n", min(fDM[fDM>0])) , file = "", sep = " ")
      fDM[fDM<=0] <- min(fDM[fDM>0])
    }
    
    Data$PepDesc   <- Data$PepDesc[imputeFeatures,];
    Data$Matrix     <- fDM;
  } 
  
  return(Data)
}

###############################################################################
imputeDM<- function(Data,how2impute,...){
  # impute NAs in Data. Data is either list with slot Matrix or numerical matrix
  # This is simplified version of old function imputeData. To extend functionality
  # MSnbase::impute function is added
  # Author: Vadim Farztdinov
  # Version: 1.0.0, 14-Nov-2022
  # Version: 1.0.1, 15-Nov-2022  Added methods from pcaMethods package: bpca and llsImpute
  # Version: 1.0.2, 22-Nov-2022  Added diagnostic of missingness in matrix
  # Version: 1.0.3, 19-Jul-2023  Updated default parameters for knn
  #
  # Arguments
  # Data - data object with slot Data$Matrix containing numerical matrix
  # how2impute - list with slots:
  # how2impute[["method"]]   <- "iPCA";          # "iPCA" - recommended, "knn", "mean", "median","MLE","QRILC", "MinDet", "MinProb", "min"
  # how2impute[["ncpMax"]]   <- 5;               # if (method=="iPCA") default = 5
  # how2impute[["nn"]]       <- 3;               # if method=="knn"
  # how2impute[["groupBased"]]<- TRUE;           # TRUE|FALSE, imputation over group members only or over all samples
  # how2impute[["group"]]    <- DataDesc$row$SubProject;  # default is DataDesc$row$Class
  # how2impute[["plotAfter"]]<- TRUE             # TRUE- default|FALSE
  #
  # Output Value
  # Data - data object with slot Data$Matrix containing imputed numerical matrix; 
  #        the observed values are kept for the non-missing entries and the 
  #        missing values are replaced by the predicted ones. 
  
  if (is.list(Data)) {
    if (!is.double(Data$Matrix)) {
      NR=nrow(Data$Matrix); DM<- matrix(as.numeric(as.matrix(Data$Matrix)), nrow=NR)  
    } else {DM <- Data$Matrix}
  } else DM<-Data;
  
  cat(sprintf("%.1f %s of NAs in DM. Imputation using %s method\n",100*sum(is.na(DM))/length(DM), "%", how2impute$method) , file = "", sep = " ")

  if (how2impute$method == "knn") {
    iD <- impute::impute.knn(data=DM, k=how2impute$nn, rowmax = 0.99, colmax = 0.99, maxp = 4500)
    DM<-iD$data;
  } else if (how2impute$method == "iPCA") {
    suppressPackageStartupMessages(library(missMDA))
    if (is.null(how2impute$ncpMax)) {how2impute$ncpMax<-min(dim(DM)-1,5)}
    if (is.null(how2impute$scaleDM)) {how2impute$scaleDM<-TRUE}
    if (is.null(how2impute$threshold)) {how2impute$threshold<-1e-4}
    ncpMax<-how2impute$ncpMax
    scaleDM<-how2impute$scaleDM
    thresh<-how2impute$threshold
    res <- pcaMethods::prep(t(DM), scale="uv", center=TRUE, simple=FALSE)
    ncomp <- estim_ncpPCA(as.data.frame(res$data),ncp.min=0, ncp.max=ncpMax,  
                          scale = scaleDM, threshold = thresh)
    iD <- imputePCA(as.data.frame(res$data), ncp = ncomp$ncp, scale = scaleDM)
    backscaled <- pcaMethods::prep(iD$completeObs, scale=res$scale, center=res$center, rev=TRUE)
    DM <- data.matrix(t(backscaled))
    detach("package:missMDA", character.only = TRUE)
    
  } else if (how2impute$method == "bpca") {
    suppressPackageStartupMessages(library(pcaMethods))
    res <- prep(t(DM), scale="uv", center=TRUE, simple=FALSE)
    resultBPCA <- pca(res$data, method = "bpca", 
                      nPcs = how2impute$ncpMax, verbose = FALSE)
    cmb <- completeObs(resultBPCA)
    DM <- t(prep(cmb, scale=res$scale, center=res$center, rev=TRUE))
    detach("package:pcaMethods", character.only = TRUE)
    
  } else if (how2impute$method == "llsImpute") {
    suppressPackageStartupMessages(library(pcaMethods))
    result <- llsImpute(t(DM), k = how2impute$nn, correlation="pearson", allVariables=TRUE)
    DM <- t(completeObs(result))
    # alternaive way using wrapper function nni (recommended by developer)
    # res2 <- nni(t(DM), method = "llsImpute", k=how2impute$nn, allVariables=TRUE)
    # DM <- t(completeObs(res2))
    detach("package:pcaMethods", character.only = TRUE)
    
  } else if (how2impute$method %in% c("QRILC", "MinDet", "MinProb", "min")) { # "bpca" is called directly from pcaMethods package 
    suppressPackageStartupMessages(library(MSnbase))
    naset <- new("MSnSet", exprs=DM)
    iD <- impute(naset, method=how2impute$method)
    DM<-exprs(iD);
    detach("package:MSnbase", character.only = TRUE)
  } else {
    if (how2impute$method == "mean") {
      DM <- ifelse(is.na(DM), rowMeans(DM, na.rm=TRUE), unlist(DM))
    } else if (how2impute$method == "median") {
      suppressPackageStartupMessages(library(matrixStats))
      DM <- ifelse(is.na(DM), rowMedians(DM, na.rm=TRUE), unlist(DM))      
      detach("package:matrixStats", character.only = TRUE)
    } else {
      cat(sprintf("%s method is not used for total set based NA imputation. Imputation skipped\n",how2impute$method) , file = "", sep = " ")
    }
  }
  if (is.list(Data)) Data$Matrix <-(DM) else Data<-(DM);
  return(Data)
}

###############################################################################
getMissingDistPar <- function(DM,impute) {
  # get propeties of missing peptides distribution
  # Author: Vadim Farztdinov
  # Version: 1.0.0, 14-Nov-2022
  # Version: 1.0.1, 02-Dec-2022 Update of mAbs, madA, and mBac definition 
  Th_pre = 34; # features % of missing values threshold below which feature is considered to be present
  Th_abs = 66; # features % of missing values threshold above which feature is considered to be absent
  if (!is.null(impute[["Th_pre"]])) {Th_pre = impute[["Th_pre"]]}
  if (!is.null(impute[["Th_abs"]])) {Th_abs = impute[["Th_abs"]]}
  NSamples<- ncol(DM);
  # analyze features by absence in the total set
  if (is.data.frame(DM)) {NApct<- 100*rowSums(is.na.data.frame(DM))/NSamples;  
  } else {NApct<- 100*rowSums(is.na(DM))/NSamples}
  
  MedianE <- matrixStats::rowMedians(DM, na.rm=T)
  mPre <- median(MedianE[NApct<=Th_pre], na.rm=T)                       # median exp of present
  madP <- mad(MedianE[NApct<=Th_pre], na.rm=T)                       # mad of present
  mAbs <- median(MedianE[(NApct >Th_abs) & (MedianE < mPre-madP)], na.rm=T) # median exp of absent
  madA <- mad(MedianE[(NApct >Th_abs) & (MedianE < mPre-madP)], na.rm=T)    # mad of not absent
  madB <- madA/2                                                      # mad of backgroung
  mBac <- mAbs-2*madB                                                   # median exp of backgroung
  return(list(medBg=mBac,madBg=madB,NApct=NApct))
}

###############################################################################
imputeMNAR <- function (DM, center = NULL, sigma = NULL) 
{
  # Performs the imputation of left-censored missing data by random draws from a
  # Gaussian distribution centered in a minimal value. Considering a
  # peptide/protein expression data matrix with n columns corresponding to
  # biological samples and p lines corresponding to peptides/proteins, for each
  # sample (column), the mean value of the Gaussian distribution is set to a
  # minimal value observed in that sample. The minimal value observed is estimated
  # as being the q-th quantile (e.g. q = 0.01) of the observed values in that
  # sample. The standard deviation is estimated as the median of the
  # peptide/protein-wise standard deviations. Note that when estimating the
  # standard deviation of the Gaussian distribution, only the peptides/proteins
  # which present more than 50% recorded values are considered.
  #
  # Author: Cosmin Lazar
  # Modified: Vadim Farztdinov
  # Version: 1.0.0, 12-Nov-2022
  # Version: 1.0.1, 01-Dec-2022  Added diagnostic of number of NÁs in matrix
  # Version: 1.0.2, 02-Dec-2022  Added set.seed for reproducibility
  # Version: 1.0.3, 19-Jul-2023  Added printing of parameters and evalMethod into log file
  
  nSamples = dim(DM)[2]
  nFeatures = dim(DM)[1]
  dataSet.imputed = DM
  cat(sprintf("%.1f %s of NAs in DM. Imputation using %s function\n",100*sum(is.na(DM))/length(DM), "%", "imputeMNAR") , file = "", sep = " ")
  
  evalMethod = "VF getMissingDistPar"
  # get center of missing values distribution
  if (is.null(center) || is.na(center)) {
    center = median(apply(dataSet.imputed, 2, quantile, prob = 0.01, na.rm = T))
    evalMethod = "Cosmin Lazar"
  }
  # get the width on missing values distribution
  if (is.null(sigma) || is.na(sigma)) {
    count.NAs = apply(!is.na(DM), 1, sum)
    count.NAs = count.NAs/nSamples
    dataSet.filtered = DM[which(count.NAs > 0.5), ]
    protSD = apply(dataSet.filtered, 1, sd)
    sigma = median(protSD, na.rm = T)
  }
  
  cat(sprintf("Centre of mnar distribution = %4.2f and sigma = %4.2f\n", center, sigma) , file = "", sep = " ")
  cat("defined by", evalMethod, "evaluation method\n", file = "", sep = " ")
  
  set.seed(12345)
  for (i in 1:(nSamples)) {
    # impute for each sample separately
    dataSet.to.impute.temp = abs(rnorm(nFeatures, mean = center, sd = sigma))
    dataSet.imputed[which(is.na(DM[, i])), i] = 
      dataSet.to.impute.temp[which(is.na(DM[,i]))]
  }
  return(dataSet.imputed)
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