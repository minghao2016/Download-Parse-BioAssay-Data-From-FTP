
require(data.table)
require(parallel)
require(stringi)
require(RCurl)


.ls.objects <- function (pos = 1, pattern, order.by,
                         decreasing = FALSE, 
						 head = FALSE, n = 5) {
  # improved list of objects
  napply <- function(names, fn) sapply(names, function(x)
    fn(get(x, pos = pos)))
  names <- ls(pos = pos, pattern = pattern)
  obj.class <- napply(names, function(x) as.character(class(x))[1])
  obj.mode <- napply(names, mode)
  obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
  obj.prettysize <- napply(names, function(x) {
    format(utils::object.size(x), units = "auto") })
  obj.size <- napply(names, object.size)
  obj.dim <- t(napply(names, function(x)
    as.numeric(dim(x))[1:2]))
  vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
  obj.dim[vec, 1] <- napply(names, length)[vec]
  out <- data.frame(obj.type, obj.size, obj.prettysize, obj.dim)
  names(out) <- c("Type", "Size", "PrettySize", "Length/Rows", "Columns")
  if (!missing(order.by))
    out <- out[order(out[[order.by]], decreasing=decreasing), ]
  if (head)
    out <- head(out, n)
  out
}


lsos <- function(..., n = 10) {
  # shorthand
  .ls.objects(..., order.by = "Size", decreasing = TRUE, head = TRUE, n = n)
}


lsdiskfiles <- function() {
	# purpose 
	#   list disk files, like `ll -h`
	# input
	#   NULL
	# output
	#   list files with size
	require(data.table)
	require(gdata)
	FI <- file.info(list.files())
	sizeInfo <- FI$size
	sizeInfoOk <- humanReadable(sizeInfo, width = 4)
	FIshow <- cbind(sizeInfoOk, FI)
	print(FIshow)
	return(invisible(NULL))
}


getPubChemInfo <- function(curAidFnFp) {
  # purpose 
  #   parse pubchem bioassay csv files from FTP
  # input 
  #   curAidFnFp: string, current file name full path, such as: "/mnt/LHRI/Bioinformatics/Projects/DAVIDUpdate/pubchem3/data/1160654.csv.gz"
  # output 
  #   data.table, parsed result

  require(data.table)
  require(stringi)


  # fnDir: file name directory
  # fileFullPath <- paste0(fnDir, "/", curFn)
  aid <- gsub(".csv.gz", "", basename(curAidFnFp))
  ###########################
  # showProgress = TRUE, default TRUE
  dat <- fread(curAidFnFp, showProgress = TRUE)
  ###########################
  colNames <- colnames(dat)
  
  if (all(is.na(dat$PUBCHEM_SID))) {
    # no compound information at all
	dtFnl <- data.table(AID = aid, SID = "e", CID = "e", OUTCOME = "e", GI = "e", GENEID = "e", TAG = "e")
	
	return(dtFnl)
	# need [savedFn] in the beginning if turn on below
	#fwrite(dtFnl, file = savedFn, append = TRUE)
    #return(invisible(NULL))
  }
  
  panelAssay <- grep("panel", dat[, PUBCHEM_RESULT_TAG], value = TRUE, ignore.case = TRUE)
  
  if (length(panelAssay) != 0) {
    # "p", panel assay
    tag <- "p"
  } else {
    resultTypeValue <- dat[PUBCHEM_RESULT_TAG %in% "RESULT_TYPE", ]
	multiTargetIdx <- which(grepl("ncbi_gene_id|ncbi_protein_id", resultTypeValue, ignore.case = TRUE))
	if (length(multiTargetIdx) != 0) {
	  # "m", multi-target assay
	  tag <- "m"
	} else {
	  # "s", single-target assay
	  tag <- "s"
	}
  }
  
  # single taget assay case
  if (tag == "s") {
    baseColNames <- c("PUBCHEM_SID", "PUBCHEM_CID", "PUBCHEM_ACTIVITY_OUTCOME")
	sidcidoutcome <- dat[which(PUBCHEM_RESULT_TAG == 1):.N, baseColNames, with = FALSE]
	dtFnl <- data.table(AID = aid, sidcidoutcome, GI = "e", GENEID = "e", TAG = "s")
	setnames(dtFnl, c("AID", "SID", "CID", "OUTCOME", "GI", "GENEID", "TAG"))
	
	return(dtFnl)
    #fwrite(dtFnl, file = savedFn, append = TRUE)
	#return(invisible(NULL))
  }
  
  
  # multiple target assay case
  if (tag == "m") {
    tmpTargetId <- dat[PUBCHEM_RESULT_TAG %in% "RESULT_TYPE", multiTargetIdx, with = FALSE]
	tmpTargetIdVec <- unname(unlist(tmpTargetId))
	isGi <- grepl("protein", tmpTargetIdVec, ignore.case = TRUE)
	if (isGi) {
	  # gi case
	  # aid = 1811
	  baseColNames <- c("PUBCHEM_SID", "PUBCHEM_CID", "PUBCHEM_ACTIVITY_OUTCOME")
	  sidcidoutcome <- dat[which(PUBCHEM_RESULT_TAG == 1):.N, baseColNames, with = FALSE]
	  targetId <- dat[which(PUBCHEM_RESULT_TAG == 1):.N, multiTargetIdx, with = FALSE]
	  dtFnl <- data.table(AID = aid, sidcidoutcome, GI = targetId, GENEID = "e", TAG = "m")
	  setnames(dtFnl, c("AID", "SID", "CID", "OUTCOME", "GI", "GENEID", "TAG"))
	  
	  return(dtFnl)
	  #fwrite(dtFnl, file = savedFn, append = TRUE)
	  #return(invisible(NULL))
	} else {
	  # gene_id case
	  # aid = 651811
	  baseColNames <- c("PUBCHEM_SID", "PUBCHEM_CID", "PUBCHEM_ACTIVITY_OUTCOME")
	  sidcidoutcome <- dat[which(PUBCHEM_RESULT_TAG == 1):.N, baseColNames, with = FALSE]
	  targetId <- dat[which(PUBCHEM_RESULT_TAG == 1):.N, multiTargetIdx, with = FALSE]
	  dtFnl <- data.table(AID = aid, sidcidoutcome, GI = "e", GENEID = targetId, TAG = "m")
	  setnames(dtFnl, c("AID", "SID", "CID", "OUTCOME", "GI", "GENEID", "TAG"))
	  
	  return(dtFnl)
      #fwrite(dtFnl, file = savedFn, append = TRUE)
	  #return(invisible(NULL))
	}
  }

  
  # panel assay case
  if (tag == "p") {
    # select rows by index
    if (any(grepl("target", panelAssay, ignore.case = TRUE))) {
      # gi or geneid
	  panelAssayTarget <- grep("target", panelAssay, value = TRUE, ignore.case = TRUE)
      panelAssayTargetId <- grep("id", panelAssayTarget, value = TRUE, ignore.case = TRUE)
      panelRow1stIdx <- dat[, which(PUBCHEM_RESULT_TAG %in% panelAssayTargetId)]
      compdStartIdx <- dat[, which(PUBCHEM_RESULT_TAG == 1)]
      selRowIdx <- c(panelRow1stIdx, compdStartIdx:nrow(dat))
    } else {
      # panel assay without target, but [cell line] or others
      # here use name, instead of id
      panelAssayName <- grep("name", panelAssay, value = TRUE, ignore.case = TRUE)
      if (length(panelAssayName) != 0) {
        # first row is panel_name
        panelRow1stIdx <- dat[, which(PUBCHEM_RESULT_TAG %in% panelAssayName)]
        compdStartIdx <- dat[, which(PUBCHEM_RESULT_TAG == 1)]
        selRowIdx <- c(panelRow1stIdx, compdStartIdx:nrow(dat))
      } else {
        # aid = 743407
        # no panel_name case, just have panel_id
        panelAssayId <- grep("id", panelAssay, value = TRUE, ignore.case = TRUE)
        panelRow1stIdx <- dat[, which(PUBCHEM_RESULT_TAG %in% panelAssayId)]
        compdStartIdx <- dat[, which(PUBCHEM_RESULT_TAG == 1)]
        selRowIdx <- c(panelRow1stIdx, compdStartIdx:nrow(dat))
      }
    }
    
    # select columns by names    
    panelOutcomePat <- "outcome"
    panelOutcomePatValue <- grep(panelOutcomePat, colNames, value = TRUE, ignore.case = TRUE)
    panelOutcomePatValue <- setdiff(panelOutcomePatValue, "PUBCHEM_ACTIVITY_OUTCOME")
    # normally, target_id and outcome in the same column
    selColNames <- c("PUBCHEM_SID", "PUBCHEM_CID", panelOutcomePatValue)
    
    
    # rare case
    if (length(panelOutcomePatValue) == 0) {
      # no panel outcome column at all
      # aid = 743028
      # cur_fn = "743028.csv.gz"
      # cur_fn = "652108.csv.gz"
      # cur_fn = "652181.csv.gz"
      # cur_fn = "434993.csv.gz"
      dtFnl <- data.table(AID = aid, 
                          SID = dat[which(PUBCHEM_RESULT_TAG == 1):.N, PUBCHEM_SID], 
                          CID = dat[which(PUBCHEM_RESULT_TAG == 1):.N, PUBCHEM_CID],
                          OUTCOME = "e", 
						  GI = "e", 
						  GENEID = "e", 
						  TAG = "p")
      
      return(dtFnl)
      #fwrite(dtFnl, file = savedFn, append = TRUE)
	  #return(invisible(NULL))
    }
    
    # sometimes, target_id is not in the same [outcome] column, but the first column with same pre-column without [outcome] part
    # then find the first part and extract the target_id
    # like: abc_ki, abc_outcome
    # like: abc ki, abc (outcome); cur_fn = "434974.csv.gz"
    # like: outcome_abc, ki_abc; cur_fn = "651549.csv.gz", header of csv is wrong, just manually add target_id
        
    outcomePatRm <- c("(outcome)", "_outcome", "outcome ", " outcome")
    # escape special character
    outcomePatRm <- gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", outcomePatRm)
    outcomePatRmComb <- paste(outcomePatRm, collapse = "|")
    # remove [outcome] part
    patNameNoOutcome <- gsub(outcomePatRmComb, "", panelOutcomePatValue, ignore.case = TRUE)
    patNameNoOutcome <- gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", patNameNoOutcome)
    patNameNoOutcomeComb <- paste(patNameNoOutcome, collapse = "|")
    allOutcomeRelatedNames <- grep(patNameNoOutcomeComb, colNames, value = TRUE, ignore.case = TRUE)
    # in case, panel outcome name is just [outcome] itself
    allOutcomeRelatedNames <- setdiff(allOutcomeRelatedNames, "PUBCHEM_ACTIVITY_OUTCOME")
    
    if (length(allOutcomeRelatedNames) == 0) {
      # just assign original [outcome] column index
      targetIdColNameIdx <- which(colNames %in% panelOutcomePatValue)
    }
    
    if (length(allOutcomeRelatedNames) == 1) {
      # only one panel outcome column exists
      targetIdColNameIdx <- which(colNames %in% panelOutcomePatValue)
    }
    
    if (length(allOutcomeRelatedNames) > 1) {
      # aid = 2393
      outcomeIdx <- grep("outcome", allOutcomeRelatedNames, ignore.case = TRUE)
      
      if (any(outcomeIdx == 1)) {
        # [outcome] in the 1st place
        # suppose, element with outcome should be at the end, not in the front
        targetIdColNameIdx <- which(colNames %in% panelOutcomePatValue)
      }
      
      # 1st method: confirm by position, but failed if [outcome] in the middle
      # like: a_ki, a_outcome, a_zz
      # aid = 2393
      # so comment it for new
      # if (all(outcome_idx != 1)) {
      #   # suppose, it is the normal case
      #   # remove the last idx, the use this as basis, and the next
      #   # element should be the 1st element of corresponding [outcome]
      #   outcome_idx_rm_end <- outcome_idx[1:(length(outcome_idx) - 1)]
      #   target_id_col_names <- all_outcome_related_names[c(1, outcome_idx_rm_end + 1)]  
      #   target_id_col_name_idx <- which(col_names %in% target_id_col_names)
      # }
      
      # 2nd method: good for [outcome] in the middle and in the end case
      # aid = 2393
      if (all(outcomeIdx != 1)) {
        numChar <- nchar(patNameNoOutcome)
        # descending order
        # longest to shortes for pattern, in order to extract correct pattern in the long string
        patNameNoOutcomeLtoS <- patNameNoOutcome[order(-numChar)]
        patNameNoOutcomeLtoScomb <- paste(patNameNoOutcomeLtoS, collapse = "|")
        # temporary patterns
        tmpPat <- stringi::stri_extract_first_regex(allOutcomeRelatedNames, patNameNoOutcomeLtoScomb)
        tmpdt <- data.table(allOutcomeRelatedNames, tmpPat, num = 1:length(allOutcomeRelatedNames))
        tmpGoodOrd <- tmpdt[, min(num), by = tmpPat]
        goodIdx <- tmpGoodOrd[, V1]
        targetIdColNames <- allOutcomeRelatedNames[goodIdx]
        targetIdColNameIdx <- which(colNames %in% targetIdColNames)
      }
    }
    
   
    # check panel assay type is gi or geneid, if not gi, then all put into geneid
    # special aid = 624278, have gene_id and gi simutaneously
    if (any(grepl("target", panelAssay, ignore.case = TRUE))) {
      panelTargetType <- grep("type", panelAssay, value = TRUE, ignore.case = TRUE)
      panelTargetTypeValue <- unique(unlist(dat[PUBCHEM_RESULT_TAG == panelTargetType, ]))
      panelTargetTypeValue <- setdiff(panelTargetTypeValue, panelTargetType)
      panelTargetTypeValue <- panelTargetTypeValue[panelTargetTypeValue != ""]
      panelTargetTypeValue <- panelTargetTypeValue[!is.na(panelTargetTypeValue)]
      # sometimes, it includes both gi and geneid, like below
      isGi <- grepl("protein|gi", panelTargetTypeValue, ignore.case = TRUE)
    } else {
      isGi <- FALSE
    }
    
    if (length(isGi) > 1) {
      # include both geneid and gi
      # cur_fn <- "624278.csv.gz"
      #stop("is_gi should be length of 1 \n")
      dtFnl <- NULL
      for (j in 1:length(panelTargetTypeValue)) {
        curTargetType <- panelTargetTypeValue[j]
        curTargetTypeVec <- unname(unlist(dat[PUBCHEM_RESULT_TAG == panelTargetType, ]))
        curTargetTypeIdx <- which(curTargetTypeVec %in% curTargetType)
        curSelColNames <- c("PUBCHEM_SID", "PUBCHEM_CID", colNames[curTargetTypeIdx])
        
		############################################################
        curDt <- dat[selRowIdx, curSelColNames, with = FALSE]
        ############################################################
        
		isGiCur <- grepl("protein|gi", curTargetType, ignore.case = TRUE)
        
        targetId <- unname(unlist(curDt[1, 3:ncol(curDt)]))
        
        sidcid <- curDt[2:.N, 1:2]
        panelRes <- curDt[2:.N, 3:ncol(curDt)]
        panelRes <- c(t(panelRes))
        
        sidcidstack <- sidcid[rep(seq_len(nrow(sidcid)), each = length(targetId))]
        
        # if proteins are used for panel assay, then gi, otherwise gene id
        if (isGiCur) {
          curDtFnl <- data.table(AID = aid, sidcidstack, OUTCOME = panelRes, GI = targetId, GENEID = "e", TAG = "p")
        } else {
          curDtFnl <- data.table(AID = aid, sidcidstack, OUTCOME = panelRes, GI = "e", GENEID = targetId, TAG = "p")
        }
        dtFnl <- rbind(curDtFnl, dtFnl)
      }
	  setnames(dtFnl, c("AID", "SID", "CID", "OUTCOME", "GI", "GENEID", "TAG"))
      
      return(dtFnl)
      #fwrite(dtFnl, file = savedFn, append = TRUE)
	  #return(invisible(NULL))
    }
    
    ####################################################
    dt <- dat[selRowIdx, selColNames, with = FALSE]
    ####################################################
      
    # important, originally it should be in the same columns as [outcome]
    # target_id can be either gi, geneid, or others (such as cell line names) 
	# if not the previous two 	
    targetId <- unname(unlist(dt[1, 3:ncol(dt)]))
      
    # 1st judge, some target_id was put in the first section of the corresponding target
    if (all(targetId == "")) {
      # aid = 434974
      targetIdValue <- dat[selRowIdx[1], targetIdColNameIdx, with = FALSE]
      targetId <- unname(unlist(targetIdValue))
    } 
      
    # 2nd judge
    if (all(targetId == "")) {
      # select target_id row index, and select panel_outcome column index - 1 becasue sometime, target_id is put into the ahead column of [outcome]
      # aid = 588327
      targetIdValue <- dat[selRowIdx[1], which(colNames %in% panelOutcomePatValue) - 1, with = FALSE]
      targetId <- unname(unlist(targetIdValue))
    }
      
    # 3rd: manual judge
    if (aid == "686975") {
      targetId <- c("38327039", "83318444")
    }
    if (aid == "624091") {
      targetId <- c("CDC2/CycB1", "CDK2/CycA2", "CDK2/CycE1", "CDK5/p25", "GSK3a", "GSK3b")
    }
    if (aid == "651573") {
      # outcome_col_names as the taget_id
      outcomeColNames <- grep("outcome", colNames, ignore.case = TRUE, value = TRUE)
      targetId <- setdiff(outcomeColNames, "PUBCHEM_ACTIVITY_OUTCOME")
    }
    if (aid == "651549") {
      targetId <- c("120649", "120649", "78486550", "78486550")
    }

    
    # pn: panel_name
    # in case no panel target id, but just panel name with number, we add the prefix, [pn]
    # which is used for distinguising the [real gene_id]
    # aid = 651838, where panel name is number, so add the "pn" on the front
    if (!any(grepl("target", panelAssay, ignore.case = TRUE))) {
      targetId <- paste("pn", targetId, sep = "-")
    }
    
    # length should be different
    if (length(targetId) != length(panelOutcomePatValue)) {
      stop("Length of targetId is different from length of panelOutcomePatValue \n")
    }
    
    
    # final result
    # just extract sid and cid column
    sidcid <- dt[2:.N, 1:2]
    # panel outcome, convert the outcome matrix into long vector
    # Active or Inactive, etc
    panelRes <- dt[2:.N, 3:ncol(dt)]
    panelRes <- c(t(panelRes))
    # length(panel_res)
    sidcidstack <- sidcid[rep(seq_len(nrow(sidcid)), each = length(targetId))]
      
    # if proteins are used for panel assay, then gi, otherwise gene id
    if (isGi) {
      dtFnl <- data.table(AID = aid, sidcidstack, OUTCOME = panelRes, GI = targetId, GENEID = "e", TAG = "p")
    } else {
      dtFnl <- data.table(AID = aid, sidcidstack, OUTCOME = panelRes, GI = "e", GENEID = targetId, TAG = "p")
    }
	setnames(dtFnl, c("AID", "SID", "CID", "OUTCOME", "GI", "GENEID", "TAG"))
	
	return(dtFnl)
	#fwrite(dtFnl, file = savedFn, append = TRUE)
	#return(invisible(NULL))
  }
  
  cat("No information was parsed at all \n")
  return(invisible(NULL))
}
##############################################


createFolders <- function(parentFolder) {
  # purpose
  #   make all folders parentFolder/data, parentFolder/mapping_file, and parentFolder/output
  # input
  #   parentFolder: string, like: /mnt/LHRI/Bioinformatics/Projects/DAVIDUpdate/pubchem
  # output
  #   created folders for pubchem ftp bioassay parsing project
  
  # example
  # parentFolder <- "/mnt/LHRI/Bioinformatics/Projects/DAVIDUpdate/pubchem4/"
  
  
  if (!dir.exists(parentFolder)) {
    dir.create(parentFolder)
  } else {
    stop("parentFolder has existed. please change folder name \n")
  }
  
  lastChar <- substr(parentFolder, nchar(parentFolder), nchar(parentFolder))
  if (lastChar == "/") {
    parentFolder <- substr(parentFolder, 1, nchar(parentFolder) - 1)
  }
  
  cat("creating [mapping_file] folder \n")
  dir.create(paste0(parentFolder, "/mapping_file"))
  cat("creating [data] folder \n")
  dir.create(paste0(parentFolder, "/data"))
  cat("creating [output] folder \n")
  dir.create(paste0(parentFolder, "/output"))
  cat("Done \n")
  
  return(invisible(NULL))
}
# nullRes <- createFolders(parentFolder = "/mnt/LHRI/Bioinformatics/Projects/DAVIDUpdate/pubchem4/")


dlPcMappingFileFromFtp <- function(parentFolder) {
  # purpose
  #   download useful mapping files from pubchem FTP (BioAssay and Compound)
  # input
  #   savedFolder: string, like savedFolderTmp <- "/mnt/LHRI/Bioinformatics/Projects/DAVIDUpdate/pubchem2"
  # output
  #   downloaded file in the $savedFolder/mapping_file directory
  
  
  lastChar <- substr(parentFolder, nchar(parentFolder), nchar(parentFolder))
  if (lastChar == "/") {
    parentFolder <- substr(parentFolder, 1, nchar(parentFolder) - 1)
  }
 
  mappingFileFolder <- paste0(parentFolder, "/mapping_file")
  if (!dir.exists(mappingFileFolder)) {
    stop("[mappingFileFolder] does not exist. please check if created in the 1st step \n")
  }
  
  # mapping files
  # $savedFolder/mapping_file
  aid2ActOutcomeMethod <- "ftp://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay/Extras/Aid2ActivityOutcomeMethod.gz"
  download.file(url = aid2ActOutcomeMethod, destfile = paste0(mappingFileFolder, "/Aid2ActivityOutcomeMethod.gz"), quiet = TRUE)
  cat("Aid2ActivityOutcomeMethod.gz was downloaded \n")
  
  aid2GiGeneidAccessionUniprot <- "ftp://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay/Extras/Aid2GiGeneidAccessionUniprot.gz"
  download.file(url = aid2GiGeneidAccessionUniprot, destfile = paste0(mappingFileFolder, "/Aid2GiGeneidAccessionUniprot.gz"), quiet = TRUE)
  cat("Aid2GiGeneidAccessionUniprot.gz was downloaded \n")
  
  cid2Synonym <- "ftp://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-Synonym-filtered.gz"
  download.file(url = cid2Synonym, destfile = paste0(mappingFileFolder, "/CID-Synonym-filtered.gz"), quiet = TRUE)
  cat("CID-Synonym-filtered.gz was downloaded \n")
  
  cid2Title <- "ftp://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-Title.gz"
  download.file(url = cid2Title, destfile = paste0(mappingFileFolder, "/CID-Title.gz"), quiet = TRUE)
  cat("CID-Title.gz was downloaded \n")
  
  cat("all mapping files were saved in:", mappingFileFolder, "\n")
  
  return(invisible((NULL)))
}
# nullRes <- dlPcMappingFileFromFtp(parentFolder = "/mnt/LHRI/Bioinformatics/Projects/DAVIDUpdate/pubchem4/")


extractGoodAid <- function(parentFolder) {
  # purpose
  #   extract good AID for parsing (conformatory and having GeneId or Gi
  # input
  #   savedDir: string, parent dir for mapping_file
  # output
  #   returned good aid for parsing

  require(data.table)
  
  lastChar <- substr(parentFolder, nchar(parentFolder), nchar(parentFolder))
  if (lastChar == "/") {
    parentFolder <- substr(parentFolder, 1, nchar(parentFolder) - 1)
  }
  
  mappingFileFolder <- paste0(parentFolder, "/mapping_file")
  if (!dir.exists(mappingFileFolder)) {
    stop("[mappingFileFolder] does not exist. please check if created in the 1st step \n")
  }
  
  aidGeneidGiMap <- fread(paste0(mappingFileFolder, "/Aid2GiGeneidAccessionUniprot.gz"))
  # both Gi and Geneid are not NA at all
  # Geneid is NA but Gi is not NA
  # Gi is NA but Geneid is not NA
  s1 <- aidGeneidGiMap[!is.na(Gi), ][!is.na(Geneid), ]
  s2 <- aidGeneidGiMap[is.na(Geneid), ][!is.na(Gi), ]
  s3 <- aidGeneidGiMap[is.na(Gi), ][!is.na(Geneid), ]
  aidGiMap <- rbind(s1, s2)
  aidGeneidMap <- copy(s3)
  aidVecGi <- aidGiMap[, unique(AID)]
  aidVecGeneid <- aidGeneidMap[, unique(AID)]
  
  aidOutputMethodMap <- fread(paste0(mappingFileFolder, "/Aid2ActivityOutcomeMethod.gz"))
  aidVecConf <- aidOutputMethodMap[Aid2ActivityOutcomeMethod == "Confirmatory", unique(AID)]
  
  aidVecGiConf <- intersect(aidVecGi, aidVecConf)  
  aidVecGeneidConf <- intersect(aidVecGeneid, aidVecConf)
  #length(c(aidVecGiConf, aidVecGeneidConf)) # 140447
  aidVecGiGeneidConf <- unique(c(aidVecGiConf, aidVecGeneidConf))
  #length(aidVecGiGeneidConf) # 140421
  # some bioassays have both Genei and Gi
  
  aidVecGiGeneidConfDt <- data.table(aidVecGiGeneidConf = aidVecGiGeneidConf)
  fwrite(aidVecGiGeneidConfDt, file = paste0(mappingFileFolder, "/aidVecGiGeneidConf.csv"))
  cat("[aidVecGiGeneidConfDt] was saved in:", mappingFileFolder, "\n")
  
  return(aidVecGiGeneidConfDt)
}
# aidDt <- extractGoodAid(parentFolder = "/mnt/LHRI/Bioinformatics/Projects/DAVIDUpdate/pubchem4/")


# dlPcAidCsvFromFtp <- function(parentFolder, aidDt, numDict = 1400, ncores = 30) {
dlPcAidCsvFromFtp <- function(parentFolder, aidDt, ncores = 30) {
  # purpose
  #   download bioassay csv file and unzip to *.csv.gz format
  # input
  #   parentFolder: string, parent path
  #   aidDt: data.table with one column for good bioassay ids
  #   ncores: integer, how many cores to be used, default = 30
  # output
  #   downloaded *.csv.gz files and put them into $savedFolder/data folder
  
  
  # update:
  # 2020-07-12
  #   numDict: this argument is removed, because this one can be automatically determined from pubchem bioassay ftp content
  
  require(data.table)
  require(parallel)
  require(RCurl)
  
  
  
  lastChar <- substr(parentFolder, nchar(parentFolder), nchar(parentFolder))
  if (lastChar == "/") {
    parentFolder <- substr(parentFolder, 1, nchar(parentFolder) - 1)
  }
  
  dataFolder <- paste0(parentFolder, "/data")
  if (!dir.exists(dataFolder)) {
    stop("[dataFolder] does not exist. please check if created in the 1st step \n")
  }
  
  mappingFileFolder <- paste0(parentFolder, "/mapping_file")
  if (!dir.exists(mappingFileFolder)) {
    stop("[mappingFileFolder] does not exist. please check if created in the 1st step \n")
  }
  
  # 2020-07-12
  # automatically determine the [numDict]
  ftpPcBioassay <- "ftp://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay/CSV/Data/"
  webContent <- getURL(ftpPcBioassay, verbose = FALSE, ftp.use.epsv = TRUE, dirlistonly = FALSE)
  # write to disk
  ftpWebContentFn <- paste0(mappingFileFolder, "/ftpWebContent.txt")
  cat(webContent, file = ftpWebContentFn)
  #
  tmp1 <- fread(ftpWebContentFn)
  tmp1 <- unlist(strsplit(tmp1[.N][[ncol(tmp1)]], split = "\\."))[1]
  tmp1 <- as.integer(unlist(strsplit(tmp1, split = "_"))[2])
  # integer division
  tmp1 <- tmp1 %/% 1000
  
  
  if (tmp1 > 1000 & tmp1 <= 1100) { 
    numDict <- 1100
  } else if (tmp1 > 1100 & tmp1 <= 1200) {
    numDict <- 1200
  } else if (tmp1 > 1200 & tmp1 <= 1300) {
    numDict <- 1300
  } else if (tmp1 > 1300 & tmp1 <= 1400) {
    numDict <- 1400
  } else if (tmp1 > 1400 & tmp1 <= 1500) {
    numDict <- 1500
  } else if (tmp1 > 1500 & tmp1 <= 1600) {
    numDict <- 1600
  } else if (tmp1 > 1600 & tmp1 <= 1700) {
    numDict <- 1700
  } else if (tmp1 > 1700 & tmp1 <= 1800) {
    numDict <- 1800
  } else if (tmp1 > 1800 & tmp1 <= 1900) {
    numDict <- 1900
  } else if (tmp1 > 1900 & tmp1 <= 2000) {
    numDict <- 2000
  } else if (tmp1 > 2000 & tmp1 <= 2100) {
    numDict <- 2100
  } else if (tmp1 > 2100 & tmp1 <= 2200) {
    numDict <- 2200
  } else if (tmp1 > 2200 & tmp1 <= 2300) {
    numDict <- 2300
  } else if (tmp1 > 2300 & tmp1 <= 2400) {
    numDict <- 2400
  } else if (tmp1 > 2400 & tmp1 <= 2500) {
    numDict <- 2500
  } else if (tmp1 > 2500 & tmp1 <= 2600) {
    numDict <- 2600
  } else if (tmp1 > 2600 & tmp1 <= 2700) {
    numDict <- 2700
  } else if (tmp1 > 2700 & tmp1 <= 2800) {
    numDict <- 2800
  } else if (tmp1 > 2800 & tmp1 <= 2900) {
    numDict <- 2900
  } else if (tmp1 > 2900 & tmp1 <= 3000) {
    numDict <- 3000
  } else {
    stop("[numDict] doese not set correctly \n")
  }
  
  
  # mimic ftp file names
  leftNum <- seq(1, numDict * 1000, by = 1000)
  rightNum <- (1:numDict) * 1000
  leftCharVec <- formatC(leftNum, width = 7, format = "d", flag = "0")
  rightCharVec <- formatC(rightNum, width = 7, format = "d", flag = "0")
  aidList <- paste0(leftCharVec, "_", rightCharVec, ".zip")
  #length(aidList) # 1400
  aidListDt <- data.table(aidList, leftNum = leftNum)
  
  # all good bioassay ids
  aidSorted <- aidDt[, sort(unique(aidVecGiGeneidConf))]
  # check the maximum AID number
  # aidDt[order(aidVecGiGeneidConf)]
  
  getInterval <- function(curAid) {
    #curAid <- aidSorted[1]
	aidRange <- data.table(aid = curAid, fileRange = aidListDt[max(which(leftNum <= curAid)), aidList])
	return(aidRange)
  }
  
  ##############################################################################
  aidRangeList <- parallel::mclapply(aidSorted, getInterval, mc.cores = ncores)
  # length(aidRangeList)
  ##############################################################################

  aidRangeDt <- rbindlist(aidRangeList)
  # aidVecPerRange <- aidRangeDt[, aid, by = fileRange]
  # aidVecPerRange[1]
  
  dlFileNames <- aidRangeDt[, unique(fileRange)] # 1012
  
  rm(aidRangeList)
  cat("aidRangeDt was finished \n")
  flush.console()
  
  cat("downloading and unzipping files. taking time this step... \n")  
  dlUnzipBioassay <- function(fn) {
    # purpose
	#   download pubhcem bioassay csv file with zip format
	# input
	#   fn: string, file name with full URL
	
	aidVecPerRange <- aidRangeDt[fileRange == fn, aid]
	aidVecPerRangeFn <- paste0(aidVecPerRange, ".csv.gz")
	# full path
	# wrong way
	#aidVecPerRangeFnFp <- paste0(dataFolder, "/", gsub(".zip", "", fn), "/", aidVecPerRangeFn)
	# good way
	aidVecPerRangeFnFp <- paste0(gsub(".zip", "", fn), "/", aidVecPerRangeFn)
	
	
	# download
	fnUrl <- paste0("ftp://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay/CSV/Data/", fn)
	destfile <- paste0(dataFolder, "/", fn)
	# set method = "wget"
	download.file(url = fnUrl, destfile = destfile, quiet = TRUE, method = "wget")
	
	# unzip specified file
	unzip(zipfile = destfile, files = aidVecPerRangeFnFp, exdir = dataFolder, junkpaths = TRUE)
	
    # delete zip file
	unlink(destfile)
	###################################
	
	return(invisible(NULL))
  }
  #################################################################################
  #dlBioassay(fn)
  nullRes <- parallel::mclapply(dlFileNames, dlUnzipBioassay, mc.cores = ncores)
  #################################################################################
  flush.console()
  cat("downloading and unziping were finished... \n")
  cat("pubchem bioassay [*.csv.gz] files were saved in:", dataFolder, "\n")  
  flush.console()
    
  return(invisible(NULL))
}
# system.time(nullRes <- dlPcAidCsvFromFtp(parentFolder = "/mnt/LHRI/Bioinformatics/Projects/DAVIDUpdate/pubchem4/", aidDt, numDict = 1400, ncores = 30))
# 5.5 min
# length(list.files("/mnt/LHRI/Bioinformatics/Projects/DAVIDUpdate/pubchem4/data"))
# 140421
# 2020-07-12 blow
# system.time(nullRes <- dlPcAidCsvFromFtp(parentFolder = "/data/Members/haom/NIAID_Projects/DAVID_MingHAO/pubchemNew20200712", aidDt, ncores = 30))


parseBioassay <- function(parentFolder, aidDt, ncores = 32) {
  # purpose
  # input
  #   fnVec: character vector, file names, such as c("280841.csv.gz", "738352.csv.gz")
  #   fnDir: string, directory for files
  #   ncores: integer, how many cores are used
  # output
  
  
  require(data.table)
  require(parallel)
  
  lastChar <- substr(parentFolder, nchar(parentFolder), nchar(parentFolder))
  if (lastChar == "/") {
    parentFolder <- substr(parentFolder, 1, nchar(parentFolder) - 1)
  }
  
  dataFolder <- paste0(parentFolder, "/data")
  if (!dir.exists(dataFolder)) {
    stop("[dataFolder] does not exist. please check if created in the 1st step \n")
  }
  
  outputFolder <- paste0(parentFolder, "/output")
  if (!dir.exists(outputFolder)) {
    stop("[outputFolder] does not exists. please check if created in the 1st step \n")
  }
  
  
  # bioassay file name full path
  aidFnFp <- paste0(dataFolder, "/", aidDt[, aidVecGiGeneidConf], ".csv.gz")
  # print(aidFnFp[1])
  cat("starting to parse all bioassay csv.gz files.please wait... \n")
  ################################################################
  resList <- mclapply(aidFnFp, getPubChemInfo, mc.cores = ncores)
  ################################################################
  dtFnl <- data.table::rbindlist(resList)
  
  fwrite(dtFnl, file = paste0(outputFolder, "/dtFnl.csv"))
  flush.console()
  cat("[dtFnl.csv] was saved in:",  outputFolder, "\n")
  flush.console()
  
  return(dtFnl)
}
# system.time(dtFnl <- parseBioassay(parentFolder = "/mnt/LHRI/Bioinformatics/Projects/DAVIDUpdate/pubchem4/", aidDt, ncores = 32))
# 5.5 min


analyzeParsedData <- function(parentFolder, dtFnl) {
  # purpose
  #   analyze the parsed bioassay data. one need geneid vs cid (all active)
  # input
  #   parentFolder:, string
  #   dtFnl: data.table, parsed data from last step
  # output
  #   tab-delimited format data such as: GeneId   CID~CompoundName   AID1,AID2,AID3
  
  # example
  # parentFolder <- "/mnt/LHRI/Bioinformatics/Projects/DAVIDUpdate/pubchem5/"
  # outputFolder <- paste0(parentFolder, "/output")
  # dtFnl <- data.table::fread(paste0(outputFolder, "/dtFnl.csv"))
  
  lastChar <- substr(parentFolder, nchar(parentFolder), nchar(parentFolder))
  if (lastChar == "/") {
    parentFolder <- substr(parentFolder, 1, nchar(parentFolder) - 1)
  }
  
  dtAct <- dtFnl[OUTCOME %in% c("Active", "Probe"), ]
  
  rm(dtFnl)
  invisible(gc())
  
  dtAct$OUTCOME[dtAct$OUTCOME == "Probe"] <- "Active"
  
  dtAct <- dtAct[SID != "e", ]
  dtAct <- dtAct[SID != "", ]
  dtAct <- dtAct[CID != "e", ]
  dtAct <- dtAct[CID != "", ]
  
  if (!all(dtAct[, unique(TAG)] %in% c("s", "p"))) {
    stop("check the result, which may include [m] case \n")
  }
  
  # required for later integer join
  # dtAct$AID <- as.integer(dtAct$AID)
  # dtAct$SID <- as.integer(dtAct$SID)
  # dtAct$CID <- as.integer(dtAct$CID)
  
  for (jcol in c("AID", "SID", "CID")) {
    set(dtAct, j = jcol, value = as.integer(dtAct[[jcol]]))
  }
  
  # panel assay result
  dtActConfP <- dtAct[TAG == "p", ]
  # here, some Gis may need map to Geneids
  # 2020-07-04, already mapped to Geneid from Gi below
  dtActConfP <- dtActConfP[GI != "", ]
  dtActConfP <- unique(dtActConfP)
  
  
  # the most common case, s, which means no target information in the bioassay csv file
  dtActConfS <- dtAct[TAG == "s", ]
  dtActConfS <- unique(dtActConfS)
  
  # used for mapping to geneid
  aidS <- dtActConfS[, unique(AID)]
  
  aidGiGeneidMapFn <- paste0(parentFolder, "/mapping_file/Aid2GiGeneidAccessionUniprot.gz")
  aidGiGeneidMapDt <- fread(aidGiGeneidMapFn)
  
  ###############
  # 2020-07-04
  # convert Gi to Geneid for panel cases
  giGeneidConvMap <- unique(aidGiGeneidMapDt[, list(Gi, Geneid)][!is.na(Gi)][!is.na(Geneid)])
  # dtActConfP[, GI][1:3]
  # char
  for (jcol in c("Gi", "Geneid")) {
    set(giGeneidConvMap, j = jcol, value = as.character(giGeneidConvMap[[jcol]]))
  }
  setnames(giGeneidConvMap, c("GI", "Geneid"))
  # giGeneidConvMap[, GI][1:3]
  # char
  # right join
  tmpRes <- giGeneidConvMap[dtActConfP, on = "GI"]
  ###
  dtActConfPGeneid <- tmpRes[!is.na(Geneid)]
  dtActConfPGi <- tmpRes[is.na(Geneid)]
  # tmpRes[, Geneid][1:3]
  ##################
  
  # only those bioassays for 's' case, which exclude 'p' case
  aidGiGeneidMapDtS <- aidGiGeneidMapDt[AID %in% aidS, ]
  
  # geneid is the main interesting one
  aidGeneidMap <- aidGiGeneidMapDtS[!is.na(Geneid), ]
  aidGeneidMap <- unique(aidGeneidMap[, list(AID, Geneid)])
  
  # some bioassays only have Gi (no Geneid, by check protein DB, these Gis may not
  # have corresponding Geneids, but need double check later
  aidGiMap <- aidGiGeneidMapDtS[is.na(Geneid), ]
  aidGiMap <- unique(aidGiMap[, list(AID, Gi)])
  
  ####################
  # inner join, geneid
  aidCidGeneidPrep <- dtActConfS[aidGeneidMap, on = "AID", nomatch = 0]
  ############
  # 2020-07-04
  aidCidGeneidPrep <- aidCidGeneidPrep[, list(AID, SID, CID, Geneid)]
  # rbind Geneid in panal case: dtActConfPGeneid
  aidCidGeneidPrep <- rbind(aidCidGeneidPrep, dtActConfPGeneid[, list(AID, SID, CID, Geneid)])
  ############
  # unique(aidCidGeneidPrep)
  cidGeneid <- unique(aidCidGeneidPrep[, list(CID, Geneid)])
  aidCid <- unique(aidCidGeneidPrep[, list(AID, CID)]) 
  cid2aidList <- aidCid[, paste(unique(AID), collapse = ","), by = CID]
  setnames(cid2aidList, c("CID", "aidList"))
  # right join
  cidGeneidAidlist <- cid2aidList[cidGeneid, on = "CID"]
  # cidGeneidAidlist[is.na(Geneid)]
  # cidGeneidAidlist[30:35]
  
  ################
  # inner join, gi
  aidCidGiPrep <- dtActConfS[aidGiMap, on = "AID", nomatch = 0]
  aidCidGiPrep <- aidCidGiPrep[, list(AID, SID, CID, Gi)]
  setnames(aidCidGiPrep, c("AID", "SID", "CID", "GI"))
  ##############
  # 2020-07-04
  # rbind panal bioassay: dtActConfPGi
  aidCidGiPrep <- rbind(aidCidGiPrep, dtActConfPGi[, list(AID, SID, CID, GI)])
  ##############
  # unique(aidCidGiPrep)
  cidGi <- unique(aidCidGiPrep[, list(CID, GI)])
  aidCid <- unique(aidCidGiPrep[, list(AID, CID)])
  cid2aidList <- aidCid[, paste(unique(AID), collapse = ","), by = CID]
  setnames(cid2aidList, c("CID", "aidList"))
  cidGiAidList <- cid2aidList[cidGi, on = "CID"]
  # cidGiAidList[is.na(GI)]
  # cidGiAidList[1:2]
  
  # used for reducing CID-Mapping file size
  allCid <- unique(c(cidGeneidAidlist[, unique(CID)], cidGiAidList[, unique(CID)]))
  
  cat("mapping compound names. please wait... \n")
  cid2nameMapFn <- paste0(parentFolder, "/mapping_file/CID-Synonym-filtered.gz")
  
  #################################
  # taking time
  compdNamePrep <- fread(cid2nameMapFn)
  setnames(compdNamePrep, c("CID", "compdName"))
  # compdName[, list(unique(CID))]
  compdName <- compdNamePrep[CID %in% allCid, ]
  rm(compdNamePrep)
  invisible(gc())
  #################################
  # cid number, small to large
  compdName <- compdName[order(CID), ]
  compdName[, ncharName := nchar(compdName)]
  compdName[, idx := 1:.N]
  cidShortNameIdx <- compdName[, idx[which.min(ncharName)], by = CID]
  cidShortNamePrep <- compdName[cidShortNameIdx$V1, ]
  rm(compdName)
  invisible(gc())
  cidShortName <- cidShortNamePrep[, list(CID, compdName)]
  
  rm(cidShortNamePrep)
  invisible(gc())
  rm(cidShortNameIdx)
  invisible(gc())
  
  # right join
  cidNameAidlistGeneid <- cidShortName[cidGeneidAidlist, on = "CID"]
  cidNameAidlistGi <- cidShortName[cidGiAidList, on = "CID"]
  
  rm(cidShortName)
  invisible(gc())
  
  cidGeneidFnl <- cidNameAidlistGeneid[, list(Geneid, paste(CID, compdName, sep = "~"), aidList)]
  setnames(cidGeneidFnl, c("Geneid", "CidCompdName", "AidList"))
  
  cidGiFnl <- cidNameAidlistGi[, list(GI, paste(CID, compdName, sep = "~"), aidList)]
  setnames(cidGiFnl, c("Gi", "CidCompdName", "AidList"))
  
  data.table::fwrite(cidGeneidFnl, file = paste0(parentFolder, "/output/cidGeneidFnl.tsv"), sep = "\t")
  data.table::fwrite(cidGiFnl, file = paste0(parentFolder, "/output/cidGiFnl.tsv"), sep = "\t")
  
  if (TRUE) {
    # cid vs. title
	flush.console()
	cat("processing CID vs. title. please wait... \n")
	##########################################
	cid2titleMapFn <- paste0(parentFolder, "/mapping_file/CID-Title.gz")
    cidTitleMap <- fread(cid2titleMapFn)
    #cidTitleMap[, list(unique(V1))]
    setnames(cidTitleMap, c("CID", "Title"))
	##########################################
	
	cidTitleAidlistGeneid <- cidTitleMap[cidGeneidAidlist, on = "CID"]
    # cidTitleAidlistGeneid[is.na(Geneid), ]
	cidTitleAidlistGi <- cidTitleMap[cidGiAidList, on = "CID"]
	
	rm(cidTitleMap)
	invisible(gc())
	
	cidTitleGeneidFnl <- cidTitleAidlistGeneid[, list(Geneid, paste(CID, Title, sep = "~"), aidList)]
    setnames(cidTitleGeneidFnl, c("Geneid", "CidTitle", "AidList"))
  
    cidTitleGiFnl <- cidTitleAidlistGi[, list(GI, paste(CID, Title, sep = "~"), aidList)]
    setnames(cidTitleGiFnl, c("Gi", "CidTitle", "AidList"))
	
	data.table::fwrite(cidTitleGeneidFnl, file = paste0(parentFolder, "/output/cidTitleGeneidFnl.tsv"), sep = "\t")
    data.table::fwrite(cidTitleGiFnl, file = paste0(parentFolder, "/output/cidTitleGiFnl.tsv"), sep = "\t")
  }
  
  flush.console()
  cat("final results were saved in:", paste0(parentFolder, "/output"), "\n")
  
  return(invisible(NULL))
}
# system.time(nullRes <- analyzeParsedData(parentFolder, dtFnl)) # 11 min


runPipeline <- function(parentFolder = "full_path_folder", ncores = 30) {
  # purpose
  #   run a pipeline for getting pubchem bioassay data
  # input
  #   parentFolder: string, full path folder
  #   ncores: integer, how many cores to be used

  st <- Sys.time()
  
  nullRes <- createFolders(parentFolder)
  nullRes <- dlPcMappingFileFromFtp(parentFolder)
  aidDt <- extractGoodAid(parentFolder)
  nullRes <- dlPcAidCsvFromFtp(parentFolder, aidDt, ncores = ncores)
  dtFnl <- parseBioassay(parentFolder, aidDt, ncores = ncores)
  nullRes <- analyzeParsedData(parentFolder, dtFnl)
  
  et <- Sys.time()
  tcost <- et - st
  tnum <- as.numeric(tcost)
  tunit <- attributes(tcost)$units
  tcostStr <- paste0("The whole pipeline was done: Time Used: ", tnum, " ", tunit)
  print(tcostStr)

  return(invisible(NULL))
}



