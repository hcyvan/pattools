library(dplyr)


configDefaultPath <- './config.default.yaml'
configPath <- './config.yaml'
if (!file.exists(configPath)) {
  if(file.exists(configDefaultPath)) {
    configPath <- configDefaultPath
  } else {
    stop('please add configure file [config.yaml] to this project root directory')
  }
}

CONFIG <- yaml::read_yaml(configPath)
CONFIG['DataRaw']<-file.path(CONFIG$DataDir, 'raw')
CONFIG['DataInter']<-file.path(CONFIG$DataDir, 'intermediate')
CONFIG['DataResult']<-file.path(CONFIG$DataDir, 'result')

############################################################################################
getExt <- function(file_path) {
  file_name <- basename(file_path)
  parts <- strsplit(file_name, "\\.")[[1]]
  if (length(parts) > 1) {
    tolower(tail(parts, 1))
  } else {
    NULL
  }
}

loadData <- function(name, ext=NULL, header=FALSE, force.refresh=FALSE){
  if (is.null(ext)) {
    ext<-getExt(name)
  }
  if (is.null(ext)) {
    stop(sprintf("Unknown file type: %s", name))
  }
  ext.path <- name
  rds.path <- paste(ext.path,'rds',sep = '.')
  if(file.exists(rds.path) && !force.refresh){
    readRDS(rds.path)
  }else if(file.exists(ext.path)) {
    if(ext=='csv'){
      data<-read.csv(ext.path)
      saveRDS(data, rds.path)
      data
    }else if(ext=='tsv'){
      data<-read.csv(ext.path,sep = '\t',check.names = FALSE)
      saveRDS(data, rds.path)
      data
    }else if(ext=='bed'){
      data<-read.csv(ext.path,sep = '\t',check.names = FALSE, header=header)
      saveRDS(data, rds.path)
      data
    } else {
      stop(sprintf("The file format <%s> is not supported", ext)); 
    }
  } else{
    stop(sprintf("File not exist: %s and %s", rds.path, ext.path))
  }
}

saveImage <- function(file,...){
  file.path=file.path(CONFIG$DataResult, file)
  if (endsWith(file, '.pdf')){
    pdf(file=file.path, ...)
  }
}
#########################################################################
Clincal <- setRefClass(
  "Clincal",
  fields = list(data = "data.frame"),
  methods = list(
    initialize = function() {
      data<<-loadData(file.path(CONFIG$DataRaw,"CellInfo.csv"),force.refresh=TRUE)
      data<<-data[data$SubjectType!="",]
      data<<-data[is.na(data$QCFailed2),]
      
    },
    pickSubjectIdByDataID=function(df){
      df<-df[,match(data$DataID, colnames(df))]
      colnames(df)<-data$SubjectID
      df
    },
    pickSubjectTypeByDataID=function(df){
      df<-df[,match(data$DataID, colnames(df))]
      colnames(df)<-data$SubjectType
      df
    },
    pickCellTypeByDataID=function(df){
      df<-df[,match(data$DataID, colnames(df))]
      colnames(df)<-data$CellType
      df
    }
  )
)

