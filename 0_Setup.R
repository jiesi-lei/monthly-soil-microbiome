library("vegan")
library("ieggr")

#clustering
library(ade4)
library(vegan)
library(gclus)
library(cluster)
library(RColorBrewer)
library(labdsv)
library(ape)
library(MASS)
library(ellipse)

#ordination
library(ape)
library(ggplot2)
library(plyr)

#dissimilarity tests
library(agricolae)
library(BiodiversityR)
library(ieggr)
library(dplyr)

#response ratio
library(plotrix)
library(openxlsx)
library(showtext)
font_add("Arial", "/Library/Fonts/Arial.ttf")  # Use the actual file path
showtext_auto()

#network analysis
library(dplyr)
library(MENA)#not publicly accessible for now, users can turn to http://ieg4.rccc.ou.edu/MENA/ instead.

#Env
library(ggsignif)

#Mantel test
library(Imap)

#CCA
library(ggvegan)

#Env variables
env <- read.csv("monthly_ENV.csv",row.names = 1)


#functions
sig_label <- function(df,col){
  df[df[,col]<0.05,which(colnames(df) == 'sig')] = '*'
  df[df[,col]<0.01,which(colnames(df) == 'sig')] = '**'
  df[df[,col]<0.001,which(colnames(df) == 'sig')] = '***'
  return(df)
}

log_scale <- function(x){
  if (x > 0){
    x <- log10(x)
  } else if (x == 0){
    x <- 0
  } else {
    x <- -log10(abs(x))}}


data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE),
      se = plotrix::std.error(x[[col]], na.rm=TRUE),
      n = length(x[[col]]))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

data_summary_varname_as_mean <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE),
      se = plotrix::std.error(x[[col]], na.rm=TRUE),
      n = length(x[[col]]))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  return(data_sum)
}

del_ghost <- function(data, n_start, n_end)
{
  data[is.na(data)] = 0
  data[is.null(data)] = 0
  
  ind_del = which(colSums(data[,n_start:n_end])==0)
  
  if(length(ind_del) == 0)
  {
    return(data)
  }
  else
  {
    return(data[,-ind_del])
  }
}

# theme
theme1 <- theme_linedraw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position="right")


