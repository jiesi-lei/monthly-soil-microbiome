# within-group dispersion
library(ieggr)
library(dplyr)
library(ggpubr)
library(plotrix)
Group<-read.csv("monthly_SampleList.csv",row.names = 1,stringsAsFactors = TRUE, check.names = FALSE) 

# 1 calculate distance matrix ----
library(ieggr)
beta.binary=beta.g(log(dat+1),dist.method = c("jaccard","bray"),abundance.weighted = FALSE)
beta.weight=beta.g(log(dat+1),dist.method = c("jaccard","bray"),abundance.weighted = TRUE)

betaTD=list(jaccard=beta.binary$jaccard.uw,ruzicka=beta.weight$ruzicka.wt,sorensen=beta.binary$sorensen.uw,bray=beta.weight$bray.wt)
save(betaTD,file="datasets/dissimilarity/16S_dist.rdata")

beta.binary=beta.g(log(gene+1),dist.method = c("jaccard","bray"),abundance.weighted = FALSE)
beta.weight=beta.g(log(gene+1),dist.method = c("jaccard","bray"),abundance.weighted = TRUE)

betaTD=list(jaccard=beta.binary$jaccard.uw,ruzicka=beta.weight$ruzicka.wt,sorensen=beta.binary$sorensen.uw,bray=beta.weight$bray.wt)
save(betaTD,file="datasets/dissimilarity/gene_dist.rdata")

# 16S ----
load(file="datasets/dissimilarity/16S_dist.rdata")
dist <- lapply(betaTD, function(x){
  dist.each <- as.matrix(x)
  id <- colnames(dist.each) %in% rownames(Group)
  dist.each <- dist.each[id,id]
  sum(colnames(dist.each) == rownames(Group))
  dist.each
})

dist.name = c("jaccard","ruzicka","sorensen","bray")
i = 4
dist = dist.3col(betaTD[[i]])
names(dist)[3] = dist.name[i]
dist %>% filter(name1 %in% rownames(Group) & name2 %in% rownames(Group)) -> dist
Group[match(dist$name1,rownames(Group)),] %>% `colnames<-`(paste0(names(Group),".1")) -> trt_name1
Group[match(dist$name2,rownames(Group)),] %>% `colnames<-`(paste0(names(Group),".2")) -> trt_name2

data.frame(dist,trt_name1,trt_name2) %>% 
  filter(Month.1 == Month.2) %>% 
  filter(Warming.1 == Warming.2) %>%
  mutate(Month.1 = as.factor(Month.1))-> 
  dist.same.month

ggline(dist.same.month, x = "Month.1", y = dist.name[i], color = "Warming.1",add = c("mean_se"))+ stat_compare_means(aes(group=Warming.1), label = "p.signif",method="wilcox.test", hide.ns=T, paired=F,vjust = 2)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

data.frame(dist,trt_name1,trt_name2) %>% 
  filter(Month.1 == Month.2) %>% 
  filter(Warming.1 != Warming.2) %>%
  filter(Block.1 == Block.2) %>%
  mutate(Month.1 = as.factor(Month.1))-> 
  dist.same.month
boxplot(dist.same.month$bray)
ggline(dist.same.month, x = "Month.1", y = dist.name[i],add = c("mean_se")) +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

xx = lm(dist.same.month$bray~poly(dist.same.month$Month.1,2))
summary(xx)

# GeoChip ----
load(file="datasets/dissimilarity/gene_dist.rdata")
dist <- lapply(betaTD, function(x){
  dist.each <- as.matrix(x)
  id <- colnames(dist.each) %in% rownames(Group)
  dist.each <- dist.each[id,id]
  sum(colnames(dist.each) == rownames(Group))
  dist.each
})

dist.name = c("jaccard","ruzicka","sorensen","bray")
i = 4
dist = dist.3col(betaTD[[i]])
names(dist)[3] = dist.name[i]
dist %>% filter(name1 %in% rownames(Group) & name2 %in% rownames(Group)) -> dist
Group[match(dist$name1,rownames(Group)),] %>% `colnames<-`(paste0(names(Group),".1")) -> trt_name1
Group[match(dist$name2,rownames(Group)),] %>% `colnames<-`(paste0(names(Group),".2")) -> trt_name2

data.frame(dist,trt_name1,trt_name2) %>% 
  filter(Month.1 == Month.2) %>% 
  filter(Warming.1 == Warming.2) %>%
  mutate(Month.1 = as.factor(Month.1))-> 
  dist.same.month
dist.same.month = filter(dist.same.month, name1 != "M04_C4")

ggline(dist.same.month, x = "Month.1", y = dist.name[i], color = "Warming.1",add = c("mean_se"))+ stat_compare_means(aes(group=Warming.1), label = "p.signif",method="wilcox.test", hide.ns=T, paired=F,vjust = 2)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

data.frame(dist,trt_name1,trt_name2) %>% 
  filter(Month.1 == Month.2) %>% 
  filter(Warming.1 != Warming.2) %>%
  filter(Block.1 == Block.2) %>%
  mutate(Month.1 = as.factor(Month.1))-> 
  dist.same.month
boxplot(dist.same.month$bray)
dist.same.month = filter(dist.same.month, name1 != "M04_H1")
dist.same.month = filter(dist.same.month, name1 != "M04_H3")
dist.same.month = filter(dist.same.month, name1 != "M04_H2")
dist.same.month = filter(dist.same.month, name1 != "M07_H4")
dist.same.month = filter(dist.same.month, name1 != "M12_H4")
#dist.same.month = filter(dist.same.month, name1 != "M04_H4")
ggline(dist.same.month, x = "Month.1", y = dist.name[i],add = c("mean_se")) +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ag <- aggregate(bray ~ Month.1, data = dist.same.month, mean)
ag$Month <- 1:12
xx = lm(ag$bray~poly(ag$Month,2))
summary(xx)

# plot ----
# 16S
load(file="datasets/dissimilarity/16S_dist.rdata")
dist <- lapply(betaTD, function(x){
  dist.each <- as.matrix(x)
  id <- colnames(dist.each) %in% rownames(Group)
  dist.each <- dist.each[id,id]
  sum(colnames(dist.each) == rownames(Group))
  dist.each
})

dist.name = c("jaccard","ruzicka","sorensen","bray")
i = 4
dist = dist.3col(betaTD[[i]])
names(dist)[3] = dist.name[i]
dist %>% filter(name1 %in% rownames(Group) & name2 %in% rownames(Group)) -> dist
Group[match(dist$name1,rownames(Group)),] %>% `colnames<-`(paste0(names(Group),".1")) -> trt_name1
Group[match(dist$name2,rownames(Group)),] %>% `colnames<-`(paste0(names(Group),".2")) -> trt_name2

data.frame(dist,trt_name1,trt_name2) %>% 
  filter(Month.1 == Month.2) %>% 
  filter(Warming.1 != Warming.2) %>%
  filter(Block.1 == Block.2) %>%
  mutate(Month.1 = as.numeric(Month.1))-> 
  dist.same.month
boxplot(dist.same.month$bray)
ggline(dist.same.month, x = "Month.1", y = dist.name[i],add = c("mean_se"),linetype = NULL) +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+theme_test()

ag <- aggregate(bray ~ Month.1, data = dist.same.month, mean)
ag$Month <- as.numeric(1:12)
ag$se <- aggregate(bray ~ Month.1, data = dist.same.month, std.error)$bray
xx = lm(ag$bray ~ ag$Month + I((ag$Month)^2))
summary(xx)
# R-squared:  0.3937;  p-value: 0.1053

ag$Month <- as.factor(ag$Month)

ggplot(ag, aes_string(x = "Month", y = "bray"))+
  geom_errorbar(aes(ymin = bray - se, ymax = bray + se ), width = 0.2, size = 0.3)+
  geom_point(shape = 21, size = 2.5, fill = "black", stroke = 0.2)+
  #scale_fill_manual(values = c("royalblue4","firebrick4"))+
  #scale_color_manual(values = c("royalblue4","firebrick4"))+
  #scale_x_continuous(limits = c(0.75,12.25), breaks = c(1:12), expand = c(0,0))+
  #scale_y_continuous(expand = c(0,0), limits = c(0,21)) +  
  expand_limits(x = 0) +
  ylab("Distance between W and C") +
  theme_test()+
  theme(axis.text = element_text(color = "black", size = 10),
        axis.title.y = element_text(color = "black", size = 10),
        legend.position = "none",
        axis.title.x = element_text(color = "black", size = 10),
        panel.border = element_rect(size = 1)) +
  stat_function(fun=function(x) 0.0003772*x + 0.0001577*x^2 + 0.5839264, size = 0.15, linetype = "dashed")

ggsave("figures/dissimilarity/16S_quadratic.pdf", width = 5, height = 3)

# GeoChip
load(file="datasets/dissimilarity/gene_dist.rdata")
dist <- lapply(betaTD, function(x){
  dist.each <- as.matrix(x)
  id <- colnames(dist.each) %in% rownames(Group)
  dist.each <- dist.each[id,id]
  sum(colnames(dist.each) == rownames(Group))
  dist.each
})

dist.name = c("jaccard","ruzicka","sorensen","bray")
i = 4
dist = dist.3col(betaTD[[i]])
names(dist)[3] = dist.name[i]
dist %>% filter(name1 %in% rownames(Group) & name2 %in% rownames(Group)) -> dist
Group[match(dist$name1,rownames(Group)),] %>% `colnames<-`(paste0(names(Group),".1")) -> trt_name1
Group[match(dist$name2,rownames(Group)),] %>% `colnames<-`(paste0(names(Group),".2")) -> trt_name2

data.frame(dist,trt_name1,trt_name2) %>% 
  filter(Month.1 == Month.2) %>% 
  filter(Warming.1 != Warming.2) %>%
  filter(Block.1 == Block.2) %>%
  mutate(Month.1 = as.numeric(Month.1))-> 
  dist.same.month
boxplot(dist.same.month$bray)
dist.same.month = filter(dist.same.month, name1 != "M04_H1")
dist.same.month = filter(dist.same.month, name1 != "M04_H3")
dist.same.month = filter(dist.same.month, name1 != "M04_H2")
dist.same.month = filter(dist.same.month, name1 != "M07_H4")
dist.same.month = filter(dist.same.month, name1 != "M12_H4")
ggline(dist.same.month, x = "Month.1", y = dist.name[i],add = c("mean_se"),linetype = NULL) +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+theme_test()

ag <- aggregate(bray ~ Month.1, data = dist.same.month, mean)
ag$Month <- as.numeric(1:12)
ag$se <- aggregate(bray ~ Month.1, data = dist.same.month, std.error)$bray
xx = lm(ag$bray ~ ag$Month + I((ag$Month)^2))
summary(xx)
# R2 = 0.64550;  p-value: 0.009398
ag$Month <- as.factor(ag$Month)

ggplot(ag, aes_string(x = "Month", y = "bray"))+
  geom_errorbar(aes(ymin = bray - se, ymax = bray + se ), width = 0.2, size = 0.3)+
  geom_point(shape = 21, size = 2.5, fill = "black", stroke = 0.2)+
  #scale_fill_manual(values = c("royalblue4","firebrick4"))+
  #scale_color_manual(values = c("royalblue4","firebrick4"))+
  #scale_x_continuous(limits = c(0.75,12.25), breaks = c(1:12), expand = c(0,0))+
  #scale_y_continuous(expand = c(0,0), limits = c(0,21)) +  
  expand_limits(x = 0) +
  ylab("Distance between W and C") +
  theme_test()+
  theme(axis.text = element_text(color = "black", size = 13),
        axis.title.y = element_text(color = "black", size = 13),
        legend.position = "none",
        axis.title.x = element_text(color = "black", size = 13),
        panel.border = element_rect(size = 1)) +
  stat_function(fun=function(x) -0.0069769*x +  0.0007512*x^2 + 0.0788570, size = 0.15)

ggsave("figures/dissimilarity/gene_quadratic.pdf", width = 5, height = 3)
