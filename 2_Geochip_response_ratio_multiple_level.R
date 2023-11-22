source("0_setup.R")
source("0_setup_gene.R")

dim(gene)

table(info$geneCategory)
table(info$subcategory1)
table(info$subcategory2)

##Select category of interest
#ategory<-c("Carbon Cycling","Nitrogen","Phosphorus","Stress","Electron transfer","Metal Homeostasis","Sulphur","Organic Remediation","Other","secondary metabolism","Secondary metabolism","virulence","Virulence","virus","Virus")
cat.select=c("Carbon cycling","Nitrogen","Phosphorus","Sulphur")
subcat.select=c("Starch","Pectin","Hemicellulose","Cellulose","Chitin","Lignin","Cutin","Cyanide","Glyoxylate cycle","Phospholipids","Tannins","Terpenes","Vanillin/Lignin")
c.deg.select = c("amyA","cda","glucoamylase","nplT","pula",#Starch
                 "ara","mannanase","xyla","xylanase",#Hemicellulose
                 "cellobiase","endoglucanase","exoglucanase",#Cellulous
                 "chitinase",#Chitin
                 "pectinase (pectate_lyase)",#Pectin
                  "glx","mnp","phenol_oxidase",#lignin
                 "cutinase",#cutin
                 "cdh")#Terpene

###By three seasons ================================================================
control = list(c(1:4,9:12),c(1:4,9:12,17:20,25:28,33:36,41:44,49:52),c(1:4,9:12,17:20))
warming = list(c(5:8,13:16),c(5:8,13:16,21:24,29:32,37:40,45:48,53:56),c(5:8,13:16,21:24))

rownames(gene1)[control[[1]]]
rownames(gene1)[warming[[1]]]
rownames(gene2)[control[[2]]]
rownames(gene2)[warming[[2]]]
rownames(gene3)[control[[3]]]
rownames(gene3)[warming[[3]]]

#c.deg.genes
rr = list()
rr.probe = list()
for (i in 1:3){
  print(c("processing cold_season1 samples -----",
          "processing warm_season samples -----",
          "processing cold_season2 samples -----")[i])
rr[[i]] = data.frame()
rr.probe[[i]] = data.frame()
for(j in 1:length(c.deg.select)){
  print(c.deg.select[j])
  subselect=which(info[,"gene"]==c.deg.select[j])
  subinformation<-info[subselect,]
  subgene.C<-get(paste0("gene",i))[control[[i]],subselect]
  subgene.W<-get(paste0("gene",i))[warming[[i]],subselect]
  rr[[i]] = rbind(rr[[i]],as.data.frame(respr(subgene.W,subgene.C,fill.zero = c("no"), mean.based = c("mean"))))
  
  for (k in 1:length(subselect)){
    ms = c()
    ms = as.data.frame(respr(subgene.W[,k],subgene.C[,k],fill.zero = c("no"), mean.based = c("mean")))
    ms = data.frame(ms,colnames(gene)[subselect][k], c.deg.select[j])
    colnames(ms)[16:17] = c("probe","gene")
    rr.probe[[i]] = rbind(rr.probe[[i]],ms)
  }
  }
row.names(rr[[i]]) = c.deg.select
}
save(rr,file = "datasets/probe/c.deg.rr.3season.Rdata")
save(rr.probe,file = "datasets/probe/c.deg.rr.probe.3season.Rdata")

load("datasets/probe/c.deg.rr.3season.Rdata")
season = c("Cold_season_1","Warm_season","Cold_season_2")
for (i in 1:3){
  print(season[i])
  print(mean(rr[[i]]$RR))
  print(std.error(rr[[i]]$RR))
rr[[i]]$gene = c.deg.select
rr[[i]]$sig = c("")
rr[[i]] = sig_label(rr[[i]],col = 3)

rr[[i]]$gene =  c("amyA","cda","glucoamylase","nplT","pula",#Starch
              "ara","mannanase","xyla","xylanase",#Hemicellulose
              "cellobiase","endoglucanase","exoglucanase",#Cellulous
              "chitinase",#Chitin
              "pectate lyase",#Pectin
              "glx","mnp","phenol oxidase",#lignin
              "cutinase",#cutin
              "cdh")#Terpene
rr[[i]]$gene =  factor(rr[[i]]$gene, levels = c("cdh",#Terpene
                                        "cutinase",#cutin
                                        "glx","mnp","phenol oxidase",#lignin
                                        "pectate lyase",#Pectin
                                        "chitinase",#Chitin
                                        "cellobiase","endoglucanase","exoglucanase",#Cellulous
                                        "ara","mannanase","xyla","xylanase",#Hemicellulose
                                        "amyA","cda","glucoamylase","nplT","pula"))#Starch

#plot response ratios
xx = ggplot(data = rr[[i]], aes(x = gene, y = RR))+ 
  geom_errorbar(aes(ymin=CI95.down, ymax=CI95.up),position = position_dodge(1), width = 0.5, size = 0.2)+#, width=.2, position=position_dodge(0.01)
  geom_point(color = "black",size = 3, pch =21, fill = "brown", stroke = 0.2)+
  geom_text(aes(y = CI95.up + 0.03, label = sig), size = 5)+
  scale_y_continuous(limits = c(-0.2,0.3))+
  theme_test()+
  geom_hline(yintercept =  0, linetype="dashed")+
  ylab("Response ratio")+
  xlab("C degradation genes")+
  coord_flip()+
  theme(axis.title.y = element_blank(),axis.text = element_text(color = "black",family = "Arial", size = 14),
        axis.title.x = element_text(color = "black",family = "Arial", size = 16),
        axis.title.x.bottom  = element_text(margin = margin(5, 0, 15, 0)))
xx
ggsave(paste0("figures/response_ratio/c.deg.response.ratio.95.",season[i],".no.mean.pdf"),xx, dpi = 500, width = 5, height = 7)
}

###By 12 months ================================================================
control = c(1:4)
warming = c(5:8)
rownames(gene)[control]
rownames(gene)[warming]

#c.deg.genes
rr = list()
rr.probe = list()
for (i in 1:12){
  print(c("processing Jan samples -----",
          "processing Feb samples -----",
          "processing Mar samples -----",
          "processing Apr samples -----",
          "processing May samples -----",
          "processing Jun samples -----",
          "processing Jul samples -----",
          "processing Aug samples -----",
          "processing Sep samples -----",
          "processing Oct samples -----",
          "processing Nov samples -----",
          "processing Dec samples -----"
          )[i])
  month.gene = gene[c((8*i-7):(8*i)),]
  rr[[i]] = data.frame()
  rr.probe[[i]] = data.frame()
  for(k in 1:length(c.deg.select)){
    print(c.deg.select[k])
    subselect=which(info[,"gene"]==c.deg.select[k])
    subinformation<-info[subselect,]
    subgene.C<-month.gene[control,subselect]
    subgene.W<-month.gene[warming,subselect]
    rr[[i]] = rbind(rr[[i]],as.data.frame(respr(subgene.W,subgene.C,fill.zero = "no", mean.based = "mean")))
    
    for (l in 1:length(subselect)){
      ms = c()
      ms = as.data.frame(respr(subgene.W[,l],subgene.C[,l],fill.zero = "no", mean.based = "mean"))
      ms = data.frame(ms,colnames(gene)[subselect][l], c.deg.select[k])
      colnames(ms)[16:17] = c("probe","gene")
      rr.probe[[i]] = rbind(rr.probe[[i]],ms)
    }
  }
  row.names(rr[[i]]) = c.deg.select
  }
save(rr, file = "datasets/probe/c.deg.rr.month.Rdata")
save(rr.probe,file = "datasets/probe/c.deg.rr.probe.month.Rdata")

load("datasets/probe/c.deg.rr.month.Rdata")
month = as.character(c(1:12))
for (i in 1:12){
  print(month[i])
  print(weighted.mean(rr[[i]]$RR,rr[[i]]$size.x))
  print(std.error(rr[[i]]$RR))
  rr[[i]]$gene = c.deg.select
  rr[[i]]$sig = c("")
  rr[[i]] = sig_label(rr[[i]],col = 3)
  
  rr[[i]]$gene =  c("amyA","cda","glucoamylase","nplT","pula",#Starch
                    "ara","mannanase","xyla","xylanase",#Hemicellulose
                    "cellobiase","endoglucanase","exoglucanase",#Cellulous
                    "chitinase",#Chitin
                    "pectate lyase",#Pectin
                    "glx","mnp","phenol oxidase",#lignin
                    "cutinase",#cutin
                    "cdh")#Terpene
  rr[[i]]$gene =  factor(rr[[i]]$gene, levels = c("cdh",#Terpene
                                                  "cutinase",#cutin
                                                  "glx","mnp","phenol oxidase",#lignin
                                                  "pectate lyase",#Pectin
                                                  "chitinase",#Chitin
                                                  "cellobiase","endoglucanase","exoglucanase",#Cellulous
                                                  "ara","mannanase","xyla","xylanase",#Hemicellulose
                                                  "amyA","cda","glucoamylase","nplT","pula"))#Starch
  
  xx = ggplot(data = rr[[i]], aes(x = gene, y = RR))+ 
    geom_errorbar(aes(ymin=CI95.down, ymax=CI95.up),position = position_dodge(1), width = 0.5, size = 0.2)+#, width=.2, position=position_dodge(0.01)
    geom_point(color = "black",size = 3, pch =21, fill = "brown", stroke = 0.2)+
    geom_text(aes(y = CI95.up + 0.03, label = sig), size = 5)+
    scale_y_continuous(limits = c(-0.3,0.4))+
    theme_test()+
    geom_hline(yintercept =  0, linetype="dashed")+
    ylab("Response ratio")+
    xlab("C degradation genes")+
    coord_flip()+
    theme(axis.title.y = element_blank(),axis.text = element_text(color = "black",family = "Arial", size = 14),
          axis.title.x = element_text(color = "black",family = "Arial", size = 18))
  xx
  ggsave(paste0("figures/response_ratio/month/c.deg.response.ratio.95.",month[i],".no.mean.pdf"),xx, dpi = 500, width = 5, height = 7)
}


#categorize RR into 4 types
##3 seasons ===================================================================
#C.deg.genes
load("datasets/probe/c.deg.rr.probe.3season.Rdata")
c.deg.probe.cold.1 = rr.probe[[1]]
c.deg.probe.hot = rr.probe[[2]]
c.deg.probe.cold.2 = rr.probe[[3]]

#summary generation
c.deg.probe.cold.1.sig = filter(c.deg.probe.cold.1,P <= 0.05 | is.na(P))
c.deg.probe.cold.1.sig = filter(c.deg.probe.cold.1.sig, size.x + size.y > 0)
c.deg.probe.cold.1.sig.stimulated = filter(c.deg.probe.cold.1.sig, RR>0)
c.deg.probe.cold.1.sig.suppressed = filter(c.deg.probe.cold.1.sig, RR<0)
c.deg.probe.cold.1.sig.w.uniq = filter(c.deg.probe.cold.1.sig, (size.x > 0) & (size.y == 0))
c.deg.probe.cold.1.sig.c.uniq = filter(c.deg.probe.cold.1.sig, (size.x == 0) & (size.y > 0))

c.deg.probe.hot.sig = filter(c.deg.probe.hot,P <= 0.05 | is.na(P))
c.deg.probe.hot.sig = filter(c.deg.probe.hot.sig, size.x + size.y > 0)
c.deg.probe.hot.sig.stimulated = filter(c.deg.probe.hot.sig, RR>0)
c.deg.probe.hot.sig.suppressed = filter(c.deg.probe.hot.sig, RR<0)
c.deg.probe.hot.sig.w.uniq = filter(c.deg.probe.hot.sig, (size.x > 0) & (size.y == 0))
c.deg.probe.hot.sig.c.uniq = filter(c.deg.probe.hot.sig, (size.x == 0) & (size.y > 0))

c.deg.probe.cold.2.sig = filter(c.deg.probe.cold.2,P <= 0.05 | is.na(P))
c.deg.probe.cold.2.sig = filter(c.deg.probe.cold.2.sig, size.x + size.y > 0)
c.deg.probe.cold.2.sig.stimulated = filter(c.deg.probe.cold.2.sig, RR>0)
c.deg.probe.cold.2.sig.suppressed = filter(c.deg.probe.cold.2.sig, RR<0)
c.deg.probe.cold.2.sig.w.uniq = filter(c.deg.probe.cold.2.sig, (size.x > 0) & (size.y == 0))
c.deg.probe.cold.2.sig.c.uniq = filter(c.deg.probe.cold.2.sig, (size.x == 0) & (size.y > 0))

df = data.frame(Group = c("Unique under warming","Increased by warming","Decreased by warming", "Unique under unwarming"),
                Early_cool_season = c(nrow(c.deg.probe.cold.1.sig.w.uniq)/nrow(c.deg.probe.cold.1.sig),
                                      nrow(c.deg.probe.cold.1.sig.stimulated)/nrow(c.deg.probe.cold.1.sig),
                                      nrow(c.deg.probe.cold.1.sig.suppressed)/nrow(c.deg.probe.cold.1.sig),
                                      nrow(c.deg.probe.cold.1.sig.c.uniq)/nrow(c.deg.probe.cold.1.sig)),
                Warm_season = c(nrow(c.deg.probe.hot.sig.w.uniq)/nrow(c.deg.probe.hot.sig),
                                nrow(c.deg.probe.hot.sig.stimulated)/nrow(c.deg.probe.hot.sig),
                                nrow(c.deg.probe.hot.sig.suppressed)/nrow(c.deg.probe.hot.sig),
                                nrow(c.deg.probe.hot.sig.c.uniq)/nrow(c.deg.probe.hot.sig)),
                Late_cool_season = c(nrow(c.deg.probe.cold.2.sig.w.uniq)/nrow(c.deg.probe.cold.2.sig),
                                     nrow(c.deg.probe.cold.2.sig.stimulated)/nrow(c.deg.probe.cold.2.sig),
                                     nrow(c.deg.probe.cold.2.sig.suppressed)/nrow(c.deg.probe.cold.2.sig),
                                     nrow(c.deg.probe.cold.2.sig.c.uniq)/nrow(c.deg.probe.cold.2.sig)))
df[,2:4] = df[,2:4]*100


df2 = melt(df)
df2[,3] = df2[,3]*100
df2$Group = factor(df2$Group, levels = c("Unique under warming","Increased by warming","Decreased by warming", "Unique under unwarming"))


link_dat <- df %>% 
  arrange(by=desc(Group)) %>% 
  mutate_if(is.numeric, cumsum) 
bar.width <- 0.5
link_dat <- link_dat[, c(1,2,rep(3:(ncol(link_dat)-1),each=2), ncol(link_dat))]
link_dat <- data.frame(y=t(matrix(t(link_dat[,-1]), nrow=2)))
link_dat$x.1 <- 1:(ncol(df)-2)+bar.width/2
link_dat$x.2 <- 1:(ncol(df)-2)+(1-bar.width/2)


ggplot(df2,aes(x=variable,y=value,fill=Group))+
  geom_bar(stat = "identity", position = "stack", width = bar.width, color = "black", alpha = 0.9, size = 0.25)+
  theme_test() +
  theme(axis.text = element_text(color = "black",family = "Arial", size = 10),
        axis.title.y = element_text(color = "black",family = "Arial", size = 10),
        panel.border = element_rect(size = 0.5),
        axis.text.x = element_text(vjust = -1.6)) +
  # scale_fill_manual(values=colors)+
  labs(x=" ",y="Significant C decomposition genes(%)",fill="Group") +
  scale_fill_manual(values = c("indianred4","indianred3","lightskyblue1","lightskyblue3")) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(labels=c("Cold season 1","Warm season","Cold season 2"))
# +geom_segment(data=link_dat,aes(x=x.1, xend=x.2, y=y.1, yend=y.2), inherit.aes = F)

ggsave("figures/response_ratio/percentage_stimulated_4_cat.pdf", dpi = 500, width = 5, height = 2.5)
