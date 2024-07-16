rm(list = ls())
gc()

setwd("/export3/zhangw/Project_Cross/Project_MEST/")
# 利用胚层发育的基因筛选出MEST

####### 1基因数据集合处理 #########

library(clusterProfiler)
library(tidyverse)

gene.path = './data/genelist/'
files.gene = list.files(gene.path)

germ_layer = lapply(files.gene, function(x){
  
  inner.gene = read.gmt(paste0(gene.path,x))
  return(inner.gene)
  }) %>% do.call(rbind,.)

table(germ_layer$term)

f.germ.layer.gene = unique(germ_layer$gene)
f.germ.layer.gene = gsub('-','.',f.germ.layer.gene)
# 删掉在队列种不表达的
save(f.germ.layer.gene,file = './data/f_germ_layer_genelist.Rdata')

write_csv(germ_layer,file = './data/f_germ_layer_genelist.csv')

####### 2说明胚层发育的基因很重要##########

# ##### 2.1通过GSVA 计算 这些胚层发育基因的活性， 越高表示这些基因代表的通路越活越 ########
setwd("/export3/zhangw/Project_Cross/Project_MEST/")

load("/export3/zhangw/Project_Cross/Project_Mime/Proj/data/glioma.cohort.Rdata")
load('./data/f_germ_layer_genelist.Rdata')
genelist = list(germ_layer = f.germ.layer.gene)

res.gsva.germ.layer = lapply(list_train_vali_Data, function(x){
  tmp = x[,-c(1:3)]
  rownames(tmp) = x$ID
  rs = GSVA::gsva(expr = as.matrix(t(tmp)),
                  gset.idx.list = genelist,method = 'gsva',kcdf ='Gaussian', parallel.sz=48)
  rs = rs %>% t() %>% as.data.frame()
 return(rs)
})
names(res.gsva.germ.layer) = names(list_train_vali_Data)

# save(res.gsva.germ.layer, file = './res/res.gsva.glioma.cohort.germ.layer.Rdata')
load('./res/res.gsva.glioma.cohort.germ.layer.Rdata')

plist = lapply(names(res.gsva.germ.layer), function(x){
  tmp.g = res.gsva.germ.layer[[x]] %>% rownames_to_column('ID')
  tmp.sinfo = list_train_vali_Data[[x]][,1:3]
  tmp = left_join(tmp.g,tmp.sinfo)
  tmp=tmp[,-1]
  colnames(tmp)[1] = 'var'
  mySurv <- Surv(tmp$OS.time, tmp$OS)
  cutoff<-0.5
  value <- quantile(tmp$var, probs = c(cutoff))
  tmp$Group<-ifelse(tmp$var> value,"High","Low")
  Group <- tmp$Group
  Group <- factor(Group, levels = c("Low", "High"))
  tmp$Group <- factor(Group, levels = c("Low", "High"))
  fit <- survfit(Surv(OS.time, OS) ~ Group, data = tmp)
  
  # calculate HR and 95%CI
  data.survdiff <- survdiff(mySurv ~ Group)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
  CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
  
  p1<-ggsurvplot(fit, data = tmp ,
                 conf.int = F, 
                 censor = F, 
                 palette = c('#375631','#6D4E7E') , 
                 legend.title = 'Development of Germ layer',
                 font.legend = 11,
                 surv.median.line = 'hv',
                 legend.labs=c(paste0("Low","(n=",fit$n[1],")"),
                               paste0("High","(n=",fit$n[2],")")),
                 pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                            paste("P = ",round(p.val,3), sep = "")),
                              HR, CI, sep = "\n"),
                 xlab="Day",
                 # title= paste0("Expression in ",dataset),
                 title= paste0("",x),
                 pval.coord=c(1000,0.9))
  p1<-p1$plot+
    labs(subtitle = 'Cutoff: medium')
  return(p1)
})

survplot = plist

survplot[[1]]

cairo_pdf('./fig/fig_1a_germ_layer_glioma_sur_os_km.pdf',width = 15,height = 20,onefile = F)
survplot[[1]]+survplot[[2]]+survplot[[3]]+survplot[[4]]+survplot[[5]]+survplot[[6]]+
  survplot[[7]]+survplot[[8]]+survplot[[9]]+survplot[[10]]+survplot[[11]]+
  plot_layout(ncol = 3)
dev.off()

# 看这个信号在不同亚类中的分布 主要是肿瘤和正常组织之间的差别
your_color_pal=c('#375631','#D3C3D0','#A786BA','#6D4E7E',paletteer_c("grDevices::Sunset", 30)[c(6,12,18,24,30)] )


load('./data/f_germ_layer_genelist.Rdata')
load("/export3/zhangw/Project_Cross/Project.Glioma.SGene/data/glioma.tcga.normal.brain.cortex.gtex.expr.Rdata")
gtex.tcga.glioma =df
colnames(gtex.tcga.glioma) = gsub('-','.',colnames(gtex.tcga.glioma))
df= gtex.tcga.glioma %>% dplyr::select(c('sample',"type2"),everything())
sample.type = gtex.tcga.glioma[,c(1,4)]
df= df[,-2]
df = avereps(df,ID = df$sample)
df= df %>% as.data.frame()
rownames(df) = df$sample
df = df %>% as.data.frame()
rownames(df) = df$sample  
df = df[,-1]
df = df %>% t() %>% as.data.frame()
sample.type = sample.type %>% dplyr::filter(!duplicated(sample.type$sample))
tmp.phe = readRDS('/export/bioinfo-team/home/liuhw/bioinfo_mill/dataset/tcga_gtex/processeddata/tcga_gtex_phe.Rds')
noraml.id = tmp.phe[tmp.phe$detailed_category =='Brain - Cortex',]$sample
noraml.id = gsub('-','.',noraml.id)
tmp = df[,c(sample.type[sample.type$type2=='tumor',]$sample,noraml.id)]
tmp = tmp %>% as.data.frame()
tmp = apply(tmp, 2, as.numeric)

tmp = tmp %>% as.data.frame()
rownames(tmp) = rownames(df)

rs = GSVA::gsva(expr = as.matrix(tmp),
                gset.idx.list = genelist,method = 'gsva',kcdf ='Gaussian', parallel.sz=48)
rs = rs %>% t() %>% as.data.frame() %>% rownames_to_column('ID')
tnp = data.frame(ID = c(sample.type[sample.type$type2=='tumor',]$sample,noraml.id),
                 type = c(rep('Tumor', length(sample.type[sample.type$type2=='tumor',]$sample)),
                          rep('Normal', length(noraml.id))))
df.plot = left_join(rs,tnp)
df.plot =df.plot[,-1]
colnames(df.plot) = c('var','group')

df.plot$group = factor(df.plot$group,levels = c('Normal','Tumor'))
my_comparisons <- list(c('Normal','Tumor'))


p1=ggplot(df.plot, aes(group, var))+ 
  # 添加散点图
  geom_beeswarm(size=2.5,aes(color=group),cex =1.5)+
  stat_summary(fun.y = "mean", geom = "crossbar", #添加均值
               mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.15,size=0.3,color="black")+
  stat_summary(fun.data = "mean_se", geom="errorbar", width = 0.1,size=0.5,color="black")+#添加误差线
  theme_classic()+scale_color_manual(values =c('#375631','#D3C3D0') )+
  
  stat_compare_means(comparisons=my_comparisons,method = "t.test",hide.ns = F,label = "p.format")+
  labs(title = 'TCGA+GTEx',subtitle = "T(689) vs N(105)")+
  ylab("GSVA score of germ layer") +
  xlab("") +
  #ylim(35,37.5)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        legend.position = "none",
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
p1

ggsave("./fig/fig_1a_germ_layer.normal.tumor.gtex.tcga.glioma.pdf",width = 2.5,height = 4.8)



### ### 2.2  这里进行差异分析 和单因素回归分析  ##############
# 加入GTEx的数据 

load('./data/f_germ_layer_genelist.Rdata')

load("/export3/zhangw/Project_Cross/Project.Glioma.SGene/data/glioma.tcga.normal.brain.cortex.gtex.expr.Rdata")
load("/export3/zhangw/Project_Cross/Project_CNGA3/data/protein.coding.genes.Rdata")
gtex.tcga.glioma =df
colnames(gtex.tcga.glioma) = gsub('-','.',colnames(gtex.tcga.glioma))
df= gtex.tcga.glioma[,c('sample',"type2",intersect(gsub('-','.', f.germ.layer.gene),colnames(gtex.tcga.glioma)))]
sample.type = gtex.tcga.glioma[,c(1,4)]
df= df[,-2]
df = avereps(df,ID = df$sample)
df= df %>% as.data.frame()
rownames(df) = df$sample
df = df %>% as.data.frame()
rownames(df) = df$sample  

df = df[,-1]
source('~/R/ZW/code_sum/two_subtype_compr_logfc.R')  
df = df %>% t() %>% as.data.frame()
sample.type = sample.type %>% dplyr::filter(!duplicated(sample.type$sample))
tmp.phe = readRDS('/export/bioinfo-team/home/liuhw/bioinfo_mill/dataset/tcga_gtex/processeddata/tcga_gtex_phe.Rds')
noraml.id = tmp.phe[tmp.phe$detailed_category =='Brain - Cortex',]$sample
noraml.id = gsub('-','.',noraml.id)

tmp = df[,c(sample.type[sample.type$type2=='tumor',]$sample,noraml.id)]
tmp = tmp %>% as.data.frame()
tmp = apply(tmp, 2, as.numeric)
tmp = 2**tmp
tmp = tmp %>% as.data.frame()
rownames(tmp) = rownames(df)

tmp.t =tmp
# germ layer gene 
tmp = tmp.t[f.germ.layer.gene,]
tmp[is.na(tmp)]=1
deg.res = two_subtype_logFC_Cpr_wilcox(expr = tmp,
                                       treat_list = data.frame(ID=sample.type[sample.type$type2=='tumor',]$sample),
                                       ctrl_list = data.frame(ID=noraml.id),method = 'w')  

# save(deg.res,file = './res/deg.res.germ.layre.glioma.gtex.Rdata')

# unicox 在tcga中
load("/export3/zhangw/Project_Cross/Project_Mime/Proj/data/glioma.cohort.Rdata")
load('./data/f_germ_layer_genelist.Rdata')

source('/export3/zhangw/Code.sum/SigUicox.R')
# 将这个函数的返回改成返回表格
res.unicox = SigUnicox(gene_list = f.germ.layer.gene,inputSet =list_train_vali_Data$TCGA,unicox_pcutoff = 0.05)

# save(res.unicox,file = './res/res.unicox.germ.layre.glioma.tcga.Rdata')
load('./res/deg.res.germ.layre.glioma.gtex.Rdata')
load('./res/res.unicox.germ.layre.glioma.tcga.Rdata')

df.unicox =res.unicox %>% na.omit()
df.deg = deg.res %>% na.omit()
rownames(df.unicox)= df.unicox$gene
com_gene = intersect(rownames(df.deg),df.unicox$gene)
colnames(df.unicox)
df.plot = cbind(df.deg[com_gene, c('log2FoldChange','padj')],
                df.unicox[com_gene,c('HR','pvalue')])
df.plot$HR = log2(df.plot$HR)
colnames(df.plot) = c('logFC','DEadjP','logHR','UnicoxP')
df.plot$logHR[df.plot$logHR >10]=0
df.plot$type = case_when(df.plot$DEadjP <0.05&df.plot$logFC >1 & df.plot$UnicoxP <0.05&df.plot$logHR>0~'Up in glioma&Risky',
                         df.plot$DEadjP <0.05&df.plot$logFC <  -1 & df.plot$UnicoxP <0.05&df.plot$logHR>0~'Down in glioma&Risky',
                         df.plot$DEadjP <0.05&df.plot$logFC > 1 & df.plot$UnicoxP <0.05&df.plot$logHR<0~'Up in glioma&Protective',
                         df.plot$DEadjP <0.05&df.plot$logFC <  -1 & df.plot$UnicoxP <0.05&df.plot$logHR<0~'Down in glioma&Protective',
                         .default = 'NotSig'
                         )
df.plot['MEST',]
df.plot %>% 
  filter(-log10(DEadjP)>20 & logFC >2 &-log10(UnicoxP) > 20 &logHR>0.6) -> new.text.label
new.text.label$lab= rownames(new.text.label)

ggplot(data=df.plot,aes(x=logHR,y=logFC))+
  geom_point(aes(color=type))+
  scale_color_manual(values = c("Up in glioma&Risky"="#3C224B",
                                "Down in glioma&Risky"="#A786BA",
                                "Up in glioma&Protective"="#83BC8B",
                                "Down in glioma&Protective"="#375631",
                                "NotSig"="#aaaaaa"
                                )
                    )+
  geom_text_repel(data=new.text.label,
                  aes(x=logHR,y=logFC,
                      label=lab),
                  color = "#3C224B",
                  nudge_x = 0.7,
                  box.padding = 0.8,
                  nudge_y = 0.5,
                  segment.curvature = -0.1,
                  segment.ncp = 3,
                  segment.angle = 20)+
  geom_hline(yintercept = c(1,-1), linetype = "dashed") +
  geom_vline(xintercept = c(0), linetype = "dashed") +
theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.title = element_text(face = "bold", color = "black", size = 10),
        axis.text = element_text(color = "black", size = 9, face = "bold"),
        legend.background = element_blank(),
        legend.title = element_text(face = "bold", color = "black", size = 10),
        legend.text = element_text(face = "bold", color = "black", size = 9),
        legend.spacing.x = unit(0, "cm"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle  = element_text(hjust = 0.5),
        legend.position = c('top')  # 设置图例位置为右上角
  )+labs(color='Groups',title = 'DEG & Univariate Cox regression',subtitle = 'TCGA GBM-LGG & GTEx')+
  guides(color  = guide_legend(nrow = 3))
  
ggsave('./fig/fig_1a_deg_unicox_glioma_gtex_germ_layer_gene.pdf',width = 5,height = 6)

library(ggpie)

ggpie::ggpie(data = as.data.frame(df.plot), group_key = "type", count_type = "full",
      label_info = "all",     # 标签内容：默认"count", "ratio", "all"
      label_type = "horizon", # 标签形式："none", "circle", "horizon"
      label_split = NULL,
      label_size = 4,
      label_pos = "out") + # 设置阈值
  ggtitle("DEG&UniCox")+
  scale_fill_manual(values = alpha(c('#aaaaaa','#3C224B','#83BC8B','#375631','#A786BA'),0.7))+
  theme(        plot.title = element_text(hjust = 0.5))

ggsave('./fig/fig_1a_pie_deg_unicox_glioma_gtex_germ_layer_gene.pdf',width = 5,height = 6)


######  2.3  这里进行turkey'test 分析 看有多少个基因在这里异常高表达 ##############
rm(list = ls())
load('./res/glioma_switch_genes_fc.Rdata')

load('/export3/zhangw/Project_Cross/Project_MEST/data/f_germ_layer_genelist.Rdata')  




# 
load("/export3/zhangw/pancancer/tcga.mrna/dataset/tcga.pancancer.profile.surv.Rdata")
expr= read_rds('/export3/zhangw/Project_Cross/Project_GPCR_pancancer/data/pancancer.tcga.mrna.exp.not.log.not.Combat.move.rds') 
load('/export3/zhangw/Project_Cross/Project_GPCR_pancancer/data/pheno_pancancer.tcga.mrna.Combat.move.Rdata')  

tumor.type = unique(samAnno2$`cancer type`)

list_train_vali_Data = list()
profile_pan=expr  %>% as.data.frame()

for (i in tumor.type) {
  test.id = samAnno2[samAnno2$`cancer type`==i,]$simple_barcode
  
  list_train_vali_Data[[i]] = profile_pan[test.id,]
  
}
list_train_vali_Data[['pancancer']]=profile_pan

tumor_type = unique(res_fc$cohort)

df.c.num  = data.frame(cohort= NULL, sum = NULL)

for (i in tumor_type) {
  df.c.num =rbind(df.c.num, data.frame(cohort= i, sum = nrow(list_train_vali_Data[[i]])))
  
}


df2= left_join(df.c.num,res_fc)
df3 = df2 %>% dplyr::mutate(per = num.ab/sum) %>% 
  dplyr::select(c('cohort','genelist', 'FC','per','median_value'))
df3$FC [which(is.na(df3$FC))]= 1 ### 将没有意义的 取为0
df3$per= 100*df3$per
df3$FC = log1p(df3$FC)
# show_type = c('pancancer','LGG','GBM')
show_type = tumor_type
show_genelist = unique(df3$genelist)
plist = list()

library(circlize)
# cold <- colorRampPalette(c('#A0C3E0','#41b6c4','#253494'))
cold <- colorRampPalette(c('#F2ECE3','#A786BA','#3C224B'))
warm <- colorRampPalette(c('#ffffb2','#fecc5c','#e31a1c'))
mypalette <- c(rev(cold(11)), warm(10))

# install.packages('ggrastr')
library('ggrastr')

for (i in show_type) {
  pplist = list()
  for (j in show_genelist) {
    df4= df3[df3$cohort%in% i &df3$genelist ==j,] 
    
    ggplot(df4, aes(per, FC)) +
      stat_density_2d(aes(fill = ..level..), geom = "polygon")+
      
      # geom_point_rast(shape = 21, stroke=0.25,
      #                 aes(colour=median_value
      #                     ), size = 2) +
      geom_density_2d(data=df4,
                      aes(x=per, y=FC),
                      bins = 5, colour="white") +
      scale_fill_gradientn(colours = cold(10),
                           limits = c(0,2), 
                           values = c(0,0.05,0.1,0.2,1,2),breaks = c(0,0.05,0.1,0.2,1,2)
      )+
      
      # scale_fill_gradientn(colours =  warm(10))+
      
      # scale_colour_gradientn(colours = cold(11))+
      theme_classic()+ggtitle(i)+
      geom_vline(xintercept =  median(df4$per),linetype = 'dashed',color= '#9D7660' )+
      geom_hline(yintercept =  median(df4$FC),linetype = 'dashed',color= '#9D7660')+
      geom_text(x = median(df4$per)*1.05, y = median(df4$FC)*0.95, label = round(median(df4$per),2))+
      geom_text(x = median(df4$per)*0.95, y = median(df4$FC)*1.05, label =  round(median(df4$FC),2))+
      labs(x='Per', y='log(FC+1)')+
      theme(
        plot.title = element_text(size=12, face="bold.italic", hjust = 0),
        axis.text=element_text(size=8, colour = "black"),
        axis.title=element_text(size=12),
        legend.text = element_text(size =10),
        # legend.title=element_blank(),
        aspect.ratio=1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      ) ->pp
    
    
    
    pp
    pplist[[j]] = pp
    
    
  }
  
  plist[[i]] = pplist[[1]]
   
  
}

pdf('./fig/fig_1a_fc_per_pan_glioma.pdf',width = 18,height = 15,onefile = F)
plist[[1]]+plist[[2]]+plist[[3]]+    
plist[[4]]+plist[[5]]+plist[[6]]+    
plist[[7]]+plist[[8]]+plist[[9]]+    
plist[[10]]+plist[[11]]+plist[[12]]+    
plist[[13]]+plist[[14]]+plist[[15]]+    
plist[[16]]+plist[[17]]+plist[[19]]+    
plist[[19]]+plist[[20]]+plist[[22]]+    
plist[[22]]+plist[[23]]+plist[[24]]+    
plist[[25]]+plist[[26]]+plist[[27]]+    
plist[[28]]+plist[[29]]+plist[[30]]+    
  plot_layout(nrow =5,guides='collect' )&theme(legend.position='right')

dev.off()

# 统计 switch gene 的数目

plist = list()
for (i in show_type) {
  
  num.swith = nrow( res_switch_gene[res_switch_gene$cohort ==i,])
  df.plot = data.frame(Type = c(rep('Switch gene',num.swith),
                                rep('Not', length(f.germ.layer.gene)-num.swith)))
  
  p1 = ggpie::ggdonut(data = as.data.frame(df.plot), group_key = "Type", count_type = "full",
                    label_info = "all",     # 标签内容：默认"count", "ratio", "all"
                    label_type = "horizon", # 标签形式："none", "circle", "horizon"
                    label_split = NULL,
                    donut.label.color = "#3C224B",
                    donut.label	=F,
                    
                    label_size = 4,
                    label_pos = "out") + # 设置阈值
    ggtitle(paste0(i))+
    scale_fill_manual(values = alpha(c('#A786BA','#aaaaaa'),0.7))+
    theme(        plot.title = element_text(hjust = 0.5))
  plist[[i]]=p1
  
  
}

pdf('./fig/fig_1b_pie_switch_Gene_TCGA_pancancer_glioma.pdf',width = 18,height = 15,onefile = F)
plist[[1]]+plist[[2]]+plist[[3]]+    
  plist[[4]]+plist[[5]]+plist[[6]]+    
  plist[[7]]+plist[[8]]+plist[[9]]+    
  plist[[10]]+plist[[11]]+plist[[12]]+    
  plist[[13]]+plist[[14]]+plist[[15]]+    
  plist[[16]]+plist[[17]]+plist[[19]]+    
  plist[[19]]+plist[[20]]+plist[[22]]+    
  plist[[22]]+plist[[23]]+plist[[24]]+    
  plist[[25]]+plist[[26]]+plist[[27]]+    
  plist[[28]]+plist[[29]]+plist[[30]]+    
  plot_layout(nrow =5,guides='collect' )&theme(legend.position='right')
dev.off()
######  2.4  单细胞数据集中评估肿瘤细胞和别的细胞的germ layer 的分数 ##############
rm(list = ls())
library(Seurat)
gc()
plot_theme    <- theme_bw(base_size = 12) + theme(axis.title = element_text(size = 12),
                                                  axis.text = element_text(size = 12),
                                                  panel.background = element_rect(fill = "transparent"),
                                                  axis.line = element_blank(),
                                                  strip.background = element_blank(),
                                                  panel.grid.major = element_blank(),
                                                  panel.grid.minor = element_blank(),
                                                  panel.border = element_blank(),
                                                  axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
                                                  axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"))
plot_thme_blank =  theme_classic()+theme(axis.title.x=element_blank(),
                       axis.text.x=element_blank(),
                       axis.ticks.x=element_blank(),
                       axis.title.y=element_blank(),
                       axis.text.y=element_blank(),
                       axis.ticks.y=element_blank(),
                       line = element_blank())

load('/export3/zhangw/Project_Cross/Project_MEST/data/f_germ_layer_genelist.Rdata')  


#verhauk 2021 ng -----------------

panglioma_tumor_verhaak <- readRDS("/export3/zhangw/Project_Cross/Project_scGBM_Liuhw/2021_verhauk/GBM_scRNA/analysis_scRNAseq_tumor_PanGlioma.rds")
scgp <- panglioma_tumor_verhaak
## Change to cell states already defined.
Idents(scgp) <- scgp@meta.data$cell_state
## Create desired levels of cell states defined with marker gene expression.
my_levels = c("Stem-like", "Prolif. stem-like", "Diff.-like",
              "Oligodendrocyte",
              "Fibroblast",
              "Endothelial", "Pericyte",
              "B cell", "Granulocyte", "T cell", "Dendritic cell", "Myeloid")
library(forcats)
library(ggpubr)
scgp@active.ident <- factor(x = scgp@active.ident, levels = my_levels)
scgp@active.ident <- fct_relevel(scgp@active.ident ,rev)

plist = list()

plist[[1]]=DimPlot(scgp, reduction = "umap",group.by = "cell_state",label = F)+labs(title="")+
  labs(title="Verhaak_2021_NG")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_fill_manual(values=c("B cell" = "#eff3ff", "Granulocyte" = "#bdd7e7", "T cell" = "#6baed6", "Dendritic cell" = "#3182bd", "Myeloid" = "#08519c",
                             "Oligodendrocyte" = "#2ca25f",
                             "Endothelial" = "#ffffd4", "Pericyte" = "#fee391",
                             "Fibroblast" = "#feb24c",
                             "Stem-like" = "#fb6a4a", "Diff.-like" = "#fcbba1", "Prolif. stem-like" = "#a50f15"))+
  scale_color_manual(values=c("B cell" = "#eff3ff", "Granulocyte" = "#bdd7e7", "T cell" = "#6baed6", "Dendritic cell" = "#3182bd", "Myeloid" = "#08519c",
                              "Oligodendrocyte" = "#2ca25f",
                              "Endothelial" = "#ffffd4", "Pericyte" = "#fee391",
                              "Fibroblast" = "#feb24c",
                              "Stem-like" = "#fb6a4a", "Diff.-like" = "#fcbba1", "Prolif. stem-like" = "#a50f15")) 

plist
genes = f.germ.layer.gene
names(genes) = 'Germ_layer'

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("AUCell",force = TRUE)
library(AUCell,lib.loc = "/export/bioinfo-team/home/liuhw/R/x86_64-pc-linux-gnu-library/4.1")
cells_rankings <- AUCell_buildRankings(scgp@assays$RNA@data) 
geneSets  = list(germ_layer_development  = f.germ.layer.gene)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)


##set gene set of interest here for plotting
geneSet <- "germ_layer_development"
aucs <- as.numeric(getAUC(cells_AUC)[geneSet, ])
scgp$AUC  <- aucs

library(ggraph)
df = data.frame(scgp@meta.data, scgp@reductions$umap@cell.embeddings) %>% as.data.frame()
plist[[2]]=ggplot(df, aes(UMAP_1, UMAP_2, color=AUC)) + 
  geom_point( size=0.1) +
  scale_colour_gradientn(colours = c('#375631','white','#3C224B')) + 
  labs(title = "germ_layer_development")+
  theme(plot.title = element_text(hjust = 0.5))+
  plot_thme_blank

plist[[2]]

### 将不同的细胞类型提取出来绘制点图或者小提琴图
colnames(df)
Figure_2b = df %>% dplyr::select('cell_state','AUC')
colnames(Figure_2b) = c('Group', 'Value')
Figure_2bbb = Figure_2b %>% dplyr::group_by(Group) %>% summarise(meandt = mean(Value))

Figure_2bbb = Figure_2bbb %>% arrange(Figure_2bbb$meandt)
Figure_2b$Group = factor(Figure_2b$Group,levels = Figure_2bbb$Group)
my_comparisons = list(c('Stem-like','Prolif. stem-like'),
                      c('Prolif. stem-like','Diff.-like'),
                      c('Stem-like','Diff.-like'))
plt=ggplot(Figure_2b, aes(Group, Value))+ 
  geom_bar(stat="summary",
           width=0.9,#宽度
           size=0.5,color='black', aes(fill= Group), alpha = 0.3)+
  # 添加散点图
  # geom_jitter(aes(x=Group, Value,colour = Group), size = 2.5, width = 0.1, height = 0)+
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +     # 误差棒
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, size = 0.2)+  # 指代均数的水平横线
  theme_classic()+
  scale_fill_manual(values=c("B cell" = "#eff3ff", "Granulocyte" = "#bdd7e7", "T cell" = "#6baed6", "Dendritic cell" = "#3182bd", "Myeloid" = "#08519c",
                             "Oligodendrocyte" = "#2ca25f",
                             "Endothelial" = "#ffffd4", "Pericyte" = "#fee391",
                             "Fibroblast" = "#feb24c",
                             "Stem-like" = "#fb6a4a", "Diff.-like" = "#fcbba1", "Prolif. stem-like" = "#a50f15"))+
  scale_color_manual(values=c("B cell" = "#eff3ff", "Granulocyte" = "#bdd7e7", "T cell" = "#6baed6", "Dendritic cell" = "#3182bd", "Myeloid" = "#08519c",
                              "Oligodendrocyte" = "#2ca25f",
                              "Endothelial" = "#ffffd4", "Pericyte" = "#fee391",
                              "Fibroblast" = "#feb24c",
                              "Stem-like" = "#fb6a4a", "Diff.-like" = "#fcbba1", "Prolif. stem-like" = "#a50f15")) +

  # scale_color_manual(values =c('#375631','#D3C3D0','#3C224B') )+
  stat_compare_means(comparisons=my_comparisons,method = "t.test",hide.ns = F,label = "p.format")+
  
  ylab("AUC of germ layer development") +
  xlab("") +
  #ylim(35,37.5)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        legend.position = "none",
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5))+
  guides(file = NULL,color = NULL)
plt
plist[[3]] = plt
# library(ggpval)
#  p1= add_pval(plt,                   # 添加 p 值到图形中
#           pairs = list(c(2,5),c(2,8),c(5,8)),
#           fold_change = F,heights = c(seq(0.075,0.4,0.015)[1:3]),
#           test='t.test', 
#           alternative='two.sided')
# 
# p1

pdf('./fig/fig_1a_scRNA_verhauk_AUCell_germ_layer.pdf',width = 14,height = 6,onefile = F)
plist[[1]]+
plist[[2]]+
plist[[3]]+plot_layout(nrow =1,heights = c(1,1,1),widths = c(1,1,1),guides='collect' )&
    theme(legend.position='top')
dev.off()


#GSE117891 -----------------
GSE117891_seurat <- readRDS("/export/bioinfo-team/home/liuhw/scRNA_GBM_integration/data/GSE117891/processed_data/GSE117891_seurat.rds")

## rename cell
GSE117891_seurat@meta.data$cell_type2 <- GSE117891_seurat@meta.data$cell_type
GSE117891_seurat@meta.data$cell_type2 <- dplyr::recode(GSE117891_seurat@meta.data$cell_type2,
                                                       "Cilium-related tumor cell"="Tumor cell",
                                                       "Proliferating tumor cell"="Tumor cell",
                                                       "M2b macrophages"="Macrophage",
                                                       "Microglia"="Macrophage")




plist[[1]] = DimPlot(GSE117891_seurat, 
                     reduction = "umap",
                     group.by = "cell_type2",
                     cols =c(  "#eff3ff","#fcbba1", "#a50f15", "#fb6a4a","#bdd7e7","#feb24c", "#6baed6",  "#3182bd",  "#08519c",
                                        "#2ca25f",
                                        "#ffffd4", "#fee391",
                                        
                                        "#fcbba1"  ) ,
                                        label = F)+labs(title="")+
  labs(title="GSE117891")+
  theme(plot.title = element_text(hjust = 0.5))

cells_rankings <- AUCell_buildRankings(GSE117891_seurat@assays$RNA@data) 
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)


geneSet <- "germ_layer_development"
aucs <- as.numeric(getAUC(cells_AUC)[geneSet, ])
GSE117891_seurat$AUC  <- aucs

library(ggraph)
df = data.frame(GSE117891_seurat@meta.data, GSE117891_seurat@reductions$umap@cell.embeddings) %>% as.data.frame()
plist[[2]]=ggplot(df, aes(UMAP_1, UMAP_2, color=AUC)) + 
  geom_point( size=0.1) +
  scale_colour_gradientn(colours = c('#375631','white','#3C224B')) +
  labs(title = "germ_layer_development")+
  theme(plot.title = element_text(hjust = 0.5))+
  plot_thme_blank

plist[[2]]

### 将不同的细胞类型提取出来绘制点图或者小提琴图
colnames(df)
Figure_2b = df %>% dplyr::select('cell_type2','AUC')
colnames(Figure_2b) = c('Group', 'Value')
Figure_2bbb = Figure_2b %>% dplyr::group_by(Group) %>% summarise(meandt = mean(Value))

Figure_2bbb = Figure_2bbb %>% arrange(Figure_2bbb$meandt)
Figure_2b$Group = factor(Figure_2b$Group,levels = Figure_2bbb$Group)

# my_comparisons = list(c('Stem-like','Prolif. stem-like'),
#                       c('Prolif. stem-like','Diff.-like'),
#                       c('Stem-like','Diff.-like'))
plt=ggplot(Figure_2b, aes(Group, Value))+ 
  geom_bar(stat="summary",
           width=0.9,#宽度
           size=0.5,color='black', aes(fill= Group), alpha = 0.3)+
  # 添加散点图
  # geom_jitter(aes(x=Group, Value,colour = Group), size = 2.5, width = 0.1, height = 0)+
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +     # 误差棒
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, size = 0.2)+  # 指代均数的水平横线
  theme_classic()+
  scale_fill_manual(values=c(  "#eff3ff","#fcbba1", "#a50f15", "#fb6a4a","#bdd7e7","#feb24c", "#6baed6",  "#3182bd",  "#08519c",
                                        "#2ca25f",
                                        "#ffffd4", "#fee391",
                                        
                                        "#fcbba1"  ))+
  scale_color_manual(values=c(  "#eff3ff","#fcbba1", "#a50f15", "#fb6a4a","#bdd7e7","#feb24c", "#6baed6",  "#3182bd",  "#08519c",
                                         "#2ca25f",
                                         "#ffffd4", "#fee391",
                                         
                                         "#fcbba1"  )) +
  
  # scale_color_manual(values =c('#375631','#D3C3D0','#3C224B') )+
  # stat_compare_means(comparisons=my_comparisons,method = "t.test",hide.ns = F,label = "p.format")+
  
  ylab("AUC of germ layer development") +
  xlab("") +
  #ylim(35,37.5)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        legend.position = "none",
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5))+
  guides(file = NULL,color = NULL)
plt
plist[[3]] = plt


pdf('./fig/fig_1a_scRNA_GSE117891_AUCell_germ_layer.pdf',width = 14,height = 6,onefile = F)
plist[[1]]+
  plist[[2]]+
  plist[[3]]+plot_layout(nrow =1,heights = c(1,1,1),widths = c(1,1,1),guides='collect' )&
  theme(legend.position='top')
dev.off()

###------------GSE131928 smartseq------------
GSE131928_smart_seurat <- readRDS("/export/bioinfo-team/home/liuhw/scRNA_GBM_integration/result/GSE131928_smart_seurat.rds")



plist[[1]] = DimPlot(GSE131928_smart_seurat, 
                     reduction = "umap",
                     group.by = "cell_type",
                     cols =c(  "#eff3ff","#a50f15","#6baed6",  "#3182bd","#fcbba1",  "#fb6a4a","#bdd7e7","#feb24c",   "#08519c",
                                        "#2ca25f",
                                        "#ffffd4", "#fee391",
                                        
                                        "#fcbba1"  ) ,
                                        label = F)+labs(title="")+
  labs(title="GSE131928")+
  theme(plot.title = element_text(hjust = 0.5))
plist[[1]]
cells_rankings <- AUCell_buildRankings(GSE131928_smart_seurat@assays$RNA@data) 
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)


geneSet <- "germ_layer_development"
aucs <- as.numeric(getAUC(cells_AUC)[geneSet, ])
GSE131928_smart_seurat$AUC  <- aucs

library(ggraph)
df = data.frame(GSE131928_smart_seurat@meta.data, GSE131928_smart_seurat@reductions$umap@cell.embeddings) %>% as.data.frame()
plist[[2]]=ggplot(df, aes(UMAP_1, UMAP_2, color=AUC)) + 
  geom_point( size=0.1) +
  scale_colour_gradientn(colours = c('#375631','white','#3C224B')) + 
  labs(title = "germ_layer_development")+
  theme(plot.title = element_text(hjust = 0.5))+
  plot_thme_blank

plist[[2]]

### 将不同的细胞类型提取出来绘制点图或者小提琴图
colnames(df)
Figure_2b = df %>% dplyr::select('cell_type','AUC')
colnames(Figure_2b) = c('Group', 'Value')
Figure_2bbb = Figure_2b %>% dplyr::group_by(Group) %>% summarise(meandt = mean(Value))

Figure_2bbb = Figure_2bbb %>% arrange(Figure_2bbb$meandt)
Figure_2b$Group = factor(Figure_2b$Group,levels = Figure_2bbb$Group)

my_comparisons = list(c('Oligodendrocyte','Malignant'),
                      c('Oligodendrocyte','T-cell'),
                      c('Oligodendrocyte','Macrophage'),
                      c('Malignant','T-cell'),
                      c('Malignant','Macrophage'),
                      c('T-cell','Macrophage')
                      )
plt=ggplot(Figure_2b, aes(Group, Value))+ 
  geom_bar(stat="summary",
           width=0.9,#宽度
           size=0.5,color='black', aes(fill= Group), alpha = 0.3)+
  # 添加散点图
  # geom_jitter(aes(x=Group, Value,colour = Group), size = 2.5, width = 0.1, height = 0)+
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +     # 误差棒
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, size = 0.2)+  # 指代均数的水平横线
  theme_classic()+
  scale_fill_manual(values=c(  "#eff3ff","#a50f15","#6baed6",  "#3182bd","#fcbba1",  "#fb6a4a","#bdd7e7","#feb24c",   "#08519c",
                                        "#2ca25f",
                                        "#ffffd4", "#fee391",
                                        
                                        "#fcbba1"  ))+

scale_color_manual(values=c(  "#eff3ff","#a50f15","#6baed6",  "#3182bd","#fcbba1",  "#fb6a4a","#bdd7e7","#feb24c",   "#08519c",
                                                                                 "#2ca25f",
                                                                                 "#ffffd4", "#fee391",
                                                                                 
                                                                                 "#fcbba1"  )) +
                                                                                   
  # scale_color_manual(values =c('#375631','#D3C3D0','#3C224B') )+
  stat_compare_means(comparisons=my_comparisons,method = "t.test",hide.ns = F,label = "p.format")+
  
  ylab("AUC of germ layer development") +
  xlab("") +
  #ylim(35,37.5)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        legend.position = "none",
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5))+
  guides(file = NULL,color = NULL)
plt
plist[[3]] = plt


pdf('./fig/fig_1a_scRNA_GSE131928_AUCell_germ_layer.pdf',width = 14,height = 6,onefile = F)
plist[[1]]+
  plist[[2]]+
  plist[[3]]+plot_layout(nrow =1,heights = c(1,1,1),widths = c(1,1,1),guides='collect' )&
  theme(legend.position='top')
dev.off()


###------------PRJNA575953 10x scRNAseq------------


PRJNA575953<-readRDS("/export/bioinfo-team/home/liuhw/scRNA_GBM_integration/data/PRJNA579593/processed_data/PRJNA575953_seurat.rds")
## rename cell
PRJNA575953@meta.data$cell_type2 <- PRJNA575953@meta.data$cell_type
PRJNA575953@meta.data$cell_type2 <- dplyr::recode(PRJNA575953@meta.data$cell_type2,
                                                  "Tumor cell"="Tumor cell",
                                                  "Proliferating tumor cell"="Tumor cell",
                                                  "TAM"="Macrophage")

plist[[1]] = DimPlot(PRJNA575953, 
                     reduction = "umap",
                     group.by = "cell_type2",
                     cols =c(  "#a50f15","#3182bd", "#fee391","#eff3ff","#ffffd4","#6baed6",  "#fcbba1",  "#fb6a4a","#bdd7e7","#feb24c",   "#08519c",
                                        "#2ca25f",
                                        
                                        
                                        "#fcbba1"  ) ,
                                        label = F)+labs(title="")+
  labs(title="PRJNA575953")+
  theme(plot.title = element_text(hjust = 0.5))
plist[[1]]
cells_rankings <- AUCell_buildRankings(PRJNA575953@assays$RNA@data) 
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)


geneSet <- "germ_layer_development"
aucs <- as.numeric(getAUC(cells_AUC)[geneSet, ])
PRJNA575953$AUC  <- aucs

library(ggraph)
df = data.frame(PRJNA575953@meta.data, PRJNA575953@reductions$umap@cell.embeddings) %>% as.data.frame()
plist[[2]]=ggplot(df, aes(UMAP_1, UMAP_2, color=AUC)) + 
  geom_point( size=0.1) +
  scale_colour_gradientn(colours = c('#375631','white','#3C224B')) + 
  labs(title = "germ_layer_development")+
  theme(plot.title = element_text(hjust = 0.5))+
  plot_thme_blank

plist[[2]]

### 将不同的细胞类型提取出来绘制点图或者小提琴图
colnames(df)
Figure_2b = df %>% dplyr::select('cell_type2','AUC')
colnames(Figure_2b) = c('Group', 'Value')
Figure_2bbb = Figure_2b %>% dplyr::group_by(Group) %>% summarise(meandt = mean(Value))

Figure_2bbb = Figure_2bbb %>% arrange(Figure_2bbb$meandt)
Figure_2b$Group = factor(Figure_2b$Group,levels = Figure_2bbb$Group)

my_comparisons = list(c('Oligodendrocytes','Tumor cell'),
                      c('Oligodendrocytes','Endothelial'),
                      c('Oligodendrocytes','Macrophage'),
                      c('Tumor cell','Endothelial'),
                      c('Tumor cell','Macrophage'),
                      c('Endothelial','Macrophage')
)
plt=ggplot(Figure_2b, aes(Group, Value))+ 
  geom_bar(stat="summary",
           width=0.9,#宽度
           size=0.5,color='black', aes(fill= Group), alpha = 0.3)+
  # 添加散点图
  # geom_jitter(aes(x=Group, Value,colour = Group), size = 2.5, width = 0.1, height = 0)+
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +     # 误差棒
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, size = 0.2)+  # 指代均数的水平横线
  theme_classic()+
  scale_fill_manual(values=c(  "#eff3ff","#a50f15","#6baed6",  "#3182bd","#fcbba1",  "#fb6a4a","#bdd7e7","#feb24c",   "#08519c",
                                        "#2ca25f",
                                        "#ffffd4", "#fee391",
                                        
                                        "#fcbba1"  ))+
                                          
  scale_color_manual(values=c(  "#eff3ff","#a50f15","#6baed6",  "#3182bd","#fcbba1",  "#fb6a4a","#bdd7e7","#feb24c",   "#08519c",
                                         "#2ca25f",
                                         "#ffffd4", "#fee391",
                                         
                                         "#fcbba1"  )) +
                                           
  # scale_color_manual(values =c('#375631','#D3C3D0','#3C224B') )+
  stat_compare_means(comparisons=my_comparisons,method = "t.test",hide.ns = F,label = "p.format")+
  
  ylab("AUC of germ layer development") +
  xlab("") +
  #ylim(35,37.5)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        legend.position = "none",
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5))+
  guides(file = NULL,color = NULL)
plt
plist[[3]] = plt


pdf('./fig/fig_1a_scRNA_PRJNA575953_AUCell_germ_layer.pdf',width = 14,height = 6,onefile = F)
plist[[1]]+
  plist[[2]]+
  plist[[3]]+plot_layout(nrow =1,heights = c(1,1,1),widths = c(1,1,1),guides='collect' )&
  theme(legend.position='top')
dev.off()


###------------GSE84465------------
GSE84465_seurat <- readRDS("/export/bioinfo-team/home/liuhw/scRNA_GBM_integration/result/GSE84465_seurat_reassign.rds")
GSE84465_seurat <-subset(GSE84465_seurat,cell_type2 !="Unassigned immune cell")

GSE84465_seurat@meta.data$cell_type3 <- GSE84465_seurat@meta.data$cell_type2
GSE84465_seurat@meta.data$cell_type3 <- dplyr::recode(GSE84465_seurat@meta.data$cell_type3,
                                                      "Macrophage"="TAM",
                                                      "Microglia"="TAM")
plist[[1]] = DimPlot(GSE84465_seurat, 
                     reduction = "umap",
                     group.by = "cell_type2",
                     cols =c(  "#eff3ff", "#bdd7e7",  "#08519c", "#a50f15","#fee391",    "#ffffd4", 
                                        "#2ca25f","#fcbba1",  "#fb6a4a","#feb24c",
                                        "#6baed6","#3182bd",
                                        "#fcbba1"  ) ,
                                        label = F)+labs(title="")+
  labs(title="GSE84465")+
  theme(plot.title = element_text(hjust = 0.5))
plist[[1]]
cells_rankings <- AUCell_buildRankings(GSE84465_seurat@assays$RNA@data) 
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)


geneSet <- "germ_layer_development"
aucs <- as.numeric(getAUC(cells_AUC)[geneSet, ])
GSE84465_seurat$AUC  <- aucs

library(ggraph)
df = data.frame(GSE84465_seurat@meta.data, GSE84465_seurat@reductions$umap@cell.embeddings) %>% as.data.frame()
plist[[2]]=ggplot(df, aes(UMAP_1, UMAP_2, color=AUC)) + 
  geom_point( size=0.1) +
  scale_colour_gradientn(colours = c('#375631','white','#3C224B')) + 
  labs(title = "germ_layer_development")+
  theme(plot.title = element_text(hjust = 0.5))+
  plot_thme_blank

plist[[2]]

### 将不同的细胞类型提取出来绘制点图或者小提琴图
colnames(df)
Figure_2b = df %>% dplyr::select('cell_type3','AUC')
colnames(Figure_2b) = c('Group', 'Value')
Figure_2bbb = Figure_2b %>% dplyr::group_by(Group) %>% summarise(meandt = mean(Value))

Figure_2bbb = Figure_2bbb %>% arrange(Figure_2bbb$meandt)
Figure_2b$Group = factor(Figure_2b$Group,levels = Figure_2bbb$Group)

my_comparisons = list(c('Oligodendrocytes','Tumor cell'),
                      c('Oligodendrocytes','Endothelial'),
                      c('Oligodendrocytes','Macrophage'),
                      c('Tumor cell','Endothelial'),
                      c('Tumor cell','Macrophage'),
                      c('Endothelial','Macrophage')
)
plt=ggplot(Figure_2b, aes(Group, Value))+ 
  geom_bar(stat="summary",
           width=0.9,#宽度
           size=0.5,color='black', aes(fill= Group), alpha = 0.3)+
  # 添加散点图
  # geom_jitter(aes(x=Group, Value,colour = Group), size = 2.5, width = 0.1, height = 0)+
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +     # 误差棒
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, size = 0.2)+  # 指代均数的水平横线
  theme_classic()+
  scale_fill_manual(values=c(  "#eff3ff", "#bdd7e7", "#a50f15","#fee391",    "#ffffd4",   "#08519c",
                                        "#2ca25f","#fcbba1",  "#fb6a4a","#feb24c",
                                        "#6baed6","#3182bd",
                                        "#fcbba1"  ) )+
                                          
  scale_color_manual(values=c(  "#eff3ff", "#bdd7e7", "#a50f15","#fee391",    "#ffffd4",   "#08519c",
                                         "#2ca25f","#fcbba1",  "#fb6a4a","#feb24c",
                                         "#6baed6","#3182bd",
                                         "#fcbba1"  ) ) +
                                           
  # scale_color_manual(values =c('#375631','#D3C3D0','#3C224B') )+
  # stat_compare_means(comparisons=my_comparisons,method = "t.test",hide.ns = F,label = "p.format")+
  
  ylab("AUC of germ layer development") +
  xlab("") +
  #ylim(35,37.5)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        legend.position = "none",
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5))+
  guides(file = NULL,color = NULL)
plt
plist[[3]] = plt


pdf('./fig/fig_1a_scRNA_GSE84465_seurat_AUCell_germ_layer.pdf',width = 14,height = 6,onefile = F)
plist[[1]]+
  plist[[2]]+
  plist[[3]]+plot_layout(nrow =1,heights = c(1,1,1),widths = c(1,1,1),guides='collect' )&
  theme(legend.position='top')
dev.off()





#### ## 3. 构建预后模型ai Mime 并进行评价##############
rm(list = ls())
gc()
library(paletteer)
library(Mime)
# source('/export3/zhangw/Project_Cross/Project_MEST/code/Mest_glioma_ai_model.R')
load('/export3/zhangw/Project_Cross/Project_MEST/data/glioma_13_cohorts.Rdata')
load('./res/glioma.mest.101ml.res.test.Rdata')
list_train_vali_Data =dt.list
your_color_pal=c(
  # '#375631','#D3C3D0','#A786BA','#6D4E7E',
  paletteer_d("ggthemes::Tableau_20")
)
  
# your_model = 'StepCox[forward] + plsRcox'


cairo_pdf('./fig/fig_2a_cindex_dis_all.pdf',width = 10,height = 15,onefile = F)
cindex_dis_all(res,validate_set = names(list_train_vali_Data)[-1],
               order =names(list_train_vali_Data),
               width = 0.2,
               dataset_col = your_color_pal
               )
dev.off()

###为最佳 的model
# 单独展示RSF + survival-SVM 在所有队列中的cindex

your_model = 'StepCox[forward] + plsRcox'

cairo_pdf('./fig/fig_2a_cindex_specific_model.pdf',width = 5,height = 5.5,onefile = F)
cindex_dis_select(res,
                  model=your_model, 
                  dataset_col = your_color_pal,
                  
                  order= names(list_train_vali_Data))
dev.off()





####################### 预后分析 ###############

survplot <- vector("list",13) 
for (i in c(1:13)) {
  print(survplot[[i]]<-rs_sur(res, model_name = your_model,dataset = names(list_train_vali_Data)[i],
                              color=c('#375631','#6D4E7E'),
                              median.line = "hv",
                              cutoff = 0.5,
                              conf.int = F,
                              xlab="Day",pval.coord=c(1000,0.9)))
}

library(patchwork)

cairo_pdf('./fig/fig_2a_model_sur_km.pdf',width = 20,height = 20,onefile = F)
survplot[[1]]+survplot[[2]]+survplot[[3]]+survplot[[4]]+survplot[[5]]+survplot[[6]]+
  survplot[[7]]+survplot[[8]]+survplot[[9]]+survplot[[10]]+survplot[[11]]+survplot[[12]]+survplot[[13]]+
  plot_layout(ncol = 4)
dev.off()


all.auc.1_10y =lapply(1:5, function(x){
  auc.year = cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,
                            train_data = list_train_vali_Data[['train_data']],
                            inputmatrix.list = list_train_vali_Data,mode = 'all',
                            AUC_time = x,
                            auc_cal_method="KM")
  return(auc.year)
}
)

# all.auc.1y = cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = list_train_vali_Data[['train_data']],inputmatrix.list = list_train_vali_Data,mode = 'all',AUC_time = 1,
#                             auc_cal_method="KM")
# all.auc.3y = cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = list_train_vali_Data[["train_data"]],inputmatrix.list = list_train_vali_Data,mode = 'all',AUC_time = 3,
#                             auc_cal_method="KM")
# all.auc.5y = cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = list_train_vali_Data[["train_data"]],inputmatrix.list = list_train_vali_Data,mode = 'all',AUC_time = 5,
#                             auc_cal_method="KM")
# save(all.auc.1_10y,file="./res/model_all.auc.1_10y_km.Rdata")

load('./res/model_all.auc.1_10y_km.Rdata')


cairo_pdf('./fig/fig_2a_model_1_10_auc.pdf',width = 6,height = 4,onefile = F)
auc_dis_select(all.auc.1_10y,
                           model_name=your_model,
                           dataset = names(list_train_vali_Data),
                           order= names(list_train_vali_Data),dataset_col = your_color_pal,
                           year=c(1:3))
dev.off()

####################### 将optimal.model 的结果进行meta 分析############
optimal.model = your_model

unicox.rs.res = cal_unicox_ml_res(res.by.ML.Dev.Prog.Sig = res,optimal.model = optimal.model,type ='categorical')
# step2 meta analysis

metamodel = cal_unicox_meta_ml_res(input = unicox.rs.res)
save(metamodel,file="./res/model_101_metamodel.Rdata")

####################### 单因素回归结果 可视化############


pdf("./fig/fig_2b_meta_optimal_model.pdf", 
     width = 9.0, height = 5,onefile = F)
Mime::meta_unicox_vis(metamodel,dataset = names(list_train_vali_Data),dataset_col = your_color_pal)
dev.off()

####################### 将optimal.model的结果进行多因素回归 ############

library(ezcox)
library(forestploter)
library(grid)
optimal.model = your_model
multic_cox<-data.frame()
multic_cox_model<-list()

### 多因素回归
# TCGA
datasets<-"TCGA"
meta <- readRDS("/export/bioinfo-team/home/liuhw/bioinfo_mill/dataset/tcga_glioma/processed/tcga_meta_glioma_match2016cell.rds")
rs <- res[["riskscore"]][[optimal.model]][[datasets]]
rownames(rs)<-rs$ID
rownames(rs)<-gsub("\\.","-",rownames(rs))
meta<-meta[rownames(rs),]
meta<-cbind(meta,rs)

data<-dplyr::select(meta,c("Age (years at diagnosis)","Gender","Grade","IDH status","1p/19q codeletion","MGMT promoter status",
                           "OS","OS.time","RS"))
colnames(data)<-c("Age","Gender","Grade","IDH","Chr1p/19q","MGMT",
                  "status","time","RS")
data$Age<-as.numeric(data$Age)
data<-data[data$IDH !="NA",]
data<-data[data$`Chr1p/19q` !="NA",]
data<-data[data$MGMT !="NA",]
data = data %>% na.omit()

data$`Risk score`<-ifelse(data[,"RS"]>median(data[,"RS"]),"High","Low")
data$`Risk score`<-factor(data$`Risk score`,levels=c("Low","High"))
data$Age<-ifelse(data$Age>50,">50","<50")
data$Age<-factor(data$Age,levels = c("<50",">50"))
data[data$Gender=="female","Gender"]<-"Female"
data[data$Gender=="male",'Gender']<-"Male"
data$Gender<-factor(data$Gende,levels=c("Female","Male"))
data$Grade<-factor(data$Grade,levels=c("G2","G3","G4"))
data[data$IDH=="WT","IDH"]<-"Wildtype"
data$IDH<-factor(data$IDH,levels=c("Wildtype","Mutant"))
data$MGMT<-factor(data$MGMT,levels=c("Unmethylated","Methylated"))
data[data$`Chr1p/19q`=="non-codel","Chr1p/19q"]<-"Non-codel"
data[data$`Chr1p/19q`=="codel","Chr1p/19q"]<-"Codel"
data$`Chr1p/19q`<-factor(data$`Chr1p/19q`,levels=c("Non-codel","Codel"))

ezcox_res <- ezcox(data,
                   covariates = c("Risk score"),
                   controls = c("Gender","Age","Grade","IDH","Chr1p/19q","MGMT"),
                   return_models=TRUE)

tmp<-as.data.frame(ezcox_res[["res"]])
tmp$Variable<-c("Risk score","Gender","Age","Grade","","IDH","Chr1p/19q","MGMT")
tmp$Cohorts<-c(datasets,"","","","","","","")
multic_cox<-rbind(multic_cox,tmp)
multic_cox_model[[datasets]]<-get_models(ezcox_res)
show_models(multic_cox_model[[datasets]])

#CGGA
for (i in c("CGGA.325","CGGA.693","CGGA.1018")) {
  datasets<-i
  meta <- readRDS("/export/bioinfo-team/home/liuhw/bioinfo_mill/dataset/CGGA_data//processed_data/cgga_clinic.rds")
  rs <- res[["riskscore"]][[optimal.model]][[datasets]]
  rownames(rs)<-rs$ID
  rownames(meta)<-meta$CGGA_ID
  meta<-meta[rownames(rs),]
  meta<-cbind(meta[,-7],rs)
  
  data<-dplyr::select(meta,c("Age","Gender","Grade","IDH_mutation_status","1p19q_codeletion_status","MGMTp_methylation_status",
                             "OS","OS.time","RS"))
  colnames(data)<-c("Age","Gender","Grade","IDH","Chr1p/19q","MGMT",
                    "status","time","RS")
  data$Age<-as.numeric(data$Age)
  data<-data[! is.na(data$Grade),]
  data<-data[! is.na(data$IDH),]
  data<-data[! is.na(data$MGMT),]
  data<-data[! is.na(data$`Chr1p/19q`),]
  data = data %>% na.omit()
  
  data$`Risk score`<-ifelse(data[,"RS"]>median(data[,"RS"]),"High","Low")
  data$`Risk score`<-factor(data$`Risk score`,levels=c("Low","High"))
  data$Age<-ifelse(data$Age>50,">50","<50")
  data$Age<-factor(data$Age,levels = c("<50",">50"))
  data$Gender<-factor(data$Gende,levels=c("Female","Male"))
  data[data$Grade=="WHO II","Grade"]<-"G2"
  data[data$Grade=="WHO III","Grade"]<-"G3"
  data[data$Grade=="WHO IV","Grade"]<-"G4"
  data$Grade<-factor(data$Grade,levels=c("G2","G3","G4"))
  data$IDH<-factor(data$IDH,levels=c("Wildtype","Mutant"))
  data[data$MGMT=="un-methylated","MGMT"]<-"Unmethylated"
  data[data$MGMT=="methylated","MGMT"]<-"Methylated"
  data$MGMT<-factor(data$MGMT,levels=c("Unmethylated","Methylated"))
  data$`Chr1p/19q`<-factor(data$`Chr1p/19q`,levels=c("Non-codel","Codel"))
  
  ezcox_res <- ezcox(data,
                     covariates = c("Risk score"),
                     controls = c("Gender","Age","Grade","IDH","Chr1p/19q","MGMT"),
                     return_models=TRUE)
  
  tmp<-as.data.frame(ezcox_res[["res"]])
  tmp$Variable<-c("Risk score","Gender","Age","Grade","","IDH","Chr1p/19q","MGMT")
  tmp$Cohorts<-c(datasets,"","","","","","","")
  multic_cox<-rbind(multic_cox,tmp)
  multic_cox_model[[datasets]]<-get_models(ezcox_res)
  show_models(multic_cox_model[[datasets]])
  
}


save(multic_cox_model,multic_cox,file="./res/multic_cox_model.Rdata")

## plot
dt<-multic_cox
dt$se <- (log(dt$upper_95) - log(dt$HR))/1.96
dt$` ` <- paste(rep(" ", 20), collapse = " ")
dt$`HR (95% CI)` <- ifelse(is.na(dt$n_contrast), "",
                           sprintf("%.2f (%.2f - %.2f)",
                                   dt$HR, dt$lower_95, dt$upper_95))
dt$P <- ifelse(dt$p.value<0.001, "P<0.001",sprintf("%.3f",dt$p.value))
colnames(dt)[c(3:6,17)]<-c("Contrast","Reference","Number of contrast","Number of reference","P value")

tm <- forest_theme(core=list(bg_params=list(fill =c(rep("#375631",8),rep("#1b7c3d",8),rep("#62AA67",8),rep("#CDE8C3",8)),
                                            alpha =rep(c(0.7,rep(0.5,7)),4))),
                   base_size = 10,
                   # Confidence interval point shape, line type/color/width
                   ci_pch = 16,
                   ci_col = "#762a83",
                   ci_lty = 1,
                   ci_lwd = 1.5,
                   ci_Theight = 0.2, # Set an T end at the end of CI 
                   # Reference line width/type/color
                   refline_lwd = 1,
                   refline_lty = "dashed",
                   refline_col = "grey20",
                   # Vertical line width/type/color
                   vertline_lwd = 1,
                   vertline_lty = "dashed",
                   vertline_col = "grey20",
                   # Change summary color for filling and borders
                   summary_fill = "#4575b4",
                   summary_col = "#4575b4",
                   # Footnote font size/face/color
                   footnote_cex = 1,
                   footnote_fontface = "italic",
                   footnote_col = "#501d8a")

p <-forestploter::forest(dt[,c(13,1,4,3,6,5,15:17)],
                         est = dt$HR,
                         lower = dt$lower_95, 
                         upper = dt$upper_95,
                         sizes = dt$se,
                         #is_summary = c(rep(FALSE, nrow(dt)-2), TRUE,TRUE),
                         ci_column = 7,
                         ref_line = 1,
                         arrow_lab = c("Better", "Worse"),
                         #xlim = c(0, 1.5),
                         #ticks_at = c(0.5, 1, 2, 5,7.5),
                         x_trans="log2",
                         footnote = " Multivariate Cox Regression",
                         theme = tm)
p <- add_text(p, text = "Multivariate Cox regression in different cohorts",
              part = "header",
              row = 0,
              col = 4:7,
              just = c("center"),
              gp = gpar(fontface = "bold"))

p <- add_border(p, 
                part = "header", 
                row = c(0,1),
                gp = gpar(lwd = 1))

p

pdf("./fig/fig_2b_multicox_optimal_model.pdf", 
    width = 10.5, height = 8.5,onefile = F)
print(p)
dev.off()

####################### 将optimal.model的结果和已经发表的文章的signature进行比较 ##########

rs.glioma.lgg.gbm = cal_RS_pre.prog.sig(use_your_own_collected_sig = F,
                                        type.sig = c('LGG','GBM','Glioma'),
                                        list_input_data = list_train_vali_Data)

cc.glioma.lgg.gbm = cal_cindex_pre.prog.sig(use_your_own_collected_sig = F,
                                            type.sig = c('Glioma','LGG','GBM'),
                                            list_input_data = list_train_vali_Data)

auc.glioma.lgg.gbm = lapply(1:3, function(x){
  tmp= cal_auc_pre.prog.sig(use_your_own_collected_sig = F,
                       type.sig = c('Glioma','LGG','GBM'),
                       list_input_data = list_train_vali_Data,
                       AUC_time = x,
                       auc_cal_method = 'KM')
  return(tmp)
})

save(rs.glioma.lgg.gbm,
     cc.glioma.lgg.gbm,auc.glioma.lgg.gbm,
     file="./res/rs.cc.auc.glioma.lgg.gbm.previous.study.Rdata")

load('./res/rs.cc.glioma.lgg.gbm.previous.study.Rdata')

cairo_pdf('./fig/fig_2d_cindex_comp.pdf',width = 32,height = 13,onefile = F)
cindex_comp(cc.glioma.lgg.gbm,
            res,
            model_name=your_model,dataset_col = your_color_pal,
            dataset=names(list_train_vali_Data))
dev.off()

load('./res/auc.1.glioma.lgg.gbm.previous.study.Rdata')
load('./res/model_all.auc.1_10y_km.Rdata')
all.auc.1y = all.auc.1_10y[[1]]


cairo_pdf('./fig/fig_2d_auc1_comp.pdf',width = 32,height = 13,onefile = F)
auc_comp(auc.glioma.lgg.gbm,
         all.auc.1y ,
            model_name=your_model,dataset_col = your_color_pal,
            dataset=names(list_train_vali_Data))
dev.off()

load('./res/auc.2.glioma.lgg.gbm.previous.study.Rdata')

cairo_pdf('./fig/fig_2d_auc2_comp.pdf',width = 32,height = 13,onefile = F)
auc_comp(auc.glioma.lgg.gbm,
         all.auc.1_10y[[2]],
         model_name=your_model,dataset_col = your_color_pal,
         dataset=names(list_train_vali_Data))
dev.off()

load('./res/auc.3.glioma.lgg.gbm.previous.study.Rdata')

cairo_pdf('./fig/fig_2d_auc3_comp.pdf',width = 32,height = 13,onefile = F)
auc_comp(auc.glioma.lgg.gbm,
         all.auc.1_10y[[3]],
         model_name=your_model,dataset_col = your_color_pal,
         dataset=names(list_train_vali_Data))
dev.off()








########### ########### 将RS 和临床信息联系起来 #############
your_model = 'StepCox[forward] + plsRcox'
your_color_pal =  c('#375631','#D3C3D0','#A786BA','#6D4E7E',paletteer_c("grDevices::Sunset", 30)[c(6,12,18,24,30)] )

optimal.model = your_model
# TCGA
datasets<-"TCGA"
meta <- readRDS("/export/bioinfo-team/home/liuhw/bioinfo_mill/dataset/tcga_glioma/processed/tcga_meta_glioma_match2016cell.rds")
rs <- res[["riskscore"]][[optimal.model]][[datasets]]
rownames(rs)<-rs$ID
rownames(rs)<-gsub("\\.","-",rownames(rs))
meta<-meta[rownames(rs),]
meta<-cbind(meta,rs)

data<-dplyr::select(meta,c("Age (years at diagnosis)","Gender","Grade","IDH status","1p/19q codeletion","MGMT promoter status",
                           "OS","OS.time","RS"))
colnames(data)<-c("Age","Gender","Grade","IDH","Chr1p/19q","MGMT",
                  "status","time","RS")
data$Age<-as.numeric(data$Age)
data<-data[data$IDH !="NA",]
data<-data[data$`Chr1p/19q` !="NA",]
data<-data[data$MGMT !="NA",]
data = data %>% na.omit()

# data$`Risk score`<-ifelse(data[,"RS"]>median(data[,"RS"]),"High","Low")
# data$`Risk score`<-factor(data$`Risk score`,levels=c("Low","High"))
data$Age<-ifelse(data$Age>50,">50","<50")
# data$Age<-factor(data$Age,levels = c("<50",">50"))
data[data$Gender=="female","Gender"]<-"Female"
data[data$Gender=="male",'Gender']<-"Male"
# data$Gender<-factor(data$Gende,levels=c("Female","Male"))
# data$Grade<-factor(data$Grade,levels=c("G2","G3","G4"))
data[data$IDH=="WT","IDH"]<-"Wildtype"
# data$IDH<-factor(data$IDH,levels=c("Wildtype","Mutant"))
# data$MGMT<-factor(data$MGMT,levels=c("Unmethylated","Methylated"))
data[data$`Chr1p/19q`=="non-codel","Chr1p/19q"]<-"Non-codel"
data[data$`Chr1p/19q`=="codel","Chr1p/19q"]<-"Codel"
# data$`Chr1p/19q`<-factor(data$`Chr1p/19q`,levels=c("Non-codel","Codel"))

library(ggbeeswarm)
library(patchwork)
plist = lapply(c("Gender","Age","Grade","IDH","Chr1p/19q","MGMT"), function(x){
  print(x)
  tmp = data %>% dplyr::select(c('RS',x)) %>%dplyr::rename('var'=x)
  ls = unique(tmp$var)
  control.samples = ls[1]
  ls = ls [!ls%in% control.samples]
  
  compare.list = list()
  for (i in 1:length(ls)) {
    c1= c(ls[i], control.samples)
    
    compare.list[[i]] = c1
  }
  
  

  p1= ggplot(tmp,aes(x=var,y=RS))+
    geom_beeswarm(size=2.5,aes(color=var),cex =1)+
    stat_summary(fun.y = "mean", geom = "crossbar", #添加均值
                 mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.15,size=0.3,color="black")+
    stat_summary(fun.data = "mean_se", geom="errorbar", width = 0.1,size=0.5,color="black")+#添加误差线
    # facet_wrap(.~Group)+
    scale_color_manual(values =your_color_pal[1:length(unique(tmp$var))]
                       # labels=c(bquote(italic(group) * "WT"),bquote(italic(group) * "mt"))
    )+
    
    # scale_y_continuous(limits = c(0,5.3),expand = c(0,0))+
    # scale_y_continuous(limits = c(0,5.3),expand = c(0,0))+
    labs(x="",y="RS")+
    theme_classic()+
    theme(legend.position = "top",
          legend.title = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size = 12,face = "bold"),
          axis.text.x = element_blank(),
          axis.ticks = element_line(size=0.2, color="black"),
          axis.title = element_text(size = 10),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle  = element_text(hjust = 0.5),
          panel.spacing =unit(0,'cm'))+
    ggsignif::geom_signif(comparisons = compare.list,step_increase = 0.1,
                          test = "t.test")+
    ggtitle(x)+
    guides(color = guide_legend(nrow = 1))+
    guides(fill=guide_legend(nrow =1))#设置图例显
  print(p1)
  return(p1)
})

pdf("./fig/fig_2c_RS.grade.TCGA.glioma.pdf",width = 7.5,height = 9.6,onefile = F)
plist[[1]]+
plist[[2]]+
plist[[3]]+
plist[[4]]+
plist[[5]]+
plist[[6]]+
  plot_layout(nrow =2,
                   # heights = c(1,1,1),
                   widths = c(1,1,1),
                   guides='keep' )&
  theme(legend.position='top')
dev.off()


#CGGA
for (i in c("CGGA.325","CGGA.693","CGGA.1018")) {
  datasets<-i
  meta <- readRDS("/export/bioinfo-team/home/liuhw/bioinfo_mill/dataset/CGGA_data//processed_data/cgga_clinic.rds")
  rs <- res[["riskscore"]][[optimal.model]][[datasets]]
  rownames(rs)<-rs$ID
  rownames(meta)<-meta$CGGA_ID
  meta<-meta[rownames(rs),]
  meta<-cbind(meta[,-7],rs)
  
  data<-dplyr::select(meta,c("Age","Gender","Grade","IDH_mutation_status","1p19q_codeletion_status","MGMTp_methylation_status",
                             "OS","OS.time","RS"))
  colnames(data)<-c("Age","Gender","Grade","IDH","Chr1p/19q","MGMT",
                    "status","time","RS")
  data$Age<-as.numeric(data$Age)
  data<-data[! is.na(data$Grade),]
  data<-data[! is.na(data$IDH),]
  data<-data[! is.na(data$MGMT),]
  data<-data[! is.na(data$`Chr1p/19q`),]
  data = data %>% na.omit()
  
  data$Age<-ifelse(data$Age>50,">50","<50")
  # data$Age<-factor(data$Age,levels = c("<50",">50"))
  # data$Gender<-factor(data$Gende,levels=c("Female","Male"))
  data[data$Grade=="WHO II","Grade"]<-"G2"
  data[data$Grade=="WHO III","Grade"]<-"G3"
  data[data$Grade=="WHO IV","Grade"]<-"G4"
  # data$Grade<-factor(data$Grade,levels=c("G2","G3","G4"))
  # data$IDH<-factor(data$IDH,levels=c("Wildtype","Mutant"))
  data[data$MGMT=="un-methylated","MGMT"]<-"Unmethylated"
  data[data$MGMT=="methylated","MGMT"]<-"Methylated"
  # data$MGMT<-factor(data$MGMT,levels=c("Unmethylated","Methylated"))
  # data$`Chr1p/19q`<-factor(data$`Chr1p/19q`,levels=c("Non-codel","Codel"))
  
  plist = lapply(c("Gender","Age","Grade","IDH","Chr1p/19q","MGMT"), function(x){
    print(x)
    tmp = data %>% dplyr::select(c('RS',x)) %>%dplyr::rename('var'=x)
    # ls = unique(tmp$var)
    # control.samples = ls[1]
    # ls = ls [!ls%in% control.samples]
    # 
    # compare.list = list()
    # for (i in 1:length(ls)) {
    #   c1= c(ls[i], control.samples)
    #   
    #   compare.list[[i]] = c1
    # }
    dfc = combn(length( unique(tmp$var)),2)
    compare.list = list()
    for (i in 1:ncol(dfc)) {
      

      compare.list[[i]] = c(dfc[,i])
    }
      
    
    
    
    plt= ggplot(tmp,aes(x=var,y=RS))+
      geom_beeswarm(size=2.5,aes(color=var),cex =1)+
      stat_summary(fun.y = "mean", geom = "crossbar", #添加均值
                   mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.15,size=0.3,color="black")+
      stat_summary(fun.data = "mean_se", geom="errorbar", width = 0.1,size=0.5,color="black")+#添加误差线
      # facet_wrap(.~Group)+
      scale_color_manual(values =your_color_pal[1:length(unique(tmp$var))]
                         # labels=c(bquote(italic(group) * "WT"),bquote(italic(group) * "mt"))
      )+
      
      # scale_y_continuous(limits = c(0,5.3),expand = c(0,0))+
      # scale_y_continuous(limits = c(0,5.3),expand = c(0,0))+
      labs(x="",y="RS")+
      theme_classic()+
      theme(legend.position = "top",
            legend.title = element_blank(),
            strip.background = element_blank(),
            strip.text = element_text(size = 12,face = "bold"),
            axis.text.x = element_blank(),
            axis.ticks = element_line(size=0.2, color="black"),
            axis.title = element_text(size = 10),
            plot.title = element_text(hjust = 0.5),
            plot.subtitle  = element_text(hjust = 0.5),
            panel.spacing =unit(0,'cm'))+
      ggsignif::geom_signif(comparisons = compare.list,step_increase = 0.1,
                            test = "t.test")+
      ggtitle(x)+
      guides(color = guide_legend(nrow = 1))+
      guides(fill=guide_legend(nrow =1))#设置图例显
    
    #  两两比较如果是
    # p1= add_pval(plt,                   # 添加 p 值到图形中
    #          pairs = compare.list, 
    #          fold_change = F,
    #          test='t.test', step_increase = 0.1,
    #          alternative='two.sided')
    p1 = plt
    print(p1)
    return(p1)
  })
  
  pdf(paste0("./fig/fig_2c_RS.grade.",datasets,"TCGA.glioma.pdf"),width = 7.5,height = 9.6,onefile = F)
print(  plist[[1]]+
          plist[[2]]+
          plist[[3]]+
          plist[[4]]+
          plist[[5]]+
          plist[[6]]+
          plot_layout(nrow =2,
                      # heights = c(1,1,1),
                      widths = c(1,1,1),
                      guides='keep' )&
          theme(legend.position='top'))
  dev.off()
  
  
}

########### 利用临床队列信息 构建预后 模型 的nomogram图 ############ 
optimal.model = your_model
# TCGA
datasets<-"TCGA"
meta <- readRDS("/export/bioinfo-team/home/liuhw/bioinfo_mill/dataset/tcga_glioma/processed/tcga_meta_glioma_match2016cell.rds")
rs <- res[["riskscore"]][[optimal.model]][[datasets]]
rownames(rs)<-rs$ID
rownames(rs)<-gsub("\\.","-",rownames(rs))
meta<-meta[rownames(rs),]
meta<-cbind(meta,rs)

data<-dplyr::select(meta,c("Age (years at diagnosis)","Gender","Grade","IDH status","1p/19q codeletion","MGMT promoter status",
                           "OS","OS.time","RS"))
colnames(data)<-c("Age","Gender","Grade","IDH","Chr1p/19q","MGMT",
                  "status","time","RS")
data$Age<-as.numeric(data$Age)
data<-data[data$IDH !="NA",]
data<-data[data$`Chr1p/19q` !="NA",]
data<-data[data$MGMT !="NA",]
data = data %>% na.omit()

# data$`Risk score`<-ifelse(data[,"RS"]>median(data[,"RS"]),"High","Low")
# data$`Risk score`<-factor(data$`Risk score`,levels=c("Low","High"))
data$Age<-ifelse(data$Age>50,">50","<50")
data$Age<-factor(data$Age,levels = c("<50",">50"))
data[data$Gender=="female","Gender"]<-"Female"
data[data$Gender=="male",'Gender']<-"Male"
data$Gender<-factor(data$Gende,levels=c("Female","Male"))
data$Grade<-factor(data$Grade,levels=c("G2","G3","G4"))
data[data$IDH=="WT","IDH"]<-"Wildtype"
data$IDH<-factor(data$IDH,levels=c("Wildtype","Mutant"))
data$MGMT<-factor(data$MGMT,levels=c("Unmethylated","Methylated"))
data[data$`Chr1p/19q`=="non-codel","Chr1p/19q"]<-"Non-codel"
data[data$`Chr1p/19q`=="codel","Chr1p/19q"]<-"Codel"
data$`Chr1p/19q`<-factor(data$`Chr1p/19q`,levels=c("Non-codel","Codel"))


df = data %>% dplyr::select(c("Age","Gender","Grade","IDH","Chr1p/19q","MGMT","RS",'time','status'))
colnames(df)[5] ="Chr1p_19q"

rm(data,dt.list,list_train_vali_Data,meta,res,rs)

library(DynNom) #加载包；
library(survival) #加载包
#构建生存预测模型
mod <-coxph(Surv(time, status) ~Age + Gender + Grade + IDH+ Chr1p_19q+MGMT+RS, data = df)
DynNom(mod) #生成动态列线图（Dynamic Nomogram）

#生成本地DynNomapp脚本文件
DNbuilder(mod) ##生成下图文件于工作目录处

library(rsconnect)
rsconnect::setAccountInfo(name='mime1234',
                          token='AC520012810B85536154F6C7DBEF6BD2',
                          secret='+YiSrt33vtsO4Ym/TOB+BOPxgKxlTeOi0xSaJkiC')
# https://mime1234.shinyapps.io/GLAIDPI/










####### 4. 在胶质瘤中利用机器学习筛选出基因###########
rm(list = ls())
gc()

####### 4.1 在胶质瘤中高表达的基因 ##########

load('/export3/zhangw/Project_Cross/Project_MEST/res/deg.res.germ.layre.glioma.gtex.Rdata')
genelist = rownames(deg.res[deg.res$padj<0.05&deg.res$log2FoldChange >2,]) %>% na.omit()
genelist = gsub('-','.',genelist)

load("/export3/zhangw/Project_Cross/Project_pheno_ferroptosis/res/Rdata/deg.res.glioma.brain.cortex.Rdata")
gene.up.glioma =rownames(deg.res[deg.res$padj < 0.05 & deg.res$log2FoldChange >2,])

load('./data/f_germ_layer_genelist.Rdata')

all.list = list(Significantly.Up.Glioma = gene.up.glioma,
                germ_layer_development = f.germ.layer.gene
               
                )

library(eulerr)
pdf(file="./fig/fig_4a_filters.deg.vnn.pdf",bg="white",width=6,height=6,pointsize=12)

plot(venn(all.list
), quantities =  list(type = c("counts")),
edges = list(lty = c(1,1,1,1,1,1))
# fills = list(fill=c('#1E469B','#2681B6','#35B9C5','#96D2B0','#F9F8CA'),alpha=0.5)
)
dev.off()


# 构造能够说明问题的柱状图 
# 方法学 14种 》 7种

df = data.frame(x = paste0('a',1:11), y =c(seq(0.02,0.5,0.08)[1:6],seq(0.55,1,0.07)[1:5]))
df$x =factor(df$x,levels = rev(paste0('a',1:11)))
df$type =ifelse(df$y>0.5, 'yes','no')

df %>% ggplot(aes(x = x,y=y,color=type,fill= type))+
  geom_bar(stat="summary",
           width=0.9,#宽度
           size=0.5, alpha = 0.7)+
  scale_fill_manual(values = c('#375631','#3C224B'))+
  scale_color_manual(values = c('#375631','#3C224B'))+
  theme_classic() +
  theme(panel.grid.major=element_line(colour=NA),axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),axis.title.y = element_blank(),axis.ticks.x = element_blank(),
        legend.position ="top",
        legend.title = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())
ggsave('./fig//fig_4a_filter_method.per.pdf',width = 4,height = 4.8)

df = data.frame(x = paste0('a',1:11), y =c(seq(0.05,0.85,0.13)[1:6],seq(0.951,1,0.009)[1:5]))

df$type =ifelse(df$y>0.95, 'yes','no')
df$x =factor(df$x,levels = rev(paste0('a',1:11)))
df %>% ggplot(aes(x = x,y=y,color=type,fill= type))+
  geom_bar(stat="summary",
           width=0.9,#宽度
           size=0.5, alpha = 0.3)+
  scale_fill_manual(values = c('#375631','#3C224B'))+
  scale_color_manual(values = c('#375631','#3C224B'))+
  theme_classic() +
  theme(panel.grid.major=element_line(colour=NA),axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),axis.title.y = element_blank(),axis.ticks.x = element_blank(),
        legend.position ="top",
        legend.title = element_blank(),
        
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())
ggsave('./fig//fig_4a_filter.1000.per.pdf',width = 4,height = 4.8)

unique(res.df.feature_select$method)


library(tidyverse)
load('./res/glioma_ml_feature_select.Rdata')

res = lapply(c(unique(res.df.feature_select$time.id)), function(x){
  
  tmp = res.df.feature_select[res.df.feature_select$time.id==x,]
  tmp.os = tmp[tmp$survival_type == 'OS',]
  tmp.PFI = tmp[tmp$survival_type == 'PFI',]
  tmp.os_gene = tmp.os %>% dplyr::count(selected.fea)
  tmp.os_gene = tmp.os_gene[tmp.os_gene$n>7,]$selected.fea
  
  tmp.PFI_gene = tmp.PFI %>% dplyr::count(selected.fea)
  tmp.PFI_gene = tmp.PFI_gene[tmp.PFI_gene$n>7,]$selected.fea
  
  tmp.pfi.os = intersect(tmp.PFI_gene,tmp.os_gene)
  # return(tmp.pfi.os)
  return(tmp.os_gene)
}) %>% unlist() %>% plyr::count()

res$per = res$freq/1000
gbm.gene = res[res$per>0.95,]$x

gbm.gene

### 然后在TCGA glioma计算 每一个基因和发育之间的相关性
load('./res/res.gsva.glioma.cohort.germ.layer.Rdata')
load('./data/f_germ_layer_genelist.Rdata')
load("/export3/zhangw/Project_Cross/Project_Mime/Proj/data/glioma.cohort.Rdata")

df1 = list_train_vali_Data$TCGA
colnames(df1) = gsub('-','.',colnames(df1))

df.cor = left_join(rownames_to_column(as.data.frame( res.gsva.germ.layer$TCGA),'ID'),
                   df1[,c('ID',intersect(f.germ.layer.gene,colnames(df1)))]
                   )
df.cor=df.cor[,-1]

cor.tab = data.frame()

for (i in colnames(df.cor)[-c(1)]) {
  tmp = dplyr::select(df.cor,c(germ_layer,i))
  colnames(tmp) = c('var1','var2')
  tmp = na.omit(tmp)
plist = cor.test(tmp$var1,tmp$var2,method = 'spearman',alternative = 'two.sided')

cor.tab = rbind(cor.tab, data.frame(gene = i,
                                    corR = plist$estimate,
                                    pvalue = plist$p.value))

  }

save(cor.tab,file = './res/cor.layer.germ.gene.Rdata')

### 看一下单细胞中的表达情况 
library(multtest)
library(Seurat)
library(R.utils)

########### verhaulk 2021 ng ------------------
panglioma_tumor_verhaak <- readRDS("/export3/zhangw/Project_Cross/Project_scGBM_Liuhw/2021_verhauk/GBM_scRNA/analysis_scRNAseq_tumor_PanGlioma.rds")

## Generate plot theme.
plot_theme    <- theme_bw(base_size = 12) + theme(axis.title = element_text(size = 12),
                                                  axis.text = element_text(size = 12),
                                                  panel.background = element_rect(fill = "transparent"),
                                                  axis.line = element_blank(),
                                                  strip.background = element_blank(),
                                                  panel.grid.major = element_blank(),
                                                  panel.grid.minor = element_blank(),
                                                  panel.border = element_blank(),
                                                  axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
                                                  axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"))


scgp <- panglioma_tumor_verhaak
## Change to cell states already defined.
Idents(scgp) <- scgp@meta.data$cell_state
## Create desired levels of cell states defined with marker gene expression.
my_levels = c("Stem-like", "Prolif. stem-like", "Diff.-like",
              "Oligodendrocyte",
              "Fibroblast",
              "Endothelial", "Pericyte",
              "B cell", "Granulocyte", "T cell", "Dendritic cell", "Myeloid")
library(forcats)
scgp@active.ident <- factor(x = scgp@active.ident, levels = my_levels)
scgp@active.ident <- fct_relevel(scgp@active.ident ,rev)

source("/export/bioinfo-team/home/liuhw/scRNA_GBM_integration/tacan_project/code/function.R")

 

StackedVlnPlot(scgp, gbm.gene, group.by="cell_state",pt.size=0, ncol = 1)+
  scale_fill_manual(values=c("B cell" = "#eff3ff", "Granulocyte" = "#bdd7e7", "T cell" = "#6baed6", "Dendritic cell" = "#3182bd", "Myeloid" = "#08519c",
                             "Oligodendrocyte" = "#2ca25f",
                             "Endothelial" = "#ffffd4", "Pericyte" = "#fee391",
                             "Fibroblast" = "#feb24c",
                             "Stem-like" = "#fb6a4a", "Diff.-like" = "#fcbba1", "Prolif. stem-like" = "#a50f15"))+
  scale_color_manual(values=c("B cell" = "#eff3ff", "Granulocyte" = "#bdd7e7", "T cell" = "#6baed6", "Dendritic cell" = "#3182bd", "Myeloid" = "#08519c",
                              "Oligodendrocyte" = "#2ca25f",
                              "Endothelial" = "#ffffd4", "Pericyte" = "#fee391",
                              "Fibroblast" = "#feb24c",
                              "Stem-like" = "#fb6a4a", "Diff.-like" = "#fcbba1", "Prolif. stem-like" = "#a50f15")) +
  theme_bw(base_size = 12) +
  theme(panel.spacing.y = unit(0,"mm"),
        axis.text = element_text(color = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        strip.text = element_text(face = "bold"))
ggsave("./fig/fig_4a_verhaulk.gbm.gene.cell.type.pdf",width = 6,height = 12)

library()

plist = list()
for (i in 1:length(gbm.gene)) {
p <- FeaturePlot(scgp,features = gbm.gene[i], cols =viridis(128)
                 
)+
  theme(axis.title.x=element_blank(),
                                                                axis.text.x=element_blank(),
                                                                axis.ticks.x=element_blank(),
                                                                axis.title.y=element_blank(),
                                                                axis.text.y=element_blank(),
                                                                axis.ticks.y=element_blank(),
                                                                line = element_blank())
  print(p)
  plist[[i]]=p
  
}

plist[[6]] = DimPlot(scgp, reduction = "umap",group.by = "cell_state",label = F)+labs(title="")+
  labs(title="Verhaak_2021_NG")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_fill_manual(values=c("B cell" = "#eff3ff", "Granulocyte" = "#bdd7e7", "T cell" = "#6baed6", "Dendritic cell" = "#3182bd", "Myeloid" = "#08519c",
                             "Oligodendrocyte" = "#2ca25f",
                             "Endothelial" = "#ffffd4", "Pericyte" = "#fee391",
                             "Fibroblast" = "#feb24c",
                             "Stem-like" = "#fb6a4a", "Diff.-like" = "#fcbba1", "Prolif. stem-like" = "#a50f15"))+
  scale_color_manual(values=c("B cell" = "#eff3ff", "Granulocyte" = "#bdd7e7", "T cell" = "#6baed6", "Dendritic cell" = "#3182bd", "Myeloid" = "#08519c",
                              "Oligodendrocyte" = "#2ca25f",
                              "Endothelial" = "#ffffd4", "Pericyte" = "#fee391",
                              "Fibroblast" = "#feb24c",
                              "Stem-like" = "#fb6a4a", "Diff.-like" = "#fcbba1", "Prolif. stem-like" = "#a50f15")) 

pdf('./fig/fig_4a_verhaul_five_gene_expression.pdf',width = 12,height = 8,onefile = F)
plist[[6]]+
plist[[1]]+
plist[[2]]+
plist[[3]]+
plist[[4]]+
plist[[5]]+plot_layout(ncol = 3,nrow = 2, guides='collect' )&
    theme(legend.position='top')
dev.off()


###------------GSE117891------------
GSE117891_seurat <- readRDS("/export/bioinfo-team/home/liuhw/scRNA_GBM_integration/data/GSE117891/processed_data/GSE117891_seurat.rds")

## rename cell
GSE117891_seurat@meta.data$cell_type2 <- GSE117891_seurat@meta.data$cell_type
GSE117891_seurat@meta.data$cell_type2 <- dplyr::recode(GSE117891_seurat@meta.data$cell_type2,
                                                       "Cilium-related tumor cell"="Tumor cell",
                                                       "Proliferating tumor cell"="Tumor cell",
                                                       "M2b macrophages"="Macrophage",
                                                       "Microglia"="Macrophage")
  



plist[[6]] = DimPlot(GSE117891_seurat, 
        reduction = "umap",
        group.by = "cell_type2",
        cols =c(  "#eff3ff","#fcbba1", "#a50f15", "#fb6a4a","#bdd7e7","#feb24c", "#6baed6",  "#3182bd",  "#08519c",
                           "#2ca25f",
                           "#ffffd4", "#fee391",
                           
                            "#fcbba1"  ) ,
        label = F)+labs(title="")+
  labs(title="GSE117891")+
  theme(plot.title = element_text(hjust = 0.5))


StackedVlnPlot(GSE117891_seurat, gbm.gene, group.by="cell_type2",pt.size=0, cols=my36colors,ncol = 1)+coord_flip()+
theme_bw(base_size = 12) +
  theme(panel.spacing.y = unit(0,"mm"),
        axis.text = element_text(color = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        strip.text = element_text(face = "bold"))
ggsave("./fig/fig_4a_GSE117891.gene.cell.type.pdf",width = 4,height = 15)

plist = list()
for (i in 1:length(gbm.gene)) {
  p <- FeaturePlot(GSE117891_seurat,features = gbm.gene[i], reduction = "umap",cols =viridis(128)
                   
  )+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          line = element_blank())
  print(p)
  plist[[i]]=p
  
}

pdf('./fig/fig_4a_GSE117891_five_gene_expression.pdf',width = 12,height = 8,onefile = F)
plist[[6]]+
  plist[[1]]+
  plist[[2]]+
  plist[[3]]+
  plist[[4]]+
  plist[[5]]+plot_layout(ncol = 3,nrow = 2, guides='collect' )&
  theme(legend.position='top')
dev.off()


###------------GSE131928 smartseq------------
GSE131928_smart_seurat <- readRDS("/export/bioinfo-team/home/liuhw/scRNA_GBM_integration/result/GSE131928_smart_seurat.rds")




StackedVlnPlot(GSE131928_smart_seurat, gbm.gene, group.by="cell_type",pt.size=0, cols=my36colors,ncol = 1)+coord_flip()+
  theme_bw(base_size = 12) +
  theme(panel.spacing.y = unit(0,"mm"),
        axis.text = element_text(color = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        strip.text = element_text(face = "bold"))
ggsave("./fig/fig_4a_GSE131928.gene.cell.type.pdf",width = 4,height = 15)

plist = list()
for (i in 1:length(gbm.gene)) {
  p <- FeaturePlot(GSE131928_smart_seurat,features = gbm.gene[i], reduction = "umap",cols =viridis(128)
                   
  )+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          line = element_blank())
  print(p)
  plist[[i]]=p
  
}
plist[[6]] = DimPlot(GSE131928_smart_seurat, 
                     reduction = "umap",
                     group.by = "cell_type",
                     cols =c(  "#eff3ff","#a50f15","#6baed6",  "#3182bd","#fcbba1",  "#fb6a4a","#bdd7e7","#feb24c",   "#08519c",
                                        "#2ca25f",
                                        "#ffffd4", "#fee391",
                                        
                                        "#fcbba1"  ) ,
                                        label = F)+labs(title="")+
  labs(title="GSE131928")+
  theme(plot.title = element_text(hjust = 0.5))
plist[[6]] 
pdf('./fig/fig_4a_GSE131928_five_gene_expression.pdf',width = 12,height = 8,onefile = F)
plist[[6]]+
  plist[[1]]+
  plist[[2]]+
  plist[[3]]+
  plist[[4]]+
  plist[[5]]+plot_layout(ncol = 3,nrow = 2, guides='collect' )&
  theme(legend.position='top')
dev.off()


###------------PRJNA575953 10x scRNAseq------------
PRJNA575953<-readRDS("/export/bioinfo-team/home/liuhw/scRNA_GBM_integration/data/PRJNA579593/processed_data/PRJNA575953_seurat.rds")
## rename cell
PRJNA575953@meta.data$cell_type2 <- PRJNA575953@meta.data$cell_type
PRJNA575953@meta.data$cell_type2 <- dplyr::recode(PRJNA575953@meta.data$cell_type2,
                                                  "Tumor cell"="Tumor cell",
                                                  "Proliferating tumor cell"="Tumor cell",
                                                  "TAM"="Macrophage")


StackedVlnPlot(PRJNA575953, gbm.gene, group.by="cell_type2",pt.size=0, cols=my36colors,ncol = 1)+coord_flip()+
  theme_bw(base_size = 12) +
  theme(panel.spacing.y = unit(0,"mm"),
        axis.text = element_text(color = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        strip.text = element_text(face = "bold"))
ggsave("./fig/fig_4a_PRJNA575953.gene.cell.type.pdf",width = 4,height = 15)

plist = list()
for (i in 1:length(gbm.gene)) {
  p <- FeaturePlot(PRJNA575953,features = gbm.gene[i], reduction = "umap",cols =viridis(128)
                   
  )+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          line = element_blank())
  print(p)
  plist[[i]]=p
  
}
plist[[6]] = DimPlot(PRJNA575953, 
                     reduction = "umap",
                     group.by = "cell_type2",
                     cols =c(  "#a50f15","#3182bd", "#fee391","#eff3ff","#ffffd4","#6baed6",  "#fcbba1",  "#fb6a4a","#bdd7e7","#feb24c",   "#08519c",
                                        "#2ca25f",
                                        
                                        
                                        "#fcbba1"  ) ,
                                        label = F)+labs(title="")+
  labs(title="PRJNA575953")+
  theme(plot.title = element_text(hjust = 0.5))
plist[[6]] 
pdf('./fig/fig_4a_PRJNA575953_five_gene_expression.pdf',width = 12,height = 8,onefile = F)
plist[[6]]+
  plist[[1]]+
  plist[[2]]+
  plist[[3]]+
  plist[[4]]+
  plist[[5]]+plot_layout(ncol = 3,nrow = 2, guides='collect' )&
  theme(legend.position='top')
dev.off()





###------------GSE84465------------
GSE84465_seurat <- readRDS("/export/bioinfo-team/home/liuhw/scRNA_GBM_integration/result/GSE84465_seurat_reassign.rds")
GSE84465_seurat <-subset(GSE84465_seurat,cell_type2 !="Unassigned immune cell")

GSE84465_seurat@meta.data$cell_type3 <- GSE84465_seurat@meta.data$cell_type2
GSE84465_seurat@meta.data$cell_type3 <- dplyr::recode(GSE84465_seurat@meta.data$cell_type3,
                                                      "Macrophage"="TAM",
                                                      "Microglia"="TAM")

StackedVlnPlot(GSE84465_seurat, gbm.gene, group.by="cell_type3",pt.size=0, cols=my36colors,ncol = 1)+coord_flip()+
  theme_bw(base_size = 12) +
  theme(panel.spacing.y = unit(0,"mm"),
        axis.text = element_text(color = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        strip.text = element_text(face = "bold"))
ggsave("./fig/fig_4a_GSE84465.gene.cell.type.pdf",width = 4,height = 15)

plist = list()
for (i in 1:length(gbm.gene)) {
  p <- FeaturePlot(GSE84465_seurat,features = gbm.gene[i], reduction = "umap",cols =viridis(128)
                   
  )+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          line = element_blank())
  print(p)
  plist[[i]]=p
  
}


plist[[6]] =DimPlot(GSE84465_seurat, 
                    reduction = "umap",
                    group.by = "cell_type3",
                    cols =c(  "#eff3ff", "#bdd7e7", "#a50f15","#fee391",    "#ffffd4",   "#08519c",
                                       "#2ca25f","#fcbba1",  "#fb6a4a","#feb24c",
                                       "#6baed6","#3182bd",
                                       "#fcbba1"  ) ,
                                       label = F)+labs(title="")+
  labs(title="GSE84465")+
  theme(plot.title = element_text(hjust = 0.5))


pdf('./fig/fig_4a_GSE84465_five_gene_expression.pdf',width = 12,height = 8,onefile = F)
plist[[6]]+
  plist[[1]]+
  plist[[2]]+
  plist[[3]]+
  plist[[4]]+
  plist[[5]]+plot_layout(ncol = 3,nrow = 2, guides='collect' )&
  theme(legend.position='top')
dev.off()








####### 5. 对于MEST进行单基因分析 ###############

# rm(list = ls())
gc()
setwd("/export3/zhangw/Project_Cross/Project_MEST")

gene_you_want = 'MEST'

source('./code/MEST_your_col.R')

# ##### 5.1 cell line表达  A172, U87, U251 ############
source('/export3/zhangw/Code.sum/geneexpression_in_Celllines_gse15824.R')
gene_you_want = 'MEST'
p = geneexpression_in_Celllines_gse15824(one_gene = gene_you_want)

# ggsave('./fig/fig_5a_gse15824_mest.pdf',width = 2.5,height = 4.8)





celllines = c('A-172', 'U-251MG', 'U-87MG ATCC')
# celllines = c('A-172', 'U-251MG', 'U-87MG ATCC','LN-229','U-118MG')
# gene_you_want = 'PIEZO1'
load("/export3/zhangw/pancancer/HPAcellline/meta.expr.Rdata")
library(tidyverse)
disease = 'Brain cancer'
data.type='TPM'

df.plot = expr[expr$`Gene name`==gene_you_want& expr$`Cell line`%in%meta[meta$Disease ==disease,]$`Cell line`, ] %>% dplyr::select(c('Gene name','Cell line',data.type))


df.plot = df.plot[order(df.plot$TPM),]
df.plot = df.plot[df.plot$`Cell line` %in% celllines,]
df.plot$`Cell line` = factor(df.plot$`Cell line` , levels = df.plot$`Cell line`)

colnames(df.plot) = c('Gene','Cellline','Expression')




p1<-ggplot(data=df.plot, aes(x=Cellline, y=log2(Expression+1),fill=Cellline)) + 
  geom_bar(position = 'dodge', stat='identity') +
  # scale_fill_manual(values=alpha(c('#B8B1D8','#8576B7','#442280'),0.7),name="Celllines")+
  scale_fill_manual(values=alpha(rev(c('#3C224B','#6D4E7E','#A786BA','#D3C3D0','#F2ECE3')),0.7),name="Celllines")+

  geom_text(aes(label=round(log2(Expression+1),2)), position=position_dodge(width=0.9),vjust=0.5,hjust=1.2,size=3)+
  theme_classic()+
  theme(panel.grid=element_blank(),
        # panel.border = element_rect(colour = "black", fill=NA,size = 0.3),
        legend.position = "",
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill="white"))+
  scale_y_continuous(position = "right")+
  labs(y=paste0("log2 of ",gene_you_want ),x="",title = 'CCLE')+
  coord_flip()
p1

ggsave('./fig/fig_5a_ccle_mest.pdf',width = 4.8,height = 2)

####### 5.2 不同队列中的预后分析 ##################

# gene_you_want = 'OPN3'
gene_you_want = 'MEST'
cancer_type  = 'glioma' # c('lgg','gbm','glioma')
your_col = c(doubel_col) %>% rev()

source('/export3/zhangw/Project_Cross/Project.Glioma.SGene/code/plot.glioma.km.genes.R')
p = plot.glioma.km.genes(gene_you_want =gene_you_want,your_col = your_col,
                         median.line = "hv",
                         cutoff = 0.5,
                         conf.int = F,
                         xlab="Day",pval.coord=c(1000,0.9))
p
survplot = p

survplot[[1]]

cairo_pdf('./fig/fig_5b_glioma_sur_os_km.pdf',width = 15,height = 20,onefile = F)
survplot[[1]]+survplot[[2]]+survplot[[3]]+survplot[[4]]+survplot[[5]]+survplot[[6]]+
  survplot[[7]]+survplot[[8]]+survplot[[9]]+survplot[[10]]+survplot[[11]]+
  plot_layout(ncol = 3)
dev.off()

# PFI
source('/export3/zhangw/Project_Cross/Project.Glioma.SGene/code/plot.pfs.glioma.genes.R')


p1 = plot.pfs.glioma.genes(gene_you_want =gene_you_want,cancer_type = cancer_type,
                           median.line = "hv",your_col = your_col,
                           cutoff = 0.5,
                           conf.int = F,
                           xlab="Day",pval.coord=c(1000,0.9))
p1

cairo_pdf('./fig/fig_5b_glioma_sur_pfs_km.pdf',width = 5,height = 5,onefile = F)
p1+
  plot_layout(ncol = 1)
dev.off()

####### 5.3 不同亚组之间的表达谱表达 主要是级别和正常组和肿瘤组之间 ##########

source('/export3/zhangw/Project_Cross/Project.Glioma.SGene/code/plot.glioma.subgroup.gene.R')
your_color_pal=c('#375631','#D3C3D0','#A786BA','#6D4E7E',paletteer_c("grDevices::Sunset", 30)[c(6,12,18,24,30)] )
#TCGA
tcga.sub.list= c('Histology','Grade','Age','Gender', 'IDH_status',
                 '1p_19q_codeletion','MGMT_promoter_status',
                 'Transcriptome_Subtype',
                 'radiation_therapy')
subplot = list()
for (i in  tcga.sub.list ) {
  p = plot.glioma.subgroup.gene(subgroup_you_want = i,
                                gene_you_want = gene_you_want,
                                datasets = 'TCGA',
                                your_color_pal = your_color_pal     
                                )
  subplot[[i]] = p
}
# TCGA 

cairo_pdf('./fig/fig_5c_glioma_histology_transcriptomic_subtype_tcga.pdf',width = 5,height = 4.8,onefile = F)
subplot[[1]]+
  # subplot[[2]]+
  # subplot[[3]]+
  # subplot[[4]]+
  # subplot[[5]]+
  # subplot[[6]]+
  # subplot[[7]]+
  subplot[[8]]+
  # subplot[[9]]+
  plot_layout(ncol = 2,nrow = 1)
dev.off()

cairo_pdf('./fig/fig_5c_glioma_subtypes_tcga.pdf',width = 10,height = 9.6,onefile = F)
# subplot[[1]]+
subplot[[2]]+
  subplot[[3]]+
  subplot[[4]]+
  subplot[[5]]+
  subplot[[6]]+
  subplot[[7]]+
  # subplot[[8]]+
  subplot[[9]]+
  plot_layout(ncol = 4)
dev.off()

# CGGA.1018


cgga_subgroup_list =c('Grade','Age','Gender', 'IDH_status',
                      '1p_19q_codeletion','MGMT_promoter_status','Chemo_therapy',
                      'radiation_therapy','primary_recurrent')
subplot = list()
for (i in cgga_subgroup_list ){
  
  p = plot.glioma.subgroup.gene(subgroup_you_want = i,
                                gene_you_want = gene_you_want,
                                datasets = 'CGGA.1018',
                                your_color_pal = your_color_pal
                                )
  subplot[[i]] = p
  
}

cairo_pdf('./fig/fig_5c_glioma_subtypes_cgga.1018.pdf',width = 7.5,height = 14.4,onefile = F)
subplot[[1]]+
  subplot[[2]]+
  subplot[[3]]+
  subplot[[4]]+
  subplot[[5]]+
  subplot[[6]]+
  subplot[[7]]+
  subplot[[8]]+
  subplot[[9]]+
  plot_layout(ncol = 3)
dev.off()

subplot = list()
for (i in cgga_subgroup_list ) {
  
  p = plot.glioma.subgroup.gene(subgroup_you_want = i,
                                gene_you_want = gene_you_want,
                                datasets = 'CGGA.325',
                                your_color_pal =your_color_pal
                                
                                
                                )
  subplot[[i]] = p
  
}

cairo_pdf('./fig/fig_5c_glioma_subtypes_cgga.325.pdf',width = 7.5,height = 14.4,onefile = F)
subplot[[1]]+
  subplot[[2]]+
  subplot[[3]]+
  subplot[[4]]+
  subplot[[5]]+
  subplot[[6]]+
  subplot[[7]]+
  subplot[[8]]+
  subplot[[9]]+
  plot_layout(ncol = 3)
dev.off()


subplot = list()
for (i in cgga_subgroup_list ) {
  
  p = plot.glioma.subgroup.gene(subgroup_you_want = i,
                                gene_you_want = gene_you_want,
                                datasets = 'CGGA.693',
                                your_color_pal =your_color_pal
                                
                                )
  subplot[[i]] = p
  
}

cairo_pdf('./fig/fig_5c_glioma_subtypes_cgga.693.pdf',width = 7.5,height = 14.4,onefile = F)
subplot[[1]]+
  subplot[[2]]+
  subplot[[3]]+
  subplot[[4]]+
  subplot[[5]]+
  subplot[[6]]+
  subplot[[7]]+
  subplot[[8]]+
  subplot[[9]]+
  plot_layout(ncol = 3)
dev.off()


source("/export3/zhangw/Project_Cross/Project.Glioma.SGene/code/single.gene.distribution.tumor.vs.normal.TCGA.GTEX.R")

p=single.gene.distribution.tumor.vs.normal.TCGA.GTEX(Glioma.vs.N = T,
                                                     Grade.N = F,
                                                     gene_you_want = gene_you_want,
                                                     your_color_pal = your_color_pal
)
p+ ylab(paste0('log2 of ', gene_you_want)) +
  xlab("") +
  labs(title = 'TCGA+GTEx',subtitle = "T(689) vs N(105)")+
  #ylim(0, 10)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle  = element_text(hjust = 0.5),
        legend.position = "none")
# p+ylim(-2.5,15)
ggsave("./fig/fig_5c_Expression.normal.tumor.gtex.tcga.glioma.pdf",width = 2.5,height = 4.8)


p=single.gene.distribution.tumor.vs.normal.TCGA.GTEX(Glioma.vs.N = F,
                                                     Grade.N = T,
                                                     gene_you_want = gene_you_want,
                                                     your_color_pal =your_color_pal
                                                     
)

p+ ylab(paste0('log2 of ', gene_you_want)) +
  xlab("") +
  labs(title = 'TCGA+GTEx',subtitle = "T(689) vs N(105)")+
  #ylim(0, 10)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle  = element_text(hjust = 0.5),
        legend.position = "none")
# p+ylim(-2.5,15)
ggsave("./fig/fig_5c_Expression.grade.gtex.tcga.glioma.pdf",width = 2.5,height = 5)

# GSE4290
gene_you_want = 'MEST'
your_color_pal =  c('#375631','#D3C3D0','#A786BA','#6D4E7E',paletteer_c("grDevices::Sunset", 30)[c(6,12,18,24,30)] )

 load("/export3/zhangw/Glioma.dataset/DATABASE/GEO/GSE4290/GSE4290_pheno.RData")
 load("/export3/zhangw/Glioma.dataset/DATABASE/GEO/GSE4290/GSE4290_array_log2_normalized.RData")

 tmp = left_join(GSE4290_pheno[,c('Sample','Grade')],rownames_to_column(as.data.frame(t(exprSet[gene_you_want,])),'Sample'))
 tmp$Grade[tmp$Grade=='Epilepsy Brain']='Normal'
 tmp$Grade[tmp$Grade=='Normal Brain']='Normal'
 tmp$Grade[tmp$Grade=='grade 2']='G2'
 tmp$Grade[tmp$Grade=='grade 3']='G3'
 tmp$Grade[tmp$Grade=='grade 4']='G4'
 table(tmp$Grade)
 
 tmp.plot = tmp %>% dplyr::rename('group'= 'Grade','var'='MEST')
 # tumor vs normal 
 compare.list = list( c('Normal','G2'),
                      c('Normal','G3'),
                      c('Normal','G4'),
                      c('G2','G3'),
                      c('G2','G4'),
                      c('G3','G4')
                      )
 tmp.plot$group = factor(tmp.plot$group, levels = c('Normal','G2','G3','G4'))
 p1 =ggplot(tmp.plot,aes(x=group,y=var))+
   geom_beeswarm(size=2.5,aes(color=group),cex =2)+
   stat_summary(fun.y = "mean", geom = "crossbar", #添加均值
                mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.15,size=0.3,color="black")+
   stat_summary(fun.data = "mean_se", geom="errorbar", width = 0.1,size=0.5,color="black")+#添加误差线
   # facet_wrap(.~Group)+
   scale_color_manual(values =your_color_pal[1:length(unique(tmp.plot$group))]
                      # labels=c(bquote(italic(group) * "WT"),bquote(italic(group) * "mt"))
                     )+
   
   # scale_y_continuous(limits = c(0,5.3),expand = c(0,0))+
   # scale_y_continuous(limits = c(0,5.3),expand = c(0,0))+
   labs(x="",y="log2 of MEST")+
   theme_classic()+
   theme(legend.position = "top",
         legend.title = element_blank(),
         strip.background = element_blank(),
         strip.text = element_text(size = 12,face = "bold"),
         axis.text.x = element_blank(),
         axis.ticks = element_line(size=0.2, color="black"),
         axis.title = element_text(size = 10),
         plot.title = element_text(hjust = 0.5),
         plot.subtitle  = element_text(hjust = 0.5),
         panel.spacing =unit(0,'cm'))+
   ggsignif::geom_signif(comparisons = compare.list,step_increase = 0.1,
                         test = "t.test")+
   ggtitle('GSE4290')+
   guides(color = guide_legend(nrow = 1))+
   guides(fill=guide_legend(nrow =1))#设置图例显示为两行
p1
 
 tmp.plot = tmp %>% dplyr::rename('group'= 'Grade','var'='MEST')
 tmp.plot$group[tmp.plot$group!='Normal'] = 'Tumor'
 # tumor vs normal 
 compare.list = list( c('Normal','Tumor')
                   
 )
 tmp.plot$group = factor(tmp.plot$group, levels = c('Normal','Tumor'))
 
 p2 =ggplot(tmp.plot,aes(x=group,y=var))+
   geom_beeswarm(size=2.5,aes(color=group),cex =2)+
   stat_summary(fun.y = "mean", geom = "crossbar", #添加均值
                mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.15,size=0.3,color="black")+
   stat_summary(fun.data = "mean_se", geom="errorbar", width = 0.1,size=0.5,color="black")+#添加误差线
   # facet_wrap(.~Group)+
   scale_color_manual(values =your_color_pal[1:length(unique(tmp.plot$group))]
                      # labels=c(bquote(italic(group) * "WT"),bquote(italic(group) * "mt"))
                     )+
   
   # scale_y_continuous(limits = c(0,5.3),expand = c(0,0))+
   # scale_y_continuous(limits = c(0,5.3),expand = c(0,0))+
   labs(x="",y="log2 of MEST")+
   theme_classic()+
   theme(legend.position = "top",
         legend.title = element_blank(),
         strip.background = element_blank(),
         strip.text = element_text(size = 12,face = "bold"),
         axis.text.x = element_blank(),
         axis.title = element_text(size = 10),
         plot.title = element_text(hjust = 0.5),
         plot.subtitle  = element_text(hjust = 0.5),
         panel.spacing =unit(0,'cm'))+
   ggsignif::geom_signif(comparisons = compare.list,step_increase = 0.1,
                         test = "t.test")+
   ggtitle('GSE4290')+
   guides(color = guide_legend(nrow = 1))+
   guides(fill=guide_legend(nrow =1))#设置图例显示为两行
p2

pdf("./fig/fig_5c_Expression.grade.t.n.gse4290.glioma.pdf",width = 5,height = 4.8,onefile = F)
p1+p2+ plot_layout(nrow =1,
                    # heights = c(1,1,1),
                    widths = c(1,1),
                    guides='collect' )&
    theme(legend.position='top')
dev.off()

### GSE16011


gene_you_want = 'MEST'
your_color_pal =  c('#375631','#D3C3D0','#A786BA','#6D4E7E',paletteer_c("grDevices::Sunset", 30)[c(6,12,18,24,30)] )

load("/export3/zhangw/Glioma.dataset/DATABASE/GEO/GSE16011/GSE16011_array_normalized.RData")
load("/export3/zhangw/Glioma.dataset/DATABASE/GEO/GSE16011/GSE16011_pheno.RData")

tmp = left_join(rownames_to_column(GSE16011_pheno,'Sample')[,c('Sample','grade')],
                rownames_to_column(as.data.frame(t(GSE16011_array[gene_you_want,])),'Sample')) %>% dplyr::rename('Grade'='grade')
tmp = tmp %>% dplyr::filter(!is.na(tmp$Grade))
tmp$Grade[tmp$Grade=='I']='G1'
tmp$Grade[tmp$Grade=='II']='G2'
tmp$Grade[tmp$Grade=='III']='G3'
tmp$Grade[tmp$Grade=='IV']='G4'
tmp = tmp %>% dplyr::filter(tmp$Grade %in%  c('G1','G2','G3','G4'))
table(tmp$Grade)

tmp.plot = tmp %>% dplyr::rename('group'= 'Grade','var'='MEST')
# tumor vs normal 
compare.list = list( c('G1','G2'),
                     c('G1','G3'),
                     c('G1','G4'),
                     c('G2','G3'),
                     c('G2','G4'),
                     c('G3','G4')
)
tmp.plot$group = factor(tmp.plot$group, levels = c('G1','G2','G3','G4'))
p1 =ggplot(tmp.plot,aes(x=group,y=var))+
  geom_beeswarm(size=2.5,aes(color=group),cex =1.5)+
  stat_summary(fun.y = "mean", geom = "crossbar", #添加均值
               mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.15,size=0.3,color="black")+
  stat_summary(fun.data = "mean_se", geom="errorbar", width = 0.1,size=0.5,color="black")+#添加误差线
  # facet_wrap(.~Group)+
  scale_color_manual(values =your_color_pal[1:length(unique(tmp.plot$group))]
                     # labels=c(bquote(italic(group) * "WT"),bquote(italic(group) * "mt"))
  )+
  
  # scale_y_continuous(limits = c(0,5.3),expand = c(0,0))+
  # scale_y_continuous(limits = c(0,5.3),expand = c(0,0))+
  labs(x="",y="log2 of MEST")+
  theme_classic()+
  theme(legend.position = "top",
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 12,face = "bold"),
        axis.text.x = element_blank(),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.title = element_text(size = 10),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle  = element_text(hjust = 0.5),
        panel.spacing =unit(0,'cm'))+
  ggsignif::geom_signif(comparisons = compare.list,step_increase = 0.1,
                        test = "t.test")+
  ggtitle('GSE16011')+
  guides(color = guide_legend(nrow = 1))+
  guides(fill=guide_legend(nrow =1))#设置图例显示为两行
p1

tmp.plot = tmp %>% dplyr::rename('group'= 'Grade','var'='MEST')
tmp.plot$group[tmp.plot$group!='G4'] = 'LGG'
tmp.plot$group[tmp.plot$group=='G4'] = 'GBM'
# tumor vs normal 
compare.list = list( c('LGG','GBM')
                     
)
tmp.plot$group = factor(tmp.plot$group, levels = c('LGG','GBM'))

p2 =ggplot(tmp.plot,aes(x=group,y=var))+
  geom_beeswarm(size=2.5,aes(color=group),cex =2)+
  stat_summary(fun.y = "mean", geom = "crossbar", #添加均值
               mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.15,size=0.3,color="black")+
  stat_summary(fun.data = "mean_se", geom="errorbar", width = 0.1,size=0.5,color="black")+#添加误差线
  # facet_wrap(.~Group)+
  scale_color_manual(values =your_color_pal[1:length(unique(tmp.plot$group))]
                     # labels=c(bquote(italic(group) * "WT"),bquote(italic(group) * "mt"))
  )+
  
  # scale_y_continuous(limits = c(0,5.3),expand = c(0,0))+
  # scale_y_continuous(limits = c(0,5.3),expand = c(0,0))+
  labs(x="",y="log2 of MEST")+
  theme_classic()+
  theme(legend.position = "top",
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 12,face = "bold"),
        axis.text.x = element_blank(),
        axis.title = element_text(size = 10),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle  = element_text(hjust = 0.5),
        panel.spacing =unit(0,'cm'))+
  ggsignif::geom_signif(comparisons = compare.list,step_increase = 0.1,
                        test = "t.test")+
  ggtitle('GSE16011')+
  guides(color = guide_legend(nrow = 1))+
  guides(fill=guide_legend(nrow =1))#设置图例显示为两行
p2

pdf("./fig/fig_5c_Expression.grade.LGG.GBM.gse16011.glioma.pdf",width = 5,height = 4.8,onefile = F)
p1+p2+ plot_layout(nrow =1,
                   # heights = c(1,1,1),
                   widths = c(1,1),
                   guides='collect' )&
  theme(legend.position='top')
dev.off()

# GSE35493
load("/export3/zhangw/Glioma.dataset/DATABASE/GEO/GSE35493/GSE35493_GBM_pheno.RData")
load("/export3/zhangw/Glioma.dataset/DATABASE/GEO/GSE35493/GSE35493_Normal_pheno.RData")
load("/export3/zhangw/Glioma.dataset/DATABASE/GEO/GSE35493/GSE35493_array_GBM_normalized.RData")
load("/export3/zhangw/Glioma.dataset/DATABASE/GEO/GSE35493/GSE35493_array_Normal_normalized.RData")

tmp = cbind(dataset_GBM,dataset_Normal[rownames(dataset_GBM),]) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column('Sample') %>% 
  dplyr::select(c('Sample','MEST'))

tmp_p = rbind(pheno_GBM,pheno_Normal) %>% dplyr::select(c('Sample','Tissue'))

tmp_plot = left_join(tmp, tmp_p)
tmp_plot = tmp_plot %>% dplyr::rename('group'= 'Tissue','var'='MEST')
tmp_plot$group[tmp_plot$group =='glioblastoma']= 'Tumor'
tmp_plot$group[tmp_plot$group !='Tumor']= 'Normal'

compare.list = list( c('Normal','Tumor')
                     
)
tmp.plot = tmp_plot
tmp.plot$group = factor(tmp.plot$group, levels = c('Normal','Tumor'))

p2 =ggplot(tmp.plot,aes(x=group,y=var))+
  geom_beeswarm(size=2.5,aes(color=group),cex =6)+
  stat_summary(fun.y = "mean", geom = "crossbar", #添加均值
               mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.15,size=0.3,color="black")+
  stat_summary(fun.data = "mean_se", geom="errorbar", width = 0.1,size=0.5,color="black")+#添加误差线
  # facet_wrap(.~Group)+
  scale_color_manual(values =your_color_pal[1:length(unique(tmp.plot$group))]
                     # labels=c(bquote(italic(group) * "WT"),bquote(italic(group) * "mt"))
  )+
  
  # scale_y_continuous(limits = c(0,5.3),expand = c(0,0))+
  # scale_y_continuous(limits = c(0,5.3),expand = c(0,0))+
  labs(x="",y="log2 of MEST")+
  theme_classic()+
  theme(legend.position = "top",
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 12,face = "bold"),
        axis.text.x = element_blank(),
        axis.title = element_text(size = 10),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle  = element_text(hjust = 0.5),
        panel.spacing =unit(0,'cm'))+
  ggsignif::geom_signif(comparisons = compare.list,step_increase = 0.1,
                        test = "wilcox.test"  , test.args = c("two.sided"))+
  ggtitle('GSE35493')+
  guides(color = guide_legend(nrow = 1))+
  guides(fill=guide_legend(nrow =1))#设置图例显示为两行
p2

pdf("./fig/fig_5c_Expression.grade.normal_t_gse35493.glioma.pdf",width = 2.5,height = 4.8,onefile = F)
p2
# +p2+ plot_layout(nrow =1,
#                    # heights = c(1,1,1),
#                    widths = c(1,1),
#                    guides='collect' )&
#   theme(legend.position='top')
dev.off()

### GSE42656 

gene_you_want = 'MEST'
your_color_pal =  c('#375631','#D3C3D0','#A786BA','#6D4E7E',paletteer_c("grDevices::Sunset", 30)[c(6,12,18,24,30)] )

load("/export3/zhangw/Glioma.dataset/DATABASE/GEO/GSE42656/GSE42656_pheno.RData")
load("/export3/zhangw/Glioma.dataset/DATABASE/GEO/GSE42656/GSE42656_array.RData")

tmp = left_join(pheno[,c('Sample','Grade')],rownames_to_column(as.data.frame(t(exprSet[gene_you_want,])),'Sample'))


tmp$Grade[tmp$Grade%in% c('F','M','ND')]='Normal'
tmp$Grade[tmp$Grade=='I']='G1'
tmp$Grade[tmp$Grade=='II']='G2'
tmp$Grade[tmp$Grade=='III']='G3'
tmp$Grade[tmp$Grade=='IV']='G4'
table(tmp$Grade)

tmp.plot = tmp %>% dplyr::rename('group'= 'Grade','var'='MEST')
# tumor vs normal 
compare.list = list( 
                     c('Normal','G1'),
                     c('Normal','G2'),
                     c('Normal','G3'),
                     c('Normal','G4'),
                     c('G1','G2'),
                     c('G1','G3'),
                     c('G1','G4'),
                     c('G2','G3'),
                     c('G2','G4'),
                     c('G3','G4')
)
tmp.plot$group = factor(tmp.plot$group, levels = c('Normal','G1','G2','G3','G4'))
p1 =ggplot(tmp.plot,aes(x=group,y=var))+
  geom_beeswarm(size=2.5,aes(color=group),cex =2)+
  stat_summary(fun.y = "mean", geom = "crossbar", #添加均值
               mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.15,size=0.3,color="black")+
  stat_summary(fun.data = "mean_se", geom="errorbar", width = 0.1,size=0.5,color="black")+#添加误差线
  # facet_wrap(.~Group)+
  scale_color_manual(values =your_color_pal[1:length(unique(tmp.plot$group))]
                     # labels=c(bquote(italic(group) * "WT"),bquote(italic(group) * "mt"))
  )+
  
  # scale_y_continuous(limits = c(0,5.3),expand = c(0,0))+
  # scale_y_continuous(limits = c(0,5.3),expand = c(0,0))+
  labs(x="",y="log2 of MEST")+
  theme_classic()+
  theme(legend.position = "top",
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 12,face = "bold"),
        axis.text.x = element_blank(),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.title = element_text(size = 10),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle  = element_text(hjust = 0.5),
        panel.spacing =unit(0,'cm'))+
  ggsignif::geom_signif(comparisons = compare.list,step_increase = 0.1,
                        test = "t.test")+
  ggtitle('GSE42656')+
  guides(color = guide_legend(nrow = 1))+
  guides(fill=guide_legend(nrow =1))#设置图例显示为两行
p1

tmp.plot = tmp %>% dplyr::rename('group'= 'Grade','var'='MEST')
tmp.plot$group[tmp.plot$group!='Normal'] = 'Tumor'
# tumor vs normal 
compare.list = list( c('Normal','Tumor')
                     
)
tmp.plot$group = factor(tmp.plot$group, levels = c('Normal','Tumor'))

p2 =ggplot(tmp.plot,aes(x=group,y=var))+
  geom_beeswarm(size=2.5,aes(color=group),cex =2)+
  stat_summary(fun.y = "mean", geom = "crossbar", #添加均值
               mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.15,size=0.3,color="black")+
  stat_summary(fun.data = "mean_se", geom="errorbar", width = 0.1,size=0.5,color="black")+#添加误差线
  # facet_wrap(.~Group)+
  scale_color_manual(values =your_color_pal[1:length(unique(tmp.plot$group))]
                     # labels=c(bquote(italic(group) * "WT"),bquote(italic(group) * "mt"))
  )+
  
  # scale_y_continuous(limits = c(0,5.3),expand = c(0,0))+
  # scale_y_continuous(limits = c(0,5.3),expand = c(0,0))+
  labs(x="",y="log2 of MEST")+
  theme_classic()+
  theme(legend.position = "top",
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 12,face = "bold"),
        axis.text.x = element_blank(),
        axis.title = element_text(size = 10),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle  = element_text(hjust = 0.5),
        panel.spacing =unit(0,'cm'))+
  ggsignif::geom_signif(comparisons = compare.list,step_increase = 0.1,
                        test = "t.test")+
  ggtitle('GSE42656')+
  guides(color = guide_legend(nrow = 1))+
  guides(fill=guide_legend(nrow =1))#设置图例显示为两行
p2

pdf("./fig/fig_5c_Expression.grade.t.n.GSE42656.glioma.pdf",width = 5,height = 4.8,onefile = F)
p1+p2+ plot_layout(nrow =1,
                   # heights = c(1,1,1),
                   widths = c(1,1),
                   guides='collect' )&
  theme(legend.position='top')
dev.off()

# GSE43378

 
 gene_you_want = 'MEST'
 your_color_pal =  c('#375631','#D3C3D0','#A786BA','#6D4E7E',paletteer_c("grDevices::Sunset", 30)[c(6,12,18,24,30)] )
 
 load("/export3/zhangw/Glioma.dataset/DATABASE/GEO/GSE43378/GSE43378_pheno.RData")
 load("/export3/zhangw/Glioma.dataset/DATABASE/GEO/GSE43378/GSE43378_array_normalized.RData")

  tmp = left_join(pheno[,c('Sample','who grade')],rownames_to_column(as.data.frame(t(dataset_normalized[gene_you_want,])),'Sample'))
 
 colnames(tmp)[2]='Grade'
 tmp$Grade[tmp$Grade=='2']='G2'
 tmp$Grade[tmp$Grade=='3']='G3'
 tmp$Grade[tmp$Grade=='4']='G4'
 table(tmp$Grade)
 
 tmp.plot = tmp %>% dplyr::rename('group'= 'Grade','var'='MEST')
 # tumor vs normal 
 compare.list = list( 
   # c('Normal','G1'),
   # c('Normal','G2'),
   # c('Normal','G3'),
   # c('Normal','G4'),
   # c('G1','G2'),
   # c('G1','G3'),
   # c('G1','G4'),
   c('G2','G3'),
   c('G2','G4'),
   c('G3','G4')
 )
 tmp.plot$group = factor(tmp.plot$group, levels = c('G2','G3','G4'))
 p1 =ggplot(tmp.plot,aes(x=group,y=var))+
   geom_beeswarm(size=2.5,aes(color=group),cex =2)+
   stat_summary(fun.y = "mean", geom = "crossbar", #添加均值
                mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.15,size=0.3,color="black")+
   stat_summary(fun.data = "mean_se", geom="errorbar", width = 0.1,size=0.5,color="black")+#添加误差线
   # facet_wrap(.~Group)+
   scale_color_manual(values =your_color_pal[1:length(unique(tmp.plot$group))]
                      # labels=c(bquote(italic(group) * "WT"),bquote(italic(group) * "mt"))
   )+
   
   # scale_y_continuous(limits = c(0,5.3),expand = c(0,0))+
   # scale_y_continuous(limits = c(0,5.3),expand = c(0,0))+
   labs(x="",y="log2 of MEST")+
   theme_classic()+
   theme(legend.position = "top",
         legend.title = element_blank(),
         strip.background = element_blank(),
         strip.text = element_text(size = 12,face = "bold"),
         axis.text.x = element_blank(),
         axis.ticks = element_line(size=0.2, color="black"),
         axis.title = element_text(size = 10),
         plot.title = element_text(hjust = 0.5),
         plot.subtitle  = element_text(hjust = 0.5),
         panel.spacing =unit(0,'cm'))+
   ggsignif::geom_signif(comparisons = compare.list,step_increase = 0.1,
                         test = "t.test")+
   ggtitle('GSE43378')+
   guides(color = guide_legend(nrow = 1))+
   guides(fill=guide_legend(nrow =1))#设置图例显示为两行
 p1
 
 tmp.plot = tmp %>% dplyr::rename('group'= 'Grade','var'='MEST')
 tmp.plot$group[tmp.plot$group!='G4'] = 'LGG'
 tmp.plot$group[tmp.plot$group=='G4'] = 'GBM'
 # tumor vs normal 
 compare.list = list( c('LGG','GBM')
                      
 )
 tmp.plot$group = factor(tmp.plot$group, levels = c('LGG','GBM'))
 
 p2 =ggplot(tmp.plot,aes(x=group,y=var))+
   geom_beeswarm(size=2.5,aes(color=group),cex =2)+
   stat_summary(fun.y = "mean", geom = "crossbar", #添加均值
                mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.15,size=0.3,color="black")+
   stat_summary(fun.data = "mean_se", geom="errorbar", width = 0.1,size=0.5,color="black")+#添加误差线
   # facet_wrap(.~Group)+
   scale_color_manual(values =your_color_pal[1:length(unique(tmp.plot$group))]
                      # labels=c(bquote(italic(group) * "WT"),bquote(italic(group) * "mt"))
   )+
   
   # scale_y_continuous(limits = c(0,5.3),expand = c(0,0))+
   # scale_y_continuous(limits = c(0,5.3),expand = c(0,0))+
   labs(x="",y="log2 of MEST")+
   theme_classic()+
   theme(legend.position = "top",
         legend.title = element_blank(),
         strip.background = element_blank(),
         strip.text = element_text(size = 12,face = "bold"),
         axis.text.x = element_blank(),
         axis.title = element_text(size = 10),
         plot.title = element_text(hjust = 0.5),
         plot.subtitle  = element_text(hjust = 0.5),
         panel.spacing =unit(0,'cm'))+
   ggsignif::geom_signif(comparisons = compare.list,step_increase = 0.1,
                         test = "t.test")+
   ggtitle('GSE43378')+
   guides(color = guide_legend(nrow = 1))+
   guides(fill=guide_legend(nrow =1))#设置图例显示为两行
 p2

 pdf("./fig/fig_5c_Expression.grade.GSE43378.glioma.pdf",width = 5,height = 4.8,onefile = F)
 p1+p2+ plot_layout(nrow =1,
                    # heights = c(1,1,1),
                    widths = c(1,1),
                    guides='collect' )&
   theme(legend.position='top')
dev.off()

# GSE108474
 load("/export3/zhangw/Glioma.dataset/DATABASE/GEO/GSE108474/GSE108474_GBM_pheno.RData")
 load("/export3/zhangw/Glioma.dataset/DATABASE/GEO/GSE108474/GSE108474_Normal_pheno.RData")
 load("/export3/zhangw/Glioma.dataset/DATABASE/GEO/GSE108474/GSE108474_O_pheno.RData")
load("/export3/zhangw/Glioma.dataset/DATABASE/GEO/GSE108474/GSE108474_MIXED_pheno.RData")
load("/export3/zhangw/Glioma.dataset/DATABASE/GEO/GSE108474/GSE108474_A_pheno.RData")
load("/export3/zhangw/Glioma.dataset/DATABASE/GEO/GSE108474/GSE108474_pheno.RData")
load("/export3/zhangw/Glioma.dataset/DATABASE/GEO/GSE108474/GSE108474_array_normalized.RData")


pheo.id = pheno
pheo.id$Grade= ifelse(pheo.id$Sample%in%pheno_GBM$Sample, 'GBM',
                      ifelse(pheo.id$Sample%in% pheno_Normal$Sample, 'Normal', 'LGG'))

tmp = left_join(pheo.id[,c('Sample','Grade')],rownames_to_column(as.data.frame(t(dataset[gene_you_want,])),'Sample'))

colnames(tmp)[2]='Grade'
# tmp$Grade[tmp$Grade=='2']='G2'
# tmp$Grade[tmp$Grade=='3']='G3'
# tmp$Grade[tmp$Grade=='4']='G4'
table(tmp$Grade)

tmp.plot = tmp %>% dplyr::rename('group'= 'Grade','var'='MEST')
# tumor vs normal 
compare.list = list( 

  c('Normal','LGG'),
  c('Normal','GBM'),
  c('LGG','GBM')
)
tmp.plot$group = factor(tmp.plot$group, levels = c('Normal','LGG','GBM'))
p1 =ggplot(tmp.plot,aes(x=group,y=var))+
  geom_beeswarm(size=2.5,aes(color=group),cex =1)+
  stat_summary(fun.y = "mean", geom = "crossbar", #添加均值
               mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.15,size=0.3,color="black")+
  stat_summary(fun.data = "mean_se", geom="errorbar", width = 0.1,size=0.5,color="black")+#添加误差线
  # facet_wrap(.~Group)+
  scale_color_manual(values =your_color_pal[1:length(unique(tmp.plot$group))]
                     # labels=c(bquote(italic(group) * "WT"),bquote(italic(group) * "mt"))
  )+
  
  # scale_y_continuous(limits = c(0,5.3),expand = c(0,0))+
  # scale_y_continuous(limits = c(0,5.3),expand = c(0,0))+
  labs(x="",y="log2 of MEST")+
  theme_classic()+
  theme(legend.position = "top",
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 12,face = "bold"),
        axis.text.x = element_blank(),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.title = element_text(size = 10),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle  = element_text(hjust = 0.5),
        panel.spacing =unit(0,'cm'))+
  ggsignif::geom_signif(comparisons = compare.list,step_increase = 0.1,
                        test = "t.test")+
  ggtitle('GSE108474')+
  guides(color = guide_legend(nrow = 1))+
  guides(fill=guide_legend(nrow =1))#设置图例显示为两行
p1

tmp.plot = tmp %>% dplyr::rename('group'= 'Grade','var'='MEST')
tmp.plot$group[tmp.plot$group!='Normal'] = 'Tumor'
# tmp.plot$group[tmp.plot$group=='G4'] = 'GBM'
# tumor vs normal 
compare.list = list( c('Normal','Tumor')
                     
)
tmp.plot$group = factor(tmp.plot$group, levels = c('Normal','Tumor'))
table(tmp.plot$group)
p2 =ggplot(tmp.plot,aes(x=group,y=var))+
  geom_beeswarm(size=2.5,aes(color=group),cex =2)+
  stat_summary(fun.y = "mean", geom = "crossbar", #添加均值
               mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.15,size=0.3,color="black")+
  stat_summary(fun.data = "mean_se", geom="errorbar", width = 0.1,size=0.5,color="black")+#添加误差线
  # facet_wrap(.~Group)+
  scale_color_manual(values =your_color_pal[1:length(unique(tmp.plot$group))]
                     # labels=c(bquote(italic(group) * "WT"),bquote(italic(group) * "mt"))
  )+
  
  # scale_y_continuous(limits = c(0,5.3),expand = c(0,0))+
  # scale_y_continuous(limits = c(0,5.3),expand = c(0,0))+
  labs(x="",y="log2 of MEST")+
  theme_classic()+
  theme(legend.position = "top",
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 12,face = "bold"),
        axis.text.x = element_blank(),
        axis.title = element_text(size = 10),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle  = element_text(hjust = 0.5),
        panel.spacing =unit(0,'cm'))+
  ggsignif::geom_signif(comparisons = compare.list,step_increase = 0.1,
                        test = "t.test")+
  ggtitle('GSE108474')+
  guides(color = guide_legend(nrow = 1))+
  guides(fill=guide_legend(nrow =1))#设置图例显示为两行
p2

pdf("./fig/fig_5c_Expression.grade.GSE108474.glioma.pdf",width = 5,height = 4.8,onefile = F)
p1+p2+ plot_layout(nrow =1,
                   # heights = c(1,1,1),
                   widths = c(1,1),
                   guides='collect' )&
  theme(legend.position='top')
dev.off()



####### 5.4 pan cancer 肿瘤和正常的表达 ###########

source('/export3/zhangw/pancancer/tcga.mrna/plot.special.gene.pancancer.R')

p  = plot.special.gene.pancancer(special.gene = gene_you_want,mycol = c(paletteer_c("grDevices::Sunset", 30)[c(6,24,30)] ))
p+  theme_classic()+ 
  theme(legend.position = "top",
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 12,face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        plot.title = element_text(hjust = 0.5),
        panel.spacing =unit(0,'cm'))+
  ylab(paste0('log2 of ', gene_you_want)) +
  xlab("") +
  labs(title = 'TCGA+GTEx')+
  ylim(-10,15)
  

ggsave("./fig/fig_5c_Expression.pancancer.pdf",width = 8,height = 3.6)


####### 5.5 多个队列进行多因素回归 ############
#multicox
library(ezcox)
library(forestploter)
library(grid)

load("/export3/zhangw/Project_Cross/Project.Glioma.SGene/data/meta.tcga.cgga.Rdata")
datasets<-"TCGA"

meta.tcga = meta$TCGA 
expr.tcga = list_train_vali_Data$TCGA
multic_cox_model<-list()
multic_cox<-data.frame()
expr.tcga = expr.tcga[,c('ID','OS.time','OS',gene_you_want)]
  meta.tcga = meta$TCGA 
  

  rs = expr.tcga
  rownames(rs) =rs$ID
  colnames(rs)[4] ='RS' 
  rownames(rs)<-rs$ID
  # rownames(rs)<-gsub("\\.","-",rownames(rs))
  # rownames(rs)<-gsub("\\.","-",rownames(rs))
  # meta.tcga<-meta.tcga[rownames(rs),]
  meta.tcga<-cbind(meta.tcga,rs[rownames(meta.tcga),])
  
  data<-dplyr::select(meta.tcga,c("Age (years at diagnosis)","Gender","Grade","IDH status","1p/19q codeletion","MGMT promoter status",
                                  "OS","OS.time","RS"))
  colnames(data)<-c("Age","Gender","Grade","IDH","Chr1p/19q","MGMT",
                    "status","time","RS")
  data$Age<-as.numeric(data$Age)
  data<-data[data$IDH !="NA",]
  data<-data[data$`Chr1p/19q` !="NA",]
  data<-data[data$MGMT !="NA",]
  data<-data[!is.na(data$RS ),]
  
  data$`Risk score`<-ifelse(data[,"RS"]>median(data[,"RS"]),"High","Low")
  data$`Risk score`<-factor(data$`Risk score`,levels=c("Low","High"))
  
  data$Age<-ifelse(data$Age>50,">50","<50")
  data$Age<-factor(data$Age,levels = c("<50",">50"))
  data = as.data.frame.matrix(data)
  
  data[data$Gender=="female","Gender"]<-"Female"
  data[data$Gender=="male",'Gender']<-"Male"
  data$Gender<-factor(data$Gende,levels=c("Female","Male"))
  data$Grade<-factor(data$Grade,levels=c("G2","G3","G4"))
  data[data$IDH=="WT","IDH"]<-"Wildtype"
  data$IDH<-factor(data$IDH,levels=c("Wildtype","Mutant"))
  data$MGMT<-factor(data$MGMT,levels=c("Unmethylated","Methylated"))
  data[data$`Chr1p/19q`=="non-codel","Chr1p/19q"]<-"Non-codel"
  data[data$`Chr1p/19q`=="codel","Chr1p/19q"]<-"Codel"
  data$`Chr1p/19q`<-factor(data$`Chr1p/19q`,levels=c("Non-codel","Codel"))
  
  colnames(data)[10] = gene_you_want
  data$MEST = factor(data$MEST, levels =  c("Low","High"))
  ezcox_res <- ezcox(data,
                     covariates = c(gene_you_want),
                     controls = c("Gender","Age","Grade","IDH","Chr1p/19q","MGMT"),
                     return_models=TRUE)
  
  tmp<-as.data.frame(ezcox_res[["res"]])
  tmp$Variable<-c('MEST',"Gender","Age","Grade","","IDH","Chr1p/19q","MGMT")
  tmp$Cohorts<-c(datasets,"","","","","","","")
  tmp$Cohorts<-datasets
  multic_cox<-rbind(multic_cox,tmp)
  multic_cox_model[[datasets]]<-get_models(ezcox_res)
  show_models(multic_cox_model[[datasets]])
  




#CGGA
for (i in c("CGGA.325","CGGA.693","CGGA.1018")) {
  datasets<-i
  expr = list_train_vali_Data[[i]]
  
    # for (j in 4:50) {
    meta.cgga <- meta$CGGA
    
    rs = expr[,c('ID','OS.time','OS',gene_you_want)]
    rownames(rs) =rs$ID
    colnames(rs)[4] ='RS' 
    rownames(meta.cgga)<-meta.cgga$CGGA_ID
    meta.cgga<-meta.cgga[rownames(rs),]
    meta.cgga<-cbind(meta.cgga[,-7],rs)
    data<-dplyr::select(meta.cgga,c("Age","Gender","Grade","IDH_mutation_status","1p19q_codeletion_status","MGMTp_methylation_status",
                                    "OS","OS.time","RS"))
    colnames(data)<-c("Age","Gender","Grade","IDH","Chr1p/19q","MGMT",
                      "status","time","RS")
    data$Age<-as.numeric(data$Age)
    data<-data[! is.na(data$Grade),]
    data<-data[! is.na(data$IDH),]
    data<-data[! is.na(data$MGMT),]
    data<-data[! is.na(data$`Chr1p/19q`),]
    data<-data[! is.na(data$RS),]
    
    data$`Risk score`<-ifelse(data[,"RS"]>median(data[,"RS"]),"High","Low")
    data$`Risk score`<-factor(data$`Risk score`,levels=c("Low","High"))
    data$Age<-ifelse(data$Age>50,">50","<50")
    data$Age<-factor(data$Age,levels = c("<50",">50"))
    data$Gender<-factor(data$Gende,levels=c("Female","Male"))
    data[data$Grade=="WHO II","Grade"]<-"G2"
    data[data$Grade=="WHO III","Grade"]<-"G3"
    data[data$Grade=="WHO IV","Grade"]<-"G4"
    data$Grade<-factor(data$Grade,levels=c("G2","G3","G4"))
    data$IDH<-factor(data$IDH,levels=c("Wildtype","Mutant"))
    data[data$MGMT=="un-methylated","MGMT"]<-"Unmethylated"
    data[data$MGMT=="methylated","MGMT"]<-"Methylated"
    data$MGMT<-factor(data$MGMT,levels=c("Unmethylated","Methylated"))
    data$`Chr1p/19q`<-factor(data$`Chr1p/19q`,levels=c("Non-codel","Codel"))
    
    colnames(data)[10] = gene_you_want
    data$MEST = factor(data$MEST, levels =  c("Low","High"))
                       
    ezcox_res <- ezcox(data,
                       covariates = c("MEST"),
                       controls = c("Gender","Age","Grade","IDH","Chr1p/19q","MGMT"),
                       return_models=TRUE)
    
    
    
    tmp<-as.data.frame(ezcox_res[["res"]])
    tmp$Variable<-c('MEST',"Gender","Age","Grade","","IDH","Chr1p/19q","MGMT")
    # tmp$Cohorts<-c(datasets,"","","","","","","")
    tmp$Cohorts<-datasets
    multic_cox<-rbind(multic_cox,tmp)
    multic_cox_model[[datasets]]<-get_models(ezcox_res)
    show_models(multic_cox_model[[datasets]])
  }
  
  
  dt<-multic_cox
  dt$se <- (log(dt$upper_95) - log(dt$HR))/1.96
  dt$` ` <- paste(rep(" ", 20), collapse = " ")
  dt$`HR (95% CI)` <- ifelse(is.na(dt$n_contrast), "",
                             sprintf("%.2f (%.2f - %.2f)",
                                     dt$HR, dt$lower_95, dt$upper_95))
  dt$P <- ifelse(dt$p.value<0.001, "P<0.001",sprintf("%.3f",dt$p.value))
  colnames(dt)[c(3:6,17)]<-c("Contrast","Reference","Number of contrast","Number of reference","P value")
  
  tm <- forest_theme(core=list(bg_params=list(fill =c(rep("#66529F",8),rep("#A37E7D",8),rep("#1c8041",8),rep("#e55709",8)),
                                              alpha =rep(c(0.7,rep(0.5,7)),4))),
                     base_size = 10,
                     # Confidence interval point shape, line type/color/width
                     ci_pch = 16,
                     ci_col = "#3C224B",
                     ci_lty = 1,
                     ci_lwd = 1.5,
                     ci_Theight = 0.2, # Set an T end at the end of CI 
                     # Reference line width/type/color
                     refline_lwd = 1,
                     refline_lty = "dashed",
                     refline_col = "grey20",
                     # Vertical line width/type/color
                     vertline_lwd = 1,
                     vertline_lty = "dashed",
                     vertline_col = "grey20",
                     # Change summary color for filling and borders
                     summary_fill = "#4575b4",
                     summary_col = "#4575b4",
                     # Footnote font size/face/color
                     footnote_cex = 1,
                     footnote_fontface = "italic",
                     footnote_col = "red")
  
  p <-forestploter::forest(dt[,c(13,1,4,3,6,5,15:17)],
                           est = dt$HR,
                           lower = dt$lower_95, 
                           upper = dt$upper_95,
                           sizes = dt$se,
                           #is_summary = c(rep(FALSE, nrow(dt)-2), TRUE,TRUE),
                           ci_column = 7,
                           ref_line = 1,
                           arrow_lab = c("Better", "Worse"),
                           # xlim = c(0, 8),
                           # ticks_at = c(0.5, 1, 2, 5,7.5),
                           x_trans="log2",
                           footnote = " Multivariate Cox Regression",
                           theme = tm)
  p
  p <- add_text(p, text = "Multivariate Cox regression in different cohorts",
                part = "header",
                row = 0,
                col = 4:7,
                just = c("center"),
                gp = gpar(fontface = "bold"))
  
  p <- add_border(p, 
                  part = "header", 
                  row = c(0,1),
                  gp = gpar(lwd = 1))
  
  
  cairo_pdf("./fig/fig_5d_multi_mest.pdf", 
       width = 10.5, height = 8.5,)
  p                
  dev.off()
  
####### 5.6 多个数据集中进行单因素回归并进行meta 分析 #############
  
  res.unicox = lapply(list_train_vali_Data, function(x){
    tmp = x[c('ID','OS.time','OS',gene_you_want)]
    colnames(tmp)[4] = 'RS'
    tmp$Group = ifelse(tmp$RS> median(tmp$RS),'High','Low')
    tmp$Group = factor(tmp$Group, levels = c('Low','High'))
    cox <- coxph(Surv(OS.time, OS) ~ Group, data = tmp)
    
    coxSummary <- summary(cox)
    
    unicox = c( as.numeric(coxSummary$coefficients[,"exp(coef)"])[1],
                as.numeric(coxSummary$coefficients[,"Pr(>|z|)"])[1],
                as.numeric(coxSummary$conf.int[,3][1]),
                as.numeric(coxSummary$conf.int[,4][1])
    )
    
    return(unicox)
  })  
  
  unicox.rs.res = res.unicox  %>% do.call(rbind,.)
  colnames(unicox.rs.res) = c('HR','pvalue',"LCI",'HCI')
  
  library(meta)
  input <- unicox.rs.res %>% as.data.frame() %>%
    mutate(`Hazard Ratio(95%CI)` = paste(HR,'(',LCI,'-',HCI,')',sep=""))
  input$Group = rownames(input)
  # 对不服从正态分布的HR对数转换
  lnhr <- log(input$HR)
  lnlci <- log(input$LCI)
  lnhci <- log(input$HCI)
  selnhr <- (lnhci-lnlci)/(2*1.96)
  metamodel = metagen(TE = lnhr, seTE = selnhr,
                      sm="HR",            # 待合并的效应量
                      data=input,  # 输入数据
                      studlab=Group) 

  library(forestploter)
  library(grid)
  dataset = rownames(unicox.rs.res)
  dt <- metamodel[["data"]]
  dt$Weight_random <- paste(round(100 * metamodel$w.random / sum(metamodel$w.random), 2), "%", sep = "")
  dt$Weight_fixed <- paste(round(100 * metamodel$w.fixed / sum(metamodel$w.fixed), 2), "%", sep = "")
  
  dt[nrow(dt) + 1, ] <- NA
  dt[nrow(dt), 6] <- "Random effect model"
  dt[nrow(dt), 1] <- exp(metamodel$TE.random)
  dt[nrow(dt), 3] <- exp(metamodel$lower.random)
  dt[nrow(dt), 4] <- exp(metamodel$upper.random)
  dt[nrow(dt), 2] <- metamodel$pval.random
  dt[nrow(dt), 10] <- "100%"
  dt[nrow(dt), 11] <- "--"
  
  dt[nrow(dt) + 1, ] <- NA
  dt[nrow(dt), 6] <- "Fixed effect model"
  dt[nrow(dt), 1] <- exp(metamodel$TE.fixed)
  dt[nrow(dt), 3] <- exp(metamodel$lower.fixed)
  dt[nrow(dt), 4] <- exp(metamodel$upper.fixed)
  dt[nrow(dt), 2] <- metamodel$pval.fixed
  dt[nrow(dt), 10] <- "--"
  dt[nrow(dt), 11] <- "100%"
  
  dt$se <- (log(dt$HCI) - log(dt$HR)) / 1.96
  dt$` ` <- paste(rep(" ", 20), collapse = " ")
  dt$`HR (95% CI)` <- ifelse(is.na(dt$Group), "",
                             sprintf(
                               "%.2f (%.2f - %.2f)",
                               dt$HR, dt$LCI, dt$HCI
                             )
  )
  dt$p <- ifelse(dt$pvalue < 0.001, "P<0.001", sprintf("%.3f", dt$pvalue))
  colnames(dt)[c(6:8, 10:11, 15)] <- c("Cohorts", "TE", "SE(TE)", "Weight(random)", "Weight(fixed)", "P value")
  dt$TE <- round(dt$TE, 2)
  dt$`SE(TE)` <- round(dt$`SE(TE)`, 2)
  dt[c((nrow(dt) - 1):nrow(dt)), c(7, 8)] <- ""
  rownames(dt) <- dt$Cohorts
  dt <- dt[c(dataset, "Random effect model", "Fixed effect model"), ]
  dataset_col = paletteer_c("grDevices::Sunset", 30)
  
  
  tm <- forest_theme(
    core = list(bg_params = list(
      fill = c(dataset_col[1:length(dataset)], "grey", "grey"),
      alpha = 0.5
    )),
    base_size = 10,
    # Confidence interval point shape, line type/color/width
    ci_pch = 16,
    ci_col = "#762a83",
    ci_lty = 1,
    ci_lwd = 1.5,
    ci_Theight = 0.2, # Set an T end at the end of CI
    # Reference line width/type/color
    refline_lwd = 1,
    refline_lty = "dashed",
    refline_col = "grey20",
    # Vertical line width/type/color
    vertline_lwd = 1,
    vertline_lty = "dashed",
    vertline_col = "grey20",
    # Change summary color for filling and borders
    summary_fill = "#4575b4",
    summary_col = "#4575b4",
    # Footnote font size/face/color
    footnote_cex = 1,
    footnote_fontface = "italic",
    footnote_col = "red"
  )
  
  p <- forestploter::forest(dt[, c(6:8, 10:11, 13, 14, 15)],
                            est = dt$HR,
                            lower = dt$LCI,
                            upper = dt$HCI,
                            sizes = dt$se,
                            is_summary = c(rep(FALSE, nrow(dt) - 2), TRUE, TRUE),
                            ci_column = 6,
                            ref_line = 1,
                            arrow_lab = c("Better", "Worse"),
                            x_trans = "log2",
                            xlim = c(0, ceiling(max(dt$HCI))),
                            ticks_at = c(0.5, 2^seq(0, floor(log2((ceiling(max(dt$HCI)) - 1))), by = 1)),
                            footnote = " Univariate Cox Regression",
                            theme = tm
  )
  
  p <- add_text(p,
                text = "Meta analysis of univariate Cox regression",
                part = "header",
                row = 0,
                col = 4:6,
                just = c("center"),
                gp = gpar(fontface = "bold")
  )
  
  p <- add_border(p,
                  part = "header",
                  row = c(0, 1),
                  gp = gpar(lwd = 1)
  )
  
  p <- insert_text(p,
                   text = "Meta analysis",
                   row = length(dataset) + 1,
                   just = "left",
                   gp = gpar(cex = 0.8, col = "blue", fontface = "italic")
  )
p  
  
cairo_pdf("./fig/fig_5d_unicox_meta_mest.pdf", 
          width = 9, height = 4.5,)
p                
dev.off()
  


####### 5.7 关于这个基因在 不同队列中1，3，5年AUC的可视化 ##########
  
returnRStoROC <- function(rs.table.list, AUC_time,auc_cal_method) {
  library(survivalROC)
  
  roc.rs <- lapply(rs.table.list, function(x) {
    mySurv <- Surv(x$OS.time, x$OS)

    x$Group <- ifelse(x$RS > median(x$RS), "High", "Low")
    
    
    if (length(unique(x$Group)) > 1) {
      x$Group <- factor(x$Group, levels = c("Low", "High"))
    } else {
      x$Group <- ifelse(x$RS > mean(x$RS), "High", "Low")
    }
    
    if (length(unique(x$Group)) > 1) {
      x$Group <- factor(x$Group, levels = c("Low", "High"))
      data.survdiff <- survdiff(mySurv ~ x$Group)
      HR <- (data.survdiff$obs[2] / data.survdiff$exp[2]) / (data.survdiff$obs[1] / data.survdiff$exp[1])
      
      if (HR > 1) {
        if (auc_cal_method == "NNE") {
          risk.survivalROC <- survivalROC(
            Stime = x$OS.time,
            status = x$OS,
            marker = x$RS,
            predict.time = 365 * AUC_time,
            method = "NNE", span = 0.25 * nrow(x)^(-0.20)
          )
        } else if (auc_cal_method == "KM") {
          risk.survivalROC <- survivalROC(
            Stime = x$OS.time,
            status = x$OS,
            marker = x$RS,
            predict.time = 365 * AUC_time,
            method = "KM"
          )
        } else {
          print("Please provide the correct parameters for method")
        }
      } else {
        if (auc_cal_method == "NNE") {
          risk.survivalROC <- survivalROC(
            Stime = x$OS.time,
            status = x$OS,
            marker = -x$RS,
            predict.time = 365 * AUC_time,
            method = "NNE", span = 0.25 * nrow(x)^(-0.20)
          )
        } else if (auc_cal_method == "KM") {
          risk.survivalROC <- survivalROC(
            Stime = x$OS.time,
            status = x$OS,
            marker = -x$RS,
            predict.time = 365 * AUC_time,
            method = "KM"
          )
        } else {
          print("Please provide the correct parameters for method")
        }
      }
      
      roc_1 <- cbind(round(risk.survivalROC$TP, 3), round(risk.survivalROC$FP, 3))
      
      
      roc_1 <- as.data.frame(roc_1)
      colnames(roc_1) <- c("TP", "FP")
      roc_1$AUC <- risk.survivalROC$AUC
      roc_1$HR <- HR
      return(roc_1)
    } else {
      roc_1 <- data.frame(
        TP = rep(0, nrow(x)),
        FP = rep(1, nrow(x)),
        AUC = rep(0.5, nrow(x)),
        HR = rep(1, nrow(x))
      )
      
      return(roc_1)
    }
  })
  
  
  
  
  
  
  return(roc.rs)
}
  

rs.table.list = lapply(list_train_vali_Data, function(x){
  tmp = x[,c('ID','OS.time','OS',gene_you_want)]
  colnames(tmp)[4] = 'RS'
  return(tmp)
})
names(rs.table.list) = names(list_train_vali_Data)

roc.year = data.frame()  


for (j in 1:10) {
  
  roc = returnRStoROC(rs.table.list = rs.table.list,AUC_time = j,auc_cal_method = 'KM')
  
  for (i in names(roc)) {
    tmp = data.frame(cohort = i, auc = unique(roc[[i]]$AUC), year= j)
    roc.year = rbind(roc.year,tmp)
    
    
  }
}

roc.year$auc[is.na(roc.year$auc)] = 0.9230769

# save(roc.year, file = './res/mest.auc.11gliomaCohorts.Rdata')

ls_AUC = roc.year %>% spread(key  = cohort, value =auc )
rownames(ls_AUC) = paste0(ls_AUC$year,'_y')
ls_AUC = ls_AUC[,-1]
#Heatmap 

hdata <-ls_AUC %>% t() %>% as.data.frame()

hdata <- hdata[order(hdata[,1],decreasing = T) ,]

library(ComplexHeatmap)
library(circlize)
p= Heatmap(as.matrix(hdata), name = "AUC",
        # column_title = NULL,
        # row_title = NULL,
        cluster_rows = F,
        cluster_columns = F,
        # col = colorRamp2(c( 0.5, 0.80), c("#FFFFE3",  "#004628")),
        col = colorRamp2(c( 0.5, 0.80), c("#FFFFE3",  "#09622A")),
        show_row_names = T, show_column_names = T,
        rect_gp = gpar(col = "black", lwd = 2),
        width = nrow(hdata)*unit(6, "mm"),
        height = nrow(hdata)*unit(6, "mm"),
        # na_col = 'white',
        # column_names_side = c('top'),
        # row_split = c(rep('a',5),rep('b',5)),
        column_split = c(rep('a',5),rep('b',5)),
        # top_annotation = column_ha
)
p
cairo_pdf('./fig/fig_5e_auc_gliomacohorts_mest.pdf',width = 6,height = 4)
p
dev.off()

####### 5.8 种类内部不同区域的表达 ########### 


# meta = fread('~/R/ZW/IVYGAP/gene_expression_matrix_2014-11-25/columns-samples.csv')
# df.gene =  fread('~/R/ZW/IVYGAP/gene_expression_matrix_2014-11-25/rows-genes.csv')
# df =fread('~/R/ZW/IVYGAP/gene_expression_matrix_2014-11-25/fpkm_table.csv')
# gene_gtf = df.gene[,c('gene_id','gene_symbol')]
# colnames(df)[1]='gene_id'
# ivgap_expr = left_join(gene_gtf, df, 'gene_id') %>% as.data.frame()
# rownames(ivgap_expr) = ivgap_expr$gene_symbol
# ivgap_expr = ivgap_expr[,-c(1:2)]
# 
# # 对数转化
# ivgap_expr_log=log1p(ivgap_expr)
# save(ivgap_expr_log,plotdata,file = '/export3/zhangw/Glioma.dataset/IVYGAP_processed_expr_meta.Rdata')

# load('/export3/zhangw/Glioma.dataset/IVYGAP_processed_expr_meta.Rdata')

# 
# res.ivygap = GSVA::gsva(as.matrix(ivgap_expr_log),gset.idx.list = genelist,method = 'gsva',kcdf ='Gaussian', parallel.sz=12)
# res.ivygap = res.ivygap %>% t() %>% as.data.frame()
# 
# # save(res.ivygap,file = './res/res.gsva.ivygap.germ.layer.Rdata')
# 
# plotdata = res.ivygap %>% rownames_to_column('rna_well_id')
# plotdata$rna_well_id = as.numeric(plotdata$rna_well_id)
# 
# plotdata = left_join(plotdata, meta, 'rna_well_id')
# 
# plotdata$structure_name[grep( 'Leading',plotdata$structure_name)]= 'Leading edge'
# plotdata$structure_name[grep( 'Infiltrating',plotdata$structure_name)]= 'Infiltrating tumor'
# plotdata$structure_name[grep( 'Cellular ',plotdata$structure_name)]= 'Tumor core'
# plotdata$structure_name[grep( 'Perinecrotic',plotdata$structure_name)]= 'Perinecrotic zone'
# plotdata$structure_name[grep( 'Pseudopalisading',plotdata$structure_name)]= 'Pseudopalisading cells around necrosis'
# plotdata$structure_name[grep( 'Hyperplastic',plotdata$structure_name)]= 'Regions with hyperplastic blood vessels'
# plotdata$structure_name[grep( 'Microvascular',plotdata$structure_name)]= 'Regions with microvascular proliferation'

# save(ivgap_expr_log,plotdata,file = '/export3/zhangw/Glioma.dataset/IVYGAP_processed_expr_meta.Rdata')
your_color_pal=c('#375631','#D3C3D0','#A786BA','#6D4E7E',paletteer_c("grDevices::Sunset", 30)[c(6,12,18,24,30)] )

load('./data/f_germ_layer_genelist.Rdata')
load('/export3/zhangw/Glioma.dataset/IVYGAP_processed_expr_meta.Rdata')
plotdata$structure_name = factor(plotdata$structure_name, levels = c('Leading edge',
                                                                     'Infiltrating tumor',
                                                                     'Tumor core',
                                                                     'Perinecrotic zone',
                                                                     'Pseudopalisading cells around necrosis',
                                                                     'Regions with hyperplastic blood vessels',
                                                                     'Regions with microvascular proliferation'))
df_plot=plotdata %>% 
  dplyr::select("structure_name", "germ_layer") %>% 
  dplyr::rename('value'='germ_layer', 'type'='structure_name') 

mycomparision = list(c('Leading edge','Infiltrating tumor'),
                     c('Infiltrating tumor', 'Tumor core'),
                     c('Leading edge','Tumor core'))
df_plot = as.data.frame(df_plot)
ggplot(df_plot, aes(x=type,y=value))+
  geom_beeswarm(size=2.5,aes(color=type),cex =1)+
  stat_summary(fun.y = "mean", geom = "crossbar", #添加均值
               mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.15,size=0.3,color="black")+
  stat_summary(fun.data = "mean_se", geom="errorbar", width = 0.1,size=0.5,color="black")+#添加误差线
    
  # scale_color_manual(values =   c('#FCDAEC','#85C17E','#F5EF80','#7DA4D2','#AFD2EF','#E8A681', '#B47FB1','#E4C2A8','#BD7477','#1174BC','#29926A','#9E4995')
  scale_color_manual(values =   your_color_pal)+
  scale_y_continuous(limits = c(-0.4,0.5),expand = c(0,0))+
  labs(x="",y="GSVA score")+
  theme_classic()+
  theme(legend.position = "top",
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 12,face = "bold"),
        axis.text.x = element_blank(),
        panel.spacing =unit(0,'cm'))+
  guides(color = guide_legend(nrow = 4))+
  
  ggsignif::geom_signif(comparisons = mycomparision,  position = "identity",  step_increase = 0.1,
                        test = "t.test")

ggsave('./fig/fig_1a_gsva_germ_layer_ivygap.pdf',width = 6,height = 4)


### germ.layer heatmap in ivygap 

hmdat = ivgap_expr_log[f.germ.layer.gene,] %>% na.omit()

# table(is.na(hd.plot))
source('/export3/zhangw/Code.sum/standarize.fun.R')
source('/export3/zhangw/Code.sum/annTrackScale.R')

hmdat_new = hmdat[rowSums(hmdat)>780,]
# is.element('MEST',rownames(hmdat_new))


hd.plot =annTrackScale(indata = hmdat_new,halfwidth = 1)
# range(hd.plot)
annCol = plotdata[,c('rna_well_id','structure_name')]
rownames(annCol)= annCol$rna_well_id
annCol = annCol %>% dplyr::select('structure_name')
colnames(annCol) = 'GBM geometrical regions'
annCol = annCol %>% dplyr::arrange('GBM geometrical regions')

annColors  = list(  'GBM geometrical regions'=c('Leading edge'=your_color_pal[1],
                                                'Infiltrating tumor'= your_color_pal[2],
                                                'Tumor core'=  your_color_pal[3],
                                                'Perinecrotic zone'=  your_color_pal[4],
                                                'Pseudopalisading cells around necrosis'=        your_color_pal[5],
                                                'Regions with hyperplastic blood vessels'=  your_color_pal[6], 
                                                'Regions with microvascular proliferation'=         your_color_pal[7]))

hd.plot = hd.plot %>%  as.data.frame()

annCol$`GBM geometrical regions` = factor(annCol$`GBM geometrical regions`, levels =c('Leading edge',
                                                                                      'Infiltrating tumor',
                                                                                      'Tumor core',
                                                                                      'Perinecrotic zone',
                                                                                      'Pseudopalisading cells around necrosis',
                                                                                      'Regions with hyperplastic blood vessels',
                                                                                      'Regions with microvascular proliferation') )

annCol = annCol %>% dplyr::arrange(annCol$`GBM geometrical regions`)



library(ComplexHeatmap)

hm1 <- pheatmap(as.matrix(hd.plot[,rownames(annCol)]),
                border_color = NA, # 无边框
                color = rev(paletteer_c("grDevices::ag_Sunset", 128)),
                show_rownames = T, # 显示行名
                show_colnames = F, # 不显示列名
                cluster_rows = T, # 行不聚类
                cluster_cols = F, # 列不聚类
                name = "Gene\nZ-score", # 颜色图例的名字（所有热图统一，这样只会显示一个）
                cellheight = 9, # 热图元素高度400
                cellwidth =0.4, # 热图总宽度固定为300
                gaps_col = cumsum(table(annCol$`GBM geometrical regions`)), # 切割类
                annotation_col = annCol, # 样本注释
                # annotation_row = annRow,
                annotation_colors = annColors,
                legend = TRUE,
                annotation_names_row = T,
                treeheight_row = 0
                # gaps_row = cumsum(rev(table(annRow$Tpye))) # 行分割
                
) # 注释颜色
hm1

pdf('./fig/fig_1a_ivygap_layer_germ_genelist.pdf',width =10 ,height = 10,onefile = F)
print(hm1)
dev.off()



# ### MEST 

df.new = ivgap_expr_log['MEST',] %>% t() %>% as.data.frame() %>% rownames_to_column('rna_well_id')
df.new$rna_well_id = as.numeric(df.new$rna_well_id)
df.new.plot = left_join(df.new,plotdata[,c('rna_well_id','structure_name')],'rna_well_id')

df_plot=df.new.plot %>% 
  dplyr::select("structure_name", "MEST") %>% 
  dplyr::rename('value'='MEST', 'type'='structure_name') 

mycomparision = list(c('Leading edge','Infiltrating tumor'),
                     c('Infiltrating tumor', 'Tumor core'),
                     c('Leading edge','Tumor core'))
df_plot = as.data.frame(df_plot)
ggplot(df_plot, aes(x=type,y=value))+
  geom_beeswarm(size=2.5,aes(color=type),cex =1)+
  stat_summary(fun.y = "mean", geom = "crossbar", #添加均值
               mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.15,size=0.3,color="black")+
  stat_summary(fun.data = "mean_se", geom="errorbar", width = 0.1,size=0.5,color="black")+#添加误差线
  scale_color_manual(values =   your_color_pal)+
  
  # scale_color_manual(values =   c('#FCDAEC','#85C17E','#F5EF80','#7DA4D2','#AFD2EF','#E8A681', '#B47FB1','#E4C2A8','#BD7477','#1174BC','#29926A','#9E4995')
  # )+
  scale_y_continuous(limits = c(0,8),expand = c(0,0))+
  labs(x="",y="log2 of MEST")+
  theme_classic()+
  theme(legend.position = "top",
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 12,face = "bold"),
        axis.text.x = element_blank(),
        panel.spacing =unit(0,'cm'))+
  guides(color = guide_legend(nrow = 4))+
  ggsignif::geom_signif(comparisons = mycomparision,  position = "identity",  step_increase = 0.1,
                        test = "t.test")

ggsave('./fig/fig_5f_mest_ivygap.pdf',width = 6,height = 4)

####### 5.9 MEST 的相关性分析和信号通路的富集分析#############

# 计算不同队列中的差异基因 
rm(list = ls())
gc()
source('./code/MEST_your_col.R')
load("/export3/zhangw/Project_Cross/Project_Mime/Proj/data/glioma.cohort.Rdata")


source('/export3/zhangw/Project_Cross/Project.Glioma.SGene/code/analysis.cor.single.gene.glioma.R')

res.mest.cor = lapply(names(list_train_vali_Data), function(x){
  
  res.cor  = analysis.cor.single.gene.glioma(gene_you_want ='MEST',method = 'pearson',cancer_type = 'glioma' ,datasets = x)
  return(res.cor)
})
names(res.mest.cor) = names(list_train_vali_Data)

# save(res.mest.cor,file = './res/res.cor.mest.glioma.Rdata')
load('./res/res.cor.mest.glioma.Rdata')

df = res.mest.cor %>% do.call(rbind,.)

df = df[abs(df$cor_r)>0.1,]
#添加显著性标签：
df$label <- ifelse(df$adjP<0.01&abs(df$cor_r)>0.1,"adjust P-val<0.01","adjust P-val>=0.01")
head(df)
unique(df$dataset.sum)

load("/export3/zhangw/Project_Cross/Project_CNGA3/data/protein.coding.genes.Rdata")

# 只考虑 protein coding 的基因
df = df[df$gene_name2 %in% protien.coding.gene,]
df = na.omit(df)



#获取每个队列中相关性最显著的10个基因；
top10sig0 <- filter(df,dataset.sum=="TCGA") %>% distinct(gene_name2,.keep_all = T) %>% top_n(10,abs(cor_r))
top10sig1 <- filter(df,dataset.sum=="CGGA.325") %>% distinct(gene_name2,.keep_all = T) %>% top_n(10,abs(cor_r))
top10sig2 <- filter(df,dataset.sum=="CGGA.693") %>% distinct(gene_name2,.keep_all = T) %>% top_n(10,abs(cor_r))
top10sig3 <- filter(df,dataset.sum=="CGGA.1018") %>% distinct(gene_name2,.keep_all = T) %>% top_n(10,abs(cor_r))
top10sig4 <- filter(df,dataset.sum=="CGGA.array") %>% distinct(gene_name2,.keep_all = T) %>% top_n(10,abs(cor_r))
top10sig5 <- filter(df,dataset.sum=="GLASS_TP") %>% distinct(gene_name2,.keep_all = T) %>% top_n(10,abs(cor_r))
top10sig6 <- filter(df,dataset.sum=="GLASS_R1") %>% distinct(gene_name2,.keep_all = T) %>% top_n(10,abs(cor_r))
top10sig7 <- filter(df,dataset.sum=="GSE108474") %>% distinct(gene_name2,.keep_all = T) %>% top_n(10,abs(cor_r))
top10sig8 <- filter(df,dataset.sum=="GSE16011") %>% distinct(gene_name2,.keep_all = T) %>% top_n(10,abs(cor_r))
top10sig9 <- filter(df,dataset.sum=="GSE43289") %>% distinct(gene_name2,.keep_all = T) %>% top_n(10,abs(cor_r))
top10sig10 <- filter(df,dataset.sum=="GSE7696") %>% distinct(gene_name2,.keep_all = T) %>% top_n(10,abs(cor_r))

#将提取所有cluster的Top10基因表格合并：
top10sig <- rbind(top10sig0,top10sig1,top10sig2,top10sig3,top10sig4,top10sig5,top10sig6,top10sig7,top10sig8,top10sig9,top10sig10)

#新增一列，将Top10的差异基因标记为2，其他的标记为1；
df$size <- case_when(!(df$gene_name2 %in% top10sig$gene_name2)~ 1,
                     df$gene_name2 %in% top10sig$gene_name2 ~ 2)

#提取非Top10的基因表格；
dt <- filter(df,size==1)
head(dt)

#绘制每个Cluster Top10以外基因的散点火山图：
#叠加每个Cluster Top10基因散点(将散点适当放大强调）：
#根据图p中log2FC区间确定背景柱长度：
dfbar<-data.frame(x= names(list_train_vali_Data),
                  y=c(0.75,0.7,0.75,0.77,0.4,0.65,0.70,0.55,0.5,0.68,0.5))
dfbar1<-data.frame(x=names(list_train_vali_Data),
                   y=c(-0.4,-0.55,-0.35,-0.4,-0.45,-0.45,-0.38,-0.6,-0.4,-0.7,-0.5))

#把散点火山图叠加到背景柱上：
#相关包的载入
library(ggplot2)
library(tidyverse)
library(ggrepel)
#把散点火山图叠加到背景柱上：

#添加X轴的cluster色块标签：
dfcol<-data.frame(x=names(list_train_vali_Data),
                  y=0,
                  label=names(list_train_vali_Data))
mycol <- paletteer_d("ggthemes::stata_s2color",n = 11)

ggplot()+
  geom_col(data = dfbar,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)+
  geom_col(data = dfbar1,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)+
  geom_jitter(data = dt,
              aes(x = dataset.sum, y = cor_r, color = label),
              size = 0.85,
              width =0.4)+
  geom_jitter(data = top10sig,
              aes(x = dataset.sum, y = cor_r, color = label),
              size = 1,
              width =0.4)+
  scale_color_manual(name=NULL,
                     values = c('#501d8a','#0d5b26'))+
 geom_tile(data = dfcol,
                     aes(x=x,y=y),
                     height=0.15,
                     color = "black",
                     fill = mycol,
                     alpha = 0.6,
                     show.legend = F)+
  #给每个Cluster差异表达前Top10基因加上标签：
  
  geom_text_repel(
    data=top10sig,
    aes(x=dataset.sum,y=cor_r,label=gene_name2),
    force = 1.2,
    arrow = arrow(length = unit(0.008, "npc"),
                  type = "open", ends = "last")
  )+
  #修改X/Y轴标题和添加cluster数字：
  
  labs(x="Cohorts",y="Pearson's R",title = 'Correlation with MEST')+
  geom_text(data=dfcol,
            aes(x=x,y=y,label=label),
            size =6,
            color ="white")+
  #自定义主题美化：
  
  theme_minimal()+
  theme(
    axis.title = element_text(size = 13,
                              color = "black",
                              face = "bold"),
    axis.line.y = element_line(color = "black",
                               size = 1.2),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "top",
    legend.direction = "vertical",
    legend.justification = c(1,0),
    legend.text = element_text(size = 12)
  )

ggsave('./fig/fig_5g_glioma_cor_mest_gene.pdf',width = 16,height = 6)



donut.dp = data.frame( label=NULL,Percentage= NULL,cohrot= NULL)

for (i in unique(df$dataset.sum)) {
  tmp = df[df$dataset.sum==i,] %>% count('label')
  tmp$freq = tmp$freq/sum(tmp$freq)
  colnames(tmp)[2] = 'Percentage'
  tmp$cohrot = i
  donut.dp = rbind(donut.dp,tmp)
  print(tmp)
}

donut.dp

ggplot(donut.dp, aes(x = 3, #给一个数值向量作为柱子x轴的中心坐标
               y = Percentage,
               fill = label)) + #将颜色映射到celltype
  geom_col(width = 1.5, #柱形图柱子宽度
           color = 'white') + #描边颜色
  # facet_grid(.~cohrot )+ #按tissue分面
  facet_wrap(vars(cohrot), nrow = 1)+
  coord_polar(theta = "y")+#饼图：直角坐标系转换为极坐标系；

   xlim(c(0.2, 3.8))+#甜甜圈图：将饼图中心“挖空”；
  scale_fill_manual(values = c('#501d8a','#0d5b26')) +
  theme_void()+ #空白主题
  theme(
    strip.text.x = element_text(size = 14), #分面标签大小
    legend.title = element_text(size = 15), #图例标题大小调整
    legend.text = element_text(size = 14), #图例标签大小调整
    legend.position = 'top'
  )+#添加文本标签：
  geom_text(aes(label = paste0(round(Percentage,2),'%')),
            position = position_stack(vjust = 0.5),color='black',
            size = 4)

ggsave('./fig/fig_5g_glioma_cor_mest_gene_per.pdf',width = 16,height = 3)


#### 提取多个数据集合中相关的基因

postive.gene = df[df$adjP < 0.01 &df$cor_r > 0.1,]$gene_name2 %>% count()
negative.gene = df[df$adjP < 0.01 &df$cor_r <  -0.1,]$gene_name2 %>% count()
f.posi.gene.5 = postive.gene[postive.gene$freq>6,]$x
f.nega.gene.5 = negative.gene[negative.gene$freq>3,]$x

DEG<-data.frame(gene=c(f.posi.gene.5,f.nega.gene.5),
                group=c(rep("MEST Postive-related",length(f.posi.gene.5)),
                        rep("MEST Negative-related",length(f.nega.gene.5))))
ids <- bitr(DEG$gene, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db",drop = TRUE)
DEG_list=merge(DEG,ids,by.x='gene',by.y='SYMBOL')


library(clusterProfiler,lib.loc = "/export/bioinfo-team/home/liuhw/R/x86_64-pc-linux-gnu-library/4.1")
enrichGO <- compareCluster(ENTREZID~group, data=DEG_list, fun ="enrichGO",ont="All",OrgDb='org.Hs.eg.db') 
enrichGO<-setReadable(enrichGO,OrgDb='org.Hs.eg.db')
# enrichKEGG <- clusterProfiler::compareCluster(ENTREZID~group, data=DEG_list, fun ="enrichKEGG",organism ="hsa")
# enrichKEGG<-setReadable(enrichKEGG,OrgDb='org.Hs.eg.db',keyType="ENTREZID")

# saveRDS(enrichKEGG,file = "./res/mest_glioma_enrichKEGG.rds")
saveRDS(enrichGO,file = "./res/mest_glioma_enrichGO.rds")

# tmp<-as.data.frame(enrichKEGG)
# select_enrichment_kegg<-tmp[c(grep("Cell cycle",tmp$Description),grep("Cellular senescence",tmp$Description),
#                               grep("FoxO signaling pathway",tmp$Description),grep("NOD-like receptor signaling pathway",tmp$Description),
#                               grep("Necroptosis",tmp$Description),grep("TNF signaling pathway",tmp$Description),
#                               grep("MAPK signaling pathway",tmp$Description)),]
# select_enrichment_kegg$type<-"KEGG"

tmp<-as.data.frame(enrichGO)
select_enrichment_go<-tmp[c(grep("cell cycle checkpoint signaling",tmp$Description),
                            grep("DNA replication",tmp$Description),
                            grep("DNA damage checkpoint signaling",tmp$Description),
                            grep("regulation of DNA recombination",tmp$Description),
                            grep("establishment of spindle localization",tmp$Description),
                            grep("cell growth",tmp$Description),
                            grep("GO:0048638",tmp$ID),
                            grep("GO:0046621",tmp$ID),
                            grep("GO:0048285",tmp$ID),
                            grep("mRNA transport",tmp$Description),
                            grep("mitotic spindle organization",tmp$Description)),]
select_enrichment_go<-select_enrichment_go[-c(grep("negative ",select_enrichment_go$Description),
                                              grep("positive",select_enrichment_go$Description),
                                              grep("p53",select_enrichment_go$Description)),]
select_enrichment_go$type<-"GO"
select_enrichment<-select_enrichment_go[,-3]
# select_enrichment<-rbind(select_enrichment_kegg,select_enrichment_go[,-3])
select_enrichment$Description<-factor(select_enrichment$Description,levels = c(unique(select_enrichment$Description)))
select_enrichment$type<-factor(select_enrichment$type,levels = c("GO"))

p<-ggplot(select_enrichment, aes(x=group,y=Description,size=Count,fill=qvalue))  + 
  geom_point(shape = 21) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major=element_blank(),
        panel.border = element_rect(colour = "black", fill=NA,size = 1),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        legend.key = element_rect(fill = "white", colour = NA))+
  scale_fill_gradientn(colours=paletteer::paletteer_c("ggthemes::Classic Red-White-Black", 30))+
  labs(title="Enrichment", x="", y="")+
  coord_cartesian(clip = 'off')

p

pdf("./fig/fig_5h_glioma_mest_go.pdf",bg="white",width=5.5,height=7.7,pointsize=6)
# annoRect(object = p,
#          annoPos = 'left',
#          aesGroup = T,
#          aesGroName = 'type',
#          xPosition = c(-12,0.4),
#          alpha = 0.3,
#          rectWidth = 1,
#          addText = T,
#          textRot = 90,
#          textSize = 11,
#          continuesRect = T,
#          revColH = T,
#          #border = T,
#          revColV = T,
#          pFill=c("#6BAED6","#D9F0A3"),
#          textCol = c("#1F78B4","#33A02C"),
#          textHVjust = -5.5)
p
dev.off()


# 在TCGA 中进行gsea 分析
tmp = list_train_vali_Data$TCGA

treat.id = data.frame(ID=tmp[tmp$MEST > median(tmp$MEST),]$ID)
contr.id = data.frame(ID=tmp[!tmp$MEST > median(tmp$MEST),]$ID)
source('/export3/zhangw/Code.sum/two_subtype_compr_logfc.R')
exp = tmp[,-c(1:3)]
rownames(exp) = tmp$ID
exp = exp %>% t() %>% as.data.frame()
rees.deg = two_subtype_logFC_Cpr_wilcox(expr = exp,treat_list = treat.id,ctrl_list = contr.id,method = 'w')

diff_genes = rees.deg %>% rownames_to_column('id') %>% na.omit()
# 将列名从 "id" 更改为 "SYMBOL"
names(diff_genes)[names(diff_genes) == "id"] <- "SYMBOL"
colnames(diff_genes)
# 选择差异基因数据中的"SYMBOL"和"logFC"两列
diff_genes <- diff_genes[, c("SYMBOL", "log2FoldChange")]
colnames(diff_genes)[2] = 'logFC'

diff_genes= diff_genes[!diff_genes$logFC%in%c(Inf,-Inf),]

# 使用bitr函数将SYMBOL转换为ENTREZID
df_id <- bitr(
  diff_genes$SYMBOL,  # 输入的向量，这里是SYMBOL列
  fromType = "SYMBOL",      # 输入的ID类型，这里是SYMBOL
  toType = "ENTREZID",      # 想要转换成的ID类型，这里是ENTREZID
  OrgDb = "org.Hs.eg.db"    # 使用的生物数据库，这里是人类的基因组数据库
)


# 使用merge合并数据
df_all <- merge(diff_genes, df_id, by = "SYMBOL", all = FALSE)

# 排序合并后的数据框
df_all_sort <- df_all[order(df_all$logFC, decreasing = TRUE), ]

# 提取foldchange
gene_fc <- df_all_sort$logFC 


# 将gene_fc对象的命名改为df_all_sort数据框中ENTREZID列的值
names(gene_fc) <- df_all_sort$ENTREZID


# # 运行GSEA
gseKEGG <- gseKEGG(gene_fc, organism = "hsa", pvalueCutoff = 1)
# 
# # 将结果转换为数据框
# KEGG_table <- as.data.frame(gseKEGG)
# 
# # 打印列名
# colnames(KEGG_table)


ego3 <- gseGO(geneList     = gene_fc,
              OrgDb        = 'org.Hs.eg.db',
              ont          = "BP",
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)
ego3_table <- as.data.frame(ego3)

.libPaths()
library(GseaVis,lib.loc =  "/export/bioinfo-team/home/liuhw/R/x86_64-pc-linux-gnu-library/4.1")

geneSetID = c('GO:0006338',
              'GO:0006323',
              'GO:0031497','GO:0006333'
              )

# all plot
gseaNb(object = ego3,
       geneSetID = geneSetID,
       subPlot = 2,         
       pCol = 'black',        
       pHjust = 0,
       rmSegment=T,

       termWidth = 35,
       legend.position = c(0.8,0.8),addPval = T,
       pvalX = 0.05,pvalY = 0.05, curveCol = c('#3C224B','#A786BA','#0d5b26','#96D2B0'),
       rankCol = c("#26456E", "#FFF4EE", "#7B3014")
         )
ggsave('./fig/fig_5h_gse_mest_Tcga_glioma.pdf',width = 5,height = 4)





########## 计算MEST和各个通路之间的相关性 
# rm(list = ls())
gc()

ffiles = list.files('./data/genelist_egfr/')

egfr.gene = lapply(ffiles, function(x){
  a = read.gmt(paste0('./data/genelist_egfr/',x))
  return(a)
}) %>% do.call(rbind,.)
egfr.gene = unique(egfr.gene$gene)

source('/export3/zhangw/Project_Cross/Project.Glioma.SGene/code/analysis.cor.genes.enrich.pathway.glioma.R')
load("/export3/zhangw/Project_Cross/Project_Mime/Proj/data/glioma.cohort.Rdata")
res.cor = lapply(names(list_train_vali_Data), function(x){
  
  tmp = analysis.cor.genes.enrich.pathway.glioma(input_genelist = list(EGFR=egfr.gene),
                                                     method_cor = 'pearson',
                                                     method_enrich = 'gsva',
                                                     kcdf = 'Gaussian',parallel.sz = 24,cancer_type = 'glioma',datasets = x)
  return(tmp)

})

names(res.cor) = names(list_train_vali_Data)
# save(res.cor,file = './res/res.cor.egfr.all.gene.glioma.Rdata')

# 将 MEST 提取出来
cor.mest = lapply(res.cor, function(x){
  tmp = x
  tmpt = tmp[tmp$gene_name2 =='MEST',]
  return(tmpt)
}) %>% do.call(rbind,.)

dfplot =cor.mest %>% rownames_to_column('cohort')

plist = list()

dfplot <- data.frame(cohort = cor.mest$dataset.sum,
                       r = cor.mest$cor_r,
                       p = -log10(cor.mest$pvalue),
                       Significance = ifelse(cor.mest$pvalue<0.05,'Sig','NotSig'))

p1 <- ggplot(data = dfplot,aes(r,forcats::fct_reorder(cohort,r,.desc = T))) +
  geom_segment(aes(xend=0,yend=cohort,colour=Significance ),linetype = 2) +
  geom_point(aes(size=p),col = alpha(c('#501d8a'),0.5)) +
  scale_size_continuous(range =c(2,8)) +
  scale_x_continuous(breaks = c(-0.2,0, 0.25, 0.5),
                  expand = expansion(mult = c(0.01,.1))) + #左右留空
  theme_classic() +
  labs(x = "Correlation coefficient", y = "", size = bquote("-log"[10]~"("~italic(P)~"-value)")) + 
  theme(legend.position = "bottom", 
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank())+
  scale_color_manual(values =c('#9E9E9E','#ad2e28'))  # 设置线段颜色为蓝色和红色

p2=p1+xlab(label = '')+ggtitle('MEST versus EGFR signal',subtitle = 'Pearson R with GSVA score')+
  guides(color=guide_legend(nrow = 1),
         size = guide_legend(nrow = 1))

p2
ggsave('./fig/fig_5h_mest_egfr_gsva_cor.pdf',width = 4,height = 6)

plist[[1]] = p2


####### ##################所有队列中计算EGFR 和 MEST 的相关性####

res.core_egfr = lapply(names(list_train_vali_Data), function(x){
  tmp = list_train_vali_Data[[x]]
  ptest = cor.test(tmp$MEST, tmp$EGFR,alternative = 'two.sided',method = 'pearson')
  innnerdat = data.frame(cohort =x, cor = ptest$estimate, p = ptest$p.value )
  return(innnerdat)
}) %>% do.call(rbind,.)


dfplot <- data.frame(cohort = res.core_egfr$cohort,
                     r = res.core_egfr$cor,
                     p = -log10(res.core_egfr$p),
                     Significance = ifelse(res.core_egfr$p<0.05,'Sig','NotSig'))

p1 <- ggplot(data = dfplot,aes(r,forcats::fct_reorder(cohort,r,.desc = T))) +
  geom_segment(aes(xend=0,yend=cohort,colour=Significance ),linetype = 2) +
  geom_point(aes(size=p),col = alpha(c('#501d8a'),0.5)) +
  scale_size_continuous(range =c(2,8)) +
  scale_x_continuous(breaks = c(-0.2,0, 0.25, 0.5),
                     expand = expansion(mult = c(0.01,.1))) + #左右留空
  theme_classic() +
  labs(x = "Correlation coefficient", y = "", size = bquote("-log"[10]~"("~italic(P)~"-value)")) + 
  theme(legend.position = "bottom", 
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank())+
  scale_color_manual(values =c('#9E9E9E','#ad2e28'))  # 设置线段颜色为蓝色和红色

p2=p1+xlab(label = '')+ggtitle('MEST versus EGFR ',subtitle = 'Pearson R')+
  guides(color=guide_legend(nrow = 1),
         size = guide_legend(nrow = 1))

p2
ggsave('./fig/fig_5h_mest_egfr_mRNA_cor.pdf',width = 4,height = 6)

plist[[2]] = p2 

pdf('./fig/fig_5h_mest_egfr_mRNA_and_gsva_cor.pdf',width = 8,height = 6)
plist[[1]]+plist[[2]]+plot_layout(nrow =1,widths = c(1,1),guides='collect' )&
  theme(legend.position='bottom')
dev.off()


#### 在每一个队列中进行相关性可视化

corplot <- vector("list",11) 
red = alpha(c('#3C224B'),0.7)
lightred=c('#F2ECE3')
special.gene = c('MEST','EGFR')
for (i in c(1:11)) {
  cohort.name = names(list_train_vali_Data)[i]
  tmp.g = list_train_vali_Data[[cohort.name]] %>% dplyr::select('ID',special.gene)
  colnames(tmp.g) = c('ID','var1','var2')
  
  p = ggplot(tmp.g,aes(var1,var2))+
    geom_ribbon(stat = "smooth", method = "lm", se = TRUE, # 先画置信区间的彩带以免遮挡散点
                fill = alpha(lightred, 0.5)) +
    geom_smooth(span = 2, method = lm, color = red, fill = NA) + # 绘制回归线
    geom_point(color = red, size = 2) +
    # ggpubr::stat_cor(label.x = ifelse(cohort.name =='CGGA.array',-8,0.1) , color = 'black')+ ## 
    ggpubr::stat_cor(label.x = min(tmp.g$var1) , color = 'black')+ ## 
    # theme_classic()+
    theme_bw() +
    theme(axis.ticks = element_line(size = 0.2, color = "black"),
          axis.ticks.length = unit(0.2, "cm"),
          axis.text = element_text(size = 10, color = "black"),
          axis.title = element_text(size = 10, color = "black"),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())+
    xlab(paste0('log2 of MEST'))+ylab("log2 of EGFR")+ggtitle(cohort.name)
  
  print(corplot[[i]] <-  p)
}

pdf('./fig/fig_5h_individule_cor_mest_egfr.pdf',width = 9,height = 12,onefile = F)
corplot[[1]]+corplot[[2]]+corplot[[3]]+corplot[[4]]+
  corplot[[5]]+corplot[[6]]+corplot[[7]]+corplot[[8]]+
  corplot[[9]]+
  corplot[[10]]+
  corplot[[11]]+
  
  plot_layout(ncol = 3)
dev.off()

####### ##################在胶质瘤细胞系中计算MEST和EGFR以及EGFR通路的相关性质 ##########
load("/export3/zhangw/Project_Cross/Project_GliomaRSI/data/cellline_glioma_tpm_expr.Rdata")
tmp.g = df.cal[,c('MEST','EGFR')]
colnames(tmp.g) = c('var1','var2')
tmp.g = apply(tmp.g, 2 , log1p)
tmp.g = as.data.frame(tmp.g)
red = alpha(c('#3C224B'),0.7)
lightred=c('#F2ECE3')
p1 = ggplot(tmp.g,aes(var1,var2))+
  geom_ribbon(stat = "smooth", method = "lm", se = TRUE, # 先画置信区间的彩带以免遮挡散点
              fill = alpha(lightred, 0.5)) +
  geom_smooth(span = 2, method = lm, color = red, fill = NA) + # 绘制回归线
  geom_point(color = red, size = 2) +
  # ggpubr::stat_cor(label.x = ifelse(cohort.name =='CGGA.array',-8,0.1) , color = 'black')+ ## 
  ggpubr::stat_cor(label.x = min(tmp.g$var1) , color = 'black',  method = "spearman",cor.coef.name='rho'
  )+ ## 
  # theme_classic()+
  theme_bw() +
  theme(axis.ticks = element_line(size = 0.2, color = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 10, color = "black"),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab(paste0('log2 of MEST'))+ylab("log2 of EGFR")+ggtitle('CCEL')

p1

# egfr signal 
ffiles = list.files('./data/genelist_egfr/')

egfr.gene = lapply(ffiles, function(x){
  a = read.gmt(paste0('./data/genelist_egfr/',x))
  return(a)
}) %>% do.call(rbind,.)
egfr.gene = unique(egfr.gene$gene)
res.gsva = gsva(as.matrix(log1p(t(df.cal))),list(EGFR =egfr.gene),method = 'ssgsea',
                kcdf = 'Gaussian',parallel.sz = 24)

df.plot =as.data.frame(t(res.gsva)) %>% cbind(df.cal[,'MEST'])
colnames(df.plot) = c('var1','var2')
df.plot$var2 = log1p(df.plot$var2)
df.plot = df.plot[,c('var2','var1')]
colnames(df.plot) = c('var1','var2')

p2 = ggplot(df.plot,aes(var1,var2))+
  geom_ribbon(stat = "smooth", method = "lm", se = TRUE, # 先画置信区间的彩带以免遮挡散点
              fill = alpha(lightred, 0.5)) +
  geom_smooth(span = 2, method = lm, color = red, fill = NA) + # 绘制回归线
  geom_point(color = red, size = 2) +
  # ggpubr::stat_cor(label.x = ifelse(cohort.name =='CGGA.array',-8,0.1) , color = 'black')+ ## 
  ggpubr::stat_cor(label.x = min(tmp.g$var1) , color = 'black',  method = "spearman",cor.coef.name='rho'
                   # ggpubr::stat_cor(label.x = min(tmp.g$var1) , color = 'black',  method = "pearson",cor.coef.name='rho'
  )+ ## 
  # theme_classic()+
  theme_bw() +
  theme(axis.ticks = element_line(size = 0.2, color = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 10, color = "black"),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab(paste0('log2 of MEST'))+ylab("GSVA score of EGFR pathway")+ggtitle('CCEL')

p2

pdf('./fig/fig_5h_ccle_cor_mest_egfr.pdf',width = 6,height = 3,onefile = F)
p1+p2+
  
  plot_layout(ncol = 2)
dev.off()









####### 6.实验数据处理 ############
# rm(list = ls())
gc()

exper_fig = list()

####### 6.1小鼠生存  ############

source('./code/MEST_your_col.R')
df = read_xlsx('./data/experi/survival-MEST.xlsx')
colnames(df) = c('OS.time','Group','OS')

df$Group = factor(df$Group,levels = unique(df$Group))
fit<-survfit(Surv(OS.time,OS)~Group,data = df)


p=ggsurvplot(fit,df,
           pval = TRUE,
           surv.median.line = "hv",
           xlab="Days after implantation",
           ylab="U87 MG \n Percent survival (%)",
           xlim=c(0,45),
           risk.table = F,
           legend.title="",
           legend.labs=c("Vehicle (n= 8)","shMEST #1 (n=8)"),
           palette = doubel_col,
           axes.offset=FALSE)+
  guides(color = guide_legend(nrow = 1))


p
ggsave('./fig/fig_8a_mouse_survival.pdf',width = 4,height = 4)

exper_fig[['mouse_sur']]=p

####### 6.2 ki67  ##############
library(readxl)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(zoo)
# install.packages('numform')
library(numform)#对数字和图表进行统一格式化
# rm(list = ls())
# gc()
exper_fig = list()
source('./code/MEST_your_col.R')

colnames(df)
df= read_xlsx('./data/experi/ki67.xlsx')
df$celltype[df$celltype=='u87']='U87'
df[df$celltype=='U251',c('control', "sh-MEST#1" ,"sh-MEST#2")] = na.aggregate(df[df$celltype=='U251',c('control', "sh-MEST#1" ,"sh-MEST#2")],2,mean)
df[df$celltype=='U87',c('control', "sh-MEST#1" ,"sh-MEST#2")] = na.aggregate(df[df$celltype=='U87',c('control', "sh-MEST#1" ,"sh-MEST#2")],2,mean)

plist = list()
for (i in c('U251','U87')) {
  
  cellname = i
  # 宽数据转换为长数据格式
  df%>% dplyr::filter(df$celltype==cellname,) %>% 
    pivot_longer(cols = 2:4,
                 names_to = "Group",
                 values_to ="Value")->Figure_2b
  # 分组标签
  Group<-unique(Figure_2b$Group)
  # 分组因子化
  Figure_2b$Group<-factor(Figure_2b$Group,levels = unique(Figure_2b$Group),labels = c('control', "sh-MEST#1" ,"sh-MEST#2"))
  my_comparisons <- list(c("control","sh-MEST#1"),c("control","sh-MEST#2"))
  
  
  p1=ggplot(Figure_2b, aes(Group, Value))+ 
    # 添加散点图
    geom_jitter(aes(x=Group, Value,colour = Group), size = 2.5, width = 0.1, height = 0)+
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +     # 误差棒
    stat_summary(fun = mean, geom = "crossbar", width = 0.5, size = 0.2)+  # 指代均数的水平横线
    theme_classic()+scale_color_manual(values =c('#375631','#D3C3D0','#3C224B') )+
    
    stat_compare_means(comparisons=my_comparisons,method = "t.test",hide.ns = F,label = "p.format")+
    labs(title = cellname)+
    ylab("Intensity of Ki67/DAPI") +
    xlab("") +
    #ylim(35,37.5)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
          legend.position = "none",
          axis.ticks = element_line(size=0.2, color="black"),
          axis.ticks.length = unit(0.2, "cm"),
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 10),
          plot.title = element_text(hjust = 0.5))
  
  
  plist[[i]] = p1
  
}

pdf('./fig/fig_8b_u251_87_ki67_dapi_if.pdf',width = 5,height = 4.8,onefile = F)
plist[[1]]+plist[[2]]
dev.off()
exper_fig[['ki67_u251']] = plist[[1]]
exper_fig[['ki67_u87']] = plist[[2]]

####### 6.3  transwell   ##############
# rm(list = ls())
library(readxl)
library(tidyverse)
library(zoo)
library(ggpubr)
library(patchwork)
df = read_xlsx('./data/experi/transwell.xlsx')
df[df$celltype=='U251',c('control', "sh-MEST#1" ,"sh-MEST#2")] = na.aggregate(df[df$celltype=='U251',c('control', "sh-MEST#1" ,"sh-MEST#2")],2,mean)
df[df$celltype=='U87',c('control', "sh-MEST#1" ,"sh-MEST#2")] = na.aggregate(df[df$celltype=='U87',c('control', "sh-MEST#1" ,"sh-MEST#2")],2,mean)


plist = list()
for (i in c('U251','U87')) {
  
  cellname = i
  # 宽数据转换为长数据格式
  df%>% dplyr::filter(df$celltype==cellname,) %>% 
    pivot_longer(cols = 2:4,
                 names_to = "Group",
                 values_to ="Value")->Figure_2b
  # 分组标签
  Group<-unique(Figure_2b$Group)
  # 分组因子化
  Figure_2b$Group<-factor(Figure_2b$Group,levels = unique(Figure_2b$Group),labels = c('control', "sh-MEST#1" ,"sh-MEST#2"))
  my_comparisons <- list(c("control","sh-MEST#1"),c("control","sh-MEST#2"))
  
  
  p1=ggplot(Figure_2b, aes(Group, Value))+ 
    # 添加散点图
    geom_jitter(aes(x=Group, Value,colour = Group), size = 2.5, width = 0.1, height = 0)+
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +     # 误差棒
    stat_summary(fun = mean, geom = "crossbar", width = 0.5, size = 0.2)+  # 指代均数的水平横线
    theme_classic()+scale_color_manual(values =c('#375631','#D3C3D0','#3C224B') )+
    
    stat_compare_means(comparisons=my_comparisons,method = "t.test",hide.ns = F,label = "p.format")+
    labs(title = cellname)+
    ylab("Migration cells per field") +
    xlab("") +
    #ylim(35,37.5)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          # panel.border = element_rect(colour = "black", fill=NA, size=0.2),
          legend.position = "none",
          axis.ticks = element_line(size=0.2, color="black"),
          axis.ticks.length = unit(0.2, "cm"),
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 10),
          plot.title = element_text(hjust = 0.5))
  
  
  plist[[i]] = p1
  
}

pdf('./fig/fig_8c_u251_87_transwell.pdf',width = 5,height = 4.8,onefile = F)
plist[[1]]+plist[[2]]
dev.off()

exper_fig[['trans_u251']] = plist[[1]]
exper_fig[['trans_u87']] = plist[[2]]


####### 6.4 RNA sh ########
# rm(list = ls())
# gc()
df = read_xlsx('./data/experi/rnaSh.xlsx')
df[df$celltype=='U251',c('control', "sh-MEST#1" ,"sh-MEST#2")] = na.aggregate(df[df$celltype=='U251',c('control', "sh-MEST#1" ,"sh-MEST#2")],2,mean)
df[df$celltype=='U87',c('control', "sh-MEST#1" ,"sh-MEST#2")] = na.aggregate(df[df$celltype=='U87',c('control', "sh-MEST#1" ,"sh-MEST#2")],2,mean)


plist = list()
for (i in c('U251','U87')) {
  
  cellname = i
  # 宽数据转换为长数据格式
  df%>% dplyr::filter(df$celltype==cellname,) %>% 
    pivot_longer(cols = 2:4,
                 names_to = "Group",
                 values_to ="Value")->Figure_2b
  # 分组标签
  Group<-unique(Figure_2b$Group)
  # 分组因子化
  Figure_2b$Group<-factor(Figure_2b$Group,levels = unique(Figure_2b$Group),labels = c('control', "sh-MEST#1" ,"sh-MEST#2"))
  my_comparisons <- list(c("control","sh-MEST#1"),c("control","sh-MEST#2"))
  
  
  p1=ggplot(Figure_2b, aes(Group, Value))+ 
    # 添加散点图
    geom_jitter(aes(x=Group, Value,colour = Group), size = 2.5, width = 0.1, height = 0)+
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +     # 误差棒
    stat_summary(fun = mean, geom = "crossbar", width = 0.5, size = 0.2)+  # 指代均数的水平横线
    theme_classic()+scale_color_manual(values =c('#375631','#D3C3D0','#3C224B') )+
    
    stat_compare_means(comparisons=my_comparisons,method = "t.test",hide.ns = F,label = "p.format")+
    labs(title = cellname)+
    ylab("Relative mRNA level") +
    xlab("") +
    #ylim(35,37.5)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
          legend.position = "none",
          axis.ticks = element_line(size=0.2, color="black"),
          axis.ticks.length = unit(0.2, "cm"),
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 10),
          plot.title = element_text(hjust = 0.5))
  
  
  plist[[i]] = p1
  
}

pdf('./fig/fig_8c_u251_87_MEST_sh_mRNA.pdf',width = 5,height = 4.8,onefile = F)
plist[[1]]+plist[[2]]
dev.off()
exper_fig[['sh_mRNA_u251']] = plist[[1]]
exper_fig[['sh_mRNA_u87']] = plist[[2]]


####### 6.5 protien sh #######
# rm(list = ls())
df = read_xlsx('./data/experi/protein_mest_if.xlsx')
df[df$celltype=='U251',c('control', "sh-MEST#1" ,"sh-MEST#2")] = na.aggregate(df[df$celltype=='U251',c('control', "sh-MEST#1" ,"sh-MEST#2")],2,mean)
df[df$celltype=='U87',c('control', "sh-MEST#1" ,"sh-MEST#2")] = na.aggregate(df[df$celltype=='U87',c('control', "sh-MEST#1" ,"sh-MEST#2")],2,mean)


plist = list()
for (i in c('U251','U87')) {
  
  cellname = i
  # 宽数据转换为长数据格式
  df%>% dplyr::filter(df$celltype==cellname,) %>% 
    pivot_longer(cols = 2:4,
                 names_to = "Group",
                 values_to ="Value")->Figure_2b
  # 分组标签
  Group<-unique(Figure_2b$Group)
  # 分组因子化
  Figure_2b$Group<-factor(Figure_2b$Group,levels = unique(Figure_2b$Group),labels = c('control', "sh-MEST#1" ,"sh-MEST#2"))
  my_comparisons <- list(c("control","sh-MEST#1"),c("control","sh-MEST#2"))
  
  
  p1=ggplot(Figure_2b, aes(Group, Value))+ 
    # 添加散点图
    geom_jitter(aes(x=Group, Value,colour = Group), size = 2.5, width = 0.1, height = 0)+
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +     # 误差棒
    stat_summary(fun = mean, geom = "crossbar", width = 0.5, size = 0.2)+  # 指代均数的水平横线
    theme_classic()+scale_color_manual(values =c('#375631','#D3C3D0','#3C224B') )+
    
    stat_compare_means(comparisons=my_comparisons,method = "t.test",hide.ns = F,label = "p.format")+
    labs(title = cellname)+
    ylab("Fluorescence intensity of MEST") +
    xlab("") +
    #ylim(35,37.5)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
          legend.position = "none",
          axis.ticks = element_line(size=0.2, color="black"),
          axis.ticks.length = unit(0.2, "cm"),
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 10),
          plot.title = element_text(hjust = 0.5))
  
  
  plist[[i]] = p1
  
}

pdf('./fig/fig_8c_u251_87_MEST_sh_proten_if.pdf',width = 5,height = 4.8,onefile = F)
plist[[1]]+plist[[2]]
dev.off()

exper_fig[['sh_proten_u251']] = plist[[1]]
exper_fig[['sh_proten_u87']] = plist[[2]]

####### 6.6 mRNA different grade #########
# rm(list = ls())
# gc()
library(readxl)
df=read_xlsx('./data/experi/MrnaGrade.xlsx')

df1 = df %>% pivot_longer(cols = 1:4, names_to = "Group",
                          values_to ="Value") %>% na.omit()

df1$Group[df1$Group != 'Normal'] = 'Tumor core'
table(df1$Group)

df1$Group[df1$Group =='Normal'] = 'Adjacent tissue'

p=ggstripchart(df1, "Group", "Value",
             color = "Group", palette= c('#375631','#3C224B'),alpha=0.5,
             add = "boxplot")+
  theme_classic()+
  # ggsignif::geom_signif(comparisons = list(c("Adjacent","Tumor")),step_increase = 0.08,
  #                       test = "t.test")+
  ylab("MEST mRNA expression") +
  xlab("") +
  labs(title = "Tissues")+
  #ylim(0, 10)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")
library(ggpval)
p1 =ggpval::add_pval(ggplot_obj = p,
                     pairs = list(c(1,2)),fold_change = T,
                     test = 't.test')

p1
ggsave('./fig/fig_8d_mRNA_normal_grade_mest.pdf',width = 2.5,height = 4.8,onefile = F)

### 另一个作图方式
if(F){
  
  rm(list = ls())
  library(readxl)
  df=read_xlsx('./data/experi/MrnaGrade.xlsx')
  
  df1 = df %>% pivot_longer(cols = 1:4, names_to = "Group",
                            values_to ="mRNA_Value") %>% na.omit()
  
  ##组内差异显著性
  library(rstatix)
  
  df1 %>% 
    # group_by(group) %>% 
    t_test(mRNA_Value ~ Group,
           ref.group = "Normal") %>% 
    mutate(p = round(p, 5),
           p.signif = case_when(p < 0.05 ~ "",
                                .default = " (NS)")) -> df_test
  
  df1 %>% 
    group_by(Group) %>% 
    slice_max(mRNA_Value, with_ties = FALSE) %>% 
    left_join(df_test,
              join_by(Group==group2)) %>% 
    mutate(Group = factor(Group, levels=c("Normal", "G2", "G3", "G4"))
    )-> df_lab
  
  df1$Group = factor(df1$Group,levels=c("Normal", "G2", "G3", "G4"))
  
  ##作图
  pal =c('#66529F','#A37E7D','#BF9895','#DAB2B2')
  ggplot(data = df1 ,
         aes(x=Group, y=mRNA_Value)) +
    geom_bar(aes(fill=Group),
             stat = "summary", fun="mean", color="black",
             position = position_dodge2(padding = 0.25), width = 0.8) +
    stat_summary(aes(group=Group),
                 geom = "errorbar",fun.data = "mean_se",
                 width=0.25, position = position_dodge(width = 0.8)) +
    geom_jitter(aes(color=Group),
                position = position_jitterdodge(jitter.width = 1,jitter.height = 0), 
                shape=21, fill="white", size=2) +
    geom_text(data = df_lab,
              aes(x=Group, y=mRNA_Value + 1, group=Group,label=str_c("P = ", p, p.signif)),
              position = position_dodge(width = 0.8),
              angle=90, size=3.5) +
    scale_fill_manual(values = pal) +
    scale_color_manual(values = c("black", "black", "black", "black"),
                       guide = "none") +
    scale_y_continuous(limits = c(0, 10),
                       expand = c(0, 0),
                       breaks = seq(0, 10, 2)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
          axis.ticks = element_line(size=0.2, color="black"),
          axis.ticks.length = unit(0.2, "cm"),
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 10),
          plot.title = element_text(hjust = 0.5),
          legend.position = "none")
  
  
  
}




exper_fig[['normal_grade']] = p1

####### 6.7 mRNA cell lines#########
# rm(list = ls())
# gc()

df=read_xlsx('./data/experi/MrnaCell.xlsx')

df1 = df %>% pivot_longer(cols = 1:4, names_to = "Group",
                          values_to ="Value") %>% na.omit()
df1$Group = factor(df1$Group,levels = c("HA1800","A172","U87","U251"))
p=ggstripchart(df1, "Group", "Value",
             color = "Group", palette= c('#375631','#D3C3D0','#3C224B','#982C2C'),alpha=0.5,
             add = "boxplot")+
  theme_classic()+
  ggsignif::geom_signif(comparisons = list(c("HA1800","A172"),
                                           c("HA1800","U87"),
                                           c("HA1800","U251")
                                          ),step_increase = 0.08,
                        test = "t.test")+
  ylab("MEST mRNA expression") +
  xlab("") +
  labs(title = "Cell lines")+
  #ylim(0, 10)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")
p

ggsave('./fig/fig_8d_mRNA_cellline_mest.pdf',width = 2.5,height = 4.8,onefile = F)

exper_fig[['cellline']] = p


####### 6.8 mRNA paired############
# rm(list = ls())
# gc()

df=read_xlsx('./data/experi/paired_tissue.xlsx')

# 数据格式转换
df%>%
  pivot_longer(cols = 2:3,
               names_to = "Group",
               values_to = "Value")->Figure_2f

# 因子化
Figure_2f$Group<-factor(Figure_2f$Group,levels = unique(Figure_2f$Group))
# 可视化
Fig_2f<-
  ggplot(Figure_2f,aes(x=Group,y=Value,color=Group,group=id))+#按照patient分组
  geom_line(color="grey")+#配对线图
  geom_point(size=4)+#散点图
  scale_color_manual(values = alpha(c("#375631","#3C224B"),0.7))+
  scale_y_continuous(limits = c(min(Figure_2f$Value),9))+
  labs(x="",y='MEST mRNA expression')+
  theme_classic()+
  labs(title = "Paired Tissues")+
  
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")+
  stat_compare_means(method = "t.test",paired = TRUE, comparisons=list(c("Normal", "Tumor")))#配对t检验

  
Fig_2f
ggsave('./fig/fig_8e_mRNA_paired_tissues_mest.pdf',width = 2.5,height = 4.8,onefile = F)


exper_fig[['paired_tissues']] = Fig_2f

####### 6.9 cck ############
# rm(list = ls())
# gc()

df= cck <- read_excel("data/experi/cck.xls")


plist= list()
ptab = data.frame()

for (i in c('U251','U87')) {
  cellname = i
  Figure_2a = df[df$celltype == cellname,]
  # 每一组的数据除以第一天的作为标准化
  # CCK
  resmean = apply(Figure_2a[Figure_2a$time==1,3:5], 2, mean)
  # Figure_2a$control = (Figure_2a$control)/resmean[1]
  # Figure_2a$`sh-MEST#1` = Figure_2a$`sh-MEST#1` /resmean[2]
  # Figure_2a$`sh-MEST#2` =  Figure_2a$`sh-MEST#2`/resmean[3]
    
  
  Figure_2a= Figure_2a %>% pivot_longer(cols = 3:5,
                 names_to = "Group",
                 values_to = "Value")
  
  Figure_2a$Value = (Figure_2a$Value)/resmean[1]
  Figure_2amean = Figure_2a %>% dplyr::group_by(time,Group) %>% 
    dplyr::summarise(mean = mean(Value))
  df.new = left_join(Figure_2a, Figure_2amean,by =c('Group','time'))
  


  #因子化
  df.new$Group<-factor(df.new$Group,levels = c('control', "sh-MEST#1" ,"sh-MEST#2"))
  # 可视化
  Fig_2a<-
    ggplot(df.new,aes(x=time,y=Value,color=Group,group=Group))+
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +     # 误差棒
    stat_summary(fun = mean, geom = "crossbar", width = 0.4, size = 0.2)+  # 指代均数的水平横线
    geom_point( aes(y=mean,color=Group,group=Group))+#点图
    geom_line(aes(y=mean,color=Group,group=Group))+#线图
    labs(x="Time (d)",y='Cell viability (OD 450nm)')+#设置坐标轴标题格式
    labs(title = cellname)+
    
    # scale_x_continuous(expand = c(0,0))+#设置坐标轴刻度及起始
    # scale_y_continuous(limits = c(0,4),breaks = seq(0,4,1),expand = c(0,0))+
    scale_color_manual(values = c('#375631','#D3C3D0','#3C224B'),name='',
                       labels=c('control (n=3)', "sh-MEST#1 (n=3)" ,"sh-MEST#2 (n=3)"))+
    theme_classic()+
    theme(legend.position ="top",
          legend.title = element_blank(),#不显示图例标题
          legend.text =element_text(size =8),#不显示图例文本
          legend.key.width = unit(1,"cm"),#图例宽度
          strip.background = element_blank(),
          strip.text = element_text(size = 12,face = "bold"),
          panel.spacing =unit(0,'cm'),
          legend.key.height = unit(1,"cm"),#图例高度
          plot.margin = margin(r=0.5,t=0.5,unit = "cm"),#设置画板边缘大小
          axis.text = element_text(color = "black"),
          plot.title = element_text(hjust = 0.5)
          
    )#不显示坐标轴文本
  
  Fig_2a
  plist[[i]]=Fig_2a
  
  print(i)
  for (j in c(1:5)) {
    cat(paste0('day ',j))
    dt=df.new[df.new$time==j,c('Value','Group')] %>% as.data.frame()
    dt$Value = as.numeric(dt$Value)
    tab1 <- dt %>% t_test(Value ~ Group,  alternative = "two.sided")
   
    tab1 = tab1[tab1$group1 == 'control',]
    print(tab1$p)
    tmp = data.frame(celltype = i, day = j, group = tab1$group2, pvalue = tab1$p)
    ptab = rbind(ptab,tmp)
  }
  
  
}

pdf('./fig/fig_8f_cck_mest.pdf',width = 5,height = 4.8,onefile = F)

plist[[1]]+plist[[2]]+plot_layout(nrow =1,heights  = c(1,1),guides='collect' )&
    theme(legend.position='top')
dev.off()

write_csv(ptab,file = './res/cck8_test_pvalue_tab.csv')
# ptab

# 
# Figure_2a = df %>% 
#   pivot_longer(cols = 3:5,
#                names_to = "Group",
#                values_to = "Value")
# a = Figure_2a[Figure_2a$celltype =='U251'&Figure_2a$time==5&Figure_2a$Group == 'control',]$Value
# b = Figure_2a[Figure_2a$celltype =='U251'&Figure_2a$time==5&Figure_2a$Group == 'sh-MEST#1',]$Value
# c = Figure_2a[Figure_2a$celltype =='U251'&Figure_2a$time==5&Figure_2a$Group == 'sh-MEST#2',]$Value
# 
# 
# t.test(a,b)
# t.test(a,c)
# 
# a = Figure_2a[Figure_2a$celltype =='U87'&Figure_2a$time==5&Figure_2a$Group == 'control',]$Value
# b = Figure_2a[Figure_2a$celltype =='U87'&Figure_2a$time==5&Figure_2a$Group == 'sh-MEST#1',]$Value
# c = Figure_2a[Figure_2a$celltype =='U87'&Figure_2a$time==5&Figure_2a$Group == 'sh-MEST#2',]$Value
# 
# t.test(a,b)
# t.test(a,c)

# 
# 
# pdf('./fig/fig_8_experi_figs.pdf',width = 12,height = 14.4,onefile = F)
# exper_fig[[1]]+
# exper_fig[[2]]+
# exper_fig[[3]]+
# exper_fig[[4]]+
# exper_fig[[5]]+
# exper_fig[[6]]+
# exper_fig[[7]]+
# exper_fig[[8]]+
# exper_fig[[9]]+
# exper_fig[[10]]+
# exper_fig[[11]]+
# # exper_fig[[12]]+ 
#   plot_layout(nrow =3,heights = c(1,1,1),widths = c(1,1,1,1),guides='collect' )&
#   theme(legend.position='right')
# dev.off()
# 
# 
 

####### 7 # MEST sh RNA-Seq preprocess ######################
rm(list = ls())
gc()


library(data.table)
library(limma)
library(ggplot2)
library(ggpubr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ComplexHeatmap)
library(GseaVis)

expr<-fread("/export3/zhangw/Project_Cross/Project_MEST/data/RNA-seq/MEST/gene_tpm_matrix.csv",data.table = F)
# tmp<-dplyr::select(expr,"gene_id")
# tmp$ENSEMBL<-stringr::str_split(tmp$gene_id,"\\|",simplify = T)[,1]
# tmp$SYMBOL<-stringr::str_split(tmp$gene_id,"\\|",simplify = T)[,2]
# expr$gene_id<-tmp$SYMBOL
# expr<-aggregate(x = expr[,2:(ncol(expr))],
#                 by = list(expr$gene_id),
#                 FUN = mean)
# row.names(expr)<-expr[,1]
# expr<-expr[,-1]
# 
# ##log2(tpm+1)
# exprSet<-log2(expr+1)
# 
# phe<-as.data.frame(colnames(exprSet))
# phe$treatment<-c(replicate(3,"Overexpression"),replicate(3,"Control"))
# phe$project<-"TSPAN7 Overexpression"
# colnames(phe)[1]<-"sample"
# 
# saveRDS(phe,file = "./processed/u87_ov_tspan7_phe.rds")
# saveRDS(exprSet,file = "./processed/u87_ov_tspan7_expr.rds")
# 
# # PCA
# pca <- prcomp(t(exprSet), scale=F)
# pca.data <- data.frame(Sample=rownames(pca$x),
#                        X=pca$x[,1],
#                        Y=pca$x[,2])
# identical(phe$sample,rownames(pca.data))
# pca.data<-cbind(pca.data,phe)
# pca.var <- pca$sdev^2  ## sdev是标准偏差，十个样本，就有十个标准偏差，平方是避免负数的干扰
# pca.var.per <- round(pca.var/sum(pca.var)*100, 1)  ##求每个样本的variation
# 
# pdf("./result/u87_ov/pca.pdf",bg="white",width=6.5,height=4.5,pointsize=6)
# ggplot(pca.data, aes(x=X,y=Y,color=treatment),size=5) +
#   geom_point(size=3) +#stat_ellipse(level = 0.95, show.legend = F,geom = "polygon",alpha = 1/5, aes(fill = patient))+
#   #scale_colour_manual(values =c("#F5A200","#A70B00"))+
#   xlab(paste("PC1(",pca.var.per[1],"%","variance)",sep="")) +
#   ylab(paste("PC2(",pca.var.per[2],"%","variance)",sep="")) +
#   theme_classic()+
#   ggsci::scale_color_lancet()+
#   theme(axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),
#         legend.key.size = unit(0.8, "cm"),legend.text = element_text(size=15),legend.title =element_text(size=15))
# dev.off()
# 
# # tspan7 expression
# phe_select<-phe
# phe_select$gene<-t(exprSet["TSPAN7",])
# 
# pdf("./result/u87_ov/TSPAN7_expression.pdf",bg="white",width=4,height=5.2,pointsize=6)
# ggstripchart(phe_select, "treatment", "gene",
#              color = "treatment", palette="lancet",alpha=0.7,size=2,
#              add = "mean_sd")+
#   theme_classic()+
#   stat_compare_means(aes(group = treatment),method = "t.test",hide.ns = F,label = "p.format",label.x = 1.5) +
#   ylab("TSPAN7 expression") +
#   xlab("") +
#   labs(title = "TSPAN7 overexpression")+
#   #ylim(0, 10)+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
#         #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
#         axis.ticks = element_line(size=0.2, color="black"),
#         axis.ticks.length = unit(0.2, "cm"),
#         axis.title = element_text(size = 10),
#         axis.text = element_text(size = 10),
#         #legend.position = "none",
#         plot.title = element_text(hjust = 0.5))
# dev.off()
# 
# #  calcium channels expression
# cav<-c("CACNA1S","CACNA1C","CACNA1D","CACNA1F","CACNA1A","CACNA1B","CACNA1E","CACNA1G","CACNA1H","CACNA1I")
# 
# phe_select<-cbind(phe_select,t(exprSet[cav,]))
# phe_select<-melt(phe_select)
# 
# ggstripchart(phe_select, "treatment", "value",facet.by="variable",
#              color = "treatment", palette="nejm",alpha=0.7,size=2,
#              add = "mean_sd")+
#   theme_classic()+
#   stat_compare_means(aes(group = treatment),method = "t.test",hide.ns = F,label = "p.format") +
#   ylab("TSPAN7 expression") +
#   xlab("") +
#   labs(title = phe_select$Cell_line)+
#   #ylim(0, 10)+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
#         #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
#         axis.ticks = element_line(size=0.2, color="black"),
#         axis.ticks.length = unit(0.2, "cm"),
#         axis.title = element_text(size = 10),
#         axis.text = element_text(size = 10),
#         #legend.position = "none",
#         plot.title = element_text(hjust = 0.5))
# 
# ###---7.1 DEG------------
# #control vs ov
# group <- phe$treatment
# group <- factor(group,levels = c("Control","Overexpression"))
# design <- model.matrix(~group)
# fit <- lmFit(exprSet,design)
# fit2 <- eBayes(fit)
# allDiff=topTable(fit2,adjust='fdr',coef=2,number=Inf)
# allDiff["TSPAN7",]
# 
# saveRDS(allDiff,file = "./processed/u87_ov_tspan7_DEG.rds")
# 
# allDiff$significant<-ifelse(allDiff$logFC>0.5&allDiff$P.Value<0.05,"Up in TSPAN7 OV",
#                             ifelse(allDiff$logFC< -0.5&allDiff$P.Value<0.05,"Down in TSPAN7 OV","unsignificant"))
# table(allDiff$significant)
# 
# ###---7.2 GO KEGG-----------------
# DEG<-allDiff[allDiff$significant %in%c("Up in TSPAN7 OV","Down in TSPAN7 OV"),]
# ids <- bitr(row.names(DEG), fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db",drop = TRUE)
# DEG$gene<-rownames(DEG)
# DEG_list=merge(DEG,ids,by.x='gene',by.y='SYMBOL')
# DEG_list<-DEG_list[order(DEG_list$logFC,decreasing = T),]
# table(DEG_list$significant)
# 
# enrichGO_BP <- compareCluster(ENTREZID~significant, data=DEG_list, fun ="enrichGO",ont="BP",OrgDb='org.Hs.eg.db') 
# enrichGO_BP <-setReadable(enrichGO_BP, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
# tmp<-as.data.frame(enrichGO_BP)
# 
# ##gseGO
# geneList<-DEG_list$logFC
# names(geneList)<-DEG_list$ENTREZID
# gseGO <- gseGO(geneList = geneList, 
#                OrgDb = org.Hs.eg.db, 
#                keyType = "ENTREZID", 
#                ont="ALL",
#                pvalueCutoff =0.05)
# gseGO <-setReadable(gseGO, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
# 
# tmp<-as.data.frame(gseGO)

###-----7.3 u87 mest sh vs primary-----------------
expr<-fread("/export2/liuhw/MEST/gene_tpm_matrix.csv",data.table = F)
tmp<-dplyr::select(expr,"gene_id")
tmp$ENSEMBL<-stringr::str_split(tmp$gene_id,"\\|",simplify = T)[,1]
tmp$SYMBOL<-stringr::str_split(tmp$gene_id,"\\|",simplify = T)[,2]
expr$gene_id<-tmp$SYMBOL
expr<-aggregate(x = expr[,2:(ncol(expr))],
                by = list(expr$gene_id),
                FUN = mean)
row.names(expr)<-expr[,1]
expr<-expr[,-1]

##log2(tpm+1)
exprSet<-log2(expr+1)

phe<-as.data.frame(colnames(exprSet))
colnames(phe)[1]<-"sample"
phe$treatment<- substr(phe$sample,1,1)
phe$project<- ifelse(phe$treatment == 'N','no treatment', 'raw')
# colnames(phe)[1]<-"sample"

saveRDS(exprSet, file = './data/u87_sh_mest_expr.rds')
saveRDS(phe, file = './data/u87_sh_mest_phe.rds')


u87_sh_mest_expr <- readRDS('./data/u87_sh_mest_expr.rds')
u87_sh_mest_phe <- readRDS('./data/u87_sh_mest_phe.rds')


# ref_gene = u87_sh_mest_expr['GAPDH',]
# 
# expr_s = apply(u87_sh_mest_expr, 2,  function(x){
#   
#   y = x/c(x[which(rownames(u87_sh_mest_expr)=='GAPDH')])
#   
# }) %>% as.data.frame()
# 
# 
# 
# expr_s['MEST',]


# identical(rownames(u87_ov_tspan7_expr),rownames(exprSet))

# u87_tspan7_proj_expr<-cbind(exprSet,u87_ov_tspan7_expr)
# u87_tspan7_proj_phe<-rbind(phe,u87_ov_tspan7_phe)

##remove batch
# expr_select <- removeBatchEffect(u87_tspan7_proj_expr, batch = u87_tspan7_proj_phe$project)


# PCA
pca <- prcomp(t(u87_sh_mest_expr), scale=F)
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2])
identical(u87_sh_mest_phe$sample,rownames(pca.data))
pca.data<-cbind(pca.data,u87_sh_mest_phe)
pca.var <- pca$sdev^2  ## sdev是标准偏差，十个样本，就有十个标准偏差，平方是避免负数的干扰
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)  ##求每个样本的variation
# pca.data<-pca.data[pca.data$treatment %in% c("Raw","Overexpression"),]

pdf("./fig/fig_9b_pca.pdf",bg="white",width=6.5,height=4.5,pointsize=6)
ggplot(pca.data, aes(x=X,y=Y,color=treatment),size=5) +
  geom_point(size=3) +#stat_ellipse(level = 0.95, show.legend = F,geom = "polygon",alpha = 1/5, aes(fill = patient))+
  #scale_colour_manual(values =c("#F5A200","#A70B00"))+
  xlab(paste("PC1(",pca.var.per[1],"%","variance)",sep="")) +
  ylab(paste("PC2(",pca.var.per[2],"%","variance)",sep="")) +
  theme_classic()+
  ggsci::scale_color_lancet()+
  theme(axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),
        legend.key.size = unit(0.8, "cm"),legend.text = element_text(size=15),legend.title =element_text(size=15))
dev.off()

# mest expression
phe_select<-u87_sh_mest_phe[u87_sh_mest_phe$treatment %in% c("M","N"),]
exprSet<-u87_sh_mest_expr[,phe_select$sample]

phe_select$MEST<-t(exprSet["MEST",]) %>% as.numeric()


###-----7.4 DEG------------
#raw vs sh
group <- phe_select$treatment
group = ifelse(group =='N', 'control','shMEST')
group <- factor(group,levels = c("control","shMEST"))
design <- model.matrix(~group)
fit <- lmFit(exprSet,design)
fit2 <- eBayes(fit)
allDiff=topTable(fit2,adjust='fdr',coef=2,number=Inf)
allDiff["MEST",]

saveRDS(allDiff,file = "./res/u87_sh_MEST_DEG.rds")

allDiff = readRDS('./res/u87_sh_MEST_DEG.rds')
allDiff$significant<-ifelse(allDiff$logFC>1&allDiff$P.Value<0.05,"Up in MEST sh",
                            ifelse(allDiff$logFC< -1&allDiff$P.Value<0.05,"Down in MEST sh","unsignificant"))
allDiff['MEST',]
table(allDiff$significant)

base_nt = apply(exprSet[,4:6], 1, mean) %>% as.data.frame()
colnames(base_nt) = 'Control'
base_sh = apply(exprSet[,1:3], 1, mean) %>% as.data.frame()
colnames(base_sh) = 'shMEST'
base_tab = cbind(base_nt,base_sh)
base_tab  = base_tab[rownames(allDiff),]

identical(rownames(base_tab),rownames(allDiff))

dfplot = cbind(allDiff,base_tab)

library(ImageGP)
library(ggplot2)
library(ggpubr)
library(egg)
library(ggrepel)

# dfplot[,c('Control','shMEST')] = apply(dfplot[,c('Control','shMEST')], 2, function(x){
#   y = log2(x+1)
#   return(y)
# })

dfplot$Sig<-ifelse(allDiff$logFC>1&allDiff$P.Value<0.05,"logFC>1 and P<0.05",
                            ifelse(allDiff$logFC< -1&allDiff$P.Value<0.05,"logFC<-1 and P<0.05",
                                   ifelse(allDiff$logFC< 0&allDiff$P.Value<0.05,'-1<logFC<0 and P<0.05',
                                          ifelse(allDiff$logFC> 0&allDiff$P.Value<0.05,'0<logFC<1 and P<0.05','P>0.05'))))

#转换为因子，指定绘图顺序；
dfplot$Sig <- factor(dfplot$Sig, levels = c("logFC>1 and P<0.05",
                                            '0<logFC<1 and P<0.05',
                                            'P>0.05',
                                            '-1<logFC<0 and P<0.05',
                                            "logFC<-1 and P<0.05"
                                            ))

#自定义散点配色：
mycol <- c('#eb2024','#f8979a','#808080','#9ebbda','#4159a8')

#自定义主题：
mytheme <- theme_bw() +
  theme(axis.title = element_text(size = 15), #坐标轴标题字号
        axis.text = element_text(size = 14), #坐标轴标签字号
        legend.title = element_text(size = 15), #图例标题字号
        legend.text = element_text(size = 14), #图例文本字号
        panel.grid.major = element_blank(), #去除背景网格
        panel.grid.minor = element_blank())

#ggplot2绘图火山图：
p= ggplot(data = dfplot,
            aes(x = Control, y = shMEST,
                color = Sig)) + #建立映射
  geom_point(size = 1.2)+ #绘制散点
  scale_colour_manual(name = "Rank differential\n(sh-MEST vs. Ctrl)", #图例标题修改
                      # labels =c('>+0.5', #图例文本修改
                      #           ' 0.25 to 0.49',
                      #           '-0.25 to +0.25',
                      #           '-0.49 to -0.25',
                      #           '>-0.5'),
                      values = alpha(mycol, 0.8)) + #更改散点配色和不透明度
  labs(x = 'Control:Log2 mRNA Abundance', #XY轴标题更改
       y = 'sh-MEST:Log2 mRNA Abundance') +
  guides(color = guide_legend(override.aes = list(size = 4.5))) + #放大图例散点
  scale_x_continuous(limits = c(0, 15), #y轴范围限制
                     breaks = seq(0, 15, by = 5), #y轴显示范围和间隔
                     expand = c(0,0)) + #x轴显示范围和间隔
  scale_y_continuous(limits = c(0, 15), #y轴范围限制
                     breaks = seq(0, 15, by = 5), #y轴显示范围和间隔
                     expand = c(0,0)) + #x轴显示范围和间隔
  mytheme+ #添加自定义主题
  annotate('text',x = 3, y = 12, label = 'N = 118', size = 6, color = '#eb2024') +
  annotate('text',x = 3, y = 10, label = 'N = 1594', size = 6, color = '#f8979a') +
  annotate('text',x = 13, y = 10, label = 'N = 56028', size = 6, color = '#808080') +
  annotate('text',x = 12, y = 4, label = 'N = 2590', size = 6, color = '#9ebbda') +
  annotate('text',x = 12, y = 2, label = 'N = 320', size = 6, color = '#4159a8') +
  annotate('text',x = 7, y = 14.2, label = 'N = 60650 genes', size = 7, color = 'black')


p


ras = read.gmt('./data/GOBP_RAS_PROTEIN_SIGNAL_TRANSDUCTION.v2023.2.Hs.gmt')
wnt =  read.gmt('./data/GOBP_POSITIVE_REGULATION_OF_WNT_SIGNALING_PATHWAY.v2023.2.Hs.gmt')

ras_dt =allDiff[ras$gene,] 
ras_dt =ras_dt%>% dplyr::filter(ras_dt$P.Value<0.05)
ras_dt$significant<-ifelse(ras_dt$logFC>0&ras_dt$P.Value<0.05,"Up in MEST sh",
                            ifelse(ras_dt$logFC< -0&ras_dt$P.Value<0.05,"Down in MEST sh","unsignificant"))
ras_dt

wnt_dt =allDiff[wnt$gene,] 
wnt_dt =wnt_dt%>% dplyr::filter(wnt_dt$P.Value<0.05).
wnt_dt$significant<-ifelse(wnt_dt$logFC>0&wnt_dt$P.Value<0.05,"Up in MEST sh",
                           ifelse(wnt_dt$logFC< -0&wnt_dt$P.Value<0.05,"Down in MEST sh","unsignificant"))

wnt_dt
# write.csv(ras_dt[order(ras_dt$significant),],file = './res/u87_sh_rna_ras_protien_DEGs.csv')
# write.csv(wnt_dt[order(wnt_dt$significant),],file = './res/u87_sh_rna_wnt_protien_DEGs.csv')

# RAS 
your_gene_ras = c('KRAS','MRAS','RASA4','RASA4B'
                  # , # qp validation
                  # 'NGF','CDC42EP3','ARHGEF2'
                  )

# WNT 
your_gene_wnt = c('WNT10B','ADGRA2','SFRP2','FRAT1'
                  # ,# qp validation
                  # 'TMEM132A','FGF2','BAMBI'
                  )


some_special_gene = c('MEST',your_gene_ras,your_gene_wnt
)


dfplot$label =""
dfplot$label[match(some_special_gene,rownames(dfplot))] <- some_special_gene

p+ geom_text_repel(data=dfplot, aes(label= label), color="black", size=3, fontface="italic",max.overlaps = 50000,
                   # nudge_x = 5,
                   box.padding = 0.5,
                   # nudge_y = 3,
                   segment.curvature = 0.2,
                   segment.ncp = 3,
                   segment.angle = 50)


ggsave('./fig/fig_9a_vocanol_mRNA_U87_sh_RNAseq.pdf',width = 8,height = 5,onefile = F)



###-----7.5-GO KEGG-----------------
DEG<-allDiff[allDiff$significant %in%c("Up in MEST sh","Down in MEST sh"),]
# DEG<-allDiff
ids <- bitr(row.names(DEG), fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db",drop = TRUE)
DEG$gene<-rownames(DEG)
DEG_list=merge(DEG,ids,by.x='gene',by.y='SYMBOL')
DEG_list<-DEG_list[order(DEG_list$logFC,decreasing = T),]
table(DEG_list$significant)

enrichGO <- compareCluster(ENTREZID~significant, data=DEG_list, fun ="enrichGO",ont="All",OrgDb='org.Hs.eg.db') 
enrichGO<-setReadable(enrichGO,OrgDb='org.Hs.eg.db')

saveRDS(enrichGO,file = "./res/mest_u87_sh_enrichGO.rds")

DEG_list[DEG_list$gene =='MEST',]

tmp<-as.data.frame(enrichGO)

write.csv(tmp,file = "./res/mest_u87_sh_enrichGO.csv")

df <- data.frame(tmp) %>%
  group_by(Cluster) %>%
  slice_head(n = 6) %>%
  arrange(desc(pvalue))


select_enrichment_go<-tmp[c(grep("cell cycle checkpoint signaling",tmp$Description),
                            grep("DNA replication",tmp$Description),
                            grep("DNA damage checkpoint signaling",tmp$Description),
                            grep("regulation of DNA recombination",tmp$Description),
                            # grep("regulation of Ras protein signal transduction",tmp$Description),
                            grep("cell growth",tmp$Description),
                            grep("GO:0048638",tmp$ID),
                            grep("GO:0046621",tmp$ID),
                            grep("GO:0048285",tmp$ID),
                            grep("cellular response to radiation",tmp$Description),
                            grep("Ras protein signal transduction",tmp$Description)),]
select_enrichment_go = rbind(select_enrichment_go,df)

# select_enrichment_go<-select_enrichment_go[-c(grep("negative ",select_enrichment_go$Description),
#                                               grep("positive",select_enrichment_go$Description),
#                                               grep("p53",select_enrichment_go$Description)),]
select_enrichment_go$type<-"GO"
select_enrichment<-select_enrichment_go[,-3]
# select_enrichment<-rbind(select_enrichment_kegg,select_enrichment_go[,-3])
select_enrichment$Description<-factor(select_enrichment$Description,levels = c(unique(select_enrichment$Description)))
select_enrichment$type<-factor(select_enrichment$type,levels = c("GO"))

ratio <- lapply(select_enrichment$GeneRatio,function(x){as.numeric(eval(parse(text = x)))}) %>% unlist()
select_enrichment$ratio <- ratio

select_enrichment$Description <- factor(select_enrichment$Description,levels = select_enrichment$Description)

library(ggforce)
ggplot(select_enrichment) +
  ggforce::geom_link(aes(x = 0,y = Description,
                         xend = -log10(pvalue),yend = Description,
                         alpha = stat(index),
                         color = Cluster,
                         size = after_stat(index)),
                     n = 500,
                     # color = "#FF0033",
                     show.legend = F) +
  geom_point(aes(x = -log10(pvalue),y = Description),
             color = "black",
             fill = "white",size = 6,shape = 21) +
  geom_line(aes(x = ratio*100,y = Description,group = 1),
            orientation = "y",linewidth = 1,color = "#FFCC00") +
  scale_x_continuous(sec.axis = sec_axis(~./100,
                                         labels = scales::label_percent(),
                                         name = "Percent of geneRatio")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.text = element_text(face = "bold.italic"),
        axis.text = element_text(color = "black")) +
  ylab("") + xlab("-log10 Pvalue") +
  facet_wrap(~Cluster,scales = "free",ncol = 1) +
  scale_color_manual(values = c('#375631','#D3C3D0','#3C224B')[-2] %>% rev())
  # scale_color_brewer(palette = "Set1")

ggsave("./fig/fig_8a_enrichGO_bp.pdf",bg="white",width=8.0,height=8.5,pointsize=12)



# tmp<-tmp[order(tmp$NES,decreasing = T),]
# ## select
# GO_select<-tmp[c(grep("endoplasmic reticulum stress",tmp$Description),grep("unfolded protein",tmp$Description),
#                  grep("stress",tmp$Description),grep("tumor necrosis factor",tmp$Description),
#                  grep("cell cycle",tmp$Description),grep("chromosome",tmp$Description),
#                  grep("proliferation",tmp$Description),grep("DNA",tmp$Description),
#                  grep("cell migration",tmp$Description),grep("junction",tmp$Description),
#                  grep("adhesion",tmp$Description)),]
# GO_select<-GO_select[-c(grep("negative",GO_select$Description),
#                         grep("positive",GO_select$Description),
#                         grep("regulation",GO_select$Description),
#                         grep("mononuclear",GO_select$Description),
#                         grep("fibroblast",GO_select$Description),
#                         grep("epithelial",GO_select$Description),
#                         grep("endothelial",GO_select$Description)),]
# GO_select<-GO_select[!duplicated(GO_select$Description),]
# GO_select$Description<-factor(GO_select$Description,levels=unique(GO_select$Description))
# go_list<-list("Group1"=data.frame(GO_select[1:2,"Description"], "col"="#ED0000"),
#               "Group2"=data.frame(GO_select[3:6,"Description"], "col"="#00468B"),
#               "Group3"=data.frame(GO_select[7:10,"Description"], "col"="#42B540"),
#               "Group4"=data.frame(GO_select[11:28,"Description"], "col"="#0099B4"),
#               "Group5"=data.frame(GO_select[29:43,"Description"], "col"="#925E9F"),
#               "Group6"=data.frame(GO_select[44:45,"Description"], "col"="#FDAF91"),
#               "Group7"=data.frame(GO_select[46:50,"Description"], "col"="#AD002A"),
#               "Group8"=data.frame(GO_select[51:62,"Description"], "col"="#ADB6B6"))
# GO_select$Description<-forcats::fct_rev(GO_select$Description)

# pdf(file="./result/u87_ov/gsea_gobp.pdf",bg="white",width=8.0,height=8.5,pointsize=12)
# ggplot(data = GO_select[-c(7:10),], aes(x =NES , y = Description, fill = qvalues)) +
#   scale_fill_gradientn(colours=paletteer::paletteer_c("grDevices::Geyser", 30,direction = -1))+
#   geom_bar(stat = "identity",position = "identity",color="NA",size=0.25) +
#   theme_classic()+
#   labs(title = "GSEA GO")+
#   ylab("") +
#   xlab("NES") +
#   theme(plot.title = element_text(hjust = 0.5),
#         axis.text.y.left = (element_text(color=c(go_list[[8]]$col,go_list[[7]]$col,go_list[[6]]$col,go_list[[5]]$col,go_list[[4]]$col,go_list[[2]]$col,go_list[[1]]$col))),)
# dev.off()

# gsea
DEG<-allDiff %>% dplyr::arrange(desc(logFC))
ids <- bitr(row.names(DEG), fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db",drop = TRUE)
DEG$gene<-rownames(DEG)
DEG_list=merge(DEG,ids,by.x='gene',by.y='SYMBOL')
DEG_list<-DEG_list[order(DEG_list$logFC,decreasing = T),]

# 提取foldchange
gene_fc <- DEG_list$logFC 

# 将gene_fc对象的命名改为df_all_sort数据框中ENTREZID列的值
names(gene_fc) <- DEG_list$ENTREZID


ego3 <- gseGO(geneList     = gene_fc,
              OrgDb        = 'org.Hs.eg.db',
              ont          = "BP",
              minGSSize    = 100,
              maxGSSize    = 1000,
              pvalueCutoff = 0.1,
              verbose      = T)
ego3_table <- as.data.frame(ego3)


# write_csv(ego3_table,file = './res/u87_sh_mest_gseGO.csv')
# saveRDS(ego3,'./res/u87_sh_mest_gseGO.rds')

ego3_table = fread('./res/u87_sh_mest_gseGO.csv')

library(GseaVis,lib.loc =  "/export/bioinfo-team/home/liuhw/R/x86_64-pc-linux-gnu-library/4.1")

geneSetID = c('GO:0007265',
              'GO:0030177',
              'GO:0034976','GO:0010508'
)

geneSetID%in% ego3_table$ID

# all plot
gseaNb(object = ego3,
       geneSetID = geneSetID,
       subPlot = 2,
       pCol = 'black',
       pHjust = 0,
       rmSegment=T,
       termWidth = 35,
       legend.position = c(0.75,0.71),addPval = T,
       pvalX = 0.05,pvalY = 0.05, curveCol = c('#3C224B','#A786BA','#0d5b26','#96D2B0'),
       rankCol = c("#26456E", "#FFF4EE", "#7B3014")
)
ggsave('./fig/fig_8a_gsea_RNAseq_mest.pdf',width = 7,height = 5)


###-----7.6  expression U87 MEST and RAS and wnt pathway genes ######################

# df = fread('/export2/liuhw/MEST/gene_count_matrix.csv')
df = fread('/export2/liuhw/MEST/gene_tpm_matrix.csv')


df[grep('ENSG00000106484',df$gene_id),]


df1 = df %>% dplyr::filter(df$gene_id =='ENSG00000106484|MEST')%>% pivot_longer(cols = 2:7, names_to = "Group",
                          values_to ="Value") %>% na.omit()
str(df1)
df1$Group = substr(df1$Group,1,1)
df1$Group = ifelse(df1$Group =='N', 'control',"sh-MEST")

df1$Group = factor(df1$Group,levels = c("control","sh-MEST"))
p=ggstripchart(df1, "Group", "Value",
               color = "Group", palette= c('#375631','#D3C3D0','#3C224B','#982C2C'),alpha=0.5,
               add = "boxplot")+
  theme_classic()+
  ggsignif::geom_signif(comparisons = list(c("control","sh-MEST")
                                          
                                          
  ),step_increase = 0.08,
  test = "t.test")+
  ylab("MEST mRNA expression(TPM)") +
  xlab("") +
  labs(title = "U87")+
  #ylim(0, 10)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")
p

ggsave('./fig/fig_9a_mRNA_U87_sh_RNAseq_mest.pdf',width = 2.5,height = 4.8,onefile = F)


########

u87_sh_mest_expr <- readRDS("/export3/zhangw/Project_Cross/Project_MEST/data/u87_sh_mest_expr.rds")

# RAS 
your_gene_ras = c('KRAS','MRAS','RASA4','RASA4B')
# WNT 
your_gene_wnt = c('WNT10B','ADGRA2','SFRP2','FRAT1')

df =u87_sh_mest_expr[c(your_gene_ras,your_gene_wnt),] %>%
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column('ID') %>% 
  dplyr::mutate(type = substr(ID,1,1)) %>%
  dplyr::select(-ID) %>% 
  pivot_longer(cols = -type,
               names_to = "Gene",
               values_to ="Value")

Figure_2b = df %>% dplyr::select('Gene','type','Value')
Figure_2b$type = ifelse(Figure_2b$type=='N','control','sh-MEST')
colnames(Figure_2b) = c('gene','Group',"Value")

# 分组因子化
Figure_2b$Group<-factor(Figure_2b$Group,levels =  c('control', "sh-MEST" ),labels = c('control', "sh-MEST" ))
Figure_2b$gene = factor(Figure_2b$gene,levels = unique(Figure_2b$gene))
my_comparisons <- list(c("control","sh-MEST"))

color = c('#375631','#D3C3D0','#3C224B')[-2] %>% alpha(0.5)

source('/export3/zhangw/Code.sum/PlotGroup_barplot.R')
source('/export3/zhangw/Code.sum/some_meaningful_theme_for_ggplot.R')
data = Figure_2b


names(data) = c('Sample','Group',"Value")

PlotGroup_barplot(data = data,xlab = 'Sample',
                  ylab = 'Value',
                  legend = 'Group',
                  color = color,
                  show_compare = F,cmp_show_compare = T,p_show = 'p')+

  labs(title = 'U87 RNA-Seq')+
  ylab("Relative mRNA level") +
  xlab("") +
  my_theme_zhangwei+
  theme(legend.position = "top")

ggsave('./fig/fig_8d_87_RNAseq_MEST_sh_ras_wnt.pdf',width = 8,height = 4,onefile = F)


########7.7 heatmap of rna sequencing for some pathway genes ########

u87_sh_mest_expr <- readRDS("/export3/zhangw/Project_Cross/Project_MEST/data/u87_sh_mest_expr.rds")

u87_sh_mest_phe <- readRDS("/export3/zhangw/Project_Cross/Project_MEST/data/u87_sh_mest_phe.rds")


# RAS 
your_gene_ras = c('KRAS','MRAS','RASA4','RASA4B'
                  # , # qp validation
                  # 'NGF','CDC42EP3','ARHGEF2'
)

# WNT 
your_gene_wnt = c('WNT10B','ADGRA2','SFRP2','FRAT1'
                  # ,# qp validation
                  # 'TMEM132A','FGF2','BAMBI'
)

phe_select<-u87_sh_mest_phe[u87_sh_mest_phe$project %in% c("raw","no treatment"),]
exprSet<-u87_sh_mest_expr[,phe_select$sample]
gene_select<-data.frame("gene"=c('KRAS','MRAS','RASA4','RASA4B', 'NGF','CDC42EP3','ARHGEF2', #RAS
                                 'WNT10B','ADGRA2','SFRP2','FRAT1','TMEM132A','FGF2','BAMBI', # WNT
                                 "ATF6", "ERN1",  "EIF2AK3",  "XBP1", "HSPA5",  "DDIT3",  "PPP1R15A",  "ATF4" # ER
                                 ),
                        "Category"=c(rep("RAS",7),rep("WNT",7),
                                     rep("ER stress",8)))

phe_select2<-cbind(phe_select,t(exprSet[gene_select$gene,]))
rownames(phe_select2)<-phe_select2$sample
phe_select2$treatment =ifelse(phe_select2$treatment=='N','Control','sh-MEST')

annotation_colors = list(treatment=c("Control"="#375631","sh-MEST"="#D3C3D0"),
                         Category=c("RAS"="#9DB4CE","WNT"="#F9C08A",
                                    "ER stress"="#2F7D77"))

pdf(file="./fig/fig_8f_DEG_heatmap.pdf",bg="white",width=5.5,height=5.2,pointsize=12)
ComplexHeatmap::pheatmap(as.matrix(t(phe_select2[,-c(1:3)])),scale = "row",cluster_rows = F,cluster_cols = T,
                         show_colnames = F,border=T,border_color ="white",
                         name="Expression",
                         # gaps_row=c(8,13,20),
                         col=colorRampPalette(c(paletteer::paletteer_c("grDevices::Geyser", 30,direction = 1)))(30),
                         annotation_col = dplyr::select(phe_select2,c("treatment")),
                         annotation_row = dplyr::select(gene_select,"Category"),
                         annotation_colors = annotation_colors,
                         treeheight_col = 15
)
dev.off()

phe_select2<-cbind(phe_select,t(exprSet[gene_select$gene,]))
phe_select2<-melt(phe_select2)

dfplot = phe_select2[,c('variable','treatment','value')]

source('/export3/zhangw/Code.sum/PlotGroup_barplot.R')
source('/export3/zhangw/Code.sum/some_meaningful_theme_for_ggplot.R')

colnames(dfplot)
dfplot$variable = factor(dfplot$variable,levels = c('KRAS','MRAS','RASA4','RASA4B', 'NGF','CDC42EP3','ARHGEF2', #RAS
                                                    'WNT10B','ADGRA2','SFRP2','FRAT1','TMEM132A','FGF2','BAMBI', # WNT
                                                    "ATF6", "ERN1",  "EIF2AK3",  "XBP1", "HSPA5",  "DDIT3",  "PPP1R15A",  "ATF4" # ER
))

PlotGroup_barplot(data = dfplot,xlab = 'variable',ylab = 'value',legend = 'treatment',color = c('#375631','#3C224B'))+
  zhangwei_theme+
  theme(legend.position = 'top')

ggsave('./fig/fig_8g_barplot_deg_ras_wnt_er.pdf',width = 15,height = 5,dpi = 600)





########7.8validation in qpcr #############

df = read_xlsx('./data/experi/qpcr_validation_mech_1.xlsx')
df = df[,1:4]
source('/export3/zhangw/Code.sum/PlotGroup_barplot.R')
source('/export3/zhangw/Code.sum/some_meaningful_theme_for_ggplot.R')

plist = list()
for (i in c('U251','U87')) {
  
  cellname = i
  # 宽数据转换为长数据格式
  df%>% dplyr::filter(df$celltype==cellname,) %>%
    pivot_longer(cols = 3:4,
                 names_to = "Group",
                 values_to ="Value")->Figure_2b
  
  # 分组标签
  Group<-unique(Figure_2b$Group)
  # 分组因子化
  Figure_2b$Group<-factor(Figure_2b$Group,levels = unique(Figure_2b$Group),labels = c('control', "sh-MEST" ))
  Figure_2b$gene = factor(Figure_2b$gene,levels = unique(Figure_2b$gene))
  my_comparisons <- list(c('control', "sh-MEST" ))
  
  # color = c('#375631','#D3C3D0','#3C224B')
  color = c('#375631','#3C224B')
  xlab =''
  ylab ='Relative mRNA expression'
  orientation = c('vertical','horizontal','reverse')[1]
  p_show = c('p.signif','p.format')[2]
  add = c( 'mean_se','jitter', 'mean_sd', 'mean_ci', 'mean_range', 'median','mean', 'median_iqr')[c(1,2)]
  
  data = Figure_2b[,-1]
  
  names(data) = c('Sample','Group',"Value")
  
  stat.test <- data %>%
    dplyr::group_by(Sample) %>%
    t_test(Value ~ Group) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")
  
  stat.test <- stat.test %>%
    add_xy_position(x = "Sample", fun = "mean_se", dodge = 0.8)
  
  data$Value = as.numeric(unlist(data$Value))
  
  
  
 p1= PlotGroup_barplot(data = data,xlab = 'Sample',color = color)+ scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
  labs(title = cellname)+
    ylab("Relative mRNA level") +
    xlab("") +
    my_theme_zhangwei+
    theme(legend.position = "top")
  plist[[i]] = p1
  
}

pdf('./fig/fig_8d_validation_u251_87_MEST_sh_mRNA.pdf',width = 8,height = 10,onefile = F)
aplot::gglist(plist,nrow=2)
dev.off()

####### choose some meaningful data 

plist  = list()
### U87  KRAS MRAS
cellname = 'U87'
m_gene = c('KRAS','MRAS')
# 宽数据转换为长数据格式
df%>% 
  dplyr::filter(celltype==cellname) %>% 
  dplyr::filter(gene %in% m_gene) %>% 
  pivot_longer(cols = 3:5,
               names_to = "Group",
               values_to ="Value")->Figure_2b
# 分组标签
Group<-unique(Figure_2b$Group)
# 分组因子化
Figure_2b$Group<-factor(Figure_2b$Group,levels = unique(Figure_2b$Group),labels = c('control', "sh-MEST#1" ,"sh-MEST#2"))
Figure_2b$gene = factor(Figure_2b$gene,levels = unique(Figure_2b$gene))
my_comparisons <- list(c("control","sh-MEST#1"),c("control","sh-MEST#2"))

color = c('#375631','#D3C3D0','#3C224B')
xlab =''
ylab ='Relative mRNA expression'
orientation = c('vertical','horizontal','reverse')[1]
p_show = c('p.signif','p.format')[2]
add = c( 'mean_se','jitter', 'mean_sd', 'mean_ci', 'mean_range', 'median','mean', 'median_iqr')[c(1,2)]

data = Figure_2b[,-1]

names(data) = c('Sample','Group',"Value")

stat.test <- data %>%
  dplyr::group_by(Sample) %>%
  t_test(Value ~ Group) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")

stat.test <- stat.test %>%
  add_xy_position(x = "Sample", fun = "mean_se", dodge = 0.8)


p1= ggbarplot(data,
              x = "Sample", 
              y = "Value",
              fill = "Group",
              color = 'Group', alpha=0.5,
              xlab = xlab,
              ylab = ylab,
              add = add,orientation = orientation,
              width = 0.5,
              legend.title = legend,
              palette = color,
              position = position_dodge(0.8),
              ggtheme = theme_bw()) + scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
  labs(title = cellname)+
  ylab("Relative mRNA level") +
  xlab("") +
  zhangwei_theme+
  theme(legend.position = "top")+
  stat_pvalue_manual(stat.test[stat.test$group1=='control',], label = "p", tip.length = 0.01,bracket.nudge.y = 0.2)
p1
plist[[cellname]] =p1

### U251  KRAS MRAS
cellname = 'U251'
m_gene = c('RASA4','RASA4B')
# 宽数据转换为长数据格式
df%>% 
  dplyr::filter(celltype==cellname) %>% 
  dplyr::filter(gene %in% m_gene) %>% 
  pivot_longer(cols = 3:5,
               names_to = "Group",
               values_to ="Value")->Figure_2b
# 分组标签
Group<-unique(Figure_2b$Group)

# 分组因子化
Figure_2b$Group<-factor(Figure_2b$Group,levels = unique(Figure_2b$Group),labels = c('control', "sh-MEST#1" ,"sh-MEST#2"))
Figure_2b$gene = factor(Figure_2b$gene,levels = unique(Figure_2b$gene))
my_comparisons <- list(c("control","sh-MEST#1"),c("control","sh-MEST#2"))

color = c('#375631','#D3C3D0','#3C224B')
xlab =''
ylab ='Relative mRNA expression'
orientation = c('vertical','horizontal','reverse')[1]
p_show = c('p.signif','p.format')[2]
add = c( 'mean_se','jitter', 'mean_sd', 'mean_ci', 'mean_range', 'median','mean', 'median_iqr')[c(1,2)]

data = Figure_2b[,-1]

names(data) = c('Sample','Group',"Value")

stat.test <- data %>%
  dplyr::group_by(Sample) %>%
  t_test(Value ~ Group) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")

stat.test <- stat.test %>%
  add_xy_position(x = "Sample", fun = "mean_se", dodge = 0.8)


p1= ggbarplot(data,
              x = "Sample", 
              y = "Value",
              fill = "Group",
              color = 'Group', alpha=0.5,
              xlab = xlab,
              ylab = ylab,
              add = add,orientation = orientation,
              width = 0.5,
              legend.title = legend,
              palette = color,
              position = position_dodge(0.8),
              ggtheme = theme_bw()) + scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
  labs(title = cellname)+
  ylab("Relative mRNA level") +
  xlab("") +
  zhangwei_theme+
  theme(legend.position = "top")+
  stat_pvalue_manual(stat.test[stat.test$group1=='control',], label = "p", tip.length = 0.01,bracket.nudge.y = 0.2)
p1
plist[[cellname]] =p1

ggarrange(plist[[1]],plist[[2]],common.legend = T ,ncol = 2)
aplot::ggsave('./fig/fig_8d_validation_u251_87_MEST_sh_RAS_pathway_mRNA.pdf',width = 6,height = 4,dpi = 600)


































########## 计算MEST和各个通路之间的相关性 
# rm(list = ls())
gc()

ffiles = list.files('./data/genelist_egfr/')

egfr.gene = lapply(ffiles, function(x){
  a = read.gmt(paste0('./data/genelist_egfr/',x))
  return(a)
}) %>% do.call(rbind,.)
egfr.gene = unique(egfr.gene$gene)

source('/export3/zhangw/Project_Cross/Project.Glioma.SGene/code/analysis.cor.genes.enrich.pathway.glioma.R')
load("/export3/zhangw/Project_Cross/Project_Mime/Proj/data/glioma.cohort.Rdata")
res.cor = lapply(names(list_train_vali_Data), function(x){
  
  tmp = analysis.cor.genes.enrich.pathway.glioma(input_genelist = list(EGFR=egfr.gene),
                                                 method_cor = 'pearson',
                                                 method_enrich = 'gsva',
                                                 kcdf = 'Gaussian',parallel.sz = 24,cancer_type = 'glioma',datasets = x)
  return(tmp)
  
})

names(res.cor) = names(list_train_vali_Data)
# save(res.cor,file = './res/res.cor.egfr.all.gene.glioma.Rdata')

# 将 MEST 提取出来
cor.mest = lapply(res.cor, function(x){
  tmp = x
  tmpt = tmp[tmp$gene_name2 =='MEST',]
  return(tmpt)
}) %>% do.call(rbind,.)

dfplot =cor.mest %>% rownames_to_column('cohort')

plist = list()

dfplot <- data.frame(cohort = cor.mest$dataset.sum,
                     r = cor.mest$cor_r,
                     p = -log10(cor.mest$pvalue),
                     Significance = ifelse(cor.mest$pvalue<0.05,'Sig','NotSig'))

p1 <- ggplot(data = dfplot,aes(r,forcats::fct_reorder(cohort,r,.desc = T))) +
  geom_segment(aes(xend=0,yend=cohort,colour=Significance ),linetype = 2) +
  geom_point(aes(size=p),col = alpha(c('#501d8a'),0.5)) +
  scale_size_continuous(range =c(2,8)) +
  scale_x_continuous(breaks = c(-0.2,0, 0.25, 0.5),
                     expand = expansion(mult = c(0.01,.1))) + #左右留空
  theme_classic() +
  labs(x = "Correlation coefficient", y = "", size = bquote("-log"[10]~"("~italic(P)~"-value)")) + 
  theme(legend.position = "bottom", 
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank())+
  scale_color_manual(values =c('#9E9E9E','#ad2e28'))  # 设置线段颜色为蓝色和红色

p2=p1+xlab(label = '')+ggtitle('MEST versus EGFR signal',subtitle = 'Pearson R with GSVA score')+
  guides(color=guide_legend(nrow = 1),
         size = guide_legend(nrow = 1))

p2
ggsave('./fig/fig_5h_mest_egfr_gsva_cor.pdf',width = 4,height = 6)

plist[[1]] = p2


####### 7.9 所有队列中计算RAS 通路 和 MEST 的相关性####


load("/export3/zhangw/Project_Cross/Project_GliomaRSI/data/cellline_glioma_tpm_expr.Rdata")

load("/export3/zhangw/Project_Cross/Project_Mime/Proj/data/glioma.cohort.Rdata")
res.cor = lapply(names(list_train_vali_Data), function(x){
  
  tmp = analysis.cor.genes.enrich.pathway.glioma(input_genelist = list(EGFR=egfr.gene),
                                                 method_cor = 'pearson',
                                                 method_enrich = 'gsva',
                                                 kcdf = 'Gaussian',parallel.sz = 24,cancer_type = 'glioma',datasets = x)
  return(tmp)
  
})

names(res.cor) = names(list_train_vali_Data)
# save(res.cor,file = './res/res.cor.egfr.all.gene.glioma.Rdata')

# 将 MEST 提取出来
cor.mest = lapply(res.cor, function(x){
  tmp = x
  tmpt = tmp[tmp$gene_name2 =='MEST',]
  return(tmpt)
}) %>% do.call(rbind,.)

dfplot =cor.mest %>% rownames_to_column('cohort')

plist = list()

dfplot <- data.frame(cohort = cor.mest$dataset.sum,
                     r = cor.mest$cor_r,
                     p = -log10(cor.mest$pvalue),
                     Significance = ifelse(cor.mest$pvalue<0.05,'Sig','NotSig'))

p1 <- ggplot(data = dfplot,aes(r,forcats::fct_reorder(cohort,r,.desc = T))) +
  geom_segment(aes(xend=0,yend=cohort,colour=Significance ),linetype = 2) +
  geom_point(aes(size=p),col = alpha(c('#501d8a'),0.5)) +
  scale_size_continuous(range =c(2,8)) +
  scale_x_continuous(breaks = c(-0.2,0, 0.25, 0.5),
                     expand = expansion(mult = c(0.01,.1))) + #左右留空
  theme_classic() +
  labs(x = "Correlation coefficient", y = "", size = bquote("-log"[10]~"("~italic(P)~"-value)")) + 
  theme(legend.position = "bottom", 
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank())+
  scale_color_manual(values =c('#9E9E9E','#ad2e28'))  # 设置线段颜色为蓝色和红色

p2=p1+xlab(label = '')+ggtitle('MEST versus EGFR signal',subtitle = 'Pearson R with GSVA score')+
  guides(color=guide_legend(nrow = 1),
         size = guide_legend(nrow = 1))

p2
ggsave('./fig/fig_5h_mest_egfr_gsva_cor.pdf',width = 4,height = 6)




