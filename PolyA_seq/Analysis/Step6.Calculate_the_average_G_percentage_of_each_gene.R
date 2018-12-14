library("stringr")

# This inputs the file and discards scaffold mapping reads 
all<-read.table("Col_1_mapping_polyA_G_information.txt",header=T)
all<-all[substr(all$gene,1,2)== "AT",]

head(all)
#gene    index   start   mapping end     cDNA_length     distance        cDNA_per        polyA_length    A_per   T_per   C_per   G_per   G_count NotA_count
#AT1G20620|AT1G20620.4   ST-E00243:521:HNFCLCCXY:3:1101:15514:1784       2712    67M     2778    2885    107     0.962911611785095       40      0.95    0       0.05    0       0       2
#AT3G62250|AT3G62250.1   ST-E00243:521:HNFCLCCXY:3:1101:15432:1819       651     57M     707     992     285     0.712701612903226       33      1       0       0       0       0       0
#AT2G01190|AT2G01190.1   ST-E00243:521:HNFCLCCXY:3:1101:10145:1836       2535    49M     2583    2714    131     0.951731761238025       25      1       0       0       0       0       0
#AT1G29930|AT1G29930.1   ST-E00243:521:HNFCLCCXY:3:1101:12824:1836       1129    27M     1155    1230     75     0.9390244               45      1       0       0       0       0       0
#AT3G08610|AT3G08610.1   ST-E00243:521:HNFCLCCXY:3:1101:8704:1907         608    26M      633     653     20     0.9693721               33      1       0       0       0       0       0
#AT1G29930|AT1G29930.1   ST-E00243:521:HNFCLCCXY:3:1101:12895:1924       1129    27M     1155    1230     75     0.9390244               43      1       0       0       0       0       0

# Paste gene_id and mapping site into a new index
all$index<-paste(all$gene,all$end,sep="-")

# This defines the major alternative polyadenylation site windown (+/- 4bp) and filters the associated poly(A) tail of each gene
data<-data.frame(tapply(all$end, all$gene, median))
data$tapply.all.end..all.gene..median.<-round(data$tapply.all.end..all.gene..median.)
data$gene<-rownames(data)
data$index<-paste(data$gene,data$tapply.all.end..all.gene..median.,sep="-")
index<-c()
for (i in c(1:nrow(data))){
  for (l in c(1:4)){
    index_sub1<-paste(data[i,2],(as.numeric(as.character(data[i,1]))+l),sep="-")
    index_sub2<-paste(data[i,2],(as.numeric(as.character(data[i,1]))-l),sep="-")
    index<-c(index,index_sub1,index_sub2)
  }
  index<-c(index,data[i,3])
}
all2<-all[which(all$index %in% index),]
all2$cal<-as.numeric("1")

# This defines the major poly(A) tail length of each gene
data2<-data.frame(tapply(all2$polyA_length, all2$gene, median))
data2$gene<-rownames(data2)
data2$gene_id<-substr(data2$gene,1,9)

# This calculates poly(A) tail and G content information of each gene
data2$polyA_readsnum<-tapply(all2$cal, all2$gene,sum)
data2$APA_site<-tapply(all2$end, all2$gene,median,na.rm=T)
data2$cDNA_percentage<-tapply(all2$cDNA_per, all2$gene,median,na.rm=T)
data2$polyA_length<-tapply(all2$polyA_length, all2$gene,median,na.rm=T)
data2$polyA_length_sum<-tapply(all2$polyA_length, all2$gene,sum,na.rm=T)
all2_sub<-all2[all2$G_count>0,]
data2$G_readsnum<-tapply(all2_sub$cal, all2_sub$gene,sum)
data2$G_reads_percentage<-data2$G_readsnum/data2$polyA_readsnum
data2$G_percentage_includingG_reads<-tapply(all2_sub$G_per, all2_sub$gene,median,na.rm=T)
data2$G_count_sum<-tapply(all2_sub$G_count,all2_sub$gene,sum,na.rm=T)
data2$average_G_content<-data2$G_count_sum/data2$polyA_length_sum

# This discards genes with less than 3 poly(A) tail reads
data2<-data2[data2$polyA_readsnum>=3,]
data2<-data2[!is.na(data2$gene),]
data2<-data2[c(2,3,7,9,13)]
data2<-rename(data2,
write.table(data2,"Col_1_mapping_polyA_G_information_readscut3_Ginfor.txt",sep="\t",row.names = F,quote = F)
