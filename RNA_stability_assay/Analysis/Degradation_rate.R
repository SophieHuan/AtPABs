library(dplyr)
all<-read.table("FPKM_Col_Cordycepin_0_to 8h.txt",header = T)
all[,2:ncol(all)]<-all[,2:ncol(all)]+0.0001
#set comparation time points with 0h
t<-as.data.frame(c(0.5,1,2,4,8))
#Calculate the slope of the regression line in the linear model ln(mRNAt/mRNA0) ~ t.
for (m in c(1:3)){
    l<-20+(m-1)*6
    all[,l]<-log((all[,3+(m-1)*6]/all[,2+(m-1)*6]),exp(1))*(-1)
    all[,l+1]<-log((all[,4+(m-1)*6]/all[,2+(m-1)*6]),exp(1))*(-1)
    all[,l+2]<-log((all[,5+(m-1)*6]/all[,2+(m-1)*6]),exp(1))*(-1)
    all[,l+3]<-log((all[,6+(m-1)*6]/all[,2+(m-1)*6]),exp(1))*(-1)
    all[,l+4]<-log((all[,7+(m-1)*6]/all[,2+(m-1)*6]),exp(1))*(-1)
    all[,l+5]<-c()
    for (i in c(1:nrow(all))){
        data<-c()
        data<-as.data.frame(t(all[i,l:(l+4)]))
        data<-cbind(data,t)
        r<-lm(data[,1]~data$`c(0.5, 1, 2, 4, 8)`)
        all[i,(l+5)]<-r$coefficients[2]
    }
}
all<-rename(all,"C1_rate"="V25")
all<-rename(all,"C2_rate"="V31")
all<-rename(all,"C3_rate"="V37")
degradation<-all[c(1,25,31,37)]
degradation$mean_of_rates<-apply(degradation[,2:4],1,mean)
write.table(degradation,"FPKM_Col_degradation_rate.txt",sep="\t",quote = F,row.names = F)

