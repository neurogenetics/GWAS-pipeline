#!/usr/bin/env Rscript

all_data <- read.table("pca.eigenvec2", header = F)

all_data$V1[all_data$V1 == 0] <- 4
all_data$V1[all_data$V1 == 1] <- 5

attach(all_data)
all_data <- all_data[order(V1),] 

pdf("raw_hapmap_plot.pdf")
plot(all_data$V4,all_data$V5,pch=16,col=all_data$V1,xlab="PC1",ylab="PC2")
legend("bottomleft",col=c("green3","red","blue","cyan"),pch=16,c("Africa","Asia","New samples","Europe"))
dev.off()


euros <- read.table("euro.txt", header = F)
asians <- read.table("asiao.txt", header = F)
africans <- read.table("afrio.txt", header = F)

eurosLowC1 <- mean(euros$V4) - (6*sd(euros$V4))
eurosHighC1 <- mean(euros$V4) + (6*sd(euros$V4))
eurosLowC2 <- mean(euros$V5) - (6*sd(euros$V5))
eurosHighC2 <- mean(euros$V5) + (6*sd(euros$V5))

asiansLowC1 <- mean(asians$V4) - (6*sd(asians$V4))
asiansHighC1 <- mean(asians$V4) + (6*sd(asians$V4))
asiansLowC2 <- mean(asians$V5) - (6*sd(asians$V5))
asiansHighC2 <- mean(asians$V5) + (6*sd(asians$V5))

africansLowC1 <- mean(africans$V4) - (6*sd(africans$V4))
africansHighC1 <- mean(africans$V4) + (6*sd(africans$V4))
africansLowC2 <- mean(africans$V5) - (6*sd(africans$V5))
africansHighC2 <- mean(africans$V5) + (6*sd(africans$V5))


temp2 = all_data[(all_data$V4>=eurosLowC1),]
temp3 = temp2[(temp2$V4<=eurosHighC1),]
temp4 = temp3[(temp3$V5>=eurosLowC2),]
EURO = temp4[(temp4$V5<=eurosHighC2),]

pdf("EURO_confirmation_plot.pdf")
plot(EURO$V4,EURO$V5,pch=16,col=EURO$V1)
legend("topright",col=c("blue","cyan"),pch=16,c("New samples","Europe"))
dev.off()


temp2 = all_data[(all_data$V4>=asiansLowC1),]
temp3 = temp2[(temp2$V4<=asiansHighC1),]
temp4 = temp3[(temp3$V5>=asiansLowC2),]
ASIA = temp4[(temp4$V5<=asiansHighC2),]


pdf("ASIA_confirmation_plot.pdf")
plot(ASIA$V4,ASIA$V5,pch=16,col=ASIA$V1)
legend("topright",col=c("blue","red"),pch=16,c("New samples","Asia"))
dev.off()


temp2 = all_data[(all_data$V4>=africansLowC1),]
temp3 = temp2[(temp2$V4<=africansHighC1),]
temp4 = temp3[(temp3$V5>=africansLowC2),]
AFRICA = temp4[(temp4$V5<=africansHighC2),]


pdf("AFRICA_confirmation_plot.pdf")
plot(AFRICA$V4,AFRICA$V5,pch=16,col=AFRICA$V1)
legend("topleft",col=c("blue","green3"),pch=16,c("New samples","Africa"))
dev.off()

library(dplyr)


main_data2 <- all_data[ ! all_data %in% EURO ]


mixed_race1 <- anti_join(all_data, EURO)
mixed_race2 <- anti_join(mixed_race1, ASIA)
mixed_race3 <- anti_join(mixed_race2, AFRICA)

total <- rbind(euros,asians,africans,mixed_race3)

pdf("mixed_race_confirmation_plot.pdf")
plot(total$V4,total$V5,pch=16,col=total$V1)
legend("topleft",col=c("green3","red","blue","black"),pch=16,c("Africa","Asia","Mixed race samples","Europe"))
dev.off()

final_euro = subset(EURO, EURO$V1 != 5)
final_euro2 =  data.frame(final_euro$V2,final_euro$V3)

final_asia = subset(ASIA, ASIA$V1 != 2)
final_asia2 =  data.frame(final_asia$V2,final_asia$V3)

final_africa = subset(AFRICA, AFRICA$V1 != 3)
final_africa2 =  data.frame(final_africa$V2,final_africa$V3)

mixed_race4 =  data.frame(mixed_race3$V2,mixed_race3$V3)

write.table(final_euro2,file = "PCA_filtered_europeans.txt", quote = FALSE,row.names=F,col.names = F)
write.table(final_asia2,file = "PCA_filtered_asians.txt", quote = FALSE,row.names=F,col.names = F)
write.table(final_africa2,file = "PCA_filtered_africans.txt", quote = FALSE,row.names=F,col.names = F)
write.table(mixed_race4,file = "PCA_filtered_mixed_race.txt", quote = FALSE,row.names=F,col.names = F)


q()
no







