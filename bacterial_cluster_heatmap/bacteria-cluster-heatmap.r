# setwd("Documents/Summer_Job_2013/")

#libraries
library(gplots) #heatmap.2
library(GMD) 	#heatmap.3
library(lattice)
library(shape)

#import data, all columns all rows
WIHS <- read.table(file="neg slow rapid long WIHS analytic 4-4-13.csv",header=TRUE,nrows=582, fill=TRUE,sep=",")

# Change NA values to zeros
WIHS[is.na(WIHS)] <- 0

#separate data
bacteria <-WIHS[,7:106]
id <- as.matrix(WIHS$id, ncol=1)
colnames(id) <- "id"
hivstat <- as.matrix(WIHS$hivstat, ncol=1)
colnames(hivstat) <- "hivstat"

#cluster and dendrogram 
#clustering with euclidean distance and ward method linkage
# write .csv file with groupid, patientid, HIVstatus, bacteria
d <- dist(bacteria, method = "euclidean") # distance matrix
fit <- hclust(d, method="ward") 
plot(fit)
# To display dendogram us command plot(fit) 
rect.hclust(fit, k=6, border="red")
cut_groups_ward<- cutree(fit, k=6)
cut_groups_ward6_id_hiv <- cbind(cut_groups_ward, id, hivstat,bacteria) 

# alternative clustering method
# clustering with euclidean distance and complete linkage
# command already run: d <- dist(bacteria, method = "euclidean") # distance matrix
# write .csv file with groupid, patientid, HIVstatus, bacteria
fit_complete <- hclust(d, method="complete") 
plot(fit_complete)
rect.hclust(fit_complete, k=6, border="red")
cut_groups_complete<- cutree(fit_complete, k=6)
length(cut_groups_ward)
cut_groups_complete6_id_hiv <- cbind(cut_groups_complete, id, hivstat,bacteria) 

#CODE TO CREATE HEATMAP WITH TWO COLOR-CODED COLUMNS ON TOP
cut_groups_ward6_id_hiv <- cbind(cut_groups_ward, id, hivstat,bacteria) 

# Normalize abundance data
abundance <- matrix(numeric(0), nrow=100, ncol=2) 
for(i in 1:100) {
     ammount <- sum(bacteria[colnames(bacteria[i])])
     abundance[i,1] <- colnames(bacteria[i])
     abundance[i,2] <- ammount  
    }
colnames(abundance) <- c("bacteria", "ammount")
newdata <- abundance[,order(ammount)]
top20<- newdata[1:20]

newheat <- cbind(cut_groups_ward6_id_hiv, bacteria)
test <- order(newheat$cut_groups_ward, newheat$hivstat)
test2<- bacteria[test,]
test3 <- order(newheat$cut_groups_ward)
cut_groups <-cut_groups_ward


# Assign color per group for image plotting
testgroupcolor <- matrix(numeric(0), nrow=581, ncol=1) 
for(i in 1:581) {
     if (cut_groups[test3][i]=="1") {
		testgroupcolor[i] <- "red" } 
	 else {
		if (cut_groups[test3][i]=="2") {
			testgroupcolor[i] <- "blue"}  
	 else {
		if (cut_groups[test3][i]=="3") {
			testgroupcolor[i] <- "cyan"} 
	 else {
	 	if (cut_groups[test3][i]=="4") {
			testgroupcolor[i] <- "darkorchid" } 
	 else {
	 	if (cut_groups[test3][i]=="5") {
			testgroupcolor[i] <- "yellow" } 
	 else {
	 	if (cut_groups[test3][i]=="6") {
			testgroupcolor[i] <- "green" } 
}}}}}}

#Assign color by hivstat for image plotting
hivcolor <- matrix(numeric(0), nrow=581, ncol=1) 
for(i in 1:581) { 
	if (newheat$hivstat[test][i]=="0") {
		hivcolor[i] <- "black" }
	else {
		if (newheat$hivstat[test][i]==1){
		hivcolor[i] <-"#EBEBEB"} #light grey
	else {
		if (newheat$hivstat[test][i]=="2"){ 
		hivcolor[i] <-"grey"}
		}
	}
}

#set side colors for heatmap
sidecolors <- cbind(testgroupcolor, hivcolor)
colnames(sidecolors)<-c("group","HIVstat")

mydistance=function(c) {dist(bacteria,method="euclidian")}
mycluster=function(c) {hclust(bacteria,method="ward")}

#Create heatmap using custom heatmap.3 source code
######FINAL IMAGE####
#modified 8-9-13
colnames(sidecolors)<-c("Group","HIV status")
tiff("bacterial_communities.tiff", width = 8, height = 8, units = 'in', res = 300, antialias='gray')
main_title="Bacterial Communities"
heatmap.3(t(test2[,1:20]), hclustfun=mycluster, distfun=mydistance, na.rm = TRUE, scale="none",dendrogram="none", margins=c(5,10),Rowv=FALSE, Colv=FALSE,ColSideColors=sidecolors, symbreaks=FALSE, key=TRUE, symkey=FALSE, density.info="none", trace="none", main=main_title, labCol=FALSE, cexRow=1, col=intpalette(c("red","#FF6000", "#FFBF00", "#DFFF00", "#80FF00", "#20FF00", "#00FF40", "#00FF9F", "cyan", "#009FFF", "#0040FF", "#2000FF", "#8000FF", "#DF00FF"),numcol = 128), NumColSideColors=2, KeyValueName="Relative Abundance")
legend("topright",legend=c("Group1","Group2","Group3","Group4","Group5","Group6","","Negative","Stable","Progressive"),fill=c("red","blue","cyan","darkorchid","yellow","green","white","black","#EBEBEB","grey"), border=FALSE, bty="n", y.intersp = 0.7, cex=0.7)
dev.off()
# now in 'preview' convert to jpeg