#!/usr/bin/Rscript
###修改，count 矫正变为染色体中值的N 倍矫正。
### 函数为 norm_chr_fold
###按嘉宝要求输出图片和结果
#install.packages("extrafont")
#library(extrafont)
##library(hash)
####增加二次矫正模块，直接从合并点开始分析---###
###修改：之前二次校正模块选择染色体名称时有误，现修正
###修改：判断性别时计算的mean_Y时过滤头尾5%数据，且需要与Y染色体seg.mean作比较后再确定（为保持一致性）。
###备注：应客户要求，使用三行的点图作为报告图，且在样本result文件夹下输出三种图以及合并点后bin copyNum的文件。注：因报告图修改，输出图片的名称与默认不一致。时间：20220328。
###修改：(1)binCount为0时赋值1；(2)GC_loess校正系数为负值时赋值1;Map_loess校正系数为负值时赋值mappability; --bin_MAP_correct2()、bin_GC_correct2()、Reduce_segment_correct_GC_MAP_raw()

plotpng1 <- function(fi){
			for (chrom_s in unique(fi$data$chrom))
				{
				band_sele=subset(cytoBand,V1==chrom_s)
				chrom_sele=data.frame(fi$data[fi$data$chrom==chrom_s,])
				name2=paste(name1,chrom_s,"unsele.png",sep="_")
				png(name2,width=3800,height=1200)
				pg=ggplot(chrom_sele,aes(x=maploc,y=Lograte))+geom_point()
				pg=pg+theme(text=element_text(face='bold'),axis.title=element_blank(),axis.text.x=element_blank(),axis.text.y=element_text(size=38,colour="black"),panel.background = element_rect(fill='white', colour='black'))+geom_hline(yintercept=c(log2(1.5),0,log2(0.5)),linetype="dashed",linewidth=1.5)
				segment_sele=subset(fi$output,chrom==chrom_s)
				for (i in 1:nrow(segment_sele))
					{
					pg=pg+annotate("segment",x=segment_sele[i,]$loc.start,xend=segment_sele[i,]$loc.end,y=segment_sele[i,]$seg.mean,yend=segment_sele[i,]$seg.mean,colour="red",size=2)
					}
				for ( i in 1:nrow(band_sele))
				{
				pg=pg+annotate("segment",x=band_sele[i,2],xend=band_sele[i,3],y=-3,yend=-3,colour=band_sele[i,6],size=3)+annotate("text",x=round((band_sele[i,2]+band_sele[i,3])/2,0),y=band_sele[i,7],label=band_sele[i,4],size=6,colour=band_sele[i,6])+annotate("rect",xmin=band_sele[i,2],xmax=band_sele[i,3],ymin=-3,ymax=3,alpha=.1,fill=band_sele[i,6])
				}
				pg=pg+scale_y_continuous(breaks=c(-3,-2,-1,0,1,2,3),labels=c(-3,-2,-1,0,1,2,3))
				print(pg)
				dev.off()
				}
			}

plotp1 <- function(fi){
			for (chrom_s in unique(fi$data$chrom))
				{
				band_sele=subset(cytoBand,V1==chrom_s)
				chrom_sele=data.frame(fi$data[fi$data$chrom==chrom_s,])
				chrom_sele$Lograte=2^chrom_sele$Lograte*2
				name2=paste(name1,chrom_s,"unsele.png",sep="_")
				png(name2,width=3800,height=1000)
				pg=ggplot(chrom_sele,aes(x=maploc,y=Lograte))+geom_point()
				pg=pg+theme(text=element_text(face='bold'),axis.title=element_blank(),axis.text.x=element_blank(),axis.text.y=element_text(size=38,colour="black"),panel.background = element_rect(fill='white', colour='black'))+geom_hline(yintercept=c(1,2,3,4),linetype="dashed",linewidth=1.5)
				segment_sele=subset(fi$output,chrom==chrom_s)
				segment_sele$seg.mean=2^segment_sele$seg.mean*2
				for (i in 1:nrow(segment_sele))
					{
					pg=pg+annotate("segment",x=segment_sele[i,]$loc.start,xend=segment_sele[i,]$loc.end,y=segment_sele[i,]$seg.mean,yend=segment_sele[i,]$seg.mean,colour="red",size=2)
					}
				for ( i in 1:nrow(band_sele))
				{
				pg=pg+annotate("segment",x=band_sele[i,2],xend=band_sele[i,3],y=0,yend=0,colour=band_sele[i,6],size=3)+annotate("text",x=round((band_sele[i,2]+band_sele[i,3])/2,0),y=band_sele[i,7],label=band_sele[i,4],size=8,colour=band_sele[i,6])+annotate("rect",xmin=band_sele[i,2],xmax=band_sele[i,3],ymin=0,ymax=3.8,alpha=.1,fill=band_sele[i,6])
				}
				pg=pg+scale_y_continuous(breaks=c(0,1,2,3,4),labels=c(0,1,2,3,4))
				print(pg)
				dev.off()
				}
			}


plotpdf1 <- function(fi){
			for (chrom_s in unique(fi$data$chrom))
				{
				band_sele=subset(cytoBand,V1==chrom_s)
				chrom_sele=data.frame(fi$data[fi$data$chrom==chrom_s,])
				name2=paste(name1,chrom_s,"new.pdf",sep="_")
				pdf(name2,width=38,height=12)
				pg=ggplot(chrom_sele,aes(x=maploc,y=Lograte))+geom_point()
				pg=pg+theme(text=element_text(face='bold'),axis.title=element_blank(),axis.text.x=element_blank(),axis.text.y=element_text(size=38,colour="black"),panel.background = element_rect(fill='white', colour='black'))+geom_hline(yintercept=c(log2(1.5),0,log2(0.5)),linetype="dashed",linewidth=1.5)
				segment_sele=subset(fi$output,chrom==chrom_s)
				for (i in 1:nrow(segment_sele))
					{
					pg=pg+annotate("segment",x=segment_sele[i,]$loc.start,xend=segment_sele[i,]$loc.end,y=segment_sele[i,]$seg.mean,yend=segment_sele[i,]$seg.mean,colour="red",size=2)
					}
				for ( i in 1:nrow(band_sele))
				{
				pg=pg+annotate("segment",x=band_sele[i,2],xend=band_sele[i,3],y=-3.1,yend=-3.1,colour=band_sele[i,6],size=3)+annotate("text",x=round((band_sele[i,2]+band_sele[i,3])/2,0),y=band_sele[i,7],label=band_sele[i,4],size=6,colour=band_sele[i,6])+annotate("rect",xmin=band_sele[i,2],xmax=band_sele[i,3],ymin=-3,ymax=3,alpha=.1,fill=band_sele[i,6])
				}
				pg=pg+scale_y_continuous(breaks=c(-3,-2,-1,0,1,2,3),labels=c(-3,-2,-1,0,1,2,3))
				print(pg)
				dev.off()
				}
			}

plotpdf2 <- function(fi){
			for (chrom_s in unique(fi$data$chrom))
				{
				band_sele=subset(cytoBand,V1==chrom_s)
				chrom_sele=data.frame(fi$data[fi$data$chrom==chrom_s,])
				name3=paste(name1,chrom_s,"sd2.pdf",sep="_")
				pdf(name3,width=38,height=12)
				pg=ggplot(chrom_sele,aes(x=maploc,y=Lograte))+geom_point()
				pg=pg+theme(text=element_text(face='bold'),axis.title=element_blank(),axis.text.x=element_blank(),axis.text.y=element_text(size=38,colour="black"),panel.background = element_rect(fill='white', colour='black'))+geom_hline(yintercept=c(log2(1.5),0,log2(0.5)),linetype="dashed",linewidth=1.5)
				segment_sele=subset(fi$output,chrom==chrom_s)
				for (i in 1:nrow(segment_sele))
					{
					pg=pg+annotate("segment",x=segment_sele[i,]$loc.start,xend=segment_sele[i,]$loc.end,y=segment_sele[i,]$seg.mean,yend=segment_sele[i,]$seg.mean,colour="red",size=2)
					}
				for ( i in 1:nrow(band_sele))
				{
				pg=pg+annotate("segment",x=band_sele[i,2],xend=band_sele[i,3],y=-3.1,yend=-3.1,colour=band_sele[i,6],size=3)+annotate("text",x=round((band_sele[i,2]+band_sele[i,3])/2,0),y=band_sele[i,7],label=band_sele[i,4],size=6,colour=band_sele[i,6])+annotate("rect",xmin=band_sele[i,2],xmax=band_sele[i,3],ymin=-3,ymax=3,alpha=.1,fill=band_sele[i,6])
				}
				pg=pg+scale_y_continuous(breaks=c(-3,-2,-1,0,1,2,3),labels=c(-3,-2,-1,0,1,2,3))
				print(pg)
				dev.off()
				}
			}

plotpng2 <- function(fi){
			for (chrom_s in unique(fi$data$chrom))
				{
				band_sele=subset(cytoBand,V1==chrom_s)
				chrom_sele=data.frame(fi$data[fi$data$chrom==chrom_s,])
				name3=paste(name1,chrom_s,"cutsd2.png",sep="_")
				png(name3,width=3800,height=1200)
				pg=ggplot(chrom_sele,aes(x=maploc,y=Lograte))+geom_point()
				pg=pg+theme(text=element_text(face='bold'),axis.title=element_blank(),axis.text.x=element_blank(),axis.text.y=element_text(size=38,colour="black"),panel.background = element_rect(fill='white', colour='black'))+geom_hline(yintercept=c(log2(1.5),0,log2(0.5)),linetype="dashed",linewidth=1.5)
				segment_sele=subset(fi$output,chrom==chrom_s)
				for (i in 1:nrow(segment_sele))
					{
					pg=pg+annotate("segment",x=segment_sele[i,]$loc.start,xend=segment_sele[i,]$loc.end,y=segment_sele[i,]$seg.mean,yend=segment_sele[i,]$seg.mean,colour="red",size=2)
					}
				for ( i in 1:nrow(band_sele))
				{
				pg=pg+annotate("segment",x=band_sele[i,2],xend=band_sele[i,3],y=-3,yend=-3,colour=band_sele[i,6],size=3)+annotate("text",x=round((band_sele[i,2]+band_sele[i,3])/2,0),y=band_sele[i,7],label=band_sele[i,4],size=6,colour=band_sele[i,6])+annotate("rect",xmin=band_sele[i,2],xmax=band_sele[i,3],ymin=-3,ymax=3,alpha=.1,fill=band_sele[i,6])
				}
				pg=pg+scale_y_continuous(breaks=c(-3,-2,-1,0,1,2,3),labels=c(-3,-2,-1,0,1,2,3))
				print(pg)
				dev.off()
				}
			}

plotp2<- function(fi,plot_type){
  for (chrom_s in unique(fi$data$chrom))
  {
	band_sele=subset(cytoBand,V1==chrom_s)
	chrom_sele=data.frame(fi$data[fi$data$chrom==chrom_s,])
	chrom_sele$Lograte=2^chrom_sele$Lograte*2
	name3=paste(name1,chrom_s,plot_type,sep="_")
	png(name3,width=3800,height=1000)
	pg=ggplot(chrom_sele,aes(x=maploc,y=Lograte))+geom_point()
	pg=pg+theme(text=element_text(face='bold'),axis.title=element_blank(),axis.text.x=element_blank(),axis.text.y=element_text(size=38,colour="black"),panel.background = element_rect(fill='white', colour='black'))+geom_hline(yintercept=c(1,2,3,4),linetype="dashed",linewidth=1.5)
	segment_sele=subset(fi$output,chrom==chrom_s)
	segment_sele$seg.mean=2^segment_sele$seg.mean*2
	for (i in 1:nrow(segment_sele))
	{
	pg=pg+annotate("segment",x=segment_sele[i,]$loc.start,xend=segment_sele[i,]$loc.end,y=segment_sele[i,]$seg.mean,yend=segment_sele[i,]$seg.mean,colour="red",size=2)
	}
	for ( i in 1:nrow(band_sele))
	{
	  pg=pg+annotate("segment",x=band_sele[i,2],xend=band_sele[i,3],y=0,yend=0,colour=band_sele[i,6],size=3)+annotate("text",x=round((band_sele[i,2]+band_sele[i,3])/2,0),y=band_sele[i,7],label=band_sele[i,4],size=8,colour=band_sele[i,6])+annotate("rect",xmin=band_sele[i,2],xmax=band_sele[i,3],ymin=0,ymax=3.8,alpha=.1,fill=band_sele[i,6])
	}
	pg=pg+scale_y_continuous(breaks=c(0,1,2,3,4),labels=c(0,1,2,3,4))
	print(pg)
	dev.off()
  }
}



plotpng3 <- function(fi){
			for (chrom_s in unique(fi$data$chrom))
				{
				band_sele=subset(cytoBand,V1==chrom_s)
				chrom_sele=data.frame(fi$data[fi$data$chrom==chrom_s,])
				name3=paste(name1,chrom_s,"sd3n.png",sep="_")
				png(name3,width=3800,height=1200)
				pg=ggplot(chrom_sele,aes(x=maploc,y=Lograte))+geom_point()
				pg=pg+theme(text=element_text(face='bold'),axis.title=element_blank(),axis.text.x=element_blank(),axis.text.y=element_text(size=38,colour="black"),panel.background = element_rect(fill='white', colour='black'))+geom_hline(yintercept=c(log2(1.5),0,log2(0.5)),linetype="dashed",linewidth=1.5)
				segment_sele=subset(fi$output,chrom==chrom_s)
				for (i in 1:nrow(segment_sele))
					{
					pg=pg+annotate("segment",x=segment_sele[i,]$loc.start,xend=segment_sele[i,]$loc.end,y=segment_sele[i,]$seg.mean,yend=segment_sele[i,]$seg.mean,colour="red",size=2)
					}
				for ( i in 1:nrow(band_sele))
				{
				pg=pg+annotate("segment",x=band_sele[i,2],xend=band_sele[i,3],y=-3,yend=-3,colour=band_sele[i,6],size=3)+annotate("text",x=round((band_sele[i,2]+band_sele[i,3])/2,0),y=band_sele[i,7],label=band_sele[i,4],size=6,colour=band_sele[i,6])+annotate("rect",xmin=band_sele[i,2],xmax=band_sele[i,3],ymin=-3,ymax=3,alpha=.1,fill=band_sele[i,6])
				}
				pg=pg+scale_y_continuous(breaks=c(-3,-2,-1,0,1,2,3),labels=c(-3,-2,-1,0,1,2,3))
				print(pg)
				dev.off()
				}
			}

plot_merge_chr <- function(fi){
	library(gridExtra)
	fi=data.frame(fi) ###chrom maploc Lograte
	fi$Lograte=2^fi$Lograte
	Group=c(paste("chr",seq(1,22),sep=""),"chrX","chrY")
	color_type=as.character(rep(c('#009900','#FFCC00','#990099','#3333CC'),6))
	color_type1=as.character(seq(0.01,1,0.01))
	
	ALL_chr1=c()
	star_ADD=c()
	X_lab=c()
	X_lab_pos=c()
	for (i in 1:5)
		{
		if(length(ALL_chr1)<2)
			{star_add=0}
		else{star_add=ALL_chr1[nrow(ALL_chr1),5]}
		star_ADD=c(star_ADD,star_add)
		sele_chr=subset(fi,chrom==Group[i])
		sele_chr$COL=color_type1[i]
		x_lab=c(seq(30000000,sele_chr[nrow(sele_chr),2],30000000))
		X_lab=c(X_lab,x_lab)
		X_lab_pos=c(X_lab_pos,x_lab+star_add)
		sele_chr$Xpos=sele_chr$maploc+star_add
		ALL_chr1=rbind(ALL_chr1,sele_chr)
		}

	star_ADD=c(star_ADD,ALL_chr1[nrow(ALL_chr1),5]+10000000)
	pg1=ggplot(ALL_chr1,aes(x=Xpos,y=Lograte))+geom_point(aes(color=COL,size=1))+scale_x_continuous(expand = c(0,0),breaks = X_lab_pos,labels = X_lab/1000000)+scale_y_continuous(expand = c(0,0),limits = c(0,4.5),breaks = c(0,1,2,3,4),labels = paste(seq(0,4),"X",sep=" "))+labs(x="",y="")+geom_hline(yintercept=c(0,1,2,3,4),linewidth=0.5)
	pg1=pg1+theme_classic()+theme(legend.position = "none")+theme(axis.text=element_text(size=38,colour="black"))+scale_color_manual(values =color_type[1:5])+annotate("segment",x=0,xend=0,y=0,yend=4.5,colour="red",size=3)
	for (i in 1:5)
		{
		pg1=pg1+annotate("segment",x=star_ADD[i],xend=star_ADD[i+1],y=0,yend=0,colour=color_type[i],size=3)+annotate("text",x=star_ADD[i]+15000000,y=4.2,label=Group[i],size=15,colour=color_type[i])
		}
	ALL_chr1=c()
	star_ADD=c()
	X_lab=c()
	X_lab_pos=c()
	for (i in 6:12)
		{
		if(length(ALL_chr1)<2)
			{star_add=0}
		else{star_add=ALL_chr1[nrow(ALL_chr1),5]}
		star_ADD=c(star_ADD,star_add)
		sele_chr=subset(fi,chrom==Group[i])
		sele_chr$COL=color_type1[i]
		
		x_lab=c(seq(30000000,sele_chr[nrow(sele_chr),2],30000000))
		X_lab=c(X_lab,x_lab)
		X_lab_pos=c(X_lab_pos,x_lab+star_add)
		sele_chr$Xpos=sele_chr$maploc+star_add
		ALL_chr1=rbind(ALL_chr1,sele_chr)
		}
	star_ADD=c(star_ADD,ALL_chr1[nrow(ALL_chr1),5])
	pg2=ggplot(ALL_chr1,aes(x=Xpos,y=Lograte))+geom_point(aes(color=COL,size=1))+scale_x_continuous(expand = c(0,0),breaks = X_lab_pos,labels = X_lab/1000000)+scale_y_continuous(expand = c(0,0),limits = c(0,4.5),breaks = c(0,1,2,3,4),labels = paste(seq(0,4),"X",sep=" "))+labs(x="",y="")+geom_hline(yintercept=c(0,1,2,3,4),linewidth=0.5)
	pg2=pg2+theme_classic()+theme(legend.position = "none")+theme(axis.text=element_text(size=38,colour="black"))+scale_color_manual(values =color_type[6:12])+annotate("segment",x=0,xend=0,y=0,yend=4.5,colour="red",size=3)
	for (i in 1:7)
		{
		pg2=pg2+annotate("segment",x=star_ADD[i],xend=star_ADD[i+1],y=0,yend=0,colour=color_type[i+5],size=3)+annotate("text",x=star_ADD[i]+15000000,y=4.2,label=Group[i+5],size=15,colour=color_type[i+5])
		}
		
	ALL_chr1=c()
	star_ADD=c()
	X_lab=c()
	X_lab_pos=c()
	for (i in 13:22)
		{
		if(length(ALL_chr1)<2)
			{star_add=0}
		else{star_add=ALL_chr1[nrow(ALL_chr1),5]}
		star_ADD=c(star_ADD,star_add)
		sele_chr=subset(fi,chrom==Group[i])
		sele_chr$COL=color_type1[i]
		x_lab=c(seq(30000000,sele_chr[nrow(sele_chr),2],30000000))
		X_lab=c(X_lab,x_lab)
		X_lab_pos=c(X_lab_pos,x_lab+star_add)
		sele_chr$Xpos=sele_chr$maploc+star_add
		ALL_chr1=rbind(ALL_chr1,sele_chr)
		}
	star_ADD=c(star_ADD,ALL_chr1[nrow(ALL_chr1),5])
	star_ADD=c(star_ADD,star_ADD[length(star_ADD)]+155270560,star_ADD[length(star_ADD)]+155270560+59373566)
	pg3=ggplot(ALL_chr1,aes(x=Xpos,y=Lograte))+geom_point(aes(color=COL,size=1))+scale_x_continuous(expand = c(0,0),breaks = X_lab_pos,labels = X_lab/1000000)+scale_y_continuous(expand = c(0,0),limits = c(0,4.5),breaks = c(0,1,2,3,4),labels = paste(seq(0,4),"X",sep=" "))+labs(x="",y="")+geom_hline(yintercept=c(0,1,2,3,4),linewidth=0.5)
	pg3=pg3+theme_classic()+theme(legend.position = "none")+theme(axis.text=element_text(size=38,colour="black"))+scale_color_manual(values =color_type[13:22])+annotate("segment",x=0,xend=0,y=0,yend=4.5,colour="red",size=3)
	for (i in 1:12)
		{
		pg3=pg3+annotate("segment",x=star_ADD[i],xend=star_ADD[i+1],y=0,yend=0,colour=color_type[i+12],size=3)+annotate("text",x=star_ADD[i]+15000000,y=4.2,label=Group[i+12],size=15,colour=color_type[i+12])
		}
	png("merge_all_chrom.png",width=4000,height=1800)
	pa=grid.arrange(pg1,pg2,pg3, ncol=1, nrow=3)
	print(pa)
	dev.off()
	}


####

plot_merge_chr2 <- function(fi,name1){
	library(gridExtra)
	fi=data.frame(fi) ###chrom maploc Lograte
	fi$Lograte=2^fi$Lograte*2
	Group=c(paste("chr",seq(1,22),sep=""),"chrX","chrY")
	color_type=as.character(rep(c('#009900','#FFCC00','#990099','#3333CC'),6))
	# color_type=as.character(rep(c("red","blue"),12))
	color_type1=as.character(seq(0.01,1,0.01))
	
	ALL_chr1=c()
	star_ADD=c()
	X_lab=c()
	X_lab_pos=c()
	for (i in 1:5)
		{
		if(length(ALL_chr1)<2)
			{star_add=0}
		else{star_add=ALL_chr1[nrow(ALL_chr1),5]}
		star_ADD=c(star_ADD,star_add)
		sele_chr=subset(fi,chrom==Group[i])
		sele_chr$COL=color_type1[i]
		x_lab=c(seq(30000000,sele_chr[nrow(sele_chr),2],30000000))
		X_lab=c(X_lab,x_lab)
		X_lab_pos=c(X_lab_pos,x_lab+star_add)
		sele_chr$Xpos=sele_chr$maploc+star_add
		ALL_chr1=rbind(ALL_chr1,sele_chr)
		}

	star_ADD=c(star_ADD,ALL_chr1[nrow(ALL_chr1),5]+10000000)
	pg1=ggplot(ALL_chr1,aes(x=Xpos,y=Lograte))+geom_point(aes(color=COL,size=1))+scale_x_continuous(expand = c(0,0),breaks = X_lab_pos,labels = X_lab/1000000)+scale_y_continuous(expand = c(0,0),limits = c(0,4.5),breaks = c(0,1,2,3,4),labels = paste(seq(0,4),"X",sep=" "))+labs(x="",y="")+geom_hline(yintercept=c(0,1,2,3,4),linewidth=0.5)
	pg1=pg1+theme_classic()+theme(legend.position = "none")+theme(axis.text=element_text(size=38,colour="black"))+scale_color_manual(values =color_type[1:5])+annotate("segment",x=0,xend=0,y=0,yend=4.5,colour="red",size=3)
	for (i in 1:5)
		{
		pg1=pg1+annotate("segment",x=star_ADD[i],xend=star_ADD[i+1],y=0,yend=0,colour=color_type[i],size=3)+annotate("text",x=star_ADD[i]+15000000,y=4.2,label=Group[i],size=15,colour=color_type[i])
		}
	ALL_chr1=c()
	star_ADD=c()
	X_lab=c()
	X_lab_pos=c()
	for (i in 6:12)
		{
		if(length(ALL_chr1)<2)
			{star_add=0}
		else{star_add=ALL_chr1[nrow(ALL_chr1),5]}
		star_ADD=c(star_ADD,star_add)
		sele_chr=subset(fi,chrom==Group[i])
		sele_chr$COL=color_type1[i]
		
		x_lab=c(seq(30000000,sele_chr[nrow(sele_chr),2],30000000))
		X_lab=c(X_lab,x_lab)
		X_lab_pos=c(X_lab_pos,x_lab+star_add)
		sele_chr$Xpos=sele_chr$maploc+star_add
		ALL_chr1=rbind(ALL_chr1,sele_chr)
		}
	star_ADD=c(star_ADD,ALL_chr1[nrow(ALL_chr1),5])
	pg2=ggplot(ALL_chr1,aes(x=Xpos,y=Lograte))+geom_point(aes(color=COL,size=1))+scale_x_continuous(expand = c(0,0),breaks = X_lab_pos,labels = X_lab/1000000)+scale_y_continuous(expand = c(0,0),limits = c(0,4.5),breaks = c(0,1,2,3,4),labels = paste(seq(0,4),"X",sep=" "))+labs(x="",y="")+geom_hline(yintercept=c(0,1,2,3,4),linewidth=0.5)
	pg2=pg2+theme_classic()+theme(legend.position = "none")+theme(axis.text=element_text(size=38,colour="black"))+scale_color_manual(values =color_type[6:12])+annotate("segment",x=0,xend=0,y=0,yend=4.5,colour="red",size=3)
	for (i in 1:7)
		{
		pg2=pg2+annotate("segment",x=star_ADD[i],xend=star_ADD[i+1],y=0,yend=0,colour=color_type[i+5],size=3)+annotate("text",x=star_ADD[i]+15000000,y=4.2,label=Group[i+5],size=15,colour=color_type[i+5])
		}

	ALL_chr1=c()
	star_ADD=c()
	X_lab=c()
	X_lab_pos=c()
	for (i in 13:22)
		{
		if(length(ALL_chr1)<2)
			{star_add=0}
		else{star_add=ALL_chr1[nrow(ALL_chr1),5]}
		star_ADD=c(star_ADD,star_add)
		sele_chr=subset(fi,chrom==Group[i])
		sele_chr$COL=color_type1[i]
		x_lab=c(seq(30000000,sele_chr[nrow(sele_chr),2],30000000))
		X_lab=c(X_lab,x_lab)
		X_lab_pos=c(X_lab_pos,x_lab+star_add)
		sele_chr$Xpos=sele_chr$maploc+star_add
		ALL_chr1=rbind(ALL_chr1,sele_chr)
		}
	star_ADD=c(star_ADD,ALL_chr1[nrow(ALL_chr1),5])
	star_ADD=c(star_ADD,star_ADD[length(star_ADD)]+155270560,star_ADD[length(star_ADD)]+155270560+59373566)
	pg3=ggplot(ALL_chr1,aes(x=Xpos,y=Lograte))+geom_point(aes(color=COL,size=1))+scale_x_continuous(expand = c(0,0),breaks = X_lab_pos,labels = X_lab/1000000)+scale_y_continuous(expand = c(0,0),limits = c(0,4.5),breaks = c(0,1,2,3,4),labels = paste(seq(0,4),"X",sep=" "))+labs(x="",y="")+geom_hline(yintercept=c(0,1,2,3,4),linewidth=0.5)
	pg3=pg3+theme_classic()+theme(legend.position = "none")+theme(axis.text=element_text(size=38,colour="black"))+scale_color_manual(values =color_type[13:22])+annotate("segment",x=0,xend=0,y=0,yend=4.5,colour="red",size=3)
	for (i in 1:12)
		{
		pg3=pg3+annotate("segment",x=star_ADD[i],xend=star_ADD[i+1],y=0,yend=0,colour=color_type[i+12],size=3)+annotate("text",x=star_ADD[i]+15000000,y=4.2,label=Group[i+12],size=15,colour=color_type[i+12])
		}

####---
	ALL_chr1=c()
	star_ADD=c()
	X_lab=c()
	X_lab_pos=c()
	for (i in 13:23)
		{
		if(length(ALL_chr1)<2)
			{star_add=0}
		else{star_add=ALL_chr1[nrow(ALL_chr1),5]}
		star_ADD=c(star_ADD,star_add)
		sele_chr=subset(fi,chrom==Group[i])
		sele_chr$COL=color_type1[i]
		x_lab=c(seq(30000000,sele_chr[nrow(sele_chr),2],30000000))
		X_lab=c(X_lab,x_lab)
		X_lab_pos=c(X_lab_pos,x_lab+star_add)
		sele_chr$Xpos=sele_chr$maploc+star_add
		ALL_chr1=rbind(ALL_chr1,sele_chr)
		}
	
	chrY=sele_chr=subset(fi,chrom=="chrY")
	if (mean(chrY$Lograte)<0.3) ### 修改，chrY小于30% 就不画了，认为缺失。
		{
		star_ADD=c(star_ADD,ALL_chr1[nrow(ALL_chr1),5])
		star_ADD=c(star_ADD,star_ADD[length(star_ADD)]+59373566)
		}else{
		chrY$COL=color_type1[24]
		star_ADD=c(star_ADD,ALL_chr1[nrow(ALL_chr1),5])
		chrY$Xpos=chrY$maploc+ALL_chr1[nrow(ALL_chr1),5]
		ALL_chr1=rbind(ALL_chr1,chrY)
		star_ADD=c(star_ADD,star_ADD[length(star_ADD)]+59373566)
			}

	pg4=ggplot(ALL_chr1,aes(x=Xpos,y=Lograte))+geom_point(aes(color=COL,size=1))+scale_x_continuous(expand = c(0,0),breaks = X_lab_pos,labels = X_lab/1000000)+scale_y_continuous(expand = c(0,0),limits = c(0,4.5),breaks = c(0,1,2,3,4),labels = paste(seq(0,4),"X",sep=" "))+labs(x="",y="")+geom_hline(yintercept=c(0,1,2,3,4),linewidth=0.5)
	pg4=pg4+theme_classic()+theme(legend.position = "none")+theme(axis.text=element_text(size=38,colour="black"))+scale_color_manual(values =color_type[13:24])+annotate("segment",x=0,xend=0,y=0,yend=4.5,colour="red",size=3)
	for (i in 1:14)
		{
		pg4=pg4+annotate("segment",x=star_ADD[i],xend=star_ADD[i+1],y=0,yend=0,colour=color_type[i+12],size=3)+annotate("text",x=star_ADD[i]+15000000,y=4.2,label=Group[i+12],size=15,colour=color_type[i+12])
		}
	
	# png(paste(name1,"merge_all_chrom.png",sep="_"),width=4000,height=1800)
	# pa=grid.arrange(pg1,pg2,pg4, ncol=1, nrow=3)
	# print(pa)
	# dev.off()
	
	# png(paste(name1,"merge_autosome.png",sep="_"),width=4000,height=1800)
	# pa2=grid.arrange(pg1,pg2,pg3, ncol=1, nrow=3)
	# print(pa2)
	# dev.off()
	
	}


####

plot_merge_chr2_new <- function(fi,name1){
	library(gridExtra)
	fi=data.frame(fi) ###chrom maploc Lograte
	fi$Lograte=2^fi$Lograte*2
	Group=c(paste("chr",seq(1,22),sep=""),"chrX","chrY")
	color_type=as.character(rep(c('#009900','#FFCC00','#990099','#3333CC'),6))
	color_type1=as.character(seq(0.01,1,0.01))
	
	ALL_chr1=c()
	star_ADD=c()
	X_lab=c()
	X_lab_pos=c()
	for (i in 1:5)
		{
		if(length(ALL_chr1)<2)
			{star_add=0}
		else{star_add=ALL_chr1[nrow(ALL_chr1),5]}
		star_ADD=c(star_ADD,star_add)
		sele_chr=subset(fi,chrom==Group[i])
		sele_chr$COL=color_type1[i]
		x_lab=c(seq(30000000,sele_chr[nrow(sele_chr),2],30000000))
		X_lab=c(X_lab,x_lab)
		X_lab_pos=c(X_lab_pos,x_lab+star_add)
		sele_chr$Xpos=sele_chr$maploc+star_add
		ALL_chr1=rbind(ALL_chr1,sele_chr)
		}

	star_ADD=c(star_ADD,ALL_chr1[nrow(ALL_chr1),5]+10000000)
	pg1=ggplot(ALL_chr1,aes(x=Xpos,y=Lograte))+geom_point(aes(color=COL,size=1))+scale_x_continuous(expand = c(0,0),breaks = X_lab_pos,labels = X_lab/1000000)+scale_y_continuous(expand = c(0,0),limits = c(0,5),breaks = c(0,1,2,3,4),labels = paste(seq(0,4),"X",sep=" "))+labs(x="",y="")+geom_hline(yintercept=c(0,1,2,3,4),linewidth=0.5)
	pg1=pg1+theme_classic()+theme(legend.position = "none")+theme(axis.text=element_text(size=38,colour="black"))+scale_color_manual(values =color_type[1:5])+annotate("segment",x=0,xend=0,y=0,yend=5,colour="red",size=3)
	for (i in 1:5)
		{
		pg1=pg1+annotate("segment",x=star_ADD[i],xend=star_ADD[i+1],y=0,yend=0,colour=color_type[i],size=3)+annotate("text",x=star_ADD[i]+15000000,y=4.7,label=Group[i],size=15,colour=color_type[i])
		}
	ALL_chr1=c()
	star_ADD=c()
	X_lab=c()
	X_lab_pos=c()
	for (i in 6:12)
		{
		if(length(ALL_chr1)<2)
			{star_add=0}
		else{star_add=ALL_chr1[nrow(ALL_chr1),5]}
		star_ADD=c(star_ADD,star_add)
		sele_chr=subset(fi,chrom==Group[i])
		sele_chr$COL=color_type1[i]
		
		x_lab=c(seq(30000000,sele_chr[nrow(sele_chr),2],30000000))
		X_lab=c(X_lab,x_lab)
		X_lab_pos=c(X_lab_pos,x_lab+star_add)
		sele_chr$Xpos=sele_chr$maploc+star_add
		ALL_chr1=rbind(ALL_chr1,sele_chr)
		}
	star_ADD=c(star_ADD,ALL_chr1[nrow(ALL_chr1),5])
	pg2=ggplot(ALL_chr1,aes(x=Xpos,y=Lograte))+geom_point(aes(color=COL,size=1))+scale_x_continuous(expand = c(0,0),breaks = X_lab_pos,labels = X_lab/1000000)+scale_y_continuous(expand = c(0,0),limits = c(0,5),breaks = c(0,1,2,3,4),labels = paste(seq(0,4),"X",sep=" "))+labs(x="",y="")+geom_hline(yintercept=c(0,1,2,3,4),linewidth=0.5)
	pg2=pg2+theme_classic()+theme(legend.position = "none")+theme(axis.text=element_text(size=38,colour="black"))+scale_color_manual(values =color_type[6:12])+annotate("segment",x=0,xend=0,y=0,yend=5,colour="red",size=3)
	for (i in 1:7)
		{
		pg2=pg2+annotate("segment",x=star_ADD[i],xend=star_ADD[i+1],y=0,yend=0,colour=color_type[i+5],size=3)+annotate("text",x=star_ADD[i]+15000000,y=4.7,label=Group[i+5],size=15,colour=color_type[i+5])
		}

	ALL_chr1=c()
	star_ADD=c()
	X_lab=c()
	X_lab_pos=c()
	for (i in 13:22)
		{
		if(length(ALL_chr1)<2)
			{star_add=0}
		else{star_add=ALL_chr1[nrow(ALL_chr1),5]}
		star_ADD=c(star_ADD,star_add)
		sele_chr=subset(fi,chrom==Group[i])
		sele_chr$COL=color_type1[i]
		x_lab=c(seq(30000000,sele_chr[nrow(sele_chr),2],30000000))
		X_lab=c(X_lab,x_lab)
		X_lab_pos=c(X_lab_pos,x_lab+star_add)
		sele_chr$Xpos=sele_chr$maploc+star_add
		ALL_chr1=rbind(ALL_chr1,sele_chr)
		}
	star_ADD=c(star_ADD,ALL_chr1[nrow(ALL_chr1),5])
	star_ADD=c(star_ADD,star_ADD[length(star_ADD)]+155270560,star_ADD[length(star_ADD)]+155270560+59373566)
	pg3=ggplot(ALL_chr1,aes(x=Xpos,y=Lograte))+geom_point(aes(color=COL,size=1))+scale_x_continuous(expand = c(0,0),breaks = X_lab_pos,labels = X_lab/1000000)+scale_y_continuous(expand = c(0,0),limits = c(0,5),breaks = c(0,1,2,3,4),labels = paste(seq(0,4),"X",sep=" "))+labs(x="",y="")+geom_hline(yintercept=c(0,1,2,3,4),linewidth=0.5)
	pg3=pg3+theme_classic()+theme(legend.position = "none")+theme(axis.text=element_text(size=38,colour="black"))+scale_color_manual(values =color_type[13:22])+annotate("segment",x=0,xend=0,y=0,yend=5,colour="red",size=3)
	for (i in 1:12)
		{
		pg3=pg3+annotate("segment",x=star_ADD[i],xend=star_ADD[i+1],y=0,yend=0,colour=color_type[i+12],size=3)+annotate("text",x=star_ADD[i]+15000000,y=4.7,label=Group[i+12],size=15,colour=color_type[i+12])
		}

####---
	ALL_chr1=c()
	star_ADD=c()
	X_lab=c()
	X_lab_pos=c()
	for (i in 13:23)
		{
		if(length(ALL_chr1)<2)
			{star_add=0}
		else{star_add=ALL_chr1[nrow(ALL_chr1),5]}
		star_ADD=c(star_ADD,star_add)
		sele_chr=subset(fi,chrom==Group[i])
		sele_chr$COL=color_type1[i]
		x_lab=c(seq(30000000,sele_chr[nrow(sele_chr),2],30000000))
		X_lab=c(X_lab,x_lab)
		X_lab_pos=c(X_lab_pos,x_lab+star_add)
		sele_chr$Xpos=sele_chr$maploc+star_add
		ALL_chr1=rbind(ALL_chr1,sele_chr)
		}
	
	chrY=sele_chr=subset(fi,chrom=="chrY")
	if (mean(chrY$Lograte)<0.3) ### 修改，chrY小于30% 就不画了，认为缺失。
		{
		star_ADD=c(star_ADD,ALL_chr1[nrow(ALL_chr1),5])
		star_ADD=c(star_ADD,star_ADD[length(star_ADD)]+59373566)
		}else{
		chrY$COL=color_type1[24]
		star_ADD=c(star_ADD,ALL_chr1[nrow(ALL_chr1),5])
		chrY$Xpos=chrY$maploc+ALL_chr1[nrow(ALL_chr1),5]
		ALL_chr1=rbind(ALL_chr1,chrY)
		star_ADD=c(star_ADD,star_ADD[length(star_ADD)]+59373566)
			}

	pg4=ggplot(ALL_chr1,aes(x=Xpos,y=Lograte))+geom_point(aes(color=COL,size=1))+scale_x_continuous(expand = c(0,0),breaks = X_lab_pos,labels = X_lab/1000000)+scale_y_continuous(expand = c(0,0),limits = c(0,5),breaks = c(0,1,2,3,4),labels = paste(seq(0,4),"X",sep=" "))+labs(x="",y="")+geom_hline(yintercept=c(0,1,2,3,4),linewidth=0.5)
	pg4=pg4+theme_classic()+theme(legend.position = "none")+theme(axis.text=element_text(size=38,colour="black"))+scale_color_manual(values =color_type[13:24])+annotate("segment",x=0,xend=0,y=0,yend=5,colour="red",size=3)
	for (i in 1:14)
		{
		pg4=pg4+annotate("segment",x=star_ADD[i],xend=star_ADD[i+1],y=0,yend=0,colour=color_type[i+12],size=3)+annotate("text",x=star_ADD[i]+15000000,y=4.7,label=Group[i+12],size=15,colour=color_type[i+12])
		}
	
	# png(paste(name1,"merge_all_chrom.png",sep="_"),width=4000,height=1800)
	# pa=grid.arrange(pg1,pg2,pg4, ncol=1, nrow=3)
	# print(pa)
	# dev.off()
	
	# png(paste(name1,"merge_autosome.png",sep="_"),width=4000,height=1800)
	# pa2=grid.arrange(pg1,pg2,pg3, ncol=1, nrow=3)
	# print(pa2)
	# dev.off()
	
	png(paste(name1,"merge_all_chrom_new.png",sep="_"),width=4000,height=1800)
	pa=grid.arrange(pg1,pg2,pg4, ncol=1, nrow=3)
	print(pa)
	dev.off()
	
	}


####

plot_merge_chr2_new2 <- function(fi,name1){
	library(gridExtra)
	fi=data.frame(fi) ###chrom maploc Lograte
	fi$Lograte=2^fi$Lograte*2
	Group=c(paste("chr",seq(1,22),sep=""),"chrX","chrY")
	# color_type=as.character(rep(c('#009900','#FFCC00','#990099','#3333CC'),6))
	color_type=as.character(rep(c("red","blue"),12))
	color_type1=as.character(seq(0.01,1,0.01))
	
	ALL_chr1=c()
	star_ADD=c()
	X_lab=c()
	X_lab_pos=c()
	for (i in 1:5)
		{
		if(length(ALL_chr1)<2)
			{star_add=0}
		else{star_add=ALL_chr1[nrow(ALL_chr1),5]}
		star_ADD=c(star_ADD,star_add)
		sele_chr=subset(fi,chrom==Group[i])
		sele_chr$COL=color_type1[i]
		x_lab=c(seq(30000000,sele_chr[nrow(sele_chr),2],30000000))
		X_lab=c(X_lab,x_lab)
		X_lab_pos=c(X_lab_pos,x_lab+star_add)
		sele_chr$Xpos=sele_chr$maploc+star_add
		ALL_chr1=rbind(ALL_chr1,sele_chr)
		}

	star_ADD=c(star_ADD,ALL_chr1[nrow(ALL_chr1),5]+10000000)
	pg1=ggplot(ALL_chr1,aes(x=Xpos,y=Lograte))+geom_point(aes(color=COL,size=1))+scale_x_continuous(expand = c(0,0),breaks = X_lab_pos,labels = X_lab/1000000)+scale_y_continuous(expand = c(0,0),limits = c(0,5),breaks = c(0,1,2,3,4),labels = paste(seq(0,4),"X",sep=" "))+labs(x="",y="")+geom_hline(yintercept=c(0,1,2,3,4),linewidth=0.5)
	pg1=pg1+theme_classic()+theme(legend.position = "none")+theme(axis.text=element_text(size=38,colour="black"))+scale_color_manual(values =color_type[1:5])+annotate("segment",x=0,xend=0,y=0,yend=5,colour="red",size=3)
	for (i in 1:5)
		{
		pg1=pg1+annotate("segment",x=star_ADD[i],xend=star_ADD[i+1],y=0,yend=0,colour=color_type[i],size=3)+annotate("text",x=star_ADD[i]+15000000,y=4.7,label=Group[i],size=15,colour=color_type[i])
		}
	ALL_chr1=c()
	star_ADD=c()
	X_lab=c()
	X_lab_pos=c()
	for (i in 6:12)
		{
		if(length(ALL_chr1)<2)
			{star_add=0}
		else{star_add=ALL_chr1[nrow(ALL_chr1),5]}
		star_ADD=c(star_ADD,star_add)
		sele_chr=subset(fi,chrom==Group[i])
		sele_chr$COL=color_type1[i]
		
		x_lab=c(seq(30000000,sele_chr[nrow(sele_chr),2],30000000))
		X_lab=c(X_lab,x_lab)
		X_lab_pos=c(X_lab_pos,x_lab+star_add)
		sele_chr$Xpos=sele_chr$maploc+star_add
		ALL_chr1=rbind(ALL_chr1,sele_chr)
		}
	star_ADD=c(star_ADD,ALL_chr1[nrow(ALL_chr1),5])
	pg2=ggplot(ALL_chr1,aes(x=Xpos,y=Lograte))+geom_point(aes(color=COL,size=1))+scale_x_continuous(expand = c(0,0),breaks = X_lab_pos,labels = X_lab/1000000)+scale_y_continuous(expand = c(0,0),limits = c(0,5),breaks = c(0,1,2,3,4),labels = paste(seq(0,4),"X",sep=" "))+labs(x="",y="")+geom_hline(yintercept=c(0,1,2,3,4),linewidth=0.5)
	pg2=pg2+theme_classic()+theme(legend.position = "none")+theme(axis.text=element_text(size=38,colour="black"))+scale_color_manual(values =color_type[6:12])+annotate("segment",x=0,xend=0,y=0,yend=5,colour="red",size=3)
	for (i in 1:7)
		{
		pg2=pg2+annotate("segment",x=star_ADD[i],xend=star_ADD[i+1],y=0,yend=0,colour=color_type[i+5],size=3)+annotate("text",x=star_ADD[i]+15000000,y=4.7,label=Group[i+5],size=15,colour=color_type[i+5])
		}

	ALL_chr1=c()
	star_ADD=c()
	X_lab=c()
	X_lab_pos=c()
	for (i in 13:22)
		{
		if(length(ALL_chr1)<2)
			{star_add=0}
		else{star_add=ALL_chr1[nrow(ALL_chr1),5]}
		star_ADD=c(star_ADD,star_add)
		sele_chr=subset(fi,chrom==Group[i])
		sele_chr$COL=color_type1[i]
		x_lab=c(seq(30000000,sele_chr[nrow(sele_chr),2],30000000))
		X_lab=c(X_lab,x_lab)
		X_lab_pos=c(X_lab_pos,x_lab+star_add)
		sele_chr$Xpos=sele_chr$maploc+star_add
		ALL_chr1=rbind(ALL_chr1,sele_chr)
		}
	star_ADD=c(star_ADD,ALL_chr1[nrow(ALL_chr1),5])
	star_ADD=c(star_ADD,star_ADD[length(star_ADD)]+155270560,star_ADD[length(star_ADD)]+155270560+59373566)
	pg3=ggplot(ALL_chr1,aes(x=Xpos,y=Lograte))+geom_point(aes(color=COL,size=1))+scale_x_continuous(expand = c(0,0),breaks = X_lab_pos,labels = X_lab/1000000)+scale_y_continuous(expand = c(0,0),limits = c(0,5),breaks = c(0,1,2,3,4),labels = paste(seq(0,4),"X",sep=" "))+labs(x="",y="")+geom_hline(yintercept=c(0,1,2,3,4),linewidth=0.5)
	pg3=pg3+theme_classic()+theme(legend.position = "none")+theme(axis.text=element_text(size=38,colour="black"))+scale_color_manual(values =color_type[13:22])+annotate("segment",x=0,xend=0,y=0,yend=5,colour="red",size=3)
	for (i in 1:12)
		{
		pg3=pg3+annotate("segment",x=star_ADD[i],xend=star_ADD[i+1],y=0,yend=0,colour=color_type[i+12],size=3)+annotate("text",x=star_ADD[i]+15000000,y=4.7,label=Group[i+12],size=15,colour=color_type[i+12])
		}

####---
	ALL_chr1=c()
	star_ADD=c()
	X_lab=c()
	X_lab_pos=c()
	for (i in 13:23)
		{
		if(length(ALL_chr1)<2)
			{star_add=0}
		else{star_add=ALL_chr1[nrow(ALL_chr1),5]}
		star_ADD=c(star_ADD,star_add)
		sele_chr=subset(fi,chrom==Group[i])
		sele_chr$COL=color_type1[i]
		x_lab=c(seq(30000000,sele_chr[nrow(sele_chr),2],30000000))
		X_lab=c(X_lab,x_lab)
		X_lab_pos=c(X_lab_pos,x_lab+star_add)
		sele_chr$Xpos=sele_chr$maploc+star_add
		ALL_chr1=rbind(ALL_chr1,sele_chr)
		}
	
	chrY=sele_chr=subset(fi,chrom=="chrY")
	if (mean(chrY$Lograte)<0.3) ### 修改，chrY小于30% 就不画了，认为缺失。
		{
		star_ADD=c(star_ADD,ALL_chr1[nrow(ALL_chr1),5])
		star_ADD=c(star_ADD,star_ADD[length(star_ADD)]+59373566)
		}else{
		chrY$COL=color_type1[24]
		star_ADD=c(star_ADD,ALL_chr1[nrow(ALL_chr1),5])
		chrY$Xpos=chrY$maploc+ALL_chr1[nrow(ALL_chr1),5]
		ALL_chr1=rbind(ALL_chr1,chrY)
		star_ADD=c(star_ADD,star_ADD[length(star_ADD)]+59373566)
			}

	pg4=ggplot(ALL_chr1,aes(x=Xpos,y=Lograte))+geom_point(aes(color=COL,size=1))+scale_x_continuous(expand = c(0,0),breaks = X_lab_pos,labels = X_lab/1000000)+scale_y_continuous(expand = c(0,0),limits = c(0,5),breaks = c(0,1,2,3,4),labels = paste(seq(0,4),"X",sep=" "))+labs(x="",y="")+geom_hline(yintercept=c(0,1,2,3,4),linewidth=0.5)
	pg4=pg4+theme_classic()+theme(legend.position = "none")+theme(axis.text=element_text(size=38,colour="black"))+scale_color_manual(values =color_type[13:24])+annotate("segment",x=0,xend=0,y=0,yend=5,colour="red",size=3)
	for (i in 1:14)
		{
		pg4=pg4+annotate("segment",x=star_ADD[i],xend=star_ADD[i+1],y=0,yend=0,colour=color_type[i+12],size=3)+annotate("text",x=star_ADD[i]+15000000,y=4.7,label=Group[i+12],size=15,colour=color_type[i+12])
		}
	
	# png(paste(name1,"merge_all_chrom.png",sep="_"),width=4000,height=1800)
	# pa=grid.arrange(pg1,pg2,pg4, ncol=1, nrow=3)
	# print(pa)
	# dev.off()
	
	# png(paste(name1,"merge_autosome.png",sep="_"),width=4000,height=1800)
	# pa2=grid.arrange(pg1,pg2,pg3, ncol=1, nrow=3)
	# print(pa2)
	# dev.off()
	
	png(paste(name1,"merge_all_chrom_new_2color.png",sep="_"),width=4000,height=1800)
	pa=grid.arrange(pg1,pg2,pg4, ncol=1, nrow=3)
	print(pa)
	dev.off()
	
	}


###
norm_fold <- function(fil,fold)
	{
	median_count=median(fil)
	fil[which(fil <median_count/fold)]=median_count/fold
	fil[which(fil >median_count*fold)]=median_count*fold
	return(fil)
	}

###
norm_chr_fold = function(fil,up_fold,down_fold) ###correct chr median N fold
	{
	new_da=c()
	chr_in=c(paste("chr",seq(1,22),sep=""),"chrX","chrY")
	for (chr_sele in chr_in)
		{
		da_sele=subset(fil,chr==chr_sele)
		median_count=median(da_sele$count)
		da_sele[which(da_sele$count <median_count*down_fold),6]=median_count*down_fold
		da_sele[which(da_sele$count >median_count*up_fold),6]=median_count*up_fold
		new_da=rbind(new_da,da_sele)
		}
	return(new_da)
	}



##normalize up down percent
norm_percent <- function(new,percent)
	{
	rang=quantile(new,prob = c(percent, 1 - percent))
	new[which(new<rang[1])]=rang[1]
	new[which(new>rang[2])]=rang[2]
	return(new)
	}

####select bin and normalize up down percent
PGS_bin_sele_normalize <-  function(in_file,GC_cut,map_cut)
	{
	in_file=na.omit(in_file)
	colnames(in_file)=c("chr","site","GC","map","count")
	in_file=subset(in_file,GC>GC_cut & map>map_cut)
	return(in_file)
	}

###GC correct
bin_GC_correct <- function(indata,Span)
	{
	message("Correcting for GC bias...")
	indata$GC=round(indata$GC,3)
	M_count=mean(subset(indata,chr!= "chrX" & chr!="chrY")$count)
	indata$GC_loess <- predict(loess((indata$count/M_count) ~ indata$GC,span=Span,degree=2))
	indata$GC_loess <- indata$count/indata$GC_loess
	return(indata)
	}

bin_GC_correct2 <- function(indata,Span)
	{
	library(hash)
	message("Correcting for GC bias...")
	indata$GC=round(indata$GC,3)
	auto_da=subset(indata,chr!= "chrX" & chr!="chrY")
	XY_da=subset(indata,chr== "chrX" |  chr=="chrY")
	M_count=mean(auto_da$count)
	auto_da$GC_loess <- predict(loess((auto_da$count/M_count) ~ auto_da$GC,span=0.35,degree=2))
	h=hash(keys=auto_da$GC,values=auto_da$GC_loess)
	
	# XY_loess=c()
	# for (i in as.character(XY_da$GC) )
		# {
		# #print(i)
		# XY_loess=c(XY_loess,h[[i]])
		# }
	# XY_da$GC_loess=XY_loess
	# indata=rbind(auto_da,XY_da)
	# indata$GC_loess <- indata$count/indata$GC_loess
	
	GC_loess=c()
	for (i in as.character(indata$GC) )
		{
		#print(i)
		if (length(h[[i]]) == 0 || h[[i]]<0 || h[[i]]>10)
			{
			h[[i]]=1
			}
		GC_loess=c(GC_loess,h[[i]])
		}
	indata$GC_loess=GC_loess
	indata$GC_loess <- indata$count/indata$GC_loess
	
	return(indata)
	}


########
GC_bias <- function(da_sele)
	{
	message("Get GC bias which LOESS")
	M_count=mean(da_sele$count)
	GL2=sum((da_sele$count-M_count*(da_sele$count/da_sele$GC_loess))^2)/nrow(da_sele)
	GM2=sum((da_sele$count-mean(da_sele$count))^2)/nrow(da_sele)
	gc_bias=1-GL2/GM2
	return(gc_bias)
	}
########


#### mappability correct
bin_MAP_correct <- function(fil,Span)
	{
	 message("Correcting for mappability bias...")
	fil$map=round(fil$map,3)
	fil$map_loess <- predict(loess((fil$GC_loess/mean(fil$GC_loess)) ~ fil$map,span=Span,degreen=2))
	
	fil$map_loess <- fil$GC_loess/fil$map_loess
	return(fil)
	}

bin_MAP_correct2 <- function(indata,Span)
	{
	 message("Correcting for mappability bias...")
	indata$map=round(indata$map,3)
	auto_da=subset(indata,chr!= "chrX" & chr!="chrY")
	XY_da=subset(indata,chr== "chrX" |  chr=="chrY")
	M_count=mean(auto_da$GC_loess)
	auto_da$map_loess <- predict(loess((auto_da$GC_loess/M_count) ~ auto_da$map,span=Span,degreen=2))
	h=hash(keys=auto_da$map,values=auto_da$map_loess)
	
	# XY_loess=c()
	# for (i in as.character(XY_da$map) )
		# {
		# #print(h[[i]])
		# if (length(h[[i]]) == 0)
			# {
			# h[[i]]=as.numeric(i)
			# }
		# XY_loess=c(XY_loess,h[[i]])
		# }
	# XY_da$map_loess=XY_loess
	# indata=rbind(auto_da,XY_da)
	# indata$map_loess <- indata$GC_loess/indata$map_loess
	
	map_loess=c()
	for (i in as.character(indata$map) )
		{
		#print(h[[i]])
		if (length(h[[i]]) == 0 || h[[i]]<0 || h[[i]]>5)
			{
			h[[i]]=as.numeric(i)
			}
			
		map_loess=c(map_loess,h[[i]])
		}
	indata$map_loess=map_loess
	indata$map_loess <- indata$GC_loess/indata$map_loess
	
	return(indata)
	}

###消除片段差异之后再进行GC和 mappability 矫正，二次矫正###
Reduce_segment_correct_GC_MAP <- function(segment.smoothed.CNA.object2,indata,Span)
	{
	smooth_data=data.frame(segment.smoothed.CNA.object2$data)
	segment=data.frame(segment.smoothed.CNA.object2$output)
	indata=indata[,1:5]
	auto_data=indata[,1:5]
	
	smooth_data_bf=smooth_data
	colnames(smooth_data_bf)=c("chr","Position","count")
	indata=merge(indata,smooth_data_bf,by=c("chr","Position"),all.x = TRUE,sort=FALSE)
	indata$count=2^indata$count
	
	segment_auto=subset(segment,chrom !="chrX" & chrom!="chrY")
	###把片段均值都矫正到 0
	for (si in 1:nrow(segment_auto))
		{
		smooth_data[which(smooth_data$chrom==segment_auto[si,]$chrom & smooth_data$maploc >= segment_auto[si,]$loc.start  & smooth_data$maploc <= segment_auto[si,]$loc.end),]$Lograte=smooth_data[which(smooth_data$chrom==segment_auto[si,]$chrom & smooth_data$maploc >= segment_auto[si,]$loc.start  & smooth_data$maploc <= segment_auto[si,]$loc.end),]$Lograte-segment_auto[si,]$seg.mean
		}
	smooth_data_at=smooth_data
	
	colnames(smooth_data_at)=c("chr","Position","count")
	auto_data=merge(auto_data,smooth_data_at,by=c("chr","Position"),all.x = TRUE,sort=FALSE)
	auto_data$count=2^auto_data$count
	auto_data=subset(auto_data,chr!="chrX"  & chr !="chrY")
	
	M_count=mean(auto_data$count)
	auto_data$GC_loess <- predict(loess((auto_data$count/M_count) ~ auto_data$GC,span=Span,degree=2))
	H_GC=hash(keys=auto_data$GC,values=auto_data$GC_loess)
	#得到矫正哈希，对所有值进行矫正
	
	COL_loess=c()
	for (i in as.character(indata$GC) )
		{##GC的值，性染色体有的常染色体都有，所以不需要进行判断复制
		COL_loess=c(COL_loess,H_GC[[i]])
		}
	indata$GC_loess=COL_loess
	indata$GC_loess=indata$count/indata$GC_loess ###完成GC矫正
	auto_data$GC_loess=auto_data$count/auto_data$GC_loess ##完成GC矫正
	###完成GC矫正，进行mappability 矫正
	
	M_count=mean(auto_data$GC_loess)
	auto_data$map_loess <- predict(loess((auto_data$GC_loess/M_count) ~ auto_data$map,span=Span,degreen=2))
	H_MP=hash(keys=auto_data$map,values=auto_data$map_loess)
	
	MAP_loess=c()
	for (i in as.character(indata$map) )
		{##map chrY有一些特别低的，之前已经做过矫正，这里不进行二次处理，直接赋值为1
		if(length(H_MP[[i]])== 0)
			{
			H_MP[[i]]=1
			}
		MAP_loess=c(MAP_loess,H_MP[[i]])
		}
		indata$map_loess=MAP_loess
		indata$map_loess=indata$GC_loess/indata$map_loess
		##完成矫正，返回矫正后的结果，用于CBS再分析
	return(indata)
	}


Reduce_segment_correct_GC_MAP_raw <- function(segment.smoothed.CNA.object2,indata,Span)
	{
	segment=data.frame(segment.smoothed.CNA.object2$output)
	indata=indata[,1:6]
	auto_data=indata
	auto_data=subset(auto_data,chr!="chrX"  & chr !="chrY")
	M_auto=mean(auto_data$count)
	segment_auto=subset(segment,chrom !="chrX" & chrom!="chrY")
	###把片段均值都矫正到 M_auto
	for (si in 1:nrow(segment_auto))
		{
		# mean_seg=mean(auto_data[which(auto_data$chr==segment_auto[si,]$chrom & auto_data$Position >= segment_auto[si,]$loc.start  & auto_data$Position <= segment_auto[si,]$loc.end),]$count)
		# auto_data[which(auto_data$chr==segment_auto[si,]$chrom & auto_data$Position >= segment_auto[si,]$loc.start  & auto_data$Position <= segment_auto[si,]$loc.end),]$count=
		# auto_data[which(auto_data$chr==segment_auto[si,]$chrom & auto_data$Position >= segment_auto[si,]$loc.start  & auto_data$Position <= segment_auto[si,]$loc.end),]$count-mean_seg+M_auto
		mean_seg=mean(auto_data[which(auto_data$chrom==segment_auto[si,]$chrom & auto_data$Position >= segment_auto[si,]$loc.start  & auto_data$Position <= segment_auto[si,]$loc.end),]$count)
		auto_data[which(auto_data$chrom==segment_auto[si,]$chrom & auto_data$Position >= segment_auto[si,]$loc.start  & auto_data$Position <= segment_auto[si,]$loc.end),]$count=
		auto_data[which(auto_data$chrom==segment_auto[si,]$chrom & auto_data$Position >= segment_auto[si,]$loc.start  & auto_data$Position <= segment_auto[si,]$loc.end),]$count-mean_seg+M_auto
		}
	
	M_count=mean(auto_data$count)
	auto_data$GC_loess <- predict(loess((auto_data$count/M_count) ~ auto_data$GC,span=Span,degree=2))
	H_GC=hash(keys=auto_data$GC,values=auto_data$GC_loess)
	#得到矫正哈希，对所有值进行矫正
	
	COL_loess=c()
	for (i in as.character(indata$GC) )
		{##GC的值，性染色体有的常染色体都有，所以不需要进行判断复制
		if(length(H_GC[[i]])== 0 || H_GC[[i]]<0.02 || H_GC[[i]]>10)
			{
			H_GC[[i]]=1
			}
		COL_loess=c(COL_loess,H_GC[[i]])
		}
	indata$GC_loess=COL_loess
	indata$GC_loess=indata$count/indata$GC_loess ###完成GC矫正
	auto_data$GC_loess=auto_data$count/auto_data$GC_loess ##完成GC矫正
	###完成GC矫正，进行mappability 矫正
	
	M_count=mean(auto_data$GC_loess)
	auto_data$map_loess <- predict(loess((auto_data$GC_loess/M_count) ~ auto_data$map,span=Span,degreen=2))
	H_MP=hash(keys=auto_data$map,values=auto_data$map_loess)
	
	MAP_loess=c()
	for (i in as.character(indata$map) )
		{##map chrY有一些特别低的，之前已经做过矫正，这里不进行二次处理，直接赋值为1
		if(length(H_MP[[i]])== 0 || H_MP[[i]]<0.1 || H_MP[[i]]>5)
			{
			H_MP[[i]]=as.numeric(i)
			}
		MAP_loess=c(MAP_loess,H_MP[[i]])
		}
		indata$map_loess=MAP_loess
		indata$map_loess=indata$GC_loess/indata$map_loess
		##完成矫正，返回矫正后的结果，用于CBS再分析
	return(indata)
	}




#####plot_before after correct
correct_plot <- function(fil,gc_cor_n)
	{
	library(gridExtra)
	p1=ggplot(fil,aes(x=GC,y=count))+geom_point(color="green")+geom_smooth(color="red")+ggtitle("GC Bias in Uncorrected Readcount")+theme(axis.title=element_text(colour = "black",size = 25),axis.text=element_text(color="black",size=20))+scale_x_continuous(breaks = seq(0.32,0.65,0.05))+theme_set(theme_set(theme_bw(base_size=30)))+theme(plot.title = element_text(hjust = 0.5))+annotate("text",x=0.54,y=mean(fil$count),label=paste("GC_bias=",round(GC_b,3),sep=""),size=12,color="red")
	p2=ggplot(fil,aes(x=GC,y=GC_loess))+geom_point(color="green")+geom_smooth(color="red")+ggtitle("GC Bias in GC corrected Readcount")+theme(axis.title=element_text(colour = "black",size = 25),axis.text=element_text(color="black",size=20))+scale_x_continuous(breaks = seq(0.32,0.65,0.05))+theme_set(theme_set(theme_bw(base_size=30)))+theme(plot.title = element_text(hjust = 0.5))
	p3=ggplot(fil,aes(x=map,y=GC_loess))+geom_point(color="green")+geom_smooth(color="red")+ggtitle("mappability Bias in GC corrected Readcount")+theme(axis.title=element_text(colour = "black",size = 25),axis.text=element_text(color="black",size=20))+scale_x_continuous(breaks = seq(0.6,1,0.1))+theme_set(theme_set(theme_bw(base_size=30)))+theme(plot.title = element_text(hjust = 0.5))
	p4=ggplot(fil,aes(x=map,y=map_loess))+geom_point(color="green")+geom_smooth(color="red")+ggtitle("GC and mappability corrected Readcount")+theme(axis.title=element_text(colour = "black",size = 25),axis.text=element_text(color="black",size=20))+scale_x_continuous(breaks = seq(0.6,1,0.1))+theme_set(theme_set(theme_bw(base_size=30)))+theme(plot.title = element_text(hjust = 0.5))
	png(gc_cor_n,width=4200,height=2000)
	pa=grid.arrange(p1,p2,p3,p4, ncol=2, nrow=2)
	print(pa)
	dev.off()
	}

median_cor <- function(col_in,span,diff_mad)
	{
	Median=c()
	Len=length(col_in)
	new_col=c(rev(col_in[(Len-span):Len]),col_in,rev(col_in[1:span]))
	for (i in (span+1):(Len+span))
		{
		median_va=median(new_col[(i-span):(i+span)])
		if ((new_col[i]-median_va)>diff_mad)
			{
			n_median=median_va+diff_mad
			}
		else if ((new_col[i]-median_va) <= -diff_mad)
			{
			n_median=median_va-diff_mad
			}
		else{n_median=new_col[i]}
		Median=c(Median,n_median)
		}
	return(Median)
	}

median_cor2 <- function(col_in,span)
	{
	Median=c()
	Len=length(col_in)
	new_col=c(rev(col_in[(Len-span):Len]),col_in,rev(col_in[1:span]))
	for (i in (span+1):(Len+span))
		{
		median_va=median(new_col[(i-span):(i+span)])
		Median=c(Median,median_va)
		}
	return(Median)
	}


cytoband_out <-function(da1,name)
	{
	FW=c()
	AW=c()
	ratio=c()
	da1=na.omit(da1)
	for (i in 1:nrow(da1))
		{
		FW=c(FW,as.character(subset(cytoBand,V1==da1[i,2] & V2 < da1[i,3] & V3 >= da1[i,3])[,4]))
		AW=c(AW,as.character(subset(cytoBand,V1==da1[i,2] & V2 < da1[i,4] & V3 >= da1[i,4])[,4]))
		ratio=c(ratio,2*(2^da1[i,6]-1))
		}
	da1=data.frame(da1,FW,AW,ratio)
	da1$ratio=round(da1$ratio,2)
	da1$CNV_len.M=(da1$loc.end-da1$loc.start)/1000000+0.6
	da1$ID=name1
	da1=da1[order(-da1$ratio),]
	da1[which(da1$chrom=="chrX" | da1$chrom=="chrY"),]$ratio="NA"
	write.table(da1,file=name,sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)
	}

###

cytoband_selet <-function(da1,lenbin_up,lenbin_dn,minimarkbin,cut_ratio)###先卡物理距离，然后物理距离里面的点数。
	{
	FW=c()
	AW=c()
	ratio=c()
	da1=na.omit(da1)
	for (i in 1:nrow(da1))
		{
		FW=c(FW,as.character(subset(cytoBand,V1==da1[i,2] & V2 < da1[i,3] & V3 >= da1[i,3])[,4]))
		AW=c(AW,as.character(subset(cytoBand,V1==da1[i,2] & V2 < da1[i,4] & V3 >= da1[i,4])[,4]))
		
		ratio=c(ratio,2*(2^da1[i,6]-1))
		}
	da1=data.frame(da1,FW,AW,ratio)
	da1$ratio=round(da1$ratio,3)
	da1$CNV_len.M=(da1$loc.end-da1$loc.start)/1000000+0.6
	da_sele=subset(da1,CNV_len.M >= lenbin_up & CNV_len.M < lenbin_dn  & num.mark >= minimarkbin & abs(ratio) >= cut_ratio & chrom !="chrX" & chrom !="chrY")
	return(da_sele)
	}


cytoband_selet_new <-function(da1,lenbin_up,lenbin_dn,minimarkbin,cut_ratio)###先卡物理距离，然后物理距离里面的点数。
	{
	FW=c()
	AW=c()
	ratio=c()
	da1=na.omit(da1)
	for (i in 1:nrow(da1))
		{
		FW=c(FW,as.character(subset(cytoBand,V1==da1[i,2] & V2 < da1[i,3] & V3 >= da1[i,3])[,4]))
		AW=c(AW,as.character(subset(cytoBand,V1==da1[i,2] & V2 < da1[i,4] & V3 >= da1[i,4])[,4]))
		
		ratio=c(ratio,2*(2^da1[i,6]-1))
		}
	da1=data.frame(da1,FW,AW,ratio)
	da1$ratio=round(da1$ratio,3)
	da1$CNV_len.M=(da1$loc.end-da1$loc.start+1)/1000000
	# da1$CNV_len.M=(da1$loc.end-da1$loc.start+bin_size)/1000000
	da_sele=subset(da1,CNV_len.M >= lenbin_up & CNV_len.M < lenbin_dn  & num.mark >= minimarkbin & abs(ratio) >= cut_ratio & chrom !="chrX" & chrom !="chrY")
	return(da_sele)
	}



####修改为chrXY_filter2  根据chrY的比例是否达到30%来判断是否有Y，如果有Y ，按chrXY为基础判断性染色体组成。
####还是按70% 的比例来判断是否是整条的还是部分或嵌合的


chrXY_filter2 <-function(da)
	{
	#da1=subset(da$output,chrom=="chrX" | chrom=="chrY")
	da1=na.omit(da$output)
	FW=c()
	AW=c()
	ratio=c()
	for (i in 1:nrow(da1))
		{
		FW=c(FW,as.character(subset(cytoBand,V1==da1[i,2] & V2 < da1[i,3] & V3 >= da1[i,3])[,4]))
		AW=c(AW,as.character(subset(cytoBand,V1==da1[i,2] & V2 < da1[i,4] & V3 >= da1[i,4])[,4]))
		}
	
	da1=data.frame(da1,FW,AW)
	#da1$ratio=round(abs(1-2^da1[,6])*2,2)
	da1$ratio=round((2^da1[,6]-1)*2,2)
	mean_X=mean(2^subset(data.frame(da$data),chrom=="chrX")$Lograte)
	mean_Y=mean(2^subset(data.frame(da$data),chrom=="chrY")$Lograte)
	
	###mean_Y 值是否超过30% 作为是否有Y，如果有Y，则按XY为基础来判断，如果没有则按XX来判断
	###修改前则按X,Y的mean值是否达到50% 来判断。
	
	sele_X=subset(da1,chrom=="chrX")
	sele_Y=subset(da1,chrom=="chrY")
	
	if (mean_Y<0.15)
		{
		##以XX为基准 70% 嵌合为阈值
		for (i in 1:nrow(sele_X))
			{
			sele_X[i,9]=(2^sele_X[i,6]-1)*2
			}
		##Y的比例全部为 NA
		for (i in 1:nrow(sele_Y))
			{
			sele_Y[i,9]="NA"
			}
		
		if (mean_X >= 1.35 & mean_X< 2)
			{
			sexchr="XXX"
			}
		else if (mean_X<0.65 & mean_X> 0.15)
			{
			sexchr="X"
			}
		else if (mean_X>=0.65 & mean_X < 1.35)
			{ 
			sexchr="XX"
			}
		else{ sexchr="Undetermined"}
		}
	else{
		#以XY为基准 70% 为阈值
		for (i in 1:nrow(sele_X))
			{
			sele_X[i,9]=(2^sele_X[i,6]-0.5)*2
			}
		
		for (i in 1:nrow(sele_Y))
			{
			sele_Y[i,9]=(2^sele_Y[i,6]-0.5)*2
			}
		
		if ( mean_Y>0.85 & mean_X <0.85 & mean_X>0.15)
			{
			sexchr="XYY"
			}
		else if (mean_X>0.85 & mean_Y <0.85 & mean_Y>0.15)
			{
			sexchr="XXY"
			}
		else if ( mean_X <0.85 & mean_X>0.15 & mean_Y <0.85 & mean_Y>0.15)
			{
			sexchr="XY"
			}
		else{sexchr="Undetermined"}
		
		}
	
	da_sele=subset(da1,chrom!="chrX" & chrom!="chrY")
	da_sele=rbind(da_sele,sele_X,sele_Y)
	sele_XY=rbind(sele_X,sele_Y)
	da_sele$CNV_len.M=(da_sele$loc.end-da_sele$loc.start)/1000000+0.6
	sele_XY=rbind(sele_X,sele_Y)
	sele_XY$CNV_len.M=(sele_XY$loc.end-sele_XY$loc.start)/1000000+0.6
	out_result=list("sexchr"=sexchr,"da_sele"=da_sele,"sele_XY"=sele_XY)
	return(out_result)
	}

sele_sign_sexchr <- function(sele_XY)
	{
	sele_XY=subset(sele_XY,ratio!="NA")
	sign_xy=c()
	sele_XY[,9] = as.numeric(as.character(sele_XY[,9]))
	
	for (i in 1:nrow(sele_XY))
		{
		if (sele_XY[i,5]<20 & sele_XY[i,5] >=14  & abs(sele_XY[i,9])>=0.50)
			{
			sign_xy=rbind(sign_xy,sele_XY[i,])
			}
		else if (sele_XY[i,5]>=20 & sele_XY[i,5]<80 & abs(sele_XY[i,9])>=0.45)
			{
			sign_xy=rbind(sign_xy,sele_XY[i,])
			}
		else if (sele_XY[i,5]>=80 & sele_XY[i,5]<140 & abs(sele_XY[i,9])>=0.34)
			{
			sign_xy=rbind(sign_xy,sele_XY[i,])
			}
		else if (sele_XY[i,5]>140  & sele_XY[i,5]<200  & abs(sele_XY[i,9])>=0.35)
			{
			sign_xy=rbind(sign_xy,sele_XY[i,])
			}
		else if (sele_XY[i,5]>200  & abs(sele_XY[i,9])>=0.30)
			{
			sign_xy=rbind(sign_xy,sele_XY[i,])
			}
		else{
			next
			}
		}
		sign_xy[,9]=round(as.numeric(as.character(sign_xy[,9])),3)
		
	return(sign_xy)
	}
#####通过基因组距离来卡
sele_sign_sexchr2 <- function(sele_XY)
	{
	sele_XY=subset(sele_XY,ratio!="NA")
	sign_xy=c()
	sele_XY[,9] = as.numeric(as.character(sele_XY[,9]))
	for (i in 1:nrow(sele_XY))
		{
		if (sele_XY[i,10]<10 & sele_XY[i,10] >=4  & sele_XY[i,5] >=3 & abs(sele_XY[i,9])>=0.50)
			{
			sign_xy=rbind(sign_xy,sele_XY[i,])
			}
		else if (sele_XY[i,10]<20 & sele_XY[i,10] >=10  & sele_XY[i,5] >=4 & abs(sele_XY[i,9])>= 0.45)
			{
			sign_xy=rbind(sign_xy,sele_XY[i,])
			}
		else if (sele_XY[i,10]<30 & sele_XY[i,10] >=20  & sele_XY[i,5] >=5 & abs(sele_XY[i,9])>= 40)
			{
			sign_xy=rbind(sign_xy,sele_XY[i,])
			}
		else if (sele_XY[i,10]<100 & sele_XY[i,10] >=30  & sele_XY[i,5] >=7 & abs(sele_XY[i,9])>=0.35)
			{
			sign_xy=rbind(sign_xy,sele_XY[i,])
			}
		else if (sele_XY[i,10] >=100  & sele_XY[i,5] >=30 & abs(sele_XY[i,9])>=0.30)
			{
			sign_xy=rbind(sign_xy,sele_XY[i,])
			}
		else{
			next
			}
		}
		sign_xy[,9]=round(as.numeric(as.character(sign_xy[,9])),3)
		
	return(sign_xy)
	}

merge_bin <- function(da,bin_nub)
	{
	chr_m=c()
	Position_m=c()
	GC_m=c()
	map_m=c()
	Bg_m=c()
	count_m=c()
	for (chrom in unique(da$chr))
		{
		chr_sele=subset(da,chr==chrom)
		len=nrow(chr_sele)
		n_less=len%%bin_nub
		
		for (i in seq(n_less+bin_nub,len,bin_nub))
			{
			Position_m=c(Position_m,median(chr_sele[(i-bin_nub+1):i,2]))
			GC_m=c(GC_m,mean(chr_sele[(i-bin_nub+1):i,3]))
			map_m=c(map_m,mean(chr_sele[(i-bin_nub+1):i,4]))
			Bg_m=c(Bg_m,mean(chr_sele[(i-bin_nub+1):i,5]))
			count_m=c(count_m,mean(chr_sele[(i-bin_nub+1):i,6]))
			}
		chr_m=c(chr_m,rep(chrom,length(seq(n_less+bin_nub,len,bin_nub))))
		}
	new_da=data.frame(chr_m,Position_m,round(GC_m,3),round(map_m,3),round(Bg_m,4),count_m)
	colnames(new_da)=c("chr","Position","GC","map","Bg","count")
	return(new_da)
	}


plot_merge_all <- function(fi,name1){
	library(gridExtra)
	fi=data.frame(fi) ###chrom maploc Lograte
	fi$Lograte=2^fi$Lograte*2
	Group=c(paste("chr",seq(1,22),sep=""),"chrX","chrY")
	color_type=as.character(rep(c('#009900','#FFCC00','#990099','#3333CC'),6))
	color_type1=as.character(seq(0.01,1,0.01))
	
	ALL_chr1=c()
	star_ADD=c()
	X_lab=c()
	X_lab_pos=c()
	for (i in 1:23)
		{
		if(length(ALL_chr1)<2)
			{star_add=0}
		else{star_add=ALL_chr1[nrow(ALL_chr1),5]}
		
		star_ADD=c(star_ADD,star_add)
		sele_chr=subset(fi,chrom==Group[i])
		sele_chr$COL=color_type1[i]
		
		x_lab=c(seq(30000000,sele_chr[nrow(sele_chr),2],30000000))
		X_lab=c(X_lab,x_lab)
		X_lab_pos=c(X_lab_pos,x_lab+star_add)
		
		sele_chr$Xpos=sele_chr$maploc+star_add
		ALL_chr1=rbind(ALL_chr1,sele_chr)
		}
	
	chrY=sele_chr=subset(fi,chrom=="chrY")
	if (mean(chrY$Lograte,na.rm=TRUE)<0.25)
		{
		star_ADD=c(star_ADD,ALL_chr1[nrow(ALL_chr1),5])
		star_ADD=c(star_ADD,star_ADD[length(star_ADD)]+59373566)
		}else{
		chrY$COL=color_type1[24]
		star_ADD=c(star_ADD,ALL_chr1[nrow(ALL_chr1),5])
		chrY$Xpos=chrY$maploc+ALL_chr1[nrow(ALL_chr1),5]
		ALL_chr1=rbind(ALL_chr1,chrY)
		star_ADD=c(star_ADD,star_ADD[length(star_ADD)]+59373566)
			}
	
	X_lab_pos_new=c()
	for (i in 1:24){
		X_lab_pos_new_temp = star_ADD[i]+0.5*(star_ADD[i+1]-star_ADD[i])
		X_lab_pos_new=c(X_lab_pos_new,X_lab_pos_new_temp)
	}
	X_lab_new = c("Chr1",2:22,"X","Y")
	pg=ggplot(ALL_chr1,aes(x=Xpos,y=Lograte))+geom_point(aes(color=COL,size=1,alpha=0.8))+scale_x_continuous(expand = c(0,0),limits = c(0,max(star_ADD)),breaks = X_lab_pos_new,labels = X_lab_new)+scale_y_continuous(,limits = c(0,6),breaks = c(0,2,4,6),labels = c(" 0"," 2"," 4"," 6"))+labs(x="",y="")+geom_hline(linetype="dashed",yintercept=c(1,2,3,4,5),linewidth=0.5)
	pg=pg+theme_bw()+theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank())      #theme_bw背景去网格
	pg=pg+ylab("Copy Number")+xlab("\n\n")
	pg=pg+theme(legend.position = "none")+theme(axis.text=element_text(size=50,colour="black"),axis.title=element_text(size=66,colour="black"))+scale_color_manual(values =color_type[1:24])
	pg=pg+theme(axis.text.x = element_text(colour = color_type[1:24]))          #设置轴刻度的颜色
	pg=pg+theme(plot.margin=unit(c(1,1,1,1),'lines'))
	pg=pg+theme(panel.background = element_rect(colour = "black",linewidth = 2))  #将图片区的背景色设置为lightblue，边框颜色为黑色，边框宽度为3

	library(cowplot)
	X_l=c(0.1,0.3,0.5,0.6,0.75,0.9)
	Y_l=c(0.65,0.55,0.8,0.5,0.85,0.45)
	for (i in 1:6)
		{
		# pg = cowplot::ggdraw(pg)+draw_label(as.character('仅供科学研究'),x=X_l[i],y=Y_l[i],fontface = "bold.italic",size=100,angle=0,alpha = .2)
		pg = cowplot::ggdraw(pg)#+draw_label(as.character('仅供科学研究'),x=X_l[i],y=Y_l[i],fontface = "bold.italic",size=100,angle=0,alpha = .2)
		}
	
	# png(paste(name1,"merge_all_chrom_new.png",sep="_"),width=6000,height=800)
	png(paste(name1,"merge_all_chrom.png",sep="_"),width=6000,height=800)
	print(pg)
	dev.off()
	}

get_sd_mean_mad <- function(segment.smoothed.CNA.object2)
	{
	segment=data.frame(segment.smoothed.CNA.object2$output)
	segment_auto=subset(segment,chrom !="chrX" & chrom!="chrY")
	all_data=data.frame(segment.smoothed.CNA.object2$data)
	auto_data=subset(all_data,chrom!="chrX"  & chrom !="chrY")
	colnames(auto_data)=c("chrom","Position","count")
	M_auto=mean(auto_data$count)
	for (si in 1:nrow(segment_auto))
		{
		mean_seg=mean(auto_data[which(auto_data$chrom==segment_auto[si,]$chrom & auto_data$Position >= segment_auto[si,]$loc.start  & auto_data$Position <= segment_auto[si,]$loc.end),]$count)
		auto_data[which(auto_data$chrom==segment_auto[si,]$chrom & auto_data$Position >= segment_auto[si,]$loc.start  & auto_data$Position <= segment_auto[si,]$loc.end),]$count=
		auto_data[which(auto_data$chrom==segment_auto[si,]$chrom & auto_data$Position >= segment_auto[si,]$loc.start  & auto_data$Position <= segment_auto[si,]$loc.end),]$count-mean_seg+M_auto
		}
	auto_data$count2=2^auto_data$count
	nor_mean=mean(auto_data$count2)
	nor_sd=sd(auto_data$count2)
	CV=nor_sd/nor_mean
	sat_out=c(nor_sd,nor_mean,CV)
	return(sat_out)
	}


####

merge_bin_new <- function(da,bin_nub)
	{
	chr_m=c()
	Position_m=c()
	Position_s=c()
	Position_e=c()
	GC_m=c()
	map_m=c()
	Bg_m=c()
	count_m=c()
	for (chrom in unique(da$chr))
		{
		chr_sele=subset(da,chr==chrom)
		len=nrow(chr_sele)
		n_less=len%%bin_nub
		
		for (i in seq(n_less+bin_nub,len,bin_nub))
			{
			Position_m=c(Position_m,median(chr_sele[(i-bin_nub+1):i,2]))
			Position_s=c(Position_s,chr_sele[(i-bin_nub+1),2])
			Position_e=c(Position_e,chr_sele[i,2])
			GC_m=c(GC_m,mean(chr_sele[(i-bin_nub+1):i,3]))
			map_m=c(map_m,mean(chr_sele[(i-bin_nub+1):i,4]))
			Bg_m=c(Bg_m,mean(chr_sele[(i-bin_nub+1):i,5]))
			count_m=c(count_m,mean(chr_sele[(i-bin_nub+1):i,6]))
			}
		chr_m=c(chr_m,rep(chrom,length(seq(n_less+bin_nub,len,bin_nub))))
		}
	# new_da=data.frame(chr_m,Position_m,round(GC_m,3),round(map_m,3),round(Bg_m,4),count_m)
	# colnames(new_da)=c("chr","Position","GC","map","Bg","count")
	new_da=data.frame(chr_m,Position_m,round(GC_m,3),round(map_m,3),round(Bg_m,4),count_m,Position_s,Position_e)
	colnames(new_da)=c("chr","Position","GC","map","Bg","count","Position_start","Position_end")
	return(new_da)
	}


###
change_position <- function(da1,da_merge,cytoBand){
	for (i in 1:nrow(da1))
		{
		# da1[i,]$loc.end=subset(da_merge,chr==da1[i,]$chrom & Position==da1[i,]$loc.end)$Position_end
		da1[i,]$loc.start=subset(da_merge,chr==da1[i,]$chrom & Position==da1[i,]$loc.start)$Position_start
		da1[i,]$loc.end=subset(da_merge,chr==da1[i,]$chrom & Position==da1[i,]$loc.end)$Position_end+bin_size-1
		
		if(da1[i,]$loc.end>max(subset(cytoBand,V1==da1[i,]$chrom)$V3)) da1[i,]$loc.end=max(subset(cytoBand,V1==da1[i,]$chrom)$V3)
		
		}
	return(da1)
}

cytoband_out_new <-function(da1,name)
	{
	FW=c()
	AW=c()
	ratio=c()
	da1=na.omit(da1)
	for (i in 1:nrow(da1))
		{
		FW=c(FW,as.character(subset(cytoBand,V1==da1[i,2] & V2 < da1[i,3] & V3 >= da1[i,3])[,4]))
		AW=c(AW,as.character(subset(cytoBand,V1==da1[i,2] & V2 < da1[i,4] & V3 >= da1[i,4])[,4]))
		ratio=c(ratio,2*(2^da1[i,6]-1))
		}
	da1=data.frame(da1,FW,AW,ratio)
	da1$ratio=round(da1$ratio,2)
	da1$CNV_len.M=(da1$loc.end-da1$loc.start+1)/1000000
	# da1$CNV_len.M=(da1$loc.end-da1$loc.start+bin_size)/1000000
	da1$ID=name1
	da1=da1[order(-da1$ratio),]
	da1[which(da1$chrom=="chrX" | da1$chrom=="chrY"),]$ratio="NA"
	write.table(da1,file=name,sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)
	}

chrXY_filter2_new <-function(da)
	{
	#da1=subset(da$output,chrom=="chrX" | chrom=="chrY")
	da1=na.omit(da$output)
	FW=c()
	AW=c()
	ratio=c()
	for (i in 1:nrow(da1))
		{
		FW=c(FW,as.character(subset(cytoBand,V1==da1[i,2] & V2 < da1[i,3] & V3 >= da1[i,3])[,4]))
		AW=c(AW,as.character(subset(cytoBand,V1==da1[i,2] & V2 < da1[i,4] & V3 >= da1[i,4])[,4]))
		}
	
	da1=data.frame(da1,FW,AW)
	#da1$ratio=round(abs(1-2^da1[,6])*2,2)
	da1$ratio=round((2^da1[,6]-1)*2,2)
	mean_X=mean(2^subset(data.frame(da$data),chrom=="chrX")$Lograte)
	mean_Y=mean(2^subset(data.frame(da$data),chrom=="chrY")$Lograte)
	
	
	###mean_Y 值是否超过30% 作为是否有Y，如果有Y，则按XY为基础来判断，如果没有则按XX来判断
	###修改前则按X,Y的mean值是否达到50% 来判断。
	
	sele_X=subset(da1,chrom=="chrX")
	sele_Y=subset(da1,chrom=="chrY")
	
	###计算mean_Y时过滤前后5%的数据
	da_sm = data.frame(da$data)
	da_sm_Y = subset(da_sm,chrom=="chrY")
	da_sm_Y$nor = 2^da_sm_Y$Lograte
	mean_Y=mean(subset(da_sm_Y,nor>quantile(da_sm_Y$nor,0.05) & nor<quantile(da_sm_Y$nor,0.95))$nor)       #seg.mean没有用所有点计算的，会过滤头尾的部分数据，因此会出现seg.mean与mean_Y性别判断不一致（一者大于0.15，一者小于0.15）的情况，
	##需要与Y染色体seg.mean作比较后再确定（为保持一致性）
	if(nrow(sele_Y)==1) mean_Y = 2^sele_Y$seg.mean
	
	if (mean_Y<0.15)
		{
		##以XX为基准 70% 嵌合为阈值
		for (i in 1:nrow(sele_X))
			{
			sele_X[i,9]=(2^sele_X[i,6]-1)*2
			}
		##Y的比例全部为 NA
		for (i in 1:nrow(sele_Y))
			{
			sele_Y[i,9]="NA"
			}
		
		if (mean_X >= 1.35 & mean_X< 2)
			{
			sexchr="XXX"
			}
		else if (mean_X<0.65 & mean_X> 0.15)
			{
			sexchr="X"
			}
		else if (mean_X>=0.65 & mean_X < 1.35)
			{ 
			sexchr="XX"
			}
		else{ sexchr="Undetermined"}
		}
	else{
		#以XY为基准 70% 为阈值
		for (i in 1:nrow(sele_X))
			{
			sele_X[i,9]=(2^sele_X[i,6]-0.5)*2
			}
		
		for (i in 1:nrow(sele_Y))
			{
			sele_Y[i,9]=(2^sele_Y[i,6]-0.5)*2
			}
		
		if ( mean_Y>0.85 & mean_X <0.85 & mean_X>0.15)
			{
			sexchr="XYY"
			}
		else if (mean_X>0.85 & mean_Y <0.85 & mean_Y>0.15)
			{
			sexchr="XXY"
			}
		else if ( mean_X <0.85 & mean_X>0.15 & mean_Y <0.85 & mean_Y>0.15)
			{
			sexchr="XY"
			}
		else{sexchr="Undetermined"}
		
		}
	
	da_sele=subset(da1,chrom!="chrX" & chrom!="chrY")
	da_sele=rbind(da_sele,sele_X,sele_Y)
	da_sele$CNV_len.M=(da_sele$loc.end-da_sele$loc.start+1)/1000000
	sele_XY=rbind(sele_X,sele_Y)
	sele_XY$CNV_len.M=(sele_XY$loc.end-sele_XY$loc.start+1)/1000000
	out_result=list("sexchr"=sexchr,"da_sele"=da_sele,"sele_XY"=sele_XY)
	return(out_result)
	}


##########-------------------------------------------------------------------------------------------

library(ggplot2)
library(DNAcopy)
library(hash)

args=commandArgs(T)
#args[1]=path  args[2]=input bed
setwd(args[1])
fil=args[2]
bg_all=args[3]
name1 = args[4]
#####band_read
data(cytoBand)
cytoBand$chromNum=sub("chr0","chr",cytoBand$chromNum)
colnames(cytoBand)=c("V1","V2","V3","V4","V5")
cytoBand$V6=rep(c("red","blue"),nrow(cytoBand)/2)
cytoBand$V7=rep(c(0.2,0.4),nrow(cytoBand)/2)
####

# name1=substr(fil,0,nchar(fil)-11)
da=read.table(fil,header=FALSE,sep=",")
# da$V6=norm_fold(da$V6,6)
da=da[,c(1:6)]
colnames(da)=c("chr","Position","GC","map","Bg","count")
ref_info=read.table(bg_all,sep="\t")


# if (nrow(da[which(da$count==0 & da$chr=="chrY"),]) >0){
# da[which(da$count==0 & da$chr=="chrY"),]$count=1
	# }

if (nrow(da[which(da$count==0),]) >0){
da[which(da$count==0),]$count=1
	}

#da=merge_bin(da,3)




da=norm_chr_fold(da,2,0.2)				##PGS 可能会出现部分完全缺失的情况，所以会设低一点
#da=bin_GC_correct(da,0.45)
da=bin_GC_correct2(da,0.45)
da$map_loess=da$GC_loess
#da=bin_MAP_correct(da,0.35)
da=bin_MAP_correct2(da,0.35)
GC_b=GC_bias(da)
da$nor=da$map_loess/median(subset(da,chr !="chrX" & chr != "chrY")$map_loess)
#da$nor=da$nor/da$Bg   ###先不做背景矫正，等二次GC mappability 矫正之后再进行背景矫正
#gc_cor_n=paste(name1,"_GCmapcor.png",sep=""
#correct_plot(da,paste(name1,"_GCmapcor.png",sep=""))

###---局部平滑，初筛片段，都归一化到均值为0，---####
# write.table(da,file=paste(name1,"_correct_da.raw.csv",sep=""),sep=",",append=FALSE,col.names=TRUE,row.names=FALSE,quote=FALSE)
		
da$nor=log2(da$nor)
da1=da[,c(1,2,9)]
colnames(da1)=c("chr","Position","Lograte")
CNA.object2 <- CNA(cbind(da1$Lograte),da1$chr,da1$Position,data.type = "logratio",sampleid = "Lograte")
smoothed.CNA.object2 <- smooth.CNA(CNA.object2,smooth.region=10, outlier.SD.scale=3, smooth.SD.scale=2,trim=0.025)
segment.smoothed.CNA.object2 <- segment(smoothed.CNA.object2,min.width=4,verbose=1)



da_correct2=Reduce_segment_correct_GC_MAP_raw(segment.smoothed.CNA.object2,da,0.35)


da_correct2$nor=da_correct2$map_loess/mean(subset(da_correct2,chr !="chrX" & chr != "chrY")$map_loess)
da=da_correct2

# write.table(da,file=paste(name1,"_correct_da.reduce.csv",sep=""),sep=",",append=FALSE,col.names=TRUE,row.names=FALSE,quote=FALSE)



Min_SD=c()
for (i in 3:ncol(ref_info))
	{
	da$Bg=da$nor/ref_info[,i]
	
	SD=sd(subset(da,chr!="chrX" & chr!="chrY")$Bg)
	Min_SD=c(Min_SD,SD)
	}
# cat(Min_SD)
bg_line=which(Min_SD==min(Min_SD))+2
# print(bg_line)
da=da_correct2
da$nor=da$nor/ref_info[,bg_line]
da_correct2$Bg=ref_info[,bg_line]#把背景值换为CV最小的那个




###----统计各项指标---###
diff_mad=mad(diff(subset(da,chr!="chrX" & chr!='chrY')$nor))

###-----进行变异查找---###
da$nor=log2(da$nor)

##中位数矫正
#da$nor=median_cor2(da$nor,3)

da1=da[,c(1,2,9)]
colnames(da1)=c("chr","Position","Lograte")
CNA.object2 <- CNA(cbind(da1$Lograte),da1$chr,da1$Position,data.type = "logratio",sampleid = "Lograte")
smoothed.CNA.object2 <- smooth.CNA(CNA.object2,smooth.region=5, outlier.SD.scale=3, smooth.SD.scale=2,trim=0.025)

segment.smoothed.CNA.object2 <- segment(smoothed.CNA.object2,min.width=4,verbose=1)
sat_out=get_sd_mean_mad(segment.smoothed.CNA.object2)
nor_sd=sat_out[1]*0.75
nor_mean=sat_out[2]
CV=sat_out[3]



bin_size = min(abs(diff(da$Position)))
# print(paste0("bin_size : ",bin_size))

da_merge=da_correct2
da_merge$count=da_merge$nor/da_merge$Bg
da_merge=da_merge[,1:6]
# ###another plot merge_3 bin
da_merge$count=median_cor2(da_merge$count,3)
# da_merge=merge_bin(da_merge,4)
da_merge=merge_bin_new(da_merge,3)


da_merge$count=log2(da_merge$count)
da_merge1=da_merge[,c(1,2,6)]
colnames(da_merge1)=c("chr","Position","Lograte")
CNA.object2 <- CNA(cbind(da_merge1$Lograte),da_merge1$chr,da_merge1$Position,data.type = "logratio",sampleid = "Lograte")
smoothed.CNA.object2 <- smooth.CNA(CNA.object2,smooth.region=5, outlier.SD.scale=3, smooth.SD.scale=2,trim=0.025)

##输出合并点后binCount结果
da_merge_smooth=data.frame(smoothed.CNA.object2)
# head(da_merge_smooth)
da_merge_smooth$Lograte=2^da_merge_smooth$Lograte*2
colnames(da_merge_smooth)=c("chr","Position","copyNum")
write.table(da_merge_smooth,file=paste(name1,"_merge_bin_correct_da.csv",sep=""),sep=",",append=FALSE,col.names=TRUE,row.names=FALSE,quote=FALSE)


sdundo.CNA.object2 <- segment(smoothed.CNA.object2,min.width=4,verbose=1)
sdundo.CNA.object2$output=change_position(sdundo.CNA.object2$output,da_merge,cytoBand)
name2=paste(name1,"result.txt",sep="_")
# cytoband_out(sdundo.CNA.object2$output,name2)
# cytoband_out_new(sdundo.CNA.object2$output,name2)


sdundo.CNA.object3 <- segment(smoothed.CNA.object2,undo.splits="sdundo",min.width=3,undo.SD=3,alpha = 0.001,verbose=1)
sdundo.CNA.object3$output=change_position(sdundo.CNA.object3$output,da_merge,cytoBand)
name3=paste(name1,"cutsd3_result.txt",sep="_")
# cytoband_out(sdundo.CNA.object3$output,name3)
# cytoband_out_new(sdundo.CNA.object3$output,name3)

# seg_sele30_300M=cytoband_selet(sdundo.CNA.object3$output,30,300,9,0.28)
# seg_sele20_30M=cytoband_selet(sdundo.CNA.object3$output,20,30,7,0.35)
# seg_sele10_20M=cytoband_selet(sdundo.CNA.object3$output,10,20,6,0.40)
# seg_sele3_10M=cytoband_selet(sdundo.CNA.object3$output,3,10,3,0.45)

seg_sele30_300M=cytoband_selet_new(sdundo.CNA.object3$output,30,300,9,0.28)
seg_sele20_30M=cytoband_selet_new(sdundo.CNA.object3$output,20,30,7,0.35)
seg_sele10_20M=cytoband_selet_new(sdundo.CNA.object3$output,10,20,6,0.40)
seg_sele3_10M=cytoband_selet_new(sdundo.CNA.object3$output,3,10,3,0.45)

# out_result3=chrXY_filter2(sdundo.CNA.object3)
out_result3=chrXY_filter2_new(sdundo.CNA.object3)
out_result3$sele_XY=sele_sign_sexchr2(out_result3$sele_XY)
stat=paste(name1,"_summary.csv",sep="")
Out=data.frame(GC_b,diff_mad,nor_sd,nor_mean,out_result3$sexchr)
colnames(Out)=c("GC_bias","MAD(diff)","SD","mean","sexchr")
write.table(Out,file=stat,sep=",",append=TRUE,col.names=TRUE,row.names=FALSE,quote=FALSE)

##cut significant_segment ##把性染色体情况追加进去。
result_sign=rbind(seg_sele3_10M,seg_sele10_20M,seg_sele20_30M,seg_sele30_300M,out_result3$sele_XY)
resu_name=paste(name1,"_significant_segment.txt",sep="")
write.table(result_sign,file=resu_name,sep=",",col.names=TRUE,row.names=FALSE,quote=FALSE)
SEXchr=paste("sexchr:",out_result3$sexchr,sep="")
write.table(SEXchr,file=resu_name,sep=",",append=TRUE,col.names=FALSE,row.names=FALSE,quote=FALSE)


## plotp2(sdundo.CNA.object2,"merge_bin.png")
plot_merge_all(sdundo.CNA.object2$data,name1)###merge_all_chrom     #一行图 -- 修改名字
## plot_merge_chr2(sdundo.CNA.object2$data,name1)     #三行图，两个，分是否包含性染色体 
plot_merge_chr2_new(sdundo.CNA.object2$data,name1)     #三行图，命名与一行图一致，以便生成报告
plot_merge_chr2_new2(sdundo.CNA.object2$data,name1)     #三行图，修改命名 -- 修改颜色为红蓝


