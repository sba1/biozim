# assumes that all rows match the same time point

drawer<-"results"
files<-dir(drawer)
basenames<-unique(gsub("\\.\\d*\\.txt","",files,perl=T))
runs<-as.integer(unique(gsub(".*\\.","",gsub(".txt$","",files,perl=T),perl=T)))

for (name in basenames)
{
	l<-list()
	for (r in runs)
	{
		dat<-read.table(file.path(drawer,paste(name,".",r,".txt",sep="")),h=T)
		l<-append(l,list(dat))
	}

	avg<-l[[1]]

	for (i in 2:length(l))
	{
		avg <- avg + l[[i]]
	}
	avg<-avg / length(l)

	fn<-file.path(drawer,paste(name,".mean.txt",sep=""))
	print(paste("Writing",fn))
	write.table(avg,file=file.path(drawer,paste(name,".mean.txt",sep="")),row.names=T)

	if (length(avg$X)>0)
	{
		fn<-file.path(drawer,paste(name,".mean.pdf",sep=""))
		print(paste("Writing",fn))
		pdf(file=fn)
		plot(avg$Time,avg$X,type="l")
		dev.off()
	}
}
