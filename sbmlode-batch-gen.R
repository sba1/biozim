# assumes that all rows match the same time point

drawer<-"results"
files<-dir(drawer)
pattern<-"\\.\\d*\\.txt" # we are only interested in these patterns
files.filtered<-grep(pattern,files,perl=T,value=T)
basenames<-unique(gsub(pattern,"",files.filtered,perl=T))
runs<-as.integer(unique(gsub(".*\\.","",gsub(".txt$","",files.filtered,perl=T),perl=T)))

for (name in basenames)
{
	l<-list()
	for (r in runs)
	{
		dat<-read.table(file.path(drawer,paste(name,".",r,".txt",sep="")),h=T)
		l<-append(l,list(dat))
	}

	avg<-l[[1]]
	taken<-1

	for (i in 2:length(l))
	{
		if (nrow(avg) == nrow(l[[i]]))
		{
			avg <- avg + l[[i]]
			taken<-taken+1
		}
	}
	avg <- avg / taken

	fn<-file.path(drawer,paste(name,".mean.txt",sep=""))
	print(paste("Writing ",fn," (average of ",taken,")",sep=""))
	write.table(avg,file=file.path(drawer,paste(name,".mean.txt",sep="")),row.names=T)

	if (length(avg$X)>0)
	{
		fn<-file.path(drawer,paste(name,".mean.pdf",sep=""))
		print(paste("Writing",fn))
		pdf(file=fn)
		plot(avg$Time,avg$X,type="l")
		dev.off()
	} else
	{
		if (length(avg$P)>0 && length(avg$P2)>0)
		{
			fn<-file.path(drawer,paste(name,".mean.pdf",sep=""))
			print(paste("Writing",fn))
			pdf(file=fn)
			plot(avg$Time,avg$P,type="l")
			lines(avg$Time,avg$P2,type="l")
			dev.off()
		}
	}
}
