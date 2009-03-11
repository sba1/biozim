drawer<-"results"
orig.drawer<-"data/stochastic/dsmts"

files<-dir(drawer)
pattern<-"-stat\\.txt" # we are only interested in names conforming with this  patterns
files.filtered<-grep(pattern,files,perl=T,value=T)
basenames<-sort(unique(gsub(pattern,"",files.filtered,perl=T)))

cols<-2
rows<-2

pdf(width=cols*5,height=rows*5,file="all.pdf")
par(mfrow=c(rows,cols))

for (name in basenames)
{
	dat<-read.table(file.path(drawer,paste(name,"-stat.txt",sep="")),h=T)
	dat<-split(dat,dat$Run)

	avg<-dat[[1]]
	sd<-dat[[2]]
	avg.orig<-read.csv(file.path(orig.drawer,paste(name,"-mean.csv",sep="")),h=T)
	sd.orig<-read.csv(file.path(orig.drawer,paste(name,"-sd.csv",sep="")),h=T)
	
	if (length(avg$X)>0)
	{
		plot(avg$Time,avg$X,type="l",main=paste(name,"by biozim"))
		plot(avg.orig$time,avg.orig$X,type="l",main=name)

		plot(sd$Time,sd$X,type="l",main=paste(name,"by biozim"))
		plot(sd.orig$time,sd.orig$X,type="l",main=name)
	} else
	{
		if (length(avg$P)>0 && length(avg$P2)>0)
		{
			plot(avg$Time,avg$P,type="l",main=paste(name,"by biozim"))
			lines(avg$Time,avg$P2,type="l")
			plot(avg.orig$Time,avg.orig$P,type="l",main=name)
			lines(avg.orig$Time,avg.orig$P2,type="l")

			plot(sd$Time,sd$P,type="l",main=paste(name,"by biozim"),ylim=c(0,max(sd$P)))
			lines(sd$Time,sd$P2,type="l")
			plot(sd.orig$Time,sd.orig$P,type="l",main=name)
			lines(sd.orig$Time,sd.orig$P2,type="l")
		}
	}
}
dev.off()

