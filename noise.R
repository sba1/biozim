library(cairoDevice)

############################################################
# Create a data frame from a tab-separated character vector 
create.frame <- function(result.raw)
{
	result.splitted<-strsplit(result.raw,"\t")
	result.frame<-data.frame(matrix(as.numeric(unlist(result.splitted[-1])),ncol=length(result.splitted[[1]]),byrow=T))
	colnames(result.frame)<-make.names(result.splitted[[1]])
	
	return(result.frame)
}

#############################################################
# Simulate the given SBML file
sbml.sim<-function(filename,maxtime=1,stochastic=T,runs=1,sample.steps=5000,params=NULL)
{
	cmd<-"./sbml2ode";
	cmd.and.args<-paste(cmd,"--stochastic","--maxtime", maxtime, "--runs",runs,"--sample-steps",sample.steps,filename);

	if (!is.null(params))
	{
		cmd.and.args<-paste(cmd.and.args," --set-var"," \"",params,"\"",collapse=" ",sep="")
	}

	result.raw<-system(cmd.and.args,intern=T)
	
	result<-create.frame(result.raw)
	
	return(result)
}

#######################
# single gene
burst<-2
gR<-log(2)/120
kR<-0.01
gP<-log(2)/3600
kP<-burst * gR

# simulate and sort data according to time
result<-sbml.sim("data/stochastic/noise/singlegene.xml",40000,stochastic=T,runs=1000,sample.steps=250,
							params=c(paste("degradation_G",gP,sep="="),
									 paste("degradation_moG",gR,sep="="),
					                 paste("transcription_oG0",kR,sep="="),
									 paste("translation_oG",kP,sep="=")))
result<-result[order(result$Time),]

# build groups according to time
rle<-rle(result$Time)
cs<-cumsum(rle$lengths)
min.max.frame<-data.frame(min=c(1,cs[-length(cs)]+1),max=cs)

# create statistics for every group
m<-data.frame(t(apply(min.max.frame,1,function(x){ mean(result[seq(x[1],x[2]),]) })))
s<-data.frame(t(apply(min.max.frame,1,function(x){ sd(result[seq(x[1],x[2]),]) })))
co<-s/m # coefficient of variation
fano<-(s*s)/m

filename<-paste("singlegene-b",burst,".pdf",sep="")
Cairo_pdf(filename,width=10,height=10)
par(mfrow=c(2, 2))
plot(m$Time,m$G,xlab="Time",type="l",main=sprintf("Protein counts for G with b=%d",burst))
plot(m$Time,s$G,xlab="Time",type="l",main=sprintf("Standard deviation G with b=%d",burst))
plot(m$Time,co$G,xlab="Time",type="l",main=sprintf("Cooefficient of variation for G with b=%d",burst))
plot(m$Time,fano$G,xlab="Time",type="l",main=sprintf("Fano factor for G with b=%d",burst))
dev.off()

#########################
# ffl
burst<-2
result<-sbml.sim("data/stochastic/noise/ffl.xml",40000,stochastic=T,runs=50,sample.steps=250,
						params=c(paste("burst",burst,sep="=")))
result<-result[order(result$Time),]
species<-"C"

# build groups according to time
rle<-rle(result$Time)
cs<-cumsum(rle$lengths)
min.max.frame<-data.frame(min=c(1,cs[-length(cs)]+1),max=cs)

# create statistics for every group
m<-data.frame(t(apply(min.max.frame,1,function(x){ mean(result[seq(x[1],x[2]),]) })))
s<-data.frame(t(apply(min.max.frame,1,function(x){ sd(result[seq(x[1],x[2]),]) })))
co<-s/m # coefficient of variation
fano<-(s*s)/m

filename<-paste("ffl-b",burst,".pdf",sep="")
#Cairo_pdf(filename,width=10,height=10)
par(mfrow=c(2, 2))
plot(m$Time,m[,species],xlab="Time",type="l",main=sprintf("Protein counts for %s with b=%d",species,burst))
plot(m$Time,s[,species],xlab="Time",type="l",main=sprintf("Standard deviation %s with b=%d",species,burst))
plot(m$Time,co[,species],xlab="Time",type="l",main=sprintf("Cooefficient of variation for %s with b=%d",species,burst))
plot(m$Time,fano[,species],xlab="Time",type="l",main=sprintf("Fano factor for %s with b=%d",species,burst))
#dev.off()


##########################
## ffl burst 2
#burst<-2
#result<-sbml.sim("data/stochastic/noise/ffl-2.xml",40000,stochastic=T,runs=200,sample.steps=250)
#result<-result[order(result$Time),]
#species<-"C"

## build groups according to time
#rle<-rle(result$Time)
#cs<-cumsum(rle$lengths)
#min.max.frame<-data.frame(min=c(1,cs[-length(cs)]+1),max=cs)

## create statistics for every group
#m<-data.frame(t(apply(min.max.frame,1,function(x){ mean(result[seq(x[1],x[2]),]) })))
#s<-data.frame(t(apply(min.max.frame,1,function(x){ sd(result[seq(x[1],x[2]),]) })))
#co<-s/m # coefficient of variation
#fano<-(s*s)/m


#filename<-paste("ffl-b",burst,".pdf",sep="")
#Cairo_pdf(filename,width=10,height=10)
#par(mfrow=c(2, 2))
#plot(m$Time,m[,species],xlab="Time",type="l",main=sprintf("Protein counts for %s with b=%d",species,burst))
#plot(m$Time,s[,species],xlab="Time",type="l",main=sprintf("Standard deviation %s with b=%d",species,burst))
#plot(m$Time,co[,species],xlab="Time",type="l",main=sprintf("Cooefficient of variation for %s with b=%d",species,burst))
#plot(m$Time,fano[,species],xlab="Time",type="l",main=sprintf("Fano factor for %s with b=%d",species,burst))
#dev.off()

#apply(min.max.frame[c(1,2),],1,function(x){ mean(result[seq(x[1],x[2]),]) })
# clumsy loop
#for (a in l)
#{
#	print(a)
#}

