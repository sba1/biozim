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

	print("Calculate")
	#print(cmd.and.args)
	result.raw<-system(cmd.and.args,intern=T)
	print("Calculation completed!")
	
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



mean = vector(mode="numeric",length=10)
fanf = vector(mode="numeric",length=10)


for (b in 1:10) {
    kP<-b * gR  ## b is burst
    # simulate and sort data according to time
    result<-sbml.sim("data/stochastic/noise/singlegene.xml",100000,stochastic=T,runs=5000,sample.steps=100,
							params=c(paste("degradation_G",gP,sep="="),
									 paste("degradation_moG",gR,sep="="),
					                 paste("transcription_oG0",kR,sep="="),
									 paste("translation_oG",kP,sep="=")))
result<-result[order(result$Time),]

qqq<-split(result,result$Time)  
l = length(qqq)
m = mean(qqq[[l]]$G)
s = sd (qqq[[l]]$G)
f = (s*s)/m 

mean[b] = m
fanf[b] = f

}

plot(fanf,mean,xlab="<p>",ylab="fano(p)")

#filename<-paste("singlegene-b",burst,".pdf",sep="")
#Cairo_pdf(filename,width=10,height=10)
#par(mfrow=c(2, 2))
#plot(m$Time,m$G,xlab="Time",type="l",main=sprintf("Protein counts for G with b=%d",burst))
#plot(m$Time,s$G,xlab="Time",type="l",main=sprintf("Standard deviation G with b=%d",burst))
#plot(m$Time,co$G,xlab="Time",type="l",main=sprintf("Cooefficient of variation for G with b=%d",burst))
#plot(m$Time,fano$G,xlab="Time",type="l",main=sprintf("Fano factor for G with b=%d",burst))
#dev.off()
