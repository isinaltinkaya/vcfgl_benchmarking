set.seed(42)


args = commandArgs(trailingOnly = TRUE)
if(length(args) != 2){
	stop("Usage: Rscript plot_beta.R <beta_mean> <beta_var>")
}
x_mean = as.numeric(args[1])
# x_var = as.numeric(args[2])
# read x_var from scientific notation (e.g. 1e-5)
x_var = as.numeric(format(as.numeric(
	as.character(args[2]))))

# stop(as.numeric(args[2]))
# stop(format(as.character(args[2])))


betaSampler<-function(mean,var,n){
	alpha = (((1.0 - mean) / var) - (1/mean)) * (mean^ 2)
	beta = alpha * ((1/mean)-1)
	return(list(alpha=alpha,beta=beta,data=rbeta(n,shape1=alpha,shape2=beta)))
}

x<-betaSampler(mean=x_mean,var=x_var,n=1000)


# plot(density(x$data),main=paste0("Sampled values from rbeta in R (n=",length(x$data),", mean=",round(mean(x$data),8),", variance=",format(round(var(x$data),8),scientific=T),")"))


# # dbeta
x2 <- seq(0, 1, length=100000)
y<-dbeta(x2,shape1=x$alpha,shape2=x$beta)

k=sqrt(x_var)*1e4
xlims=c(x_mean-(x_mean*k), x_mean+(x_mean*k))
xlims[1]=0

fn=paste0("betaDensity","-mean_",x_mean,
			 "-var_",x_var,
			 ".png"
			 )

cat("Writing plot to file ",fn,"\n")
png(fn)
plot(x2,y,xlab="",ylab="",xlim = xlims, type="l", main=paste0("Density of Beta Distribution\nmean=",x_mean,", variance=",format(x_var,scientific=T),"\nalpha=",x$alpha,", beta=",x$beta))

abline(v=x_mean,col="red")

legend("top",legend=c("mean"), fill=c( "red"))

dev.off()



