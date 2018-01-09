#install.packages(....)
#install.packages("itsmr")
library(itsmr)
library(dplyr) #manipulation functions for dataframes
library(ggplot2) #allows for cool plots
library(zoo) #for rolling median (runmed)
library(forecast) #for decomp function 
library(TSA) 
library(stats) 
#install.packages("devtools") #################install this once then comment out
#devtools::install_github("twitter/AnomalyDetection") ##################### install this once then comment out
library(AnomalyDetection) #twitters stuff

init_datatable <- read.table("timedata.txt",sep=",") #change the directory into the file with time as miliseconds
init_datatable <- init_datatable[2:nrow(init_datatable),] ###########RUN THIS TO REMOVE FIRST LINE OR ELSE THERE WILL BE NA VALUES####

init_datatable
########################CLEAN AND SCALE DATA################################
mutate_data <- function(x) { #takes in the data w/time series with time as miliseconds, returns scaled data
datatablex <- x
datatablex[,2] <- as.numeric(as.character(datatablex[,2]))
datatablex[,3] <- as.numeric(as.character(datatablex[,3]))
colnames(datatablex) <- c("no","time","memory")
#################time starts at 0#############################
datatablex[,2] <- datatablex[,2] - min(datatablex[,2])
#feature scaling memory
datatablex$memory <- (datatablex$memory - min(datatablex$memory))/(max(datatablex$memory) - min(datatablex$memory))
return(datatablex)
}

datatable <- mutate_data(init_datatable) #actually a dataframe

#ssp <- spectrum(datatable[,3])
#per <- 1/ssp$freq[ssp$spec==max(ssp$spec)]
estimate_periodicity <- function(x){ #estimate the periodicty, input is the dataframe, returns periodicity
y <- x
ttime <- ts(y)
p <- periodogram(ttime)
dd <- data.frame(freq=p$freq, spec=p$spec)
order <- dd[order(-dd$spec),]
top2 <- head(order, 2)
# display the 2 highest "power" frequencies
per <- 1/top2$f[1]
return(per)
}

per <- estimate_periodicity(datatable[,3])

##################TIME SERIES = T + S + W##########################################
######################MEDIAN DECOMP################################################
detrend_periodicity_func <- function(dataval,period) { #inputs the dataframe and the periodicity, returns residuals
mediantrend <- runmed(dataval,period) #uses rolling median
#plot(mediantrend)
detrended <- dataval - mediantrend
#plot(detrended)
newx <- seq(1:length(detrended))
newdetrend <- data.frame(newx,detrended)
fit.lm2 <- lm(detrended ~ sin(2*pi/per*newx)+cos(2*pi/per*newx)+sin(4*pi/per*newx)+cos(4*pi/per*newx)) #fit harmonic regression curves
summary(fit.lm2) #summary
####################convert to easy equation################ #convert to beta + r*sin(x + alpha) simple sin wave reduction
b0 <- coef(fit.lm2)[1]
alpha <- coef(fit.lm2)[2]
beta <- coef(fit.lm2)[3]
r <- sqrt(alpha^2 + beta^2)
phi <- atan2(beta, alpha)
###############subtract periodicity function###
randomsub <- detrended -(b0 + r * sin(newx + phi)) #subtract
return(randomsub)
}

randomsub <- detrend_periodicity_func(datatable[,3],per) #this will change period to 811, but that doesnt affect the accuracy of the model

hist(randomsub)##########approximately normal
plot(randomsub)
qqnorm(randomsub)
qqline(randomsub)

###############################ANOMALY DETECTION FUNCTION####################
anomaly_detect <- function(randomsub) {
####################anomaly detection###########################
min <-  mean(randomsub, na.rm = T) - 4*sd(randomsub, na.rm = T) #4 stds below mean
max <-  mean(randomsub, na.rm = T) + 4*sd(randomsub, na.rm = T) #4stds above mean
plot(as.numeric(randomsub))
abline(h=max, col="#e15f3f", lwd=2)
abline(h=min, col="#e15f3f", lwd=2)

position <- data.frame(id=seq(1, length(randomsub)), value=randomsub) #assigns id to residuals

anomalyH <- position[position[,2] >= max, ] #keeps anomalous rows
anomalyH <- anomalyH[!is.na(anomalyH[,2]), ] #removes any NAS
anomalyL <- position[position[,2] <= min, ]
anomalyL <- anomalyL[!is.na(anomalyL[,2]), ]
anomaly <- data.frame(id=c(anomalyH[,1], anomalyL[,1]), #dataframe of ids and anomaly values
                     value=c(anomalyH[,2], anomalyL[,2]))
anomaly <- anomaly[!is.na(anomaly$value), ] #removes any NAS if there are any


realx <- data.frame(id=seq(1, length(datatable[,3])), value=datatable[,3])
realAnomalyx <- realx[anomaly[,1],]
medianplot <- plot((datatable[,3]),xlab="data points", ylab="memory")
###for time series plot: uncomment next line####
####medianplot <- plot(as.ts(datatable[,3]),xlab="data points", ylab="memory")
medianplot + points(x = realAnomalyx[,1], y =realAnomalyx[,2], type = "p", col="red", pch=16)
time <- datatable[realAnomalyx[,1],2]
memory <- realAnomalyx[,2]
realAnomalyx[,2] <- time
realAnomalyx <- cbind(realAnomalyx,memory)
realAnomalyx <- cbind(realAnomalyx,init_datatable[realAnomalyx[,1],3])
colnames(realAnomalyx) <- c("id","time","memoryscaled","memory")
return(realAnomalyx)
}

final_anom <- anomaly_detect(randomsub) #plots
final_anom #check in external window for plot, this also is a dataframe of all the anomaly values




#####################DIFFERENCING METHOD#########################################
seasonal_diff <- diff(datatable[,3], differences = 810)
plot(seasonal_diff)
length(seasonal_diff)
deseasoned <- datatable[811:nrow(datatable),3]
deseasoned <- deseasoned - seasonal_diff
plot(deseasoned)
deseasoned_detrended <- diff(deseasoned,differences=1)
plot(deseasoned_detrended)
feature_scaled_detrend <- (deseasoned_detrended - min(deseasoned_detrended))/max(deseasoned_detrended) - min(deseasoned_detrended)
plot(feature_scaled_detrend)
##############EXTRA STUFF HERE#########################################
#######################TWITTER#########################################
res = AnomalyDetectionVec(as.numeric(datatable[,3]), max_anoms=.49, direction='both', period=810, only_last=FALSE, plot=TRUE) 
res$plot


#############################using decomp###################################
decomposed <- decompose(ts(datatable$memory, frequency = 810))
plot(decomposed)
random <- decomposed$random
min = mean(random, na.rm = T) - 4*sd(random, na.rm = T)
max = mean(random, na.rm = T) + 4*sd(random, na.rm = T)
plot(as.ts(as.vector(random)), ylim = c(-0.5,2.5))
abline(h=max, col="red", lwd=2)
abline(h=min, col="red", lwd=2)

position = data.frame(id=seq(1, length(random)), value=random)
anomalyH = position[position$value > max, ]
anomalyH = anomalyH[!is.na(anomalyH$value), ]
anomalyL = position[position$value < min, ]
anomalyL = anomalyL[!is.na(anomalyL$value)]
anomaly = data.frame(id=c(anomalyH$id, anomalyL$id),
                     value=c(anomalyH$value, anomalyL$value))
points(x = anomaly$id, y =anomaly$value, col="#e15f3f")

plot(datatable[,3])
real = data.frame(id=seq(1, length(datatable[,3])), value=datatable[,3])
realAnomaly = real[anomaly[,1], ]
points(x = realAnomaly[,1], y =realAnomaly[,2], col="#e15f3f")
random_na <- random[!is.na(random)]
density(random_na)
hist(random_na)
qqnorm(random_na)
qqline(random_na)
summary(random_na)
#plot(decomposed)
#hist(random)



##################ALT APPROACH###################
det <- as.vector(detrended)
kx <- c(1:length(det))
idnumbers <-  split(kx,ceiling(seq_along(kx)/811))
#det <- data.frame(kx,det)
detsplit <- split(det,ceiling(seq_along(det)/811))
length(detsplit)

splitlist <- list()
for(i in 1:length(detsplit)) { #add id numbers
  splitlist[[length(splitlist)+1]] <- data.frame(idnumbers[[i]],detsplit[[i]])
}
xx <- plot(splitlist[[1]]) 
xx
