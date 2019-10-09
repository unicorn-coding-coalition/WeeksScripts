####find and plot distribution of data

library(ggplot2)
library(fitdistrplus)
library(logspline)

data <- scan("D:/Weeks/Data/NetworkAnalysis/resCorrDistFiles/distListAll.txt", 
             what = numeric(), sep = ",",
             quote = NULL)

sum(data)
hist(data)
histInfo <- hist(data)
descdist(data, discrete = FALSE)
# distAll from ringmap follows a beta distribution with alpha = 2 and beta = 5
#dbeta(data,2,5)
#fitBeta <- fitdist(data, "beta")

sum = 0
probList = list()
probList <- histInfo$counts
probList <- lapply(probList, data= data, function(x,data) {
        x/length(data)
})

for (count in histInfo$counts){
        sum = sum + count/length(data)
        print(count/length(data))
}


if (sum != 1){
        print("sum of probabilities does not equal 1")
}



