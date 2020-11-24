#setwd("/Users//Documents/R/PPAR pathway")
library("readxl")
library("infotheo")
library("repr")

Merged <- read_excel("data.xlsx")
Merged1 <- t(Merged)
colnames(Merged1) <- as.character(unlist(Merged1[1,]))
Merged1 <- Merged1[-1,]
Merged1 <- t(Merged1)
Merged1 <- type.convert(Merged1) # Merged1 has gene(row) by data(column) 


a1 <- nrow(Merged1) # Number of rows (genes)
b1 <- ncol(Merged1) # Number of columns (samples)
# Create an empty matrix that will have the discretized values
Merged2 <- matrix(, nrow = a1, ncol = b1)
for (i in 1:a1)
{
  mean1 = sum(Merged1[i,])/b1
  for (j in 1:b1)
  {
    if(Merged1[i,j] > mean1){
      Merged2[i,j] <- 1
    } else {
      Merged2[i,j] <- 0
    }
  }
}


Chow1 <- Merged2[,1:40]
HFD1 <- Merged2[,41:80]


comp1 <- matrix(, nrow = a1, ncol = 3)

for (i in 1:a1)
{
  comp1[i,1] = sum(Chow1[i,])
  comp1[i,2] = sum(HFD1[i,])
  comp1[i,3] = comp1[i,2]/comp1[i,1]
}



## Probability of (all) genes upregulated Chow vs HFD
prob1_1 <- matrix(, nrow = nrow(Chow1), ncol = 1)
prob1_2 <- matrix(, nrow = nrow(HFD1), ncol = 1)
prob1 <- matrix(, nrow = nrow(Chow1), ncol = 3)

for (i in 1:nrow(Chow1)){
  prob1_1[i,1] <- (length(which(Chow1[i,] == 1)+1)/(ncol(Chow1)+2))
  prob1_2[i,1] <- (length(which(HFD1[i,] == 1)+1)/(ncol(HFD1)+2))
}

prob1[,1] <- prob1_1[,1]
prob1[,2] <- prob1_2[,1]
prob1[,3] <- prob1_2[,1]/prob1_1[,1]


options(repr.plot.width=8, repr.plot.height=5)
#Plot of absolute probabilities
#plot1_1 <- barplot(t(prob1[,1:2]), legend=c("Chow", "HFD"), col=c("darkblue","red"), beside=TRUE) 
#Plot of Ratios
cols1 <- c("blue", "red")[(prob1[,3] > 1) + 1] 
plot1_2 <- barplot(prob1[,3],  col=cols1, main="Ratio of up-regulated probabilities of genes HFD/Chow (Red = Ratio > 1)"
                   ,xlab="Genes",ylab="Ratio", ylim=c(0, 10), cex.main=0.85) #Plot of ratios
abline(h=1,col=1,lty=20)


## Probability of only output processes upregulated Chow vs HFD
prob2_1 <- matrix(, nrow = 10, ncol = 1)
prob2_2 <- matrix(, nrow = 10, ncol = 1)
prob2 <- matrix(, nrow = 10, ncol = 3)

prob2_1[1,1] <- (length(which(Chow1[23,] == 1))+length(which(Chow1[24,] == 1)))/(2*ncol(Chow1))
prob2_2[1,1] <- (length(which(HFD1[23,] == 1))+length(which(HFD1[24,] == 1)))/(2*ncol(HFD1))

prob2_1[2,1] <- (length(which(Chow1[25,] == 1))+length(which(Chow1[26,] == 1))+length(which(Chow1[27,] == 1))
                 +length(which(Chow1[28,] == 1)))/(4*ncol(Chow1))
prob2_2[2,1] <- (length(which(HFD1[25,] == 1))+length(which(HFD1[26,] == 1))+length(which(HFD1[27,] == 1))
                 +length(which(HFD1[28,] == 1)))/(4*ncol(HFD1))

prob2_1[3,1] <- (length(which(Chow1[29,] == 1)))/(ncol(Chow1))
prob2_2[3,1] <- (length(which(HFD1[29,] == 1)))/(ncol(HFD1))

prob2_1[4,1] <- (length(which(Chow1[30,] == 1))+length(which(Chow1[31,] == 1))+length(which(Chow1[32,] == 1)))/(3*ncol(Chow1))
prob2_2[4,1] <- (length(which(HFD1[30,] == 1))+length(which(HFD1[31,] == 1))+length(which(HFD1[32,] == 1)))/(3*ncol(HFD1))

prob2_1[5,1] <- (length(which(Chow1[33,] == 1)))/(ncol(Chow1))
prob2_2[5,1] <- (length(which(HFD1[33,] == 1)))/(ncol(HFD1))

prob2_1[6,1] <- (length(which(Chow1[34,] == 1)))/(ncol(Chow1))
prob2_2[6,1] <- (length(which(HFD1[34,] == 1)))/(ncol(HFD1))

prob2_1[7,1] <- (length(which(Chow1[35,] == 1))+length(which(Chow1[36,] == 1))+length(which(Chow1[37,] == 1))
                 +length(which(Chow1[38,] == 1))+length(which(Chow1[39,] == 1)))/(5*ncol(Chow1))
prob2_2[7,1] <- (length(which(HFD1[35,] == 1))+length(which(HFD1[36,] == 1))+length(which(HFD1[37,] == 1))
                 +length(which(HFD1[38,] == 1))+length(which(HFD1[39,] == 1)))/(5*ncol(HFD1))

prob2_1[8,1] <- (length(which(Chow1[40,] == 1))+length(which(Chow1[41,] == 1))
                 +length(which(Chow1[42,] == 1)))/(3*ncol(Chow1))
prob2_2[8,1] <- (length(which(HFD1[40,] == 1))+length(which(HFD1[41,] == 1))
                 +length(which(HFD1[42,] == 1)))/(3*ncol(HFD1))

prob2_1[9,1] <- (length(which(Chow1[43,] == 1))+length(which(Chow1[44,] == 1))+length(which(Chow1[45,] == 1))
                 +length(which(Chow1[46,] == 1)))/(4*ncol(Chow1))
prob2_2[9,1] <- (length(which(HFD1[43,] == 1))+length(which(HFD1[44,] == 1))+length(which(HFD1[45,] == 1))
                 +length(which(HFD1[46,] == 1)))/(4*ncol(HFD1))

prob2_1[10,1] <- (length(which(Chow1[47,] == 1))+length(which(Chow1[48,] == 1))
                  +length(which(Chow1[49,] == 1))+length(which(Chow1[50,] == 1))+length(which(Chow1[51,] == 1))
                  +length(which(Chow1[52,] == 1))+length(which(Chow1[53,] == 1))+length(which(Chow1[54,] == 1))
                  +length(which(Chow1[55,] == 1)))/(9*ncol(Chow1))
prob2_2[10,1] <- (length(which(HFD1[47,] == 1))+length(which(HFD1[48,] == 1))
                  +length(which(HFD1[49,] == 1))+length(which(HFD1[50,] == 1))+length(which(HFD1[51,] == 1))
                  +length(which(HFD1[52,] == 1))+length(which(HFD1[53,] == 1))+length(which(HFD1[54,] == 1))
                  +length(which(HFD1[55,] == 1)))/(9*ncol(HFD1))

prob2[,1] <- prob2_1[,1]
prob2[,2] <- prob2_2[,1]
prob2[,3] <- prob2_2[,1]/prob2_1[,1]


options(repr.plot.width=6, repr.plot.height=5)
cols2 <- c("blue", "red")[(prob2[,3] >= 1) + 1]
#plot2_1 <- barplot(t(prob2[,1:2]), col=c("darkblue","red"), beside=TRUE) #Plot of absolute probabilities
plot2_2 <- barplot(prob2[,3], col=cols2, main="Ratio of up-regulated probabilities of processes HFD/Chow (Red = Ratio > 1)"
                   ,xlab="Output Processes",ylab="Ratios", ylim=c(0, 4), cex.main=0.85, cex.names=0.8) #Plot of ratios
abline(h=1,col=1,lty=20)



g1 = matrix(c(0,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),nrow=1)
g2 = matrix(c(0,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),nrow=1)
g3 = matrix(c(0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),nrow=1)
g4 = matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),nrow=1)
g5 = matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1),nrow=1)
g6 = matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),nrow=1)
g7 = matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),nrow=1)
g8 = matrix(c(0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),nrow=1)
g9 = matrix(c(0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),nrow=1)
g10 = matrix(c(0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),nrow=1)
g11 = matrix(c(0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),nrow=1)
g12 = matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),nrow=1)
g13 = matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),nrow=1)
g14 = matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),nrow=1)
g15 = matrix(c(0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),nrow=1)
g16 = matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),nrow=1)
g17 = matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),nrow=1)
g18 = matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),nrow=1)
g19 = matrix(c(0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),nrow=1)
g20 = matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),nrow=1)
g21 = matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),nrow=1)
g22 = matrix(c(0,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),nrow=1)

g = rbind(g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,g17,g18,g19,g20,g21,g22)

ProbChow1 = matrix(, nrow = nrow(g), ncol = ncol(g))
ProbHFD1 = matrix(, nrow = nrow(g), ncol = ncol(g))

for (x in 1:nrow(g)){
  for (y in 1:ncol(g))
    if (g[x,y] == 1){
      m00 <- 0
      m0 <- 0
      n11 <- 0
      n1 <- 0
      for (z in 1:40){
        if(Chow1[x,z] == 0 & Chow1[y,z] == 0){
          m00 <- m00 + 1}
        if(Chow1[x,z] == 0){
          m0 <- m0 + 1}
        if(HFD1[x,z] == 1 & HFD1[y,z] == 1){
          n11 <- n11 + 1}
        if(HFD1[x,z] == 1){
          n1 <- n1 + 1}
      }
      p00 <- (m00+1)/(m0+2) # A->B, P(B=0|A=0)
      q11 <- (n11+1)/(n1+2) # A->B, P(B=1|A=1)
      
      ProbChow1[x,y] <- p00
      ProbHFD1[x,y] <- q11
    }
}

ProbChow1 <- round(ProbChow1,2)
ProbHFD1 <- round(ProbHFD1,2)


set1 = list(1,3,4,6,7,12,13,14,15)
set2 = list(25,26,27,28,40,41,42,29)

cond_prob <- matrix(, nrow = length(set1), ncol = length(set2))

i <- 1
j <- 1

for (iteration1 in set1){
  for (iteration2 in set2){
    x <- iteration1
    y <- iteration2
    n00 <- 0
    n0 <- 0
    for (z in 1:40){
      if(Merged2[x,z] == 0 & Merged2[y,z] == 0){
        n00 <- n00 + 1}
      if(Merged2[x,z] == 0){
        n0 <- n0 + 1}
    }
    r00 <- (n00+1)/(n0+2) # A->B, P(B=0|A=0)
    cond_prob[i,j] <- r00
    j <- j+1
  }
  j <- 1
  i <- i+1
}

cond_prob <- round(cond_prob,2)
cond_prob <- t(cond_prob)

#write.csv(cond_prob, file="cond_prob.csv")