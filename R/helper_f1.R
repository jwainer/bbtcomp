ss = read.csv(system.file("testdata","ss.csv", package = "bbtcomp"), header=TRUE)
ssmean = read.csv(system.file("testdata","ssmean.csv", package = "bbtcomp"),header = TRUE)
ssmeanx = ssmean[-c(3,5,7,8,12, 20),]

amean = as.matrix(as.data.frame(ssmean)[,-1])

sssd = read.csv(system.file("testdata","sssd.csv", package = "bbtcomp"), header = TRUE)
