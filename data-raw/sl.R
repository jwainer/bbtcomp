## code to prepare `DATASET` dataset goes here

ll <- read.csv("ll.csv", header=TRUE)


dbs = c("biomed","breast","breast_w","buggyCrx", "clean1", "cmc" , "colic", "corral" ,  "credit_g",
        "diabetes"  ,"ionosphere", "irish" , "molecular_b...y_promoters",  "monk3", "prnn_crabs",
        "prnn_synth" , "saheart", "threeOf9", "tokyo1", "vote" )

algs = c("dt","lda","lgbm","xgb", "svm")

sl <-  ll[, c("db",algs)]

ss <- sl[sl$db %in% dbs, ]

usethis::use_data(ll, sl, ss, overwrite = TRUE)

