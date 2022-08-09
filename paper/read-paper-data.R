ss = read.csv("ss.csv")

ssmean = ss %>% group_by(db) %>% summarise(across(.fns = mean)) 
ssalgnames = colnames(ss)[-1]
am = ssmean[,-1]


sl = read.csv("sl.csv")
sl = sl[complete.cases(sl),]

slmean = sl %>% group_by(db) %>% summarise(across(.fns = mean)) 



ll = read.csv("ll.csv")
ll = ll[complete.cases(ll),]
llmean = ll %>% group_by(db) %>% summarise(across(.fns = mean)) 


ss10 = read.csv("ss10.csv")
ll10 = read.csv("ll10.csv")
