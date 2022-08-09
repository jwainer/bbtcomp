
library(conflicted)
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
library(scmamp)
library(dplyr)
#library(R.cache)
devtools::load_all("../bbtcomp/")
library(bbtcomp)

options(mc.cores = parallel::detectCores(logical = FALSE))

source("read-paper-data.R")

colMedians <- function(m) apply(m,2,median)


pairbayes <- function(xx){
  ord = names(sort(colMedians(xx), decreasing = T))
  yy = xx[, ord]
  n = length(ord)
  mat = NULL
  for (i in 1 : (n - 1))
    for (j in (i + 1) : n) {
      nam = paste(ord[i], ">", ord[j])
      r0=bSignedRankTest(yy[[i]], y=yy[[j]], rope=c(-0.0, 0.0))
      above0 = r0$posterior.probabilities[3]
      res=bSignedRankTest(yy[[i]], y=yy[[j]], rope=c(-0.01, 0.01))
      inrope = res$posterior.probabilities[2]
      aboverope = res$posterior.probabilities[3]
      mat = rbind(mat, data.frame(name = nam,
                                  above0 = above0,
                                  in.rope = inrope,
                                  above.rope = aboverope))
      
    }
  rownames(mat)=NULL
  return(list(tab = mat, order=ord))
}  

pairwilcox <- function(xx,adj="hochberg"){
  n1 = sort(colnames(xx))
  zz = sort(apply(xx, 2, median), decreasing = T)
  ord = names(zz)
  xx = xx[, ord]
  ndb = dim(xx)[1]
  x = as.vector(as.matrix(xx))
  g = rep(colnames(xx), each = ndb)
  defaultW <- getOption("warn") 
  options(warn = -1) 
      res = pairwise.wilcox.test(x, g, paired = T, p.adjust.method = adj)
  options(warn = defaultW)
  pv = res$p.value
  nn = length(ord)
  m = matrix(NA, nn, nn)
  m[lower.tri(m)] = pv[lower.tri(pv,diag = T)]
  m[upper.tri(m)] = t(pv)[upper.tri(pv,diag = T)]
  colnames(m) = n1
  rownames(m) = n1
  n = length(ord)
  out = data.frame(name = character(0), p.value = numeric(0))
  for (i in 1 : (n - 1)){
    xi = ord[i]
    for (j in (i + 1) : n){
      xj = ord[j]
      out = rbind(out, list(name = paste(xi, ">", xj),
                            p.value = round(m[xi, xj], 2)))
    }
  }
  order = data.frame(algorithm = ord, "median"= as.numeric(zz))
  rownames(order)=NULL
  return(list(tab = out,order = order))
}


procNemenyi <- function(dat){
  mm = rankMatrix(dat)
  ord = names(sort(colMeans(mm)))
  nn = nemenyiTest(dat)
  mat = ifelse(abs(nn$diff.matrix) > nn$statistic, "yes", " ")
  colnames(mat) = colnames(nn$diff.matrix)
  rownames(mat) = colnames(mat)
  n = length(ord)
  out = data.frame(name = character(0), "signif?" = character(0))
  for (i in 1 : (n - 1)){
    xi = ord[i]
    for (j in (i + 1) : n){
      xj = ord[j]
      out = rbind(out, list(name = paste(xi, ">", xj),
                           dif = mat[xi, xj]))
    }
  }
  order = data.frame(algorithm = ord, "mean rank"= as.numeric(sort(colMeans(mm))))
  rownames(order)=NULL
  return(list(tab = out, order = order))
}




#bbtcomp=addMemoization("bbtcomp")

#m_ssmean = bbtcomp(ssmean) 
#m_slmean = bbtcomp(slmean) 
#m_llmean = bbtcomp(llmean) 

splitdf <- function(df,split=3){
  n = nrow(df)
  nc = ncol(df)
  y = rep(NA, nc + 1)
  ndiv = ceiling(n / split)
  out = NULL
  for (i in 1 : split){
    j = (i - 1) * ndiv
    x1 = df[(j + 1) : (j + ndiv), ]
    x1$n = (j + 1) : (j + ndiv)
    while (nrow(x1) < ndiv) x1 = rbind(x1, y)
    colnames(x1) = c(colnames(df), "n")
    if (is.null(out)) 
      out = x1[,c(nc + 1, 1 : nc)]
    else
      out = cbind(out, x1[,c(nc + 1, 1 : nc)])
  }
  return(out)
}  
      
 
round_df <- function(df, digits) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  
  df[,nums] <- round(df[,nums], digits = digits)
  
  (df)
}
## logNormal(0,1) LogNormal (0,0.25)
## cauchy(0,1)
## normal(0,1)

## https://github.com/stan-dev/example-models/tree/master/knitr/bradley-terry


###
###
###  Predictive tests
###

run_alternative_predictive <- function(debug=F){
  if (debug) n=1 else n=5
  ss1 = proc_all_predictive(predictive_all_experiment(n=2*n, nalg=5, ndb1 = 20, ndb2 = 10))
  mm1 = proc_all_predictive(predictive_all_experiment(n=2*n, nalg=10, ndb1 = 50, ndb2 = 25))
  sl1 = proc_all_predictive(predictive_all_experiment(n=n, nalg=5, ndb1 = 100, ndb2 = 30))
  return(list(ss=ss1, mm=mm1, sl= sl1))
}


run_predictive_test <- function(debug=F){
  if (debug) n=1 else n=5
  ss1 = proc_predictive(predictive_experiment(n=2*n, nalg=5, ndb1 = 20, ndb2 = 10, lrope1=TRUE,ties1="s",pair=T))
  mm1 = proc_predictive(predictive_experiment(n=2*n, nalg=10, ndb1 = 50, ndb2 = 25, lrope1=TRUE,ties1="s",pair=T))
  sl1 = proc_predictive(predictive_experiment(n=n, nalg=5, ndb1 = 100, ndb2 = 30, lrope1=TRUE,ties1="s",pair=T))
  return(list(ss=ss1, mm=mm1, sl= sl1))
}

matches1 <- function(ord, tablex){
  tab = print_wintable(tablex)
  pairs = strsplit(ord, " > ")
  n = length(pairs)
  out = data.frame(name = character(n) , win1 = numeric(n), win2 = numeric(n), ties = numeric(n))
  for (i in seq_along(ord)){
    x = pairs[[i]]
    if (any(tab$alg1 == x[1] & tab$alg2 == x[2])) {
      p = which(tab$alg1 == x[1] & tab$alg2 == x[2])
      out[i,1] = ord[i]
      out[i,2] = tab$win1[p]
      out[i,3] = tab$win2[p]
      out[i,4] = tab$ties[p]
    } else {
      p = which(tab$alg1 == x[2] & tab$alg2 == x[1])
      out[i,1] = ord[i]
      out[i,2] = tab$win2[p]
      out[i,3] = tab$win1[p]
      out[i,4] = tab$ties[p]
    }
  }
  return(out)
}


predictive_experiment <- function(n=10, nalg=5, ndb1 = 20, ndb2 = 20, lrope1=TRUE,ties1="s",pair=T){
    algs = colnames(ll)[-1]
    dbs = unique(ll$db)
    out = NULL
    for (v in 1:n){
        inx = sample(dbs,ndb1)
        noinx = sample(setdiff(dbs,inx), ndb2)
        thisalg = c("db",sample(algs, nalg))
        mod = bbtcomp(ll[ll$db %in% inx, thisalg], lrope = lrope1, deal_with_ties = ties1, paired = pair)
        pred = summary_pwin(mod, columns = c("high","low","above.50","above.rope", "in.rope"), hdi=0.90)
        pred2 = summary_pwin(mod, columns = c("high","low"),hdi = 0.70)[,-1]
        colnames(pred2) = c("h70","l70")
        pred3 = summary_pwin(mod, columns = c("high","low"),hdi = 0.50)[,-1]
        colnames(pred3) = c("h50","l50")
        
        res = make_wintable(ll[ll$db %in% noinx, thisalg],
                               lrope = F, deal_with_ties = "d" ) 
        res2 = matches1(pred$name, res)
        res2$pwin = res2$win1/(res2$win1 + res2$win2)
        
        aa = cbind(pred, pred2, pred3, res2[,-1])
        if (is.null(out)) out = aa else out = rbind(out,aa)
    }
  return(out)
}

predictive_all_experiment <- function(n=10, nalg=5, ndb1 = 20, ndb2 = 20){
  algs = colnames(ll)[-1]
  dbs = unique(ll$db)
  out = NULL
  for (v in 1:n){
    inx = sample(dbs,ndb1)
    noinx = sample(setdiff(dbs,inx), ndb2)
    thisalg = c("db",sample(algs, nalg))
    train = ll[ll$db %in% inx, thisalg]
    test = ll[ll$db %in% noinx, thisalg]
    res = make_wintable(test, lrope = F, deal_with_ties = "d" )
    resx = make_wintable(test, lrope = T, deal_with_ties = "d" ) 
    for (ties in  c("s", "f", "d","a"))
      for (lropex in list(c(FALSE, FALSE), c(TRUE, FALSE), c(TRUE,TRUE))) {
        lrope = lropex[1]
        pair = lropex[2]
        mod = bbtcomp(train, lrope = lrope, deal_with_ties = ties, paired = pair)
        pred = summary_pwin(mod, columns = c("mean", "low","high","above.50","above.rope", "in.rope"), hdi=0.90)
        pred2 = summary_pwin(mod, columns = c("low","high"),hdi = 0.70)[,-1]
        colnames(pred2) = c("l70","h70")
        pred3 = summary_pwin(mod, columns = c("low","high"),hdi = 0.50)[,-1]
        colnames(pred3) = c("l50","h50")
        res2 = matches1(pred$pair, res)
        res2x = matches1(pred$pair, resx)
        res2$pwin = res2$win1/(res2$win1 + res2$win2)
        res2x$pwin = res2x$win1/(res2x$win1 + res2x$win2)
        colnames(res2x) = paste0(colnames(res2x),"x")
        aa = cbind(data.frame(policy = ties, lrope = lrope, paired = pair), pred, pred2, pred3, res2[,-1],  res2x[,-1])
        if (is.null(out)) out = aa else out = rbind(out,aa)
  }}
  return(out)
}


proc_predictive <- function(res,zzz){
  r90 = sum(res$pwin >= res$low & res$pwin <= res$high)/nrow(res)
  r70 = sum(res$pwin >= res$l70 & res$pwin <= res$h70)/nrow(res)
  r50 = sum(res$pwin >= res$l50 & res$pwin <= res$h50)/nrow(res)
  dh = res$high-res$mean
  dl = res$mean-res$low
  r2x90 = sum(res$pwin >= (res$mean-3*dl) & res$pwin <= (res$mean+3*dh))/nrow(res)
  above = sum(res$pwin > res$high)/nrow(res)
  below = sum(res$pwin < res$low)/nrow(res)
  err = sum(res$mean - res$pwin)/nrow(res)
  mad = median(abs(res$mean-res$pwin))
  res = mutate(res,x7 = above.50 < 0.7,  x9 = above.50 < 0.9 & !x7, x10 = above.50 >= 0.9)
  r7 = filter(res, x7) %>% summarise(prob = 0.7, pred7 = sum(above.50), real7 = sum(win1>win2))
  #r8 = filter(res, x8) %>% summarise(prob = 0.8,pred = sum(above.50), real = n(win1>win2))
  r9 = filter(res, x9) %>% summarise(prob = 0.9,pred9 = sum(above.50), real9 = sum(win1>win2))
  r10 = filter(res, x10) %>% summarise(prob = 1.0,pred10 = sum(above.50), real10 = sum(win1>win2))
  return(data.frame(within.90 = r90, within.70 = r70, within.50= r50, above = above, below = below,
                    err = err, mad = mad,
              pred70 =r7$pred7, real70=r7$real7,
              pred90 =r9$pred9, real90=r9$real9,
              pred100 =r10$pred10, real100=r10$real10,
              within.2x90=r2x90))
}

proc_predictivex <- function(res,zzz){
  r90 = sum(res$pwinx >= res$low & res$pwinx <= res$high)/nrow(res)
  r70 = sum(res$pwinx >= res$l70 & res$pwinx <= res$h70)/nrow(res)
  r50 = sum(res$pwinx >= res$l50 & res$pwinx <= res$h50)/nrow(res)
  above = sum(res$pwinx > res$high)/nrow(res)
  below = sum(res$pwinx < res$low)/nrow(res)
  err = sum(res$mean - res$pwinx)/nrow(res)
  mad = median(abs(res$mean-res$pwinx))
  res = mutate(res,x7 = above.50 < 0.7,  x9 = above.50 < 0.9 & !x7, x10 = above.50 >= 0.9)
  r7 = filter(res, x7) %>% summarise(prob = 0.7, pred7 = sum(above.50), real7 = sum(win1>win2))
  #r8 = filter(res, x8) %>% summarise(prob = 0.8,pred = sum(above.50), real = n(win1>win2))
  r9 = filter(res, x9) %>% summarise(prob = 0.9,pred9 = sum(above.50), real9 = sum(win1>win2))
  r10 = filter(res, x10) %>% summarise(prob = 1.0,pred10 = sum(above.50), real10 = sum(win1>win2))
  return(data.frame(within.90 = r90, within.70 = r70, within.50= r50, 
                    above = above, below = below,
                    err = err, mad = mad,
                    pred70 =r7$pred7, real70=r7$real7,
                    pred90 =r9$pred9, real90=r9$real9,
                    pred100 =r10$pred10, real100=r10$real10))
}

proc_all_predictive <- function(res){
  r1 = res  %>% group_by(policy, lrope, paired) %>% group_modify(proc_predictive) %>% ungroup()
#  r2 = res  %>% group_by(policy, lrope, paired) %>% group_modify(proc_predictivex) %>% ungroup()
  #return(list(nolrope=r1, yeslrope=r2))
  return(list(nolrope=r1))
  }

###
### Compare test approaches
###
###

run_compare_tests <- function(debug=F){
  if (debug) nn=1 else nn=5
  ss1 = compare_tests(n=2*nn, nalg=5, ndb = 20,  lrope=TRUE,ties = "s", paired = T)
  mm1 = compare_tests(n=2*nn, nalg=10, ndb = 50,  lrope=TRUE,ties = "s", paired = T)
  sl1 = compare_tests(n=nn, nalg=5, ndb = 100,  lrope=TRUE,ties = "s", paired = T)
  ll1 = compare_tests(n=1, lrope=TRUE,ties = "s", paired = T, sample=F)
  return(list(ss=ss1, mm=mm1,sl=sl1,ll=ll1))
}

run_compare_tests80 <- function(debug=F){
  if (debug) nn=1 else nn=5
  ss1 = compare_tests(n=2*nn, nalg=5, ndb = 20,  lrope=TRUE,ties = "s", paired = T)
  mm1 = compare_tests(n=2*nn, nalg=10, ndb = 50,  lrope=TRUE,ties = "s", paired = T)
  sl1 = compare_tests(n=nn, nalg=5, ndb = 100,  lrope=TRUE,ties = "s", paired = T)
  ll1 = compare_tests(n=1, lrope=TRUE,ties = "s", paired = T, sample=F)
  return(list(ss=ss1, mm=mm1,sl=sl1,ll=ll1))
}

splitname <- function(tab){
  aa = strsplit(tab$name," > ")
  alg1 = unlist(lapply(aa, function(x) x[1]))
  alg2 = unlist(lapply(aa, function(x) x[2]))
  tob = dplyr::select(tab,-name)
  tob$alg1 = alg1
  tob$alg2 = alg2
  tob
}

ordnem <- function(bbt, nem) {
  nem2 = splitname(nem$tab)
  if (all(bbt$alg1 == nem2$alg1)) return(nem2)
  out = nem2
  for (i in 1:nrow(bbt)){
    a1 = bbt$alg1[i]
    a2 = bbt$alg2[i]
    out$alg1[i] = a1
    out$alg2[i] = a2
    if (any(nem2$alg1 == a1 & nem2$alg2 == a2)){
      p = which(nem2$alg1 == a1 & nem2$alg2 == a2)
      out$dif[i] = nem2$dif[p]
    } else {
      p = which(nem2$alg1 == a2 & nem2$alg2 == a1)
      out$dif[i] = nem2$dif[p]
      if (out$dif[i] == "yes") out$dif[i] = "no"
    }
  }
  return(out)
}

ordpairwilc <- function(bbt, nem) {
  nem2 = splitname(nem$tab)
  if (all(bbt$alg1 == nem2$alg1)) return(nem2)
  out = nem2
  for (i in 1:nrow(bbt)){
    a1 = bbt$alg1[i]
    a2 = bbt$alg2[i]
    out$alg1[i] = a1
    out$alg2[i] = a2
    if (any(nem2$alg1 == a1 & nem2$alg2 == a2)){
      p = which(nem2$alg1 == a1 & nem2$alg2 == a2)
      out$p.value[i] = nem2$p.value[p]
    } else {
      p = which(nem2$alg1 == a2 & nem2$alg2 == a1)
      out$p.value[i] = - nem2$p.value[p]
    }
  }
  return(out)
}

ordbayeswilc <- function(bbt, nem) {
  nem2 = splitname(nem$tab)
  if (all(bbt$alg1 == nem2$alg1)) return(nem2)
  out = nem2
  for (i in 1:nrow(bbt)){
    a1 = bbt$alg1[i]
    a2 = bbt$alg2[i]
    out$alg1[i] = a1
    out$alg2[i] = a2
    if (any(nem2$alg1 == a1 & nem2$alg2 == a2)){
      p = which(nem2$alg1 == a1 & nem2$alg2 == a2)
      out$above0[i] = nem2$above0[p]
      out$in.rope[i] = nem2$in.rope[p]
      out$above.rope[i] = nem2$above.rope[p]
    } else {
      p = which(nem2$alg1 == a2 & nem2$alg2 == a1)
      out$above0[i] = 1-nem2$above0[p]
      out$in.rope[i] = nem2$in.rope[p]
      out$above.rope[i] = 1 - (nem2$above.rope[p] + nem2$in.rope[p])
    }
  }
  return(out)
}

join_compare_tests <- function(tab, meanx = FALSE, lrope, ties, paired){
  if (meanx) 
    m = bbtcomp(tab, lrope=F, paired = F)
  else
    m = bbtcomp(tab, lrope=lrope, deal_with_ties = ties, paired=paired)
  aa = m$model$summary("beta",c("mean"))
  thisord = m$wintable$alg_names[order(aa$mean, decreasing = T)]
  tab = tab %>% group_by(db) %>% 
    summarise(across(everything(), mean, na.rm = T), .groups = "drop") %>%
    dplyr::select(-db)
  nm = procNemenyi(tab)
  if (all(thisord == nm$order$algorithm)) onem = T else onem = F
  pairw = pairwilcox(tab)
  if (all(thisord == pairw$order$algorithm)) opw = T else opw = F
  bwil = pairbayes(tab)
  if (all(thisord == bwil$order)) obw = T else obw = F
  tn = ordnem(m, nm)
  tp = ordpairwilc(m,pairw)
  tb = ordbayeswilc(m,bwil)
  oo = cbind(summary_pwin(m, columns=c("above.50","in.rope","above.rope","mean")),
             tn[,1,drop=FALSE],  tp[,1,drop=FALSE],
             round(tb[,1:3],2),  order.nem = onem, order.pw = opw, order.bpw = obw)
  colnames(oo)[8] = "in.rope2"
  colnames(oo)[9] = "above.rope2"
  oo
}


proc_compare_tests <- function(res){
  oo = res %>% summarize(npairs = n(),
    nembetter = sum(above.50>0.95 & dif != "yes"),
    nemworse =  sum(above.50< 0.95 & dif == "yes"),
    wilbetter = sum(above.50>=.95 &  p.value > 0.05),
    wilworse =  sum(above.50<.95 &  p.value <= 0.05),
    bayshigher = sum(above.50 > above0),
    bayslower =  sum(above.50 < above0)
  )
  g1 = ggplot(res,aes(x=above.50, y=above0))+geom_point(size=0.4)+
    geom_abline(intercept = 0, slope = 1)+ coord_cartesian(xlim = c(0.25,1), ylim=c(0.25,1))
  return(list(oo, g1))
}

proc_compare_tests80 <- function(res, tresh = 0.8){
  oo = res %>% summarize(npairs = n(),
                         nem.bbtbetter = sum(mean>=tresh & dif != "yes"),
                         nem.bbtworse =  sum(mean< tresh & dif == "yes"),
                         nem.bbtsame = sum(mean>=tresh & dif == "yes"),
                         wil.bbtbetter = sum(mean>= tresh &  p.value > 0.05),
                         wil.bbtworse =  sum(mean< tresh &  p.value <= 0.05),
                         wil.bbtsame = sum(mean>= tresh &  p.value <= 0.05),
                         bsr.bbthigher = sum(mean > above.rope2),
                         bsr.bbtslower =  sum(mean < above.rope2)
  )
  # g1 = ggplot(res,aes(x=mean, y=above0))+geom_point(size=0.4)+
  #   geom_abline(intercept = 0, slope = 1)+ coord_cartesian(xlim = c(0.25,1), ylim=c(0.25,1))
  return(list(oo ))
}



## compare BBT with other tests (nemenyi, pair wilcoxon, bayesian wilcoxon)
compare_tests <-function(n=10, nalg=5, ndb = 20,  lrope=TRUE,ties = "s", paired = T, sample=T){
  algs = colnames(ll)[-1]
  dbs = unique(ll$db)
  out = NULL
  if (!sample) n=1
  for (v in 1:n){
    inx = sample(dbs,ndb)
    thisalg = c("db",sample(algs, nalg))
    if (sample) 
      res2 = join_compare_tests(ll[ll$db %in% inx, thisalg], F, lrope,  ties, paired)
    else
      res2 = join_compare_tests(ll, F, lrope,  ties, paired)

    if (is.null(out)) out = res2 else out = rbind(out,res2)
  }
  return(out)
}


###
###
###  Test ties
###
###

run_alternatives_ppc <- function(debug=F){
  if (debug) n=1 else n=5
  outss = data.frame()
  outmm = data.frame()
  outsl = data.frame()
  outll = data.frame()
  for (cc in list(c(FALSE,FALSE), c(TRUE,FALSE), c(TRUE,TRUE))) {
    lrope = cc[1]
    paired = cc[2]
    ss1 = test_ties(n=2*n, nalg=5, ndb=20, lrope=lrope, paired=paired)
    outss = rbind(outss,ss1)
    mm1 = test_ties(n=2*n, nalg=10, ndb=50, lrope=lrope, paired=paired)
    outmm = rbind(outmm,mm1)
    sl1 = test_ties(n=n, nalg=5, ndb=100, lrope=lrope, paired=paired)
    outsl = rbind(outsl,sl1)
    ll1 = test_ties(sample=F, lrope=lrope, paired=paired)
    outll = rbind(outll,ll1)
  }
  outss = outss %>% group_by(lrope, paired, policy) %>% summarise(across(everything(), mean), .groups = "drop") %>% arrange(policy,lrope)
  outmm = outmm %>% group_by(lrope, paired, policy) %>% summarise(across(everything(), mean), .groups = "drop") %>% arrange(policy,lrope)
  outsl = outsl %>% group_by(lrope, paired, policy) %>% summarise(across(everything(), mean), .groups = "drop") %>% arrange(policy,lrope)
  outll = outll %>% arrange(policy,lrope)
  return(list(ss = outss, mm= outmm, sl=outsl, ll = outll))
}



test_ties <- function(n=10, nalg=5, ndb = 20, lrope=T , paired = T, sample=T){
  algs = colnames(ll)[-1]
  dbs = unique(ll$db)
  out = data.frame()
  if (!sample) n=1
  for (v in 1:n){
    inx = sample(dbs,ndb)
    thisalg = c("db",sample(algs, nalg))
    for (dd in c("s","f","a","d")) {
      if (sample)
        m = bbtcomp(ll[ll$db %in% inx, thisalg],deal_with_ties = dd, use_log_lik = T, lrope=lrope,
                     paired = paired)
      else 
        m = bbtcomp(ll,deal_with_ties = dd, use_log_lik = T, lrope=lrope,
                    paired = paired)
      w = get_waic(m)
      aw = ppc_summary(m)
      ww = aw$proportion
      if (dd=="d") ww = (ww +aw$ties)/2
      #lo = get_loo(m)
      res2 = data.frame(lrope = lrope, paired= paired, policy = dd, 
                        waic = w$estimates[3,1], h50 = ww[1], 
                        h90 = ww[2], h95 = ww[3], h100 = ww[4])
      out = rbind(out,res2)
    }
  }
  out = arrange(out,policy)
  return(out)
}


###
###
###  Test hyper
###
###

run_test_hyper <- function(debug=F){
  chname <- function(df) {
    df$hyper = case_when(df$hyper==1 ~ "lognormal",
                                              df$hyper==2 ~"cauchy",
                                              df$hyper==3 ~"normal")
    df}
  if (debug) n=1 else n=5
  ss1 = chname(test_hyper(n=2*n,nalg=5,ndb=20))
  ss1 = ss1 %>% group_by(hyper, scale) %>% summarise(across(everything(),mean),.groups = "drop")
  mm1 = chname(test_hyper(n=2*n,nalg=10,ndb=50))
  mm1 = mm1 %>% group_by(hyper, scale) %>% summarise(across(everything(),mean),.groups = "drop")
  sl1 = chname(test_hyper(n=n,nalg=5,ndb=100))
  sl1 = sl1 %>% group_by(hyper, scale) %>% summarise(across(everything(),mean),.groups = "drop")
  out = NULL
  for (dd in list(c(1,0.5), c(1,0.25), c(1,1), c(2,1), c(2,0.5), c(3,2), c(3,5))) {
    hy = dd[1]
    sc = dd[2]
    m = bbtcomp(ll, use_log_lik = T, lrope=T, hyper_prior = hy, scale = sc)
    w = get_waic(m)
    aw = ppc_summary(m)
    ww = aw$proportion
    res2 = c(hyper = hy, scale = sc,  waic = w$estimates[3,1],  h50 = ww[1], 
             h90 = ww[2], h95 = ww[3], h100 = ww[4])
    if (is.null(out)) out = res2 else out = rbind(out,res2)
  }
  out = as.data.frame(out)
  rownames(out) = NULL
  ll1 = chname(out)
  return(list(ss= ss1, mm = mm1, sl = sl1, ll = ll1))
}


test_hyper <- function(n=10, nalg=5, ndb = 20){
  algs = colnames(ll)[-1]
  dbs = unique(ll$db)
  out = data.frame()
  for (v in 1:n){
    inx = sample(dbs,ndb)
    thisalg = c("db",sample(algs, nalg)) 
    for (dd in list(c(1,0.5), c(1,0.25), c(1,1), c(2,1), c(2,0.5), c(3,2), c(3,5))) {
      hy = dd[1]
      sc = dd[2]
      m = bbtcomp(ll[ll$db %in% inx, thisalg], use_log_lik = T, lrope=T, hyper_prior = hy, scale = sc)
      w = get_waic(m)
      aw = ppc_summary(m)
      ww = aw$proportion
      res2 = data.frame(hyper = hy, scale = sc,  waic = w$estimates[3,1],  h50 = ww[1], 
                        h90 = ww[2], h95 = ww[3], h100 = ww[4])
      out = rbind(out,res2)
    }}
  rownames(out) = NULL
  return(out)
}


###
###  Tests PPC
### 

run_test_ppc <- function(debug=F){
  if (debug) n = 1 else n = 5
  ss1 = test_ppc(n=2*n,nalg=5,ndb=20,lrope=F,paired=F)
  ss2 = test_ppc(n=2*n,nalg=5,ndb=20,lrope=T,paired=F)
  ss3 = test_ppc(n=2*n,nalg=5,ndb=20,lrope=F,paired=T)
  mm1 = test_ppc(n=2*n,nalg=10,ndb=50,lrope=F,paired=F)
  mm2 = test_ppc(n=2*n,nalg=10,ndb=50,lrope=T,paired=F)
  mm3 = test_ppc(n=2*n,nalg=10,ndb=50,lrope=F,paired=T)
  sl1 = test_ppc(n=n,nalg=5,ndb=100,lrope=F,paired=F)
  sl2 = test_ppc(n=n,nalg=5,ndb=100,lrope=T,paired=F)
  sl3 = test_ppc(n=n,nalg=5,ndb=100,lrope=F,paired=T)
  xss = cbind(ss1,ss2[,-1],ss3[,-1])
  colnames(xss) = c("hdi","m.nolrope","sd.nolrope","m.nopair","sd.nopair","m.pair","sd.pair")
  xmm = cbind(mm1,mm2[,-1],mm3[,-1])
  colnames(xmm) = c("hdi","m.nolrope","sd.nolrope","m.nopair","sd.nopair","m.pair","sd.pair")
  xsl = cbind(sl1,sl2[,-1],sl3[,-1])
  colnames(xsl) = c("hdi","m.nolrope","sd.nolrope","m.nopair","sd.nopair","m.pair","sd.pair")
  ll1 = ppc_summary(bbtcomp(ll,lrope=F,paired=F))
  ll2 = ppc_summary(bbtcomp(ll,lrope=T,paired=F))
  ll3 = ppc_summary(bbtcomp(ll,lrope=T,paired=T))
  xll = cbind(ll1,ll2[,-1,drop=T],ll3[,-1,drop=T])
  colnames(xll) = c("hdi","p.nolrope","p.nopair","m.pair")
  return(list(ss=xss, mm = xmm, sl = xsl, ll = xll))
}

test_ppc <- function(n=10, nalg=5, ndb = 20, lrope = T, paired = T, ties = "s", sample = T){
  algs = colnames(ll)[-1]
  dbs = unique(ll$db)
  ww = NULL
  out = NULL
  if (!sample) n=1
  for (v in 1:n){
    inx = sample(dbs,ndb)
    thisalg = c("db",sample(algs, nalg))
    if (sample) 
      m = bbtcomp(ll[ll$db %in% inx, thisalg],  lrope=lrope, paired = paired, use_log_lik = T, 
                deal_with_ties = ties)
    else
      m = bbtcomp(ll, lrope=lrope, paired = paired, use_log_lik = T, 
                  deal_with_ties = ties)
    res2 = ppc_summary(m)
    w = get_waic(m)$estimates[3,1]
    if (is.null(ww)) ww = w else ww = c(ww,w)
    if (is.null(out)) out = res2 else out = cbind(out,res2[,-1,drop=FALSE])
  }
  x1 = data.frame(hdi = as.character(out$hdi), mean = round(apply(out[,-1, drop=F],1,mean),2), 
                    sd = round(apply(out[,-1, drop=F],1,sd),2))
  x2 = rbind(x1, list(hdi="waic", mean = round(mean(ww),2), sd = round(sd(ww),2)))
  return(x2)
}

plottests <- function(tests){
  library(ggplot2)
  library(gridExtra)
  xss = tests$ss %>% select(above.50, above0) %>% rename(BBT= above.50, BSR = above0)
  xmm = tests$mm %>% select(above.50, above0) %>% rename(BBT= above.50, BSR = above0)
  xsl = tests$sl %>% select(above.50, above0) %>% rename(BBT= above.50, BSR = above0)
  xll = tests$ll %>% select(above.50, above0) %>% rename(BBT= above.50, BSR = above0)
  c1 = cor(xss$BBT,xss$BSR)
  c2 = cor(xmm$BBT,xmm$BSR)
  c3 = cor(xsl$BBT,xsl$BSR)
  c4 = cor(xll$BBT,xll$BSR)
  g1 = ggplot(xss, aes (x=BBT, y = BSR)) + geom_point(size=0.4) + geom_abline()+ theme_bw()+ggtitle("ss")+
    annotate("text", x = 0.7, y =0.90,  label=paste("R =",round(c1,2)))
  g2 = ggplot(xmm, aes (x=BBT, y = BSR)) + geom_point(size=0.4) + geom_abline()+ theme_bw()+ggtitle("mm")+
    annotate("text", x = 0.7, y =0.90,  label=paste("R =",round(c2,2)))
  g3 = ggplot(xsl, aes (x=BBT, y = BSR)) + geom_point(size=0.4) + geom_abline()+ theme_bw()+ggtitle("sl")+
    annotate("text", x = 0.7, y =0.90,  label=paste("R =",round(c3,2)))
  g4 = ggplot(xll, aes (x=BBT, y = BSR)) + geom_point(size=0.4) + geom_abline()+ theme_bw()+ggtitle("ll")+
    annotate("text", x = 0.7, y =0.90,  label=paste("R =",round(c4,2)))
  arrangeGrob(g1,g2,g3,g4, ncol=2)
}

