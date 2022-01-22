#### run once ####
#install.packages("lpSolve")
#install.packages("lpSolveAPI")
#install.packages("plyr")
##############################
# library(lpSolve)
# library(lpSolveAPI)
# library(plyr)priny

library(linprog)
library(summarytools)
library(rpart)
library(rpart.plot)
library(formattable)
library(crayon)
#library(boot)

rm(list = ls()) # clear all parameters from memory

compares = data.frame(comparing = logical(),
                      ClassicSolved = integer(),
                      NewSolved = integer(),
                      n = integer(),
                      m = integer(),
                      used.constraints = integer(),
                      n.binding = integer(),
                      prop.binding = numeric(),
                      rnd1 = integer(),
                      rnd2 = integer(),
                      rnd3 = integer(),
                      rnd4 = integer(),
                      rnd5 = integer(),
                      rnd6 = integer(),
                      rnd7 = integer(),
                      rnd8 = integer(),
                      rnd9 = integer(),
                      rnd10 = integer())

sel.per.rand.total = data.frame(Total = numeric(),
                                Srnd1 = numeric(),
                                Srnd2 = numeric(),
                                Srnd3 = numeric(),
                                Srnd4 = numeric(),
                                Srnd5 = numeric(),
                                Srnd6 = numeric(),
                                Srnd7 = numeric(),
                                Srnd8 = numeric(),
                                Srnd9 = numeric(),
                                Srnd10 = numeric())

Max.selection = data.frame(n.desmeytikoi = numeric(),
                           m = numeric(),
                                Nrnd1 = numeric(),
                                Nrnd2 = numeric(),
                                Nrnd3 = numeric(),
                                Nrnd4 = numeric(),
                                Nrnd5 = numeric(),
                                Nrnd6 = numeric(),
                                Nrnd7 = numeric(),
                                Nrnd8 = numeric(),
                                Nrnd9 = numeric(),
                           Nrnd10 = numeric())
# to do ...
# καλύτερο ήταν το 6o > 1o > 3o & 4o > 8o
# 1o 77% , 
# 2o 57% , 
# 3o 40%, 
# 4o 40%, 
# 3o-4o 81%
# 5o 89%, 
# 6o 88%
# 7o, 68%
# 8o, 74%

randsTF = c(F, F, F, F, F, F, F, F, T, F)
SimplexSuccess = 0
h = 0
TotalProblems = 0

rangeN = 3:21
#rangeN = 40:100
mysuccess = 0

set.seed(1000) # mesw ayths ths entolhs eksasflizetai oti tha exoume to idio provlhma kathe fora pou trexoume ton kwdika se periptwsh otan dinoume ton idio arithmo sto seed...
while ( SimplexSuccess < 100 ){
#  h = h+1
  TotalProblems = TotalProblems + 1
  n = sample(rangeN, 1, T) # dhmiourgia enos plithous metavlhtwn apo 2 ews 10 (tyxaio to poses akrivws tha einai)...
  m1 = sample(seq(10*n, 15*n,1),1,T)  # posoi tha einai oi periorismoi (o arithmos aytos tha kymainetai apo n+1 mexri 3*n...)
  # ### calculate c , coefficients for objective function
  c = runif(n, min = -1, max = 1) # edw dhmiourgoume tyxaious syntelestes ths obj.function pou einai osoi kai to plhthos twn metavlhtwn kai oi times tous kymainonde apo 0 ews 1... 
  # 
  # ### calculate a, coefficients of m x n matrix, m number
  # ### of constraints and n number of variables with the same
  # ### slope with objective function
  a1 = runif(m1*n, min = 0, max = 1)
  a1 = matrix(a1, nrow = m1, ncol = n)
  
  # 
  # 
  # 
  # ### calculate random solution x
  xsol = runif(n, min = 0, max = 1)
  # 

  # my.k = sort(sample(seq(2, m1-1), n-1)) # briskoume tixaious ari8mous kai tous diatasoume
  # my.k = c(my.k,m1)
  # for (i in seq(1,n-1)){
  #   a1[seq(my.k[i]+1,my.k[i+1]), n-i] = -a1[seq(my.k[i]+1,my.k[i+1]), n-i]
  # }
  
  for (k in 1:m1) {
    a1[k, ] = sample(c(-1,1), n, replace = T) *  a1[k, ]
  }

  # my.signs = matrix(1,nrow = 2^n, ncol = n)
  # for(j in 1:n){
  #   my.signs[,j] = rep(rep(c(1,-1),each=2^(n-j)),2^(j-1))
  #   
  # }
  # 
  # my.k = sort(sample(seq(2, m1-1), 2^n-1, replace = T)) # briskoume tyxaious ari8mous kai tous diatasoume
  # my.k = c(1,my.k,m1+1)
  # j = 0
  # for (i in 1:2^n){
  #   j = j+1
  #   multi = my.signs[j,]
  #   istart = my.k[i]
  #   iend = my.k[i+1]-1
  #   for(k in seq(istart,iend,1)){
  #     a1[k,]=multi * a1[k,]
  #   }
  # }
  
  # ### find b for constraints with the same slope with objective function
  b1 = a1 %*% xsol
  for (i in 1:m1) {   # elafra diataraxi sta b
    b1[i] = jitter(b1[i], amount = 0.1*abs(b1[i]))
  }
  
#  b1 = matrix(sort(b1,decreasing = T),nrow = m1)  # f8inousa sera oi sintelestes tou b
  
  pinakas.const = a1 # o pinakas me tous syntelestes twn x1,x2,x3 gia kathe ena periorismo...
  sta8eroi.oroi  = b1
  sintelestes.obj = c
  dimension = n
  Initial.no.const = m1
  
  monadiaios = diag(n) # gia na valoume kai tous periorismous xi>=0
  
  statheroi.oroi.xi = matrix(0,n,1)
  Classic = solveLP(sintelestes.obj, sta8eroi.oroi, pinakas.const, maximum = T)

#    Classic = boot:: simplex(a = sintelestes.obj, A1 = pinakas.const, b1 = sta8eroi.oroi,  
#                           A2 = monadiaios, b2 = statheroi.oroi.xi, A3 = NULL,b3 = NULL, maxi = TRUE, 
#                           n.iter = n + 2 * m1, eps = 1e-6)
  
#  SimplexSuccess = SimplexSuccess + as.numeric(Classic$solved==1)
  SimplexSuccess = SimplexSuccess + as.numeric(Classic$status==0)
  
#  print(paste(h, SimplexSuccess))
  
  if (Classic$status==0) {
    h = SimplexSuccess
    desmeytikoi = abs(Classic$con$actual - b1) <= 0.00000001
  ############################################################################################################ 
  #for (j in combn(3,2))
  syndiasmoi = combn(n,2) 
  frow = c(syndiasmoi[1, ], syndiasmoi[2,]) 
  lrow = c(syndiasmoi[2, ], syndiasmoi[1,])
  syndiasmoi = rbind(frow, lrow)
  mymaxobj = which(c==max(c))
 tohold = syndiasmoi[1,]==mymaxobj
#  tohold = syndiasmoi[1,]==mymaxobj | syndiasmoi[2,]==mymaxobj
 syndiasmoi = syndiasmoi[, tohold]
  myl = length(syndiasmoi[1,])
  k1.const = matrix(, nrow = m1, ncol = myl) # dhmiourgia enos pinaka me tous syntelestes twn xi kathe periorismou
  k1.obj = 1:1
  for (i in 1:myl){ # gia osoi einai oi syndiasmoi
    sel = syndiasmoi[,i]
    
    # Evresi klisewn constraints
    
    for(k1 in 1:m1){
      k1.const[k1,i] = -a1[k1,sel[1]]/a1[k1,sel[2]]
    }
    k1.obj[i]  = -c[sel[1]]/c[sel[2]]
  }
  
  k1.ext = rbind(k1.const,k1.obj)  
  seldesm = rep(FALSE, m1) # dhmiourgei ena keno dianysma
  #### αθροιστής για κάθε rand
  sel.per.rand = data.frame(Total = seldesm,
                            sR1 = seldesm,
                            sR2 = seldesm,
                            sR3 = seldesm,
                            sR4 = seldesm,
                            sR5 = seldesm,
                            sR6 = seldesm,
                            sR7 = seldesm,
                            sR8 = seldesm,
                            sR9 = seldesm,
                            sR10 = seldesm)
  # for (k in 1:m1) {
  #   seldesm[k] = sum(sign(a1[k,])==rep(1,n))==n
  # }

#    seldesm[1:my.k[1]] = T # 

# rand1
  apo = 2
  eos = n 
  mymymin = min(apo, eos)
  mymymax = max(apo, eos)
  myrand1 = sample(seq(mymymin, mymymax), size = 1)
  Nrand1 = m1
  for(i in 1:myl){
    sumrank = k1.ext[,i]
    randchoice = rank(abs(sumrank[1:(length(sumrank)-1)]))
    rank.const.obj = randchoice <= myrand1
    if (max(randchoice[desmeytikoi])<Nrand1) Nrand1 = max(randchoice[desmeytikoi])
    if (randsTF[1]==T) sel.per.rand$sR1 = sel.per.rand$sR1 | rank.const.obj
  }

  stathmikos.mesos = abs(b1/rowSums(a1))
  # rand7
  apo = 2
  eos = n
  mymymin = min(apo, eos)
  mymymax = max(apo, eos)
  myrand7 = sample(seq(mymymin, mymymax), size = 1)
  Nrand7 = m1
  for(i in 1:myl){
    sumrank = k1.ext[,i] * c(stathmikos.mesos, median(stathmikos.mesos))
    randchoice = rank(abs(sumrank[1:(length(sumrank)-1)]))
    rank.const.obj = randchoice <= myrand7
    if (max(randchoice[desmeytikoi])<Nrand7) Nrand7 = max(randchoice[desmeytikoi])
    if (randsTF[7]==T) sel.per.rand$sR7 = sel.per.rand$sR7 | rank.const.obj
  }
  
#rand6 ### eggitita synteleston a_ij me c_j
  apo = 2
  eos = m1
  mymymin = min(apo, eos)
  mymymax = max(apo, eos)
  myrand6 = sample(seq(mymymin, mymymax), size = 1)
  Nrand6 = m1
  Aext = rbind(a1, c)
  for (nstiles in 1:n) {
    sumrank = rank(Aext[,nstiles])
    sumrank = sumrank - tail(sumrank,1)
    randchoice = rank(abs(sumrank[1:(length(sumrank)-1)]))
    rank.const.obj = randchoice <= myrand6
    if (max(randchoice[desmeytikoi])<Nrand6) Nrand6 = max(randchoice[desmeytikoi])
    if (randsTF[6]==T) sel.per.rand$sR6 = sel.per.rand$sR6 | rank.const.obj
  }

  #rand10 ### eggitita synteleston a_ij me c_j fou diaireso me bi
  apo = 2
  eos = n
  mymymin = min(apo, eos)
  mymymax = max(apo, eos)
  myrand10 = sample(seq(mymymin, mymymax), size = 1)
  Nrand10 = m1
  Abext = rbind(a1, c)
  for (nstiles in 1:n) {
    sumrank = rank(Abext[,nstiles]/ c(b1,1))
    sumrank = sumrank - tail(sumrank,1)
    randchoice = rank(abs(sumrank[1:(length(sumrank)-1)]))
    rank.const.obj = randchoice <= myrand10
    if (max(randchoice[desmeytikoi])<Nrand10) Nrand10 = max(randchoice[desmeytikoi])
    if (randsTF[10]==T) sel.per.rand$sR10 = sel.per.rand$sR10 | rank.const.obj
  }
  
  # rand9 ### synoliki eggitita synteleston a_ij me c_j
  apo = n+1
  eos = m1
  mymymin = min(apo, eos)
  mymymax = max(apo, eos)
  myrand9 = sample(seq(mymymin, mymymax), size = 1)
  Nrand9 = m1
  Aext = rbind(a1, c)
  cmult = rank(c) - 1
#  cmult[c<0] = 0
  tsumrank = rep(0, length(Aext[,1]))
  for (nstiles in 1:n) {
    sumrank = Aext[,nstiles]
    sumrank = rank(abs(sumrank - tail(sumrank,1)))
    tsumrank = tsumrank + sumrank*cmult[nstiles]
  }
  randchoice = rank(abs(tsumrank[1:(length(tsumrank)-1)]))
  rank.const.obj = randchoice <= myrand9
  if (max(randchoice[desmeytikoi])<Nrand9) Nrand9 = max(randchoice[desmeytikoi])
  if (randsTF[9]==T) sel.per.rand$sR9 = sel.per.rand$sR9 | rank.const.obj
  
  # rand5
  apo = 2 
  eos = n # %/% 2 #(n)%/%2
  mymymin = min(apo, eos)
  mymymax = max(apo, eos)
  myrand5 = sample(seq(mymymin, mymymax), size = 1)
  Nrand5 = m1
  for(i in 1:myl){
    sumrank = k1.ext[,i]
    randchoice = rank(abs(sumrank[1:(length(sumrank)-1)]))
    rank.const.obj = randchoice >= m1 - myrand5 # good
#    rank.const.obj = rank(abs(sumrank[1:(length(sumrank)-1)])) <= myrand5 # bad
    if (1+m1-min(randchoice[desmeytikoi])<Nrand5) Nrand5 = 1+m1 - min(randchoice[desmeytikoi])
    if (randsTF[5]==T) sel.per.rand$sR5 = sel.per.rand$sR5 | rank.const.obj
  }

  # rand8
  apo = 2 # %/% 2
  eos = n # %/% 4 # (n)%/%2
  mymymin = min(apo, eos)
  mymymax = max(apo, eos)
  myrand8 = sample(seq(mymymin, mymymax), size = 1)
  Nrand8 = m1
  for(i in 1:myl){
    sumrank = k1.ext[,i] * c(stathmikos.mesos, median(stathmikos.mesos))
    randchoice = rank(abs(sumrank[1:(length(sumrank)-1)]))
    rank.const.obj = randchoice >= m1 - myrand8
    if (1+m1-min(randchoice[desmeytikoi])<Nrand8) Nrand8 = 1+m1 - min(randchoice[desmeytikoi])
    if (randsTF[8]==T) sel.per.rand$sR8 = sel.per.rand$sR8 | rank.const.obj
  }
  
  k1.ext.rank = k1.ext
  for(i in 1:myl){
    toupdate = rank(k1.ext.rank[,i])
    k1.ext.rank[,i]=toupdate - toupdate[m1+1]
  }
  sumrank = rowSums(k1.ext.rank)
  sumrank = sumrank - sumrank[length(sumrank)]
  rank.const.obj = rank(abs(sumrank[1:(length(sumrank)-1)]))
  randchoice = rank.const.obj

  # rand2
  Nrand2 = max(randchoice[desmeytikoi])
  apo = n + 1 # n%/%3
  eos = m1 # (2*n)%/%3
  mymymin = min(apo, eos)
  mymymax = max(apo, eos)
  myrand2 = sample(seq(mymymin, mymymax), size = 1)
  selupddesm = rank.const.obj <= myrand2  
  if (randsTF[2]==T) sel.per.rand$sR2 = sel.per.rand$sR2 | selupddesm
    
  stathmikos.mesos = abs(b1/rowSums(a1))
  ranksm = rank(stathmikos.mesos)
  Nrand3 = min(ranksm[desmeytikoi])
  Nrand4 = max(ranksm[desmeytikoi])
  
# rand3 and rand4
  apo = n+1 # (2*n)%/%3
  eos = m1 # (3*n)%/%2
  mymymin = min(apo, eos)
  mymymax = max(apo, eos)
  myrand3 = sample(seq(mymymin, mymymax), size = 1)
  apo = n+1 #(5*n)%/%4
  eos = m1 # (8*n)%/%4
  mymymin = min(apo, eos)
  mymymax = max(apo, eos)
  myrand4 = sample(seq(mymymin, mymymax), size = 1)
  st.mes.mean = median(stathmikos.mesos, na.rm = T)
  randchoice = 
  myselection = c(head(sort(ranksm[stathmikos.mesos-st.mes.mean>=0],decreasing = F),myrand3))
  theseis = integer()
  for(i in 1:length(myselection)){
    theseis = c(theseis,
                which(ranksm == myselection[i]))
  }
  if (randsTF[3]==T) sel.per.rand$sR3[theseis]= T 
  
  myselection = c(head(sort(ranksm[stathmikos.mesos-st.mes.mean<=0],decreasing = T),myrand4))
  theseis = integer()
  for(i in 1:length(myselection)){
    theseis = c(theseis,
                which(ranksm == myselection[i]))
  }
  if (randsTF[4]==T) sel.per.rand$sR4[theseis]= T 
  
  # άθροισμα εφαρμοζόμενων περιορισμών
  seldesm = rowSums(sel.per.rand[,-1]) > 0
  sel.per.rand$Total = seldesm
    
  a1.desm = a1[seldesm, ]
  b1.desm = b1[seldesm]


  new.pinakas.const = a1.desm
  new.sta8eroi.oroi  = b1.desm
  New.no.const = length(b1.desm)
  
  New = solveLP(sintelestes.obj, b1.desm, a1.desm, maximum = T)
#  Classic = solveLP(sintelestes.obj, sta8eroi.oroi, pinakas.const, maximum = T)
  #  New = boot::simplex(a = sintelestes.obj, A1 = a1.desm, b1 = b1.desm, A2 = monadiaios, 
#                      b2 = statheroi.oroi.xi, A3 = NULL,
#                      b3 = NULL, maxi = TRUE, n.iter = n + 2 * New.no.const, eps = 1e-6)
  # print("Simplex with New Method:")
  # print(New)
  
#  Classic = boot:: simplex(a = sintelestes.obj, A1 = pinakas.const, b1 = sta8eroi.oroi,  
#                           A2 = monadiaios, b2 = statheroi.oroi.xi, A3 = NULL,b3 = NULL, maxi = TRUE, 
#                           n.iter = n + 2 * m1, eps = 1e-6)
  
  # cat("Dimension:", n)
  # 
  # cat("No. Used Constraints:", New.no.const)
  # cat("No. Initial Constraints:", Initial.no.const)
  # print(Classic$soln)
  # print(New$soln)
  
  my.norm = norm(as.matrix(Classic$solution- New$solution), type = c("F"))
  compare.sol = my.norm <= 0.0001
  No.Used.const = New.no.const
  
  desmeytikoi = abs(Classic$con$actual - b1) <= 0.00000001
  binding = sum(desmeytikoi & seldesm)
  prop.binding = binding / sum(desmeytikoi)
  nstiles = length(sel.per.rand.total)
  for (i in 1:nstiles) sel.per.rand.total[h, i] = sum(desmeytikoi & sel.per.rand[,i]) / sum(desmeytikoi)
  
  compares[h,] = c(comparing = compare.sol,
                   ClassicSolved = Classic$status,
                   NewSolved = New$status,
                   n = n,
                   m = m1,
                   used.constraints = New.no.const,
                   n.binding = binding,
                   prop.binding = prop.binding,
                   rnd1 = myrand1,
                   rnd2 = myrand2,
                   rnd3 = myrand3,
                   rnd4 = myrand4,
                   rnd5 = myrand5,
                   rnd6 = myrand6,
                   rnd7 = myrand7,
                   rnd8 = myrand8,
                   rnd9 = myrand9,
                   rnd10 = myrand10)
  
  Max.selection[h,] = c(sum(desmeytikoi), m1, 
                             Nrand1,
                             Nrand2,
                             Nrand3,
                             Nrand4,
                             Nrand5,
                             Nrand6,
                             Nrand7,
                             Nrand8, 
                             Nrand9,
                        Nrand10)
  
  mysuccess = mysuccess + as.numeric(compare.sol)
  toprint = paste(h, "(",TotalProblems,")",
                  "suc =", round(100*mysuccess/h,2), 
                  "(",
                  "n=",n,
                  "m=",m1,
                  "economy=",m1 - New.no.const,
                  compare.sol,
                  ")")
  if (compare.sol==T) print(cat(blue(toprint))) else print(cat(red(toprint)))
  
} # close if Simplex succeed 

} # 

boxplot(Max.selection$m - Max.selection[, -c(1:2,5:6)],horizontal = T)


#print(table(compares$comparing))
#print(table(compares$used.constraints[compares$comparing==TRUE]))
#print(summary(compares$used.constraints))

toremove = which(randsTF == F)
toremove = seq(9,18,1)[toremove]
fit = rpart(comparing ~ ., data=compares[, -c(2:8, toremove)])
fit.prune = prune(fit,cp = fit$cptable[which.min(fit$cptable[,"xerror"]),"CP"])
#formattable(as.data.frame(round(100*prop.table(fit.prune$variable.importance),2)))
as.data.frame(round(100*prop.table(fit.prune$variable.importance),2))
fit.prune
rpart.plot(fit.prune, cex =1, tweak = 1, box.palette = "BuRd", fallen.leaves = F)
title("Αύξηση πιθανότητας εύρεσης λύσης")
fit1.lm = lm(comparing ~ ., data=compares[, -c(2:8, toremove)])
fit1.lm = step(fit1.lm)
summary(fit1.lm)
economy = compares$m - compares$used.constraints
compares = data.frame(compares, economy)

fit2 = rpart(economy ~ ., data=compares[, -c(1:8, toremove)])
fit2.prune = prune(fit2,cp = fit2$cptable[which.min(fit2$cptable[,"xerror"]),"CP"])
#formattable(as.data.frame(round(100*prop.table(fit2.prune$variable.importance),2)))
as.data.frame(round(100*prop.table(fit2.prune$variable.importance),2))
fit2.prune
rpart.plot(fit2.prune, cex =1, tweak = 1, box.palette = "BuRd", fallen.leaves = F)
title("Μείωση περιορισμών")
fit2.lm = lm(economy ~ ., data=compares[, -c(1:8, toremove)])
fit2.lm = step(fit2.lm)
summary(fit2.lm)

fit3 = rpart(prop.binding ~ ., data=compares[, -c(1:7,toremove, 19)])
fit3 = rpart(prop.binding ~ ., data = data.frame(prop.binding = compares$prop.binding,
                                                 k = compares$rnd9/compares$m))
fit3.prune = prune(fit3,cp = fit3$cptable[which.min(fit3$cptable[,"xerror"]),"CP"])
#formattable(as.data.frame(round(100*prop.table(fit3.prune$variable.importance),2)))
as.data.frame(round(100*prop.table(fit3.prune$variable.importance),2))
fit3.prune
rpart.plot(fit3.prune, cex =1, tweak = 1, box.palette = "BuRd", fallen.leaves = F)
title("Πιθανότητα εντοπισμού δεσμευτικού περιορισμού")
fit3.lm = lm(prop.binding ~ ., data=compares[, -c(1:7,toremove, 19)])
fit3.lm = step(fit3.lm)
summary(fit3.lm)

#print(summarytools::dfSummary(compares), method = "viewer")
print(summarytools::dfSummary(compares[compares$comparing==1,]), method = "browser")
print(summarytools::dfSummary(sel.per.rand.total), method = "browser")
print(summarytools::dfSummary(Max.selection/Max.selection$m), method = "browser")
print(summarytools::dfSummary(Max.selection), method = "browser")

DescTools::Freq(compares$prop.binding, ord = "desc")
DescTools::Freq(compares$m - compares$used.constraints, ord = "desc")
DescTools::Freq(Max.selection$Nrnd6 - Max.selection$Nrnd9)

library(stargazer)
stargazer(fit3.lm)

#write.table(compares, file = "compares.small.2000.csv", row.names = F, sep=",")
#write.table(Max.selection, file = "Max.selection.small.2000.csv", row.names = F, sep=",")
write.table(compares, file = "compares.321.rand9.2000.csv", row.names = F, sep=",")
write.table(Max.selection, file = "Max.selection.321.rand9.2000.csv", row.names = F, sep=",")

# right.tree <- snip.rpart(fit.prune, toss = 2) # trim subtree at node 2
# rpart.plot(right.tree, cex =1, tweak = 1., box.palette = "BuRd",fallen.leaves = F)
# formattable(as.data.frame(round(100*prop.table(right.tree$variable.importance),2)))


# abs(pinakas.const %*% New$soln - sta8eroi.oroi) < 0.000000000001  
# new.sel = abs(pinakas.const %*% New$soln - sta8eroi.oroi) < 0.000000000001
# which(new.sel)        