rm(list = ls()) # clear all parameters from memory

MyResults = data.frame(LPprob = character(),
                      n = integer(),
                      m = integer(),
                      suc = logical(),
                      n.binding = integer(),
                      Rank0 = numeric(),
                      Rank1 = numeric(),
                      Rank2 = numeric(),
                      Rank3 = numeric(),
                      Rank4 = numeric(),
                      Rank5 = numeric(),
                      Rank6 = numeric(),
                      Rank7 = numeric(),
                      Rank8 = numeric(),
                      Rank9 = numeric(),
                      Rank10 = numeric())

rankPci2021 = function(a1, b1, C) {
  require(DescTools)
  n = ncol(a1)
  m1 = nrow(a1)
  k1.const = matrix(, nrow = n, ncol = m1)
  k1.rank = matrix(, nrow = n, ncol = m1)
  syndiasmoi = matrix(, nrow = 2, ncol = n-1)
  mymaxobj = which(C==max(C))[1]
  syndiasmoi[1,] = rep(mymaxobj, n-1)
  syndiasmoi[2,] = seq(1,n,1)[-mymaxobj]
  a1[a1==0]= 10^(-64)
  k1.const = -a1[,mymaxobj] / a1
  k1.rank = apply(k1.const, 2, rank)
#  v = apply(k1.rank, 1, FUN = Gmean, na.rm = T)
#  v = apply(k1.rank, 1, FUN = min, na.rm = T)
  v = apply(k1.rank, 1, FUN = sum, na.rm = T)
  v = rank(v)
  return(v)
}

require(R.matlab)
#require(linprog)
require(lpSolve)
require(DescTools)
require(stringr)
require(summarytools)
require(ggplot2)

MAT.names = str_sub(list.files("LP_MATLAB/."), start = 1, end = -5)
MPS.names = list.files("MPS/.")

Netlib.names = intersect(MAT.names, gsub('-', '_', MPS.names))

i=0
NumofProbs = length(Netlib.names)
for (Testing in 1:NumofProbs) {
i = i + 1
TestName = Netlib.names[i]
p1 = readMat(paste('LP_MATLAB/', TestName, '.mat', sep = ''))
if (names(p1)[1]=="Problem") {
  myA = as.matrix(p1$Problem[[3]])
  myb = as.vector(p1$Problem[[4]])
  myc = as.vector(p1$Problem[[5]])
} else {
  myA = p1$A
  myb = p1$b
  myc = -p1$c
}
n =  ncol(myA)
m = nrow(myA)
print(paste('i=', i, '/', length(Netlib.names), 'Prob:',TestName, ': n=', n, ' m=', m))


# calculate Rankings
Rankings = rankPci2021(myA, myb, myc)
#print(Rankings)

#Solve with simplex - linprog::solveLP
#Classic = solveLP(-p1$c, p1$b, p1$A, maximum = T)
#desmeytikoi = abs(Classic$con$actual - p1$b) <= 10^(-12)

mytxt = read.table(file=paste("MPS/", gsub('_', '-', TestName), sep = ""), sep="\n", quote="", comment.char="")
mytxt = mytxt[seq(3, m+3,1),1]
mytxt = substr(mytxt,2,2)
mytxt = mytxt[-which(mytxt=='N')]
myconst = mytxt
myconst[myconst=='E'] = '='
myconst[myconst=='G'] = '>='
myconst[myconst=='L'] = '<='

Classic1 = lpSolve::lp(direction = 'min', 
                      objective.in = -myc,
                      const.mat = as.matrix(myA),
                      const.dir = myconst, 
                      const.rhs = myb,
                      all.int = F)
desmeytikoi = as.vector(abs(myA %*% Classic1$solution - myb) <= 10^(-12))
print(paste('Success:', Classic1$status==0))

#sum(desmeytikoi)
#fivenum(Rankings[desmeytikoi]/max(Rankings))
#quantile(Rankings[desmeytikoi]/max(Rankings), probs = seq(0,1,0.1))

MyResults[i,1] = TestName
MyResults[i,2] = n
MyResults[i,3] = m
MyResults[i,4] = Classic1$status==0
MyResults[i,5] = sum(desmeytikoi)/m
MyResults[i,6:16] = quantile(Rankings[desmeytikoi]/max(Rankings), probs = seq(0,1,0.1))
write.table(MyResults, file = "NetLibResults.csv", sep = ',', row.names=F)
}

print(summarytools::dfSummary(MyResults), method = "browser")

Pos_of_Binding = 10*round(MyResults$n.binding/MyResults$m,1)
Score = 0
for (i in 1:nrow(MyResults)) {
  Score[i] = MyResults[i, 5+Pos_of_Binding[i]]
}

ggplot(data.frame(Pos_of_Binding, Score), aes(x = Pos_of_Binding, y = Score)) + geom_smooth() + geom_point()


