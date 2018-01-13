BNNsel <- function(X,Y,train_num = as.integer(0.8*N), hid_num = 3, lambda=0.025, total_iteration = 1000000, popN=20, nCPUs = 20 )
{
  N <- dim(X)[1]
  P <- dim(X)[2]
  test_num <- N-train_num
  Y <- as.matrix(Y)
  OUT_UNIT <- dim(Y)[2]
  DataX=DataY=NULL;
  for(i in 1:N)
  {
    for(j in 1:P) DataX[(i-1)*P+j]=X[i,j];
    for(j in 1:OUT_UNIT) DataY[(i-1)*OUT_UNIT+j]=Y[i,j];
  }
 # system("ulimit -s unlimited")
  .C("posratio", as.numeric(DataX),as.numeric(DataY),
     as.integer(train_num),as.integer(test_num),
     as.integer(OUT_UNIT),as.integer(P),as.integer(hid_num),as.numeric(lambda),
     as.integer(total_iteration),as.integer(popN),as.integer(nCPUs))
  net <- as.numeric(read.table("aa.net"))
  prob <- as.numeric(read.table("aa.BF"))[1]
  mar <- as.numeric(read.table("aa.mar"))
  fit <- as.numeric(read.table("aa.fit")[,1])
  pred <- as.numeric(read.table("aa.pred")[,1])
  result <- list(net=net,prob=prob,mar=mar,fit=fit,pred=pred)
  unlink("aa.net")
  unlink("aa.BF")
  unlink("aa.mar")
  unlink("aa.fit")
  unlink("aa.pred")
  return(result)
}


BNNprior <- function(dimX, dimY, hid_num = 3,lambda=0.025, total_iteration = 1000000, popN = 20)
{
  .C("pratio", as.integer(dimY),as.integer(dimX),as.integer(hid_num),as.numeric(lambda),as.integer(total_iteration),as.integer(popN))
  prob <- as.numeric(read.table("bb.BF"))[1]
  unlink("bb.BF")
  return(prob)
}


## in linux set C stage useage : $ulimit -s unlimited