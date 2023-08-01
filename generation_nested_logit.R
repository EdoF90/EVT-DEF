# this is the code to create the istances with dependence 

# Loading usefull Libraries
library(evd)
require(evd)


Compute_G <- function(x,beta_m,nest_id) {
  C_nest=c();
  for (n in unique(nest_id)){
    C_nest[n]=(sum((x[nesut_id==n])^beta_m[n]))^(1/beta_m[n]);
  }
  Val=sum(C_nest);
  return(Val)
}


# loop with 10 different seeds to create 10 instances per each level of oscillation 
for (seed in 1:10){

  #setting the parameters
  set.seed(seed)
  n_nodes=100
  n_alternatives=10
  n_scenarios=100
  n_nest=2

  # Choosing the level of oscillation  and  Generating the deterministic costs: 

  determistic_cost=matrix(runif(n_nodes*n_alternatives, 0,5),nrow =n_nodes,ncol = n_alternatives); # low 30% 
  #determistic_cost=matrix(runif(n_nodes*n_alternatives, 0,1.5),nrow =n_nodes,ncol = n_alternatives); # medium 50%
  #determistic_cost=matrix(runif(n_nodes*n_alternatives, 0,0.3),nrow =n_nodes,ncol = n_alternatives); # high 80%

  # Setting the laction costs
  f_i=matrix(runif(n_nodes*n_alternatives, 0/3,0.5/3),nrow =1,ncol = n_nodes); 

  # Setting the Euler constant 
  Gamma=0.5772;

  # computing the teoritecal costs h_i
  h_i=array(0,c(1,n_nodes));
  beta_m=sample(2:5, n_nest, replace = FALSE);
  fisrt_nest_length=sample(2:(n_alternatives-2), 1);
  second_nest_length=n_alternatives-fisrt_nest_length;

  nest_1=1:fisrt_nest_length;
  nest_2=(fisrt_nest_length+1):n_alternatives;

  nest_id=rep(0,n_alternatives);
  nest_id[nest_1]=rep(1,fisrt_nest_length)
  nest_id[nest_2]=rep(2,second_nest_length)

  for (id in 1:n_nodes){
    # note: beta is 1.
    h_i[1,id]=(log(Compute_G (exp(determistic_cost[id,]),beta_m,nest_id))+Gamma) 
  }

  # Generating the scenarios.
  t=list()
  pos=0
  pos_nest_of_length_n=rep(0,n_alternatives)
  for (i in 1:n_alternatives){
    pos_nest_of_length_n[i]=  pos_nest_of_length_n[i]+choose(n_alternatives,i)
      for (j in 1:choose(n_alternatives,i)){
        pos=pos+1
        t[[pos]]=rep(0,i);
      }
  }

  nest_of_first_length=combn(1:n_alternatives,fisrt_nest_length)
  nest_of_second_length=combn(1:n_alternatives,second_nest_length)

  pos_first_nest=0;
  for (id in 1:choose(n_alternatives,fisrt_nest_length)){
    pos_first_nest=pos_first_nest+1
    if (sum(nest_of_first_length[,id]==nest_1)==fisrt_nest_length){
      break
    }
  }

  pos_second_nest=0;
  for (id in 1:choose(n_alternatives,second_nest_length)){
    pos_second_nest=pos_second_nest+1
    if (sum(nest_of_second_length[,id]==nest_2)==second_nest_length){
      break
    }
  } 

  start_id_first_next=sum(pos_nest_of_length_n[1:(fisrt_nest_length-1)])
  start_id_second_next=sum(pos_nest_of_length_n[1:(second_nest_length-1)])

  t[[start_id_first_next+pos_first_nest]]=rep(1,fisrt_nest_length)
  t[[start_id_second_next+pos_second_nest]]=rep(1,second_nest_length)

  dep=c( rep ( 1/beta_m[1], sum(pos_nest_of_length_n[2:(fisrt_nest_length)]) ),  rep ( 1/beta_m[2], sum(pos_nest_of_length_n[(fisrt_nest_length+1):n_alternatives]) ) )

  cost=array(0,c(n_nodes,n_alternatives,n_scenarios))
  max_oscillation=matrix(0,nrow=n_scenarios,ncol=1)

  perc=rep(0,n_scenarios)
  for (s in 1:n_scenarios){
    d=rmvevd(n_nodes, dep = dep, asy = t, model = "alog", d = n_alternatives,mar=c(0,1,-1))
    oscillations=-log(1-d)
    cost[,,s]=determistic_cost+oscillations;
    perc[s]=mean(abs(determistic_cost)/(abs(determistic_cost)+abs(oscillations)))
  }

  # Writting the generated instance on file 
  output=c(n_nodes,n_alternatives,n_scenarios)
  instance_name=paste("instance_low_oscillation","_seed_",as.character(seed),".csv",sep = "")
  write.table(output, file = paste("./data/",instance_name,sep=""), append = FALSE, quote = FALSE, sep = " ", 
  eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, 
  qmethod = c("escape", "double"), fileEncoding = "")

  write.table(f_i, file = paste("./data/",instance_name,sep=""), append = TRUE, quote = FALSE, sep = " ", 
  eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, 
  qmethod = c("escape", "double"), fileEncoding = "")

  write.table(h_i, file = paste("./data/",instance_name,sep=""), append = TRUE, quote = FALSE, sep = " ", 
  eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, 
  qmethod = c("escape", "double"), fileEncoding = "")

  for   (id in 1:n_scenarios){
    write.table(cost[,,id], file = paste("./data/",instance_name,sep=""), append = TRUE, quote = FALSE, sep = " ", 
    eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, 
    qmethod = c("escape", "double"), fileEncoding = "")
  }
}





