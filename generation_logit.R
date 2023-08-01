# this is the code to create the istances without dependency, i.e. Logi-based uncertainty

# Loading useful packages
library(evd)
require(evd)

# loop with 10 different seeds to create 10 instances per each level of oscillation 
for (seed in 1:10){
  set.seed(seed)
  n_nodes=100
  n_alternatives=10
  n_scenarios=100

  # Choosing the level of oscillation  and  Generating the deterministic costs: 
  
  #determistic_cost=matrix(runif(n_nodes*n_alternatives, 0,10),nrow =n_nodes,ncol = n_alternatives); # low 30% 
  #determistic_cost=matrix(runif(n_nodes*n_alternatives, 0,3),nrow =n_nodes,ncol = n_alternatives); # medium 50%
  determistic_cost=matrix(runif(n_nodes*n_alternatives, 0,0.45),nrow =n_nodes,ncol = n_alternatives); # high 80%
  

  # laction costs
  f_i=matrix(runif(n_nodes*n_alternatives, 0/3,0.5/3),nrow =1,ncol = n_nodes); 
  
  
  # Setting the Euler constant 
  Gamma=0.5772;
  
  # computing the teoritecal costs h_i
  h_i=array(0,c(1,n_nodes));
  beta=runif(1,min = 1,max = 2);
  for (id in 1:n_nodes){
    h_i[1,id]= (1/beta)*(log( sum( exp (beta*determistic_cost[id,] )))+Gamma)
  }
  
  # Generating the scenarios.
  cost=array(0,c(n_nodes,n_alternatives,n_scenarios))
  max_oscillation=matrix(0,nrow=n_scenarios,ncol=1)
  perc=rep(0,n_scenarios)
  
  for (s in 1:n_scenarios){
    oscillations=matrix(rgumbel(n_nodes*n_alternatives, loc = 0, scale = beta),nrow = n_nodes,ncol = n_alternatives)
    cost[,,s]=determistic_cost+oscillations;
    perc[s]=mean(abs(determistic_cost)/(abs(determistic_cost)+abs(oscillations)))
  }

  # Saving the created instance 
  output=c(n_nodes,n_alternatives,n_scenarios)
  instance_name=paste("instance_high_oscillation","_seed_",as.character(seed),".csv",sep = "")
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





