library(CRFutil)
library(CRFutilRcppComponents)
library(rbenchmark)
library(microbenchmark)

# Make up a random graph
num.nodes <- 12
g <- erdos.renyi.game(num.nodes, 1, typ="gnp")
dev.off()
plot(g)

# Get its adjacency matrix and genrate an MRF sample
adj <- as.matrix(as_adj(g))

f0       <- function(y){ as.numeric(c((y==1),(y==2)))}
rmod     <- make.empty.field(adj.mat = adj, parameterization.typ = "standard", plotQ = F)

# "true" theta
rmod$par <- runif(rmod$n.par,-1.5,1.5)
rmod$par

# Make true pots from true theta
out.pot <- make.pots(parms = rmod$par,  crf = rmod,  rescaleQ = T, replaceQ = T)
rmod$edges
rmod$node.pot
rmod$edge.pot

# So now sample from the model as if we obtained an experimental sample:
num.samps <- 500
samps <- sample.exact(rmod, num.samps)
colnames(samps) <- 1:ncol(samps)
mrf.sample.plot(samps)

configs.and.counts <- as.data.frame(ftable(data.frame(samps)))
head(configs.and.counts)
configs <- configs.and.counts[,1:ncol(samps)]
configs <- as.matrix(sapply(configs, as.numeric))
head(configs)

plot(1:nrow(configs.and.counts),configs.and.counts[,ncol(configs.and.counts)],typ="h", ylab="config freqs",xlab="config #")

# Here we remove the extra index put in by CRF using the C function
theta.pars   <- fix_node_and_edge_par(node_par = rmod$node.par, edge_par = rmod$edge.par)

# Tests:
tesf <- make.empty.field(adj.mat = adj, parameterization.typ = "standard", plotQ = F)
tesf$par


#jridx <- sample(1:nrow(configs), size = 1)
for(i in 1:nrow(configs)) {
  jridx <- i
  
  for(j in 1:num.nodes) {
    
    # R version:
    ra <- symbolic.conditional.energy(config             = as.matrix(configs[jridx,]), 
                                      condition.element.number = 1, 
                                      #crf                      = tesf,
                                      crf                      = rmod,
                                      ff                       = f0,
                                      format                   = "conditional.phi",
                                      printQ                   = F)
    # C version:
    ca<- alpha_vector(config              = as.matrix(configs[jridx,]), 
                      condition_element_number = 1, 
                      node_par                 = theta.pars$node_par, 
                      edge_par                 = theta.pars$edge_par, 
                      edge_mat                 = rmod$edges,
                      adj_nodes                = rmod$adj.nodes)
    #print(paste(jridx, " ", j))
    if(sum(ra!=ca) != 0 ){
      stop("Ack, not the same!")
    }
  }
}



# R version:
ra <- symbolic.conditional.energy(config             = as.matrix(configs[jridx,]), 
                            condition.element.number = 1, 
                            #crf                      = tesf,
                            crf                      = rmod,
                            ff                       = f0,
                            format                   = "conditional.phi",
                            printQ                   = F)
# C version:
ca<- alpha_vector(config              = as.matrix(configs[jridx,]), 
             condition_element_number = 1, 
             node_par                 = theta.pars$node_par, 
             edge_par                 = theta.pars$edge_par, 
             edge_mat                 = rmod$edges,
             adj_nodes                = rmod$adj.nodes)

rbind(ra,ca)
sum(ra!=ca)

# Speed tests
benchmark(replications = 1000,
          # C version:
          ca<- alpha_vector(config              = as.matrix(configs[jridx,]), 
                            condition_element_number = 1, 
                            node_par                 = theta.pars$node_par, 
                            edge_par                 = theta.pars$edge_par, 
                            edge_mat                 = rmod$edges,
                            adj_nodes                = rmod$adj.nodes),
          # R version:
          ra <- symbolic.conditional.energy(config             = as.matrix(configs[jridx,]), 
                                            condition.element.number = 1, 
                                            #crf                      = tesf,
                                            crf                      = rmod,
                                            ff                       = f0,
                                            format                   = "conditional.phi",
                                            printQ                   = F)
)

# C
start_time <- Sys.time()
for(i in 1:nrow(configs)) {
  jridx <- i
  
  for(j in 1:num.nodes) {
    # C version:
    ca<- alpha_vector(config              = as.matrix(configs[jridx,]), 
                      condition_element_number = 1, 
                      node_par                 = theta.pars$node_par, 
                      edge_par                 = theta.pars$edge_par, 
                      edge_mat                 = rmod$edges,
                      adj_nodes                = rmod$adj.nodes)
    #print(paste(jridx, " ", j))
    if(sum(ra!=ca) != 0 ){
      stop("Ack, not the same!")
    }
  }
}
end_time <- Sys.time()
end_time - start_time

# R
start_time <- Sys.time()
for(i in 1:nrow(configs)) {
  jridx <- i
  
  for(j in 1:num.nodes) {
    # R version:
    ra <- symbolic.conditional.energy(config             = as.matrix(configs[jridx,]), 
                                      condition.element.number = 1, 
                                      #crf                      = tesf,
                                      crf                      = rmod,
                                      ff                       = f0,
                                      format                   = "conditional.phi",
                                      printQ                   = F)
    #print(paste(jridx, " ", j))
    if(sum(ra!=ca) != 0 ){
      stop("Ack, not the same!")
    }
  }
}
end_time <- Sys.time()
end_time - start_time


get_par_off(
  config      = as.matrix(configs[jridx,]), 
  i_in        =  1, 
  j_in        = -1, 
  node_par_in = theta.pars$node_par, 
  edge_par_in = NULL, 
  edge_mat_in = array(-1,c(1,1)), 
  printQ      = F)

array(-1,c(1,1))

get.par.idx(
  config   = as.matrix(configs[jridx,]), 
  i        = 1, 
  #j        = NULL, 
  node.par = rmod$node.par, 
  #edge.par = rmod$edge.par, 
  #edge.mat = rmod$edges, 
  ff       = f0, 
  printQ   = F)
