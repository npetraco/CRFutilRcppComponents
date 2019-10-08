library(CRFutil)
library(CRFutilRcppComponents)
library(rbenchmark)
library(microbenchmark)

# Make up a random graph
g <- erdos.renyi.game(16, 1, typ="gnp")
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

theta.pars   <- fix_node_and_edge_par(node_par = rmod$node.par, edge_par = rmod$edge.par)
jridx <- 1
jridx <- sample(1:nrow(configs), size = 1)
#configs[jridx,]
rmod$edges
as.matrix(configs[jridx,])

get_par_off(as.matrix(configs[jridx,]), i_in = 17, node_par_in = theta.pars$node_par) # This should not work

get_par_off(as.matrix(configs[jridx,]), 
            i_in = 14,
            j_in = 15,
            #node_par_in  = theta.pars$node_par, 
            edge_par_in  = theta.pars$edge_par,
            edge_mat_in = rmod$edges) + 1

get.par.idx(as.matrix(configs[jridx,]), 
            i = 14,
            j = 15,
            #node.par = rmod$node.par, 
            edge.par = rmod$edge.par,
            edge.mat = rmod$edges, ff = f0)

emptym <- array(-1,c(1,1)) # This indicates empty node_par_in 
get_par_off2(as.matrix(configs[jridx,]), 
            i_in = 14,
            j_in = 15,
            node_par_in  = emptym,
            edge_par_in  = theta.pars$edge_par,
            edge_mat_in = rmod$edges) + 1


# Check edges:
num.rand.configs <- 2000
jridxs <- sample(1:nrow(configs), size = num.rand.configs, replace = F)

for(jii in 1:length(jridxs)) {
  
  jridx <- jridxs[jii]
  a.config <- configs[jridx,]
  
  for(jedge.idx in 1:nrow(rmod$edges)) {
    
    jedg <- rmod$edges[jedge.idx, ]
    cidx <- get_par_off2(as.matrix(a.config), 
                        i_in = jedg[1],
                        j_in = jedg[2],
                        node_par_in  = emptym, 
                        edge_par_in  = theta.pars$edge_par,
                        edge_mat_in = rmod$edges) + 1
    
    ridx <- get.par.idx(a.config, 
                        i = jedg[1],
                        j = jedg[2],
                        #node.par = rmod$node.par, 
                        edge.par = rmod$edge.par,
                        edge.mat = rmod$edges, ff = f0)
    
    if(cidx != ridx) {
      print(jedg)
      print(cidx)
      print(ridx)
      print(a.config)
      stop("Ack!")
    }
    
  }
  
}


# Speed check:
jridx    <- sample(1:nrow(configs), size = 1)
a.config <- configs[jridx,]
jedg     <- rmod$edges[sample(1:nrow(rmod$edges), size = 1), ]
#a.config
#jedg

benchmark(replications = 1000,
          coff <- get_par_off2(as.matrix(a.config), 
                              i_in = jedg[1],
                              j_in = jedg[2],
                              node_par_in  = emptym, 
                              edge_par_in  = theta.pars$edge_par,
                              edge_mat_in = rmod$edges),
          
          ridx <- get.par.idx(a.config, 
                              i = jedg[1],
                              j = jedg[2],
                              #node.par = rmod$node.par, 
                              edge.par = rmod$edge.par,
                              edge.mat = rmod$edges, ff = f0)
)

