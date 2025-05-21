library(profvis)
profvis(source("code/net_n.R"))
profvis(source("code/net_n2.R"))
profvis(source("code/net_n3.R"))
profvis(source("code/motif_prob_par.R"))


profvis(source("code/grp_pred.R"))
profvis(source("code/old_py_pred.R"))



profvis(source("code/alpha_density.R"))






library(microbenchmark)
microbenchmark(source("code/net_n.R"), times = 5)


# Comparison collapse package
microbenchmark(grp_pred(5, 0, 1:10000, 10000),
               benchmark_py_pred(5, 0, 1:10000, 10000),
               times = 10000)

microbenchmark(GRPN(1:10000, expand = F),
               tabulate(1:10000),
               times = 10000)

# comparison
microbenchmark(net_n0(1000,c(3,5),c(5,5)),
               net_n1(1000,c(3,5),c(5,5)),
               net_n2(1000,c(3,5),c(5,5)),
               times = 1000)
# > microbenchmark(net_n0(1000,c(3,5),c(5,5)),
#                  +                net_n1(1000,c(3,5),c(5,5)),
#                  +                net_n2(1000,c(3,5),c(5,5)),
#                  +                times = 1000)
# Unit: milliseconds
# expr     min      lq     mean   median       uq      max neval
# net_n0(1000, c(3, 5), c(5, 5)) 30.9604 40.7896 47.50220 46.46710 50.47360 153.6788  1000
# net_n1(1000, c(3, 5), c(5, 5)) 30.4900 41.1322 47.32339 45.03855 50.61535 199.4748  1000
# net_n2(1000, c(3, 5), c(5, 5)) 27.2179 35.3380 40.74169 37.66075 44.05700 117.6281  1000

microbenchmark(net_n2(1000,c(3,5),c(5,5)),
               net_n3(1000,c(3,5),c(5,5)),
               times = 5000)



microbenchmark(sim(100,seq(1000,10000,1000) , c(3,5,10), c(5,5)),
               sim_par2(100,seq(1000,10000,1000) , c(3,5,10), c(5,5)),
               times = 1)


microbenchmark(source("code/motif_prob.R"),
               source("code/motif_prob_par.R"),
               times = 1)



microbenchmark(source("code/test.R"),
               times = 5)
