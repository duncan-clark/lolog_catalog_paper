# This script fits ERGM and LOLOG to the data

library("statnet")
library("devtools")
library("pkgbuild")
library("lolog")

# extended LOLOG package included as tar.gz for the fitting of cyclic triples
library("LologExtension")

load("processed_data.RData")

# -------- Utility Functions for fitting and GOF procedure ------- #

# We define some utilities that allowed for standardisation accross the many examples we considered.

ergm_fit <- function(net,
                     terms,
                     gofit_terms,
                     offset.coef = NULL,
                     control.ergm = NULL){
  
  formula <- as.formula(paste(c("net ~",paste(terms,collapse = "+")),sep = ""))
  gofit_formula <- as.formula(paste("net~",paste(gofit_terms,collapse = "+"),sep = ""))
  
  if(!is.null(control.ergm)){
    model <- tryCatch(ergm(formula,control = control.ergm,offset.coef = offset.coef),error = function(e){return(e)})
  }else{
    model <- tryCatch(ergm(formula,offset.coef = offset.coef),error = function(e){return(e)})
  }
  
  if(class(model)[1] != "ergm"){return(model)}
  
  summary <- summary(model)
  gofit <- gof(model,gofit_formula)
  
  return(list(model = model,
              summary = summary,
              gofit = gofit))
}

lolog_fit <- function(net,
                      terms,
                      gofit_terms,
                      gofit_name,
                      fit_terms = NULL,
                      vertex_order = NULL,
                      reverse_order = F, #nodes are added from low to high values of vetex order, if need high to low values specify as TRUE
                      cl = NULL,
                      gofit_draws = 100,
                      nsamp = 1000,
                      maxIter = 100,
                      ...){
  t <- proc.time()
  
  formula <- paste("net ~",paste(terms,collapse = "+"),sep = "")
  gofit_formula <- paste("net ~",paste(gofit_terms,collapse = "+"),sep = "")
  
  if(!is.null(vertex_order)){
    formula <- paste( formula, " | ",eval(vertex_order),sep= "")
    gofit_formula <- paste( gofit_formula, " | ", eval(vertex_order))
    
    assign(eval(vertex_order),get.vertex.attribute(net,eval(parse(text = "vertex_order"))))
    if(reverse_order){
      tmp <- eval(parse(text = vertex_order))
      tmp <- (max(tmp) - tmp + 1)
      assign(eval(vertex_order),tmp)
    }
  }
  if(!is.null(fit_terms)){
    auxFormula <- as.formula(paste("net ~",paste(fit_terms,collapse = "+"),sep=""))
  }else{
    auxFormula <- NULL
  }
  
  formula <- as.formula(formula)
  gofit_formula <- as.formula(gofit_formula)
  
  model <- tryCatch(lolog(formula,
                          auxFormula = auxFormula,
                          cl=cl,
                          nsamp = nsamp,
                          maxIter = maxIter),
                    error = function(e){return(e)})
  
  if(class(model)[2] != "lolog"){return(model)}
  
  if(is.na(model)){
    return(list(
      model = NA,
      summary = NA,
      gofit = NA
    ))
  }
  
  summary <- summary(model,term_names = terms)
  gofit <- gofit(model,gofit_formula,nsim = gofit_draws)
  
  #for each of the terms in the gofit formula do the plot:
  plots <- list()
  gofits <- list()
  for(term in gofit_terms){
    tmp <- gofit(model,
                 as.formula(paste("net ~",term)),nsim = gofit_draws)
    gofits[[length(gofits)+1]] <- tmp
    plots[[length(plots)+1]] <- plot(tmp)
  }
  t<- proc.time() -t
  
  return(list(model = model,
              summary = summary,
              gofits = gofits,
              plots = plots,
              time = t))
}

# 0)
#Start with dyad independent model to see if this fits adequately:
ergm_fit_1_0 <- ergm_fit(net,
                         terms = c("edges",
                                   "nodeicov('handicapped')",
                                   "nodeocov('handicapped')",
                                   "nodeicov('sweetgiver')",
                                   "nodeocov('sweetgiver')",
                                   "nodeicov('rank')",
                                   "nodeocov('rank')",
                                   "nodeicov('repeater')",
                                   "nodeocov('repeater')"),
                         gofit_terms = c("edges","degree(0:25)","esp(0:15)"),
                         gofit_name = "ergm_1_0")

#See that fit isn't very good:
#Seems that rank and repeater are the main covariate to have an effect
#Possible sweetgiver also has an effct but since there is only one sweetgiver in 
#the network the parameter MLE has a high s.d. and is not significant.

#Does not fit well
#Bad fit on degree and esp 

# 1) Model 1 in paper
#Fit as in paper to replicate table 1 Markov model:
ergm_fit_1_1 <- ergm_fit(net,
                         terms = c("edges",
                                   "nodeicov('repeater_or_sweets')",
                                   "nodeicov('handicapped')",
                                   "nodeicov('rank')",
                                   "nodeocov('rank')",
                                   "absdiff('rank')",
                                   "gwesp()",
                                   "mutual",
                                   "cycle(3)",
                                   "m2star"),
                         gofit_terms = c("edges","degree(0:25)","esp(0:15)"),
                         gofit_name = "ergm_1_1")

#Comment on fit :
#same as model in paper - managed to match!
#Much improved fit on esp and min distance,slight improvement on in and out degree fit.

# 2) Model 2 in paper
# Fit as in paper with uprank statistic to match results
ergm_fit_1_2 <- ergm_fit(net,
                         terms = c("edges",
                                   "nodeicov('repeater_or_sweets')",
                                   "nodeicov('handicapped')",
                                   "nodeicov('rank')",
                                   #"nodeocov('rank')",
                                   "edgecov(edges_adj_uprank)",
                                   "absdiff('rank')",
                                   "gwesp()",
                                   "mutual",
                                   "cycle(3)",
                                   "m2star"),
                         gofit_terms = c("edges","degree(0:25)","esp(0:15)"),
                         gofit_name = "ergm_1_2")

#Similar to fit of model 1


# 3) Model 3 in paper
ergm_fit_1_3 <- ergm_fit(net,
                         terms = c("edges",
                                   "nodeicov('repeater_or_sweets')",
                                   "nodeicov('handicapped')",
                                   "nodeicov('rank')",
                                   "nodeocov('rank')",
                                   "edgecov(edges_adj_uprank)",
                                   "absdiff('rank')",
                                   "gwesp()",
                                   "mutual",
                                   "cycle(3)",
                                   "m2star"),
                         gofit_terms = c("edges","degree(0:25)","esp(0:15)"),
                         gofit_name = "ergm_1_3")

#Similar to fit of model 1

# 4) Include popularity terms this time
# suspect this might make some of the other covariate estimates not significant:
ergm_fit_1_4 <- ergm_fit(net,
                         terms = c("edges",
                                   "nodeicov('repeater_or_sweets')",
                                   "nodeicov('handicapped')",
                                   "nodeicov('rank')",
                                   "nodeocov('rank')",
                                   "edgecov(edges_adj_uprank)",
                                   "absdiff('rank')",
                                   "gwesp()",
                                   "gwidegree(fixed = FALSE)",
                                   "gwodegree(fixed = FALSE)",
                                   "mutual",
                                   "cycle(3)",
                                   "m2star"),
                         gofit_terms = c("edges","degree(0:25)","esp(0:15)"),
                         gofit_name = "ergm_1_4")

# ESP and minimum distance fit similar to model 1, slight improvement on in and out degree fit
# May suggest some popularity process not captured by model 1

# 5)
#Check to see if triangle and star model is degenerative:
ergm_5 <- function(net,MCMC.samplesize,seed=1){
  tmp <- ergm(net ~ nodeicov("rank")+
                nodeicov("repeater_or_sweets")+
                nodeicov("handicapped")+
                nodeocov("rank")+
                edgecov(edges_adj_uprank)+
                absdiff("rank")+
                ttriple+
                mutual+
                cycle(3)+
                m2star+
                istar(c(2,3,4))+
                ostar(c(2,3,4))+
                edges,
              estimate = "MLE",
              control = control.ergm(seed = seed,MCMC.samplesize =MCMC.samplesize))
  return(tmp)
}


#ergm_1_5 <- ergm_5(net,2**10,seed = 2)
#summary(ergm_1_5)
#gof_plot(ergm_1_5)
#Did not converge - ergm degenerate

# 6) Model 2 in paper but with g weighted degrees instead of triangles
ergm_fit_1_6 <- ergm_fit(net,
                         terms = c("edges",
                                   "nodeicov('repeater_or_sweets')",
                                   "nodeicov('handicapped')",
                                   "nodeicov('rank')",
                                   "nodeocov('rank')",
                                   "edgecov(edges_adj_uprank)",
                                   "absdiff('rank')",
                                   #"gwesp()",
                                   "gwidegree(fixed = FALSE)",
                                   "gwodegree(fixed = FALSE)",
                                   "mutual",
                                   "cycle(3)",
                                   "m2star"),
                         gofit_terms = c("edges","degree(0:25)","esp(0:15)"),
                         gofit_name = "ergm_1_6")

# Key is to check whether this model captures the ESP distribution, since if it does
# the ESP distribution may not be due to transitive closure.

# Captures degree distribution arguabley better than models with gwesp term
# ESP distribution fits somewhat, though issues with esp 0 and 2, so perhaps does not account for these well

#--------- LOLOG Fit ---------

# 0)
#Fit inital LOLOG model with dyad independence
lolog_fit_1_0 <- lolog_fit(net,
                         terms = c("edges",
                                   "nodeCov('handicapped_in','in')",
                                   "nodeCov('handicapped_out','out')",
                                   "nodeCov('sweetgiver_in','in')",
                                   "nodeCov('sweetgiver_out','out')",
                                   "nodeCov('rank_in','in')",
                                   "nodeCov('rank_out','out')",
                                   "nodeCov('repeater_in','in')",
                                   "nodeCov('repeater_out','out')"
                                   ),
                         gofit_terms = c("edges","degree(0:25)","esp(0:15)"),
                         gofit_name = "lolog_1_0")

#same as dyad independent ergm - fit poor as expected!

# 1)
#Add parameters in paper model 1
# lolog_fit_1_1 <- lolog_fit(net,
#                            terms = c("edges",
#                                      "nodeCov('repeater_or_sweets_in','in')",
#                                      "nodeCov('handicapped_in','in')",
#                                      "nodeCov('rank_in','in')",
#                                      "nodeCov('rank_out','out')",
#                                      "edgeCov(edges_adj_rankdiff,'rankdiff')",
#                                      #"edgeCov(edges_adj_uprank,'uprank')",
#                                      "gwesp(0.5)",
#                                      "mutual",
#                                      "cTriple",
#                                      "twoPath"),
#                            gofit_terms = c("edges","degree(0:25)","esp(0:15)"),
#                            gofit_name = "lolog_1_1")
#Comments
#This did not converge

#Replace gwesp term with triangles 
lolog_fit_1_1 <- lolog_fit(net,
                           terms = c("edges",
                                     "nodeCov('repeater_or_sweets_in','in')",
                                     "nodeCov('handicapped_in','in')",
                                     "nodeCov('rank_in','in')",
                                     "nodeCov('rank_out','out')",
                                     "edgeCov(edges_adj_rankdiff,'rankdiff')",
                                     #"edgeCov(edges_adj_uprank,'uprank')",
                                     #"gwesp(0.5)",
                                     "triangles",
                                     "mutual",
                                     "cTriple",
                                     "twoPath"),
                           gofit_terms = c("edges","degree(0:15,'out')","degree(0:15,'in')","esp(0:15)"),
                           gofit_name = "lolog_1_1")

# 2)
#Add terms in paper model 2
# lolog_fit_1_2 <- lolog_fit(net,
#                            terms = c("edges",
#                                      "nodeCov('repeater_or_sweets_in','in')",
#                                      "nodeCov('handicapped_in','in')",
#                                      #"nodeCov('rank_in','in')",
#                                      "nodeCov('rank_out','out')",
#                                      "edgeCov(edges_adj_rankdiff,'rankdiff')",
#                                      "edgeCov(edges_adj_uprank,'uprank')",
#                                      "gwesp(0.5)",
#                                      "mutual",
#                                      "cTriple",
#                                      "twoPath"),
#                            gofit_terms = c("edges","degree(0:25)","esp(0:15)"),
#                            gofit_name = "lolog_1_2")

# Comments
# This did not converge

#Replace gwesp term with triangles 
lolog_fit_1_2 <- lolog_fit(net,
                           terms = c("edges",
                                     "nodeCov('repeater_or_sweets_in','in')",
                                     "nodeCov('handicapped_in','in')",
                                     #"nodeCov('rank_in','in')",
                                     "nodeCov('rank_out','out')",
                                     "edgeCov(edges_adj_rankdiff,'rankdiff')",
                                     "edgeCov(edges_adj_uprank,'uprank')",
                                     #"gwesp(0.5)",
                                     "triangles",
                                     "mutual",
                                     "cTriple",
                                     "twoPath"),
                           gofit_terms = c("edges","degree(0:15,'out')","degree(0:15,'in')","esp(0:15)"),
                           gofit_name = "lolog_1_2")


# 3)
#Add terms in paper model 3
# lolog_fit_1_3 <- lolog_fit(net,
#                            terms = c("edges",
#                                      "nodeCov('repeater_or_sweets_in','in')",
#                                      "nodeCov('handicapped_in','in')",
#                                      "nodeCov('rank_in','in')",
#                                      "nodeCov('rank_out','out')",
#                                      "edgeCov(edges_adj_rankdiff,'rankdiff')",
#                                      "edgeCov(edges_adj_uprank,'uprank')",
#                                      "gwesp(0.5)",
#                                      "mutual",
#                                      "cTriple",
#                                      "twoPath"),
#                            gofit_terms = c("edges","degree(0:25)","esp(0:15)"),
#                            gofit_name = "lolog_1_3")

#Replace gwesp term with triangles 
lolog_fit_1_3 <- lolog_fit(net,
                           terms = c("edges",
                                     "nodeCov('repeater_or_sweets_in','in')",
                                     "nodeCov('handicapped_in','in')",
                                     "nodeCov('rank_in','in')",
                                     "nodeCov('rank_out','out')",
                                     "edgeCov(edges_adj_rankdiff,'rankdiff')",
                                     "edgeCov(edges_adj_uprank,'uprank')",
                                     #"gwesp(0.5)",
                                     "triangles",
                                     "mutual",
                                     "cTriple",
                                     "twoPath"),
                           gofit_terms = c("edges","degree(0:15,'out')","degree(0:15,'in')","esp(0:15)"),
                           gofit_name = "lolog_1_3")


#Comment


# 4)
#Lolog fitting proceedure
#Keep all covariate terms in model
#Fit just with  triangles and mutual
lolog_fit_1_4 <- lolog_fit(net,
                           terms = c("edges",
                                     "nodeCov('repeater_or_sweets_in','in')",
                                     "nodeCov('handicapped_in','in')",
                                     "nodeCov('rank_in','in')",
                                     "nodeCov('rank_out','out')",
                                     "edgeCov(edges_adj_rankdiff,'rankdiff')",
                                     "edgeCov(edges_adj_uprank,'uprank')",
                                     "triangles",
                                     "mutual"),
                           gofit_terms = c("edges","degree(0:15,'out')","degree(0:15,'in')","esp(0:15)"),
                           gofit_name = "lolog_1_4")

# 5) 
#Fit with just 2,3 stars
lolog_fit_1_5 <- lolog_fit(net,
                           terms = c("edges",
                                     "nodeCov('repeater_or_sweets_in','in')",
                                     "nodeCov('handicapped_in','in')",
                                     "nodeCov('rank_in','in')",
                                     "nodeCov('rank_out','out')",
                                     "edgeCov(edges_adj_rankdiff,'rankdiff')",
                                     "edgeCov(edges_adj_uprank,'uprank')",
                                     "star(c(2,3))",
                                     "mutual"),
                           gofit_terms = c("edges","degree(0:15,'out')","degree(0:15,'in')","esp(0:15)"),
                           gofit_name = "lolog_1_5")

# 6) 
#Fit with just 2,3,4,5 stars
lolog_fit_1_6 <- lolog_fit(net,
                           terms = c("edges",
                                     "nodeCov('repeater_or_sweets_in','in')",
                                     "nodeCov('handicapped_in','in')",
                                     "nodeCov('rank_in','in')",
                                     "nodeCov('rank_out','out')",
                                     "edgeCov(edges_adj_rankdiff,'rankdiff')",
                                     "edgeCov(edges_adj_uprank,'uprank')",
                                     "star(c(2,3,4,5))",
                                     "mutual"),
                           gofit_terms = c("edges","degree(0:15,'out')","degree(0:15,'in')","esp(0:15)"),
                           gofit_name = "lolog_1_6")


# 7) 
#Fit with triangles and 2 and 3 stars
lolog_fit_1_7 <- lolog_fit(net,
                           terms = c("edges",
                                     "nodeCov('repeater_or_sweets_in','in')",
                                     "nodeCov('handicapped_in','in')",
                                     "nodeCov('rank_in','in')",
                                     "nodeCov('rank_out','out')",
                                     "edgeCov(edges_adj_rankdiff,'rankdiff')",
                                     "edgeCov(edges_adj_uprank,'uprank')",
                                     "triangles",
                                     "star(c(2,3))",
                                     "mutual"),
                           gofit_terms = c("edges","degree(0:15,'out')","degree(0:15,'in')","esp(0:15)"),
                           gofit_name = "lolog_1_7")

# 8) 
#Fit with triangles and 2,3,4 and 5 stars
lolog_fit_1_8 <- lolog_fit(net,
                           terms = c("edges",
                                     "nodeCov('repeater_or_sweets_in','in')",
                                     "nodeCov('handicapped_in','in')",
                                     "nodeCov('rank_in','in')",
                                     "nodeCov('rank_out','out')",
                                     "edgeCov(edges_adj_rankdiff,'rankdiff')",
                                     "edgeCov(edges_adj_uprank,'uprank')",
                                     "triangles",
                                     "star(c(2,3,4,5))",
                                     "mutual"),
                           gofit_terms = c("edges","degree(0:15,'out')","degree(0:15,'in')","esp(0:15)"),
                           gofit_name = "lolog_1_8")

# 9)
#Fit with triangles, 2,3 stars and cTriples
lolog_fit_1_9 <- lolog_fit(net,
                           terms = c("edges",
                                     "nodeCov('repeater_or_sweets_in','in')",
                                     "nodeCov('handicapped_in','in')",
                                     "nodeCov('rank_in','in')",
                                     "nodeCov('rank_out','out')",
                                     "edgeCov(edges_adj_rankdiff,'rankdiff')",
                                     "edgeCov(edges_adj_uprank,'uprank')",
                                     "triangles",
                                     "star(c(2,3))",
                                     "cTriple",
                                     "mutual"),
                           gofit_terms = c("edges","degree(0:15,'out')","degree(0:15,'in')","esp(0:15)"),
                           gofit_name = "lolog_1_9")

# 10)
#Fit with triangles, 2,3 stars and twoPaths
lolog_fit_1_10 <- lolog_fit(net,
                           terms = c("edges",
                                     "nodeCov('repeater_or_sweets_in','in')",
                                     "nodeCov('handicapped_in','in')",
                                     "nodeCov('rank_in','in')",
                                     "nodeCov('rank_out','out')",
                                     "edgeCov(edges_adj_rankdiff,'rankdiff')",
                                     "edgeCov(edges_adj_uprank,'uprank')",
                                     "triangles",
                                     "star(c(2,3))",
                                     "twoPaths",
                                     "mutual"),
                           gofit_terms = c("edges","degree(0:15,'out')","degree(0:15,'in')","esp(0:15)"),
                           gofit_name = "lolog_1_10")
# 11)
#Fit with triangles, 2,3 stars and twoPaths and cTriple
lolog_fit_1_11 <- lolog_fit(net,
                           terms = c("edges",
                                     "nodeCov('repeater_or_sweets_in','in')",
                                     "nodeCov('handicapped_in','in')",
                                     "nodeCov('rank_in','in')",
                                     "nodeCov('rank_out','out')",
                                     "edgeCov(edges_adj_rankdiff,'rankdiff')",
                                     "edgeCov(edges_adj_uprank,'uprank')",
                                     "triangles",
                                     "star(c(2,3))",
                                     "cTriple",
                                     "twoPaths",
                                     "mutual"),
                           gofit_terms = c("edges","degree(0:15,'out')","degree(0:15,'in')","esp(0:15)"),
                           gofit_name = "lolog_1_11")

#Cannot have both triangles two paths

# 12)
#Fit with star(2,3), twoPaths, 
lolog_fit_1_12 <- lolog_fit(net,
                            terms = c("edges",
                                      "nodeCov('repeater_or_sweets_in','in')",
                                      "nodeCov('handicapped_in','in')",
                                      "nodeCov('rank_in','in')",
                                      "nodeCov('rank_out','out')",
                                      "edgeCov(edges_adj_rankdiff,'rankdiff')",
                                      "edgeCov(edges_adj_uprank,'uprank')",
                                      "star(c(2,3))",
                                      "twoPaths",
                                      "mutual"),
                            gofit_terms = c("edges","degree(0:15,'out')","degree(0:15,'in')","esp(0:15)"),
                            gofit_name = "lolog_1_11")

# 13)
#Fit with star(c(2,3)), and twoPaths and cTriple
lolog_fit_1_13 <- lolog_fit(net,
                            terms = c("edges",
                                      "nodeCov('repeater_or_sweets_in','in')",
                                      "nodeCov('handicapped_in','in')",
                                      "nodeCov('rank_in','in')",
                                      "nodeCov('rank_out','out')",
                                      "edgeCov(edges_adj_rankdiff,'rankdiff')",
                                      "edgeCov(edges_adj_uprank,'uprank')",
                                      "star(c(2,3))",
                                      "cTriple",
                                      "twoPaths",
                                      "mutual"),
                            gofit_terms = c("edges","degree(0:15,'out')","degree(0:15,'in')","esp(0:15)"),
                            gofit_name = "lolog_1_13")

# 14)
#Try to fit with transitive triples
lolog_fit_1_14 <- lolog_fit(net,
                            terms = c("edges",
                                      "nodeCov('repeater_or_sweets_in','in')",
                                      "nodeCov('handicapped_in','in')",
                                      "nodeCov('rank_in','in')",
                                      "nodeCov('rank_out','out')",
                                      "edgeCov(edges_adj_rankdiff,'rankdiff')",
                                      "edgeCov(edges_adj_uprank,'uprank')",
                                      #"triangles",
                                      #"star(c(2,3))",
                                      #"cTriple",
                                      #"twoPaths",
                                      "tTriple",
                                      "mutual"),
                            gofit_terms = c("edges","degree(0:15,'out')","degree(0:15,'in')","esp(0:15)"),
                            gofit_name = "lolog_1_14")

# 15)
#Try to fit with transitive triples
lolog_fit_1_15 <- lolog_fit(net,
                            terms = c("edges",
                                      "nodeCov('repeater_or_sweets_in','in')",
                                      "nodeCov('handicapped_in','in')",
                                      "nodeCov('rank_in','in')",
                                      "nodeCov('rank_out','out')",
                                      "edgeCov(edges_adj_rankdiff,'rankdiff')",
                                      "edgeCov(edges_adj_uprank,'uprank')",
                                      #"triangles",
                                      "star(c(2,3))",
                                      #"cTriple",
                                      #"twoPaths",
                                      "tTriple",
                                      "mutual"),
                            gofit_terms = c("edges","degree(0:15,'out')","degree(0:15,'in')","esp(0:15)"),
                            gofit_name = "lolog_1_15")

# 16)
#Try to fit with transitive triples
lolog_fit_1_16 <- lolog_fit(net,
                            terms = c("edges",
                                      "nodeCov('repeater_or_sweets_in','in')",
                                      "nodeCov('handicapped_in','in')",
                                      "nodeCov('rank_in','in')",
                                      "nodeCov('rank_out','out')",
                                      "edgeCov(edges_adj_rankdiff,'rankdiff')",
                                      "edgeCov(edges_adj_uprank,'uprank')",
                                      #"triangles",
                                      "star(c(2,3))",
                                      "cTriple",
                                      #"twoPaths",
                                      "tTriple",
                                      "mutual"),
                            gofit_terms = c("edges","degree(0:15,'out')","degree(0:15,'in')","esp(0:15)"),
                            gofit_name = "lolog_1_16")

# 17)
#Try to fit with transitive triples
lolog_fit_1_17 <- lolog_fit(net,
                            terms = c("edges",
                                      "nodeCov('repeater_or_sweets_in','in')",
                                      "nodeCov('handicapped_in','in')",
                                      "nodeCov('rank_in','in')",
                                      "nodeCov('rank_out','out')",
                                      "edgeCov(edges_adj_rankdiff,'rankdiff')",
                                      "edgeCov(edges_adj_uprank,'uprank')",
                                      "triangles",
                                      "star(c(2,3))",
                                      "cTriple",
                                      #"twoPaths",
                                      "tTriple",
                                      "mutual"),
                            gofit_terms = c("edges","degree(0:15,'out')","degree(0:15,'in')","esp(0:15)"),
                            gofit_name = "lolog_1_17")

# 18)
#Fit with preferential attachment term, fitting with star(2,3) and triangle terms
lolog_fit_1_18 <- lolog_fit(net,
                            terms = c("edges",
                                      "nodeCov('repeater_or_sweets_in','in')",
                                      "nodeCov('handicapped_in','in')",
                                      "nodeCov('rank_in','in')",
                                      "nodeCov('rank_out','out')",
                                      "edgeCov(edges_adj_rankdiff,'rankdiff')",
                                      "edgeCov(edges_adj_uprank,'uprank')",
                                      #"star(c(2,3))",
                                      #"cTriple",
                                      #"twoPaths",
                                      "mutual",
                                      "preferentialAttachment(1)"),
                            gofit_terms = c("edges","degree(0:15,'out')","degree(0:15,'in')","esp(0:15)"),
                            fit_terms = c("star(c(2,3))","triangles"),
                            gofit_name = "lolog_1_18")

#Doens't converge

# 19)
# Fit with preferential attachment term, fitting with star(2,3), including triangle term in model
lolog_fit_1_19 <- lolog_fit(net,
                            terms = c("edges",
                                      "nodeCov('repeater_or_sweets_in','in')",
                                      "nodeCov('handicapped_in','in')",
                                      "nodeCov('rank_in','in')",
                                      "nodeCov('rank_out','out')",
                                      "edgeCov(edges_adj_rankdiff,'rankdiff')",
                                      "edgeCov(edges_adj_uprank,'uprank')",
                                      "triangles",
                                      #"star(c(2,3))",
                                      #"cTriple",
                                      #"twoPaths",
                                      "mutual",
                                      "preferentialAttachment(1)"),
                            gofit_terms = c("edges","degree(0:15,'out')","degree(0:15,'in')","esp(0:15)"),
                            fit_terms = c("star(c(2,3))"),
                            gofit_name = "lolog_1_19")

#Doens't converge

# 20)
#Fit with preferential attachment term, fitting with star(2,3), including triangle term in model
#Try with higher values of k
#1  - didn't work
#10 - didn't work
#15 - didn't work
#5  - didn't work
#2  - didn't work
lolog_fit_1_20 <- lolog_fit(net,
                            terms = c("edges",
                                      "nodeCov('repeater_or_sweets_in','in')",
                                      "nodeCov('handicapped_in','in')",
                                      "nodeCov('rank_in','in')",
                                      "nodeCov('rank_out','out')",
                                      "edgeCov(edges_adj_rankdiff,'rankdiff')",
                                      "edgeCov(edges_adj_uprank,'uprank')",
                                      "triangles",
                                      #"star(c(2,3))",
                                      #"cTriple",
                                      #"twoPaths",
                                      "mutual",
                                      "preferentialAttachment(2)"),
                            gofit_terms = c("edges","degree(0:15,'out')","degree(0:15,'in')","esp(0:15)"),
                            fit_terms = c("star(c(2,3))"),
                            gofit_name = "lolog_1_20")

# 21)
#Fit with preferential attachment term, fitting with star(2,3), and triangles.
#Remove all other structural terms
#Didn't work
#Remove nodeCov(rank_in)
#Remove nodeCov(rank_in) and uprank
lolog_fit_1_21 <- lolog_fit(net,
                            terms = c("edges",
                                      #"nodeCov('repeater_or_sweets_in','in')",
                                      #"nodeCov('handicapped_in','in')",
                                      #"nodeCov('rank_in','in')",
                                      #"nodeCov('rank_out','out')",
                                      #"edgeCov(edges_adj_rankdiff,'rankdiff')",
                                      #"edgeCov(edges_adj_uprank,'uprank')",
                                      #"triangles",
                                      #"star(c(2,3))",
                                      #"cTriple",
                                      #"twoPaths",
                                      #"mutual",
                                      "preferentialAttachment(1)"),
                            gofit_terms = c("edges","degree(0:15,'out')","degree(0:15,'in')","esp(0:15)"),
                            fit_terms = c("star(c(2,3))","triangles"),
                            gofit_name = "lolog_1_21")

# 22)
#Try to use our best LOLOG fit with some ordering.
#Suspect there is some ordering process, first idea is that repeaters and sweet givers come into social
#contact with other students first therefore their friendship ties form first
#may or may not be plausible - see if it provides a better fit.

repeater_or_sweets <- (net %v% "repeater_or_sweets")
repeater_or_sweets[repeater_or_sweets == 0] <- 2

lolog_fit_1_22 <- lolog_fit(net,
                            terms = c("edges",
                                      "nodeCov('repeater_or_sweets_in','in')",
                                      "nodeCov('handicapped_in','in')",
                                      "nodeCov('rank_in','in')",
                                      "nodeCov('rank_out','out')",
                                      "edgeCov(edges_adj_rankdiff,'rankdiff')",
                                      "edgeCov(edges_adj_uprank,'uprank')",
                                      #"triangles",
                                      "star(c(2,3))",
                                      "cTriple",
                                      #"twoPaths",
                                      "tTriple",
                                      "mutual"),
                            gofit_terms = c("edges","degree(0:15,'out')","degree(0:15,'in')","esp(0:15)"),
                            gofit_name = "lolog_1_22",
                            vertex_order = repeater_or_sweets)
#similar fit to lolog without ordering

# 23)
#Try ordering based on rank, i.e. since we have a clear heirarchy in the class see if 
#accounting for ties for high nodes first gives a better model

rank <- (net %v% "rank")
#best student is ranked as 1 so so need to flip to consider worst (i.e most social) students first
rank <- 54-rank

lolog_fit_1_23 <- lolog_fit(net,
                            terms = c("edges",
                                      "nodeCov('repeater_or_sweets_in','in')",
                                      "nodeCov('handicapped_in','in')",
                                      "nodeCov('rank_in','in')",
                                      "nodeCov('rank_out','out')",
                                      "edgeCov(edges_adj_rankdiff,'rankdiff')",
                                      "edgeCov(edges_adj_uprank,'uprank')",
                                      #"triangles",
                                      "star(c(2,3))",
                                      "cTriple",
                                      #"twoPaths",
                                      "tTriple",
                                      "mutual"),
                            gofit_terms = c("edges","degree(0:15,'out')","degree(0:15,'in')","esp(0:15)"),
                            gofit_name = "lolog_1_23",
                            vertex_order = rank)

#Doesn't seem to be an appreciabley better fit than LOLOG without ordering

#Would like to look at edge specified ordering process, 
#i.e. consider neighbours in the ranking proceedure since students sit in rank order
#not suppported by LOLOG currently.

#------------ Save all the models for later ----------
setwd("..")
#list making function
make_list <- function(networks,models,type,space=ls()){
  models <- rep(models,each = length(networks))
  tmp <- unlist(mapply(networks,models,FUN = function(x,y){return(paste(type,"_",x,"_",y,sep = ""))},SIMPLIFY = FALSE))  
  return(tmp)
}

fit_ergm_list <- make_list(c(1),seq(0,4,1),"ergm_fit")
ergm_list <- make_list(c(1),seq(0,4,1),"ergm")
summary_ergm_list <- make_list(c(1),seq(0,4,1),"summary_ergm")
gofit_ergm_list <- make_list(c(1),seq(0,4,1),"gofit_ergm")

fit_lolog_list <- make_list(c(1),seq(0,23,1),"lolog_fit")
lolog_list <- make_list(c(1),seq(0,23,1),"lolog")
summary_lolog_list <- make_list(c(1),seq(0,23,1),"summary_lolog")
gofit_lolog_list <- make_list(c(1),seq(0,23,1),"gofit_lolog")

#assign ergm variables for lists
for(i in ergm_list){
  if(is.na(eval(parse(text = paste("ergm_fit",substr(i,5,nchar(i)),sep = ""))))){
    assign(paste("ergm",substr(i,5,nchar(i)),sep = ""),NA)
    assign(paste("summary_ergm",substr(i,5,nchar(i)),sep = ""),NA)
    assign(paste("gofit_ergm",substr(i,5,nchar(i)),sep = ""),NA)
  }else{
    assign(paste("ergm",substr(i,5,nchar(i)),sep = ""),eval(parse(text = paste("ergm_fit",substr(i,5,nchar(i)),"$model",sep = ""))))
    assign(paste("summary_ergm",substr(i,5,nchar(i)),sep = ""),eval(parse(text = paste("ergm_fit",substr(i,5,nchar(i)),"$summary",sep = ""))))
    assign(paste("gofit_ergm",substr(i,5,nchar(i)),sep = ""),eval(parse(text = paste("ergm_fit",substr(i,5,nchar(i)),"$gofit",sep = ""))))
  }
}

#assign lolog variables for lists
for(i in  lolog_list){
  if(is.na(eval(parse(text = paste("lolog_fit",substr(i,6,nchar(i)),sep =""))))){
    assign(paste("lolog",substr(i,6,nchar(i)),sep = ""),NA)
    assign(paste("summary_lolog",substr(i,6,nchar(i)),sep = ""),NA)
    assign(paste("gofit_lolog",substr(i,6,nchar(i)),sep = ""),NA)
  }else{
    assign(paste("lolog",substr(i,6,nchar(i)),sep = ""),eval(parse(text = paste("lolog_fit",substr(i,6,nchar(i)),"$model",sep = ""))))
    assign(paste("summary_lolog",substr(i,6,nchar(i)),sep = ""),eval(parse(text = paste("lolog_fit",substr(i,6,nchar(i)),"$summary",sep = ""))))
    assign(paste("gofit_lolog",substr(i,6,nchar(i)),sep = ""),eval(parse(text = paste("lolog_fit",substr(i,6,nchar(i)),"$gofit",sep = ""))))
  }
}

#save lists:
save(list= c(ergm_list,summary_ergm_list,gofit_ergm_list),file ="ergm_fit.RData")
save(list= c(lolog_list,summary_lolog_list,gofit_lolog_list),file ="lolog_fit.RData")
