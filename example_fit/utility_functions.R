# -------- Utility Functions for fitting and GOF procedure ------- #

# We define some utilities that allowed for standardisation accross the many examples we considered.

ergm_fit <- function(net,
                     terms,
                     gofit_terms,
                     gofit_name,
                     offset.coef = NULL,
                     control.ergm = NULL){
  
  formula <- as.formula(paste(c("net ~",paste(terms,collapse = "+")),collapse = " "))
  gofit_formula <- as.formula(paste("net~",paste(gofit_terms,collapse = "+"),collapse = " "))
  
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
  
  formula <- paste("net ~",paste(terms,collapse = "+"),collapse = " ")
  gofit_formula <- paste("net ~",paste(gofit_terms,collapse = "+"),collapse = " ")
  
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
    auxFormula <- as.formula(paste("net ~",paste(fit_terms,collapse = "+"),collapse=" "))
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
  
  if(length(model) <= 1 & is.na(model)){
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
