# library(stats4)

mnl_likelihood <- function(beta = rep(0.00,length(mnl.x.names))){
  X_pro <- data.training[, c("PNR","Itin_No","Booked_Itin","pnr_weight")]
  mnl.x <- data.training[, ..mnl.x.names]
  X_pro[,v:= exp(as.matrix(mnl.x)%*% beta)]
  X_pro[,total_v:=sum(v),by=c("PNR")]
  booked_X = X_pro[Booked_Itin==1,]
  LogL = sum(booked_X$pnr_weight*(log(booked_X$v)-log(booked_X$total_v)))
  return(LogL)
}

gradient_mnl <- function(beta = rep(0.01,length(mnl.x.names))){
  # grad = numeric(length(mnl.x.names))
  X_pro <- data.training[, c("PNR","Itin_No","Booked_Itin","pnr_weight")]
  mnl.x <- data.training[, ..mnl.x.names]
  X_pro[,v:= exp(as.matrix(mnl.x)%*% beta)]
  X_pro[,total_v:=sum(v),by=c("PNR")]
  X_pro[,prob := v/total_v]
  evaluate_grad = function(column){
    #show(column)
    g = sum(X_pro$Booked_Itin * X_pro$pnr_weight * mnl.x[,..column] 
            - X_pro$pnr_weight * X_pro$prob * mnl.x[,..column])
    #show(g)
    return (g)
  }
  grad = do.call(c,lapply(mnl.x.names,
                          FUN = evaluate_grad))
  # show(grad)
  return(grad)
}

# mnl.formula <- paste0("Booked_Itin~",paste(mnl.x.names,collapse = "+"),"|-1")
# data_input = data.training
# # mnl.fit = mle(mnl_likelihood, method = "BFGS")
# ptm <- proc.time()
# mnl.fit = optim(par = rep(0.01,length(mnl.x.names)),
#                 fn = mnl_likelihood,
#                 gr = gradient_mnl,
#                 method = "BFGS")
# elapsed_time = proc.time() - ptm
# show(mnl.fit$par)
# show(elapsed_time[3])
# 
# ptm <- proc.time()
# mnl.fit = optim(par = rep(0.01,length(mnl.x.names)),
#                 fn = mnl_likelihood,
#                 method = "BFGS")
# elapsed_time = proc.time() - ptm
# show(mnl.fit$par)
# show(elapsed_time[3])

ptm <- proc.time()
mnl.model = maxLik(mnl_likelihood, grad = gradient_mnl,
                   start = rep(0.01,length(mnl.x.names)),
                   method = "BFGS")
elapsed_time = proc.time() - ptm
show(elapsed_time[3])
out = summary(mnl.model)
temp = data.frame(Attribute = mnl.x.names, out$estimate)
show(temp)

if(enable_result_output){
  write("********************************************************************",
        file=result.file.name,append = TRUE)
  if(is.null(mnl.model)){
    write(mnl.formula,file=result.file.name,append = TRUE)
    write("Numerical error occurs!!!",file = result.file.name,append = TRUE)
  }else{
    write(paste0("Elapsed time: ", elapsed_time[3]),
          file=result.file.name,append = TRUE)
    write("\t\t\t Estimate \t Std. Error \t t value \t Pr(>t)",
          file =result.file.name,append = TRUE)
    temp = data.frame(Attribute = mnl.x.names, out$estimate)
    write.table(temp,file =result.file.name,append = TRUE,row.names = FALSE,
                col.names = FALSE)
    write(paste0("Log-Likelihood: ",out$loglik),
          file = result.file.name,append = TRUE)
  }
}

if(enable_xlsx_output){
  if(!is.null(mnl.model)){
    temp = data.frame(Attribute = mnl.x.names, out$estimate)
    addWorksheet(wb = workbook, sheetName = paste0(mnl.model.name,"-coef"))
    writeDataTable(wb = workbook, sheet = paste0(mnl.model.name,"-coef"),
                   x = temp)
    
    temp = data.frame(logLik = out$loglik,
                      elaps.time = elapsed_time[3],
                      nb.iter = out$iterations,
                      eps = mnl.model$control@tol,
                      method = out$maximType)
    addWorksheet(wb = workbook, sheetName = paste0(mnl.model.name,"-stat"))
    writeDataTable(wb = workbook, sheet = paste0(mnl.model.name,"-stat"),
                   x = temp)
    
    # addWorksheet(wb = workbook, sheetName = paste0(mnl.model.name,"-gradient"))
    # writeDataTable(wb = workbook, sheet = paste0(mnl.model.name,"-gradient"),
    #                x = data.frame(Attribute = mnl.x.names,
    #                               Gradient = mnl.model$gradient))
    # 
    # temp = data.frame(mnl.x.names,mnl.model$hessian)
    # setnames(temp,c("Attribute",mnl.x.names))
    # addWorksheet(wb = workbook, sheetName = paste0(mnl.model.name,"-hessian"))
    # writeDataTable(wb = workbook, sheet = paste0(mnl.model.name,"-hessian"),
    #                x = as.data.frame(temp))
    
    saveWorkbook(workbook, file = result.xlsx.name, overwrite = TRUE)
  }
}


