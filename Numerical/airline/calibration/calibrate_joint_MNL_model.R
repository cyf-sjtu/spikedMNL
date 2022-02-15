# library(stats4)

likelihood_mnl <- function(beta = rep(0.01,length(mnl.x.names))){
  len = length(mnl.x.names)
  X_pro <- data.training[, c("PNR","Itin_No","Booked_Itin","pnr_weight",
                             "book_DtD","Dept_Dt")]
  mnl.x <- data.training[, ..mnl.x.names]
  X_pro[,v:= exp(as.matrix(mnl.x)%*% beta)]
  X_pro[,total_v:=sum(v),by=c("PNR")]
  # index for demand estimation
  X_pro[,lambda_idx := ifelse(book_DtD<= -29,-29,book_DtD)] 
  X_pro[,count_b := sum(pnr_weight*Booked_Itin), by=c("lambda_idx")]
  
  temp = unique(X_pro,by=c("Dept_Dt","lambda_idx"))
  temp[, U_b := sum(total_v/(total_v+1.0)), 
       by=c("lambda_idx")]
  booked_X = X_pro[Booked_Itin==1,]
  booked_X = merge(booked_X,temp[,c("Dept_Dt","lambda_idx","U_b")],
                   by = c("Dept_Dt","lambda_idx"))
  
  LogL = sum(booked_X$pnr_weight*(log(booked_X$v) - log(booked_X$total_v+1.0)
                             - log(booked_X$U_b)))
  temp = unique(temp, by=c("lambda_idx"))
  setorder(temp,lambda_idx)
  lambda_hat <<- temp$count_b/temp$U_b
  return(LogL)
}

gradient_mnl <- function(beta = rep(0.01,length(mnl.x.names))){
  len = length(mnl.x.names)
  X_pro <- data.training[, c("PNR","Itin_No","Booked_Itin","pnr_weight",
                             "book_DtD","Dept_Dt")]
  mnl.x <- data.training[, ..mnl.x.names]
  X_pro[,v:= exp(as.matrix(mnl.x)%*% beta)]
  X_pro[,total_v:=sum(v),by=c("PNR")]
  X_pro[,prob := v/(total_v+1.0)]
  X_pro[,lambda_idx := ifelse(book_DtD<= -29,-29, book_DtD)] # index for demand estimation
  X_pro[,count_b := sum(pnr_weight*Booked_Itin), by=c("lambda_idx")]
  temp = unique(X_pro,by=c("Dept_Dt","lambda_idx"))
  temp[, U_b := sum(total_v/(total_v+1.0)),
       by=c("lambda_idx")]
  X_pro = merge(X_pro,unique(temp[,c("lambda_idx","U_b")],by=c("lambda_idx")),
                   by = c("lambda_idx"))
  setorder(X_pro,"PNR","Itin_No")
  PNR_selected = temp$PNR
  X_pro_b = X_pro[PNR %in% PNR_selected,]
  setorder(X_pro_b, PNR, Itin_No)
  # mnl.x.b = data.training[PNR %in% temp$PNR, ..mnl.x.names]
  mnl.x.b = data.training[PNR %in% PNR_selected, c("PNR", "Itin_No", ..mnl.x.names)]
  setorder(mnl.x.b,PNR,Itin_No)
  evaluate_grad = function(column){
    g = sum(X_pro$Booked_Itin * X_pro$pnr_weight * mnl.x[,..column] 
            - X_pro$pnr_weight * X_pro$prob * mnl.x[,..column])
    X_pro_b$xxx = mnl.x.b[,..column]
    X_pro_b[, diff_hbk := prob/(total_v+1.0)*xxx]
    X_pro_b[, diff_hb := sum(diff_hbk), by=c("PNR")]
    X_pro_bb = X_pro_b[Booked_Itin == 1,]
    X_pro_bb[, diff_b := sum(diff_hb), by=c("lambda_idx")]
    X_pro_bbb = unique(X_pro_bb, by=c("lambda_idx"))
    g = g - sum(X_pro_bbb$count_b*X_pro_bbb$diff_b/X_pro_bbb$U_b)
    return (g)
  }
  grad = do.call(c,lapply(mnl.x.names, FUN = evaluate_grad))
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

# data.training = data.training[book_DtD>= -29, ]

setorder(data.training,"PNR","Itin_No")

lambda_hat = numeric(30)

ptm <- proc.time()
mnl.model = maxLik(likelihood_mnl, grad = gradient_mnl,
                   start = rep(0.0,length(mnl.x.names)),
                   method = "BFGS")
elapsed_time = proc.time() - ptm
show(elapsed_time[3])
out = summary(mnl.model)
temp = data.frame(Attribute = c(mnl.x.names), out$estimate)
show(temp)

plot(-29:0, lambda_hat, "l")

data.training[,lambda_idx := ifelse(book_DtD<= -29,-29, book_DtD)]
booked = data.training[Booked_Itin==1,]
demand = count(booked, vars = c("lambda_idx"),wt_var = "pnr_weight")
setorder(demand,lambda_idx)
num_sample_path = uniqueN(booked$Dept_Dt)

lines(-29:0, demand$freq/num_sample_path, "l")
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
    
    addWorksheet(wb = workbook, sheetName = paste0(mnl.model.name,"-gradient"))
    writeDataTable(wb = workbook, sheet = paste0(mnl.model.name,"-gradient"),
                   x = data.frame(Attribute = mnl.x.names,
                                  Gradient = mnl.model$gradient))
    
    temp = data.frame(mnl.x.names,mnl.model$hessian)
    setnames(temp,c("Attribute",mnl.x.names))
    addWorksheet(wb = workbook, sheetName = paste0(mnl.model.name,"-hessian"))
    writeDataTable(wb = workbook, sheet = paste0(mnl.model.name,"-hessian"),
                   x = as.data.frame(temp))
    
    saveWorkbook(workbook, file = result.xlsx.name, overwrite = TRUE)
  }
}


