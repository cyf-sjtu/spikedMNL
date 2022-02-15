# Start the clock!
if(enable_result_output){
  write("********************************************************************",
        file=result.file.name,append = TRUE)
}

ptm <- proc.time()

mnl.formula <- paste0("Booked_Itin~",paste(mnl.x.names,collapse = "+"),"|-1")
mnl.model<- tryCatch(mlogit(as.Formula(mnl.formula), 
                            data = data.training, shape = "long", 
                            alt.var = "Itin_No" , chid.var = "PNR",method = "bfgs"),
                     error=function(e) {
                       show(paste0("Error: ",e))
                       show("BFGS failed!")
                       elapsed_time = proc.time() - ptm
                       write(paste0("Elapsed time: ",elapsed_time[3]),
                             file=result.file.name,append = TRUE)
                       write("BFGS failed!",file = result.file.name,
                             append = TRUE)
                       ptm <<- proc.time()
                       mnl.model<-tryCatch(mlogit(as.Formula(mnl.formula), 
                                                  data = data.training, shape = "long", 
                                                  alt.var = "Itin_No" , chid.var = "PNR",
                                                  method = "nr"),
                                           error = function(e){
                                             show(paste0("Error: ",e))
                                             show("NR failed!")
                                             elapsed_time = proc.time() - ptm
                                             write(paste0("Elapsed time: ",elapsed_time[3]),
                                                   file=result.file.name,append = TRUE)
                                             write("NR failed!",file = result.file.name,
                                                   append = TRUE)
                                             return(NULL)
                                           })
                       return(mnl.model)
                     })
# Stop the clock


out <- summary(mnl.model)
show(out)

if(enable_result_output){
  if(is.null(mnl.model)){
    write(mnl.formula,file=result.file.name,append = TRUE)
    write("Numerical error occurs!!!",file = result.file.name,append = TRUE)
  }else{
    # elapsed_time = proc.time() - ptm
    write(paste0("Elapsed time: ",out$est.stat$elaps.time[3]),
          file=result.file.name,append = TRUE)
    write("\t\t\t Estimate \t Std. Error \t z-value \t Pr(>|z|)",
          file =result.file.name,append = TRUE)
    write.table(out$CoefTable,file =result.file.name,append = TRUE,
                col.names = FALSE)
    # (as.character(print(mnl.model)),file =result.file.name,append = TRUE)
    write(paste0("Log-Likelihood: ",out$logLik),
          file = result.file.name,append = TRUE)
  }
}

if(enable_xlsx_output){
  if(!is.null(mnl.model)){
    addWorksheet(wb = workbook, sheetName = paste0(mnl.model.name,"-coef"))
    writeDataTable(wb = workbook, sheet = paste0(mnl.model.name,"-coef"),
                   x = as.data.frame(out$CoefTable), rowNames = TRUE)
    
    addWorksheet(wb = workbook, sheetName = paste0(mnl.model.name,"-logLike"))
    writeDataTable(wb = workbook, sheet = paste0(mnl.model.name,"-logLike"),
                   x = as.data.frame(out$logLik))
    
    temp = data.frame(elaps.time = out$est.stat$elaps.time[3],
                      nb.iter = out$est.stat$nb.iter,
                      eps = out$est.stat$eps,
                      method = out$est.stat$method)
    addWorksheet(wb = workbook, sheetName = paste0(mnl.model.name,"-stat"))
    writeDataTable(wb = workbook, sheet = paste0(mnl.model.name,"-stat"),
                   x = temp)
    
    addWorksheet(wb = workbook, sheetName = paste0(mnl.model.name,"-gradient"))
    writeDataTable(wb = workbook, sheet = paste0(mnl.model.name,"-gradient"),
                   x = as.data.frame(out$gradient), rowNames = TRUE)
    
    addWorksheet(wb = workbook, sheetName = paste0(mnl.model.name,"-hessian"))
    writeDataTable(wb = workbook, sheet = paste0(mnl.model.name,"-hessian"),
                   x = as.data.frame(out$hessian), rowNames = TRUE)
    
    # write.xlsx(out$CoefTable,file = result.xlsx.name,
    #            sheetName = paste0(mnl.model.name,"-coef"),
    #            row.names = TRUE,
    #            append=TRUE)
    # write.xlsx(out$logLik,file = result.xlsx.name,
    #            sheetName = paste0(mnl.model.name,"-logLike"),
    #            col.names = FALSE,
    #            append=TRUE)
    # temp = data.frame(elaps.time = out$est.stat$elaps.time[3],
    #                   nb.iter = out$est.stat$nb.iter,
    #                   eps = out$est.stat$eps,
    #                   method = out$est.stat$method)
    # write.xlsx(temp,file = result.xlsx.name,
    #            sheetName = paste0(mnl.model.name,"-stat"),
    #            append=TRUE)
    # write.xlsx(out$gradient,file = result.xlsx.name,
    #            sheetName = paste0(mnl.model.name,"-gradient"),
    #            row.names = TRUE,
    #            append=TRUE)
    # write.xlsx(out$hessian,file = result.xlsx.name,
    #            sheetName = paste0(mnl.model.name,"-hessian"),
    #            row.names = TRUE,
    #            append = TRUE)
    saveWorkbook(workbook, file = result.xlsx.name, overwrite = TRUE)
  }
}
