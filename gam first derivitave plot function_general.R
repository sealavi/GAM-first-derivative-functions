species_interact_deriv <- function (model, term, main, eps, response = NULL, confidence=95, output){
  ###model must be a brms model object
  ###term is a character string of the smooth term, same syntax as used in the model
  ###main is a character string of the predictor variable, must not be wrapped in a smooth function
  ###eps is the amount to offset the original data, to be differenced from original to calculate slope
  ###response is an optional character string indicating the response variable to use, only relevant in the multivariate case
  ###confidence is the confidence level used to calculate the posterior intervals
  ###The desired name of the resulting ggplot object 
  
  require(dplyr)
  require(ggplot2)
  require(brms)
  
  Response=response
  if(is.null(Response)){
    Response=model$formula$resp
  }
    
  upper=(50+(confidence/2))/100
  lower=(50-(confidence/2))/100
  
  newdat=model$data
  newdat[,which(names(newdat)==main)]=newdat[,which(names(newdat)==main)]+eps
  dir=posterior_smooths(model, smooth = term, resp=response)
  dir2=posterior_smooths(model, smooth = term, resp=response, newdata = newdat)
  
  dir_model=(dir2-dir)/eps
  
  mean_der <- apply(dir_model,MARGIN = 2,FUN = mean)
  lower_der <- apply(dir_model,MARGIN = 2,FUN = quantile, prob = lower)
  upper_der <- apply(dir_model,MARGIN = 2,FUN = quantile, prob = upper)
  
  der_data=data.frame(mean_der) %>%
    cbind(lower_der) %>%
    cbind(upper_der) %>%
    cbind(model$data[,which(names(model$data)==main)]) %>%
    cbind(model$data[,which(names(model$data)=="Species")])
  colnames(der_data)=c("mean","lower","upper","main","Species")
  
  der_data$Significance=NA
  der_data$Significance[which(sign(der_data$lower)<0&sign(der_data$upper)<0)]=-1
  der_data$Significance[which(sign(der_data$lower)>0&sign(der_data$upper)>0)]=1
  der_data$Significance[which(sign(der_data$lower)!=sign(der_data$upper))]=0
  
  der_data$siglab <- with(rle(der_data$Significance), rep(cumsum(lengths >= 1),lengths))
  
  if(length(which(der_data$Significance!=0))==0){
    model_plot=plot(conditional_effects(model,spaghetti=FALSE),plot=FALSE)
    if(is.null(response)){
      index=which(names(model_plot)==paste(main,":Species",sep=""))
    }else{
      index=which(names(model_plot)==paste(response,".",response,"_",main,":Species",sep=""))
      
    }
    model_est <- as.data.frame(model_plot[[index]][[1]])
    index2=which(names(model_est)==main)
    colnames(model_est)[index2]="Main"
    
    model_plot=ggplot(model_est,aes(Main,estimate__)) + 
      geom_line(aes(color=I("black")),size=1)+
      geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = .1)+
      ylab(Response)+xlab(main)+
      theme_classic()+ guides(color="none")+facet_wrap(~Species)
    assign(output,model_plot, envir = parent.frame())
    return(model_plot)
    
  } else{
    der_data_SIG=der_data[which(der_data$Significance!=0),]
    
    sigranges=tapply(der_data_SIG$main,as.factor(der_data_SIG$siglab),range)

    model_plot=plot(conditional_effects(model,spaghetti=FALSE),plot=FALSE)
    if(is.null(response)){
      index=which(names(model_plot)==paste(main,":Species",sep=""))
    }else{
      index=which(names(model_plot)==paste(response,".",response,"_",main,":Species",sep=""))
      
    }
    model_est <- as.data.frame(model_plot[[index]][[1]])
    model_est$Sig=NA
    index2=which(names(model_est)==main)
    colnames(model_est)[index2]="Main"
    
    for(i in 1:nrow(model_est)){
      for(j in 1:length(sigranges)){
        if(model_est$Main[i]>=sigranges[[j]][1] & model_est$Main[i]<sigranges[[j]][2]){
          model_est$Sig[i]=1
        }
      }
    }
    

    
    model_est$Sig[-which(model_est$Sig==1)]=0
    if(length(which(model_est$Sig==1))==0){
      model_est$Sig=0
    }
    model_plot=ggplot(model_est,aes(Main,estimate__)) + 
      geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = .1)+
      geom_line(aes(color=(Sig)),size=1)+
      scale_color_gradient2(low="black", mid="black",high="cyan" )+
      ylab(Response)+xlab(main)+
      theme_classic()+ guides(color="none")+facet_wrap(~Species)
    assign(output,model_plot, envir = parent.frame())
    return(model_plot)
    
  }

}


deriv_plot <- function (model, term, main, eps, response = NULL, spaghetti=FALSE, rug = TRUE, confidence=95,output){
  require(dplyr)
  require(ggplot2)
  require(brms)
  
  ###model must be a brms model object
  ###term is a character string of the smooth term, same syntax as used in the model
  ###main is a character string of the predictor variable, must not be wrapped in a smooth function
  ###eps is the amount to offset the original data, to be differenced from original to calculate slope
  ###response is an optional character string indicating the response variable to use, only relevant in the multivariate case
  ###confidence is the confidence level used to calculate the posterior intervals
  ###The desired name of the resulting ggplot object 
  Response=response
  if(is.null(Response)){
    Response=model$formula$resp
  }
  
  if(length(names(model$data))>6){
    model$data=model$data[,c(1:6)]
  }
  
  upper=(50+(confidence/2))/100
  lower=(50-(confidence/2))/100
  
  newdat=model$data
  newdat[,which(names(newdat)==main)]=newdat[,which(names(newdat)==main)]+eps
  dir=posterior_smooths(model, smooth = term, resp=response)
  dir2=posterior_smooths(model, smooth = term, resp=response, newdata = newdat)
  
  dir_model=(dir2-dir)/eps
  
  mean_der <- apply(dir_model,MARGIN = 2,FUN = mean)
  lower_der <- apply(dir_model,MARGIN = 2,FUN = quantile, prob = lower)
  upper_der <- apply(dir_model,MARGIN = 2,FUN = quantile, prob = upper)
  
  der_data=data.frame(mean_der) %>%
    cbind(lower_der) %>%
    cbind(upper_der) %>%
    cbind(model$data[,which(names(model$data)==main)])
  colnames(der_data)=c("mean","lower","upper","main")
  
  der_data$Significance=NA
  der_data$Significance[which(sign(der_data$lower)<0&sign(der_data$upper)<0)]="Significant"
  der_data$Significance[which(sign(der_data$lower)>0&sign(der_data$upper)>0)]="Significant"
  der_data$Significance[which(sign(der_data$lower)!=sign(der_data$upper))]="Not Significant"
  #sigranges=tapply(der_data$main,as.factor(der_data$Significance),range)
  
  der_data$Significance=NA
  der_data$Significance[which(sign(der_data$lower)<0&sign(der_data$upper)<0)]=-1
  der_data$Significance[which(sign(der_data$lower)>0&sign(der_data$upper)>0)]=1
  der_data$Significance[which(sign(der_data$lower)!=sign(der_data$upper))]=0
  der_data=der_data[order(der_data$main),]
  der_data$siglab <- with(rle(der_data$Significance), rep(cumsum(lengths >= 1),lengths))
  
  if(length(which(der_data$Significance!=0))==0){
    model_plot=plot(conditional_effects(model,spaghetti=spaghetti),rug = rug,errorbar_args = list(alpha=0.1),plot=FALSE)
    if(is.null(response)){
      index=which(names(model_plot)==paste(main,sep=""))
      }else{
        index=which(names(model_plot)==paste(response,".",response,"_",main,sep=""))
    }    
    model_est <- as.data.frame(model_plot[[index]][[1]])
    model_plot=plot(conditional_effects(model,spaghetti=spaghetti),rug = rug,errorbar_args = list(alpha=0.1),plot=FALSE)[[index]]
    
    index2=which(names(model_est)==main)
    colnames(model_est)[index2]="Main"
    
    model_plot2=model_plot+
      geom_line(data=model_est,aes(Main,estimate__,color=I("black")),size=1)+
      ylab(Response)+xlab(main)+
      theme_classic()+ guides(color="none")
    assign(output,model_plot2, envir = parent.frame())
    return(model_plot2)
    
  } else{
    der_data_SIG=der_data[which(der_data$Significance!=0),]
    
    sigranges=tapply(der_data_SIG$main,as.factor(der_data_SIG$siglab),range, na.rm=T)
    
    model_plot=plot(conditional_effects(model,spaghetti=spaghetti),rug = rug,errorbar_args = list(alpha=0.1),plot=FALSE)
    if(is.null(response)){
      index=which(names(model_plot)==paste(main,sep=""))
      }else{
        index=which(names(model_plot)==paste(response,".",response,"_",main,sep=""))
    }    
    model_plot=plot(conditional_effects(model,spaghetti=spaghetti),rug = rug, errorbar_args = list(alpha=0.1),plot=FALSE)[[index]]
    
    model_est <- as.data.frame(model_plot[[1]])
    model_est$Sig=NA
    model_est$Sig2=NA
    model_est$Sig2[which(model_est$Sig==0)]=.8
    model_est$Sig2[which(model_est$Sig==1)]=1.5
    index2=which(names(model_est)==main)
    colnames(model_est)[index2]="Main"
    

    for(i in 1:nrow(model_est)){
      for(j in 1:length(sigranges)){
        if(model_est$Main[i]>=sigranges[[j]][1] & model_est$Main[i]<sigranges[[j]][2]){
          model_est$Sig[i]=1
        }
      }
      
      
    }
    model_est$Sig[-which(model_est$Sig==1)]=0
    if(length(which(model_est$Sig==1))==0){
      model_est$Sig=0
    }
    model_plot2=model_plot+ 
      geom_line(data=model_est,aes(Main,estimate__,color=(Sig)),size=1)+
      scale_color_gradient2(low="black", mid="black",high="cyan" )+
      ylab(Response)+xlab(main)+
      theme_classic()+ guides(color="none")

    assign(output,model_plot2, envir = parent.frame())
    output2=gsub("_plot", "", output)
    output2=paste("VOI",output2,sep="_")
    if(length(which(model_est$Sig==1))>0){
      VOIdat=model_est[which(model_est$Sig==1),]
      assign(output2,VOIdat, envir = parent.frame())
    }
    return(model_plot2)
    
  }
  
}



deriv_plot2 <- function (model, term, main, eps, response = NULL, confidence=95, output){
  require(dplyr)
  require(ggplot2)
  require(brms)
  
  ###model must be a brms model object
  ###term is a character string of the smooth term, same syntax as used in the model
  ###main is a character string of the predictor variable, must not be wrapped in a smooth function
  ###eps is the amount to offset the original data, to be differenced from original to calculate slope
  ###response is an optional character string indicating the response variable to use, only relevant in the multivariate case
  ###confidence is the confidence level used to calculate the posterior intervals
  ###The desired name of the resulting ggplot object 
  
  Response=response
  if(is.null(Response)){
    Response=model$formula$resp
  }
  if(length(names(model$data))>6){
    model$data=model$data[,c(1:6)]
  }
  
  upper=(50+(confidence/2))/100
  lower=(50-(confidence/2))/100
  
  newdat=model$data
  newdat[,which(names(newdat)==main)]=newdat[,which(names(newdat)==main)]+eps
  dir=posterior_smooths(model, smooth = term, resp=response)
  dir2=posterior_smooths(model, smooth = term, resp=response, newdata = newdat)
  
  dir_model=(dir2-dir)/eps
  
  mean_der <- apply(dir_model,MARGIN = 2,FUN = mean)
  lower_der <- apply(dir_model,MARGIN = 2,FUN = quantile, prob = lower)
  upper_der <- apply(dir_model,MARGIN = 2,FUN = quantile, prob = upper)
  
  der_data=data.frame(mean_der) %>%
    cbind(lower_der) %>%
    cbind(upper_der) %>%
    cbind(model$data[,which(names(model$data)==main)])
  colnames(der_data)=c("mean","lower","upper","main")
  
  der_data$Significance=NA
  der_data$Significance[which(sign(der_data$lower)<0&sign(der_data$upper)<0)]="Significant"
  der_data$Significance[which(sign(der_data$lower)>0&sign(der_data$upper)>0)]="Significant"
  der_data$Significance[which(sign(der_data$lower)!=sign(der_data$upper))]="Not Significant"
  sigranges=tapply(der_data$main,as.factor(der_data$Significance),range)
  der_data$Significance=as.factor(der_data$Significance)
  derivplot=ggplot(der_data,aes(main,mean,group=1)) + 
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2)+
    geom_line(aes(color=Significance)) + 
    geom_hline(yintercept = 0)+
    scale_color_manual(values = c("Significant" = "cyan",
                                  "Not Significant"="black")) +
    theme_classic()
  assign(output,derivplot, envir = parent.frame())
  
  
  
  ## Adapted from original deriv_plot function 

## function allowing for estimation of derivative at different confidence levels set by the user
## adapted to also include how much of derivative is on one side of 0
deriv_plot_zprob <- function (model, dimensions = 1, by = FALSE, term, main, eps, response = NULL, spaghetti=FALSE, rug = TRUE, confidence=95,output, meanmain, sdmain){
  require(dplyr)
  require(ggplot2)
  require(plotly)
  require(brms)
  
  ###model must be a brms model object
  ###dimensions should be the number of variables in your spline
  ###term is a character string of the smooth term, same syntax as used in the model
  ###main is a character string (or vector of characters equal to dimensions) of the predictor variable, must not be wrapped in a smooth function
  ###eps is the amount to offset the original data (or a vector of offsets equal to dimensions), to be differenced from original to calculate slope
  ###response is an optional character string indicating the response variable to use, only relevant in the multivariate case
  ###confidence is the confidence level used to calculate the posterior intervals
  ###The desired name of the resulting ggplot object 
  ## meanmain is a vector of the means of your main predictor variable(s)
  ## sdmain is a vector of the means of your main predictor variable(s)
  Response=response
  if(is.null(Response)){
    Response=model$formula$resp
  }
  
  if(length(names(model$data))>6){
    model$data=model$data[,c(1:6)]
  }
  
  upper=(50+(confidence/2))/100
  lower=(50-(confidence/2))/100
  
  newdat=model$data
  newdat_b=model$data
  newdat_c=model$data
  newdat_d=model$data
  
  ##for 2D smooth finite difference aprox something like this
  ##fxy(x,y)~(f(x+eps_h,y+eps_k)-f(x+eps_h,y-eps_k)-f(x-eps_h,y+eps_k)+f(x-eps_h,y-eps_k))/(4*eps_h*eps*k)
  if(dimensions > 1) {
    if(length(eps)==1){
      eps[2]=eps[1]
    }
    for(i in 1:dimensions) {
      newdat[,which(names(newdat)==main[i])]=newdat[,which(names(newdat)==main[i])]+eps[i] # h + K
      newdat_b[,which(names(newdat_b)==main[i])]=newdat_b[,which(names(newdat_b)==main[i])]-eps[i] # -h - K
      
    }
    
    
    #h - k
    newdat_c[,which(names(newdat_c)==main[1])]=newdat_c[,which(names(newdat_c)==main[1])]+eps[1] 
    newdat_c[,which(names(newdat_c)==main[2])]=newdat_c[,which(names(newdat_c)==main[2])]-eps[2] 
    
    #-h + k
    newdat_d[,which(names(newdat_d)==main[1])]=newdat_d[,which(names(newdat_d)==main[1])]-eps[1]  
    newdat_d[,which(names(newdat_d)==main[2])]=newdat_d[,which(names(newdat_d)==main[2])]+eps[2] 
    
    
  } else{
    newdat[,which(names(newdat)==main)]=newdat[,which(names(newdat)==main)]+eps
  } 
  
  if(dimensions > 1){
    #dir=posterior_smooths(model, smooth = term, resp=response)
    dir2=posterior_smooths(model, smooth = term, resp=response, newdata = newdat)
    dir2_b=posterior_smooths(model, smooth = term, resp=response, newdata = newdat_b)
    dir2_c=posterior_smooths(model, smooth = term, resp=response, newdata = newdat_c)
    dir2_d=posterior_smooths(model, smooth = term, resp=response, newdata = newdat_d)
    
    dir_model=(dir2-dir2_c-dir2_d+dir2_b)/(4*prod(eps))
    
    mean_der <- apply(dir_model,MARGIN = 2,FUN = mean)
    lower_der <- apply(dir_model,MARGIN = 2,FUN = quantile, prob = lower)
    upper_der <- apply(dir_model,MARGIN = 2,FUN = quantile, prob = upper)
    
    probrange <- function(x) {
      sum(x>0)/length(x)
    }
    probrange2 <- function(x) {
      sum(x < 0)/length(x)
    }
    sum(dir_model[,1] > 0 )/length(dir_model[,1])
    prob_above <- apply(dir_model, MARGIN = 2, probrange)
    prob_below <- apply(dir_model, MARGIN = 2, probrange2)
    
    der_data = cbind(mean_der, lower_der, upper_der, prob_above, prob_below)
    
    for(i in 1:length(main)) {
      der_data = cbind(der_data, model$data[,which(names(model$data)==main[i])])
    }
    der_data <- as.data.frame(der_data)
    colnames(der_data)=c("mean","lower","upper", "prob_above", "prob_below", main[1:length(main)])
    
    if(is.null(by)==TRUE) {
      interpdat <- with(der_data, akima::interp(x = der_data[,6], y = der_data[,7], z = mean, duplicate = "mean"))
      interpdat2 <- reshape2::melt(interpdat$z, na.rm = TRUE)
      names(interpdat2) <- c("x", "y", "dir")
      interpdat2$main1 <- interpdat$x[interpdat2$x]
      interpdat2$main2 <- interpdat$y[interpdat2$y]
      interpdat_low <- with(der_data, akima::interp(x = der_data[,6], y = der_data[,7], z = lower, duplicate = "mean"))
      interpdat2_low <- reshape2::melt(interpdat_low$z, na.rm = TRUE)
      names(interpdat2_low) <- c("x", "y", "dir")
      interpdat2_low$main1 <- interpdat_low$x[interpdat2_low$x]
      interpdat2_low$main2 <- interpdat_low$y[interpdat2_low$y]
      interpdat_high <- with(der_data, akima::interp(x = der_data[,6], y = der_data[,7], z = upper, duplicate = "mean"))
      interpdat2_high <- reshape2::melt(interpdat_high$z, na.rm = TRUE)
      names(interpdat2_high) <- c("x", "y", "dir")
      interpdat2_high$main1 <- interpdat_high$x[interpdat2_high$x]
      interpdat2_high$main2 <- interpdat_high$y[interpdat2_high$y]
      interpdat2$upper=interpdat2_high$dir
      interpdat2$lower=interpdat2_low$dir
      interpdat2$threshold=0
    } else {
      # add by column to der_data
      # for now only set up for by variable with TWO LEVELS and in quite explicit/roundabout way
      der_data = cbind(der_data, model$data[,which(names(model$data)==by)])
      colnames(der_data)=c("mean","lower","upper", "prob_above", "prob_below", main[1:length(main)], by)
      ## transform back to real scale for z-transformed variables
      der_data[,which(names(der_data)==main[1])] <- der_data[,which(names(der_data)==main[1])] * sdmain[1] + meanmain[1]
      der_data[,which(names(der_data)==main[2])] <- der_data[,which(names(der_data)==main[2])] * sdmain[2] + meanmain[2]
      
      # factor level 1
      der_data_by1 <- der_data[which(der_data[,8] == levels(der_data[,8])[1]),]
      interpdat_a <- with(der_data_by1, akima::interp(x = der_data_by1[,6], y = der_data_by1[,7], z = mean, duplicate = "mean"))
      interpdat_a2 <- reshape2::melt(interpdat_a$z, na.rm = TRUE)
      names(interpdat_a2) <- c("x", "y", "dir")
      interpdat_a2$main1 <- interpdat_a$x[interpdat_a2$x]
      interpdat_a2$main2 <- interpdat_a$y[interpdat_a2$y]
      interpdat_a_low <- with(der_data_by1, akima::interp(x = der_data_by1[,6], y = der_data_by1[,7], z = lower, duplicate = "mean"))
      interpdat_a2_low <- reshape2::melt(interpdat_a_low$z, na.rm = TRUE)
      names(interpdat_a2_low) <- c("x", "y", "dir")
      interpdat_a2_low$main1 <- interpdat_a_low$x[interpdat_a2_low$x]
      interpdat_a2_low$main2 <- interpdat_a_low$y[interpdat_a2_low$y]
      interpdat_a_high <- with(der_data_by1, akima::interp(x = der_data_by1[,6], y = der_data_by1[,7], z = upper, duplicate = "mean"))
      interpdat_a2_high <- reshape2::melt(interpdat_a_high$z, na.rm = TRUE)
      names(interpdat_a2_high) <- c("x", "y", "dir")
      interpdat_a2_high$main1 <- interpdat_a_high$x[interpdat_a2_high$x]
      interpdat_a2_high$main2 <- interpdat_a_high$y[interpdat_a2_high$y]
      interpdat_a_probabove <- with(der_data_by1, akima::interp(x = der_data_by1[,6], y = der_data_by1[,7], z = prob_above, duplicate = "mean"))
      interpdat_a2_probabove <- reshape2::melt(interpdat_a_probabove$z, na.rm = TRUE)
      names(interpdat_a2_probabove) <- c("x", "y", "dir")
      interpdat_a2_probabove$main1 <- interpdat_a_probabove$x[interpdat_a2_probabove$x]
      interpdat_a2_probabove$main2 <- interpdat_a_probabove$y[interpdat_a2_probabove$y]
      interpdat_a_probbelow <- with(der_data_by1, akima::interp(x = der_data_by1[,6], y = der_data_by1[,7], z = prob_below, duplicate = "mean"))
      interpdat_a2_probbelow <- reshape2::melt(interpdat_a_probbelow$z, na.rm = TRUE)
      names(interpdat_a2_probbelow) <- c("x", "y", "dir")
      interpdat_a2_probbelow$main1 <- interpdat_a_probbelow$x[interpdat_a2_probbelow$x]
      interpdat_a2_probbelow$main2 <- interpdat_a_probbelow$y[interpdat_a2_probbelow$y]
      interpdat_a2$upper=interpdat_a2_high$dir
      interpdat_a2$lower=interpdat_a2_low$dir
      interpdat_a2$probabove = interpdat_a2_probabove$dir
      interpdat_a2$probbelow = interpdat_a2_probbelow$dir
      interpdat_a2$threshold=0
      assign(paste(output, "1p", sep = "_"),interpdat_a2, envir = parent.frame())
      
      
      # factor level 2
      der_data_by2 <- der_data[which(der_data[,8] == levels(der_data[,8])[2]),]
      interpdat_b <- with(der_data_by2, akima::interp(x = der_data_by2[,6], y = der_data_by2[,7], z = mean, duplicate = "mean"))
      interpdat_b2 <- reshape2::melt(interpdat_b$z, na.rm = TRUE)
      names(interpdat_b2) <- c("x", "y", "dir")
      interpdat_b2$main1 <- interpdat_b$x[interpdat_b2$x]
      interpdat_b2$main2 <- interpdat_b$y[interpdat_b2$y]
      interpdat_b_low <- with(der_data_by2, akima::interp(x = der_data_by2[,6], y = der_data_by2[,7], z = lower, duplicate = "mean"))
      interpdat_b2_low <- reshape2::melt(interpdat_b_low$z, na.rm = TRUE)
      names(interpdat_b2_low) <- c("x", "y", "dir")
      interpdat_b2_low$main1 <- interpdat_b_low$x[interpdat_b2_low$x]
      interpdat_b2_low$main2 <- interpdat_b_low$y[interpdat_b2_low$y]
      interpdat_b_high <- with(der_data_by2, akima::interp(x = der_data_by2[,6], y = der_data_by2[,7], z = upper, duplicate = "mean"))
      interpdat_b2_high <- reshape2::melt(interpdat_b_high$z, na.rm = TRUE)
      names(interpdat_b2_high) <- c("x", "y", "dir")
      interpdat_b2_high$main1 <- interpdat_b_high$x[interpdat_b2_high$x]
      interpdat_b2_high$main2 <- interpdat_b_high$y[interpdat_b2_high$y]
      interpdat_b_probabove <- with(der_data_by2, akima::interp(x = der_data_by2[,6], y = der_data_by2[,7], z = prob_above, duplicate = "mean"))
      interpdat_b2_probabove <- reshape2::melt(interpdat_b_probabove$z, na.rm = TRUE)
      names(interpdat_b2_probabove) <- c("x", "y", "dir")
      interpdat_b2_probabove$main1 <- interpdat_b_probabove$x[interpdat_b2_probabove$x]
      interpdat_b2_probabove$main2 <- interpdat_b_probabove$y[interpdat_b2_probabove$y]
      interpdat_b_probbelow <- with(der_data_by2, akima::interp(x = der_data_by2[,6], y = der_data_by2[,7], z = prob_below, duplicate = "mean"))
      interpdat_b2_probbelow <- reshape2::melt(interpdat_b_probbelow$z, na.rm = TRUE)
      names(interpdat_b2_probbelow) <- c("x", "y", "dir")
      interpdat_b2_probbelow$main1 <- interpdat_b_probbelow$x[interpdat_b2_probbelow$x]
      interpdat_b2_probbelow$main2 <- interpdat_b_probbelow$y[interpdat_b2_probbelow$y]
      interpdat_b2$upper=interpdat_b2_high$dir
      interpdat_b2$lower=interpdat_b2_low$dir
      interpdat_b2$probabove=interpdat_b2_probabove$dir
      interpdat_b2$probbelow = interpdat_b2_probbelow$dir
      interpdat_b2$threshold=0
      assign(paste(output, "2p", sep = "_"),interpdat_b2, envir = parent.frame())
      
    }
    
    if(is.null(by)==TRUE){
      axx <- list(
        title = names(model$data)[3]
      )
      
      axy <- list(
        title = names(model$data)[4]
      )
      
      
      p <- plot_ly(interpdat2, x=~main1, y=~main2, 
                   z=~dir, intensity = ~dir,type="mesh3d") %>% 
        add_mesh(x=~main1, y=~main2, 
                 z=~upper, intensity = ~upper, opacity=0.30) %>%
        add_mesh(x=~main1, y=~main2, 
                 z=~lower, intensity = ~lower, opacity=0.30)  %>%
        add_mesh(x=~main1, y=~main2, 
                 z=~threshold, intensity = ~threshold, colorscale='Hot' )
      p=p%>% hide_colorbar()
      p <- p %>% layout(title = "Derivative",
                        scene = list(xaxis=axx, yaxis=axy,
                                     aspectmode='cube'))
      assign(output,p, envir = parent.frame())
      return(p)
    } else{ 
      axx <- list(
        title = names(model$data)[3]
      )
      
      axy <- list(
        title = names(model$data)[4]
      )
      
      p1 <- plot_ly(interpdat_a2, x=~main1, y=~main2, 
                    z=~dir, intensity = ~dir, scene= 'scene1', type="mesh3d") %>% 
        add_mesh(x=~main1, y=~main2, 
                 z=~upper, intensity = ~upper, opacity=0.30) %>%
        add_mesh(x=~main1, y=~main2, 
                 z=~lower, intensity = ~lower, opacity=0.30)  %>%
        add_mesh(x=~main1, y=~main2, 
                 z=~threshold, intensity = ~threshold, colorscale='Hot' )
      p1=p1%>% hide_colorbar()
      p1 <- p1 %>% layout(annotations = list(x = 0.2 , y = 0.95, text = paste(by, levels(der_data[,8])[1], sep = ": "),
                                             showarrow = F, xref='paper', yref='paper', font = list(size = 15)), showlegend = FALSE) 
      
      p2 <- plot_ly(interpdat_b2, x=~main1, y=~main2, 
                    z=~dir, intensity = ~dir, scene= 'scene2', type="mesh3d") %>% 
        add_mesh(x=~main1, y=~main2, 
                 z=~upper, intensity = ~upper, opacity=0.30) %>%
        add_mesh(x=~main1, y=~main2, 
                 z=~lower, intensity = ~lower, opacity=0.30)  %>%
        add_mesh(x=~main1, y=~main2, 
                 z=~threshold, intensity = ~threshold, colorscale='Hot' )
      p2=p2%>% hide_colorbar()
      p2 <- p2 %>% layout(annotations = list(x = 0.2 , y = 0.95, text = paste(by, levels(der_data[,8])[2], sep = ": "),
                                             showarrow = F, xref='paper', yref='paper', font = list(size = 15)), showlegend = FALSE) 
      
      pp <- subplot(p1, p2)
      pp <- pp %>% layout(title = paste("Derivative at confidence", confidence, sep = " "),
                          scene = list(xaxis=axx, yaxis=axy,
                                       aspectmode='cube'),
                          scene2 = list(xaxis=axx, yaxis=axy,
                                        aspectmode='cube'))
      assign(output,pp, envir = parent.frame())
      return(pp) 
    }
    
    
    
  } else{
    newdat=model$data
    newdat[,which(names(newdat)==main)]=newdat[,which(names(newdat)==main)]+eps
    dir=posterior_smooths(model, smooth = term, resp=response)
    dir2=posterior_smooths(model, smooth = term, resp=response, newdata = newdat)
    
    dir_model=(dir2-dir)/eps
    
    mean_der <- apply(dir_model,MARGIN = 2,FUN = mean)
    lower_der <- apply(dir_model,MARGIN = 2,FUN = quantile, prob = lower)
    upper_der <- apply(dir_model,MARGIN = 2,FUN = quantile, prob = upper)
    
    der_data=data.frame(mean_der) %>%
      cbind(lower_der) %>%
      cbind(upper_der) %>%
      cbind(model$data[,which(names(model$data)==main)])
    colnames(der_data)=c("mean","lower","upper","main")
    
    
    der_data$Significance=NA
    der_data$Significance[which(sign(der_data$lower)<0&sign(der_data$upper)<0)]="Significant"
    der_data$Significance[which(sign(der_data$lower)>0&sign(der_data$upper)>0)]="Significant"
    der_data$Significance[which(sign(der_data$lower)!=sign(der_data$upper))]="Not Significant"
    #sigranges=tapply(der_data$main,as.factor(der_data$Significance),range)
    
    der_data$Significance=NA
    der_data$Significance[which(sign(der_data$lower)<0&sign(der_data$upper)<0)]=-1
    der_data$Significance[which(sign(der_data$lower)>0&sign(der_data$upper)>0)]=1
    der_data$Significance[which(sign(der_data$lower)!=sign(der_data$upper))]=0
    #der_data=der_data[with(der_data, order(der_data[,4], der_data[,5])),]
    der_data$siglab <- with(rle(der_data$Significance), rep(cumsum(lengths >= 1),lengths))
    
    
    if(length(which(der_data$Significance!=0))==0){
      model_plot=plot(conditional_effects(model,spaghetti=spaghetti),rug = rug,errorbar_args = list(alpha=0.1),plot=FALSE)
      if(is.null(response)){
        index=which(names(model_plot)==paste(main,sep=""))
      }else{
        index=which(names(model_plot)==paste(response,".",response,"_",main,sep=""))
      }    
      model_est <- as.data.frame(model_plot[[index]][[1]])
      model_plot=plot(conditional_effects(model,spaghetti=spaghetti),rug = rug,errorbar_args = list(alpha=0.1),plot=FALSE)[[index]]
      
      index2=which(names(model_est)==main)
      colnames(model_est)[index2]="Main"
      
      model_plot2=model_plot+
        geom_line(data=model_est,aes(Main,estimate__,color=I("black")),size=1)+
        ylab(Response)+xlab(main)+
        theme_classic()+ guides(color="none")
      assign(output,model_plot2, envir = parent.frame())
      return(model_plot2)
      
    } else{
      der_data_SIG=der_data[which(der_data$Significance!=0),]
      
      sigranges=tapply(der_data_SIG$main,as.factor(der_data_SIG$siglab),range, na.rm=T)
      
      model_plot=plot(conditional_effects(model,spaghetti=spaghetti),rug = rug,errorbar_args = list(alpha=0.1),plot=FALSE)
      if(is.null(response)){
        index=which(names(model_plot)==paste(main,sep=""))
      }else{
        index=which(names(model_plot)==paste(response,".",response,"_",main,sep=""))
      }    
      model_plot=plot(conditional_effects(model,spaghetti=spaghetti),rug = rug, errorbar_args = list(alpha=0.1),plot=FALSE)[[index]]
      
      model_est <- as.data.frame(model_plot[[1]])
      model_est$Sig=NA
      model_est$Sig2=NA
      model_est$Sig2[which(model_est$Sig==0)]=.8
      model_est$Sig2[which(model_est$Sig==1)]=1.5
      index2=which(names(model_est)==main)
      colnames(model_est)[index2]="Main"
      
      
      for(i in 1:nrow(model_est)){
        for(j in 1:length(sigranges)){
          if(model_est$Main[i]>=sigranges[[j]][1] & model_est$Main[i]<sigranges[[j]][2]){
            model_est$Sig[i]=1
          }
        }
        
        
      }
      model_est$Sig[-which(model_est$Sig==1)]=0
      if(length(which(model_est$Sig==1))==0){
        model_est$Sig=0
      }
      model_plot2=model_plot+ 
        geom_line(data=model_est,aes(Main,estimate__,color=(Sig)),size=1)+
        scale_color_gradient2(low="black", mid="black",high="cyan" )+
        ylab(Response)+xlab(main)+
        theme_classic()+ guides(color="none")
      
      assign(output,model_plot2, envir = parent.frame())
      output2=gsub("_plot", "", output)
      output2=paste("VOI",output2,sep="_")
      if(length(which(model_est$Sig==1))>0){
        VOIdat=model_est[which(model_est$Sig==1),]
        assign(output2,VOIdat, envir = parent.frame())
      }
    }
    return(model_plot2)
    
  }
  
}

## extracting ranges of derivative at one side of 0
deriv_ranges <- function(der_data_50_1, der_data_50_2, der_data_70_1, der_data_70_2, der_data_90_1, der_data_90_2, factorlevels, modelname, seventy = TRUE, ninety = TRUE){
  # supply all derivative dataframes
  # levels of factor
  # name of model
  
  der_data_50_1$Significance <- ifelse((sign(der_data_50_1$lower) <0 & sign(der_data_50_1$upper)<0) | (sign(der_data_50_1$lower) >0 & sign(der_data_50_1$upper)>0), 1, 0)
  der_data_50_2$Significance <- ifelse((sign(der_data_50_2$lower) <0 & sign(der_data_50_2$upper)<0) | (sign(der_data_50_2$lower) >0 & sign(der_data_50_2$upper)>0), 1, 0)
  der_data_50_1$factor <- factorlevels[1]
  der_data_50_2$factor <- factorlevels[2]
  der_data_50_1$confidence <- 50
  der_data_50_2$confidence <- 50 
  
  if(seventy == TRUE){
    der_data_70_1$Significance <- ifelse((sign(der_data_70_1$lower) <0 & sign(der_data_70_1$upper)<0) | (sign(der_data_70_1$lower) >0 & sign(der_data_70_1$upper)>0), 1, 0)
    der_data_70_2$Significance <- ifelse((sign(der_data_70_2$lower) <0 & sign(der_data_70_2$upper)<0) | (sign(der_data_70_2$lower) >0 & sign(der_data_70_2$upper)>0), 1, 0)
    der_data_70_1$factor <- factorlevels[1]
    der_data_70_2$factor <- factorlevels[2]
    der_data_70_1$confidence <- 70
    der_data_70_2$confidence <- 70 
  }
  if(ninety==TRUE){
    der_data_90_1$Significance <- ifelse((sign(der_data_90_1$lower) <0 & sign(der_data_90_1$upper)<0) | (sign(der_data_90_1$lower) >0 & sign(der_data_90_1$upper)>0), 1, 0)
    der_data_90_2$Significance <- ifelse((sign(der_data_90_2$lower) <0 & sign(der_data_90_2$upper)<0) | (sign(der_data_90_2$lower) >0 & sign(der_data_90_2$upper)>0), 1, 0)
    der_data_90_1$factor <- factorlevels[1]
    der_data_90_2$factor <- factorlevels[2]
    der_data_90_1$confidence <- 90
    der_data_90_2$confidence <- 90
  }
  
  if(seventy==FALSE & ninety == FALSE){
    der_data <- rbind(der_data_50_1, der_data_50_2)
  }
  
  if(seventy == TRUE & ninety == FALSE){
    der_data <- rbind(der_data_50_1, der_data_50_2, der_data_70_1, der_data_70_2)
  }
  
  if(seventy == TRUE & ninety == TRUE){
    der_data <- rbind(der_data_50_1, der_data_50_2, der_data_70_1, der_data_70_2, der_data_90_1, der_data_90_2)
  }
  
  assign(paste(modelname, "overlay", sep = "_"), der_data, envir = parent.frame())
  
}
  return(derivplot)
}

