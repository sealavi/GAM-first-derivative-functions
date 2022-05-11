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
  return(derivplot)
}

