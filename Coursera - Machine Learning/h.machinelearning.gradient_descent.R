# Source
message("h.machinelearning.gradient_descent :: v0.2: 2023 Nov 24")

# Function
h.machinelearning.gradient_descent <- function(x, y, w_init, b_init,
                                               model = "linear", 
                                               alpha, n_iters,
                                               scale = FALSE,
                                               max_iters = 1e+10,
                                               conv_crit = 1e-50,
                                               div_crit = 1e+100,
                                               search=FALSE, searchspace=3,
                                               cost.method = NULL, 
                                               cost.lambda_w = 0,
                                               cost.lambda_b = 0,
                                               gradient.method = NULL,
                                               verbose=TRUE, plot=TRUE) {

  " Performs gradient descent to fit w,b. Updates w,b by taking 
    n_iters gradient steps with learning rate alpha.
    
    Args:
      x                     : Predictive values
      y                     : Predicted, target values
      w_init,b_init (scalar): Initial values of model parameters
      model (char)          : Regression model: 'linear', 'logistic'
      alpha (float)         : Learning rate
      n_iters (int)         : Number of iterations to run gradient descent
      scale (bool)          : Scale the features? This is recommended (default = TRUE)
      max_iters (int)       : Maximum number of iterations to run gradient descent
      conv_crit (float)     : Minimum cost criterion under which descent is considered converged.
                              So, if J(w,b) < conv_crit, descent is halted. Default = 1e-100.
      div_crit (float)      : Maximum cost criterion above which descent is considered diverged.
                              So, if J(w,b) > div_crit, descent is halted. Default = 1e+1000.
      search (bool)         : If TRUE, do a further search around the found w and b values
      searchspace (int)     : If search==TRUE, search all combination of w and b integer values
                              between w+/-searchspace*var(desc_hist of w) and b+/-searchspace*var(desc_hist of b)
                              (values > 4 can lead to long searches!)
      cost.method           : Choice of method to be called by h.machinelearning.cost [optional]
      cost.lambda_w         : Regularization parameter for feature coefficients w (lower is more overfit, higher is more underfit).
                              Value must be >= 0.
      cost.lambda_b         : Regularization parameter for b. Value must be >= 0.
      gradient.method       : Choice of method to be called by h.machinelearning.gradient [optional]
      verbose (bool)        : If TRUE, print values of w and b to the console every 2000 iterations
      plot (bool)           : If TRUE, plot the cost function
      
    Returns:
      w (scalar)        : Updated value of parameter after running gradient descent
      b (scalar)        : Updated value of parameter after running gradient descent
      desc_hist (array) : History of cost (J), w, and b (values at each iteration)
      idx (int, array)  : All local minima found with broad search"
  
  # Flex parameters (can be adjusted, if you want)
  sigdigs=6 # Number of significant digits to be used when printing to console
  search_crit=20 # When plotting the cost surface, broad search computes variance around the w,b parameters. If this value is low (var < search_crit), the final plot may look odd. This parameter adjust the cutoff.
  col.surface=NULL # Color vector for surface plot (see options below, NULL = standard viridis scale)
  #col.surface=rev(c("#F2DC5D", "#F2A359", "#DB9065", "#A4031F", "#240B36"))
  #col.surface=rev(c("#F7AEF8", "#B388EB", "#8093F1", "#72DDF7", "#F4F4ED"))
  
  
  ### UNIVARIATE GRADIENT DESCENT ###
  if (is.null(dim(x))) {
    
    # Initialize array for cost (J), w, and b at each iteration (for graphing)
    {
      desc_hist = array(NA, c(n_iters,3)); colnames(desc_hist) = c("J","w","b")
      n_examples = length(x)
      b = b_init
      w = w_init
      
      # Message to console
      message(paste0(
        "Starting univariate gradient descent :::
        model      :   ",model,"
        alpha      :   ",alpha,"
        n examples :   ",n_examples,"
        n features :   1
        initial b  :   ",b,"
        initial w  :   ",w))
    }
    
    # Scale feature
    if (scale) {
    message("NOTE: Scaling requested, but there's only one feature. No scaling applied.
      You can scale the feature yourself, if you insist.")
    }
    
    # Gradient descent
    for (i in 1:n_iters) {
      
      # Calculate the gradient and update the parameters using gradient_function
      gradient = h.machinelearning.gradient(x, y, w, b, model = model)
      dj_dw=gradient[[1]]
      dj_db=gradient[[2]]
      
      # Update parameters
      w = w - alpha * dj_dw
      b = b - alpha * dj_db
      
      # Save cost J and parameters w & b at each iteration
      if (i < max_iters) { # prevent resource exhaustion 
        desc_hist[i,"J"] = h.machinelearning.cost(x, y, w, b, model=model, lambda_w=cost.lambda_w, lambda_b=cost.lambda_b) # Cost with w & b on this iteration
        desc_hist[i,"w"] = w # w on this iteration
        desc_hist[i,"b"] = b # b on this iteration
        
        if (is.na(desc_hist[i,"J"])) { stop("h.machinelearning.cost returned a cost of NaN. Perhaps try a different alpha.") }
        
      } else { stop(paste0("Failed to converge. Search reached the maximum number of iterations (",max_iters,")!
                          Set max_iters to a higher value to continue trying with current parameters. Or, try
                          different parameters.")) }
      
      # Do we have convergence?
      if (desc_hist[i,"J"] < conv_crit) {
        
        # Remove all NA
        desc_hist<-na.omit(desc_hist)
        
        message(paste0("
    Gradient descent has converged after ",i," iteration(s)."))
        break
      }
          
      # Do we have divergence?
      if (any(desc_hist[i,] > div_crit | desc_hist[i,] == -Inf)) {
        
        # Remove all NA
        desc_hist<-na.omit(desc_hist)
  
        stop(paste0("
    Gradient descent has diverged after ",i," iteration(s).
    w (",signif(desc_hist[i,2],digits=sigdigs),") and b (",signif(desc_hist[i,3],digits=sigdigs),"), with a cost (J_wb) of ",signif(desc_hist[i,1],digits=sigdigs),"."))
      }
      
      # Print cost every at intervals 10 times or as many iterations if < 10
      if (verbose) {
        if ( i == 1 | i %% ceiling(n_iters/5) == 0) {
          message(paste0("Iteration ",i," ::"))
          message(paste0("     dj_dw: ",signif(dj_dw,digits=sigdigs),", dj_db: ",signif(dj_db,digits=sigdigs)))
          message(paste0("     w: ",signif(w,digits=sigdigs),", b: ",signif(b,digits=sigdigs)))
          message(paste0("     Cost (J_wb): ",signif(desc_hist[i,1],digits=sigdigs)))
        }
      }
      
    } # Descent end
    
    # Print final results to console
    {
    message(paste0("
                 
    Gradient descent (",model,") found w (",signif(w,digits=sigdigs),") and b (",signif(b,digits=sigdigs),"), with a cost (J_wb) of ",signif(desc_hist[i,1],digits=sigdigs),"."))
    }
      
    if (scale) {
      message("
    NOTE: This is another reminder that scaling was NOT applied!
          Scale the feature (x) yourself, if you really insist.")
    }
    
    # Compute broader cost surface
    if (search) {
      
      # Is searchspace too big?
      w_var=var(desc_hist[,2])
      b_var=var(desc_hist[,3])
      
      if (w_var < search_crit | b_var < search_crit) {
        if (w_var < search_crit) { w_var=50 }
        if (b_var < search_crit) { b_var=50 }
        message("
    Broad search requested. Note that there was low variance in the w,b history trace. Attempting to adjust plotting parameters, but it may look odd! Computing now..")
        
      } else {
        message("
    Broad search requested. Computing now..")
      }
      
      # Initialize search
      lowerbound_w=floor(w)-floor(w_var*searchspace); upperbound_w=floor(w)+floor(w_var*searchspace)
      lowerbound_b=floor(b)-floor(b_var*searchspace); upperbound_b=floor(b)+floor(b_var*searchspace)
      grad_search<-data.frame(
        w=rep(lowerbound_w:upperbound_w,each=length(lowerbound_b:upperbound_b)),
        b=rep(lowerbound_b:upperbound_b),
        J=NA
      )
      
      # Search this space..
      for (i in 1:dim(grad_search)[1]) {
        grad_search[i,3] = h.machinelearning.cost(x, y, grad_search[i,1], grad_search[i,2], model=model, lambda_w=cost.lambda_w, lambda_b=cost.lambda_b)
      }
      
      # Do the results (approximately) match?
      idx=which(grad_search$J==min(grad_search$J))
      
      message(paste0("
                 
    Broad search found w (",signif(grad_search[idx,1],digits=sigdigs),") and b (",signif(grad_search[idx,2],digits=sigdigs),"), with a cost (J_wb) of ",signif(grad_search[idx,3],digits=sigdigs),"."))
      
      # Make cost surface matrix
      surface <- matrix(grad_search$J, nrow = length(lowerbound_w:upperbound_w), ncol = length(lowerbound_b:upperbound_b))
      mat_surface <- list(as.vector(lowerbound_w:upperbound_w),
                          as.vector(lowerbound_b:upperbound_b),
                          surface)
      
    } # Search end
    
    
    # Plot cost?
    if (plot) {
      
      # With low n_iters, make one cost curve
      if (n_iters < 2000) {
        plot(1:n_iters, desc_hist[,1], type="l",
             main="Cost over iterations", xlab="Iteration", ylab="Cost")
        
        # If we have high n_iters...
      } else {
        if (search) { # ..and we ran a broad search..
          par(mfrow=c(2,1))
          splitval=1000
          
          # Early search
          plot(1:splitval, desc_hist[1:splitval,1], type="l", 
               main="Early search", xlab="Iteration", ylab="Cost")
          
          # Late search
          plot(splitval:n_iters, desc_hist[splitval:n_iters,1], type="l", 
               main="Late search", xlab="Iteration", ylab="Cost")
          
          # 2d surface
          ggplot(grad_search, aes(w, b, z=J)) + 
            geom_tile() + 
            geom_contour_filled(colour="white")
          
          # 3d surface
          descent_trace=data.frame(cost=desc_hist[,1], w=desc_hist[,2], b=desc_hist[,3])
          search_dot=data.frame(w=grad_search[idx,1], b=grad_search[idx,2], cost=grad_search[idx,3])
          J=xtabs(J~w+b, data=grad_search)
          suppressWarnings(print(
            plot_ly(z=~J, type="surface", colors=col.surface, contours = list(
            z = list(show=TRUE,usecolormap=TRUE,highlightcolor="lightblue",project=list(z=TRUE)),
            y = list(show=TRUE,usecolormap=TRUE,highlightcolor="lightblue",project=list(y=TRUE)),
            x = list(show=TRUE,usecolormap=TRUE,highlightcolor="lightblue",project=list(x=TRUE)))) %>%
              add_trace(x=~b, y=~w, z=~cost, data=descent_trace,
                        type="scatter3d", mode="lines", 
                        line=list(color="red", width=4)) %>%
              add_trace(x=~b, y=~w, z=~cost, data=search_dot,
                        type="scatter3d", mode="markers",
                        marker=list(color="red", size=10)) %>%
              layout(title = 'Cost Surface with Descent Trace',
                     scene = list(xaxis=list(title="b"), yaxis=list(title="w"))) %>%
              hide_colorbar() %>% hide_legend()))
          
        } else { # ..and we didn't run a broad search..
          par(mfrow=c(2,1))
          splitval=1000
          
          # Early search
          plot(1:splitval, desc_hist[1:splitval,1], type="l",
               main="Early search", xlab="Iteration", ylab="Cost")
          
          # Late search
          plot(splitval:n_iters, desc_hist[splitval:n_iters,1], type="l",
               main="Late search", xlab="Iteration", ylab="Cost")
          
        } # End high iters, no search
      } # End high iters, with search
    } # End plot
    
    # Finished!
    if (search) { return(list(w, b, desc_hist, idx)) } else { return(list(w, b, desc_hist)) }
    
    
  }
  
  ### MULTIVARIATE GRADIENT DESCENT ###  
  if (!is.null(dim(x))) {
    
    # Setup
    n_examples = dim(x)[1]
    n_features = dim(x)[2]
    
    # Initialize array for cost (J), w, and b at each iteration (for graphing)
    {
    w.labels=NULL; for (i in 1:n_features) { w.labels = append(w.labels, paste0("w_",i)) }
    desc_hist = array(NA, c(n_iters,n_features+2)); colnames(desc_hist) = c("J",w.labels,"b")
    b = b_init
    w = w_init
    
    # Message to console
     message(paste0(
      "Starting multivariate gradient descent :::
        model      :   ",model,"
        alpha      :   ",alpha,"
        n examples :   ",n_examples,"
        n features :   ",n_features,"
        initial b  :   ",b))
    for (i in 1:n_features) {
    message(paste0(
      "        initial w_",i,":   ",w[i]))
    }
    
    # Indices for last of all w's and the b (to improve readability in subsequent sections)
    idx_w = 2:(n_features+1)
    idx_b = (n_features+2)
    }
    
    # Scale features
    if (scale) {
      message("Scaling features..")
      for (i in 1:n_features) { 
      message("On x_",i,", unscaled |max-min| = ",diff(range(x[,i]))) }
      x = as.data.frame(scale(x)) # h.machinelearning.gradient and _.cost expect a data frame
      for (i in 1:n_features) { 
      message("On x_",i,", the newly scaled |max-min| = ",diff(range(x[,i]))) }
      message("The range of the scaled features should now be more similar. Is that right?")
    } else {
      message("NOTE: Features are not being scaled! Is that really what you want?")
      for (i in 1:n_features) { 
      message("On x_",i,", |max-min| = ",diff(range(x[,i]))) }
      message("If the range of the feature values are very different, scaling is strongly recommended!")
    }
    
    # Gradient descent
    for (i in 1:n_iters) {
      
      # Calculate the gradient and update the parameters using gradient_function
      gradient = h.machinelearning.gradient(x, y, w, b, model = model)
      dj_dw = gradient[[1]]
      dj_db = gradient[[2]]
      
      # Update parameters
      w = w - alpha * dj_dw
      b = b - alpha * dj_db
      
      # Save cost J and parameters w & b at each iteration
      if (i < max_iters) { # prevent resource exhaustion 
        desc_hist[i,"J"] = h.machinelearning.cost(x, y, w, b, model=model, lambda_w=cost.lambda_w, lambda_b=cost.lambda_b) # Cost with w & b on this iteration
        desc_hist[i,idx_w] = w # w on this iteration
        desc_hist[i,idx_b] = b # b on this iteration
        
        if (is.na(desc_hist[i,"J"])) { stop("h.machinelearning.cost returned a cost of NaN. Perhaps try a different alpha.") }
        
      } else { stop(paste0("Failed to converge. Search reached the maximum number of iterations (",max_iters,")!
                          Set max_iters to a higher value to continue trying with current parameters. Or, try
                          different parameters.")) }
      
      # Do we have convergence?
      if (desc_hist[i,"J"] < conv_crit) {
        
        # Remove all NA
        desc_hist<-na.omit(desc_hist)
        
        if (i == 1) {
          message(paste0("
    Gradient descent has converged after 1 iteration."))
        } else {
          message(paste0("
    Gradient descent has converged after ",i," iterations."))}
        break
      }
      
      # Do we have divergence?
      if (desc_hist[i,"J"] > div_crit | abs(desc_hist[i,"J"]) == Inf) {
        
        # Remove all NA
        desc_hist<-na.omit(desc_hist)
        
        # Divergence message
        message("
                
                Divergence criterion reached!
                ")
        message(paste0(" w: ",signif(desc_hist[i,idx_w],digits=sigdigs)))
        message(paste0(" b: ",signif(desc_hist[i,idx_b],digits=sigdigs)))
        message(paste0(" J: ",signif(desc_hist[i,1],digits=sigdigs)))
        
        stop(paste0("
    Gradient descent has diverged after ",i," iteration(s)."))
      }
      
      # Print cost every at intervals 10 times or as many iterations if < 10
      if (verbose) {
        if (i == 1 | i %% ceiling(n_iters/5) == 0) {
          message(paste0("Iteration ",i," ::"))
          message(paste0("     dj_dw: ",signif(dj_dw,digits=sigdigs)))
          message(paste0("     dj_db: ",signif(dj_db,digits=sigdigs)))
          message(paste0("     w: ",signif(w,digits=sigdigs)))
          message(paste0("     b: ",signif(b,digits=sigdigs)))
          message(paste0("     Cost (J_wb): ",signif(desc_hist[i,1],digits=sigdigs)))
        } # Message end
      } # Verbose end
    } # Descent end

    # Print final results to console
    {
    message(paste0("
    
    Gradient descent (",model,") found the following values:"))
    for (i in 1:n_features) {
    message(paste0("      w_",i,":  ",signif(w[i],digits=sigdigs)))
    }
    message(paste0("      b:    ",signif(b,digits=sigdigs)))
    message(paste0("      J:    ",signif(desc_hist[n_iters,1],digits=sigdigs)))
    }
    
    # Compute broader cost surface
    if (search) {
      
      # Search bounds
      if (n_features <= 5 ) { searchspace_multivar=c(3,.1) } else { searchspace_multivar=c(5,.01) }
      message(paste0("
    Broad search requested. This can take a while! Computing now with searchspace_multivar = ",searchspace_multivar[1]," and ",searchspace_multivar[2],"."))
      
      # Initialize search
      foo = list() # Create empty list for eventual data frame
      search_bounds=array(NA, c(2,n_features+1)) # Search bounds array
      highvar_bound=5
      for (i in 1:(n_features+1)) { # For all features and b..
        
        if (i < (n_features+1)) { # All features w_i
          
          # Variance on multivariate approaches can be low. If so, adjust.
          w_var=var(desc_hist[,(1+i)])
          w_mean=mean(desc_hist[,(1+i)])
          if (w_var < 1) {
            message(paste0("Low variance found on w_",i," (",signif(w_var,digits=sigdigs),"), which is not uncommon with mulivariate gradient descent. Adjusting search space accordingly.."))
            
            search_bounds[1,i] = floor(w_mean) # lower for each feature w
            search_bounds[2,i] = ceiling(w_mean) # upper for each feature w
            foo[[i]] = c(seq(search_bounds[1,i],search_bounds[2,i],searchspace_multivar[2])) # Append final data frame
          } else if (w_var < 10) {
            message(paste0("Variance found on w_",i," was < 10 (",signif(w_var,digits=sigdigs),"). Adjusting search space accordingly.."))
            
            search_bounds[1,i] = floor(w_mean)-floor(w_var*searchspace_multivar[1]) # lower for each feature w
            search_bounds[2,i] = floor(w_mean)+floor(w_var*searchspace_multivar[1]) # upper for each feature w
            foo[[i]] = c(seq(search_bounds[1,i],search_bounds[2,i])) # Append final data frame
          } else {
            message(paste0("Variance found on w_",i," was > 10 (",signif(w_var,digits=sigdigs),"). Adjusting search space accordingly.."))
            
            search_bounds[1,i] = floor(desc_hist[n_iters,(1+i)])-highvar_bound # lower for each feature w
            search_bounds[2,i] = floor(desc_hist[n_iters,(1+i)])+highvar_bound # upper for each feature w
            foo[[i]] = c(seq(search_bounds[1,i],search_bounds[2,i])) # Append final data frame
          }
        } else { # b..
          
          b_var=var(desc_hist[,(1+i)])
          b_mean=mean(desc_hist[,(1+i)])
          if (b_var < 1) {
            message(paste0("Low variance found on b (",signif(b_var,digits=sigdigs),"), which is not uncommon with mulivariate gradient descent. Adjusting search space.."))
            
            search_bounds[1,i] = floor(b_mean) # lower for b
            search_bounds[2,i] = ceiling(b_mean) # upper for b
            foo[[i]] = c(seq(search_bounds[1,i],search_bounds[2,i],searchspace_multivar[2])) # Append final data frame
          } else if (b_var < 10) {
            message(paste0("Variance found on b was < 10 (",signif(b_var,digits=sigdigs),"). Adjusting search space accordingly.."))
            
            search_bounds[1,i] = floor(b_mean)-floor(b_var*searchspace_multivar[1]) # lower b
            search_bounds[2,i] = floor(b_mean)+floor(b_var*searchspace_multivar[1]) # upper b
            foo[[i]] = c(seq(search_bounds[1,i],search_bounds[2,i])) # Append final data frame
            
          } else {
            message(paste0("Variance found on b was > 10 (",signif(b_var,digits=sigdigs),"). Adjusting search space accordingly.."))
            
            search_bounds[1,i] = floor(desc_hist[n_iters,(2+n_features)])-highvar_bound # lower for each feature w
            search_bounds[2,i] = floor(desc_hist[n_iters,(2+n_features)])+highvar_bound # upper for each feature w
            foo[[i]] = c(seq(search_bounds[1,i],search_bounds[2,i])) # Append final data frame
          }
        }
      }


      # Format final data frame for the search
      foo = h.format.expand_grid(foo) # Expand grid to form all unique combinations
      grad_search = data.frame(J = rep(NA,dim(foo)[1])) # Set up new data frame, add features+b, and relabel
      grad_search = cbind(grad_search,foo)
      colnames(grad_search) = c("J",w.labels,"b")
      
      # Search this space..
      for (i in 1:dim(grad_search)[1]) {
        w_vec = NULL
        for (j in 1:n_features) {
          w_vec = append(w_vec, grad_search[i,j+1])
        }
        grad_search$J[i] = h.machinelearning.cost(x, y, w_vec, grad_search$b[i], model=model, lambda_w=cost.lambda_w, lambda_b=cost.lambda_b)
        
        # Print cost every at intervals 10 times or as many iterations if < 10
        if (verbose) {
          if ( i == 1 | i %% ceiling(n_iters/5) == 0) {
            message(paste0("Iteration ",i,"/",dim(grad_search)[1]," (~",signif((i/dim(grad_search)[1])*100,digits=5)," %) ::"))
            message(paste0("     w: ",signif(grad_search[i,2:(1+n_features)],digits=sigdigs)))
            message(paste0("     b: ",signif(grad_search$b[i],digits=sigdigs)))
            message(paste0("     Cost (J_wb): ",signif(grad_search$J[i],digits=sigdigs)))
          } # Message end
        } # Verbose end
      }
      
      # Do the results (approximately) match?
      idx=which(grad_search$J==min(grad_search$J))
      
      # If multiple minima
      if (length(idx) > 1) {
        message(paste0("
                 
    Broad search found multiple minima, with a minimal cost (J_wb) of ",signif(min(grad_search$J),digits=sigdigs),"."))
 
      } else {
        
        message(paste0("
                 
    Broad search found the following values:"))
        for (i in 1:n_features) {
          message(paste0("      w_",i,":  ",signif(grad_search[idx,(1+i)],digits=sigdigs)))
        }
        message(paste0("      b:    ",signif(grad_search$b[idx],digits=sigdigs)))
        message(paste0("      J:    ",signif(grad_search$J[idx],digits=sigdigs)))  
        
      }
    }
    
    # Plot cost?
    if (plot) {
      
      # With low n_iters, make one cost curve
      if (n_iters < 2000) {
        plot(1:n_iters, desc_hist[,1], type="l",
             main="Cost over iterations", xlab="Iteration", ylab="Cost")
        
        # If we have high n_iters..
      } else {
        if (search) { # ..and we ran a broad search..
          par(mfrow=c(2,1))
          splitval=1000
          
          # Early search
          plot(1:splitval, desc_hist[1:splitval,1], type="l", 
               main="Early search", xlab="Iteration", ylab="Cost")
          
          # Late search
          plot(splitval:n_iters, desc_hist[splitval:n_iters,1], type="l", 
               main="Late search", xlab="Iteration", ylab="Cost")
          
          # 2d surface
          for (i in 1:n_features) {
            w_i=grad_search[,i+1]
            ggplot(grad_search, aes(w_i, b, z=J)) + 
              geom_tile() + 
              geom_contour_filled(colour="white")
          }
          
          # 3d surface
          for (i in 1:n_features) {
            descent_trace=data.frame(cost=desc_hist[,1], w=desc_hist[,(1+i)], b=desc_hist[,(2+n_features)])
            search_dot=data.frame(cost=grad_search$J[idx], w=grad_search[idx,(1+i)], b=grad_search$b[idx])
            w_i=grad_search[,(1+i)]
            J=xtabs(J~w_i+b, data=grad_search)
            suppressWarnings(print(
              plot_ly(z=~J, type="surface", colors=col.surface, contours = list(
              z = list(show=TRUE,usecolormap=TRUE,highlightcolor="lightblue",project=list(z=TRUE)),
              y = list(show=TRUE,usecolormap=TRUE,highlightcolor="lightblue",project=list(y=TRUE)),
              x = list(show=TRUE,usecolormap=TRUE,highlightcolor="lightblue",project=list(x=TRUE)))) %>%
                add_trace(x=~b, y=~w, z=~cost, data=descent_trace,
                          type="scatter3d", mode="lines", 
                          line=list(color="red", width=4)) %>%
                add_trace(x=~b, y=~w, z=~cost, data=search_dot,
                          type="scatter3d", mode="markers",
                          marker=list(color="red", size=10)) %>%
                layout(title = 'Cost Surface with Descent Trace',
                       scene = list(xaxis=list(title="b",range=c(search_bounds[1,dim(search_bounds)[2]],search_bounds[2,dim(search_bounds)[2]])),
                                    yaxis=list(title=paste0("w_",i),range=c(search_bounds[1,i],search_bounds[2,i])))) %>%
                hide_colorbar() %>% hide_legend()))
          }
        } else { # ..and we didn't run a broad search..
          par(mfrow=c(2,1))
          splitval=1000
          
          # Early search
          plot(1:splitval, desc_hist[1:splitval,1], type="l",
               main="Early search", xlab="Iteration", ylab="Cost")
          
          # Late search
          plot(splitval:n_iters, desc_hist[splitval:n_iters,1], type="l",
               main="Late search", xlab="Iteration", ylab="Cost")
          
        } # End high iters, no search
      } # End high iters, with search
    }
    
    # Finished!
    if (search) { return(list(w, b, desc_hist, idx)) } else { return(list(w, b, desc_hist)) }
    
  }
}


  