h.machinelearning.gradient <- function(x, y, w, b, method = NULL) {

  "Compute the gradient for linear regression.
  
    Args:
      x               : Predictive values (rows = examples, columns = features)
      y               : Predicted, target values (single vector/column)
      w,b (scalar)    : Model parameters  𝑓𝑤,𝑏(𝑥)=𝑤𝑥+𝑏 (with either single or multiple predictors 𝑤*𝑥)
      method (char)   : Method of computation
  
    Returns:
      dj_dw (scalar)  : The gradient of the cost w.r.t. the parameters w
      dj_db (scalar)  : The gradient of the cost w.r.t. the parameter b"
  
  # What are the number of features and number of training examples?
  if (is.null(dim(x)[1])) { # Single predictor/feature
    n_examples = length(x)
    dj_dw = 0
    dj_db = 0
    
    # Compute gradient
    for (i in 1:n_examples) {  
      f_wb = w * x[i] + b            # Regression model with provided parameters
      dj_dw_i = (f_wb - y[i]) * x[i] # Error * ith feature (part of the partial derivative)
      dj_db_i = f_wb - y[i]          # Error for b
      dj_dw = dj_dw+dj_dw_i          # Update w and b
      dj_db = dj_db+dj_db_i
    }
    
    dj_dw = dj_dw / n_examples 
    dj_db = dj_db / n_examples
    
    return(list(dj_dw, dj_db))
    
  } else { # Multiple predictors/features
    n_examples = dim(x)[1]
    n_features = dim(x)[2]
    dj_dw = array(0, c(1,n_features))
    dj_db = 0
    
    # Compute gradient
    for (i in 1:n_examples) {  
      err = (as.numeric(x[i,]) %*% w + b)[1] - y[i] # Error from the model with n_features on the ith example
      for (j in 1:n_features) {
        dj_dw[j] = dj_dw[j] + err * x[i,j]
      }
      dj_db = dj_db + err
    }
    
    dj_dw = dj_dw / n_examples 
    dj_db = dj_db / n_examples
    
    return(list(dj_dw, dj_db))
  }
    
}