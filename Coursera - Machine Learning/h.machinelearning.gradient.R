# Source
message("h.machinelearning.gradient :: v0.1: 2023 Nov 24")

# Function
h.machinelearning.gradient <- function(x, y, w, b, method = NULL, model = "linear", lambda_w = 0, lambda_b = 0) {

  "Compute the cost gradient of a regression.
  
    Args:
      x               : Predictive values (rows = examples, columns = features)
      y               : Predicted, target values (single vector/column)
      w,b (scalar)    : Model parameters  𝑓𝑤,𝑏(𝑥)=𝑤𝑥+𝑏 (with either single or multiple predictors 𝑤*𝑥)
      method (char)   : Method of computation
      model (char)    : Regression model: 'linear', 'logistic'
      lambda_w        : Regularization parameter for feature coefficients w
      lambda_b        : Regularization parameter for b
  
    Returns:
      dj_dw (scalar)  : The gradient of the cost w.r.t. the parameters w
      dj_db (scalar)  : The gradient of the cost w.r.t. the parameter b"
  
  
  # Single predictor/feature
  if (is.null(dim(x))) {
    
    # Setup
    n_examples = length(x)
    n_features = 1
    dj_dw = 0
    dj_db = 0
    
    # Compute gradient
    for (i in 1:n_examples) {
      
      # Model prediction with current weight and bias
      if (model == "linear")    { f_wb_i = w * x[i] + b }    
      if (model == "logistic")  { f_wb_i = 1 / (1 + exp(-(w * x[i] + b))) }
      
      error_i = f_wb_i - y[i]     # Error for the ith example
      dj_dw_i = error_i * x[i]    # Error times ith example, part of the partial derivative of w
      dj_db_i = error_i           # Error for b is simply error_i
      
      dj_dw = dj_dw + dj_dw_i     # Update w and b by adding the ith error to running sum
      dj_db = dj_db + dj_db_i
      
    }
  }
  
  # Multiple predictors/features
  if (!is.null(dim(x))) {
    
    # Setup
    n_examples = dim(x)[1]
    n_features = dim(x)[2]
    dj_dw = array(0, c(1,n_features))
    dj_db = 0
    
    # Compute gradient
    for (i in 1:n_examples) {
      
      # Model prediction with current weight and bias
      if (model == "linear")    { f_wb_i = sum(w * x[i,]) + b } # Dot product + b
      if (model == "logistic")  { f_wb_i = 1 / (1 + exp(- (sum(w * x[i,]) + b) )) } # Dot product + b.. fed into the sigmoid function
      
      error_i = f_wb_i - y[i]         # Error for the ith example
      for (j in 1:n_features) {
        dj_dw[j] = dj_dw[j] + error_i * x[i,j]
      }
      dj_db = dj_db + error_i
    }
  }
  
  # Compute average over n_examples and return
  dj_dw = dj_dw / n_examples
  dj_db = dj_db / n_examples
  
  # Regularization
  for (i in 1:n_features) { dj_dw[i]  = dj_dw[i] + (lambda_w/n_examples) * w[i] }
  dj_db = dj_db + (lambda_b/(n_examples)) * b
  
  return(list(dj_dw, dj_db))
    
}