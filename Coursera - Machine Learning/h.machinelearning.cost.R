# Source
message("h.machinelearning.cost :: v0.1: 2023 Nov 24")

# Function
h.machinelearning.cost <- function(x, y, w, b, model = "linear", lambda_w = 0, lambda_b = 0) {

  "Compute the cost of a regression.
    
    Args:
      x           : Predictive values (rows = examples, columns = features)
      y           : Predicted, target values (single vector/column)
      w,b (scalar): Model parameters  𝑓𝑤,𝑏(𝑥)=𝑤𝑥+𝑏 (with either single or multiple predictors 𝑤*𝑥)
      model (char): Regression model: 'linear', 'logistic'
      lambda_w    : Regularization parameter for feature coefficients w
      lambda_b    : Regularization parameter for b
    
    Returns:
      cost (float): The cost of using w,b as the parameters for linear regression
                    to fit the data points y using weighted transformations of x"
  


  # Initialize cost
  cost_sum = 0
  reg_w_sum = 0
  reg_b_sum = 0
  
  # Single predictor/feature
  if (is.null(dim(x))) {
    
    # Setup
    n_examples = length(x)
    if (!is.numeric(x) | !is.numeric(y) ) { stop("For single feature cost computation, both x and y need to be numeric vectors.") }
    
    # Compute sum of all errors
    if (model == "linear") {
      
      # Cost
      for (i in 1:n_examples) {
        f_wb_i = w * x[i] + b           # Product of ith example plus bias
        error = (f_wb_i - y[i]) ** 2    # Error of ith example's predicted score
        cost_sum = cost_sum + error     # Add ith error to total cost sum
      }
      
    }
    if (model == "logistic") {
      for (i in 1:n_examples) {
        z_i = w * x[i] + b              # Product of ith example plus bias
        f_wb_i = 1 / (1 + exp(-z_i))    # Sigmoid
        error = -y[i] * log(f_wb_i) - (1 - y[i]) * log(1 - f_wb_i)  # Loss of ith example's predicted score
        cost_sum = cost_sum + error     # Add ith loss to total cost sum
      }
    }
  }
  
  # Multiple predictors/features
  if (!is.null(dim(x))) {
    
    # Setup
    n_examples = dim(x)[1]
    n_features = dim(x)[2]
    if (!is.data.frame(x) | !is.numeric(y)) { stop("For multiple feature cost computation, x should be a data frame and y should be a numeric vector.") }
    
    # Compute sum of all errors
    if (model == "linear") {
      
      # Cost
      for (i in 1:n_examples) {
        f_wb_i = sum(x[i,] * w) + b     # Dot product of ith example plus bias
        error = (f_wb_i - y[i]) ** 2    # Error for ith example's predicted score
        cost_sum = cost_sum + error     # Add ith error to total cost sum
      }
      
    }
    if (model == "logistic") {
      for (i in 1:n_examples) {
        z_i = sum(x[i,] * w) + b        # Dot product of ith example plus bias
        f_wb_i = 1 / (1 + exp(-z_i))    # Sigmoid
        error = -y[i] * log(f_wb_i) - (1 - y[i]) * log(1 - f_wb_i)  # Loss of ith example's predicted score
        cost_sum = cost_sum + error     # Add ith loss to total cost sum
      }
    }
  }
  
  # Regularization
  reg_w_sum = sum(w^2)
  reg_b_sum = b^2
  
  # Compute average total cost and return
  if (model == "linear") { cost = cost_sum / (2*n_examples) }
  if (model == "logistic") { cost = cost_sum / n_examples }
  
  reg_w = (lambda_w/(2*n_examples)) * reg_w_sum
  reg_b = (lambda_b/(2*n_examples)) * reg_b_sum
  
  total_cost = cost + reg_w + reg_b
  return(total_cost)
  
}

