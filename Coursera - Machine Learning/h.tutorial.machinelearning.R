### h.lab tutorial: machine learning
# Based on the Coursera 'Machine Learning' course
# https://www.coursera.org/learn/machine-learning/home/

# Dependencies
library(ggplot2)
library(plotly)

# h.lab
hlab_dir='~/Desktop/Research/+h.lab/'
source(paste0(hlab_dir,'h.format.expand_grid.R'))
source(paste0(hlab_dir,'h.machinelearning.cost.R'))
source(paste0(hlab_dir,'h.machinelearning.gradient.R'))
source(paste0(hlab_dir,'h.machinelearning.gradient_descent.R'))




##################
##### Week 1 #####
##################


"Gradient descent for linear regression (C1_W1_Lab04)"
# In this step, we will learn to automate optimization of w and b using gradient descent.
# We will start with a 'simple' (single variable) linear regression.
{
  # Our training set contains information on the size of a house (1000 sqft) and its price
  # (in thousands of dollars)
  training_set <- data.frame(
    size = c(1,2),
    price = c(300,500)
  )
  
  # Run gradient descent
  grad_desc = h.machinelearning.gradient_descent(x=training_set$size, y=training_set$price,
                                                 w_in=0, b_in=0, 
                                                 alpha=1.0e-2, num_iters=1e4, 
                                                 verbose=T, search=T, plot=T)
  
}


"VECTORS ::
Vectors, as you will use them in this course, are ordered arrays of numbers. In notation, vectors 
are denoted with lower case bold letters such as 𝐱. The elements of a vector are all the same type.
A vector does not, for example, contain both characters and numbers. The number of elements in the 
array is often referred to as the 'dimension', though mathematicians may prefer 'rank'. The vector shown 
has a dimension of𝑛.

vector = c(77,2,3,10,99,...,n)

The elements of a vector can be referenced with an index. In math settings, 
indexes typically run from 1 to n. In computer science, indexing will often run from 0 to n-1.
In notation, elements of a vector, when referenced individually will indicate the 
index in a subscript, for example, the 0th element, of the vector 𝐱 is 𝑥_0. Note, the x is not bold 
in this case. R is 1-based and so is Matlab, meaning that indexing starts with 1, but Python is 0-based, for example.
This means that in the vector above, to access the first entry in the vector in R, you'd call 'vector[1]'.
In 'Matlab, vector(1)'. But, in Python, it'd be 'vector[0]' to get the number 77 back (since 77 is the corresponding entry)."



"Multiple Variable Linear Regression (C1_W2_Lab02)"
# Now, we will move to a more complex (multiple variable) linear regression.
{
  # You will use the motivating example of housing price prediction.
  # The training dataset contains three examples with four features (size, bedrooms, floors, and age of home)
  # shown in the table below. Note that, unlike the previous example, size is in sqft rather than 1000 sqft.
  # This causes an issue, which you will solve in the next section!"
  training_set <- data.frame(
    size = c(2104,1416,852),
    rooms = c(5,3,2),
    floors = c(1,2,1),
    age = c(45,40,35),
    price = c(460,232,178)
  )
  
  # What are the features and what's being predicted?
  features_set=training_set[,c("size","rooms","floors","age")]
  target_set=training_set[,"price"]
  
  # Initialize parameters
  b_init = 785.1811367994083
  w_init = c(0.39133535, 18.75376741, -53.36032453, -26.42131618) # One for each feature, in order of columns in the features set
  #b_init = 0
  #w_init = c(0,0,0,0)
  
  # Run gradient descent
  grad_desc = h.machinelearning.gradient_descent(x=features_set, y=target_set,
                                                 w_in=w_init, b_in=b_init, 
                                                 alpha=.01, num_iters=3000, 
                                                 verbose=TRUE, search=FALSE, plot=TRUE, scale=TRUE,
                                                 model_type = "linear", max_iters = 1e+100,
                                                 conv_crit = 1e-100, div_crit = 1e+1000, searchspace=3, 
                                                 cost.method = NULL, gradient.method = NULL)
  
}



















# Grad descent debug
x=training_set$size; y=training_set$price
w_init=0; b_init=0;
alpha=1.0e-2; num_iters=1e4;
verbose=T; search=T; plot=T
scale=T
searchspace=3
max_iters = 1e+100
conv_crit = 1e-100
div_crit = 1e+1000













if (TRUE) {
x=list(c(1.24e+03, 3.00e+00, 1.00e+00, 6.40e+01),
      c(1.95e+03, 3.00e+00, 2.00e+00, 1.70e+01),
      c(1.72e+03, 3.00e+00, 2.00e+00, 4.20e+01),
      c(1.96e+03, 3.00e+00, 2.00e+00, 1.50e+01),
      c(1.31e+03, 2.00e+00, 1.00e+00, 1.40e+01),
      c(8.64e+02, 2.00e+00, 1.00e+00, 6.60e+01),
      c(1.84e+03, 3.00e+00, 1.00e+00, 1.70e+01),
      c(1.03e+03, 3.00e+00, 1.00e+00, 4.30e+01),
      c(3.19e+03, 4.00e+00, 2.00e+00, 8.70e+01),
      c(7.88e+02, 2.00e+00, 1.00e+00, 8.00e+01),
      c(1.20e+03, 2.00e+00, 2.00e+00, 1.70e+01),
      c(1.56e+03, 2.00e+00, 1.00e+00, 1.80e+01),
      c(1.43e+03, 3.00e+00, 1.00e+00, 2.00e+01),
      c(1.22e+03, 2.00e+00, 1.00e+00, 1.50e+01),
      c(1.09e+03, 2.00e+00, 1.00e+00, 6.40e+01),
      c(8.48e+02, 1.00e+00, 1.00e+00, 1.70e+01),
      c(1.68e+03, 3.00e+00, 2.00e+00, 2.30e+01),
      c(1.77e+03, 3.00e+00, 2.00e+00, 1.80e+01),
      c(1.04e+03, 3.00e+00, 1.00e+00, 4.40e+01),
      c(1.65e+03, 2.00e+00, 1.00e+00, 2.10e+01),
      c(1.09e+03, 2.00e+00, 1.00e+00, 3.50e+01),
      c(1.32e+03, 3.00e+00, 1.00e+00, 1.40e+01),
      c(1.59e+03, 0.00e+00, 1.00e+00, 2.00e+01),
      c(9.72e+02, 2.00e+00, 1.00e+00, 7.30e+01),
      c(1.10e+03, 3.00e+00, 1.00e+00, 3.70e+01),
      c(1.00e+03, 2.00e+00, 1.00e+00, 5.10e+01),
      c(9.04e+02, 3.00e+00, 1.00e+00, 5.50e+01),
      c(1.69e+03, 3.00e+00, 1.00e+00, 1.30e+01),
      c(1.07e+03, 2.00e+00, 1.00e+00, 1.00e+02),
      c(1.42e+03, 3.00e+00, 2.00e+00, 1.90e+01),
      c(1.16e+03, 3.00e+00, 1.00e+00, 5.20e+01),
      c(1.94e+03, 3.00e+00, 2.00e+00, 1.20e+01),
      c(1.22e+03, 2.00e+00, 2.00e+00, 7.40e+01),
      c(2.48e+03, 4.00e+00, 2.00e+00, 1.60e+01),
      c(1.20e+03, 2.00e+00, 1.00e+00, 1.80e+01),
      c(1.84e+03, 3.00e+00, 2.00e+00, 2.00e+01),
      c(1.85e+03, 3.00e+00, 2.00e+00, 5.70e+01),
      c(1.66e+03, 3.00e+00, 2.00e+00, 1.90e+01),
      c(1.10e+03, 2.00e+00, 2.00e+00, 9.70e+01),
      c(1.78e+03, 3.00e+00, 2.00e+00, 2.80e+01),
      c(2.03e+03, 4.00e+00, 2.00e+00, 4.50e+01),
      c(1.78e+03, 4.00e+00, 2.00e+00, 1.07e+02),
      c(1.07e+03, 2.00e+00, 1.00e+00, 1.00e+02),
      c(1.55e+03, 3.00e+00, 1.00e+00, 1.60e+01),
      c(1.95e+03, 3.00e+00, 2.00e+00, 1.60e+01),
      c(1.22e+03, 2.00e+00, 2.00e+00, 1.20e+01),
      c(1.62e+03, 3.00e+00, 1.00e+00, 1.60e+01),
      c(8.16e+02, 2.00e+00, 1.00e+00, 5.80e+01),
      c(1.35e+03, 3.00e+00, 1.00e+00, 2.10e+01),
      c(1.57e+03, 3.00e+00, 1.00e+00, 1.40e+01),
      c(1.49e+03, 3.00e+00, 1.00e+00, 5.70e+01),
      c(1.51e+03, 2.00e+00, 1.00e+00, 1.60e+01),
      c(1.10e+03, 3.00e+00, 1.00e+00, 2.70e+01),
      c(1.76e+03, 3.00e+00, 2.00e+00, 2.40e+01),
      c(1.21e+03, 2.00e+00, 1.00e+00, 1.40e+01),
      c(1.47e+03, 3.00e+00, 2.00e+00, 2.40e+01),
      c(1.77e+03, 3.00e+00, 2.00e+00, 8.40e+01),
      c(1.65e+03, 3.00e+00, 1.00e+00, 1.90e+01),
      c(1.03e+03, 3.00e+00, 1.00e+00, 6.00e+01),
      c(1.12e+03, 2.00e+00, 2.00e+00, 1.60e+01),
      c(1.15e+03, 3.00e+00, 1.00e+00, 6.20e+01),
      c(8.16e+02, 2.00e+00, 1.00e+00, 3.90e+01),
      c(1.04e+03, 3.00e+00, 1.00e+00, 2.50e+01),
      c(1.39e+03, 3.00e+00, 1.00e+00, 6.40e+01),
      c(1.60e+03, 3.00e+00, 2.00e+00, 2.90e+01),
      c(1.22e+03, 3.00e+00, 1.00e+00, 6.30e+01),
      c(1.07e+03, 2.00e+00, 1.00e+00, 1.00e+02),
      c(2.60e+03, 4.00e+00, 2.00e+00, 2.20e+01),
      c(1.43e+03, 3.00e+00, 1.00e+00, 5.90e+01),
      c(2.09e+03, 3.00e+00, 2.00e+00, 2.60e+01),
      c(1.79e+03, 4.00e+00, 2.00e+00, 4.90e+01),
      c(1.48e+03, 3.00e+00, 2.00e+00, 1.60e+01),
      c(1.04e+03, 3.00e+00, 1.00e+00, 2.50e+01),
      c(1.43e+03, 3.00e+00, 1.00e+00, 2.20e+01),
      c(1.16e+03, 3.00e+00, 1.00e+00, 5.30e+01),
      c(1.55e+03, 3.00e+00, 2.00e+00, 1.20e+01),
      c(1.98e+03, 3.00e+00, 2.00e+00, 2.20e+01),
      c(1.06e+03, 3.00e+00, 1.00e+00, 5.30e+01),
      c(1.18e+03, 2.00e+00, 1.00e+00, 9.90e+01),
      c(1.36e+03, 2.00e+00, 1.00e+00, 1.70e+01),
      c(9.60e+02, 3.00e+00, 1.00e+00, 5.10e+01),
      c(1.46e+03, 3.00e+00, 2.00e+00, 1.60e+01),
      c(1.45e+03, 3.00e+00, 2.00e+00, 2.50e+01),
      c(1.21e+03, 2.00e+00, 1.00e+00, 1.50e+01),
      c(1.55e+03, 3.00e+00, 2.00e+00, 1.60e+01),
      c(8.82e+02, 3.00e+00, 1.00e+00, 4.90e+01),
      c(2.03e+03, 4.00e+00, 2.00e+00, 4.50e+01),
      c(1.04e+03, 3.00e+00, 1.00e+00, 6.20e+01),
      c(1.62e+03, 3.00e+00, 1.00e+00, 1.60e+01),
      c(8.03e+02, 2.00e+00, 1.00e+00, 8.00e+01),
      c(1.43e+03, 3.00e+00, 2.00e+00, 2.10e+01),
      c(1.66e+03, 3.00e+00, 1.00e+00, 6.10e+01),
      c(1.54e+03, 3.00e+00, 1.00e+00, 1.60e+01),
      c(9.48e+02, 3.00e+00, 1.00e+00, 5.30e+01),
      c(1.22e+03, 2.00e+00, 2.00e+00, 1.20e+01),
      c(1.43e+03, 2.00e+00, 1.00e+00, 4.30e+01),
      c(1.66e+03, 3.00e+00, 2.00e+00, 1.90e+01),
      c(1.21e+03, 3.00e+00, 1.00e+00, 2.00e+01),
      c(1.05e+03, 2.00e+00, 1.00e+00, 6.50e+01))
y=c(300.  , 509.8 , 394.  , 540.  , 415.  , 230.  , 560.  , 294.  ,
    718.2 , 200.  , 302.  , 468.  , 374.2 , 388.  , 282.  , 311.8 ,
    401.  , 449.8 , 301.  , 502.  , 340.  , 400.28, 572.  , 264.  ,
    304.  , 298.  , 219.8 , 490.7 , 216.96, 368.2 , 280.  , 526.87,
    237.  , 562.43, 369.8 , 460.  , 374.  , 390.  , 158.  , 426.  ,
    390.  , 277.77, 216.96, 425.8 , 504.  , 329.  , 464.  , 220.  ,
    358.  , 478.  , 334.  , 426.98, 290.  , 463.  , 390.8 , 354.  ,
    350.  , 460.  , 237.  , 288.3 , 282.  , 249.  , 304.  , 332.  ,
    351.8 , 310.  , 216.96, 666.34, 330.  , 480.  , 330.3 , 348.  ,
    304.  , 384.  , 316.  , 430.4 , 450.  , 284.  , 275.  , 414.  ,
    258.  , 378.  , 350.  , 412.  , 373.  , 225.  , 390.  , 267.4 ,
    464.  , 174.  , 340.  , 430.  , 440.  , 216.  , 329.  , 388.  ,
    390.  , 356.  , 257.8 )
} # Load large housing price data set

size=NULL; rooms=NULL; floors=NULL; age=NULL
for (i in 1:length(x)) {
  size = append(size, x[[i]][1])
  rooms = append(rooms, x[[i]][2])
  floors = append(floors, x[[i]][3])
  age = append(age, x[[i]][4])
}

training_set <- data.frame(
  size = size,
  rooms = rooms,
  floors = floors,
  age = age,
  price = y
)

plot(training_set$size,training_set$price)
plot(training_set$rooms,training_set$price)
plot(training_set$floors,training_set$price)
plot(training_set$age,training_set$price)

features_set=training_set[,c("size","rooms","floors","age")]
target_set=training_set[,"price"]
x=features_set
y=target_set

w_init=c(0,0,0,0)
b_init=0
tmp_alpha=1e-7
iters=10

grad_desc = h.machinelearning.gradient_descent(x=features_set, y=target_set,
                                               w_in=w_init, b_in=b_init, 
                                               alpha=tmp_alpha, num_iters=iters, 
                                               verbose=TRUE, search=FALSE, plot=TRUE, scale=FALSE)


training_set <- data.frame(
  size = scale(size),
  rooms = scale(rooms),
  floors = scale(floors),
  age = scale(age),
  price = scale(y)
)

plot(training_set$size,training_set$age)


features_set=training_set[,c("size","rooms","floors","age")]
target_set=training_set[,"price"]
x=features_set
y=target_set
tmp_alpha=1.0e-1 
grad_desc = h.machinelearning.gradient_descent(x=features_set, y=target_set,
                                               w_in=w_init, b_in=b_init, 
                                               alpha=tmp_alpha, num_iters=1000, 
                                               verbose=TRUE, search=FALSE, plot=TRUE, scale=TRUE)





x = c(6.1101,  5.5277,  8.5186,  7.0032,  5.8598,  8.3829,  7.4764,
   8.5781,  6.4862,  5.0546,  5.7107, 14.164 ,  5.734 ,  8.4084,
   5.6407,  5.3794,  6.3654,  5.1301,  6.4296,  7.0708,  6.1891,
   20.27  ,  5.4901,  6.3261,  5.5649, 18.945 , 12.828 , 10.957 ,
   13.176 , 22.203 ,  5.2524,  6.5894,  9.2482,  5.8918,  8.2111,
   7.9334,  8.0959,  5.6063, 12.836 ,  6.3534,  5.4069,  6.8825,
   11.708 ,  5.7737,  7.8247,  7.0931,  5.0702,  5.8014, 11.7   ,
   5.5416,  7.5402,  5.3077,  7.4239,  7.6031,  6.3328,  6.3589,
   6.2742,  5.6397,  9.3102,  9.4536,  8.8254,  5.1793, 21.279 ,
   14.908 , 18.959 ,  7.2182,  8.2951, 10.236 ,  5.4994, 20.341 ,
   10.136 ,  7.3345,  6.0062,  7.2259,  5.0269,  6.5479,  7.5386,
   5.0365, 10.274 ,  5.1077,  5.7292,  5.1884,  6.3557,  9.7687,
   6.5159,  8.5172,  9.1802,  6.002 ,  5.5204,  5.0594,  5.7077,
   7.6366,  5.8707,  5.3054,  8.2934, 13.394 ,  5.4369)
y = c(17.592  ,  9.1302 , 13.662  , 11.854  ,  6.8233 , 11.886  ,
      4.3483 , 12.     ,  6.5987 ,  3.8166 ,  3.2522 , 15.505  ,
      3.1551 ,  7.2258 ,  0.71618,  3.5129 ,  5.3048 ,  0.56077,
      3.6518 ,  5.3893 ,  3.1386 , 21.767  ,  4.263  ,  5.1875 ,
      3.0825 , 22.638  , 13.501  ,  7.0467 , 14.692  , 24.147  ,
      -1.22   ,  5.9966 , 12.134  ,  1.8495 ,  6.5426 ,  4.5623 ,
      4.1164 ,  3.3928 , 10.117  ,  5.4974 ,  0.55657,  3.9115 ,
      5.3854 ,  2.4406 ,  6.7318 ,  1.0463 ,  5.1337 ,  1.844  ,
      8.0043 ,  1.0179 ,  6.7504 ,  1.8396 ,  4.2885 ,  4.9981 ,
      1.4233 , -1.4211 ,  2.4756 ,  4.6042 ,  3.9624 ,  5.4141 ,
      5.1694 , -0.74279, 17.929  , 12.054  , 17.054  ,  4.8852 ,
      5.7442 ,  7.7754 ,  1.0173 , 20.992  ,  6.6799 ,  4.0259 ,
      1.2784 ,  3.3411 , -2.6807 ,  0.29678,  3.8845 ,  5.7014 ,
      6.7526 ,  2.0576 ,  0.47953,  0.20421,  0.67861,  7.5435 ,
      5.3436 ,  4.2415 ,  6.7981 ,  0.92695,  0.152  ,  2.8214 ,
      1.8451 ,  4.2959 ,  7.2029 ,  1.9869 ,  0.14454,  9.0551 ,
      0.61705)
plot(x,y)

grad_desc = h.machinelearning.gradient_descent(x=x, y=y,
                                               w_in=100, b_in=100, 
                                               alpha=.001, num_iters=5000, 
                                               verbose=TRUE, search=TRUE, plot=TRUE, scale=FALSE)



