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
{
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
                                                   w_init=0, b_init=0, 
                                                   alpha=1.0e-2, n_iters=1e4, 
                                                   verbose=T, search=F, plot=T)
    
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
}

##################
##### Week 2 #####
##################
{
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
    w_inits = c(0.39133535, 18.75376741, -53.36032453, -26.42131618) # One for each feature, in order of columns in the features set
    #w_inits = runif(length(features_set))
    b_init = .1
    w_inits = c(.1,.1,.1,.1)
    
    # Run gradient descent
    grad_desc = h.machinelearning.gradient_descent(x=features_set, y=target_set,
                                                   w_init=w_inits, b_init=b_init, 
                                                   alpha=5e-7, n_iters=1000,
                                                   verbose=TRUE, search=FALSE, plot=TRUE, scale=FALSE,
                                                   max_iters=1e+10, conv_crit=1e-50, div_crit=1e+100,
                                                   searchspace=3, cost.method = NULL, gradient.method = NULL)
    
  }
  
  
  "Polynomial Features

Out of the box, linear regression provides a means of building models of the form:
𝑓(𝐰,𝑏) = 𝑤1𝑥1 + 𝑤2𝑥2 +...+ 𝑤𝑛𝑥𝑛 +𝑏

But what if your features/data are non-linear or are combinations of features? For example, housing prices do not
tend to be linear with living area, but penalize very small or very large houses. How can we use the machinery of 
linear regression to fit this curve? Recall, the 'machinery' we have is the ability to modify the parameters 𝐰 and 𝐛
to fit the model (i.e. the equation) to the training data. However, no amount of adjusting those parameters will achieve
a good fit to a non-linear curve.

Let's try fitting a non-linear curve. We'll start with a simple quadratic:  y = 1 + x^2
"
  
  # Let's engineer some polynomial features and see how our model handles those.
  {
    
    training_set <- data.frame(
      x = 0:19
    )
    training_set$y = (training_set$x^2) + 1 # Thus, the true model of the data is 'y = 1 * x^2 + 1'
    
    # As you can see, the values we're trying to model are quadratic
    plot(training_set$x, training_set$y)
    
    # Descent parameters
    n_iters = 1e4
    alpha = 1.0e-2
    w_initial = 0
    b_initial = 0
    
    # Run gradient descent
    grad_desc = h.machinelearning.gradient_descent(x=training_set$x, y=training_set$y,
                                                   w_init=w_initial, b_init=b_initial, 
                                                   alpha=alpha, n_iters=n_iters, 
                                                   verbose=T, search=F, plot=F)
    # Look at the results! Our best fitting w and b combo is:
    "Gradient descent found w (19) and b (-56), with a cost (J_wb) of 438.9."
    w_fit = grad_desc[[1]]
    b_fit = grad_desc[[2]]
    
    # We can now plot this fitted line over the data itself.
    plot(training_set$x, training_set$y) # original data
    lines(training_set$x, w_fit*training_set$x + b_fit, col="red") # our fitted line, simply plugging our fitted w,b into the model template:  y = w*x + b
    # Okay, that doesn't look to great unfortunately and our model's cost (of ~440) isn't great either.
    
    # Because we know that our y-values are polynomial, we can engineer a feature that has
    # better chances of predicting that. So, let's square our x-values this time and see what happens.
    training_set$x_sq = training_set$x^2
    
    # Run gradient descent on our new squared feature
    alpha = 1e-5
    grad_desc = h.machinelearning.gradient_descent(x=training_set$x_sq, y=training_set$y,
                                                   w_init=w_initial, b_init=b_initial,
                                                   alpha=alpha, n_iters=n_iters, 
                                                   verbose=T, search=F, plot=F)
    # This should have done a much better job.
    w_fit = grad_desc[[1]] # w (Should be ~1.004)
    b_fit = grad_desc[[2]] # b (Should be ~0.049)
    grad_desc[[3]][n_iters,1] # J_wb (Should be ~0.207)
    # In other words, gradient descent is telling us that our model of the data should be following equation:
    # y = 1.004 * x^2 + 0.049   <- This is our computational model
    # Does that make sense? Recall, our model "template" is :   y = w * x + b
    # So, just plug in the values we found with gradient descent. Don't forget that we squared our x-feature, so
    # you need to square it in the final equation of the model too.
    lines(w_fit*training_set$x^2 + b_fit, col="blue")
    # Alright! Our algorithm found the equation 'y = 1.004 * x^2 + .049', which is pretty close to the 'true' model
    # of 'y = 1 * x^2 + 1. Not perfect, but certainly looks better than the first model fit.
    
    
    ### Selecting Features
    # In our examples above, we already knew a priori that a squared term (i.e., x^2) was necessary. But, with
    # real data, that isn't always going to be as obvious (even after visualizing data -- which you should always
    # do regardless). One could try multiple features to find the most useful. For instance, what if
    # we had a first degree term (i.e., x), a second degree term (x^2), and a third degree term (x^3)?
    # Thus, the model we'd be fitting would be: 𝑦=𝑤1𝑥1 + 𝑤2𝑥2^2 + 𝑤3𝑥3^3  + 𝑏?
    
    training_set <- data.frame(
      x1 = 0:19,
      x2 = (0:19)^2,
      x3 = (0:19)^3
    )
    training_set$y = 1 * training_set$x1^2 + 0 # True model
    
    # What are the features and what's being predicted?
    features_set=training_set[,c("x1","x2","x3")]
    target_set=training_set[,"y"]
    
    # Parameters
    n_iters=1e4
    alpha = 1e-7
    b_initial = 0
    w_initial = c(0,0,0) # Recall, we need an initial w for each feature (in this case, 3 of them)
    #w_initial = c(0.1,0.9,0.1)
    
    grad_desc = h.machinelearning.gradient_descent(x=features_set, y=target_set,
                                                   w_init=w_initial, b_init=b_initial,
                                                   alpha=alpha, n_iters=n_iters, 
                                                   verbose=T, search=F, plot=F)
    "Gradient descent found the following values:
      w_1:  0.0823753
      w_2:  0.535521
      w_3:  0.0275222
      b:    0.0105619
      J:    79.0693"
    w1_fit = grad_desc[[1]][1]
    w2_fit = grad_desc[[1]][2]
    w3_fit = grad_desc[[1]][3]
    b_fit = grad_desc[[2]]
    plot(1:20, training_set$y) # our data
    lines(1:20, w1_fit*training_set$x1 + w2_fit*training_set$x2 + w3_fit*training_set$x3 + b_fit, col="blue") # our model
    # Gradient descent is telling us that our final model should be 'y = .08*x + .54*x^2 + .03*x^3 + .0106'
    # In other words, the algorithm is pushing for an emphasis on the second degree term (x^2) and minimizing
    # the impact of the first and third terms. If we ran the algorithm longer, it would further diminish the
    # impact of those two and further emphasize the squared term (by moving the .54 closer and closer to 1.0).
    # The important point, therefore, is that gradient descent can learn to pick the 'correct' features for us.
    # In short, the importance of a feature can be inferred from the weight given to it in the fitted model.
    
    # Since our three features vary in range, let's see what happens when we scale them (by adding scale=TRUE).
    grad_desc = h.machinelearning.gradient_descent(x=features_set, y=target_set,
                                                   w_init=w_initial, b_init=b_initial,
                                                   alpha=alpha, n_iters=n_iters, 
                                                   verbose=T, search=F, plot=F,
                                                   scale=T)
    "Gradient descent found the following values:
      w_1:  0.106637
      w_2:  0.110469
      w_3:  0.108901
      b:    0.123438
      J:    14016"
    w1_fit = grad_desc[[1]][1]; w2_fit = grad_desc[[1]][2]; w3_fit = grad_desc[[1]][3]; b_fit = grad_desc[[2]]
    lines(1:20, w1_fit*scale(training_set$x1) + w2_fit*scale(training_set$x2) + w3_fit*scale(training_set$x3) + b_fit, col="red")
    # What happened here?
    
    # What if we change the learning rate to be a bigger number?
    alpha=1e-1; n_iters=1e5
    #alpha=.6; n_iters=1e3
    grad_desc = h.machinelearning.gradient_descent(x=features_set, y=target_set,
                                                   w_init=w_initial, b_init=b_initial,
                                                   alpha=alpha, n_iters=n_iters, 
                                                   verbose=T, search=F, plot=F,
                                                   scale=T)
    "Gradient descent found the following values:
      w_1:  0.000104596
      w_2:  116.443
      w_3:  0.000167483
      b:    123.5
      J:    6.90045e-11"
    w1_fit = grad_desc[[1]][1]; w2_fit = grad_desc[[1]][2]; w3_fit = grad_desc[[1]][3]; b_fit = grad_desc[[2]]
    lines(w1_fit*scale(training_set$x1) + w2_fit*scale(training_set$x2) + w3_fit*scale(training_set$x3)  + b_fit, col="green")
    # Okay, our cost has dropped enormously and both w_1 and w_3 have effectively been eliminated in their importance.
    # This shows us that we really need to pay attention to all the elements of our data and the models we run on them:
    # the learning rate, the final cost, the number of iterations, whether our data are scaled or not, and so on.
    
    # Depending on whether you're running these algorithms in R, python, or Matlab, or Windows versus Mac,
    # you might have noticed that it can take a long time to run. Try playing more with the alpha and n_iters
    # while paying attention to how quickly the cost (J) drops. The algorithm will, of course, finish faster if
    # you ask for fewer n_iters. But you can also try a higher learning rate, so see if the algorithm can bring
    # the cost down at a faster pace (e.g., from 1e-1, you can try .3, .6. or .9). At some point, you will
    # end up with a learning rate that's too high and the algorithm will diverge (instead of converge), meaning
    # that the cost kept going up at each iteration instead of going down. Scaling can also speed up processing.
    # Playing around with these (and other) parameters for a while can already give you strong suggestions as to
    # which features are helpful and which aren't (at least, relative to each other). That way, you may be able to
    # speed up your work by quickly noticing which features can be dropped and then reducing the work that gradient
    # descent has to do with a smaller set of features.
    
    
    
    ### Complex function
    # With feature engineering, we can even model quite complex functions.
    xvals = 0:19
    training_set <- data.frame(
      x1 = xvals,
      x2 = xvals^2,
      x3 = xvals^3,
      x4 = xvals^4,
      x5 = xvals^5,
      x6 = xvals^6,
      x7 = xvals^7,
      x8 = xvals^8,
      x9 = xvals^9,
      x10 = xvals^10,
      x11 = xvals^11,
      x12 = xvals^12,
      x13 = xvals^13
    )
    training_set$y = cos(xvals/2) # True model
    
    # What are the features and what's being predicted?
    features_set=training_set[,1:13] # The first 13 columns, which are our 13 features
    target_set=training_set[,"y"]
    plot(xvals, training_set$y) # our data
    
    # Parameters
    n_iters=1e5
    alpha = .1
    b_initial = 0
    w_initial = c(0,0,0,0,0,0,0,0,0,0,0,0,0)
    
    # This will take a long time in R. Very long!
    grad_desc = h.machinelearning.gradient_descent(x=features_set, y=target_set,
                                                   w_init=w_initial, b_init=b_initial,
                                                   alpha=alpha, n_iters=n_iters, 
                                                   verbose=T, search=F, plot=F,
                                                   scale=T)
    "Gradient descent found the following values:
      w_1:  -3.81903
      w_2:  3.76706
      w_3:  6.43265
      w_4:  0.298999
      w_5:  -4.22198
      w_6:  -4.99168
      w_7:  -3.21612
      w_8:  -0.570534
      w_9:  1.69477
      w_10:  2.85959
      w_11:  2.61382
      w_12:  0.91283
      w_13:  -2.13653
      b:    -0.00728719
      J:    0.0172613"
    lines(xvals, -3.81903*scale(training_set$x1) + 3.76706*scale(training_set$x2) + 6.43265*scale(training_set$x3) +
                 0.298999*scale(training_set$x4) + -4.22198*scale(training_set$x5) + -4.99168*scale(training_set$x6) + 
                 -3.21612*scale(training_set$x7) + -0.570534*scale(training_set$x8) + 1.69477*scale(training_set$x9) + 
                 2.85959*scale(training_set$x10) + 2.61382*scale(training_set$x11) + 0.91283*scale(training_set$x12) +
                -2.13653*scale(training_set$x13) + -0.00728719, col="blue")
    # Not bad!
    
    # There are ways to speed up the process in R, but that matter lies beyond the scope of these sessions.
    # R's focus is on running algebraic computations in a stable and accurate manner, rather than with a
    # trade-off in accuracy to speed things up. For many purposes, this trade-off in accuracy is probably
    # negligible, so feel free to look into R Open, GotoBlas, OpenBLAS, MKL, and related software that can
    # speed up algebraic processes in R. Do note that in, e.g., Python or Matlab, these processes are already
    # much faster. So, you may consider running some gradient descent algorithms in Python or Matlab instead.
    
  }
  
}

##################
##### Week 3 #####
##################
{
  "Classification"
  {
    # In addition to predicting values using a model, another common machine learning problem is that of
    # classifying data points into categories. For classification problems, our linear regression approach
    # is - sadly - not very suited. Instead, we'll need something called logistic regression.
    
    # Let's see how they both perform by setting up some data (x, y) that has groupings (g) and a color map (c)
    training_set1a <- data.frame(
      x = c(0, 1, 2, 3, 4, 5),
      g = c(0, 0, 0, 1, 1, 1),
      c = c("blue","blue","blue","red","red","red")
    )
    plot(training_set1a$x, training_set1a$g, col=training_set1a$c, pch=16) # One variable (x)
    
    # Run gradient descent
    grad_desc = h.machinelearning.gradient_descent(x=training_set1a$x, y=training_set1a$g,
                                                   w_init=0, b_init=0, 
                                                   alpha=.1, n_iters=1e6, 
                                                   verbose=T, search=F, plot=F)
    "Gradient descent found w (0.257143) and b (-0.142857), with a cost (J_wb) of 0.0285714."
    lines(training_set1a$x, 0.257143*training_set1a$x + -0.142857, lwd=4) # Fitted linear regression line
    # Note, lines() takes an x and y argument, so we first pass it the independent x values and then the model
    
    # We can then draw a threshold, say at g = .5
    abline(h=.5) # abline will now draw a horizontal line at .5 (hence, h=.5)
    # Next, for example, we decide to classify points that lie to one side of the intersection
    # between our threshold and our regression line will be classified as one group and points
    # on the other side as part of the other group.
    f1 <- function(x) 0.257143*x + -0.142857 # Our fitted model
    f2 <- function() .5                      # Our arbitrary threshold
    
    uniroot(function(x) f1(x)-f2(),c(0,5))$root # This asks for the intersection (root) between
    # two functions in the interval 0 to 5. And indeed, it finds x is 2.499998 -- that looks right.
    abline(v=uniroot(function(x) f1(x)-f2(),c(0,5))$root) # abline will now draw a vertical line at ~2.5
    # In this case, all data are classified correctly (blue on the left side of ~2.5, red on the right side).
    
    # But now let's add another data point
    training_set1b <- data.frame(
      x = c(0, 1, 2, 3, 4, 5, 10),
      g = c(0, 0, 0, 1, 1, 1, 1),
      c = c("blue","blue","blue","red","red","red","red")
    )
    plot(training_set1b$x, training_set1b$g, col=training_set1b$c, pch=16) # One variable (x)
    
    # Run gradient descent
    grad_desc = h.machinelearning.gradient_descent(x=training_set1b$x, y=training_set1b$g,
                                                   w_init=0, b_init=0, 
                                                   alpha=.01, n_iters=1e6, 
                                                   verbose=T, search=F, plot=F)
    "Gradient descent found w (0.117391) and b (0.152174), with a cost (J_wb) of 0.057764."
    
    # Plot the model
    lines(training_set1b$x, 0.117391*training_set1b$x + -0.152174, lwd=4) # Fitted linear regression line
    abline(h=.5) # Threshold
    f1 <- function(x) 0.117391*x + -0.152174 # Our fitted model
    f2 <- function() .5                      # Our arbitrary threshold
    abline(v=uniroot(function(x) f1(x)-f2(),c(0,10))$root) # abline will now draw a vertical line at ~5.56
    
    # Okay, that looks bad! We're now misclassifying the original 3 red data points and only capturing the
    # new data point that we just introduced. Clearly, we need a better approach than linear regression.
  }
    
  "Logistic regression"
  {
    # Enter logistic regression, one of the (if not actually the) most used classification algorithm in the world!
    # Instead of a linear model 'y(x) = w*x + b', we will use a so-called sigmoid function (a.k.a. logistic
    # function). The sigmoid function is 'g(z) = 1 / (1+e^-z)', where 0 < g(z) < 1 (so, y can only take on values
    # between 0 and 1). Importantly, the z variable in the sigmoid function is actually y(x) from earlier, the 
    # output of the linear regression function.
    
    # Logistic regression also uses a cost function that is different from linear regression:
    #
    # loss(f_wb_i, y_i) =  { -log(f_wb_i)      if y_i = 1 
    #                        -log(1 - f_wb_i)  if y_i = 0 }
    #
    # The error at the level of each individual training sample is often referred to as 'loss', but it is
    # essentially another word for 'error'. Thus, the average loss -- across all training examples -- is
    # the cost of the model. It's referred to as loss, because we're "losing (some) information" by modeling.
    # Importantly, as you can see, there are two elements to this function. It uses two 'loss curves', one
    # for each category in the data (y_i = 1 or y_i = 0).
    
    # The loss function above can be rewritten as follows and that's how it's implemented in h.lab:
    #
    # loss(f_wb_i, y_i) = (-y_i * log(f_wb_i)) - (1 - y_i) * log(1 - f_wb_i)
    #
    # Note that y_i is the category value, so, for instance, it would only take on a value of either 1 or 0.
    # If y_i is 1, you'll see that it indeed reduces to top loss function in the first example. When it's 0,
    # it indeed reduces to the bottom one.
    #
    # You might also be wondering why the cost function looks the way it does. It is derived from a statistical
    # concept called 'maximum likelihood'. No need to worry about what that is just yet. We'll return to it later.
    # For now, let's see how it works.
    
    
    # Let's see how logistic regression handles the following data
    training_set2 <- data.frame(
      x1 = c(.5, 1, 1.5, 3, 2, 1),
      x2 = c(1.5, 1, .5, .5, 2, 2.5),
      g = c(0, 0, 0, 1, 1, 1),
      c = c("blue","blue","blue","red","red","red")
    )
    plot(training_set2$x1, training_set2$x2, col=training_set2$c, pch=16, xlim=c(0, 3), ylim=c(0, 3.5))  # Two variables (x1, x2)
    
    features_set=training_set2[,c("x1","x2")]
    target_set=training_set2[,"g"]
    
    # Run logistic gradient descent
    grad_desc = h.machinelearning.gradient_descent(x=features_set, y=target_set,
                                                   w_init=c(0,0), b_init=0,
                                                   alpha=.1, n_iters=1e4, 
                                                   verbose=T, search=F, plot=F,
                                                   model="logistic") # <- Note that we need to request a logistic model
    "Gradient descent found the following values:
      w_1:  5.28123
      w_2:  5.07816
      b:    -14.2224
      J:    0.0171178"
    
    # Plot the decision boundary
    x1_intercept = -(-14.2224) / 5.28123 # Decision boundary intercepts the x1 variable's axis at this point
    x2_intercept = -(-14.2224) / 5.07816 # ..and the x2 variable's axis at this point.
    # Thus, the rule is:  - b / w  for each w we have. Note that we are flipping the sign of b.
    
    # We can use xspline to draw the boundary. This draws a line between provided coordinates, but its approach
    # may feel a bit odd at first. We first provide all the x-coordinates, then all the y-coordinates. In this case,
    # xspline draws the line from the first provided x-coordinate (x1_intercept) to the first provided y-coordinate
    # which is the SECOND 0 below. The line is then drawn to the coordinate based on the second provided x-coordinate,
    # which is our FIRST 0 below and the second provided y-coordinate, which is x2_intercept.
    xspline(c(x1_intercept,0),c(0,x2_intercept),lwd=4)
    #segments(x0=0, y0=x2_intercept, x1=x1_intercept, y1=0,lwd=4) # These are different ways to draw the line
    #lines(x=c(0,x1_intercept),y=c(x2_intercept,0),lwd=4)
    
    # Great! Now, when a data point falls on one side of the decision boundary, we categorize it as one category, and
    # when it falls on the other side, we categorize it as the other.
    # In other words, what we now need /mathematically/ (i.e. not just visually in our plot) is some way to take
    # our data (and future data points) and convert them into the categories, by computing on which side of the
    # decision boundary they fall. Let's start by looking at our model again:
    5.28123*training_set2$x1 + 5.07816*training_set2$x2 + -14.2224
    # When you run this line, you get the following output:
    "-3.964545 -3.863010 -3.761475  4.160370  6.496380  3.754230"
    # Recall we had 6 data points and something should immediately become clear: the first 3 data points are
    # negative in sign, the last 3 are positive. Think of the decision boundary as a zero point, with some data
    # falling "above" it and others "below" it. Or "right" versus "left" of the line. Or.. "positive" versus
    # "negative" /relative to/ the line. Thus, for any given data point, we compute following, for example with
    # the first data point in the data frame ( x1=.5, x2=1.5 ):
    "5.28123*.5 + 5.07816*1.5 + -14.2224 = -3.964545"
    # And since this is negative, it falls in the range of the blue group. Easy-peasy!
    # And the previous line of code executed that computation for the entire data frame. If we want, we can then
    # convert the numeric results into categories ourselves as follows:
    classification_set <- data.frame(
      true_labels = training_set2$g,   # Ground truth of our labels
      relative_to_bound = 5.28123*training_set2$x1 + 5.07816*training_set2$x2 + -14.2224 # Data relative to decision boundary
    )
    classification_set$modeled_labels[classification_set$relative_to_bound < 0] = 0
    classification_set$modeled_labels[classification_set$relative_to_bound > 0] = 1
    # Naturally, we could now compare what percentage of our modeled_labels match the ground truth,
    # but that would be pointless -- they will match 100%. In future steps, we will first set up
    # a decision boundary with /some/ portion of our data (the 'training set') and then test that
    # boundary on data for which we do know the true labels but which weren't used to train the model
    # (the 'test set'). There's no point in testing a model on the same data it was trained on.
    
    # Also look at the following. We can plot the model with our found weights and bias.
    plot(seq(.5,3,by=.5), (1 / (1 + exp(- (5.28123*training_set2$x1 + 5.07816*training_set2$x2 + -14.2224)))), col=training_set2$c, pch=16)
    # As you can see, it has labeled the first three points as belonging to category 0 and the last three as belonging to category 1.
    
    
    # We'll return to our data from the classification section, but now run a logistic model
    training_set1a <- data.frame(
      x = c(0, 1, 2, 3, 4, 5),
      g = c(0, 0, 0, 1, 1, 1),
      c = c("blue","blue","blue","red","red","red")
    )
    plot(training_set1a$x, training_set1a$g, col=training_set1a$c, pch=16) # One variable (x)
    features_set=training_set1a[,c("x")]
    target_set=training_set1a[,"g"]
    # Run gradient descent
    grad_desc = h.machinelearning.gradient_descent(x=features_set, y=target_set,
                                                   w_init=0, b_init=0, 
                                                   alpha=.0001, n_iters=1e6, 
                                                   verbose=T, search=F, plot=F,
                                                   model = "logistic")    
    "Gradient descent (logistic) found w (2.18657) and b (-5.16871), with a cost (J_wb) of 0.113415."
    
    # Plot the fitted sigmoid
    # Logistic model: 1 / (1 + exp(- (w*x + b) ))
    fitted_sigmoid = 1 / (1 + exp(- (2.18657*training_set1a$x + -5.16871) ))
    "0.005659692 0.048239281 0.310973935 0.800751839 0.972816095 0.996871724"
    lines(training_set1a$x, fitted_sigmoid, pch = 18, col = "black", type = "b")
    lines(training_set1a$x, 2.18657*training_set1a$x + -5.16871, lwd=4) # We can even use the found weights on a linear line  
    abline(h=.5)   
    
    # And also the other data set
    training_set1b <- data.frame(
      x = c(0, 1, 2, 3, 4, 5, 10),
      g = c(0, 0, 0, 1, 1, 1, 1),
      c = c("blue","blue","blue","red","red","red","red")
    )
    plot(training_set1b$x, training_set1b$g, col=training_set1b$c, pch=16) # One variable (x)
    features_set=training_set1b[,c("x")]
    target_set=training_set1b[,"g"]
    # Run gradient descent
    grad_desc = h.machinelearning.gradient_descent(x=features_set, y=target_set,
                                                   w_init=0, b_init=0, 
                                                   alpha=.0001, n_iters=1e6, 
                                                   verbose=T, search=F, plot=F,
                                                   model = "logistic")
    "Gradient descent (logistic) found w (2.05927) and b (-4.83612), with a cost (J_wb) of 0.104991."
    
    # Plot the fitted sigmoid
    # Logistic model: 1 / (1 + exp(- (w*x + b) ))
    fitted_sigmoid = 1 / (1 + exp(- (2.05927*training_set1b$x + -4.83612) ))
    "0.007875281 0.058588054 0.327926106 0.792767724 0.967734524 0.995765430 0.999999856"
    lines(training_set1b$x, fitted_sigmoid, pch = 18, col = "black", type = "b")
    lines(training_set1b$x, 2.05927*training_set1b$x + -4.83612, lwd=4) # We can even use the found weights on a linear line  
    abline(h=.5)
    
    # Thus, each sample in the data set has a corresponding sigmoid value. For instance, the third blue sample
    # has a value of 0.327926106 for the fitted sigmoid. The classification is now based on the values
    # for the fitted sigmoid: any below 0.5 belong to one group, any above 0.5 to the other group. And indeed,
    # the logistic regression has circumvented the issues we ran into with the linear regression. Now, all
    # the samples are accurately classified.
    
    
    ### Non-linear decision boundaries
    # At this point, you might be wondering "But what if my data aren't neatly clustered like this?
    # What if it's immediately obvious (from just visualizing my data) that a straight line is clearly
    # not gonna work?" Well, then I've got good news for you! Take a look at this:
    training_set3 <- data.frame(
      x1 = c(.2, 1, 1.5, .2, 1, 1.5, -.2, -1, -1.5, -.2, -1, -1.5,
             .2, .5, .6, .7, .4, .3, -.2, -.4, -.5, -.6, -.7, -.6),
      x2 = c(1.8, 1.5, .8, -1.6, -1.5, -1.8, -1.8, -1.5, -1.1, 1.8, 1.5, 1.8,
             .3, .4, .1, -.2, -.4, -.3, .3, .5, .7, -.3, -.4, -.6),
      g = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
      c = c("red","red","red","red","red","red","red","red","red","red","red","red",
            "blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue")
    )
    plot(training_set3$x1, training_set3$x2, col=training_set3$c, pch=16, xlim=c(-2, 2), ylim=c(-2, 2)) # Two variables (x1, x2)
    # Okay, it looks like this isn't easily 'solvable' with a simple decision boundary (i.e., a single straight line)
    # Recall that, earlier, we used feature engineering to add some polynomial features and then let the algorithm find
    # weights that select the features for us. Let's try that here.
    training_set3$x1_2 = training_set3$x1^2 # Squares
    training_set3$x2_2 = training_set3$x2^2
    training_set3$x1_3 = training_set3$x1^3 # Cubes
    training_set3$x2_3 = training_set3$x2^3
    
    features_set=training_set3[,c("x1","x2","x1_2","x1_3","x2_2","x2_3")]
    target_set=training_set3[,"g"]
    
    grad_desc = h.machinelearning.gradient_descent(x=features_set, y=target_set,
                                                   w_init=c(0,0,0,0,0,0), b_init=0,
                                                   alpha=.1, n_iters=1e4, 
                                                   verbose=T, search=F, plot=F,
                                                   model="logistic")
    "Gradient descent (logistic) found the following values:
      w_1:  -0.627192
      w_2:  0.0364092
      w_3:  -3.44463
      w_4:  -0.174596
      w_5:  -4.87699
      w_6:  0.160551
      b:    7.29609
      J:    0.00316245"
    # From the looks of it, w_3 and w_5 appear to be useful features (the others are, relatively speaking, much closer to 0). These two
    # weights correspond to our two squared terms: x1^2 and x2^2.
    #     features_set=training_set3[,c("x1","x2","x1_2","x1_3","x2_2","x2_3")] <- weights 3 and 5 are the squared terms
    
    # So now let's replot the data, but now using the squared terms instead.
    plot(training_set3$x1^2, training_set3$x2^2, col=training_set3$c, pch=16)
    x1_intercept = -(7.29609) / -3.44463; x2_intercept = -(7.29609) / -4.87699
    xspline(c(x1_intercept,0),c(0,x2_intercept),lwd=4)
    
    # As you can see, a little bit of feature engineering can sometimes reduce a complex problem to a simpler problem.
  }
  
  "Overfitting"
  # You may now be tempted to think "Gee, why not always engineer a ton of features and then fit my data -- easy!"
  # Well, with more and more (higher-order, polynomial) features, the algorithm will contort itself to find a highly
  # complicated decision boundary. At that point, you are 'overfitting' the data. This is often also called 'high variance'
  # in a model (because predictions on new data can highly vary). The reverse is also possible: 'underfitting' (a.k.a. bias).
  # When our model doesn't have an appropriate fit, it will not generalize to new data (it will be biased towards certain data).
  
  # How can we diagnose a wrong fit and what can we do to mitigate it? First, more data -- more data will always help. However,
  # this may not always be possible. Second, reducing the number of features (i.e., feature selection). However, it is possible,
  # in theory and even in practice, that all feature are informative to a substantial extent. At that point, reducing features
  # would mean throwing away valuable information. Third, we can use regularization. Regularization minimizes the extent to which
  # a given feature can impact the overall model fitting process. We're gonna look at this last option next.
  
  ### Regularization
  # Regularization is performed on the features (i.e., w_1, w_2, w_3, ... w_n). It's possible to regularize the b parameter too,
  # although that should have little impact in practice. And so, by convention, only the coefficient weights are regularized.
  # Note that regularization is meant to reduce overfitting (as opposed to underfitting).
  
  # An example of over/underfitting with linear regression
  training_set4 <- data.frame(
    x = 0:25
  )
  for (i in 1:length(training_set4$x)) { training_set4$y[i] = (training_set4$x[i]+runif(1)*3)^2 }
  # True model is quadratic with some noise
  plot(training_set4$x, training_set4$y, pch=16, xlim=c(0, length(training_set4$x)), ylim=c(0, 1100))
  
  # Set up squared terms, cubed terms, and so on up until a 6th degree polynomial
  training_set4$x_2 = training_set4$x^2 # Squares
  training_set4$x_3 = training_set4$x^3 # Cubes
  training_set4$x_4 = training_set4$x^4
  training_set4$x_5 = training_set4$x^5
  training_set4$x_6 = training_set4$x^6
  
  # Select however many of the features you wanna use
  features_set=training_set4[,c("x","x_2","x_3","x_4","x_5","x_6")]
  features_set=training_set4[,c("x","x_2","x_3")]
  features_set=training_set4[,c("x","x_2")]
  w_inits=0 # If 1 feature
  w_inits=runif(length(features_set)) # If >1 features
  
  # Fit
  grad_desc = h.machinelearning.gradient_descent(x=features_set, y=training_set4[,"y"],
                                                 w_init=w_inits, b_init=0, 
                                                 alpha=1e-100, n_iters=1e4, 
                                                 div_crit = 1e+2000, scale=F,
                                                 verbose=T, search=F, plot=F,
                                                 model = "linear",
                                                 cost.lambda_w = 0)
  
  lines(training_set4$x, 1.25765*training_set4$x_2 + 0.0956655, lwd=4)
  lines(training_set4$x, 3.65684*training_set4$x + 3.37048*training_set4$x_2 + 3.66208, lwd=4)
  lines(training_set4$x, 0.515225*training_set4$x + 0.285974*training_set4$x_2 + 0.0360659*training_set4$x_3 + 2.49354e-05, lwd=4)
  lines(training_set4$x, 0.898841*training_set4$x + 0.0671861*training_set4$x_2 + 0.671709*training_set4$x_3 + 0.658829*training_set4$x_4 + 0.920746*training_set4$x_5 + 0.615722*training_set4$x_6 + -2.53318e-89, lwd=4)
  
  
  # Now an example of over/underfitting with logistic regression
  training_set5 <- data.frame(
    x = c(c(-.75,-72,-.5,-.39,-.4,-.39,-.28,-.27,-.25,-.14,.6,.8,-.10,-.8,-.7,-.2,-.29,-.25,.45,.55,.23,.21,.08,.1,.12,.05,.17,.12,.22,.4),c(-.66,-.63,-.68,-.59,-.53,-.42,-.30,-.12,-.27,-.07,-.01,-.13,-.10,-.33,-.28,.25,.29,.49,.59,.67,.61,.99)),
    y = c(c(.2,.04,-.16,-.30,-.78,-.87,-.22,-.46,-.91,-.97,0,-.1,-.95,-.59,-.33,-.81,-.15,-.36,-.87,-.83,-.94,-.83,-.46,-.66,-.31,-.01,-.004,-.002,-.11,-.11),c(.76,.56,.40,.24,-.74,.03,.10,.54,.70,.90,-.23,-.15,.90,-.23,-.04,.06,.32,.31,.16,.47,.21,.93)),
    g = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    c = c("blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue",
          "red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red")
  )
  plot(training_set5$x, training_set5$y, col=training_set5$c, pch=16, xlim=c(-1, 1), ylim=c(-1, 1))
  
  # We'll try fitting a decision boundary using some nth degree polynomial and logistic regression
  features_set = training_set5$x
  grad_desc = h.machinelearning.gradient_descent(x=features_set, y=training_set5[,"g"],
                                                 w_init=0, b_init=0, 
                                                 alpha=1e-5, n_iters=1e5, 
                                                 div_crit = 1e+2000, scale=F,
                                                 verbose=T, search=F, plot=F,
                                                 model = "logistic",
                                                 cost.lambda_w = 0)
  "Gradient descent (logistic) found w (0.0632092) and b (-0.0598015), with a cost (J_wb) of 0.676335."
  plot(training_set5$x, 0.0632092*training_set5$x, col=training_set5$c, pch=16) # Does a decent job of separating data
  lines(training_set5$x, (0.0632092*training_set5$x + -0.0598015), pch = 18, col = "black", type = "b")

  training_set5$class = 0.0632092*training_set5$x;
  training_set5$class_col[training_set5$class < 0] = "blue"; training_set5$class_col[training_set5$class > 0] = "red";
  training_set5$class_g[training_set5$class < 0] = 0; training_set5$class_g[training_set5$class > 0] = 1;
  plot(training_set5$x, training_set5$y, col=training_set5$class_col, pch=16, xlim=c(-1, 1), ylim=c(-1, 1))
  
  
  #fitted_curve = 1 / (1 + exp(- (0.0632092*training_set5$x + -0.0598015) ))
  #lines(training_set5$x, fitted_curve, pch = 18, col = "green", type = "b")
  
}


###################
###### TESTS ######
###################
{
  # h.lab
  hlab_dir='~/Desktop/Research/+h.lab/'
  source(paste0(hlab_dir,'h.format.expand_grid.R'))
  source(paste0(hlab_dir,'h.machinelearning.cost.R'))
  source(paste0(hlab_dir,'h.machinelearning.gradient.R'))
  source(paste0(hlab_dir,'h.machinelearning.gradient_descent.R'))
  
  x = data.frame(
    a = c(1,2,3),
    a2 = c(4,5,6),
    b = c(7,8,9),
    b2 = c(10,11,12),
    c = c(13,14,15),
    c2 = c(16,17,18)
  )
  y=c(0,1,0,1,0)
  w=c(-.40165317,-0.07889237,.45788953,.03316528,.19187711,-.18448437)
  b=.5
  cost.lambda_w=.7
  cost.lambda_b=0
  model = "linear"
  
  
  # Test the cost and gradient scripts. Do you get the expected output?
  
  h.machinelearning.cost(x, y, w, b, model="linear", lambda_w=cost.lambda_w, lambda_b=cost.lambda_b)
  "3.406914"
  
  h.machinelearning.cost(x, y, w, b, model="logistic", lambda_w=cost.lambda_w, lambda_b=cost.lambda_b)
  "2.027153"
  
  h.machinelearning.gradient(x, y, w, b, method = NULL, model = "linear")
  " [[1]]
         [,1]     [,2]     [,3]     [,4]     [,5]     [,6]
    [1,] 5.105691 12.74632 20.38696 28.02759 35.66823 43.30886

    [[2]]
    [1] 2.546878
  "
  h.machinelearning.gradient(x, y, w, b, method = NULL, model = "logistic")
  " [[1]]
         [,1]     [,2]     [,3]     [,4]   [,5]     [,6]
    [1,] 1.227643 3.068207 4.908771 6.749336 8.5899 10.43046

    [[2]]
    [1] 0.6135214
  "
  
}

###################
##### FIGURES #####
###################

# Figure 2
{
  # Features
  training_set <- data.frame(
    x1 = 0:19,
    x2 = (0:19)^2,
    x3 = (0:19)^3
  )
  
  # Plot of data
  training_set$y = 1 * training_set$x1^2 + 0 # True model
  plot(1:20, training_set$y)
  
  # Unscaled
  w1_fit= 0.0823753;
  w2_fit= 0.535521;
  w3_fit= 0.0275222;
  b_fit= 0.0105619;
  lines(w1_fit*training_set$x1 + w2_fit*training_set$x2 + w3_fit*training_set$x3 + b_fit, col="blue")
  
  # Scaled, but same learning rate
  w1_fit= 0.106637;
  w2_fit= 0.110469;
  w3_fit= 0.108901;
  b_fit= 0.123438;
  lines(w1_fit*scale(training_set$x1) + w2_fit*scale(training_set$x2) + w3_fit*scale(training_set$x3) + b_fit, col="red")
  
  # Scaled and adjusted learning rate
  w1_fit= 0.000104596;
  w2_fit= 116.443;
  w3_fit= 0.000167483;
  b_fit= 123.5;
  lines(w1_fit*scale(training_set$x1) + w2_fit*scale(training_set$x2) + w3_fit*scale(training_set$x3)  + b_fit, col="green")
}



#######
## Normal equation to get regression coefficients
#######

if (method == "normal") {
  
  # Print to console
  message("Computation with the normal equation hasn't been validated yet! Use at your own peril!")
  message("The normal equation is computationally expensive with a large number of features!") 
  
  # Check args
  if (is.null(dim(x))) { n_examples = length(x); if (!n_examples > 1) { stop("The number of examples must exceed the number of features.") } }
  if (!is.null(dim(x))) { n_examples = dim(x)[1]; n_features = dim(x)[2]; if (!n_examples > n_features) { stop("The number of examples must exceed the number of features.") } }
  
  # Normal equation
  # Note the leading 1s in X, which are for the intercept
  X = as.matrix(cbind(1,features_set)); Y = as.matrix(target_set);
  beta = solve(t(X) %*% X) %*% (t(X) %*% Y)
  return(beta)
  
}


