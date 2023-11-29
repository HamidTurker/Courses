# Discovering Statistics using R - Andy Field


"################################################################################################### 
 ################################################################################################### 
 #############################                                   ################################### 
 #############################            CHAPTER 1              ################################### 
 #############################                                   ###################################
 ################################################################################################### 
 ###################################################################################################"

# What are the important parts of research? How are theories generated? Why do we need data to test
# those theories? This chapter will attempt to explain what, how, and why.

# There are qualitative methods and quantitative. Our concern is with the latter, which involves
# quantifying - in some fashion - our experimental data and using statistics to make inferences.
# The research process starts with an observation, from which you make a prediction (hypothesis).
# Then, you collect relevant data to test that hypothesis and analyze those data.

# A good theory allows us to make statements about the (state of the) world and, thereby, produce
# hypotheses that are testable (confirmable/falsifiable).

# When we collect data, we need to decide on two things: i. what to measure and ii. how to measure it.
# Then, some of those so-called variables will function as predictors and others as outcomes. Another
# way to refer to them, respectively, is independent versus dependent variables.

# In addition to there being different functions of variables, there are also different levels (of 
# measurement). For instance, categorical variables (i.e. variable 'animal' can have the levels 'dog', 
# 'cat', 'hamster', etc.). If a categorical variable has only two levels, it's also known as a binary
# variable (e.g., gender: 'male', 'female'). If the levels do not have any kind of inherent order,
# rank, or other such relationship, it can also be referred to as a nominal variable (e.g., the jersey
# number of players on a sports team, where player 1 and 32 and 7 are simply labeled by those numbers,
# but player 1 is not better/worse/more/less than player 7) and computing stats like the mean or deviation 
# would be meaningless (computing the 'average jersey number' would not be informative).
#           In addition to categorical variables, there are also continuous variables. Here, the
# measurement can, in principle, take on any (numeric) value and is only limited by the range
# or other limits of the measurement scale we're using. For example, a common continuous variable is the 
# so-called interval variable, where each possible response on our scale is separated by a constant
# interval (e.g., On the following scale, how much do you like X? options: 0 - 0.5 - 1 - 1.5 - 2 - 2.5 - 3).
# So-called ratio variables impose an additional criterion: a 0 on the scale must be 'meaningfully zero'
# and e.g. a doubling of the outcome measure must be 'meaningfully double'. So, in our example, if
# someone picked 0, it must mean their liking of X is 'truly' zero. Not an ounce of liking, nada. Zero.
# And a shift from going from 1 to 2 must mean a doubling of the liking. A commonly encountered ratio
# variable is response/reaction times ( 0 ms is truly 0 and going from 100 ms to 200 ms truly means a
# doubling of the time).




