""" 

https://www.youtube.com/watch?v=HP93wpG2AZI


Model 1:
-----------------------------------------

   Yi = b0 + b1 Xi1 + Xi2 + Ei

   Xi1 = cos (2 pi i/ 365)  (not random)
   Xi2 = sin (2 pi i/ 365)  (not random)

   Ei = N(0,sigma*2)  independent of each other


Model 2: Moving average model
-----------------------------------------

   Yi = b0 + b1 Xi1 + Xi2 + Ei + rho E_{i-1}

   Ei = N(0,sigma*2)  independent of each other


   Y = [Y1, ... Yn ] is a multi-variate normal vector


   E[Yi] = b0 + b1 xi1 + b2 xi2
   V[Yi] = rho sigma**2 + sigma**2 = (rho+1)*sigma**2

   Cov[Yi, Yi-1] = E[ rho Ei-1 + Ei ] [rho Ei-2 + Ei-1] 
                 = rho**2 E[Ei-1 Ei-2] +  rho E[Ei-1 Ei-1] + rho E[ Ei Ei-2] + E[Ei Ei-1] 
                 = 0                   +  rho sigma**2     +  0 + 0
   Cov[Yi, Yj] =  0, |i-j|>1

   Cov = sigma**2 [ rho**2+1  rho      0 ...
                  [ rho      rho**2+1   0
                  [ 0         rho ....


   Can use Residual maximum likelyhood (REML) or maximum likelihood to identify the b's and rho, sigma



Model 3: Autoregressive model
-----------------------------------------

   E[Yi] = b0 + b1 xi1 + b2 xi2 
 
   Yi - E[Yi] =  rho ( Yi-1 - E[Yi-1] ) + Ei
   Yi = E[Yi] +  rho ( Yi-1 - E[Yi-1] ) + Ei

   Ei = N(0,sigma*2)  independent of each other

   We look at deviation between the the response and its expected value.

   Also a multivariate normal vector


   Var[Yi] = E[ ( Yi - E[Yi] ) **2  ]
           = E[   rho**2 ( Yi-1 - E[Yi-1] )**2 + 2 *rho ( Yi-1 - E[Yi-1] )Ei  +  Ei**2  ]
           = rho**2  E[ (Yi-1 - E[Yi-1])**2 ]  +                0           + E[ Ei**2  ]
           = rho**2 Var [Yi-1]                                              + sigma**2
   
    E[ Ei Yi-1 ] = 0 


   Assume stationarity (Var[i-1] = Var[i] 
     => Var[Yi] = sigma**2 / (1-rho**2)


   Cov( Yi, Yi-1)  = E[ (Yi - E[Yi] ) (Yi-1 - E[Yi-1] ) ]
                   = rho Var[Yi] = rho sigma**2 / (1-rho**2)
 
   Cov( Yi, Yk)   = rho**|i-k| sigma**2 / (1-rho**2)
 


Model 4: Higher-order autoregressive model
-----------------------------------------

   Yi - E[Yi] = Sum_j=1..m   rho_j (Yi-j - E[Yi-j] )   + Ei


Model 5: Gaussian Process model
-------------------------------

   Responses : y1, ... yn
   Inputs    : t1, ... tn  (times, statial locations, etc.)
   Covariates: x1, ... xn  (other information observed at those time points)


   In previous models, time was discrete with equispacing.
   In Gaussian Process model, there is no equispacing assumed.
   Tries to account for the fact that correlation between 
   two inputs depends on how close they are to each other.


 Goals: 
    - infer relation between y and x
    - interpolate or extrapolate y at new time tn+1

 
  Model for yi:   
     Yi = b0 + b1 xi  + Zi

     b0, b1 is not random
     xi is observed, part of the dataset

     [Z1, Zi,  Zn]  ~ N(0_, Sigma)

     [Y1, ... Yn] ~ N( Xb_ , Sigma)
 

    How do we get teh covariance matrix Sigma?


    We need a model, a way to get the Cov(Z1, Z2) for any pair of inputs t1 and t2. 
    We need a continuous covariance function "K":

       Cov(Z1, Z2) = K(t1,  t2)

    Once we have K, we can take any pairs of t1,t2 and construct a Cov matrix

    K needs to satisfy certain criteria. In particular, positive definitiveness.


   Some Examples:

     K(t1, t2)  = sigma**2 exp ( - |t2-t1| / r )
     K(t1, t2)  = sigma**2 (1 + |t2 -t1| )  exp ( - |t2-t1| / r )



   We typically use maximum likelihood to estimate parameters in the model (e.g. b0, b1, sigma**2, r)


   We use the multivariate normal distribution to make predictions.












"""

