# ----------------------------------------------------------------------------------------------------------------------------
#                  Wave function library
#
#
#
#
# ----------------------------------------------------------------------------------------------------------------------------
# NOTE: This function are used to calculate drag cape

import numpy as np
from scipy.special import sici


# ----------------------------------------------------------------------------------------------------------------------------

def comp_cdwave(F):
   """
   Compute the cumulative distribution wave function of F, given as
   
   Cd(F) = 4*(F*sin(2/F) - (F**2)*(sin(1/F)**2) - cosint(2/F) + log(2/F) - 1 + gamma_em)
   
   This required in CAPES DRAG function

   Parameters
   ----------
   F : array_like
       The input frequency array
   
   Returns
   -------
   Cd : array_like
       The cumulative distribution wave function of F
   
   Notes
   -----
   The cumulative distribution wave function is derived from the probability
   distribution function of F, given as
   
   p(F) = 4*F*sin(2/F) - (4*F**2)*(sin(1/F)**2)
   
   The cumulative distribution wave function is then computed as
   
   Cd(F) = int[F->inf] p(x) dx
   
   """
   gamma_em = np.euler_gamma
   return 4*(F*np.sin(2/F) - (F**2)*(np.sin(1/F)**2) - cosint(2/F) + np.log(2/F) - 1 + gamma_em)

# ----------------------------------------------------------------------------------------------------------------------------

def cosint(x):
   """
   Return the cosine integral function of x.
   
   Parameters
   ----------
   x : array_like
       The input array

   Returns
   -------
   ci : array_like
       The cosine integral function of x

   Notes
   -----
   The cosine integral function, Ci(x), is defined as [1]_:
   
       Ci(x) = - \int_x^\infty \frac{\cos t}{t} dt
   
   References
   ----------
   .. [1]  https://en.wikipedia.org/wiki/Trigonometric_integral#Cosine_integral
   .. [2]  https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.sici.html
   
   """
   _, ci = sici(x)
   return ci
