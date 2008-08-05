'''
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright
    notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
    notice, this list of conditions and the following disclaimer in the
    documentation and/or other materials provided with the distribution.
    * Neither the name of the University of Auckland, New Zealand nor
    the names of its contributors may be used to endorse or promote
    products derived from this software without specific prior written
    permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL,EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.'''

import scipy.stats, numpy as n
from copy import copy
import unittest
from numpy.testing import NumpyTestCase

# data - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

winglengths = [36] + \
              [37] + \
              [38] * 2 + \
              [39] * 2 + \
              [40] * 4 + \
              [41] * 6 + \
              [42] * 7 + \
              [43] * 8 + \
              [44] * 9 + \
              [45] * 10 + \
              [46] * 10 + \
              [47] * 9 + \
              [48] * 8 + \
              [49] * 7 + \
              [50] * 6 + \
              [51] * 4 + \
              [52] * 2 + \
              [53] * 2 + \
              [54] + \
              [55]

#check winglengths mean = 45.5, var = 15.21,
#std. (pop) = 3.9, std. (sample) = 3.91

def box4_1():
    '''Birth weights of male Chinese in ounces. Returns Y, f.'''
    Y = n.arange(15) * 8 + 59.5
    f = [2, 6, 39, 385, 888, 1729, 2240, 2007, 1233, 641, 201, 74, 14, 5, 1]
    return Y, f

def box8_1():
    n_females = 10
    n_males = 10
    mean_females = 8.5 # days
    mean_males = 4.8 # days

def box9_4():
    a = [75, 67, 70, 75, 65, 71, 67, 67, 76, 68]
    b = [57, 58, 60, 59, 62, 60, 60, 57, 59, 61]
    c = [58, 61, 56, 58, 57, 56, 61, 60, 57, 58]
    d = [58, 59, 58, 61, 57, 56, 58, 57, 57, 59]
    e = [62, 66, 65, 63, 64, 62, 65, 65, 62, 67]
    return (a, b, c, d, e)

def box13_7():
    a = [104, 109, 112, 114, 116, 118, 118, 119, 121, 123, \
         125, 126, 126, 128, 128, 128]
    b = [100, 105, 107, 107, 108, 111, 116, 120, 121, 123]
    return a, b

def box15_2():
    '''body weight of 12 crabs in grams'''
    return [1.41, 2.5, 4.19, 9.52, 11.30, 14.40, \
            14.90, 15.20, 15.39, 15.81, 17.25, 22.70]
    


# functions - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    
def significance_of_difference_between_two_variances( \
    ar1=None, ar2=None, alpha=0.05, \
    v1=None, n1=None, v2=None, n2=None):
    '''from box 8.1, Sokal and Rohlf'''

    if ar1 is None:
        if v1 is None or n1 is None:
            raise ValueError("Either ar1 or s1 & n1 must be specified.")
    else:
        v1, n1 = scipy.stats.var(ar1), ar1.size
    if ar2 is None:
        if v2 is None or n2 is None:
            raise ValueError("Either ar2 or s2 & n2 must be specified.")
    else:
        v2, n2 = scipy.stats.var(ar2), ar2.size

    print "var1 = %0.2f, n1 = %d" % (v1, n1)
    print "var2 = %0.2f, n2 = %d" % (v2, n2)
    print ""
    
    if v1 > v2:
        numerator, denominator = v1, v2
        df1, df2 = n1 - 1, n2 - 1
    else:
        numerator, denominator = v2, v1
        df1, df2 = n2 - 1, n1 - 1
    Fs = numerator / denominator
    print "Fs = %0.2f / %0.2f = %0.2f" % (numerator, denominator, Fs)
    print "df1 = %d, df2 = %d" % (df1, df2)
    print ""

    F_upper = scipy.stats.f.ppf(1-alpha/2., df1, df2)
    F_lower = scipy.stats.f.ppf(alpha/2., df1, df2)
    print "F_lower = %0.3f" % (F_lower)
    print "F_upper = %0.3f" % (F_upper)

    if Fs < F_lower or Fs > F_upper:
        return True
    else:
        return False

def welchs_approximate_ttest(n1, mean1, sem1, \
                            n2, mean2, sem2, alpha):
    '''Welch''s approximate t-test for the difference of two means of
heteroscedasctic populations.

Implemented from Biometry, Sokal and Rohlf, 3rd ed., 1995, Box 13.4

:Parameters:
    n1 : int
        number of variates in sample 1
    n2 : int
        number of variates in sample 2
    mean1 : float
        mean of sample 1
    mean2 : float
        mean of sample 2
    sem1 : float
        standard error of mean1
    sem2 : float
        standard error of mean2
    alpha : float
        desired level of significance of test

:Returns:
    significant : bool
        True if means are significantly different, else False
    t_s_prime : float
        t_prime value for difference of means
    t_alpha_prime : float
        critical value of t_prime at given level of significance
    
Copyright (c) 2007, Angus McMorland, released under BSD licence.
'''
    svm1 = sem1**2 * n1
    svm2 = sem2**2 * n2
    t_s_prime = (mean1 - mean2)/n.sqrt(svm1/n1+svm2/n2)
    t_alpha_df1 = scipy.stats.t.ppf(1-alpha/2, n1 - 1)
    t_alpha_df2 = scipy.stats.t.ppf(1-alpha/2, n2 - 1)
    t_alpha_prime = (t_alpha_df1 * sem1**2 + t_alpha_df2 * sem2**2) / \
                    (sem1**2 + sem2**2)
    return abs(t_s_prime) > t_alpha_prime, t_s_prime, t_alpha_prime

def F_delta(delta, n_, F):
    '''Returns the F_delta value for delta-correction of the
    Kolmogorov-Smirnov goodness of fit test, as per Sokal and Rohlf, p711

    INPUTS:
      i = array of data indices, 1-based to n'''
    return (F - delta)/(n_ - 2 * delta + 1.)

def std_with_f(Y, f, mean=None):
    Y = n.asarray(Y)
    f = n.asarray(f)
    
    if mean==None:
        mean = (Y * f).sum() / f.sum()

    y = Y - mean
    var = (f * y**2).sum() / (f.sum() - 1)
    return n.sqrt(var)
        
def ks_samples(data, F_hat=None):
    '''Test for goodness of fit, using delta adjustment, for individual
    data points Taken from Sokal and Rohlf, Table 17.4. Requires
    comparing the returned g values with table X (extrinsic hypothesis)
    or Y (intrinsic hypothesis) from the Biometry Statistical Tables book.'''
    data = copy(n.asarray(data))
    data.sort()
    n_ = data.size
    print "n = %d" % n_
    i = n.arange(n_) + 1 # 1, 2, 3,...
    if F_hat == None:
        mean = data.mean()
        s = scipy.stats.std(data)
        stded_devs = (data - mean) / s
        F_hat = scipy.stats.norm.cdf(stded_devs)
    F_0 = F_delta(  0. , n_, i)
    F_0p5 = F_delta(0.5, n_, i)
    F_1 = F_delta(  1. , n_, i)
    g_0   = n.abs(F_hat - F_0)
    g_0p5 = n.abs(F_hat - F_0p5)
    g_1   = n.abs(F_hat - F_1)
    d_max = g_0p5.max() + 1/(2. * n_)
    print "g0 max = %f; g0.5 max = %f, g1 max = %f" % \
          (g_0.max(), g_0p5.max(), g_1.max())
    if n_ > 100:
        g0_5c = 0.89196 / n.sqrt(n_) - (1/(2. * n_))
        g0_5d = 1.0427 / n.sqrt(n_) - (1/(2. * n_))
        print "g_0.5 critical (a=0.05) = %.3f" % g0_5c
        print "g_0.5 critical (a=0.01) = %.3f" % g0_5d
    return None


def ties(data):
    data = n.asarray(data)
    return data.size <> n.unique(data).size

def sumT_j(data):
    '''copied from scipy.stats.tiecorrect
    sorted,posn = fastsort(asarray(rankvals))
    
    returns sum^m(T_j) as per box 3.6
    '''
    sorted,posn = scipy.stats.fastsort(n.asarray(data))
    num = len(sorted)
    T = 0.0
    i = 0
    while (i<num-1):
        if sorted[i] == sorted[i+1]:
            nties = 1
            while (i<num-1) and (sorted[i] == sorted[i+1]):
                nties = nties +1
                i = i +1
            T = T + nties**3 - nties
        i = i+1
    return T

def wilcoxon(da, db):
    '''Wilcoxon non-parametric test for a significant difference between the
    location (mean) of two samples - unlike scipy.stats this variant allows
    for uneven n numbers in the two groups. Requires n_max > 20 to allow use of
    normal approximation and t-estimation of probabilities.

    Returns t value for comparison and p value. You want p > level of
    significance for rejecting null hypothesis.'''
    da = n.asarray(da)
    na = da.size
    db = n.asarray(db)
    nb = db.size

    if da.size > db.size:
        d1, d2 = da, db
    else:
        d1, d2 = db, da
    n1 = d1.size
    n2 = d2.size

    alld = n.hstack((d1,d2))
    allr = scipy.stats.rankdata(alld)
    d2r = allr[n1:]
    C = n1 * n2 + n2*(n2+1)/2 - d2r.sum()

    U = max([C, n1 * n2 - C])

    if not ties(alld):
        # simple formula
        t = (U - n1 * n2 / 2.) / n.sqrt( n1 * n2 * (n1 + n2 + 1) / 12. )
    else:
        pass
    t_denom1 = n1 * n2 / ((n1 + n2) * (n1 + n2 - 1.))
    t_denom2 = ((n1 + n2)**3 - (n1 + n2) - sumT_j(alld)) / 12.
    t = (U - n1 * n2 / 2) / n.sqrt(t_denom1 * t_denom2)

    p = scipy.stats.t.sf(t, 1e20) * 2. # 1e20 simulates infinite d.f.

    return t, p

# def two_way_anova_with_replication(data):
#     '''Implements box 11.2, two-way ANOVA with replication,
#     also known as randomised blocks.
    
#     Data is presented as an array (cols x rows x replicates)

#     not finished
#     '''
#     # grand mean
#     grandmean = data.mean()
#     print "Grand mean = %.4f" % grandmean

#     # SS_A (SS of columns)
#     # mean of columns
#     colmeans = data.mean(axis=2).mean(axis=1)
#     print "column means = ", colmeans
#     SS_A = ((colmeans - grandmean)**2).sum() 
#     print "SS_A = %.4f" % SS_A

# tests - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    
def test_significance_of_difference_between_two_variances():
    '''from the values taken in Box 8.1'''
    v1, n1 = 3.6, 10
    v2, n2 = 0.9, 10
    result = significance_of_difference_between_two_variances(\
        v1=v1, n1=n1, v2=v2, n2=n2)
    return result == False

def test_kolmogorov_smirnov_test():
    data = n.array([ 1.41, 2.50, 4.19, 9.52, 11.30, 14.40, \
                         14.90, 15.20, 15.39, 15.81, 17.25, 22.70 ])
    kolmogorov_smirnov_one(data)

def test_wilcoxon():
    d1 = [104, 109, 112, 114, 116, 118, 118, 119,
          121, 123, 125, 126, 126, 128, 128, 128]
    d2 = [100, 105, 107, 107, 108, 111, 116, 120, 121, 123]
    n1, n2, U = wilcoxon(d1, d2)
    print "n1 = %d, n2 = %d, U_s = %f" % (n1, n2, U)

# def test_two_way_anova_without_replication():
#     jul29 = n.array((23.8, 22.6, 22.2, 21.2, 18.4, \
#                      13.5, 9.8, 6.0, 5.8, 5.6))
#     jul30 = n.array((24.0, 22.4, 22.1, 21.8, 19.3, \
#                      14.4, 9.9, 6.0, 5.9, 5.6))
#     jul31 = n.array((24.6, 22.9, 22.1, 21.0, 19.0, \
#                      14.2, 10.4, 6.3, 6.0, 5.5))
#     aug01 = n.array((24.8, 23.2, 22.2, 21.2, 18.8, \
#                      13.8, 9.6, 6.3, 5.8, 5.6))
#     result = two_way_anova_without_replication([jul29, jul30, jul31, aug01])

def test_two_way_anova_with_replication():
    a_scabra100 = n.array((7.16, 6.78, 13.60, 8.93, 8.26, 14.00, 16.10, 9.66))
    a_scabra75 = n.array((5.2, 13.2, 5.2, 8.39, 7.18, 10.4, 6.37, 7.18))
    a_scabra50 = n.array((11.11, 10.5, 9.74, 14.6, \
                          18.80, 11.1, 9.74, 11.8))
    a_digit100 = n.array((6.14, 6.14, 3.86, 10., 10.4, 11.6, 5.49, 5.8))
    a_digit75 = n.array((4.47, 4.95, 9.9, 6.49, 5.75, 5.44, 11.8, 9.9))
    a_digit50 = n.array((9.63, 14.5, 6.38, 10.2, 13.4, 17.7, 14.5, 12.3))
    a_scabra = n.array((a_scabra100, a_scabra75, a_scabra50))
    a_digit = n.array((a_digit100, a_digit75, a_digit50))
    data = n.array((a_scabra, a_digit))
    results = two_way_anova_with_replication(data)

class TestBiometry(NumpyTestCase):
    def test_welchs_approximate_ttest(self):
        chimpanzees = (37, 0.115, 0.017) # n, mean, sem
        gorillas = (6, 0.511, 0.144)
        case1 = welchs_approximate_ttest(chimpanzees[0], \
                                    chimpanzees[1], \
                                    chimpanzees[2], \
                                    gorillas[0], \
                                    gorillas[1], \
                                    gorillas[2], \
                                    0.05)
        self.assertTrue( case1[0] )
        self.assertAlmostEqual( case1[1], -2.73, 2 )
        self.assertAlmostEqual( case1[2], 2.564, 2 )
        
        female = (10, 8.5, n.sqrt(3.6)/n.sqrt(10))
        male = (10, 4.8, n.sqrt(0.9)/n.sqrt(10))
        case2 = welchs_approximate_ttest(female[0], \
                                female[1], \
                                female[2], \
                                male[0], \
                                male[1], \
                                male[2], 0.001)
        self.assertTrue( case2[0] )
        self.assertAlmostEqual( case2[1], 5.52, 2 )
        self.assertAlmostEqual( case2[2], 4.781, 2 )
        

def test():
    suite = unittest.TestLoader().loadTestsFromTestCase( \
        testBiometry)
    unittest.TextTestRunner(verbosity=2).run(suite)
