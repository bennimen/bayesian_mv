# -*- coding: utf-8 -*-
from scipy.stats import beta
import numpy as np
#from calc_prob import calc_prob_between
from math import lgamma

def g0(a, b, c):    
    return np.exp(lgamma(a + b) + lgamma(a + c) - (lgamma(a + b + c) + lgamma(a)))

def h(a, b, c, d):
    num = lgamma(a + c) + lgamma(b + d) + lgamma(a + b) + lgamma(c + d)
    den = lgamma(a) + lgamma(b) + lgamma(c) + lgamma(d) + lgamma(a + b + c + d)
    return np.exp(num - den)


def hiter(a, b, c, d):
    while d > 1:
        d -= 1
        yield h(a, b, c, d) / d

def g(a, b, c, d):
    return g0(a, b, c) + sum(hiter(a, b, c, d))

def calc_prob_between(beta1, beta2):
    return g(beta1.args[0], beta1.args[1], beta2.args[0], beta2.args[1])

#This is the known data: imporessions and conversions for the Control and Test set
imps_ctrl,convs_ctrl=16500, 30 
imps_test, convs_test=17000, 50

imps_ctrl,convs_ctrl=734, 43
imps_test, convs_test=704, 30

imps_ctrl,convs_ctrl=903, 54
imps_test, convs_test=886, 65

imps_ctrl,convs_ctrl=102659,56
imps_test, convs_test=17137, 9

imps_ctrl,convs_ctrl=815,16
imps_test, convs_test=796, 15

imps_ctrl,convs_ctrl=43824,3
imps_test, convs_test=43043,6

imps_ctrl,convs_ctrl=9280,5
imps_test, convs_test=7715,5

imps_ctrl,convs_ctrl=2079,10
imps_test, convs_test=2035,3

imps_ctrl,convs_ctrl=40254,23
imps_test, convs_test=40058,30

imps_ctrl,convs_ctrl=75251,44
imps_test, convs_test=74879,49

imps_ctrl,convs_ctrl=21073,6
imps_test, convs_test=15724,6

imps_ctrl,convs_ctrl=1463,6
imps_test, convs_test=1463,3

imps_ctrl,convs_ctrl=106158+13764,9+3
imps_test, convs_test=6927,1

imps_ctrl,convs_ctrl=26404,14
imps_test, convs_test=26261,22

convs_test_og = convs_test

print('\n')
print('Option A:\nImpressions: {}\nConversions: {}\n'.format(imps_ctrl, convs_ctrl))
print('Option B:\nImpressions: {}\nConversions: {}\n'.format(imps_test, convs_test))

#here we create the Beta functions for the two sets
a_C, b_C = convs_ctrl+1, imps_ctrl-convs_ctrl+1
beta_C = beta(a_C, b_C)
a_T, b_T = convs_test+1, imps_test-convs_test+1
beta_T = beta(a_T, b_T)

#calculating the lift
lift=(beta_T.mean()-beta_C.mean())/beta_C.mean()

#calculating the probability for Test to be better than Control
prob=calc_prob_between(beta_T, beta_C)

print (f"Test option lift Conversion Rates by {lift*100:2.2f}% with {prob*100:2.1f}% probability.")
#output: Test option lift Conversion Rates by 59.68% with 98.2% probability.


while prob*100 < 95:


    #imps_ctrl,convs_ctrl=734, 43
    #imps_test, convs_test=704, 59
    convs_test += 1
    #here we create the Beta functions for the two sets
    a_C, b_C = convs_ctrl+1, imps_ctrl-convs_ctrl+1
    beta_C = beta(a_C, b_C)
    a_T, b_T = convs_test+1, imps_test-convs_test+1
    beta_T = beta(a_T, b_T)

    #calculating the lift
    lift=(beta_T.mean()-beta_C.mean())/beta_C.mean()

    #calculating the probability for Test to be better than Control
    prob=calc_prob_between(beta_T, beta_C)
    convs_test += 1

print('\nConversions in Option B need for at least 95% probability: {}'.format(convs_test))
print('Additional: {}\n'.format(convs_test-convs_test_og))
print (f"SIMULATION: lift Conversion Rates by {lift*100:2.2f}% with {prob*100:2.1f}% probability.")
#output: Test option lift Conversion Rates by 59.68% with 98.2% probability.







# +
##test

# +
import math
import scipy.special
import scipy.special as sc
import numpy as np

def probability_B_beats_A(α_A, β_A, α_B, β_B):
    total = 0.0
    counter = 0-1
    for i in range(0,α_B-1+1):
        counter += 1
        total += math.exp(scipy.special.betaln(α_A+i, β_B+β_A) - math.log(β_B+i) - scipy.special.betaln(1+i, β_B) - scipy.special.betaln(α_A, β_A))
        #print('{}: {}'.format(counter,total))
    return total

def probability_B_beats_A2(alpha_a, beta_a, alpha_b, beta_b):
    ## this function results in higher that double precision numbers
    ## so returns nan without modification
    total = 0.0
    counter = 0-1
    for i in range(0,alpha_b-1+1):
        counter += 1
        num = sc.beta(alpha_a+i, beta_a + beta_b)
        den = (beta_b+i)*sc.beta(1+i,beta_b)*sc.beta(alpha_a,beta_a)
        equals = num/den
        total += equals
        #print('{}: {}\nadded: {}\nnum: {}\nden: {}'.format(counter,total,equals,num,den))
        #print(alpha_a, beta_a, alpha_b, beta_b)
        #print('\n')
    return total

def probability_C_beats_A_and_B(α_A, β_A, α_B, β_B, α_C, β_C):
    total = 0.0
    counter = 0-1
    for i in range(0,α_A-1+1):
        for j in range(0,α_B-1+1):
            counter += 1
            total += math.exp(scipy.special.betaln(α_C+i+j, β_A+β_B+β_C) - math.log(β_A+i) - math.log(β_B+j) - scipy.special.betaln(1+i, β_A) - scipy.special.betaln(1+j, β_B) - scipy.special.betaln(α_C, β_C))
            
            den = -math.log(β_A+i) - math.log(β_B+j) - scipy.special.betaln(1+i, β_A) - scipy.special.betaln(1+j, β_B) - scipy.special.betaln(α_C, β_C)
            print_this = False
            if print_this:
                print('\nalphaC: {}\ni: {}\nj: {}\n'.format(α_C,i,j))
                print('\nbetaA: {}\nbetaB: {}\nbetaC: {}\n'.format(β_A, β_B, β_C))
                print('\nn1: {}\nn2: {}\n'.format(α_C+i+j, β_A+β_B+β_C))
                print('\nnum: {}\nden: {}\n'.format(scipy.special.betaln(α_C+i+j, β_A+β_B+β_C),-den))
                print('{}: {}\n'.format(counter,total))
                print('*'*80)

    return (1 - probability_B_beats_A(α_C, β_C, α_A, β_A) - probability_B_beats_A(α_C, β_C, α_B, β_B) + total)




# +
#sc.beta(1.7, 2.4)

# +
#math.log(sc.beta(1.7, 2.4))
# -

imps_ctrl,convs_ctrl=75251,44
imps_test, convs_test=74879,49

imps_ctrl,convs_ctrl=1000,1
imps_test, convs_test=1000,3

imps_ctrl,convs_ctrl=1000,8
imps_test, convs_test=1000,3

alpha_a = 1 + convs_ctrl
#beta_a = imps_ctrl - convs_ctrl + 1
beta_a = imps_ctrl - alpha_a
alpha_b = 1 + convs_test
#beta_b = imps_test - convs_test + 1
beta_b = imps_test - alpha_b


results = probability_B_beats_A(alpha_a,beta_a,alpha_b,beta_b)
results

results2 = probability_B_beats_A2(alpha_b,beta_b,alpha_a,beta_a)
results2

# +
#math.log(2)

# +
imps_t3, convs_t3 = 5149,2

alpha_c = 1 + convs_t3
beta_c = imps_t3 - alpha_c

# +
imps_t3, convs_t3 = 1000,8
imps_t3, convs_t3 = 1000,1

alpha_c = 1 + convs_t3
beta_c = imps_t3 - alpha_c
# -

resultsC = probability_C_beats_A_and_B(alpha_a,beta_a,alpha_b,beta_b,alpha_c,beta_c)
resultsC

resultsB = probability_C_beats_A_and_B(alpha_a,beta_a,alpha_c,beta_c,alpha_b,beta_b)
resultsB

resultsA = probability_C_beats_A_and_B(alpha_b,beta_b,alpha_c,beta_c,alpha_a,beta_a)
resultsA

resultsA + resultsB + resultsC





# +
imps_t4, convs_t4 = 3000,33

alpha_d = 1 + convs_t4
beta_d = imps_t4 - convs_t4
# -



# +
def bln(one,two):
    return scipy.special.betaln(one, two)

def summation(alpha,beta):
    ## returns probabilty that first value (?) [0]
    ## is likely to be better than the rest
    sizer = len(alpha)
    summer = np.zeros((len(alpha),1))
    beta0 = bln(alpha[0],beta[0])
    
    done = False
    l2 = False

    total = 0

    counter = 0-1
    while not done:
        counter += 1
        #print('Iteration: {}\n'.format(counter))

        num = 0
        den = 0

        n1 = sum(summer[1:len(summer)])+alpha[0]
        #print('\nalpha_c: {}\ni: {}\nj: {}\nn1: {}\n'.format(alpha[0],summer[1],summer[2],n1))
        #n2 = beta0
        n2 = sum(beta)
        #print('\nbeta_a: {}\nbeta_b: {}\nbeta_c: {}\nsum: {}\n'.format(beta[0],beta[1],beta[2],sum(beta)))
        num = bln(n1,n2)
        #print('num: {}'.format(num))

        ## den; part 01/02
        c2 = 0-1
        for de in range(0,sizer):
            c2 += 1
            den += math.log(beta[de] + summer[de])

        ## den; part 02/02
        c3 = 0-1
        for de in range(0,sizer):
            c3 += 1
            den += bln(1+summer[de],beta[de])

        den += beta0
        #print('den: {}'.format(den))
        #print('\n')

        ## total; part 01/01
        total += math.exp(num-den)
        #print('thisTotal: {}\ntotal: {}'.format(math.exp(num-den),total))

        printer = True
        printer = False
        if printer:
            print('Iteration: {}\n'.format(counter))
            print('\nalpha_c: {}\ni: {}\nj: {}\nn1: {}\n'.format(alpha[0],summer[1],summer[2],n1))
            print('\nbeta_a: {}\nbeta_b: {}\nbeta_c: {}\nsum: {}\n'.format(beta[0],beta[1],beta[2],sum(beta)))
            print('num: {}'.format(num))
            print('den: {}'.format(den))
            print('\n')
            print('thisTotal: {}\ntotal: {}'.format(math.exp(num-den),total))
            print('\n')
            print('*'*80,'\n')

        
        
        #print('\n')
        ## summer manipulation
        for cd in range(sizer-1,0-1,-1):
            #print(cd)
            if not summer[cd]+1 >= alpha[cd]:
                #print('if1')
                #print('summer: {}\nalpha: {}'.format(summer[cd],alpha[cd]))
                summer[cd] += 1
                break
            else:
                #print('else')
                l2count = 0
                for cd2 in range(sizer-1,0-1,-1):
                    if summer[cd2]+1 == alpha[cd2]: 
                        l2count += 1
                if l2count+1 == sizer:
                    l2=True

                summer[cd] = 0
        #l2 = True
        if l2:
            done = True
        #print('\n')
        #print('summerArray: {}'.format(summer))




        ## temporary
        #if counter >=10: done = True
        #print('\n')

        #print('*'*80,'\n')

    total
    return total

def alpha_beta(my_array):
    alpha = []
    beta = []
    for i in my_array:
        alpha.append(i[1]+1)
        beta.append(i[0]-(i[1]+1))
    #print(alpha)
    #print(beta)
    return alpha,beta

#alpha_beta(my_array)
# -





# +
## test data
my_array = []
my_array.append([1000,8])
my_array.append([1000,3])
my_array.append([1000,1])
#my_array.append([3000,33])
#my_array.append([1000,11])

print(my_array)
print('\n')

# alpha = []
# beta = []
# for i in my_array:
#     alpha.append(i[1]+1)
#     beta.append(i[0]-(i[1]+1))

alpha, beta = alpha_beta(my_array)

print(alpha)
print(beta)


# +
p21 = probability_B_beats_A2(alpha[0],beta[0],alpha[1],beta[1])
p31 = probability_B_beats_A2(alpha[0],beta[0],alpha[2],beta[2])
#p41 = probability_B_beats_A2(alpha[0],beta[0],alpha[3],beta[3])

# p213 = probability_C_beats_A_and_B(alpha[0],beta[0],alpha[2],beta[2],alpha[1],beta[1])
# p214 = probability_C_beats_A_and_B(alpha[0],beta[0],alpha[3],beta[3],alpha[1],beta[1])

# p312 = probability_C_beats_A_and_B(alpha[0],beta[0],alpha[1],beta[1],alpha[2],beta[2])
# p314 = probability_C_beats_A_and_B(alpha[0],beta[0],alpha[3],beta[3],alpha[2],beta[2])

# p412 = probability_C_beats_A_and_B(alpha[0],beta[0],alpha[1],beta[1],alpha[3],beta[3])
# p413 = probability_C_beats_A_and_B(alpha[0],beta[0],alpha[2],beta[2],alpha[3],beta[3])


# +
print('p21: {}'.format(p21))
print('p31: {}'.format(p31))
# print('p41: {}'.format(p41))
# print('\n')

# print('p213: {}'.format(p213))
# print('p214: {}'.format(p214))
# print('\n')

# print('p312: {}'.format(p312))
# print('p314: {}'.format(p314))
# print('\n')

# print('p412: {}'.format(p412))
# print('p413: {}'.format(p413))
# -


total = summation(alpha,beta)
total


# +
# p1234 = 1 - p21 - p31 - p41 - p213 - p214 - p312 - p314 - p412 - p413 + total
# p1234

p123 = 1 - p21 - p31 + total
p123
# -








# +
## test data
my_array = []
my_array.append([1000,8])
my_array.append([1000,3])
my_array.append([1000,1])
my_array.append([3000,33])
#my_array.append([1000,11])

print(my_array)
print('\n')

# alpha = []
# beta = []
# for i in my_array:
#     alpha.append(i[1]+1)
#     beta.append(i[0]-(i[1]+1))

alpha, beta = alpha_beta(my_array)

print(alpha)
print(beta)
# -


## retuns total for first term of {alpha,beta} being focus
## alpha[0], beta[0]
total = summation(alpha,beta)
total


# +
p21 = probability_B_beats_A2(alpha[0],beta[0],alpha[1],beta[1])
p31 = probability_B_beats_A2(alpha[0],beta[0],alpha[2],beta[2])
p41 = probability_B_beats_A2(alpha[0],beta[0],alpha[3],beta[3])

p213 = probability_C_beats_A_and_B(alpha[0],beta[0],alpha[2],beta[2],alpha[1],beta[1])
p214 = probability_C_beats_A_and_B(alpha[0],beta[0],alpha[3],beta[3],alpha[1],beta[1])

p312 = probability_C_beats_A_and_B(alpha[0],beta[0],alpha[1],beta[1],alpha[2],beta[2])
p314 = probability_C_beats_A_and_B(alpha[0],beta[0],alpha[3],beta[3],alpha[2],beta[2])

p412 = probability_C_beats_A_and_B(alpha[0],beta[0],alpha[1],beta[1],alpha[3],beta[3])
p413 = probability_C_beats_A_and_B(alpha[0],beta[0],alpha[2],beta[2],alpha[3],beta[3])


# +
print('p21: {}'.format(p21))
print('p31: {}'.format(p31))
print('p41: {}'.format(p41))
print('\n')

print('p213: {}'.format(p213))
print('p214: {}'.format(p214))
print('\n')

print('p312: {}'.format(p312))
print('p314: {}'.format(p314))
print('\n')

print('p412: {}'.format(p412))
print('p413: {}'.format(p413))
print('\n')

print('total: {}'.format(total))

# -


## value for these numbers should be ~0.48
p1234 = 1 - p21 - p31 - p41 - p213 - p214 - p312 - p314 - p412 - p413 + total
p1234



# +
p12 = 1 - p21 
p13 = 1 - p31
p14 = 1 - p41

p123 = 1 - p312 - p213
p134 = 1 - p413 - p314
p124 = 1 - p214 - p412




# +
print('p12: {}'.format(p12))
print('p13: {}'.format(p13))
print('p14: {}'.format(p14))
print('\n')

print('p123: {}'.format(p123))
print('p134: {}'.format(p134))
print('p124: {}'.format(p124))
print('\n')
# -


p1234real = 1 - p12 - p13 - p14 + p123 + p124 + p134 - total
p1234real






def probability_D_beats_ABC(α_A, β_A, α_B, β_B, α_C, β_C, α_D, β_D):
    num = 0
    den = 0
    total = 0.0
    counter = 0-1
    for i in range(0,α_A-1+1):
        for j in range(0,α_B-1+1):
            for k in range(0,α_C-1+1):
                counter += 1
                
                num = scipy.special.betaln(α_D+i+j+k, β_A+β_B+β_C+β_D)
                den = -math.log(β_A+i) - math.log(β_B+j) - math.log(β_C+k) - scipy.special.betaln(1+i, β_A) - scipy.special.betaln(1+j, β_B) - scipy.special.betaln(1+k, β_C) - scipy.special.betaln(α_D, β_D)
                result = math.exp(num + den)
                total += result
                
                #print('{}: num: {} | den: {}\nresult: {} | total: {}\n'.format(counter,num,den,result,total))
                #print(result)
                #print(total)
    #print(total)
    
    print_this = False
    if print_this:
        print('d > a: {}'.format(probability_B_beats_A(α_A, β_A, α_D, β_D)))
        print('d > b: {}'.format(probability_B_beats_A(α_B, β_B, α_D, β_D)))
        print('d > c: {}'.format(probability_B_beats_A(α_C, β_C, α_D, β_D)))

        print('d > [a,b]: {}'.format(probability_C_beats_A_and_B(α_A, β_A, α_B, β_B, α_D, β_D)))
        print('d > [a,c]: {}'.format(probability_C_beats_A_and_B(α_A, β_A, α_C, β_C, α_D, β_D)))
        print('d > [b,c]: {}'.format(probability_C_beats_A_and_B(α_B, β_B, α_C, β_C, α_D, β_D)))
        
        print('total: {}\n'.format(total))
              
    return (1 - probability_B_beats_A(α_A, β_A, α_D, β_D)
              - probability_B_beats_A(α_B, β_B, α_D, β_D)
              - probability_B_beats_A(α_C, β_C, α_D, β_D)
              + probability_C_beats_A_and_B(α_A, β_A, α_B, β_B, α_D, β_D)
              + probability_C_beats_A_and_B(α_A, β_A, α_C, β_C, α_D, β_D)
              + probability_C_beats_A_and_B(α_B, β_B, α_C, β_C, α_D, β_D)
              - total)



resultsD4=probability_D_beats_ABC(alpha_a,beta_a,alpha_b,beta_b,alpha_c,beta_c,alpha_d,beta_d)
resultsD4

resultsC4=probability_D_beats_ABC(alpha_d,beta_d,alpha_a,beta_a,alpha_b,beta_b,alpha_c,beta_c)
resultsC4

resultsB4=probability_D_beats_ABC(alpha_c,beta_c,alpha_d,beta_d,alpha_a,beta_a,alpha_b,beta_b)
resultsB4

resultsA4=probability_D_beats_ABC(alpha_b,beta_b,alpha_c,beta_c,alpha_d,beta_d,alpha_a,beta_a)
resultsA4

print('a: {}'.format(resultsA4))
print('b: {}'.format(resultsB4))
print('c: {}'.format(resultsC4))
print('d: {}'.format(resultsD4))

resultsA4 + resultsB4 + resultsC4 + resultsD4





import fiverr

fiverr.probability_B_beats_A(alpha_a, beta_a, alpha_b, beta_b)

fiverr.probability_C_beats_A_and_B(alpha_a, beta_a, alpha_b, beta_b, alpha_c, beta_c)

dd = fiverr.probability_D_beats_ABC(alpha_a, beta_a, alpha_b, beta_b, alpha_c, beta_c, alpha_d, beta_d)

cc = fiverr.probability_D_beats_ABC(alpha_d, beta_d, alpha_a, beta_a, alpha_b, beta_b, alpha_c, beta_c)

bb = fiverr.probability_D_beats_ABC(alpha_c, beta_c, alpha_a, beta_a, alpha_d, beta_d, alpha_b, beta_b)

aa = fiverr.probability_D_beats_ABC(alpha_d, beta_d, alpha_b, beta_b, alpha_c, beta_c, alpha_a, beta_a)

print('a: {}'.format(aa))
print('b: {}'.format(bb))
print('c: {}'.format(cc))
print('d: {}'.format(dd))

check = aa + bb + cc + dd
check

# +
imps_t5, convs_t5 = 1000,11

alpha_e = 1 + convs_t5
beta_e = imps_t5 - convs_t5
# -

fiverr.probability_E_beats_ABC(alpha_a, beta_a, alpha_b, beta_b, alpha_c, beta_c, alpha_d, beta_d, alpha_e, beta_e)

print(alpha_a, beta_a, alpha_b, beta_b, alpha_c, beta_c, alpha_d, beta_d, alpha_e, beta_e)

fiverr.probability_E_beats_ABC2(alpha_a, beta_a, alpha_b, beta_b, alpha_c, beta_c, alpha_d, beta_d, alpha_e, beta_e)



# +
imps_t6, convs_t6 = 1000,12

alpha_f = 1 + convs_t6
beta_f = imps_t6 - convs_t6

# +
#fiverr.probability_F_beats_ABC(alpha_a, beta_a, alpha_b, beta_b, alpha_c, beta_c, alpha_d, beta_d, alpha_e, beta_e,alpha_f,beta_f)
# -



# +
print('impressions_a: {}\nleads_a: {}\n'.format(imps_ctrl,convs_ctrl))
print('impressions_b: {}\nleads_b: {}\n'.format(imps_test,convs_test))
print('impressions_c: {}\nleads_c: {}\n'.format(imps_t3,convs_t3))
print('impressions_d: {}\nleads_d: {}\n'.format(imps_t4,convs_t4))
print('impressions_e: {}\nleads_e: {}\n'.format(imps_t5,convs_t5))
print('impressions_f: {}\nleads_f: {}\n'.format(imps_t6,convs_t6))

print('*'*80)

print('alpha_a: {}\nbeta_a: {}\n'.format(alpha_a,beta_a))
print('alpha_b: {}\nbeta_b: {}\n'.format(alpha_b,beta_b))
print('alpha_c: {}\nbeta_c: {}\n'.format(alpha_c,beta_c))
print('alpha_d: {}\nbeta_d: {}\n'.format(alpha_d,beta_d))
print('alpha_e: {}\nbeta_e: {}\n'.format(alpha_e,beta_e))
print('alpha_f: {}\nbeta_f: {}\n'.format(alpha_f,beta_f))
# -



# +
import time
from time import gmtime, strftime, localtime

import logging


logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO) ## INFO, DEBUG, WARNING

print('Starting now...\n')

time_started = strftime("%d%m%Y-%H%M%S", gmtime())
# time_now(now)
now = gmtime()
now = localtime()
time_now = lambda x: strftime("%d%m%Y-%H%M%S", x)
ntime_now = lambda x: strftime("%d%m%Y", x)
today = ntime_now(now)


now = localtime()
today = ntime_now(now)


k = 0    
start = time.time()

# +
run_iter = False
if run_iter:
    my_arrays = []
    my_arrays.append([10000,4])
    my_arrays.append([10000,5])
    my_arrays.append([10000,6])
    my_arrays.append([10000,7])
    my_arrays.append([10000,8])
    my_arrays.append([10000,9])
    my_arrays.append([10000,10])
    my_arrays.append([10000,11])
    # my_arrays.append([10000,12])
    # my_arrays.append([10000,13])

    alphas, betas = alpha_beta(my_arrays)

    print(my_arrays)
    print(len(my_arrays))
    print('\n')
    print(alphas)
    print(betas)

    totals = []
    print('Starting iterations...\n')
    for i in range(0,len(my_arrays)):

        start2 = time.time()

        ttt = summation(alphas,betas)
        print(ttt)
        totals.append(ttt)

        toend = my_arrays.pop(0)
        my_arrays.append(toend)

        alphas, betas = alpha_beta(my_arrays)

        ####################################################################################################################
        ## finally

        end = time.time()
        duration = end - start2
        print('Code snippet took %4.3f seconds' % (duration))
        print('Code snippet took %4.3f minutes' % (duration/60))

        print('*'*80)

    print('#'*80)



# -

if run_iter:
    ####################################################################################################################
    ## finally

    end = time.time()
    duration = end - start
    print('Code snippet took %4.3f seconds' % (duration))
    print('Code snippet took %4.3f minutes' % (duration/60))






