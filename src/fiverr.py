# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 07:43:39 2020

@author: MR UMAIR
"""

import scipy.special
import math

def bln(one, two):
    return scipy.special.betaln(one, two)


def probability_B_beats_A(alpha_A, beta_A, alpha_B, beta_B):
    total = 0.0
    for i in range(0, (alpha_B)):
        total += math.exp(bln(alpha_A + i, beta_B + beta_A)
                     - math.log(beta_B + i) - bln(1 + i, beta_B) - bln(alpha_A, beta_A))

    return total

def probability_C_beats_A_and_B(alpha_A, beta_A, alpha_B, beta_B, alpha_C, beta_C):
    total = 0.0
    for i in range(0, (alpha_A)):
        for j in range(0, (alpha_B)):
            total += math.exp(bln(alpha_C + i + j, beta_A + beta_B + beta_C) - math.log(beta_A + i) - math.log(beta_B + j)
                         - bln(1 + i, beta_A) - bln(1 + j, beta_B) - bln(alpha_C, beta_C))

    return (1 - probability_B_beats_A(alpha_C, beta_C, alpha_A, beta_A)
            - probability_B_beats_A(alpha_C, beta_C, alpha_B, beta_B) + total)

def probability_D_beats_ABC(alpha_A, beta_A, alpha_B, beta_B, alpha_C, beta_C, alpha_D, beta_D):
    total = 0.0
    for i in range(0, (alpha_A)):
        for j in range(0, (alpha_B)):
            for k in range(0, (alpha_C)):



                total += math.exp(bln(alpha_D + i + j + k, beta_A + beta_B + beta_C + beta_D)
                             - math.log(beta_A + i) - math.log(beta_B + j) - math.log(beta_C + k)
                             - bln(1 + i, beta_A) - bln(1 + j, beta_B) - bln(1 + k, beta_C) 
                             - bln(alpha_D, beta_D))

    return (1 - probability_B_beats_A(alpha_A, beta_A, alpha_D, beta_D)
            - probability_B_beats_A(alpha_B, beta_B, alpha_D, beta_D)
            - probability_B_beats_A(alpha_C, beta_C, alpha_D, beta_D)
            + probability_C_beats_A_and_B(alpha_A, beta_A, alpha_B, beta_B, alpha_D, beta_D)
            + probability_C_beats_A_and_B(alpha_A, beta_A, alpha_C, beta_C, alpha_D, beta_D)
            + probability_C_beats_A_and_B(alpha_B, beta_B, alpha_C, beta_C, alpha_D, beta_D)
            - total)

def probability_E_beats_ABC(alpha_A, beta_A, alpha_B, beta_B, alpha_C, beta_C, alpha_D, beta_D, alpha_E, beta_E):
    total = 0.0
    for i in range(0, (alpha_A)):
        for j in range(0, (alpha_B)):
            for k in range(0, (alpha_C)):
                for l in range(0, (alpha_D)):

                    total += math.exp(bln(alpha_E + i + j + k + l, beta_A + beta_B + beta_C + beta_D + beta_E)
                                 - math.log(beta_A + i) - math.log(beta_B + j) - math.log(beta_C + k) - math.log(beta_D + l)
                                 - bln(1 + i, beta_A) - bln(1 + j, beta_B) - bln(1 + k, beta_C) - bln(1 + l, beta_D)
                                 - bln(alpha_E, beta_E))



    return (1 - probability_B_beats_A(alpha_A, beta_A, alpha_E, beta_E)
            - probability_B_beats_A(alpha_B, beta_B, alpha_E, beta_E)
            - probability_B_beats_A(alpha_C, beta_C, alpha_E, beta_E)
            + probability_C_beats_A_and_B(alpha_A, beta_A, alpha_B, beta_B, alpha_E, beta_E)
            + probability_C_beats_A_and_B(alpha_A, beta_A, alpha_C, beta_C, alpha_E, beta_E)
            + probability_C_beats_A_and_B(alpha_B, beta_B, alpha_C, beta_C, alpha_E, beta_E)
            + probability_D_beats_ABC(alpha_A, beta_A, alpha_B, beta_B, alpha_C, beta_C, alpha_E, beta_E)
            + probability_D_beats_ABC(alpha_A, beta_A, alpha_B, beta_B, alpha_D, beta_D, alpha_E, beta_E)
            + probability_D_beats_ABC(alpha_A, beta_A, alpha_D, beta_D, alpha_C, beta_C, alpha_E, beta_E)
            + probability_D_beats_ABC(alpha_B, beta_B, alpha_C, beta_C, alpha_D, beta_D, alpha_E, beta_E)
            - total)





def probability_E_beats_ABC2(alpha_A, beta_A, alpha_B, beta_B, alpha_C, beta_C, alpha_D, beta_D, alpha_E, beta_E):
    total = 0.0
    for i in range(0, (alpha_A)):
        for j in range(0, (alpha_B)):
            for k in range(0, (alpha_C)):
                for l in range(0, (alpha_D)):

                    total += math.exp(bln(alpha_E + i + j + k + l, beta_A + beta_B + beta_C + beta_D + beta_E)
                                 - math.log(beta_A + i) - math.log(beta_B + j) - math.log(beta_C + k) - math.log(beta_D + l)
                                 - bln(1 + i, beta_A) - bln(1 + j, beta_B) - bln(1 + k, beta_C) - bln(1 + l, beta_D)
                                 - bln(alpha_E, beta_E))



    return (1 - probability_B_beats_A(alpha_A, beta_A, alpha_E, beta_E)
            - probability_B_beats_A(alpha_B, beta_B, alpha_E, beta_E)
            - probability_B_beats_A(alpha_C, beta_C, alpha_E, beta_E)
            - probability_B_beats_A(alpha_D, beta_D, alpha_E, beta_E)
            + probability_C_beats_A_and_B(alpha_A, beta_A, alpha_B, beta_B, alpha_E, beta_E)
            + probability_C_beats_A_and_B(alpha_A, beta_A, alpha_C, beta_C, alpha_E, beta_E)
            + probability_C_beats_A_and_B(alpha_C, beta_C, alpha_D, beta_D, alpha_E, beta_E)
            + probability_C_beats_A_and_B(alpha_B, beta_B, alpha_C, beta_C, alpha_E, beta_E)
            + probability_C_beats_A_and_B(alpha_D, beta_D, alpha_C, beta_C, alpha_E, beta_E)
            + probability_C_beats_A_and_B(alpha_B, beta_B, alpha_D, beta_D, alpha_E, beta_E)
            + probability_C_beats_A_and_B(alpha_A, beta_A, alpha_D, beta_D, alpha_E, beta_E)
            + probability_D_beats_ABC(alpha_A, beta_A, alpha_B, beta_B, alpha_C, beta_C, alpha_E, beta_E)
            + probability_D_beats_ABC(alpha_A, beta_A, alpha_B, beta_B, alpha_D, beta_D, alpha_E, beta_E)
            + probability_D_beats_ABC(alpha_A, beta_A, alpha_D, beta_D, alpha_C, beta_C, alpha_E, beta_E)
            + probability_D_beats_ABC(alpha_B, beta_B, alpha_C, beta_C, alpha_D, beta_D, alpha_E, beta_E)
            - total)


def probability_F_beats_ABC(alpha_A, beta_A, alpha_B, beta_B, alpha_C, beta_C, alpha_D, beta_D, alpha_E, beta_E, alpha_F, beta_F):
    total = 0.0
    for i in range(0, (alpha_A)):
        for j in range(0, (alpha_B)):
            for k in range(0, (alpha_C)):
                for l in range(0, (alpha_D)):
                    for m in range(0, (alpha_E)):

                        total += math.exp(bln(alpha_F + i + j + k + l + m, beta_A + beta_B + beta_C + beta_D + beta_E + beta_F)
                                     - math.log(beta_A + i) - math.log(beta_B + j) - math.log(beta_C + k) - math.log(beta_D + l) - math.log(beta_E + m)
                                     - bln(1 + i, beta_A) - bln(1 + j, beta_B) - bln(1 + k, beta_C) - bln(1 + l, beta_D) - bln(1 + m, beta_E)
                                     - bln(alpha_F, beta_F))



    return (1 - probability_B_beats_A(alpha_A, beta_A, alpha_F, beta_F)
            - probability_B_beats_A(alpha_B, beta_B, alpha_F, beta_F)
            - probability_B_beats_A(alpha_C, beta_C, alpha_F, beta_F)
            - probability_B_beats_A(alpha_D, beta_D, alpha_F, beta_F)
            - probability_B_beats_A(alpha_E, beta_E, alpha_F, beta_F)
            + probability_C_beats_A_and_B(alpha_A, beta_A, alpha_B, beta_B, alpha_F, beta_F)
            + probability_C_beats_A_and_B(alpha_A, beta_A, alpha_C, beta_C, alpha_F, beta_F)
            + probability_C_beats_A_and_B(alpha_A, beta_A, alpha_D, beta_D, alpha_F, beta_F)
            + probability_C_beats_A_and_B(alpha_A, beta_A, alpha_E, beta_E, alpha_F, beta_F)
            + probability_C_beats_A_and_B(alpha_B, beta_B, alpha_C, beta_C, alpha_F, beta_F)
            + probability_C_beats_A_and_B(alpha_B, beta_B, alpha_D, beta_D, alpha_F, beta_F)
            + probability_C_beats_A_and_B(alpha_B, beta_B, alpha_E, beta_E, alpha_F, beta_F)
            + probability_C_beats_A_and_B(alpha_D, beta_D, alpha_C, beta_C, alpha_F, beta_F)
            + probability_C_beats_A_and_B(alpha_E, beta_E, alpha_C, beta_C, alpha_F, beta_F)
            + probability_C_beats_A_and_B(alpha_E, beta_E, alpha_D, beta_D, alpha_F, beta_F)
            + probability_D_beats_ABC(alpha_A, beta_A, alpha_B, beta_B, alpha_C, beta_C, alpha_F, beta_F)
            + probability_D_beats_ABC(alpha_A, beta_A, alpha_B, beta_B, alpha_D, beta_D, alpha_F, beta_F)
            + probability_D_beats_ABC(alpha_A, beta_A, alpha_B, beta_B, alpha_E, beta_E, alpha_F, beta_F)
            + probability_D_beats_ABC(alpha_A, beta_A, alpha_C, beta_C, alpha_E, beta_E, alpha_F, beta_F)
            + probability_D_beats_ABC(alpha_A, beta_A, alpha_D, beta_D, alpha_C, beta_C, alpha_F, beta_F)
            + probability_D_beats_ABC(alpha_B, beta_B, alpha_C, beta_C, alpha_D, beta_D, alpha_F, beta_F)
            + probability_D_beats_ABC(alpha_B, beta_B, alpha_C, beta_C, alpha_E, beta_E, alpha_F, beta_F)
            + probability_E_beats_alpha_and_B(alpha_A, beta_A, alpha_B, beta_B, alpha_C, beta_C, alpha_E, beta_E, alpha_F, beta_F)
            + probability_E_beats_alpha_and_B(alpha_A, beta_A, alpha_B, beta_B, alpha_C, beta_C, alpha_D, beta_D, alpha_F, beta_F)
            + probability_E_beats_alpha_and_B(alpha_A, beta_A, alpha_B, beta_B, alpha_D, beta_D, alpha_E, beta_E, alpha_F, beta_F)
            + probability_E_beats_alpha_and_B(alpha_A, beta_A, alpha_D, beta_D, alpha_C, beta_C, alpha_E, beta_E, alpha_F, beta_F)
            + probability_E_beats_alpha_and_B(alpha_B, beta_B, alpha_C, beta_C, alpha_D, beta_D, alpha_E, beta_E, alpha_F, beta_F)
            - total)


