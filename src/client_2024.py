from scipy.special import betaln
import numpy as np


def probability_B_beats_A(alpha_A, beta_A, alpha_B, beta_B):

    i = np.arange(0, alpha_B-1)
    total = np.sum(np.exp(betaln(alpha_A + i, beta_A + beta_B) -
                          np.log(beta_B + i) -
                          betaln(1 + i, beta_B) -
                          betaln(alpha_A, beta_A)))

    return (total)


def probability_C_beats_AB(alpha_A, beta_A, alpha_B, beta_B, alpha_C, beta_C):

    i = np.arange(0, alpha_A-1)
    j = np.arange(0, alpha_B-1)
    total = np.sum(np.exp(betaln(alpha_C + i + j, beta_A + beta_B + beta_C) -
                          np.log(beta_A + i) - np.log(beta_B + j) -
                          betaln(1 + i, beta_A) - betaln(1 + j, beta_B) -
                          betaln(alpha_C, beta_C)))

    return (1 - probability_B_beats_A(alpha_C, beta_C, alpha_A, beta_A) -
            probability_B_beats_A(alpha_C, beta_C, alpha_B, beta_B) + total)


def probability_D_beats_ABC(alpha_A, beta_A, alpha_B, beta_B, alpha_C, beta_C, alpha_D, beta_D):

    i = np.arange(0, alpha_A-1)
    j = np.arange(0, alpha_B-1)
    k = np.arange(0, alpha_C-1)
    total = np.sum(np.exp(betaln(alpha_D + i + j + k, beta_A + beta_B + beta_C + beta_D) -
                          np.log(beta_A + i) - np.log(beta_B + j) - np.log(beta_C + k) -
                          betaln(1 + i, beta_A) - betaln(1 + j, beta_B) - betaln(1 + k, beta_C) -
                          betaln(alpha_D, beta_D)))

    return (1 - probability_B_beats_A(alpha_A, beta_A, alpha_D, beta_D) -
            probability_B_beats_A(alpha_B, beta_B, alpha_D, beta_D) -
            probability_B_beats_A(alpha_C, beta_C, alpha_D, beta_D) +
            probability_C_beats_AB(alpha_A, beta_A, alpha_B, beta_B, alpha_D, beta_D) +
            probability_C_beats_AB(alpha_A, beta_A, alpha_C, beta_C, alpha_D, beta_D) +
            probability_C_beats_AB(alpha_B, beta_B, alpha_C, beta_C, alpha_D, beta_D) - total)


def probability_E_beats_ABC(alpha_A, beta_A, alpha_B, beta_B, alpha_C, beta_C, alpha_D, beta_D, alpha_E, beta_E):

    i = np.arange(0, alpha_A-1)
    j = np.arange(0, alpha_B-1)
    k = np.arange(0, alpha_C-1)
    l = np.arange(0, alpha_D-1)
    total = np.sum(np.exp(betaln(alpha_E + i + j + k + l, beta_A + beta_B + beta_C + beta_D + beta_E) -
                          np.log(beta_A + i) - np.log(beta_B + j) -
                          np.log(beta_C + k) - np.log(beta_D + l) -
                          betaln(1 + i, beta_A) - betaln(1 + j, beta_B) -
                          betaln(1 + k, beta_C) - betaln(1 + l, beta_D) -
                          betaln(alpha_E, beta_E)))

    return (1 - probability_B_beats_A(alpha_A, beta_A, alpha_E, beta_E) -
            probability_B_beats_A(alpha_B, beta_B, alpha_E, beta_E) -
            probability_B_beats_A(alpha_C, beta_C, alpha_E, beta_E) +
            probability_C_beats_AB(alpha_A, beta_A, alpha_B, beta_B, alpha_E, beta_E) +
            probability_C_beats_AB(alpha_A, beta_A, alpha_C, beta_C, alpha_E, beta_E) +
            probability_C_beats_AB(alpha_B, beta_B, alpha_C, beta_C, alpha_E, beta_E) +
            probability_D_beats_ABC(alpha_A, beta_A, alpha_B, beta_B, alpha_C, beta_C, alpha_E, beta_E) +
            probability_D_beats_ABC(alpha_A, beta_A, alpha_B, beta_B, alpha_D, beta_D, alpha_E, beta_E) +
            probability_D_beats_ABC(alpha_A, beta_A, alpha_D, beta_D, alpha_C, beta_C, alpha_E, beta_E) +
            probability_D_beats_ABC(alpha_B, beta_B, alpha_C, beta_C, alpha_D, beta_D, alpha_E, beta_E) - total)


def probability_E_beats_AB(alpha_A, beta_A, alpha_B, beta_B, alpha_C, beta_C, alpha_D, beta_D, alpha_E, beta_E):

    i = np.arange(0, alpha_A-1)
    j = np.arange(0, alpha_B-1)
    k = np.arange(0, alpha_C-1)
    l = np.arange(0, alpha_D-1)
    total = np.sum(np.exp(betaln(alpha_E + i + j + k + l, beta_A + beta_B + beta_C + beta_D + beta_E) -
                          np.log(beta_A + i) - np.log(beta_B + j) -
                          np.log(beta_C + k) - np.log(beta_D + l) -
                          betaln(1 + i, beta_A) - betaln(1 + j, beta_B) -
                          betaln(1 + k, beta_C) - betaln(1 + l, beta_D) -
                          betaln(alpha_E, beta_E)))

    return (1 - probability_B_beats_A(alpha_A, beta_A, alpha_E, beta_E) -
            probability_B_beats_A(alpha_B, beta_B, alpha_E, beta_E) -
            probability_B_beats_A(alpha_C, beta_C, alpha_E, beta_E) -
            probability_B_beats_A(alpha_D, beta_D, alpha_E, beta_E) +
            probability_C_beats_AB(alpha_A, beta_A, alpha_B, beta_B, alpha_E, beta_E) +
            probability_C_beats_AB(alpha_A, beta_A, alpha_C, beta_C, alpha_E, beta_E) +
            probability_C_beats_AB(alpha_C, beta_C, alpha_D, beta_D, alpha_E, beta_E) +
            probability_C_beats_AB(alpha_B, beta_B, alpha_C, beta_C, alpha_E, beta_E) +
            probability_C_beats_AB(alpha_D, beta_D, alpha_C, beta_C, alpha_E, beta_E) +
            probability_C_beats_AB(alpha_B, beta_B, alpha_D, beta_D, alpha_E, beta_E) +
            probability_C_beats_AB(alpha_A, beta_A, alpha_D, beta_D, alpha_E, beta_E) +
            probability_D_beats_ABC(alpha_A, beta_A, alpha_B, beta_B, alpha_C, beta_C, alpha_E, beta_E) +
            probability_D_beats_ABC(alpha_A, beta_A, alpha_B, beta_B, alpha_D, beta_D, alpha_E, beta_E) +
            probability_D_beats_ABC(alpha_A, beta_A, alpha_D, beta_D, alpha_C, beta_C, alpha_E, beta_E) +
            probability_D_beats_ABC(alpha_B, beta_B, alpha_C, beta_C, alpha_D, beta_D, alpha_E, beta_E) - total)


def probability_F_beats_ABC(alpha_A, beta_A, alpha_B, beta_B, alpha_C, beta_C, alpha_D, beta_D, alpha_E, beta_E, alpha_F, beta_F):

    i = np.arange(0, alpha_A-1)
    j = np.arange(0, alpha_B-1)
    k = np.arange(0, alpha_C-1)
    l = np.arange(0, alpha_D-1)
    m = np.arange(0, alpha_E-1)
    total = np.sum(np.exp(betaln(alpha_F + i + j + k + l + m, beta_A + beta_B + beta_C + beta_D + beta_E + beta_F) -
                          np.log(beta_A + i) - np.log(beta_B + j) - np.log(beta_C + k) -
                          np.log(beta_D + l) - np.log(beta_E + m) -
                          betaln(1 + i, beta_A) - betaln(1 + j, beta_B) - betaln(1 + k, beta_C) -
                          betaln(1 + l, beta_D) - betaln(1 + m, beta_E) -
                          betaln(alpha_F, beta_F)))

    return (1 - probability_B_beats_A(alpha_A, beta_A, alpha_F, beta_F) -
            probability_B_beats_A(alpha_B, beta_B, alpha_F, beta_F) -
            probability_B_beats_A(alpha_C, beta_C, alpha_F, beta_F) -
            probability_B_beats_A(alpha_D, beta_D, alpha_F, beta_F) -
            probability_B_beats_A(alpha_E, beta_E, alpha_F, beta_F) +
            probability_C_beats_AB(alpha_A, beta_A, alpha_B, beta_B, alpha_F, beta_F) +
            probability_C_beats_AB(alpha_A, beta_A, alpha_C, beta_C, alpha_F, beta_F) +
            probability_C_beats_AB(alpha_A, beta_A, alpha_D, beta_D, alpha_F, beta_F) +
            probability_C_beats_AB(alpha_A, beta_A, alpha_E, beta_E, alpha_F, beta_F) +
            probability_C_beats_AB(alpha_B, beta_B, alpha_C, beta_C, alpha_F, beta_F) +
            probability_C_beats_AB(alpha_B, beta_B, alpha_D, beta_D, alpha_F, beta_F) +
            probability_C_beats_AB(alpha_B, beta_B, alpha_E, beta_E, alpha_F, beta_F) +
            probability_C_beats_AB(alpha_D, beta_D, alpha_C, beta_C, alpha_F, beta_F) +
            probability_C_beats_AB(alpha_E, beta_E, alpha_C, beta_C, alpha_F, beta_F) +
            probability_C_beats_AB(alpha_E, beta_E, alpha_D, beta_D, alpha_F, beta_F) +
            probability_D_beats_ABC(alpha_A, beta_A, alpha_B, beta_B, alpha_C, beta_C, alpha_F, beta_F) +
            probability_D_beats_ABC(alpha_A, beta_A, alpha_B, beta_B, alpha_D, beta_D, alpha_F, beta_F) +
            probability_D_beats_ABC(alpha_A, beta_A, alpha_B, beta_B, alpha_E, beta_E, alpha_F, beta_F) +
            probability_D_beats_ABC(alpha_A, beta_A, alpha_C, beta_C, alpha_E, beta_E, alpha_F, beta_F) +
            probability_D_beats_ABC(alpha_A, beta_A, alpha_D, beta_D, alpha_C, beta_C, alpha_F, beta_F) +
            probability_D_beats_ABC(alpha_B, beta_B, alpha_C, beta_C, alpha_D, beta_D, alpha_F, beta_F) +
            probability_D_beats_ABC(alpha_B, beta_B, alpha_C, beta_C, alpha_E, beta_E, alpha_F, beta_F) +
            probability_E_beats_AB(alpha_A, beta_A, alpha_B, beta_B, alpha_C, beta_C, alpha_E, beta_E, alpha_F, beta_F) +
            probability_E_beats_AB(alpha_A, beta_A, alpha_B, beta_B, alpha_C, beta_C, alpha_D, beta_D, alpha_F, beta_F) +
            probability_E_beats_AB(alpha_A, beta_A, alpha_B, beta_B, alpha_D, beta_D, alpha_E, beta_E, alpha_F, beta_F) +
            probability_E_beats_AB(alpha_A, beta_A, alpha_D, beta_D, alpha_C, beta_C, alpha_E, beta_E, alpha_F, beta_F) +
            probability_E_beats_AB(alpha_B, beta_B, alpha_C, beta_C, alpha_D, beta_D, alpha_E, beta_E, alpha_F, beta_F) - total)
