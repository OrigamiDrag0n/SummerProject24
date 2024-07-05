## OrigamiDrag0n 2024, 05/07/2024

## It has been shown that, for elements of the Thompson Group F,
## it is possible to conjugate by the Minkowski question mark
## function and create an embedding of F into the collection
## of smooth automorphisms of [0, 1]. In this piece of code, I
## demonstrate this.

# Modules

import matplotlib.pyplot as plt
import numpy as np
from math import floor

# Functions

def question(x, n):

    ## Evaluates the Minkowski question mark function on [0, 1]

    ## The Minkowski function ?: [0, 1] -> [0, 1] is defined to be the
    ## function which converts a continued fraction into run-length encoding
    ## e.g, ?([0: a_1, a_2, ...]) = 0.000...011...1... with a_1 - 1 0s, then
    ## a_2 1s, then a_3 0s, etc. It satisfies the recurrence relations
    ## ?(x/(1 + x)) = ?(x)/2, ?(1 - x) = 1 - ?(x), for x in [0, 1].

    ## Inputs:
    ## x - value in [0, 1]
    ## n - number of iterations

    ## Outputs:
    ## z - the inverse Minkowski question mark function evaluated on x
    ##     to n digits. 

    if x == 0 or n == 0: # Edge cases
        return x
    else:                # Applying recurrence relation
        if x >= 1/2: 
            y = question(2 - 1/x, n - 1)
            return (y + 1)/2
        else:
            y = question(1/(1 - x) - 1, n - 1)
            return y/2

def invquestion(x, n): 

    ## Evaluates the inverse Minkowski question mark function on [0, 1]
    
    ## Inputs:
    ## x - value in [0, 1]
    ## n - number of binary digits
    
    ## Outputs:
    ## z - continued fraction from run-length encoding of x in binary
    ##     up to the nth digit

    if x == 0 or n == 0: # Edge cases
        return x
    else:                # Applying recurrence relation
        if x >= 1/2:
            y = invquestion(2*x - 1, n - 1)
            return 1/(2 - y)
        else:
            y = invquestion(2*x, n - 1)
            return y/(1 + y)

def to_bin(x, n):

    ## Converts a number to a binary string

    ## Inputs:
    ## x - value in [0, 1]
    ## n - number of binary digits

    ## Outputs:
    ## string - the first n digits of x following the binary point

    if n == 0:
        return ''
    else:
        if x >= 1/2:
            return '1' + to_bin(2*x - 1, n - 1)
        else:
            return '0' + to_bin(2*x, n - 1)

def to_float(string):

    ## Converts a binary string to a number

    ## Inputs:
    ## string - string in '0' and '1'

    ## Outputs:
    ## x - value in [0, 1]

    if string == '':
        return 0
    else:
        return (int(string[0]) + to_float(string[1:]))/2

def treemap(sigma):

    ## Tree maps used in the construction of the Thompson group.

    ## A tree map takes some prefix sets A and B such that every 
    ## infinite binary string has a unique prefix a_i in A and a
    ## unique prefix b_i, and some bijection \sigma: A -> B. The
    ## resulting map then substitutes a the prefix a in A for
    ## \sigma(a) in B. If the tree map preserves lexicographic order
    ## then we can also consider it acting on binary expansions in [0, 1]
    ## in a well-defined way (producing the Thompson group F). 
    
    ## Inputs:
    ## sigma - a dictionary from prefix sets A to B. 

    ## Outputs:
    ## f - the tree map from [0, 1] to [0, 1].

    def f(x):

        n = max(len(s) for s in sigma)
        string = to_bin(x, n)
        for string_0 in sigma:
            if string[:len(string_0)] == string_0:
                y = (x - to_float(string_0)) * (2**len(string_0))
                return to_float(sigma[string_0]) + y * (2**(-len(sigma[string_0])))
    return f

## Main routine

if __name__ == '__main__':

    ## Constants

    n = 10  # Binary depth
    m = 20  # Iterations of the (inverse) Minkowski function

    q = lambda x: question(x, m)  # Minkowski function 
    invq = lambda x: invquestion(x, m)  # Inverse Minkowski function

    sigma1 = {'00': '0', '01': '10', '1': '11'} # Dictionary for element 'a' of the Thompson group
    sigma2 = {'000': '00', '001': '010', '01': '011', '1': '1'} # Dictionaty for element 'b' of the Thompson group
    
    f1 = treemap(sigma1) # The corresponding tree maps
    f2 = treemap(sigma2) 

    npq = np.vectorize(q) # Vectorising functions
    npinvq = np.vectorize(invq)
    npf1 = np.vectorize(f1)
    npf2 = np.vectorize(f2)

    x = np.arange(0, 1 + 1/(2**n), 1/(2**n)) # Base range
    
    z1 = npf1(x) # Tree maps
    z2 = npf2(x)
    
    y1 = npinvq(npf1(npq(x))) # Conjugations of the tree maps
    y2 = npinvq(npf2(npq(x)))

    ## Plotting tree maps and conjugations of tree maps
    
    plt.plot(x, z1)
    plt.plot(x, z2)
    plt.plot(x, y1)
    plt.plot(x, y2)
    
    #plt.plot(x, npq(x)) # Question mark function
    
    plt.gca().set_aspect('equal')
    plt.show()
