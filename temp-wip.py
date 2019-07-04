#!/usr/bin/python3

from sympy import *





if __name__ == "__main__":
    alpha = Symbol('alpha')
    Re = Symbol('Re')
    A_I = Symbol('A_I')
    A_D = Symbol('A_D')
    A_T = Symbol('A_T')
    A_B = Symbol('A_B')
    A_K = Symbol('A_K')
    beta = Symbol('beta')
    S = Symbol('S');
    
    c = Symbol('c')
    z = Symbol('z')
    phi = Function('phi')

    base_system = [phi(z).diff(z,z,z,z) -
        2 * alpha ** 2 * phi(z).diff(z,z) +
        alpha ** 4 * phi(z) - I * alpha * Re * (z * (2 - z) - c) *
        (phi(z).diff(z,z) - alpha ** 2 * phi(z)) -
        2 * alpha * I * Re * phi(z),
        phi(z).diff(z,z) * (1 - c) + (alpha**2 * (1 - c) + 2) * phi(z),
        I * phi(z).diff(z,z,z) * (1 - c) + (alpha * Re * (1 - c) - 3 * I * alpha**2) * (1 - c) * phi(z).diff(z) - (2 * alpha * cot(beta) + alpha * S) * phi(z),
        c * phi(z).diff(z) + 2 * phi(z),
        I * alpha * (-c**2 * A_I - I * c * A_D + A_T + A_B + A_K - 2 * cot(beta)) * phi(z) + c * phi(z).diff(z,z,z) + 2 * alpha**2 * phi(z)]
    
    
    phi_0 = Function('phi_0')
    phi_1 = Function('phi_1')
    phi_2 = Function('phi_2')
    
    c_0 = Symbol('c_0')
    c_1 = Symbol('c_1')
    c_2 = Symbol('c_2')

    expanded_system = []

    for i in range(5):
        expanded_system.append(
            base_system[i].subs(
                [(phi(z),phi_0(z)+alpha*phi_1(z)+alpha**2*phi_2(z)),
                (c,c_0 + alpha * c_1 + alpha ** 2 * c_2)]))
        expanded_system[i].collect(alpha)

    nth_order_system = [[expanded_system[k].coeff(alpha,j) for k in range (5)] for j in range(3)]
    
    for i in range(5):
        print(nth_order_system[1][i])
    
    

