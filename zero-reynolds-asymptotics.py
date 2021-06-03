from sympy import *

def orrSommerfeld():
    eqn = diff(phi, z, 4) - 2 * k ** 2 * diff(phi, z, 2) + k ** 4
    bc = zeros(4, 1)
    bc[0] = diff(phi, z, 2) * (1 - c) + (k ** 2 * (1 - c) + 2) * phi
    bc[1] = ((I * diff(phi, z, 3) - 3 * I * k ** 2 * diff(phi, z)) * (1 - c) -
             (2 * k * cotbeta + S * k ** 3) * phi)
    bc[2] = diff(phi, z) * c + 2 * phi
    bc[3] = (diff(phi, z, 3) * c +
             (I * k) * (- I * k * c * AD +
              k ** 2 * AT + k ** 4 * AB - 2 * cotbeta -
              2 * k * I) * phi)
    return eqn, bc

def subsSys(eqn, bc, exp):
    eqn = eqn.subs(exp)
    for i in range(4):
        bc[i] = bc[i].subs(exp)
        
    return eqn, bc

def solveAsymptotics(order, variablesToSolveFor):
    bc0 = zeros(4, 1)
    if order == 0:
        diffeq = eqn.subs(k, order)
        for i in range(4):
            bc0[i] = bc[i].subs(k, 0)
    else:
        diffeq = eqn.coeff(k ** order)
        for i in range(4):
            bc0[i] = bc[i].coeff(k ** order)
    
    sol = dsolve(diffeq).rhs

    for i in range(4):
        bc0[i] = bc0[i].subs(phi0, sol).doit()
    
    for i in range(2):
        bc0[i] = simplify(bc0[i].subs(z, 1))
    
    for i in range(2, 4):
        bc0[i] = simplify(bc0[i].subs(z, 0))
    
    coeffSol = linsolve(bc0[0:2], variablesToSolveFor).args[0]
    
    sol = sol.subs([(variablesToSolveFor[i], coeffSol[i]) for i in range(4)])
    
    return sol


if __name__=="__main__":
    z, c, c0, c1, c2 = symbols('z c c0 c1 c2')
    phi = Function('phi')(z)
    phi0 = Function('phi0')(z)
    phi1 = Function('phi1')(z)
    phi2 = Function('phi2')(z)
    phi3 = Function('phi3')(z)
    
    A, B, C, D = symbols('A B C D')
    C1, C2, C3, C4 = symbols('C1, C2, C3, C4')
    B1, B2, B3, B4 = symbols('B1, B2, B3, B4')
    k, cotbeta = symbols('k cotbeta')
    AD, AT, AB, S = symbols('AD AT AB S')
    
    eqn, bc = orrSommerfeld()
    
    eqn = simplify(eqn.subs(phi, phi0 + k * phi1 + k ** 2 * phi2).subs(c, c0 + k * c1 + k ** 2 * c2))
    for i in range(4):
        bc[i] = simplify(bc[i].subs(phi, phi0 + k * phi1 + k ** 2 * phi2).subs(c, c0 + k * c1 + k ** 2 * c2))
    
    sol = solveAsymptotics(0, [C1, C2, C3, C4])

    eqn = expand(eqn.subs(phi0, sol).doit())
    for i in range(4):
        bc[i] = expand(bc[i].subs(phi0, sol).doit())


    order = 1

    bc1 = zeros(4, 1)
    if order == 0:
        diffeq = eqn.subs(k, order)
        for i in range(4):
            bc1[i] = bc[i].subs(k, 0)
    else:
        diffeq = eqn.coeff(k ** order)
        for i in range(4):
            bc1[i] = bc[i].coeff(k ** order)

    sol = dsolve(diffeq).rhs

    sol = ode.constant_renumber(sol, newconstants=symbols('B1:4'))

    for i in range(4):
        bc1[i] = bc1[i].subs(phi1, sol).doit()
    
    for i in range(2):
        bc1[i] = simplify(bc1[i].subs(z, 1))
    
    for i in range(2, 4):
        bc1[i] = simplify(bc1[i].subs(z, 0))

    print(bc1)
    print(linear_eq_to_matrix(bc1, B1, B3, B4))
    print(nonlinsolve(bc1, B1, B3, B4))

    #coeffSol = linsolve(bc1, B1, B2, B3, B4).args[0]
    #
    #print(coeffSol)
    #with open("out","w") as f:
    #    print(limit1,file=f)
