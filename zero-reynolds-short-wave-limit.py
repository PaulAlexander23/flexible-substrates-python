from sympy import *

phi, z, c = symbols('phi z c')
A, B, C, D = symbols('A B C D')
k, cotbeta = symbols('k cotbeta')
AD, AT, AB, AK, S = symbols('AD AT AB AK S')

phi = (A * cosh(k * z) + B * sinh(k * z) +
       C * k * z * cosh(k * z) + D * k * z * sinh(k * z))

bc = zeros(4, 1)
bc[0] = diff(phi, z, 2) * (1 - c) + (k ** 2 * (1 - c) + 2) * phi
bc[1] = ((I * diff(phi, z, 3) - 3 * I * k ** 2 * diff(phi, z)) * (1 - c) -
         (2 * k * cotbeta + S * k ** 3) * phi)
bc[2] = diff(phi, z) * c + 2 * phi
bc[3] = (diff(phi, z, 3) * c / (I * k) +
         (- I * k * c * AD +
          k ** 2 * AT + k ** 4 * AB + AK - 2 * cotbeta -
          2 * k * I) * phi)

for i in range(2):
    bc[i] = simplify(bc[i].subs(z, 1))

for i in range(2, 4):
    bc[i] = simplify(bc[i].subs(z, 0))

mat = Matrix([[expand(z).coeff(t) for t in [A, B, C, D]] for z in bc])
vec = [z.subs([(A, 0), (B, 0), (C, 0), (D, 0)]) for z in bc]

eqn = mat.det()
myc = solve(eqn, c)

myc = myc[2:5]

limit1 = [limit(z, k, oo) for z in myc]
limit2 = [limit(z.subs(AB, 0), k, oo) for z in myc]
limit3 = [limit(z.subs([(AB, 0), (AT, 0), (S, 0)]), k, oo) for z in myc]

with open("out","w") as f:
    print(limit1,file=f)
    print(limit2,file=f)
    print(limit3,file=f)
