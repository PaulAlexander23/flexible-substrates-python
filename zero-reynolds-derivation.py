from sympy import *

phi, z, c = symbols('phi z c')
A, B, C, D = symbols('A B C D')
alpha, cot_beta = symbols('alpha cot_beta')
A_I, A_D, A_T, A_B, A_K, S = symbols('A_I A_D A_T A_B A_K S')

phi = (A * cosh(alpha * z) + B * sinh(alpha * z) +
       C * alpha * z * cosh(alpha * z) + D * alpha * z * sinh(alpha * z))

bc = zeros(4, 1)
bc[0] = diff(phi, z, 2) * (1 - c) + (alpha ** 2 * (1 - c) + 2) * phi
bc[1] = ((I * diff(phi, z, 3) - 3 * I * alpha ** 2 * diff(phi, z)) * (1 - c) -
         (2 * alpha * cot_beta + S * alpha ** 3) * phi)
bc[2] = diff(phi, z) * c + 2 * phi
bc[3] = (diff(phi, z, 3) * c / (I * alpha) +
         (- alpha ** 2 * c ** 2 * A_I - I * alpha * c * A_D +
          alpha ** 2 * A_T + alpha ** 4 * A_B + A_K - 2 * cot_beta -
          2 * alpha * I) * phi)

for i in range(2):
    bc[i] = simplify(bc[i].subs(z, 1))

for i in range(2, 4):
    bc[i] = simplify(bc[i].subs(z, 0))

mat = Matrix([[expand(z).coeff(t) for t in [A, B, C, D]] for z in bc])
vec = [z.subs([(A, 0), (B, 0), (C, 0), (D, 0)]) for z in bc]

eqn = mat.det()
myc = solve(eqn, c)

print(myc)
