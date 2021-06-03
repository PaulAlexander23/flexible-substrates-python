from sympy import *
#restart;
#with(LinearAlgebra):
#with(PDEtools):
#with(plots):
#with(CodeGeneration):

phi, z, c = symbols('phi z c')
phi0, phi1, phi2 = functions('phi0 phi1 phi2')
c1, c2 = symbols('c1 c2')
k, cotbeta, R, S = symbols('k cotbeta R S')
AI, AD, AT, AB, AK = symbols('AI AD AT AB AK')
#declare(phi,phi__0,phi__1,phi__2):
#constants:=constants,R,A__I,A__D,A__T,A__B,A__I:

orr_sommerfeld = diff(phi,z,4) - 2 * k**2 * diff(phi,z,2) + k**4 * phi - I * k * R * (z * (2 - z) - c) * (diff(phi,z,2) - k**2 * phi) - 2 * k * I * R * phi
bc = zeros(4, 1)
bc[0] = diff(phi,z,2) * (1 - c) + (k**2 * (1 - c) + 2) * phi
bc[1] = I * diff(phi,z,3) * (1 - c) + (k * R * (1 - c) - 3 * I * k**2) * (1 - c) * diff(phi,z) - (2 * k * cot(beta) + k * S) * phi
bc[2] = c * diff(phi,z) + 2 * phi
bc[3] = I * k * (-c**2 * A__I - I * c * A__D + A__T + A__B + A__K - 2 * cot(beta)) * phi + c * diff(phi,z,3) + 2 * k**2 * phi
#OrrSommerfeld:= diff(phi,z,4) - 2 * k**2 * diff(phi,z,2) + k**4 * phi - I * k * R * (z * (2 - z) - c) * (diff(phi,z,2) - k**2 * phi) - 2 * k * I * R * phi:
#BC1:= diff(phi,z,2) * (1 - c) + (k**2 * (1 - c) + 2) * phi:
#BC2:= I * diff(phi,z,3) * (1 - c) + (k * R * (1 - c) - 3 * I * k**2) * (1 - c) * diff(phi,z) - (2 * k * cot(beta) + k * S) * phi:
#BC3:= c * diff(phi,z) + 2 * phi:
#BC4:= I * k * (-c**2 * A__I - I * c * A__D + A__T + A__B + A__K - 2 * cot(beta)) * phi + c * diff(phi,z,3) + 2 * k**2 * phi:

expansion = [(phi , phi0 + k * phi1 + k**2 * phi2), (c, 2 + k * c1 + k**2 * c2)]
orr_sommerfeld_expanded = orr_sommerfeld.subs(expansion)
bc_expanded = [z.subs(expansion) for z in bc]
#base_system:= [OrrSommerfeld, BC1, BC2, BC3, BC4]:
#expanded_system:= map(x-&gt;taylor(subs([phi = phi__0 + k * phi__1 + k**2 * phi__2, c = 2 + k * c__1 + k**2 * c__2],x),k,3),base_system):

orr_sommerfeld_collected = collect(orr_sommerfeld_expanded, k, evaluate=False)
bc_collected = [collect(z, k, evaluate=False) for z in bc_expanded]

zeroth_order = orr_sommerfeld_collected[1]
zeroth_order_bc = [z[1] for z in bc_collected]
for i in range(2):
    zeroth_order_bc[i] = zeroth_order_bc[i].subs(z, 1)
for i in range(2,4):
    zeroth_order_bc[i] = zeroth_order_bc[i].subs(z, 0)

#zeroth_order:= map(x-&gt;coeff(x,k,0),expanded_system):
#zeroth_order:= [zeroth_order[1],eval(zeroth_order[2],z=1),eval(zeroth_order[3],z=1),eval(zeroth_order[4],z=0),eval(zeroth_order[5],z=0)]:
#print~(simplify(zeroth_order))[];
#sol1:= dsolve(zeroth_order)

first_order = orr_sommerfeld_collected[k]
first_order_bc = [z[k] for z in bc_collected]
#first_order:= eval(map(x-&gt;coeff(x,k,1),expanded_system),[phi__0 = z**2 + a__3 * (z - 1), c__0 = 2]): #a__1 = 1
#first_order:= [first_order[1]=0,eval(first_order[2],z=1)=0,eval(first_order[3],z=1)=0,eval(first_order[4],z=0)=0,eval(first_order[5],z=0)=0,diff(c__1,z)=0,diff(a__3,z)=0]:
#print~(simplify(first_order))[];
#sol2:= dsolve(first_order);
#assign(sol2[2]);
#ans:=simplify(taylor(eval(c__1,[A__T=k**2*A__T,S=k**2*S,A__D=k*A__D,A__B=k**4*A__B,A__I=k**2*A__I]),k,7));
##expand(Re(ans));
##Matlab(ans)
##simplify(numer(ans));
##latex(simplify(numer(ans)));

second_order = orr_sommerfeld_collected[k**2]
second_order_bc = [z[k**2] for z in bc_collected]
#second_order:= eval(map(x-&gt;coeff(x,k,2),expanded_system),[phi__0 = z**2 + a__3 * (z - 1), c__0 = 2,phi__1 = 2 * Re * i * (2 + a__3) * (z**5/(5!) - z**4/(4!)) + b__1 * z**3 + b__2 * z**2 + b__3 * z + b__4]):
#second_order:= [second_order[1],eval(second_order[2],z=1),eval(second_order[3],z=1),eval(second_order[4],z=0),eval(second_order[5],z=0)]:
#print~(simplify(second_order))[];
#sol3:= dsolve(second_order);
