import numpy as np
import linalg
import math

R = 8.14
gamma = 1.17
mu = 0.021
Q = 5.02e6
Ea = 113e3
A = 6.85e9

eps = 1e-5
delta = 2e-7
x = 0

## начальные условия параметров перед фронтом
p0 = 1e5
r0 = 0.873  ## плотность
T0 = p0*mu/R/r0

r = []
u = []
p = []
Z = []
T = []

c = []

## параметры фон Неймана в качестве ну при х = 0

Dcj = math.sqrt((gamma * gamma - 1)*Q/2 + gamma*p0/r0) + math.sqrt((gamma*gamma - 1)*Q/2)

p.insert(0, (2*r0*Dcj*Dcj/(gamma + 1) - (gamma -1)*p0/(gamma + 1)))
u.insert(0, (Dcj - (p[0] - p0)/math.sqrt(r0*( (gamma + 1)*p[0]/2 + (gamma - 1)*p0/2) ) ) )
r.insert(0, r0*Dcj/u[0])
Z.insert(0, 1.0)
T.insert(0, p[0] * mu / r[0] / R)

c.insert(0, math.sqrt(gamma * p[0]/r[0]))

f4 = []
f2 = []
f1 = []
f3 = []

f4.insert (0, (-1)*A*r[0]*Z[0]/u[0] * math.exp((-1)*Ea*r[0]/p[0]/mu))
f2.insert(0, (-1)*Q*u[0]*(gamma - 1)/c[0]/c[0]/(1 - u[0]*u[0]/c[0]/c[0]) * f4[0])
f1.insert(0, (-1)*r[0]/u[0] * f2[0])
f3.insert(0, (-1)*r[0]*u[0]*f2[0])

f = open('text.txt', 'w')
f.write(str(Dcj))
f.write("x" + '\t' + str(x) + '\t' + "ro" + '\t' + str(r[0]) + '\t' + "u" + '\t' + str(u[0]) + '\t' + "p" + '\t' + str(p[0]) + '\t' + "Z" + '\t' + str(Z[0]) + '\t' + "T" + '\t' + str(T[0]) + '\n')

## Решение исходной сист - Гира 1 порядка, m = 1

i = 0
k = 0

p_k0 = p_kk0 = p[0]
Z_k0 = Z_kk0 = Z[0]
u_k0 = u_kk0 = u[0]
r_k0 = r_kk0 = r[0]
f4_k0 = f4[0]
f2_k0 = f2[0]
f3_k0 = f3[0]
f1_k0 = f1[0]
c_k = c[0]

while True:

    #print(A * r_k * Z_k / u_k / f4_k)    
    r_kk = r_k0 + (-1) * delta * r_kk0 / u_kk0 * f2_k0
    u_kk = u_k0 + delta * (-1) * Q * u_kk0 * (gamma - 1) / c_k / c_k / (1 - u_kk0 * u_kk0 / c_k / c_k) * f4_k0
    p_kk = p_k0 + delta * (-1) * r_kk0 * u_kk0 * f2_k0
    Z_kk = Z_k0 + delta * (-1) * A * r_kk0 * Z_kk0 / u_kk0 * math.exp( (-1)*r_kk0*Ea/p_kk0/mu )
    temp1=abs((p_kk - p_kk0)/p_kk)
    temp2 = abs((Z_kk - Z_kk0)/Z_kk)
    temp3 = abs((u_kk - u_kk0)/u_kk)
    temp4 = abs((r_kk - r_kk0)/r_kk)
    if max([temp1, temp2, temp3, temp4]) < eps:
        break
    ##print(max([temp1, temp2, temp3, temp4]))
    ##print("max1 " + str(max(p_kk - p_k0, Z_kk - Z_k0, u_kk - u_k0, r_kk - r_k0)))

    p_kk0 = p_kk
    Z_kk0 = Z_kk
    u_kk0 = u_kk
    r_kk0 = r_kk

    c_k = math.sqrt(gamma * p_kk0/r_kk0)

    f4_k0 = (-1) * A * r_kk0 * Z_kk0 / u_kk0 * math.exp((-1)*Ea*r_kk0/p_kk0/mu)
    f2_k0 = (-1) * Q * u_kk0 * (gamma - 1) / c_k / c_k / (1 - u_kk0 * u_kk0 / c_k / c_k) * f4_k0
    f1_k0 = (-1) * r_kk0 / u_kk0 * f2_k0
    f3_k0 = (-1) * r_kk0 * u_kk0 * f2_k0


p.insert(1, p_kk)
u.insert(1, u_kk)
r.insert(1, r_kk)
Z.insert(1, Z_kk)
T.insert(1, p_kk * mu / r_kk / R)

c.insert(1, math.sqrt(gamma * p[1]/r[1]))

if c[1] == u[1]: ## обращение в ноль знаменателя во 2 уравнении
    print("stop integration")

x += delta

f.write("x" + '\t' + str(x) + '\t' + "ro" + '\t' + str(r[1]) + '\t' + "u" + '\t' + str(u[1]) + '\t' + "p" + '\t' + str(p[1]) + '\t' + "Z" + '\t' + str(Z[1]) + '\t' + "T" + '\t' + str(T[1]) + '\n')

# m = 2

p_k1 = p_kk1 = p[1]
Z_k1 = Z_kk1 = Z[1]
u_k1 = u_kk1 = u[1]
r_k1 = r_kk1 = r[1]
p_k0 = p[0]
Z_k0 = Z[0]
u_k0 = u[0]
r_k0 = r[0]

while True:
    
    c_k = math.sqrt(gamma * p_kk1 / r_kk1)
    f4_k1 = (-1) * A * r_kk1 * Z_kk1 / u_kk1 * math.exp((-1)*Ea*r_kk1/p_kk1/mu)
    f2_k1 = (-1) * Q * u_kk1 * (gamma - 1) / c_k / c_k / (1 - u_kk1 * u_kk1 / c_k / c_k) * f4_k1
    f1_k1 = (-1) * r_kk1 / u_kk1 * f2_k1
    f3_k1 = (-1) * r_kk1 * u_kk1 * f2_k1
    
    r_kk = (4 * r_k1 - r_k0)/3 + (-1) * 2/3 * delta * r_kk1 / u_kk1 * f2_k1
    u_kk = (4 * u_k1 - u_k0)/3 + 2/3 * delta * (-1) * Q * u_kk1 * (gamma - 1) / c_k / c_k / (1 - u_kk1 * u_kk1 / c_k / c_k) * f4_k1
    p_kk = (4 * p_k1 - p_k0)/3 + 2/3 * delta * (-1) * r_kk1 * u_kk1 * f2_k1
    Z_kk = (4 * Z_k1 - Z_k0)/3 + 2/3 * delta * (-1) * A * r_kk1 * Z_kk1 / u_kk1 * math.exp( (-1)*r_kk1*Ea/p_kk1/mu )

    temp1 = abs((p_kk - p_kk1)/p_kk)
    temp2 = abs((Z_kk - Z_kk1)/Z_kk)
    temp3 = abs((u_kk - u_kk1)/u_kk)
    temp4 = abs((r_kk - r_kk1)/r_kk)

    if max([temp1, temp2, temp3, temp4]) < eps:
        break
    ##print(max([temp1, temp2, temp3, temp4]))

    ##print ("max2 " + str(max(p_kk - p_k1, Z_kk - Z_k1, u_kk - u_k1, r_kk - r_k1)))

    p_kk1 = p_kk
    Z_kk1 = Z_kk
    u_kk1 = u_kk
    r_kk1 = r_kk


p.insert(2, p_kk)
u.insert(2, u_kk)
r.insert(2, r_kk)
Z.insert(2, Z_kk)
T.insert(2, p_kk * mu / r_kk / R)

c.insert(2, math.sqrt(gamma * p[2]/r[2]))

x += delta

f.write("x" + '\t' + str(x) + '\t' + "ro" + '\t' + str(r[2]) + '\t' + "u" + '\t' + str(u[2]) + '\t' + "p" + '\t' + str(p[2]) + '\t' + "Z" + '\t' + str(Z[2]) + '\t' + "T" + '\t' + str(T[2]) + '\n')

if c[2] == u[2]: ## обращение в ноль знаменателя во 2 уравнении
    print("stop integration")

# m = 3

i = 3

while True:

    p_k2 = p_kk2 = p[i-1]
    Z_k2 = Z_kk2 = Z[i-1]
    u_k2 = u_kk2 = u[i-1]
    r_k2 = r_kk2 = r[i-1]

    p_k1 = p[i-2]
    Z_k1 = Z[i-2]
    u_k1 = u[i-2]
    r_k1 = r[i-2]

    p_k0 = p[i-3]
    Z_k0 = Z[i-3]
    u_k0 = u[i-3]
    r_k0 = r[i-3]

    while True:

        c_k = math.sqrt(gamma * p_kk2 / r_kk2)
        f4_k2 = (-1) * A * r_kk2 * Z_kk2 / u_kk2 * math.exp((-1)*Ea*r_kk2/p_kk2/mu)
        f2_k2 = (-1) * Q * u_kk2 * (gamma - 1) / c_k / c_k / (1 - u_kk2 * u_kk2 / c_k / c_k) * f4_k2
        f1_k2 = (-1) * r_kk2 / u_kk2 * f2_k2
        f3_k2 = (-1) * r_kk2 * u_kk2 * f2_k2
    
        r_kk = (18 * r_k2 - 9 * r_k1 + 2 * r_k0)/11 + 6/11 * delta * (-1) * r_kk2 / u_kk2 * f2_k2
        u_kk = (18 * u_k2 - 9 * u_k1 + 2 * u_k0)/11 + 6/11 * delta * (-1) * Q * u_kk2 * (gamma - 1) / c_k / c_k / (1 - u_kk2 * u_kk2 / c_k / c_k) * f4_k2
        p_kk = (18 * p_k2 - 9 * p_k1 + 2 * p_k0)/11 + 6/11 * delta * (-1) * r_kk2 * u_kk2 * f2_k2
        Z_kk = (18 * Z_k2 - 9 * Z_k1 + 2 * Z_k0)/11 + 6/11 * delta * (-1) * A * r_kk2 * Z_kk2 / u_kk2 * math.exp( (-1)*r_kk2*Ea/p_kk2/mu )

        temp1 = abs((p_kk - p_kk2)/p_kk)
        temp2 = abs((Z_kk - Z_kk2)/Z_kk)
        temp3 = abs((u_kk - u_kk2)/u_kk)
        temp4 = abs((r_kk - r_kk2)/r_kk)
        if max([temp1, temp2, temp3, temp4]) < eps:
            break
        ##print(max([temp1, temp2, temp3, temp4]))

        ##print(str(i) + " " + str(max(p_kk - p_k2, Z_kk - Z_k2, u_kk - u_k2, r_kk - r_k2)))

        p_kk2 = p_kk
        Z_kk2 = Z_kk
        u_kk2 = u_kk
        r_kk2 = r_kk

    p.insert(i, p_kk)
    u.insert(i, u_kk)
    r.insert(i, r_kk)
    Z.insert(i, Z_kk)
    T.insert(i, p_kk * mu / r_kk / R)

    x += delta

    print(str(i) + " ro " + str(r[i]) + "  u " + str(u[i]) + "p " + str(p[i]) + " Z " + str(Z[i]) + " T " + str(T[i]))
    print(x)

    f.write("x" + '\t' + str(x) + '\t' + "ro" + '\t' + str(r[i]) + '\t' + "u" + '\t' + str(u[i]) + '\t' + "p" + '\t' + str(p[i]) + '\t' + "Z" + '\t' + str(Z[i]) + '\t' + "T" + '\t' + str(T[i]) + '\n')

    c.insert(i, math.sqrt(gamma * p[i]/r[i]))

    #if x > 0.7:
    #    break
    #else:
    #    i += 1
       
    if abs(c[i] - u[i])/c[i]<1e-2: ## обращение в ноль знаменателя во 2 уравнении
        break
    else:
        i += 1

f.close()








        

