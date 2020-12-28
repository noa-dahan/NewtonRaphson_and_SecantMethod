#Noa Dahan, id:208375709
#Ofir Ben Zaken, id:208839712

from sympy import *
from sympy.utilities.lambdify import lambdify

x = Symbol('x')

def secant_method(f,start_point,end_point,epsi):
    x1=(start_point+end_point)/2
    x2=x1+0.1
    num_iter=0
    flag=True
    while f(x1)!=0 and flag:
        num_iter+=1
        temp=x1
        x1=x2
        x2=(temp*f(x2)-x2*f(temp))/(f(x2)-f(temp))
        if(x1-temp)<epsi:
            flag=False
    print("The number of iterations in secant_method is:", num_iter)
    return x1

def newton_raphson(f,f_prime,start_point,end_point,epsilon):
    Xr=(start_point+end_point)/2
    num_iter=0
    prev_flag=True
    while f(Xr)!=0 and prev_flag:
        fun=f(Xr)
        der=f_prime(Xr)
        Xr_prev=Xr
        Xr=Xr_prev-(fun/der)
        num_iter+=1
        if abs(Xr-Xr_prev)<epsilon:
            prev_flag=False
    print("The number of iterations in newton raphson is:",num_iter)
    return Xr

def MainF(f,f_prime,start_point,end_point,epsil,num_Of_Parts):
    D_jump = float((end_point - start_point) / num_Of_Parts)
    for i in range(num_Of_Parts):
        temp1 = f(start_point)
        temp2 = f(start_point + D_jump)
        if temp1 * temp2 < 0:
            print("we found a range where there is an extreme point")
            root1 = newton_raphson(f, f_prime, start_point, start_point + D_jump, epsil)
            root2=secant_method(f,start_point,start_point + D_jump,epsil)
        start_point = start_point + D_jump
    print("The root in newton raphson is: ", root1)
    print("The root in secant_method is: ", root2)

f = x**3-x-1
f_prime = f.diff(x)
start_point=1
end_point=2
num_Of_Parts=9
epsil=0.000001
print("example:")
print("The start_point is:",start_point)
print("The end_point is:",end_point)
print("f: ",f)
print("f': ",f_prime)
MainF(lambdify(x, f),lambdify(x, f_prime),start_point,end_point,epsil,num_Of_Parts)