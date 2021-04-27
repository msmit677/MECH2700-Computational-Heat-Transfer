"""
Author/s: Matthew Smith and Adam Wormington
Date: 15/10/19
Title: MECH2700 Assignment 2
"""

from numpy import zeros, hstack, cosh, sinh, sqrt, linspace
import pylab

#Define all global parameters that have been given
h = 2000.0 # W/m^2.K
k = 250.0 # W/m.K
Tw = 1710.0 # K
Ta = 23.0 # K
L = 10*10**-3 # m
t = 1*10**-3 # m
d = 300*10**-3 # m
P = 2*d + 2*t # m
Ac = t*d # m^2
qcomb = 14000 # W

#Task 4

"""
Create a function that returns an (n+1)x1 matrix 'b' in the form of:
    |           Tw            |
    |-(dx)^2(h)(P)(Ta)/(k)(Ac)|
    |-(dx)^2(h)(P)(Ta)/(k)(Ac)|
    |-(dx)^2(h)(P)(Ta)/(k)(Ac)|
    |-(dx)^2(h)(P)(Ta)/(k)(Ac)|
    |          h*Ta           |
"""

def b(n,dx):  
    b = zeros((n+1,1),float)   #Begin with a matrix full of zeros
    b[0,0] = Tw     #Add Tw to the first row
    b[n,0] = h*Ta   #Add h*Ta to the last row
    for j in range(1,n):
        b[j,0] = -dx**2*h*P*Ta/(k*Ac)   #Add the middle row entries
    return b

"""
Create a function which returns an (n+1)x(n+1) matrix 'A' of the form:
    |1             0                      0                    0   |
    |1   (2+dx^2(h)(P)/(k)(Ac))           1                    0   |
    |0             1            (2+dx^2(h)(P)/(k)(Ac))         1   |
    |0             0                    -k/dx                h+k/dx|
"""

def A(n,dx):
    A = zeros((n+1,n+1),float)  #Begin with a matrix full of zeros
    A[0,0] = 1      #Add 1 to the first entry in the first row
    A[n,n] = h+k/dx #Add the value for the last entry in the last row
    A[n,n-1] = -k/dx #Add the value for the second last entry in last row
    for j in range(1,n):
        A[j,j] = -(2+(dx**2*h*P/(k*Ac)))  #Add diagonal entry for middle rows
        A[j,j-1] = 1    #Add 1 to the left of the diagonal entry for each row
        A[j,j+1] = 1    #Add 1 to the right of the diagonal entry for each row
    return A

"""
Use a more efficient Gauss-Jordan algorithm for solving the linear system of
equations:
                               AT = b
"""

def solve(A,b):
    nrows, ncols = A.shape
    c = hstack([A,b])   #Create a matrix which adds b to the last column of A
    for j in range(1,nrows):
        c[j,:] = c[j,:] - c[j,j-1]*c[j-1,:]/c[j-1,j-1] #Make the 1s to the left of the diagonal 0
    for j in range(0,nrows):
        c[j,:] = c[j,:]/c[j,j] #Make the leading diagonal all 1s
    j = nrows - 2
    while j >= 0:
        c[j,:] = c[j,:] - c[j,j+1]*c[j+1,:]/c[j+1,j+1] #Make the values to the right of the diagonal 1s
        j += -1
    T = c[:,-1] #Return the final column of the matrix
    return T

#Task 5

"""
Create the analytic function that computes the analytic temperature along the
fin, as found in the text: Incropera and De Witt, “Fundamentals of Heat and 
Mass Transfer”
"""

def T_Analytic(x):
    theta_b = Tw-Ta
    m = sqrt((h*P)/(k*Ac))
    theta = theta_b*((cosh(m*(L-x))+(h/(m*k))*sinh(m*(L-x)))/(cosh(m*L)+(h/(m*k))*sinh(m*L)))
    Tx = theta + Ta
    return Tx

"""
Plot the temperature distribution over the fin at n = 5,10,20 and 200 on the
same graph along with the analytic temperture distribution
"""

ns = [5,10,20,200]
for n in ns:
    dx = L/n
    T = solve(A(n,dx),b(n,dx))  #Run the solve function
    xs = linspace(0,L,n+1)  #Create a list of x values along the fin
    pylab.figure(1)
    pylab.plot(xs,T, '-.', label="Numerical Solution (n = {0:0.1f})".format(n))  #Plot the Numerical Solution
    pylab.xlim(0,L)
    pylab.ylim(min(T)-100,max(T)+100)
    pylab.xlabel("x(m)")
    pylab.ylabel("T(K)")
    pylab.title("Temperature Distribution Over a Fin")   
pylab.figure(1)
xs = linspace(0,L,ns[-1]+1)
Ts = [T_Analytic(x) for x in xs]    #Run the analytic function
pylab.plot(xs,Ts,'r--', label="Analytical Solution") #Plot the analytic function on the same figure
pylab.legend()
pylab.show()

#Task 6

"""
Plot the relative error of the numerical solution compared with the analytic
function of the fin temperature along with the fin temperature over a range of 
n values 
"""

ns = [i for i in range(1,201)]
error = []  #List to store the error values at each value of n
T_num = []  #List to store the numerical temperatures at each value of n
T_real = [] #List to store the analytic temperatures at each value of n
for n in ns:
    dx = L/n
    xs = linspace(0,L,n+1)
    Ts = [T_Analytic(x) for x in xs]    #Run the Analytic function
    T = solve(A(n,dx),b(n,dx))      #Run the solve function
    error.append((abs(T[-1]-Ts[-1])/Ts[-1])*100)    #Add to the error list
    T_num.append(T[-1])     #Add to the numerical list
    T_real.append(T_Analytic(L))    #Add to the analytic list
pylab.figure(1)
pylab.plot(ns,error,'b-')   #Plot n vs. error (%) on one figure 
pylab.title("Relative Error at the Fin Tip")
pylab.xlim(0,max(ns))
pylab.ylim(0,max(error)+1)
pylab.xlabel("n")
pylab.ylabel("error (%)")
pylab.figure(2)
pylab.plot(ns,T_num,'b-',label="Numerical Solutions") #Plot n vs. numerical temp on another figure 
pylab.plot(ns,T_real,'r--',label="Analytical Solution (T = {0:0.3f}K)".format(T_Analytic(L)))
pylab.title("Temperature at the Fin Tip")
pylab.xlim(0,max(ns))
pylab.ylim(T_Analytic(L)-20,max(T_num)+5)
pylab.xlabel("n")
pylab.ylabel("T(K)")
pylab.legend()
pylab.show()

"""
Create a function which computes the number of nodes needed to produce a fin
temperature error within 0.1% of the analytic solution
"""

def T_Convergence():
    n = 1
    tol = 0.1   #Set the target error as 0.1%
    err = 1
    while err > tol:    #While the error is greater than 0.1 the function runs
        dx = L/n
        xs = linspace(0,L,n+1)
        Ts = [T_Analytic(x) for x in xs]    #Run analytic solution
        T = solve(A(n,dx),b(n,dx))  #Run numerical solution
        err = (abs(T[-1]-Ts[-1])/Ts[-1])*100    #Compute the error as a %
        n += 1
    return n-1

print("The number of nodes needed to produce an error within 0.1% =",T_Convergence())

#Task 7

"""
Using the finite-difference method, a function for the numerical solution 
for the heat transfer rate is created along with the analytical solution 
presented in the text: Incropera and De Witt, “Fundamentals of Heat and Mass 
Transfer”
"""

def q(n):
    dx = L/n
    T = solve(A(n,dx),b(n,dx))  #Run the solve function
    return k*Ac*(Tw-T[1])/dx    #q depends on the second T value

def q_Analytic(L,Ac,P):
    theta_b = Tw-Ta
    M = sqrt(h*P*k*Ac)*theta_b
    m = sqrt((h*P)/(k*Ac))
    q = M*((sinh(m*L)+(h/(m*k))*cosh(m*L))/(cosh(m*L)+(h/(m*k))*sinh(m*L)))
    return q

"""
Plot the numerical heat transfer rate and the analytical heat transfer rate 
on the same figure over a range of n values
"""

q_real = [] #Create an empty list to store the analytic q's
q_num = [] #Create an empty list to store the numerical q's
for n in ns:
    q_num.append(q(n))  #Add to the numerical list
    q_real.append(q_Analytic(L,Ac,P))   #Add to the analytic list
pylab.plot(ns,q_num,'b-', label="Numerical Solution") #Plot numerical solution
pylab.plot(ns,q_real,'r--', label="Analytical Solution (q = {0:0.3f}W)".format(q_Analytic(L,Ac,P)))
pylab.xlim(0,max(ns))
pylab.xlabel("n")
pylab.ylabel("q(W)")
pylab.title("Heat Transfer Rate from Fin")
pylab.legend()
pylab.show()

"""
Create a function which computes the number of nodes needed to produce a heat 
transfer rate error within 0.1% of the analytic solution
"""

def q_Convergence():
    n = 1
    error = (abs(q(n)-q_Analytic(L,Ac,P))/q_Analytic(L,Ac,P))*100 #Percentage error
    tol = 0.1 #Set the target error as 0.1%
    while error > tol:  #Function runs until error is lower than 0.1
        n += 1
        error = (abs(q(n)-q_Analytic(L,Ac,P))/q_Analytic(L,Ac,P))*100 #Percentage error
    return n

print("The number of nodes needed to produce an error within 0.1% =",q_Convergence())

print("Combustion Heat Transfer Rate =",q_Analytic(L,Ac,P))
if qcomb - q_Analytic(L,Ac,P) > 0:
    print("Therefore the fin doesn't provide sufficient cooling")
    print("We must change the geometry of the fin...")
else:
    print("Therefore the fin does provide sufficient cooling")

"""
Create functions which change different geometry parameters (L,d and t) so that 
qcomb is greater than 14kW. Compare the changes for each one to determine the
best geometry change
"""

L0 = L 
def Lchange():
    L = 10*10**-3
    q_Analytic(L,Ac,P)
    while qcomb > q_Analytic(L,Ac,P): #Runs while the q value is lower than qcomb
        L += 1*10**-5 #Increase by a reasonable number
        q_Analytic(L,d,t)
    change = (abs(q_Analytic(L,Ac,P)-q_Analytic(L0,Ac,P))/(q_Analytic(L0,Ac,P)))*100
    return L,change  #Return the required L value and the percentage change

Ac_0 = t*d
P_0 = 2*d + 2*t
def dchange():
    d = 300*10**-3
    Ac = t*d
    P = 2*t + 2*d
    q_Analytic(L,Ac,P)
    while qcomb > q_Analytic(L,Ac,P): #Runs while the q value is lower than qcomb
        d += 1*10**-5 #Increase by a reasonable number
        P = 2*d + 2*t
        Ac = t*d
        q_Analytic(L,Ac,P)
        change = (abs(q_Analytic(L,Ac,P)-q_Analytic(L,Ac_0,P_0))/(q_Analytic(L,Ac_0,P_0)))*100
    return d,change #Return the required d value and the percentage change

def tchange():
    t = 1*10**-3
    Ac = t*d
    P = 2*t + 2*d
    q_Analytic(L,Ac,P)
    while qcomb > q_Analytic(L,Ac,P): #Runs while the q value is lower than qcomb
        t += 1*10**-5 #Increase by a reasonable number
        P = 2*d + 2*t
        Ac = t*d
        q_Analytic(L,Ac,P)
        change = (abs(q_Analytic(L,Ac,P)-q_Analytic(L,Ac_0,P_0))/(q_Analytic(L,Ac_0,P_0)))*100
    return t,change #Return the required t value and the percentage change

print("When increasing L:")
print(Lchange())
print("When increasing d:")
print(dchange())
print("When increasing t:")
print(tchange())

  
        




