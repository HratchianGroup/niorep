#!/usr/bin/python

#######################################################################################################################
#
#  NIO  v. 1.0
#  By Hassan Harb
#  Program that reads in SCF Densities of ground and detached states and returns the eigenvalues and eigenvectors of
#  natural ionization orbitals
#
#  usage: python nio.py ground.fchk detached.fchk
#
#  For more details, check: J. Chem. Phys., 144, 204117 (2016)
#
#  Last edited by Hassan Harb, September 13, 2018
#
#######################################################################################################################

from __future__ import division
import sys
import math
import numpy as np
from numpy import genfromtxt
import csv
from decimal import Decimal 


#Function: Symmetrize: Transform packed array into NBasis x NBasis 
def symmetrize(a):
  Nbas = int((np.sqrt(8*len(a)+1)-1)/2)
  b = np.zeros((Nbas,Nbas))
  n = 0
  for i in range(0,Nbas):
    for j in range(0,i+1):
      b[i,j]=a[n]
      b[j,i]=a[n]
      n=n+1
  return b

#Function: convert output to scientific notation
def sci_notation(n):
#    n = np.real(n)
    a = '%.8E' % n
    return '%.8E' % np.real(n)
#    return a.split('E')[0].rstrip('0').rstrip('.') + 'E' + a.split('E')[1]


#Part 1: Read in the matrix files from two checkpoint files

NBasis = 0
filename1 = sys.argv[1]
filename2 = sys.argv[2]
filename3 = "NIO-"+filename1
acc = 8 #accuracy to the acc's decimal place

print "NIO: Calculates the density difference between two checkpoint files.\n"
print "Checkpoint 1:", filename1
print "Checkpoint 2:", filename2
 
with open(filename1,'r') as origin:
    for line in origin:
        if  "Number of basis functions" in line:
          words = line.split()
          for i in words:
              for letter in i:
                  if(letter.isdigit()):
                      NBasis = NBasis*10 + int(letter)   

print "Number of Basis Functions = ", NBasis, "\n"
MOElements = NBasis * NBasis
print "The code will look for ", MOElements, " elements of the MO Coeffienct matrix\n"
PElements = int(NBasis*(NBasis+1)/2)
print "The code will look for ", PElements, " elements of the Density matrix\n"
NBasis2=0

with open(filename2,'r') as origin:
    for line in origin:
        if  "Number of basis functions" in line:
          words = line.split()
          for i in words:
              for letter in i:
                  if(letter.isdigit()):
                      NBasis2 = NBasis2*10 + int(letter)

if (NBasis != NBasis2):
    print "ERROR: Files have different numbers of basis functions"
    exit()

MOlines = int(MOElements/5) + 1
Plines = int(PElements/5) + 1

if (MOElements % 5 ==0):
   MOlines = int(MOElements/5)
if (PElements % 5 ==0):
   Plines = int(PElements/5)

print "MO Lines = ", MOlines, "\n"
print "P Lines = ", Plines, "\n"

MOraw = np.zeros(MOElements)
TotalPraw1 = np.zeros(PElements) 
SpinPraw1 = np.zeros(PElements)

TotalPraw2 = np.zeros(PElements)
SpinPraw2 = np.zeros(PElements)

p = 0
r = 0
AOE = 0
with open(filename1,'r') as origin:
    for i, line  in enumerate(origin):
        if "Alpha Orbital Energies" in line:
              AOE = i
        if  "Alpha MO coefficients" in line:              
              i=i+1
              AMO=i
              print "Alpha MO coefficients starts at line :", i
              j=i+MOlines-1
              print "Alpha MO coefficients ends at line :", j
              for m in range(0,j-i+1):
                 nextline = origin.next()
                 nextline = nextline.split()
                 for p in range(p,len(nextline)):
                   MOraw[r] = nextline[p]
                   r = r+1
                 p = 0
        if "Beta MO coefficients" in line:
              i=i+1
              BMO = i
        if  "Total SCF Density" in line:
              i=i+1
              r = 0
              p = 0
              print "Total SCF Density starts at line :", i
              j=i+Plines-1
              print "Total SCF Density ends at line :", j
              for m in range(0,j-i+1):
                 nextline = origin.next()
                 nextline = nextline.split()
                 for p in range(p,len(nextline)):
                   TotalPraw1[r] = nextline[p]
                   r = r+1
                 p = 0
        if  "Spin SCF Density" in line:
              i=i+1
              r = 0
              p = 0
              print "Spin SCF Density starts at line: ", i
              j=i+Plines-1
              print "Spin SCF Density ends at line: ", j
              for m in range(0,j-i+1):
                 nextline = origin.next()
                 nextline = nextline.split()
                 for p in range(p,len(nextline)):
                   SpinPraw1[r] = nextline[p]
                   r = r+1
                 p = 0

with open(filename2,'r') as origin:
    for i, line  in enumerate(origin):
        if  "Total SCF Density" in line:
              i=i+1
              r = 0
              p = 0
              print "Total SCF Density starts at line :", i
              j=i+Plines-1
              print "Total SCF Density ends at line :", j
              for m in range(0,j-i+1):
                 nextline = origin.next()
                 nextline = nextline.split()
                 for p in range(p,len(nextline)):
                   TotalPraw2[r] = nextline[p]
                   r = r+1
                 p = 0
        if  "Spin SCF Density" in line:
              i=i+1
              r = 0
              p = 0
              print "Spin SCF Density starts at line: ", i
              j=i+Plines-1
              print "Spin SCF Density ends at line: ", j
              for m in range(0,j-i+1):
                 nextline = origin.next()
                 nextline = nextline.split()
                 for p in range(p,len(nextline)):
                   SpinPraw2[r] = nextline[p]
                   r = r+1
                 p = 0

MOraw = np.around(MOraw,decimals=acc)
TotalPraw1 = np.around(TotalPraw1,decimals=acc)
SpinPraw1 = np.around(SpinPraw1,decimals=acc)
TotalPraw2 = np.around(TotalPraw2,decimals=acc)
SpinPraw2 = np.around(SpinPraw2,decimals=acc)

print "Stored Alpha MO Coefficients matrix = \n", MOraw
print "Stored Total SCF Density matrix = \n", TotalPraw1
print "Stored Spin SCF Density matrix = \n", SpinPraw1

PalphaRaw1 = (np.add(TotalPraw1,SpinPraw1)) * 0.5
PbetaRaw1 = (np.subtract(TotalPraw1,SpinPraw1)) * 0.5 

PalphaRaw2 = (np.add(TotalPraw2,SpinPraw2)) * 0.5
PbetaRaw2 = (np.subtract(TotalPraw2,SpinPraw2)) * 0.5

print "Ground state: Alpha Density matrix =\n", PalphaRaw1
print "Ground state: Beta Density matrix =\n", PbetaRaw1

print "Detached state: Alpha Density matrix =\n", PalphaRaw2
print "Detached state: Beta Density matrix =\n", PbetaRaw2

#Part 2: The NIO Code

P1a = np.zeros((NBasis,NBasis))
P2a = np.zeros((NBasis,NBasis))
P1b = np.zeros((NBasis,NBasis))
P2b = np.zeros((NBasis,NBasis))

#Density Matrices

P1a = symmetrize(PalphaRaw1)
P2a = symmetrize(PalphaRaw2)
P1b = symmetrize(PbetaRaw1)
P2b = symmetrize(PbetaRaw2)

print "Alpha Density 1 = \n", P1a
print "Beta Density 1 = \n", P1b
print "Alpha Density 2 = \n", P2a
print "Beta Density 2 = \n", P2b

#Difference Density Matrices
Dpa = np.subtract(P2a,P1a)
Dpb = np.subtract(P2b,P1b)

print "Alpha Difference Density = \n", Dpa
print "Beta Difference Density = \n", Dpb

#MO Coefficient Matrix
C = np.zeros((NBasis,NBasis))
t=0
for i in range(0,NBasis):
   for j in range(0,NBasis):
     C[j,i]=MOraw[t]
     t=t+1

print "MO Coefficient Matrix = \n", C

#Inverse of C:
CInv = np.linalg.inv(C)
#Overlap Matrix:
# S = (C**(-1)t*(C**(-1))
S = np.dot(np.transpose(CInv),CInv)
print "Overlap Matrix = \n", S

#Checkpoint: <Dp.S> = -1

DpSa = np.dot(Dpa,S)
DpSb = np.dot(Dpb,S)

print "Alpha <DeltaP S> = ", np.trace(DpSa)
print "Beta <DeltaP S> = ", np.trace(DpSb)

#Calculate S**(0.5)
Svals, Svecs = np.linalg.eig(S)
Sval_minhalf = (np.diag(Svals**(0.5)))
Shalf = np.dot(Svecs,np.dot(Sval_minhalf,np.transpose(Svecs)))
print "S**(0.5) = \n", Shalf

#Calculae S**(0.5).Dpa.S**(0.5) and S**(0.5).Dpb.S**(0.5)
#Let these two matrices be Z1 and Z2, respectively

Z1 = np.dot(Shalf,np.dot(Dpa,Shalf))
Z2 = np.dot(Shalf,np.dot(Dpb,Shalf))

print "Z1 = \n", Z1
print "Z2 = \n", Z2

#Diagonaize Z1 and Z2, obtain eigenvalues (e1 and e2) and eigenvectors (U1 and U2)

e1, U1 = np.linalg.eig(Z1)
e2, U2 = np.linalg.eig(Z2)

#e1 = np.sort(e1)
#e2 = np.sort(e2)

print "Alpha NIO Eigenvalues = \n", e1
print "Beta NIO Eigenvalues = \n", e2
print "Alpha NIO Eigenvectors = \n", U1
print "Beta NIO Eigenvectors = \n", U2

#Backtransform U1 and U2 into the AO Basis

V1 = np.dot(np.linalg.inv(Shalf),U1)
V2 = np.dot(np.linalg.inv(Shalf),U2)

print "V1 =\n", V1, "\n"
print "V2 =\n", V2, "\n" 

# Part 3: write everything to a new fchk file:
# We need to replace Alpha Orbital Energies with e1 & Beta orbital energies with e2
# We also need to replace Alpha MO Coefficients with V1 & Beta MO Coefficients with V2

pointer=0
counter=1
print "AMO = ", AMO
print "BMO = ", BMO
with open(filename1,'r') as origin:
  data = origin.readlines()
  print "Line 5 = ", data[5]
  print "Line 6 = ", data[6]
  print "length of data = ", len(data)
  if "Alpha Orbital Energies" in line:
    AOE = i  
    BOE = AOE + int(NBasis/5) + 1
    print "AOE line is: ", data[AOE]
    print "BOE line is: ", data[BOE]
  with open(filename3,'w') as f2:

      print "Writing results to new output file: ", filename3, " ... "

#Read in data from old fchk (unchanged)
      while (pointer < AOE+1):
         f2.write(data[pointer])
         pointer = pointer+1
#Write the Alpha eigenvalues, overwrite the ALpha Orbital energies
      for j in range(0,NBasis):
          f2.write(" ")
          if (e1[j] >= 0):
             f2.write(" ")
          f2.write(str(sci_notation(e1[j])))
          if (counter%5 == 0):
              f2.write("\n")
              counter=0
          counter=counter+1
      counter =1      
#Write the Beta eigenvalues, overwrite the Beta Orbital energies
      f2.write("\n")
      BOE = AOE + (int(NBasis/5)+2)
      f2.write(data[BOE])
      for j in range(0,NBasis):
          f2.write(" ")
          if (e2[j] >= 0):
             f2.write(" ")
          f2.write(str(sci_notation(e2[j])))
          if (counter%5 ==0):
              f2.write("\n")
              counter=0
          counter = counter+1
      counter =1
      f2.write("\n")
#Write the Alpha NIO eigenvectors, overwrite the Alpha MO Coefficients
      AMO = BOE + (int(NBasis/5)+2)
      f2.write(data[AMO])
      for i in range(0,NBasis):
          for j in range(0,NBasis):
               f2.write(" ")
               if (V1[j,i] >= 0):
                  f2.write(" ")
               f2.write(str(sci_notation(V1[j,i])))
               if (counter%5 ==0):
                   f2.write("\n")
                   counter=0
               counter = counter + 1
         # counter =1
      counter = 1
      f2.write("\n")
#Write the Beta NIO eigenvectors, overwrite the Beta MO coefficients
      BMO = AMO + (int(NBasis*NBasis/5))+2
      f2.write(data[BMO])
      for i in range(0,NBasis):
             for j in range(0,NBasis):
                  f2.write(" ")
                  if (V2[j,i] >= 0):
                     f2.write(" ")
                  f2.write(str(sci_notation(V2[j,i])))
                  if (counter%5 ==0):
                      f2.write("\n")
                      counter=0
                  counter = counter + 1
           #  counter = 1
      counter = 1
      f2.write("\n")
      print "NBasis element = ", V2[0,NBasis-1]
#Copy the rest of the checkpoint file (unchanged)      
      pointer = BMO + (int(NBasis*NBasis/5))+2
      while (pointer < len(data)):
         f2.write(data[pointer])
         pointer = pointer+1
print "Copying data to new checkpoint file: DONE"
print "NIO Calculation done!"
