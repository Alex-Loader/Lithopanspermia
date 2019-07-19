# -*- coding: utf-8 -*-
"""
Created on Tue Feb  5 18:17:53 2019

@author: c1647707
"""


import numpy as np 
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


#Correct numbers - each corresponding to the number of Orders in each Class and Phylum, allowing for
#us to "tick off" each Class/Phylum as Orders are selected
PhylaOrder=(4,8,11,16,18,24,29,33,36,40,44,49,53,56,60,66,73,77,82,86,90)
ClassOrder=(4,8,11,16,18,24,29,33,36,40,44,49,53,56,60,66,73,77,82,86,90,
            94,96,103,105,110,115,118,121,127,133,138,143,149,152,158,161,
            169,173,179,183,186,190,194,200,204,208,211,216,218,224,229,
            233,236,240,244,249,256,260,266,273,277,286,290,296,303,305,
            310,315,321,325,329,333,338,343,349,352,358,361,367,369,375,
            379,383,386,390,394,400,406,410)


plt.figure()
Phyla_number=21
Class_number=90
Order_number=410

for i in range(1):
        
    #A is a masked array used to simulate the missing taxa and their arrival through uncovering masked orders as time progresses
    A4=np.full(Order_number, False)
    A4[:(Order_number-410)]=True
    np.random.shuffle(A4)
    A3=np.full(Order_number, False)
    A3[:(Order_number-410)]=True
    np.random.shuffle(A3)
    A2=np.full(Order_number, False)
    A2[:(Order_number-410)]=True
    np.random.shuffle(A2)
    A1=np.full(Order_number, False)
    A1[:(Order_number-205)]=True
    np.random.shuffle(A1)
    
    
    A3=A3+A4
    A2=A2+A3+A4
    A1=A1+A2+A3+A4
    
    base=np.arange(Order_number)
    repeats=25
    length=np.random.randint(0,11,repeats)
    maxlength=10
    b=np.zeros(maxlength)
    R=np.arange(100)
    
    Order1=np.ma.masked_array(base,A1)
    Order2=np.ma.masked_array(base,A2)
    Order3=np.ma.masked_array(base,A3)
    Order4=np.ma.masked_array(base,A4)
    
    
    #Important values for how many times the loops will repeat,and how the random selection process will
    #work. Using a random number of random Orders each step, with between 0 and 30 selected each time.
    #length shows the amount selected for each step - will change for each loop
    
    
    #The Xnew arrays will be where we tick off each taxa as it is selected, by turning the zeros to ones
    #as the loop runs. Choice is where we will put the randomly selected numbers for the Order selection
    Onew=np.tile(np.zeros(Order_number),(repeats,1))
    Cnew=np.tile(np.zeros(Class_number),(repeats,1))
    Pnew=np.tile(np.zeros(Phyla_number),(repeats,1))
    choice=np.zeros((repeats,maxlength))
    
    
    Order1=np.ma.masked_array(base,A1)
    #Section A of the loop - random number selection, where a randomly selects a number within the range of
    #the Order_number limit a random number of times, and repeats this for how many repeats are set. By then
    #adding these to the zero matrix b we can construct the newly selected orders
    for i in np.arange(repeats):
        a=(np.random.choice(Order1[~A1],length[i]))
        a.resize(b.shape)
        a=(b+a)
        choice[i]=a
        a=np.array([a], dtype=int)
        Onew[i,a]=1
    
    #For each newly selected group we need to "tick off" and prevent repeats - done with this check where if
    #any values are already a 1, the rest of the column will be set to all 0s and thus the Order has been
    #"ticked off"
    for i in np.arange(Order_number):
        for j in np.arange(repeats):
            if Onew[j,i]==1:
                Onew[j+1:,i]=0
    
    #Section B - To find the Classes these Orders correspond to, we use the ClassOrder array and match the ticked off 
    #Orders to each Class based off the numbers of Orders in each Class. The initial Class has to be done 
    #seperately due to the nature of the calculation
    for i in np.arange(repeats):
        for j in np.arange(Class_number):
            if j==0:
                if np.sum(Onew[i][ClassOrder[0]:ClassOrder[j]])==0:
                    Cnew[i][j]=0
                else:
                    Cnew[i][j]=1
            else:            
                if np.sum(Onew[i][ClassOrder[j-1]:ClassOrder[j]])==0:
                    Cnew[i][j]=0
                else:
                    Cnew[i][j]=1
    
    #Have to check for repeats again
    for i in np.arange(Class_number):
        for j in np.arange(repeats):
            if Cnew[j,i]==1:
                Cnew[j+1:,i]=0
    
    #Section C - The same process is then repeated for the Phyla, using PhylaOrder. 
    for i in np.arange(repeats):
        for j in np.arange(Phyla_number):
            if j==0:
                if np.sum(Onew[i][PhylaOrder[0]:PhylaOrder[j]])==0:
                    Pnew[i][j]=0
                else:
                    Pnew[i][j]=1
            else:            
                if np.sum(Onew[i][PhylaOrder[j-1]:PhylaOrder[j]])==0:
                    Pnew[i][j]=0
                else:
                    Pnew[i][j]=1
    
    #Checking for repeats
    for i in np.arange(Phyla_number):
        for j in np.arange(repeats):
            if Pnew[j,i]==1:
                Pnew[j+1:,i]=0
    
    #finally Section D - constructing the Xplot arrays. These are the summed values of
    #each row of the Xnew arrays, and are going to be used to plot the final results.
    Oplot=np.zeros(repeats)
    Cplot=np.zeros(repeats)
    Pplot=np.zeros(repeats)
    
    for i in np.arange(repeats):
        Oplot[i]=np.sum(Onew[i,:],dtype=int)
        Cplot[i]=np.sum(Cnew[i,:],dtype=int)
        Pplot[i]=np.sum(Pnew[i,:],dtype=int)
    
    #Finding the percentage of each taxa appearing in this period
    Opicked1=np.sum(Oplot)
    Opercentage1=Opicked1/Order_number
    Cpicked1=np.sum(Cplot)
    Cpercentage1=Cpicked1/Class_number
    Ppicked1=np.sum(Pplot)
    Ppercentage1=Ppicked1/Phyla_number
    
    print("""The percentage of taxa appearing in the first 25 timesteps are:
        Phyla %.5f, Classes %.5f, Orders %.5f""" %(Ppercentage1, Cpercentage1, Opercentage1))
    "_____________________________________________________________________________"
    #Unmasking more of the Order array
    
    Onew2=np.tile(np.zeros(Order_number),(repeats,1))
    Cnew2=np.tile(np.zeros(Class_number),(repeats,1))
    Pnew2=np.tile(np.zeros(Phyla_number),(repeats,1))
    length=np.random.randint(0,11,repeats)
    
    for i in np.arange(repeats):
        a2=(np.random.choice(Order2[~A2],length[i]))
        a2.resize(b.shape)
        a2=(b+a2)
        choice[i]=a2
        a2=np.array([a2], dtype=int)
        Onew2[i,a2]=1
    
    for i in np.arange(Order_number):
        for j in np.arange(repeats):
            if Onew2[j,i]==1:
                Onew2[j+1:,i]=0
            if np.sum(Onew[j:,i])!=0:
                Onew2[j:,i]=0
    
    for i in np.arange(repeats):
        for j in np.arange(Class_number):
            if j==0:
                if np.sum(Onew2[i][ClassOrder[0]:ClassOrder[j]])==0:
                    Cnew2[i][j]=0
                else:
                    Cnew2[i][j]=1
            else:            
                if np.sum(Onew2[i][ClassOrder[j-1]:ClassOrder[j]])==0:
                    Cnew2[i][j]=0
                else:
                    Cnew2[i][j]=1
    
    for i in np.arange(Class_number):
        for j in np.arange(repeats):
            if Cnew2[j,i]==1:
                Cnew2[j+1:,i]=0
            if np.sum(Cnew[j:,i])!=0:
                Cnew2[j:,i]=0
    
    for i in np.arange(repeats):
        for j in np.arange(Phyla_number):
            if j==0:
                if np.sum(Onew2[i][PhylaOrder[0]:PhylaOrder[j]])==0:
                    Pnew2[i][j]=0
                else:
                    Pnew2[i][j]=1
            else:            
                if np.sum(Onew2[i][PhylaOrder[j-1]:PhylaOrder[j]])==0:
                    Pnew2[i][j]=0
                else:
                    Pnew2[i][j]=1
    
    for i in np.arange(Phyla_number):
        for j in np.arange(repeats):
            if Pnew2[j,i]==1:
                Pnew2[j+1:,i]=0
            if np.sum(Pnew[j:,i])!=0:
                Pnew2[j:,i]=0
    
    
    Oplot2=np.zeros(repeats)
    Cplot2=np.zeros(repeats)
    Pplot2=np.zeros(repeats)
    
    for i in np.arange(repeats):
        Oplot2[i]=np.sum(Onew2[i,:],dtype=int)
        Cplot2[i]=np.sum(Cnew2[i,:],dtype=int)
        Pplot2[i]=np.sum(Pnew2[i,:],dtype=int)
    
    #Finding the percentage of each taxa appearing in this period
    Opicked2=np.sum(Oplot2)
    Opercentage2=Opicked2/Order_number
    Cpicked2=np.sum(Cplot2)
    Cpercentage2=Cpicked2/Class_number
    Ppicked2=np.sum(Pplot2)
    Ppercentage2=Ppicked2/Phyla_number
    
    print("""The percentage of taxa appearing in the second 25 timesteps are:
        Phyla %.5f, Classes %.5f, Orders %.5f""" %(Ppercentage2, Cpercentage2, Opercentage2))
    print("""The total percentage of taxa appearing by the 50th timestep are:
        Phyla %.5f, Classes %.5f, Orders %.5f""" %(Ppercentage2+Ppercentage1, 
        Cpercentage2+Cpercentage1, Opercentage2+Opercentage1))
    "_____________________________________________________________________________"
    
    Onew3=np.tile(np.zeros(Order_number),(repeats,1))
    Cnew3=np.tile(np.zeros(Class_number),(repeats,1))
    Pnew3=np.tile(np.zeros(Phyla_number),(repeats,1))
    length=np.random.randint(0,11,repeats)
    
    for i in np.arange(repeats):
        a3=(np.random.choice(Order3[~A3],length[i]))
        a3.resize(b.shape)
        a3=(b+a3)
        choice[i]=a3
        a3=np.array([a3], dtype=int)
        Onew3[i,a3]=1
    
    for i in np.arange(Order_number):
        for j in np.arange(repeats):
            if Onew3[j,i]==1:
                Onew3[j+1:,i]=0
            if np.sum(Onew[j:,i])!=0:
                Onew3[j:,i]=0
            if np.sum(Onew2[j:,i])!=0:
                Onew3[j:,i]=0
    
    for i in np.arange(repeats):
        for j in np.arange(Class_number):
            if j==0:
                if np.sum(Onew3[i][ClassOrder[0]:ClassOrder[j]])==0:
                    Cnew3[i][j]=0
                else:
                    Cnew3[i][j]=1
            else:            
                if np.sum(Onew3[i][ClassOrder[j-1]:ClassOrder[j]])==0:
                    Cnew3[i][j]=0
                else:
                    Cnew3[i][j]=1
    
    for i in np.arange(Class_number):
        for j in np.arange(repeats):
            if Cnew3[j,i]==1:
                Cnew3[j+1:,i]=0
            if np.sum(Cnew[j:,i])!=0:
                Cnew3[j:,i]=0
            if np.sum(Cnew2[j:,i])!=0:
                Cnew3[j:,i]=0
    
    for i in np.arange(repeats):
        for j in np.arange(Phyla_number):
            if j==0:
                if np.sum(Onew3[i][PhylaOrder[0]:PhylaOrder[j]])==0:
                    Pnew3[i][j]=0
                else:
                    Pnew3[i][j]=1
            else:            
                if np.sum(Onew3[i][PhylaOrder[j-1]:PhylaOrder[j]])==0:
                    Pnew3[i][j]=0
                else:
                    Pnew3[i][j]=1
    
    for i in np.arange(Phyla_number):
        for j in np.arange(repeats):
            if Pnew3[j,i]==1:
                Pnew3[j+1:,i]=0
            if np.sum(Pnew[j:,i])!=0:
                Pnew3[j:,i]=0
            if np.sum(Pnew2[j:,i])!=0:
                Pnew3[j:,i]=0
                
    Oplot3=np.zeros(repeats)
    Cplot3=np.zeros(repeats)
    Pplot3=np.zeros(repeats)
    
    for i in np.arange(repeats):
        Oplot3[i]=np.sum(Onew3[i,:],dtype=int)
        Cplot3[i]=np.sum(Cnew3[i,:],dtype=int)
        Pplot3[i]=np.sum(Pnew3[i,:],dtype=int)
    
    #Finding the percentage of each taxa appearing in this period
    Opicked3=np.sum(Oplot3)
    Opercentage3=Opicked3/Order_number
    Cpicked3=np.sum(Cplot3)
    Cpercentage3=Cpicked3/Class_number
    Ppicked3=np.sum(Pplot3)
    Ppercentage3=Ppicked3/Phyla_number
    
    print("""The percentage of taxa appearing in the third 25 timesteps are:
        Phyla %.5f, Classes %.5f, Orders %.5f""" %(Ppercentage3, Cpercentage3, Opercentage3))
    print("""The total percentage of taxa appearing by the 75th timestep are:
        Phyla %.5f, Classes %.5f, Orders %.5f""" %(Ppercentage3+Ppercentage2+Ppercentage1, 
        Cpercentage3+Cpercentage2+Cpercentage1, Opercentage3+Opercentage2+Opercentage1))
    "_____________________________________________________________________________"
    
    Onew4=np.tile(np.zeros(Order_number),(repeats,1))
    Cnew4=np.tile(np.zeros(Class_number),(repeats,1))
    Pnew4=np.tile(np.zeros(Phyla_number),(repeats,1))
    length=np.random.randint(0,11,repeats)
    
    for i in np.arange(repeats):
        a4=(np.random.choice(Order4[~A4],length[i]))
        a4.resize(b.shape)
        a4=(b+a4)
        choice[i]=a4
        a4=np.array([a4], dtype=int)
        Onew4[i,a4]=1
    
    for i in np.arange(Order_number):
        for j in np.arange(repeats):
            if Onew4[j,i]==1:
                Onew4[j+1:,i]=0
            if np.sum(Onew[j:,i])!=0:
                Onew4[j:,i]=0
            if np.sum(Onew2[j:,i])!=0:
                Onew4[j:,i]=0
            if np.sum(Onew3[j:,i])!=0:
                Onew4[j:,i]=0
    
    for i in np.arange(repeats):
        for j in np.arange(Class_number):
            if j==0:
                if np.sum(Onew4[i][ClassOrder[0]:ClassOrder[j]])==0:
                    Cnew4[i][j]=0
                else:
                    Cnew4[i][j]=1
            else:            
                if np.sum(Onew4[i][ClassOrder[j-1]:ClassOrder[j]])==0:
                    Cnew4[i][j]=0
                else:
                    Cnew4[i][j]=1
    
    for i in np.arange(Class_number):
        for j in np.arange(repeats):
            if Cnew4[j,i]==1:
                Cnew4[j+1:,i]=0
            if np.sum(Cnew[j:,i])!=0:
                Cnew4[j:,i]=0
            if np.sum(Cnew2[j:,i])!=0:
                Cnew4[j:,i]=0
            if np.sum(Cnew3[j:,i])!=0:
                Cnew4[j:,i]=0
    
    for i in np.arange(repeats):
        for j in np.arange(Phyla_number):
            if j==0:
                if np.sum(Onew4[i][PhylaOrder[0]:PhylaOrder[j]])==0:
                    Pnew4[i][j]=0
                else:
                    Pnew4[i][j]=1
            else:            
                if np.sum(Onew4[i][PhylaOrder[j-1]:PhylaOrder[j]])==0:
                    Pnew4[i][j]=0
                else:
                    Pnew4[i][j]=1
    
    for i in np.arange(Phyla_number):
        for j in np.arange(repeats):
            if Pnew4[j,i]==1:
                Pnew4[j+1:,i]=0
            if np.sum(Pnew[j:,i])!=0:
                Pnew4[j:,i]=0
            if np.sum(Pnew2[j:,i])!=0:
                Pnew4[j:,i]=0
            if np.sum(Pnew3[j:,i])!=0:
                Pnew4[j:,i]=0
    
                
    Oplot4=np.zeros(repeats)
    Cplot4=np.zeros(repeats)
    Pplot4=np.zeros(repeats)
    
    for i in np.arange(repeats):
        Oplot4[i]=np.sum(Onew4[i,:],dtype=int)
        Cplot4[i]=np.sum(Cnew4[i,:],dtype=int)
        Pplot4[i]=np.sum(Pnew4[i,:],dtype=int)
    
    #Finding the percentage of each taxa appearing in this period
    Opicked4=np.sum(Oplot4)
    Opercentage4=Opicked4/Order_number
    Cpicked4=np.sum(Cplot4)
    Cpercentage4=Cpicked4/Class_number
    Ppicked4=np.sum(Pplot4)
    Ppercentage4=Ppicked4/Phyla_number
    
    print("""The percentage of taxa appearing in the final 25 timesteps are:
        Phyla %.5f, Classes %.5f, Orders %.5f""" %(Ppercentage4, Cpercentage4, Opercentage4))
    print("""The total percentage of taxa appearing by the final timestep are:
        Phyla %.5f, Classes %.5f, Orders %.5f""" %(Ppercentage4+Ppercentage3+Ppercentage2+Ppercentage1, 
        Cpercentage4+Cpercentage3+Cpercentage2+Cpercentage1, Opercentage4+Opercentage3+Opercentage2+Opercentage1))
    
    #The Xplot arrays now need to be combined to create the final results, which
    #can be plotted together along the full time period of the combined loops
    Oplot=np.append(Oplot,Oplot2)
    Oplot=np.append(Oplot,Oplot3)
    Oplot=np.append(Oplot,Oplot4)
    
    Cplot=np.append(Cplot,Cplot2)
    Cplot=np.append(Cplot,Cplot3)
    Cplot=np.append(Cplot,Cplot4)
    
    Pplot=np.append(Pplot,Pplot2)
    Pplot=np.append(Pplot,Pplot3)
    Pplot=np.append(Pplot,Pplot4)
    
    
    #Xfinal arrays are combining all the Xnew arrays - the times each array was selected - into one single array.
    #Used to find the total numbers picked of each taxa as well as looking at which Phyla, Classes and Orders
    #were never selected by the code
    Ofinal=Onew+Onew2+Onew3+Onew4
    Cfinal=Cnew+Cnew2+Cnew3+Cnew4
    Pfinal=Pnew+Pnew2+Pnew3+Pnew4
    
    #Checking which Phyla, Classes and Orders were picked and which ones were missed - will be more useful when
    #changing taxa numbers and extinctions are included.
    Opicked=np.sum(Ofinal, axis=0)
    Cpicked=np.sum(Cfinal, axis=0)
    Ppicked=np.sum(Pfinal, axis=0)
    
    #Plotting the results - using R, a sum of all the times the loop has repeated, and each of the Xplot arrays
    
    
    plt.subplot(3,1,1)
    plt.title("Proportional Taxa Distribution compared to Erwin Data",fontsize='18')
    plt.bar(R,Pplot,color=(0.2, 0.4, 0.6, 1))
    plt.ylabel("Phyla",fontsize="16")
    plt.subplot(3,1,2)
    plt.bar(R,Cplot,color=(0.2, 0.4, 0.6, 1))
    plt.ylabel("Classes",fontsize="16")
    plt.subplot(3,1,3)
    plt.bar(R,Oplot,color=(0.2, 0.4, 0.6, 1))
    plt.ylabel("Orders",fontsize="16")
    plt.xlabel("Time",fontsize='18')
    
    
    #Finally looking at the numbers of each taxa selected, and comparing it as a percentage of the total available
    Ofinal=np.sum(Ofinal, axis=1)
    Onum=np.sum(Ofinal)
    Opercent4=Onum/Order_number
    Cfinal=np.sum(Cfinal, axis=1)
    Cnum=np.sum(Cfinal)
    Cpercent4=Cnum/Class_number
    Pfinal=np.sum(Pfinal, axis=1)
    Pnum=np.sum(Pfinal)
    Ppercent4=Pnum/Phyla_number

 #PLotting a best fit line
x_phyla=np.arange(1,101)
y_phyla=Pplot
z_phyla=np.polyfit(x_phyla,y_phyla,2)
f_phyla=np.poly1d(z_phyla)
x_phyla_new=np.linspace(x_phyla[0],x_phyla[-1],100)
y_phyla_new=f_phyla(x_phyla_new)
plt.subplot(3,1,1)
def func(x, a, b, c):
        return a * np.exp(-b * x) + c
poptp, pcovp = curve_fit(func, x_phyla, y_phyla)
plt.plot(x_phyla, func(x_phyla, *poptp), color=(0.2,0.4,0.6,1), label="Varying Random Selection")


x_class=np.arange(1,101)
y_class=Cplot
z_class=np.polyfit(x_class,y_class,6)
f_class=np.poly1d(z_class)
x_class_new=np.linspace(x_class[0],x_class[-1],100)
y_class_new=f_class(x_class_new)
plt.subplot(3,1,2)
poptc, pcovc = curve_fit(func, x_class, y_class)
plt.plot(x_class, func(x_class, *poptc), color=(0.2,0.4,0.6,1))
   

x_order=np.arange(1,101)
y_order=Oplot
z_order=np.polyfit(x_order,y_order,9)
f_order=np.poly1d(z_order)
x_order_new=np.linspace(x_order[0],x_order[-1],100)
y_order_new=f_order(x_order_new)
plt.subplot(3,1,3)

popto, pcovo = curve_fit(func, x_order, y_order)
plt.plot(x_order, func(x_order, *popto), color=(0.2,0.4,0.6,1))


PhylaOrder=(4,8,11,16,18,24,29,33,36,40,44,49,53,56,60,66,73,77,82,86,90)
ClassOrder=(4,8,11,16,18,24,29,33,36,40,44,49,53,56,60,66,73,77,82,86,90,
            94,96,103,105,110,115,118,121,127,133,138,143,149,152,158,161,
            169,173,179,183,186,190,194,200,204,208,211,216,218,224,229,
            233,236,240,244,249,256,260,266,273,277,286,290,296,303,305,
            310,315,321,325,329,333,338,343,349,352,358,361,367,369,375,
            379,383,386,390,394,400,406,410)



Phyla_number=21
Class_number=90
Order_number=410

base=np.arange(Order_number)
repeats=100
length=5
R=np.arange(repeats)


Onew=np.tile(np.zeros(Order_number),(repeats,1))
Cnew=np.tile(np.zeros(Class_number),(repeats,1))
Pnew=np.tile(np.zeros(Phyla_number),(repeats,1))
choice=np.zeros((repeats,length))

for i in range(1):
    for i in np.arange(repeats):
        a=(np.random.choice(base,length))
        choice[i]=a
        Onew[i,a]=1
    
    for i in np.arange(Order_number):
        for j in np.arange(repeats):
            if Onew[j,i]==1:
                Onew[j+1:,i]=0
    
    for i in np.arange(repeats):
        val=0
        for j in np.arange(Class_number):
            if np.sum(Onew[i][val:val+4])==0:
                Cnew[i][j]=0
            else:
                Cnew[i][j]=1
            val=val+4
        
    for i in np.arange(Class_number):
        for j in np.arange(repeats):
            if Cnew[j,i]==1:
                Cnew[j+1:,i]=0
    
    for i in np.arange(repeats):
        val=0
        for j in np.arange(Phyla_number):
            if np.sum(Onew[i][val:val+16])==0:
                Pnew[i][j]=0
            else:
                Pnew[i][j]=1
            val=val+16
    
    for i in np.arange(Phyla_number):
        for j in np.arange(repeats):
            if Pnew[j,i]==1:
                Pnew[j+1:,i]=0
    
    Oplot=np.zeros(repeats)
    Cplot=np.zeros(repeats)
    Pplot=np.zeros(repeats)
    
    for i in np.arange(repeats):
        Oplot[i]=np.sum(Onew[i,:])
        Cplot[i]=np.sum(Cnew[i,:])
        Pplot[i]=np.sum(Pnew[i,:])
    
    

    """plt.subplot(3,1,1)
    plt.title("Arrival of Phyla, Classes and Order",fontsize='16')
    plt.bar(R,Pplot,color=(0.2, 0.4, 0.6, 0.4))
    plt.ylabel("Phyla",fontsize="16")
    plt.subplot(3,1,2)
    plt.bar(R,Cplot,color=(0.2, 0.4, 0.6, 0.4))
    plt.ylabel("Classes",fontsize="16")
    plt.subplot(3,1,3)
    plt.bar(R,Oplot,color=(0.2, 0.4, 0.6, 0.4))
    plt.ylabel("Orders",fontsize="16")
    plt.xlabel("Time",fontsize='18')
    plt.show"""
    
#PLotting a best fit line
x_phyla2=np.arange(1,101)
y_phyla2=Pplot
z_phyla2=np.polyfit(x_phyla2,y_phyla2,2)
f_phyla2=np.poly1d(z_phyla2)
x_phyla2_new=np.linspace(x_phyla2[0],x_phyla2[-1],100)
y_phyla2_new=f_phyla2(x_phyla2_new)
plt.subplot(3,1,1)
#plt.plot(x_phyla_new, y_phyla_new, color=(0,1,0,0.6))
def func2(x, a, b, c):
    return a * np.exp(-b * x) + c

popt2p, pcov2p = curve_fit(func, x_phyla2, y_phyla2)
plt.plot(x_phyla, func2(x_phyla2, *popt2p), "r--", label="Constant Selection")
plt.legend()

x_class2=np.arange(1,101)
y_class2=Cplot
z_class2=np.polyfit(x_class2,y_class2,6)
f_class2=np.poly1d(z_class2)
x_class2_new=np.linspace(x_class2[0],x_class2[-1],100)
y_class2_new=f_class2(x_class2_new)
plt.subplot(3,1,2)
#plt.plot(x_class_new, y_class_new, color=(0,1,0,0.6))
popt2c, pcov2c= curve_fit(func2, x_class2, y_class2)

plt.plot(x_class2, func2(x_class2, *popt2c), "r--")

    

x_order2=np.arange(1,101)
y_order2=Oplot
z_order2=np.polyfit(x_order2,y_order2,9)
f_order2=np.poly1d(z_order2)
x_order2_new=np.linspace(x_order2[0],x_order2[-1],100)
y_order2_new=f_order2(x_order2_new)
plt.subplot(3,1,3)
#plt.plot(x_order_new, y_order_new, color=(0,1,0,0.6))
popt2o, pcov2o = curve_fit(func2, x_order2, y_order2)
plt.plot(x_order2, func2(x_order, *popt2o), "r--")

plt.show



"""

#Finding the difference in the lines
NormalPhyla=func2(x_phyla2, *popt2p)
NormalClass=func2(x_class2, *popt2c)
NormalOrder=func2(x_order2, *popt2o)

VaryingPhyla=func(x_phyla, *poptp)
VaryingClass=func(x_class, *poptc)
VaryingOrder=func(x_order, *popto)

Pdiff=VaryingPhyla-NormalPhyla
Cdiff=VaryingClass-NormalClass
Odiff=VaryingOrder-NormalOrder

zeros=np.zeros_like(R)

plt.figure()
plt.subplot(3,1,1)
plt.title("Difference in taxa levels for the Proportional Distribution",fontsize='18')
plt.plot(R,Pdiff,color=(0.2, 0.4, 0.6, 1))
plt.plot(R,zeros, "k--")
plt.ylabel("Phyla",fontsize="16")
plt.subplot(3,1,2)
plt.plot(R,Cdiff,color=(0.2, 0.4, 0.6, 1))
plt.plot(R,zeros, "k--")
plt.ylabel("Classes",fontsize="16")
plt.subplot(3,1,3)
plt.plot(R,Odiff,color=(0.2, 0.4, 0.6, 1))
plt.plot(R,zeros, "k--")
plt.ylabel("Orders",fontsize="16")
plt.xlabel("Time",fontsize='18')
"""