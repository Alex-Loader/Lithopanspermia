# -*- coding: utf-8 -*-
"""
Created on Wed May  1 10:57:01 2019

@author: Alex
"""

PPhyla=np.array([0,55,75,85,100])
PClasses=np.array([0,79.3,94.2,97.6,100])
POrders=np.array([0,35.2,61.7,78.8,100])

MPhyla=np.array([0,66.7,85.2,96.3,100])
MClasses=np.array([0,46.3,69.5,91.6,100])
MOrders=np.array([0,25.8,55.7,83.6,100])

EPhyla=np.array([0,41.6,71.4,90.5,100])
EClasses=np.array([0,35.6,86.7,95.6,100])
EOrders=np.array([0,23.2,64.7,81.5,100])


x_data1=np.array([0,20,45,70,95])
x_data2=np.array([0,25,50,75,100])
x_data3=np.array([0,30,55,80,105])

plt.figure()
plt.subplot(3,1,1)
plt.bar(x_data1, PPhyla,5, alpha=1, label="Phyla")
plt.bar(x_data2, PClasses,5, alpha=1, label="Classes")
plt.bar(x_data3, POrders,5, alpha=1, label="Orders")
plt.xlim(0,112.5)


plt.title("Cumulative rate of arrival - Proportional Taxonomic Levels")


plt.subplot(3,1,2)
plt.bar(x_data1, MPhyla,5, alpha=1, label="Phyla")
plt.bar(x_data2, MClasses,5, alpha=1, label="Classes")
plt.bar(x_data3, MOrders,5, alpha=1, label="Orders")
plt.xlim(0,112.5)
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

plt.ylabel("Percentage of Taxa selected")
plt.title("Modern ")


plt.subplot(3,1,3)
plt.bar(x_data1, EPhyla,5, alpha=1, label="Phyla")
plt.bar(x_data2, EClasses,5, alpha=1, label="Classes")
plt.bar(x_data3, EOrders,5, alpha=1, label="Orders")
plt.xlim(0,112.5)

plt.xlabel("Percentage of run completed")

plt.title("Erwin")


plt.figure()
plt.subplot(3,1,1)
plt.plot(x_data1, PPhyla, label="Proportional")
plt.plot(x_data1, MPhyla, label="Modern")
plt.plot(x_data1, EPhyla, label="Erwin")
plt.title("Cumulative Arrival Rates of Phyla, Classes and Orders")
plt.ylabel("Phyla")
plt.legend(loc='lower right')

plt.subplot(3,1,2)
plt.plot(x_data1, PClasses, label="Proportional")
plt.plot(x_data1, MClasses, label="Modern")
plt.plot(x_data1, EClasses, label="Erwin")
plt.ylabel("Classes")



plt.subplot(3,1,3)
plt.plot(x_data1, POrders, label="Proportional")
plt.plot(x_data1, MOrders, label="Modern")
plt.plot(x_data1, EOrders, label="Erwin")
plt.ylabel("Orders")
plt.xlabel("Percentage of run completed")