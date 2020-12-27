import math
import numpy
import matplotlib.pyplot as plt
L = 0.3;
R1 = 0.5;
R2 = 0.5;
c = 299792458;
R = 0.99995;
F = (math.pi*math.sqrt(R))/(1 - R);
FSR = c/(2*L);
linewidth=FSR/F;
g1 = 1-L/R1;
g2 = 1-L/R2;

def freq(q,m):
    return q+(1+m)/(math.pi)*numpy.arccos(math.sqrt(g1*g2));

datatext1="Cavity Length ="+str(L)+"m, R1="+str(R1)+"m, R2="+str(R2)+"m, R="+str(R)
datatext2="FSR ="+str("{0:.2f}".format(FSR/1000000))+"MHz, Finesse="+str("{0:.2f}".format(F))+", linewidth="+str("{0:.2f}".format(linewidth/1000))+"kHz"
for q in range(-20,20,1):
    for m in range(0,30):
        f=freq(q,m)-freq(0,0);
        lb="("+str(q)+","+str(m)+")"

        if (f>-0.5)&(f<1.2):
            plt.plot(f, 1-0.02*m, marker='o', markersize=5, color=str(m*0.03))
            plt.plot([f, f], [ 1-0.02*m, 0], color=str(m*0.03), linestyle='-', linewidth=2)
            plt.text(f, 1-0.02*m,lb)

plt.text(0,1.3,"Resonant frequency/FSR for mode(q,m+n)")
plt.text(0,1.2,datatext1)
plt.text(0,1.1,datatext2)


plt.axis([-0.5, 1.3, 0, 1.4])
plt.show()
