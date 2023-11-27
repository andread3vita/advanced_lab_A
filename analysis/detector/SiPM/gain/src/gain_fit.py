import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy.odr import *

Sipm=input("Which SiPM? [ex. A_SX]	")

filepath="txt_files/"+Sipm
df=pd.read_csv(filepath+".txt",sep="\t")

x_data=df["V[dV]"]
y_data=df["Gain[mV]"]

err_x=df["SigmaV[dV]"]
err_y=df["SigmaGain[mv]"]

def linear(p,x):
    m,q=p
    return m*x+q


# Create a model for fitting.
linear_model = Model(linear)

# Create a RealData object using our initiated data from above.
data = RealData(x_data, y_data, sx=err_x, sy=err_y)

# Set up ODR with the model and data.
odr = ODR(data, linear_model, beta0=[0., -50.])

# Run the regression.
out = odr.run()

# Use the in-built pprint method to give us results.
out.pprint()

fig,ax=plt.subplots(1,1,figsize=(13,12))
x_fit=np.linspace(280,x_data[len(x_data)-1],1000)
V_b=-out.beta[1]/out.beta[0]

errb0= out.sd_beta[0]
errb1=out.sd_beta[1]
sigmaVb=np.sqrt((out.beta[1]/out.beta[0]**2)**2*errb0**2+(1/out.beta[0])**2*errb1**2)

chi_squared = np.sum(out.delta**2)
ndf = len(x_data) - len(out.beta)


ax.errorbar(x=x_data,y=y_data,xerr=err_x,yerr=err_y,linestyle="none",capsize=8)
ax.axvline(V_b, color='red', linestyle='dashed', linewidth=1,label= "$V_b$ [dV]               	{}$\pm${}".format(round(V_b,2),round(sigmaVb,2)))
ax.plot(x_fit,linear(out.beta,x_fit),linestyle="dashed",label=f'y={round(out.beta[0],2)}$\cdot$x - {np.abs(round(out.beta[1],2))}')
#ax.set_title("Gain vs $V_{bias}$",fontsize=22)
ax.set_xlabel("$V_{bias}$ [dV]",fontsize=16)
ax.set_ylabel("Gain [mV]",fontsize=16)
ax.grid()
ax.legend(loc="upper left", fontsize=14,title_fontsize=16, title=f"$\chi^2$ / ndf				{chi_squared:.2f} / {ndf}")

plt.show()
image_name=Sipm+".pdf"
fig.savefig(image_name)
