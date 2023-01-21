# %%
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

# %%
class Tbdy_spektrum:
    
    Ss_range = [0.25 , 0.50 , 0.75, 1.00 , 1.25 , 1.50]
    Fs_table = {"ZA": [0.8 , 0.8 , 0.8 , 0.8 , 0.8 , 0.8], 
                "ZB": [0.9 , 0.9 , 0.9 , 0.9 , 0.9 , 0.9], 
                "ZC": [1.3 , 1.3 , 1.2 , 1.2 , 1.2 , 1.2],
                "ZD": [1.6 , 1.4 , 1.2 , 1.1 , 1.0 , 1.0],
                "ZE": [2.4 , 1.7 , 1.3 , 1.1 , 0.9 , 0.8]}

    S1_range = [0.10 , 0.20 , 0.30, 0.40 , 0.50 , 0.60]
    F1_table = {"ZA": [0.8 , 0.8 , 0.8 , 0.8 , 0.8 , 0.8], 
                "ZB": [0.8 , 0.8 , 0.8 , 0.8 , 0.8 , 0.8], 
                "ZC": [1.5 , 1.5 , 1.5 , 1.5 , 1.5 , 1.4],
                "ZD": [2.4 , 2.2 , 2.0 , 1.9 , 1.8 , 1.7],
                "ZE": [4.2 , 3.3 , 2.8 , 2.4 , 2.2 , 2.0]}
    
    def __init__(self, Ss, S1, site_class):
        
        self.Ss = Ss
        self.S1 = S1
        self.site_class = site_class
        
    def calculate_Fs(self):
        if self.Ss <= self.Ss_range[0]:
            Fs = self.Fs_table[self.site_class][0]
            
        elif self.Ss >= self.Ss_range[-1]:
            Fs = self.Fs_table[self.site_class][-1]
        
        else:
            Fs_s = interp1d(self.Ss_range, self.Fs_table[self.site_class], kind='linear')
            Fs = Fs_s(self.Ss)
        return Fs
        
    def calculate_F1(self):
        if self.S1 <= self.S1_range[0]:
            F1 = self.F1_table[self.site_class][0]
            
        elif self.S1 >= self.S1_range[-1]:
            F1 = self.F1_table[self.site_class][-1]
        
        else:
            F1_s = interp1d(self.S1_range, self.F1_table[self.site_class], kind='linear')
            F1 = F1_s(self.S1)
        return F1
    
    def get_Sds(self):
        return round(self.Ss*self.calculate_Fs(),3)
    
    def get_Sd1(self):
        return round(self.S1*self.calculate_F1(),3)
    
    def Sae_(self):
        Ta = 0.2*(self.get_Sd1()/self.get_Sds())
        Tb = self.get_Sd1()/self.get_Sds()
        Tl = 6
        
        self.T_list = list(np.arange(0,8.01,0.01))
        self.Sae_list = []

        
        for a in self.T_list:
            if a <= Ta:
                Sae = (0.4+0.6*(a/Ta))*self.get_Sds()
                self.Sae_list.append(Sae)
            elif a <= Tb:
                Sae = self.get_Sds()
                self.Sae_list.append(Sae)
            elif a <= Tl:
                Sae = self.get_Sd1()/a
                self.Sae_list.append(Sae)
            else:
                Sae = self.get_Sd1()*Tl / (a**2)
                self.Sae_list.append(Sae)
        
        return self.T_list, self.Sae_list
            
    def plot(self):
        self.Sae_()
        plt.plot(self.T_list, self.Sae_list, label="Sae")
        plt.title("Yatay Elastik İvme Spektrumu", fontsize = 15)
        plt.grid()
        plt.xlabel("Periyot(s)")
        plt.ylabel("Sae(g)")
        plt.xlim(left=0, right=8)        
        plt.ylim(bottom=0)
        plt.legend()
        
# %%

Ss = 1.68
S1 = 0.63
site_class = "ZB"

example = Tbdy_spektrum(Ss, S1, site_class)

example.plot()

