
import numpy as np


class AWTLS_Data:
# ----------------------------------------------------------------------------------------------------------------
# This is a class containing necessary data to proceed cell capacity estimate through Q_Est_AWTLS function
# ----------------------------------------------------------------------------------------------------------------
    
    def __init__(self):
        
        self.c1 = None               # c1_n = c1_n-1 + x_n^2/Vary
        self.c2 = None               # c1_n = c1_n-1 + x_n^2/Vary
        self.c3 = None               # c3_n = c2_n-1 + y_n^2/Vary
        self.c4 = None               # c4_n = c4_n-1 + x_n^2/Varx
        self.c5 = None               # c5_n = c5_n-1 + x_n*y_n/Varx
        self.c6 = None               # c6_n = c6_n-1 + y_n^2/Varx
        self.gamma = None            # fading memory factor 
        
    def setup(self,Q,Varx,Vary,gamma):
        # Initiate the data structure
        
        if Q==0:
            self.c1 = 0        
            self.c2 = 0      
            self.c3 = 0    
            self.c4 = 0               
            self.c5 = 0            
            self.c6 = 0 
        else:
            K = (Varx/Vary)**0.5
            self.c1 = 1**2/Vary/K**2         
            self.c2 = 1*Q/Vary/K      
            self.c3 = Q**2/Vary     
            self.c4 = 1**2/Varx               
            self.c5 = 1*Q/Varx*K            
            self.c6 = Q**2/Varx*K**2         

        self.gamma = gamma
        
def Q_Est_AWTLS(x,y,Varx,Vary,AWTLSData):
# ----------------------------------------------------------------------------------------------------------------
# This function estimates cell capcity through AWTLS.
# Input: SOC change(x), integral current(y), variance of SOC change(Varx), variance of integral current(Vary), 
#        AWTLSdata (AWTLSData)
# Output: cell capcity estimate(Q_hat), variance of capacity estimate(VarQ), AWTLSdata Update(AWTLSData)
# ----------------------------------------------------------------------------------------------------------------
    
    # Load AWLTSdata
    gamma = AWTLSData.gamma
    c1 = AWTLSData.c1
    c2 = AWTLSData.c2
    c3 = AWTLSData.c3
    c4 = AWTLSData.c4
    c5 = AWTLSData.c5
    c6 = AWTLSData.c6
    
    # Update AWTLS paramters 
    K = (Varx/Vary)**0.5
    c1 = c1*gamma + x**2/Vary/K**2
    c2 = c2*gamma + x*y/Vary/K
    c3 = c3*gamma + y**2/Vary
    c4 = c4*gamma + x**2/Varx
    c5 = c5*gamma + x*y/Varx*K
    c6 = c6*gamma + y**2/Varx*K**2
    
    # Find Q that minimizes cost function
    r = np.roots([c5,-c1+2*c4-c6,3*c2-3*c5,c1-2*c3+c6,-c2]) # roots of the Jacobian of cost function
    r = r[np.imag(r)==0]
    r = np.real(r)
    r = r[r>0]
    cost = (1/(r**2+1)**2) * (c4*r**4-2*c5*r**3+(c1+c6)*r**2-2*c2*r+c3) # cost corresponding to each root
    r = r[cost == min(cost)][0]
    Q_hat = r/K
    
    # Find the corresponding Hessian of cost function
    H = (2/(r**2+1)**4) * (-2*c5*r**5+(3*c1-6*c4+3*c6)*r**4+(-12*c2+16*c5)*r**3+\
                           (-8*c1+10*c3+6*c4-8*c6)*r**2+(12*c2-6*c5)*r+(c1-2*c3+c6))
    
    # Calcualte the variance of Q estimate
    VarQ = 2/H/K**2
    
    # Update AWLTSdata
    AWTLSData.c1 = c1
    AWTLSData.c2 = c1
    AWTLSData.c3 = c3
    AWTLSData.c4 = c4
    AWTLSData.c5 = c5
    AWTLSData.c6 = c6
    
    return Q_hat, VarQ, AWTLSData