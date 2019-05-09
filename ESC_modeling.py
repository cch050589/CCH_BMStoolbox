
import numpy as np
import pandas as pd


class ESCmodel:
# ----------------------------------------------------------------------------------------------------------------
# Class members are parameters and relationships betwwen SOC and OCV that describe cell behaviors.  Member 
# fuctions are methods that can retrive model parametes, calcualte OCV, calculate SOC and simulate cell voltage.
# ----------------------------------------------------------------------------------------------------------------

    def __init__(self):
        
        self.temps = None       # temperature
        self.numpoles = None    # number of parallel R-C subcircuits
        self.etaParam = None    # Coulombic efficiency
        self.QParam = None      # cell capacity
        self.GParam = None      # tuning rate of dynamic hysteresis  
        self.MParam = None      # magnitude of dynamic hysteresis
        self.M0Param = None     # magnitude of instantaneous hysteresis
        self.R0Param = None     # ohmic resistance
        self.RCParam = None     # RC time constant of parallel resistor-capacitor sub-circuit(s)
        self.RParam = None      # resistance of resistor in parallel resistor-capacitor sub-circuit(s)
        
        self.SOC = None         #
        self.OCV0 = None        # OCV-SOC relationship
        self.OCVrel = None      #

        self.OCV = None         # 
        self.SOC0 = None        # SOC-OCV relationship
        self.SOCrel = None      # 

    
    def getParam(self,param,T):
    # -----------------------------------------------------------------------------------------
    # This member function retrives a model parameter
    # Input: parameter name(param), temperature(T)
    # Output: parameter
    # -----------------------------------------------------------------------------------------
        
        x = self.temps.index(T)
        
        if param == 'numpoles':
            return self.numpoles
        elif param == 'etaParam':
            return self.etaParam[x]
        elif param == 'QParam': 
            return self.QParam[x]
        elif param == 'GParam':
            return self.GParam[x]
        elif param == 'M0Param':
            return self.M0Param[x]
        elif param == 'MParam':
            return self.MParam[x]
        elif param == 'R0Param':
            return self.R0Param[x]
        elif param =='RCParam':
            return self.RCParam[x]
        else:
            return self.RParam[x]
    
    
    def OCVfromSOCtemp(self,soc,temp):
    # -----------------------------------------------------------------------------------------
    # This member function calculates OCV
    # Input: state of charge(soc), temperature(temp)
    # Output: open circuit voltage(ocv)
    # -----------------------------------------------------------------------------------------
    
        if np.isscalar(soc):
            soc = [soc]
        if np.isscalar(temp):
            temps = [temp for x in range(len(soc))]
    
        ocv = np.zeros(len(soc))
        deltaSOC = self.SOC[1] - self.SOC[0]

        for x,T in enumerate(temps):
            z = soc[x]
            if z >= 1:
                indub = -1
                indlb = -2
            elif z <= 0:
                indub = 1
                indlb = 0
            else:
                indlb = int((z - self.SOC[0])/deltaSOC)
                indub = indlb + 1
        
            vub = (self.OCV0[indub] + T*self.OCVrel[indub])
            vlb = (self.OCV0[indlb] + T*self.OCVrel[indlb])
            deltav = vub - vlb
            ocv[x] = vlb + (z - self.SOC[indlb])*(deltav/deltaSOC)

        return ocv
    
        
    def SOCfromOCVtemp(self,ocv,temp):
    # -----------------------------------------------------------------------------------------
    # This member function calculates SOC
    # Input: open circuit voltage(ocv), temperature(temp)
    # Output: state of charge(soc)
    # -----------------------------------------------------------------------------------------
    
        if np.isscalar(ocv):
            ocv = [ocv]
        if np.isscalar(temp):
            temps = [temp for x in range(len(ocv))]
    
        soc = np.zeros(len(ocv))
        deltaOCV = self.OCV[1] - self.OCV[0]

        for x,T in enumerate(temps):
            v = ocv[x]
            if v >= self.OCV[-1]:
                indub = -1
                indlb = -2
            elif v <= self.OCV[0]:
                indub = 1
                indlb = 0
            else:
                indlb = int((v - self.OCV[0])/deltaOCV)
                indub = indlb + 1
        
            zub = (self.SOC0[indub] + T*self.SOCrel[indub])
            zlb = (self.SOC0[indlb] + T*self.SOCrel[indlb])
            deltaz = zub - zlb
            soc[x] = zlb + (v - self.OCV[indlb])*(deltaz/deltaOCV)

        return soc
    
    
    def simCell(self,i,T,delta_t,z0,iR0,h0):
    # -----------------------------------------------------------------------------------------
    # This member function estimates cell voltage
    # Input: current(i), temperature(T), sampling intergval(delta_t), initial states(z0,iR0,h0)
    # Output: cell voltage(v), state of charge(z), open circuit voltage(OCV)
    # -----------------------------------------------------------------------------------------
        
        x = self.temps.index(T)
        eta = self.etaParam[x]
        Q = self.QParam[x]
        G = self.GParam[x]
        M0= self.M0Param[x]
        M = self.MParam[x]
        R0= self.R0Param[x]
        RC = self.RCParam[x]
        R = self.RParam[x]
        
        # Compute SOC and correpslidng OCV
        etai = i
        for k in range(len(i)):
            if i[k] < 0: etai[k] = etai[k]*eta
        z = np.ones(len(etai))*z0 - np.cumsum(np.concatenate(([0],etai[:-1])))*delta_t/(Q*3600)
        OCV = self.OCVfromSOCtemp(z,T)
        
        # Compute other dyanmic states
        h = np.zeros(len(i)) 
        h[0] = h0
        si = np.zeros(len(i)) 
        AH = np.exp(-abs(G*etai*delta_t/(3600*Q)))
        iR = np.zeros((self.numpoles,len(i)))
        iR[:,0] = iR0
        RCfact = np.exp(-delta_t/RC)
        for k in np.arange(1,len(i),1):
            h[k] = AH[k-1]*h[k-1] - (1-AH[k-1])*np.sign(etai[k-1])
            si[k] = np.sign(etai[k])
            if abs(etai[k]) < Q/100: si[k] = si[k-1]
            iR[:,k] = np.diag(RCfact).dot(iR[:,k-1])+(1-RCfact).dot(etai[k-1])
        
        # Compute cell voltage
        v = OCV + M*h + M0*si - R0*etai - R.dot(iR)
        
        return v, z, OCV