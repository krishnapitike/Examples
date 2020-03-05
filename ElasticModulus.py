import numpy as np
def computeEl(c11,c12,c44):
    c   = np.array([ [c11,c12,c12,0,0,0], 
                     [c12,c11,c12,0,0,0], 
                     [c12,c12,c11,0,0,0],
                     [0,0,0,c44,0,0],
                     [0,0,0,0,c44,0],
                     [0,0,0,0,0,c44] ])
    s   = np.linalg.inv(c)
    s11 = s[0,0]
    s12 = s[0,1]
    s44 = s[3,3]
    KV  = ( c11+2*c12 )/3           #Voight average of bulk modulus
    GV  = c11/5 - c12/5 + 3*c44/5   #Voight average of shear modulus
    KR  = 1/( 3*s11+6*s12 )         #Reiss average of bulk modulus
    GR  = 5/( 4*s11 -4*s12 +3*s44 ) #Reiss average of shear modulus
    K   = (KV+KR)/2                 #Voight-Reiss-Hill average of bulk modulus
    G   = (GV+GR)/2                 #Voight-Reiss-Hill average of shear modulus
    E   = 1/(1/(3*G) + 1/(9*K))     #Voight-Reiss-Hill average of young's modulus
    nu  = (1-(3*G/(3*K+G)))/2       #Poisson's ratio
    AZ  = 2*c44/(c11-c12)           #Zenner anisotropy factor
    AU  = 5*GV/GR + KV/KR - 6       #Universal anisotropy factor
    return(K,G,E,nu,AZ,AU)

#set c11,c12,c44
#c11,c12,c44=211,121,109       # NiO exp at 0K Plessis
#c11,c12,c44=224,97,110        # NiO exp at 298K Plessis
#c11,c12,c44=271,125,105       # NiO exp at 298 K Uchida and Saito
#c11,c12,c44=450,163,163       # NiO Towler semi imperical model
#c11,c12,c44=344.6, 141.3, 40  # NiO exp at 298 K Wang
#c11,c12,c44=365.9,81.9,77.5   # NiO cinthia DFT U=0
#c11,c12,c44=279.6,231.9,126.3 # NiO cinthia DFT U=8
#c11,c12,c44=260.3,163.8,40.1  # CoO cinthia DFT U=0
#c11,c12,c44=212.8,226.7,130.3 # CoO cinthia DFT U=8
c11,c12,c44=266.4, 147.2,83.6
print("c11     c12      c44    K_VRH    G_VRH  E_VRH     nu      A_Z     A_U")
print("%7.3f %7.3f %7.3f"%(c11, c12, c44), end = ' ')
print("%7.3f %7.3f %7.3f %7.3f %7.3f %7.3f"%(computeEl(c11,c12,c44)))
