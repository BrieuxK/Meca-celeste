import numpy as np
import matplotlib.pyplot as plt

tmax = 5000 * 12 #mois
k = 24 #mois
n_k = int(tmax/k)

G = 1.15488171e-4
#G = 2.959e-4

mass1 = 0.000954786104043          #1.898*1e27   #Masse Jupiter
mass2 = 1                          #1.989*1e30   #Masse Soleil

masses = [1, 0.000954786104043, 0.0002857214681, 0.0000436450477]

def dq(pos, imp):
    zq = np.zeros(3*len(masses))
    for i in range(len(masses)):
        j = i*3
        zq[j] = imp[j]/masses[i]
        zq[j+1] = imp[j+1]/masses[i]
        zq[j+2] = imp[j+2]/masses[i]
    return zq

def dq1(param_pos, param_imp): #Explose dés la 1ère itération
    """
    param_pos = array : [qX1,qY1,qZ1,qX2,qY2,qZ2,qX3,qY3,qZ3] 3premiers : Soleil, 3suivants Jupiter, 3 derniers Saturne
    param_imp = array : [pX1,pY1,pZ1,pX2,pY2,pZ2,pX3,pY3,pZ3] 3premiers : Soleil, 3suivants Jupiter, 3 derniers Saturne

    """
    zq = np.zeros(3*len(masses))
    zq[0] = param_imp[0]/masses[0]
    zq[1] = param_imp[1]/masses[0]
    zq[2] = param_imp[2]/masses[0]
    zq[3] = param_imp[3]/masses[1]
    zq[4] = param_imp[4]/masses[1]
    zq[5] = param_imp[5]/masses[1]
    zq[6] = param_imp[6]/masses[2]
    zq[7] = param_imp[7]/masses[2]
    zq[8] = param_imp[8]/masses[2]
    zq[9] = param_imp[9]/masses[3]
    zq[10] = param_imp[10]/masses[3]
    zq[11] = param_imp[11]/masses[3]
    """
    zq : array des dq/dt longeur 6 
    """
    #print(zq)
    return zq

def dp(param_pos,param_imp):
    zp = np.zeros(len(param_pos))
    for i in range(len(masses)):
        j = i*3
        for m in range(len(masses)):
            n = m*3
            if i != m :
                zp[j]  += -2*((G*masses[i]*masses[m])/distance([param_pos[j],param_pos[j+1],param_pos[j+2]], [param_pos[n],param_pos[n+1],param_pos[n+2]])**3) * (param_pos[j]-param_pos[n])
                zp[j+1]+= -2*((G*masses[i]*masses[m])/distance([param_pos[j],param_pos[j+1],param_pos[j+2]], [param_pos[n],param_pos[n+1],param_pos[n+2]])**3) * (param_pos[j+1]-param_pos[n+1])
                zp[j+2]+= -2*((G*masses[i]*masses[m])/distance([param_pos[j],param_pos[j+1],param_pos[j+2]], [param_pos[n],param_pos[n+1],param_pos[n+2]])**3) * (param_pos[j+2]-param_pos[n+2])
    return zp

def dp1(param_pos, param_imp): #On divise par distance XYZ
    zp = np.zeros(len(param_pos))
    
    zp[0] = -2*(((G*masses[0]*masses[1])/distance(param_pos[:3],param_pos[3:6])**3) * (param_pos[0] - param_pos[3])) - 2*(((G*masses[0]*masses[2])/distance(param_pos[:3],param_pos[6:])**3) * (param_pos[0] - param_pos[6]))#Marqueur 
    zp[1] = -2*(((G*masses[0]*masses[1])/distance(param_pos[:3],param_pos[3:6])**3) * (param_pos[1] - param_pos[4])) - 2*(((G*masses[0]*masses[2])/distance(param_pos[:3],param_pos[6:])**3) * (param_pos[1] - param_pos[7]))
    zp[2] = -2*(((G*masses[0]*masses[1])/distance(param_pos[:3],param_pos[3:6])**3) * (param_pos[2] - param_pos[5])) - 2*(((G*masses[0]*masses[2])/distance(param_pos[:3],param_pos[6:])**3) * (param_pos[2] - param_pos[8]))
    
    zp[3] = -2*(((G*masses[0]*masses[1])/distance(param_pos[:3],param_pos[3:6])**3) * (param_pos[3] - param_pos[0])) - 2*(((G*masses[1]*masses[2])/distance(param_pos[3:6],param_pos[6:])**3) * (param_pos[3] - param_pos[6]))
    zp[4] = -2*(((G*masses[0]*masses[1])/distance(param_pos[:3],param_pos[3:6])**3) * (param_pos[4] - param_pos[1])) - 2*(((G*masses[1]*masses[2])/distance(param_pos[3:6],param_pos[6:])**3) * (param_pos[4] - param_pos[7]))
    zp[5] = -2*(((G*masses[0]*masses[1])/distance(param_pos[:3],param_pos[3:6])**3) * (param_pos[5] - param_pos[2])) - 2*(((G*masses[1]*masses[2])/distance(param_pos[3:6],param_pos[6:])**3) * (param_pos[5] - param_pos[8]))
    
    zp[6] = -2*(((G*masses[2]*masses[0])/distance(param_pos[6:],param_pos[:3])**3) * (param_pos[6] - param_pos[0])) - 2*(((G*masses[1]*masses[2])/distance(param_pos[3:6],param_pos[6:])**3) * (param_pos[6] - param_pos[3]))
    zp[7] = -2*(((G*masses[2]*masses[0])/distance(param_pos[6:],param_pos[:3])**3) * (param_pos[7] - param_pos[1])) - 2*(((G*masses[1]*masses[2])/distance(param_pos[3:6],param_pos[6:])**3) * (param_pos[7] - param_pos[4]))
    zp[8] = -2*(((G*masses[2]*masses[0])/distance(param_pos[6:],param_pos[:3])**3) * (param_pos[8] - param_pos[2])) - 2*(((G*masses[1]*masses[2])/distance(param_pos[3:6],param_pos[6:])**3) * (param_pos[8] - param_pos[5]))
    
    """
    zp : array des dp/dt longueur 9
    """
    #print(zp)
    return zp

def distance(pos1, pos2):
    """
    pos1,pos2 = [X,Y,Z]
    """
    dist = np.sqrt(((pos1[0]-pos2[0])**2)+((pos1[1]-pos2[1])**2)+((pos1[2]-pos2[2])**2))
    return dist

def F_Heun(position, impulsion):
    """
    position,impulsion : array [x1 y1 z1 x2 y2 z2]
    """
    q1 = k * dq(position,impulsion)
    #print("q1 :",q1)
    p1 = k * dp(position,impulsion)
    q_tilde = position + q1
    #print("q_tilde :", q_tilde)
    p_tilde = impulsion + p1
    q2 = k * dq(q_tilde, p_tilde)
    p2 = k * dp(q_tilde, p_tilde)
    return 0.5*(q1+q2), 0.5*(p1+p2)

def Heun(q_Heun, p_Heun):
    """
    q/p_Heun : matrice (n_k,6)
    CHaque colonne représente une donnée 
    Chaque ligne représente un pas de temps
    """
    for i in range(n_k-1):
        delta_q, delta_p = F_Heun(q_Heun[i],p_Heun[i])
        q_Heun[i+1] = q_Heun[i] + delta_q
        p_Heun[i+1] = p_Heun[i] + delta_p
    return q_Heun, p_Heun

def SV(q_SV,p_SV):
    for i in range(n_k-1):
        p_mid = p_SV[i] + 0.5 * k * dp(q_SV[i], p_SV[i])
        q_SV[i+1] = q_SV[i] + k * dq(q_SV[i], p_mid)
        p_SV[i+1] = p_mid + 0.5 * k * dp(q_SV[i+1],p_mid)
    return q_SV,p_SV

def Energy(result_q, result_p): 
    energies = np.zeros(n_k)
    for i in range(n_k):       
        T_Sun = (result_p[i,0]**2 + result_p[i,1]**2 + result_p[i,2]**2)/2*masses[0]
        T_Jup = (result_p[i,3]**2 + result_p[i,4]**2 + result_p[i,5]**2)/2*masses[1]
        T_Sat = (result_p[i,6]**2 + result_p[i,7]**2 + result_p[i,8]**2)/2*masses[2]
        V_SoJ = G*masses[0]*masses[1]/distance(result_q[i,:3],result_q[i,3:6])
        V_SoSat = G*masses[0]*masses[2]/distance(result_q[i,:3],result_q[i,6:])
        V_JSat = G*masses[1]*masses[2]/distance(result_q[i,3:6],result_q[i,6:])
        energies[i] = (T_Jup + T_Sun + T_Sat) - V_SoJ - V_SoSat - V_JSat
    return energies

def angular_m(result_q, result_p, coord): #Moment angulaire
    """
    result_q/_p => (n_k,6) : matrice de solutions finales
    coord => 0 pour x, 1 pour y, 2 pour z
    """
    ang_m = np.zeros(n_k)
    ang_m = result_q[:,0 + coord] * result_p[:,0 + coord] + result_q[:,3 + coord] * result_p[:,3 + coord]
    """
    return : moment angulaire par rapport à une seule coordonnée
    """
    return ang_m

def cm(result_q, coord):  #Centre de masse
    """
    coord => 0 pour x, 1 pour y, 2 pour z
    """
    cmass = np.zeros(n_k)
    cmass = 1 * result_q[:,0 + coord] + 0.000954786104043 * result_q[:,3 + coord] + 0.0002857214681 * result_q[:,6 + coord] + 0.0000436450477 * result_q[:,9 + coord]
    """
    return : centre de masse par rapport à une seule coordonnée
    """
    return cmass/(1 + 0.000954786104043 + 0.0002857214681 + 0.0000436450477)
#--------------------------------------------------------------------------
q_SV = np.zeros((n_k, len(masses) * 3))  #1st col : XSun, 2nd : YSun, 3rd : ZSun,...., 6th : ZJup
p_SV = np.zeros((n_k, len(masses) * 3))

q_SV[0] = np.array([0, 0, 0,4.6580839119399,-1.7945574704081,-0.0967627432056,6.9600794880566,-7.0666123077577,-0.1541266614081,14.3976311704513,13.4787169790607,-0.1365126314185])
p_SV[0] = np.array([0, 0, 0,0.0026255239604*masses[1],0.0074041879836*masses[1],-0.0000894925649*masses[1],0.0036673089657*masses[2],0.0039101598352*masses[2],-0.0002138877692*masses[2],-0.0027147204414*masses[3],0.0026953325910*masses[3],0.0000450465708*masses[3]])

"""
q_SV[0] = np.array([0, 0, 0,4.8320528131341,-1.1900491790446,-0.1031658265828,7.2530549979978,-6.7394166060155,-0.1714717076894])
p_SV[0] = np.array([0, 0, 0,0.0017167836330*masses[1],0.0076897387172*masses[1],-0.0000703506472*masses[1],0.0034894004389*masses[2],0.0040828028757*masses[2],-0.0002098149920*masses[2]])
"""

q_Heun = np.zeros((n_k, len(masses) * 3))  #1st col : XSun, 2nd : YSun, 3rd : ZSun,...., 6th : ZJup
p_Heun = np.zeros((n_k, len(masses) * 3))

#COnditions initiales dans la 1ère ligne de chaque matrice ##########
q_Heun[0] = np.array([0, 0, 0,4.6580839119399,-1.7945574704081,-0.0967627432056,6.9600794880566,-7.0666123077577,-0.1541266614081,14.3976311704513,13.4787169790607,-0.1365126314185])
p_Heun[0] = np.array([0, 0, 0,0.0026255239604*masses[1],0.0074041879836*masses[1],-0.0000894925649*masses[1],0.0036673089657*masses[2],0.0039101598352*masses[2],-0.0002138877692*masses[2],-0.0027147204414*masses[3],0.0026953325910*masses[3],0.0000450465708*masses[3]])


finalq, finalp = SV(q_SV,p_SV)
#finalq, finalp = Heun(q_Heun, p_Heun)
print(np.shape(finalq))

fig = plt.figure(figsize=plt.figaspect(0.5))
ax = fig.add_subplot(1, 2, 1, projection='3d')

ax.plot3D(finalq[:,9],finalq[:,10],finalq[:,11], color = "Orange") #Uranus
ax.plot3D(finalq[:,6],finalq[:,7],finalq[:,8], color = "Blue") #Saturne
ax.plot3D(finalq[:,3],finalq[:,4],finalq[:,5], color = "Red") #Jupiter
ax.plot3D(finalq[:,0],finalq[:,1],finalq[:,2], color = "Black") #Soleil
plt.title('n')
plt.show()

"""
ax.plot3D(finalq[:,9] - cm(finalq,0),finalq[:,10] - cm(finalq,1),finalq[:,11] - cm(finalq,2), color = "Orange") #Uranus
ax.plot3D(finalq[:,6] - cm(finalq,0),finalq[:,7] - cm(finalq,1),finalq[:,8] - cm(finalq,2), color = "Blue") #Saturne
ax.plot3D(finalq[:,3] - cm(finalq,0),finalq[:,4] - cm(finalq,1),finalq[:,5] - cm(finalq,2), color = "Red") #Jupiter
ax.plot3D(finalq[:,0] - cm(finalq,0),finalq[:,1] - cm(finalq,1),finalq[:,2] - cm(finalq,2), color = "Black") #Soleil
plt.show()
"""

"""
en = Energy(finalq, finalp)
plt.plot(range(n_k), en)
plt.show()
"""


