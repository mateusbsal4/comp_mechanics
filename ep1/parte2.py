import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


#Parâmetros (SI)
h = 0.15
L = 3
d = L/2
H = 2*L
V = 100/3.6
R = L/2
epsilon = 0.01
p_atm = 101325
k = 0.026
rho = 1.25
cp = 1002
gama = 1.4

def plot_psi(delta, graph): #plota a função de corrente
    fig = plt.figure()
    x = np.linspace(0, L + 2*d, int(2*L/delta)+1)
    y = np.linspace(0, H, int(2*L/delta)+1)
    X, Y = np.meshgrid(x, y)
    plt.pcolor(X, Y, graph, cmap = 'jet')
    plt.colorbar(label='Magnitude')
    plt.title('Função de corrente')
    plt.show()

def plot_temp(delta, graph): #plota a função de corrente
    fig = plt.figure()
    x = np.linspace(0, L + 2*d, int(2*L/delta)+1)
    y = np.linspace(0, H, int(2*L/delta)+1)
    X, Y = np.meshgrid(x, y)
    plt.pcolor(X, Y, graph, cmap = 'jet')
    plt.colorbar(label='Magnitude')
    plt.title('Função de temperatura')
    plt.show()


def plot_velfield(delta, u, v):        #plota o campo de velocidades
    x = np.linspace(0, L + 2 * d, int(2 * L / delta) + 1)
    y = np.linspace(0, H, int(2 * L / delta) + 1)
    X, Y = np.meshgrid(x, y)
    norm = np.sqrt(u**2 + v**2) 
    
    plt.figure()
    plt.quiver(X, Y, u, v, norm, cmap='viridis')
    plt.colorbar(label='Magnitude (m/s)') 
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Campo de velocidades')
    plt.show()

def plot_pressure(delta, graph): #plota a função de corrente
    fig = plt.figure()
    x = np.linspace(0, L + 2*d, int(2*L/delta)+1)
    y = np.linspace(0, H, int(2*L/delta)+1)
    X, Y = np.meshgrid(x, y)
    plt.pcolor(X, Y, graph, cmap = 'jet')
    plt.colorbar(label='Magnitude (Pa)')
    plt.title('Pressão no domínio')
    plt.show()




def plot_p_upper_contour(sorted_x_contour, sorted_p_contour):
    #print("X_contour", x_contour)
    #print("P_contour", p_contour)
    min_index = np.argmin(sorted_p_contour)
    min_x = sorted_x_contour[min_index]
    min_p = sorted_p_contour[min_index]
    #plt.scatter(x_contour, p_contour)
    print("sortedX_contour", sorted_x_contour)
    print("sortedP_contour", sorted_p_contour)
    plt.plot(sorted_x_contour, sorted_p_contour)
    plt.xlabel('x [m]')
    plt.ylabel('Pressão [Pa]')
    plt.title('Pressão ao longo da parte superior da carroceira')
    plt.grid(True)
    plt.annotate(f'Pressão mínima: ({min_x:.2f}, {min_p:.2f})', xy=(min_x, min_p),
             xytext=(min_x+0.4, min_p + 1), arrowprops=dict(arrowstyle='->'))
    #plt.show()

def plot_p_lower_contour(x_contour, p_contour):
    print("sortedX_contour", x_contour)
    print("sortedP_contour", p_contour)
    plt.plot(x_contour, p_contour)
    plt.xlabel('x [m]')
    plt.ylabel('Pressão [Pa]')
    plt.title('Pressão ao longo da parte superior da carroceira')
    plt.grid(True)
    #plt.show()

def isInsideCar(x, y, inclusive = 0):
    if(inclusive):
        return (x-d-L/2)**2+(y-h)**2<=R**2 and y>=h
        #return (y <= np.sqrt((L/2)**2 - (x-d-L/2)**2) + h) and y >= h
    return (x-d-L/2)**2+(y-h)**2<R**2 and y > h
    




def mdf_psi(delta, lamb):
    n_lines = int(2*L/delta)+1        
    n_columns = int(L/delta)+1 
    assert n_lines%2 != 0         #escolher n.o de linhas ímpar devido a simetria   
    psi = np.full((n_lines, n_columns), 0, dtype=float)        #precisa estar entre 0 e 1
    print(psi)
    max_error = 1       #inicialmente psi e psi velho são iguais
    while(max_error>epsilon):
        psi_old = np.copy(psi)
        #print("Psi old antes", psi_old)
        #for i in range(psi.shape[0] - 1, -1, -1):
        for i in range(psi.shape[0]): # linhas, iteração de baixo para cima
            for j in range(psi.shape[1]): # colunas, iteração da esquerda para a direita
                if (i == psi.shape[0]-1 and j == 0): # canto superior
                    psi[i][j] = 0.5*(psi[i][j+1]+psi[i-1][j]+V*delta)
                elif (i == 0 and j == 0): # canto inferior
                    psi[i][j] = 0
                elif (j == 0): # Parede (esquerda) do dominio
                    psi[i][j] = (2*psi[i][1]+psi[i+1][0]+psi[i-1][0])/4
                elif (i == psi.shape[0]-1): # Teto do dominio
                    if (j == psi.shape[1]-1): # Caso esteja no MEIO do dominio (CUIDADO)
                        psi[i][j] = 0.25*(psi[i][j-1]+psi[i][j-1]+2*psi[i-1][j]+2*V*delta)
                    else:
                        psi[i][j] = 0.25*(psi[i][j+1]+psi[i][j-1]+2*psi[i-1][j]+2*V*delta)
                elif (i == 0): # Chao do domínio
                    psi[i][j] = 0
                elif isInsideCar(j*delta, i*delta, 1): #dentro do carro
                    psi[i][j] = 0
                elif isInsideCar((j+1)*delta, i*delta) and isInsideCar(j*delta, (i-1)*delta):       #contorno a direita e embaixo
                    n_center = n_columns - j -2                    #numero de deltas (inteiros) dentro do circulo na direção horizontal até o centro
                    y_center = i*delta - h                  #altura do pto do contorno (e do pto a sua esquerda) ate o centro
                    bij = (delta - (np.sqrt(R**2-(y_center)**2)-n_center*delta))/delta
                    assert bij<1
                    n_center = n_columns - j -1 #numero de deltas (inteiros) dentro do circulo na direção horizontal até o centro
                    y_center = (i-1)*delta - h   #altura do pto do contorno ate o centro
                    aij = (delta - (np.sqrt(R**2-(n_center*delta)**2)-y_center))/delta  #pitagoras
                    assert aij<1
                    psi[i][j] = ((aij*bij)/(aij+bij))*((psi[i][j-1]/(bij+1))+(psi[i+1][j])/(aij+1))
                elif isInsideCar((j+1)*delta, i*delta): #Contorno está a direita
                    n_center = n_columns - j -2                    #numero de deltas (inteiros) dentro do circulo na direção horizontal até o centro
                    y_center = i*delta - h                  #altura do pto do contorno (e do pto a sua esquerda) ate o centro
                    bij = (delta - (np.sqrt(R**2-(y_center)**2)-n_center*delta))/delta  #pitagoras
                    assert bij<1
                    psi[i][j] = (bij/(2*(bij+1)))*(psi[i+1][j]+psi[i-1][j]+(2/(bij+1))*psi[i][j-1])
                elif isInsideCar(j*delta, (i-1)*delta): #Contorno está abaixo
                    n_center = n_columns - j - 1 #numero de deltas (inteiros) dentro do circulo na direção horizontal até o centro
                    y_center = (i-1)*delta - h   #altura do pto do contorno ate o centro
                    aij = (delta - (np.sqrt(R**2-(n_center*delta)**2)-y_center))/delta  #pitagoras
                    assert aij<1
                    if (j == psi.shape[1]-1): # Caso esteja no MEIO do dominio (CUIDADO)
                        psi[i][j] = (aij/(2*(aij+1)))*(psi[i][j-1]+psi[i][j-1]+(2/(aij+1))*psi[i+1][j])
                    else:
                        psi[i][j] = (aij/(2*(aij+1)))*(psi[i][j+1]+psi[i][j-1]+(2/(aij+1))*psi[i+1][j])
                elif isInsideCar(j*delta, (i+1)*delta, 1): #Contorno acima
                    c = (h-i*delta)/delta
                    if (j == psi.shape[1]-1): # Caso esteja no MEIO do dominio (CUIDADO)
                        psi[i][j] = (c/(2*(c+1)))*(psi[i][j-1]+psi[i][j-1]+(2/(c+1))*psi[i-1][j])
                    else:
                        psi[i][j] = (c/(2*(c+1)))*(psi[i][j+1]+psi[i][j-1]+(2/(c+1))*psi[i-1][j])
                else: # Parte central do dominio
                    if (j == psi.shape[1]-1): # Caso esteja no MEIO do dominio (CUIDADO)
                        psi[i][j] = 0.25*(psi[i+1][j]+psi[i-1][j]+psi[i][j-1]+psi[i][j-1])
                    else:
                        psi[i][j] = 0.25*(psi[i+1][j]+psi[i-1][j]+psi[i][j+1]+psi[i][j-1])
                psi[i][j] = lamb*psi[i][j]+(1-lamb)*psi_old[i][j]
        max_error = np.max(np.abs((psi-psi_old)))
        print("PSI", psi)
        print("psi old", psi_old)
        print("Max error ", max_error)
    return np.concatenate((psi, np.flip(psi, 1)[:, 1:]), axis=1)



def mdf_vel(delta, psi):
    n_lines = int(2*L/delta)+1        
    n_columns = n_lines
    u = np.zeros((n_lines, n_columns), dtype = float) 
    v = np.zeros((n_lines, n_columns), dtype = float)  
    for i in range(u.shape[0]): # linhas, iteração de baixo para cima
        for j in range(u.shape[1]): # colunas, iteração da esquerda para a direita
            if (i == psi.shape[0]-1 and (j == 0 or j == psi.shape[1]-1)): # cantos superiores 
                u[i][j] = V
                v[i][j] = 0        
            elif (i == 0 and (j == 0 or j ==psi.shape[1]-1)): # cantos inferiores
                u[i][j] = (-psi[i+2][j]+4*psi[i+1][j]-3*psi[i][j])/(2*delta)    #prog em y
                v[i][j] = 0
            elif (j == 0 or j == psi.shape[1]-1): # Paredes do dominio
                u[i][j] = (psi[i+1][j]-psi[i-1][j])/(2*delta)
                v[i][j] = 0
            elif (i == psi.shape[0]-1): # Teto do dominio
                u[i][j] = V
                v[i][j] = -(psi[i][j+1]-psi[i][j-1])/(2*delta)
            elif (i == 0): # Chao do domínio
                u[i][j] = (-psi[i+2][j]+4*psi[i+1][j]-3*psi[i][j])/(2*delta)    #prog em y
                v[i][j] = -(psi[i][j+1]-psi[i][j-1])/(2*delta)   
            elif isInsideCar(j*delta, i*delta, 1):      #dentro do carro
                u[i][j] = 0
                v[i][j] = 0
            elif isInsideCar((j+1)*delta, i*delta, 1) and isInsideCar(j*delta, (i-1)*delta, 1):       #contorno a direita e abaixo
                u[i][j] = (-psi[i+2][j]+4*psi[i+1][j]-3*psi[i][j])/(2*delta)            #progressiva em y 
                v[i][j] =  -(3*psi[i][j]-4*psi[i][j-1]+psi[i][j-2])/(2*delta)           #reg em x
            elif isInsideCar((j-1)*delta, i*delta, 1) and isInsideCar(j*delta, (i-1)*delta, 1):       #contorno a esquerda e abaixo
                u[i][j] = (-psi[i+2][j]+4*psi[i+1][j]-3*psi[i][j])/(2*delta)        #prog nas 2 direções
                v[i][j] =  -(-psi[i][j+2]+4*psi[i][j+1]-3*psi[i][j])/(2*delta)
            elif isInsideCar(j*delta, (i-1)*delta, 1): #Contorno está abaixo
                u[i][j] = (-psi[i+2][j]+4*psi[i+1][j]-3*psi[i][j])/(2*delta)    #prog em y       
                v[i][j] = -(psi[i][j+1]-psi[i][j-1])/(2*delta)  
            elif isInsideCar((j+1)*delta, i*delta, 1): #Contorno está a direita
                u[i][j] = (psi[i+1][j]-psi[i-1][j])/(2*delta)                        
                v[i][j] =  -(3*psi[i][j]-4*psi[i][j-1]+psi[i][j-2])/(2*delta)      #reg em x           
            elif isInsideCar((j-1)*delta, i*delta, 1): #Contorno está a esquerda
                u[i][j] = (psi[i+1][j]-psi[i-1][j])/(2*delta)                           
                v[i][j] = -(-psi[i+2][j]+4*psi[i+1][j]-3*psi[i][j])/(2*delta)            #prog em x
            elif isInsideCar(j*delta, (i+1)*delta, 1): #Contorno está acima
                u[i][j] = (3*psi[i][j]-4*psi[i-1][j]+psi[i-2][j])/(2*delta)      #reg em y                     
                v[i][j] = -(psi[i][j+1]-psi[i][j-1])/(2*delta) 
            else: # Parte central do dominio
                u[i][j] = (psi[i+1][j]-psi[i-1][j])/(2*delta)     
                v[i][j] = -(psi[i][j+1]-psi[i][j-1])/(2*delta)   
    return u, v

def mdf_pressure(u, v):
    p = p_atm+0.5*rho*((gama-1)/gama)*(V**2-(u**2+v**2))
    return p

def mdf_p_upper_contour(delta, p):
    n_rows, n_cols = p.shape
    inContour = np.zeros((n_rows, n_cols), dtype=bool)    
    p_contour = []
    x_contour = []   
    for i in range(n_rows):
        for j in range(n_cols):
            #if (not isInsideCar(j*delta, i*delta, 1) and (isInsideCar(j*delta, (i-1)*delta) or isInsideCar((j+1)*delta, i*delta) or isInsideCar((j-1)*delta, i*delta))):
            if (not isInsideCar(j*delta, i*delta, 1) and (isInsideCar(j*delta, (i-1)*delta) or isInsideCar((j+1)*delta, i*delta) or isInsideCar((j-1)*delta, i*delta))) or (((j*delta)-d-L/2)**2+(i*delta-h)**2==R**2 and i*delta > h):
            #if (not isInsideCar(j*delta, i*delta, 1) and (isInsideCar(j*delta, (i-1)*delta) or isInsideCar((j+1)*delta, i*delta) or isInsideCar((j-1)*delta, i*delta))) or (((j*delta)-d-L/2)**2+(i*delta-h)**2==R**2 and i*delta > h):
                #print("i e j que entraram: ", i, j)
                inContour[i][j] = 1
                p_contour.append(p[i][j])
                x_contour.append(j*delta)
    #print("INCONTOUR UPPER", inContour)
    x_contour = np.array(x_contour)
    p_contour = np.array(p_contour)    
    sorted_indices = np.argsort(x_contour)
    sorted_x_contour = x_contour[sorted_indices]
    sorted_p_contour = p_contour[sorted_indices]
    return sorted_x_contour, sorted_p_contour

def mdf_p_lower_contour(delta, p):
    n_rows, n_cols = p.shape
    inContour = np.zeros((n_rows, n_cols), dtype=bool)    
    p_contour = []
    x_contour = []   
    for i in range(n_rows):
        for j in range(n_cols):
            #if (not isInsideCar(j*delta, i*delta, 1) and isInsideCar(j*delta, (i+1)*delta, 1)):
            #if ((j*delta >= d and j*delta <= d+L) and i*delta==h) or (not isInsideCar(j*delta, i*delta, 1) and isInsideCar(j*delta, (i+1)*delta)):
            if not isInsideCar(j*delta, i*delta, 1) and isInsideCar(j*delta, (i+1)*delta):
                inContour[i][j] = 1
                p_contour.append(p[i][j])
                x_contour.append(j*delta)
    #print("INCONTOUR LOWER", inContour)
    x_contour = np.array(x_contour)
    p_contour = np.array(p_contour)  
    return x_contour, p_contour

def calculate_lift_force(x_upper, p_upper, x_lower, p_lower):
    f = 0
    unique_x_positions, unique_indices = np.unique(x_upper, return_index=True)
    x_contour = x_upper[np.sort(unique_indices)]
    p_contour = np.array([np.mean(p_upper[x_upper == x]) for x in x_contour])

    if len(x_contour) > len(x_lower):
        indices = np.isin(x_contour, x_lower)
        x_filtered = x_contour[indices]
        p_filtered = p_contour[indices]
        assert len(x_filtered) == len(x_lower)
        for j in range(1, len(x_filtered)):
            if x_filtered[j] < d+R:
                dtheta = np.arccos((d+R-x_filtered[j])/R) - np.arccos((d+R-x_filtered[j-1])/R)
                theta = np.arccos((d+R-x_filtered[j])/R)
            elif x_filtered[j] > d+R and x_filtered[j-1] < d+R:
                dtheta = np.pi - np.arccos((x_filtered[j]-d-R)/R) - np.arccos((d+R-x_filtered[j-1])/R)
                theta = np.arccos((x_filtered[j]-d-R)/R)
            elif  x_filtered[j-1] > d+R:
                dtheta = np.arccos((x_filtered[j-1]-d-R)/R) - np.arccos((x_filtered[j]-d-R)/R)
                theta = np.arccos((x_filtered[j]-d-R)/R)
            f -= p_filtered[j]*dtheta*np.sin(theta)                #aproximando fndl = frdthetasin(theta) na integral
        f *= R
        for j in range(1, len(x_lower)):
            f += p_lower[j]*(x_lower[j]-x_lower[j-1])
    else:
        indices = np.isin(x_lower, x_contour)
        x_filtered = x_lower[indices]
        p_filtered = p_lower[indices]
        assert len(x_contour) == len(x_filtered)
        for j in range(1, len(x_filtered)):
            if x_contour[j] < d+R:
                dtheta = np.arccos((d+R-x_contour[j])/R) - np.arccos((d+R-x_contour[j-1])/R)
                theta = np.arccos((d+R-x_contour[j])/R)
            elif x_contour[j] > d+R and x_contour[j-1] < d+R:
                dtheta = np.pi - np.arccos((x_contour[j]-d-R)/R) - np.arccos((d+R-x_contour[j-1])/R)
                theta = np.arccos((x_contour[j]-d-R)/R)
            elif  x_contour[j-1] > d+R:
                dtheta = np.arccos((x_contour[j-1]-d-R)/R) - np.arccos((x_contour[j]-d-R)/R)
                theta = np.arccos((x_contour[j]-d-R)/R)
            f -= p_contour[j]*dtheta*np.sin(theta)                 #aproximando fndl = frdthetasin(theta) na integral
        f *= R
        for j in range(1, len(x_filtered)):
            f += p_filtered[j]*(x_filtered[j]-x_filtered[j-1])
    return f

#def calculate_lift_force(x_contour, p_contour, x_lower, p_lower):
#    f = 0
#    print("len(x_contour)", len(x_contour))
#    print("len(x_lower)", len(x_lower))
#    if len(x_contour) > len(x_lower):
#        indices = np.isin(x_contour, x_lower)
#        x_filtered = x_contour[indices]
#        p_filtered = p_contour[indices]
#        assert len(x_filtered) == len(x_lower)
#        for j in range(1, len(x_filtered)):
#            if x_filtered[j] < d+R:
#                dtheta = np.arccos((d+R-x_filtered[j])/R) - np.arccos((d+R-x_filtered[j-1])/R)
#                theta = np.arccos((d+R-x_filtered[j])/R)
#            elif x_filtered[j] > d+R and x_filtered[j-1] < d+R:
#                dtheta = np.pi - np.arccos((x_filtered[j]-d-R)/R) - np.arccos((d+R-x_filtered[j-1])/R)
#                theta = np.arccos((x_filtered[j]-d-R)/R)
#            elif  x_filtered[j-1] > d+R:
#                dtheta = np.arccos((x_filtered[j-1]-d-R)/R) - np.arccos((x_filtered[j]-d-R)/R)
#                theta = np.arccos((x_filtered[j]-d-R)/R)
#            f -= p_filtered[j]*dtheta*np.sin(theta)                #aproximando fndl = frdthetasin(theta) na integral
#        f *= R
#        for j in range(1, len(x_lower)):
#            f += p_lower[j]*(x_lower[j]-x_lower[j-1])
#    else:
#        indices = np.isin(x_lower, x_contour)
#        x_filtered = x_lower[indices]
#        p_filtered = p_lower[indices]
#        assert len(x_contour) == len(x_filtered)
#        for j in range(1, len(x_filtered)):
#            if x_contour[j] < d+R:
#                dtheta = np.arccos((d+R-x_contour[j])/R) - np.arccos((d+R-x_contour[j-1])/R)
#                theta = np.arccos((d+R-x_contour[j])/R)
#            elif x_contour[j] > d+R and x_contour[j-1] < d+R:
#                dtheta = np.pi - np.arccos((x_contour[j]-d-R)/R) - np.arccos((d+R-x_contour[j-1])/R)
#                theta = np.arccos((x_contour[j]-d-R)/R)
#            elif  x_contour[j-1] > d+R:
#                dtheta = np.arccos((x_contour[j-1]-d-R)/R) - np.arccos((x_contour[j]-d-R)/R)
#                theta = np.arccos((x_contour[j]-d-R)/R)
#            f -= p_contour[j]*dtheta*np.sin(theta)                 #aproximando fndl = frdthetasin(theta) na integral
#        f *= R
#        for j in range(1, len(x_filtered)):
#            f += p_filtered[j]*(x_filtered[j]-x_filtered[j-1])
#    return f

def mdf_temp(delta, lamb, u, v):
    n_lines = int(2*L/delta)+1        
    n_columns = n_lines
    temp = np.full((n_lines, n_columns), 0, dtype=float)        #precisa estar entre 0 e 1
    max_error = 1       #inicialmente temp e temp velho são iguais
    while(max_error>epsilon):
        temp_old = np.copy(temp) 
        for i in range(temp.shape[0]): # linhas, iteração de baixo para cima
            for j in range(temp.shape[1]): # colunas, iteração da esquerda para a direita
                if (i == temp.shape[0]-1 and j == temp.shape[1]-1): # canto superior direito 
                    temp[i][j] = (temp[i][j-1] + temp[i-1][j])/2      
                elif (i == 0 and j ==temp.shape[1]-1): # canto inferior direito
                    temp[i][j] = (temp[i][j-1] + temp[i+1][j])/2
                elif (j == 0): # Parede esquerda do dominio, incluindo os cantos superior e inferior
                    temp[i][j] = 20
                elif (j == temp.shape[1]-1): # Parede direita do dominio
                    if v[i][j] >= 0:
                        temp[i][j] = k*(2*temp[i][j-1] + temp[i-1][j] + temp[i+1][j])/(delta**2) + rho*cp*u[i][j]*(temp[i-1][j])/delta
                        temp[i][j] /= 4*k/(delta**2) + rho*cp*u[i][j]/delta
                    else:
                        temp[i][j] = k*(2*temp[i][j-1] + temp[i-1][j] + temp[i+1][j])/(delta**2) - rho*cp*u[i][j]*(temp[i+1][j])/delta
                        temp[i][j] /= 4*k/(delta**2) - rho*cp*u[i][j]/delta
                elif (i == temp.shape[0]-1): # Teto do dominio
                    temp[i][j] = k*(temp[i][j+1] + temp[i][j-1] + 2*temp[i-1][j])/(delta**2) + rho*cp*u[i][j]*(temp[i][j-1])/delta
                    temp[i][j] /= 4*k/(delta**2) + rho*cp*u[i][j]/delta
                elif (i == 0): # Chao do domínio
                    temp[i][j] = k*(temp[i][j+1] + temp[i][j-1] + 2*temp[i+1][j])/(delta**2) + rho*cp*u[i][j]*(temp[i][j-1])/delta
                    temp[i][j] /= 4*k/(delta**2) + rho*cp*u[i][j]/delta
                elif isInsideCar(j*delta, i*delta, 1) and i*delta <= np.tan(np.pi/3) * j*delta + h - 3*np.tan(np.pi/3): # Motor do carro
                    temp[i][j] = 80
                elif isInsideCar(j*delta, i*delta, 1):      #dentro do carro
                    temp[i][j] = 25
                elif isInsideCar((j+1)*delta, i*delta, 1) and isInsideCar(j*delta, (i-1)*delta, 1):       #contorno a direita e abaixo
                    # n_center = n_columns - j -2                    #numero de deltas (inteiros) dentro do circulo na direção horizontal até o centro
                    # y_center = i*delta - h                  #altura do pto do contorno (e do pto a sua esquerda) ate o centro
                    # bij = (delta - (np.sqrt(R**2-(y_center)**2)-n_center*delta))/delta
                    # assert bij<1
                    # n_center = n_columns - j -1 #numero de deltas (inteiros) dentro do circulo na direção horizontal até o centro
                    # y_center = (i-1)*delta - h   #altura do pto do contorno ate o centro
                    # aij = (delta - (np.sqrt(R**2-(n_center*delta)**2)-y_center))/delta  #pitagoras
                    # assert aij<1
                    # temp[i][j] = k*(2*temp[i+1][j]/(aij+1) + 2*temp[i-1][j]/(aij*(aij+1)) + temp[i][j+1]/(bij*(bij+1)) + temp[i][j-1]/(bij + 1))/(delta**2) + rho*cp*u[i][j]*(temp[i][j-1] + temp[i-1][j])/delta
                    # temp[i][j] /= 2*k/(aij*delta**2) + 2*k/(bij*delta**2) + 2*rho*cp*u[i][j]/delta
                    pass
                elif isInsideCar((j-1)*delta, i*delta, 1) and isInsideCar(j*delta, (i-1)*delta, 1):       #contorno a esquerda e abaixo
                    # n_center = n_columns - j -2                    #numero de deltas (inteiros) dentro do circulo na direção horizontal até o centro
                    # y_center = i*delta - h                  #altura do pto do contorno (e do pto a sua esquerda) ate o centro
                    # bij = (delta - (np.sqrt(R**2-(y_center)**2)-n_center*delta))/delta
                    # assert bij<1
                    # n_center = n_columns - j -1 #numero de deltas (inteiros) dentro do circulo na direção horizontal até o centro
                    # y_center = (i-1)*delta - h   #altura do pto do contorno ate o centro
                    # aij = (delta - (np.sqrt(R**2-(n_center*delta)**2)-y_center))/delta  #pitagoras
                    # assert aij<1
                    # temp[i][j] = k*(2*temp[i+1][j]/(aij+1) + 2*temp[i-1][j]/(aij*(aij+1)) + temp[i][j-1]/(bij*(bij+1)) + temp[i][j+1]/(bij + 1))/(delta**2) - rho*cp*u[i][j]*(temp[i+1][j] - temp[i][j-1])/delta
                    # temp[i][j] /= 2*k/(aij*delta**2) + 2*k/(bij*delta**2)
                    pass
                elif isInsideCar(j*delta, (i-1)*delta, 1): #Contorno está abaixo
                    # n_center = n_columns/2 - j - 1 #numero de deltas (inteiros) dentro do circulo na direção horizontal até o centro
                    # y_center = (i-1)*delta - h   #altura do pto do contorno ate o centro
                    # aij = (delta - (np.sqrt(R**2-(n_center*delta)**2)-y_center))/delta  #pitagoras
                    # print(aij)
                    # assert aij<1
                    # if v[i][j]>=0:
                    #     temp[i][j] = k*(2*temp[i+1][j]/(aij+1) + 2*temp[i-1][j]/aij*(aij+1) + temp[i][j+1] + temp[i][j-1])/(delta**2) + rho*cp*u[i][j]*(temp[i][j-1] + temp[i-1][j])/delta
                    #     temp[i][j] /= 2*k/(aij*delta**2) + 2*k/(delta**2) + 2*rho*cp*u[i][j]/delta
                    # else:
                    #     temp[i][j] = k*(2*temp[i+1][j]/(aij+1) + 2*temp[i-1][j]/aij*(aij+1) + temp[i][j+1] + temp[i][j-1])/(delta**2) - rho*cp*u[i][j]*(temp[i+1][j] - temp[i][j-1])/delta
                    #     temp[i][j] /= 2*k/(aij*delta**2) + 2*k/(delta**2)
                    pass
                elif isInsideCar((j+1)*delta, i*delta, 1): #Contorno está a direita
                    n_center = n_columns/2 - j -2                    #numero de deltas (inteiros) dentro do circulo na direção horizontal até o centro
                    y_center = i*delta - h                  #altura do pto do contorno (e do pto a sua esquerda) ate o centro
                    bij = (delta - (np.sqrt(R**2-(y_center)**2)-n_center*delta))/delta  #pitagoras
                    assert bij<1    
                    temp[i][j] = k*(2*temp[i][j-1]/(bij+1) + 2*temp[i][j+1]/bij*(bij+1) + temp[i+1][j] + temp[i-1][j])/(delta**2) + rho*cp*u[i][j]*(temp[i][j-1] + temp[i-1][j])/delta
                    temp[i][j] /= 2*k/(bij*delta**2) + 2*k/(delta**2) + 2*rho*cp*u[i][j]/delta
                elif isInsideCar((j-1)*delta, i*delta, 1): #Contorno está a esquerda
                    # n_center = n_columns/2 - j -2                    #numero de deltas (inteiros) dentro do circulo na direção horizontal até o centro
                    # y_center = i*delta - h                  #altura do pto do contorno (e do pto a sua esquerda) ate o centro
                    # bij = (delta - (np.sqrt(R**2-(y_center)**2)-n_center*delta))/delta  #pitagoras
                    # assert bij<1    
                    # temp[i][j] = k*(2*temp[i][j+1]/(bij+1) + 2*temp[i][j-1]/bij*(bij+1) + temp[i+1][j] + temp[i-1][j])/(delta**2) - rho*cp*u[i][j]*(temp[i+1][j] - temp[i][j-1])/delta
                    # temp[i][j] /= 2*k/(bij*delta**2) + 2*k/(delta**2)
                    pass
                elif isInsideCar(j*delta, (i+1)*delta, 1): #Contorno está acima
                    # n_center = n_columns/2 - j - 1 #numero de deltas (inteiros) dentro do circulo na direção horizontal até o centro
                    # y_center = (i-1)*delta - h   #altura do pto do contorno ate o centro
                    # aij = (delta - (np.sqrt(R**2-(n_center*delta)**2)-y_center))/delta  #pitagoras
                    # print(aij)
                    # assert aij<1
                    # if v[i][j]>=0:
                    #     temp[i][j] = k*(2*temp[i-1][j]/(aij+1) + 2*temp[i+1][j]/aij*(aij+1) + temp[i][j+1] + temp[i][j-1])/(delta**2) + rho*cp*u[i][j]*(temp[i][j-1] + temp[i-1][j])/delta
                    #     temp[i][j] /= 2*k/(aij*delta**2) + 2*k/(delta**2) + 2*rho*cp*u[i][j]/delta
                    # else:
                    #     temp[i][j] = k*(2*temp[i+1][j]/(aij+1) + 2*temp[i-1][j]/aij*(aij+1) + temp[i][j+1] + temp[i][j-1])/(delta**2) - rho*cp*u[i][j]*(temp[i+1][j] - temp[i][j-1])/delta
                    #     temp[i][j] /= 2*k/(aij*delta**2) + 2*k/(delta**2)
                    pass
                else: # Parte central do dominio
                    if v[i][j] >= 0:
                        temp[i][j] = k*(temp[i][j+1] + temp[i][j-1] + temp[i-1][j] + temp[i+1][j])/(delta**2) + rho*cp*u[i][j]*(temp[i-1][j] + temp[i][j-1])/delta
                        temp[i][j] /= 4*k/(delta**2) + 2*rho*cp*u[i][j]/delta
                    # else:
                    #     temp[i][j] = k*(temp[i][j+1] + temp[i][j-1] + temp[i-1][j] + temp[i+1][j])/(delta**2) - rho*cp*u[i][j]*(temp[i+1][j] - temp[i][j-1])/delta
                    #     temp[i][j] /= 4*k/(delta**2)
                    pass
                temp[i][j] = lamb*temp[i][j]+(1-lamb)*temp_old[i][j]
        max_error = np.max(np.abs((temp-temp_old)))
        print("T", temp)
        print("T old", temp_old)
        print("Max error ", max_error)
    return temp

def main():
    # print(mdf_psi(0.4))
    #mdf_psi(0.045)
    psi = mdf_psi(L/8, 1.85)         #NÃO usar divisor de 0.15 ou 0.2 com os filtros implementados 
    u, v = mdf_vel(L/8, psi)   #cálculo de F_lift fica errado
    temp = mdf_temp(L/8, 1.15, u, v)
    # p = mdf_pressure(u, v)
    print(temp)
    plot_temp(L/8, temp)
    # print("U:", u)
    # print("V", v)
    plot_velfield(L/8, u, v)
    #plot_pressure(0.1, p)
    # x_upper, p_upper = mdf_p_upper_contour(0.06, p)
    # x_lower, p_lower = mdf_p_lower_contour(0.06, p)
    #assert len(x_lower) ==len(x_upper)
    # plot_p_upper_contour(x_upper, p_upper)
    # plot_p_lower_contour(x_lower, p_lower)
    #print(x_contour)s
    #print(p_contour)
    #print("X UPPER ", x_upper)
    #print("X LOWER", x_lower)
    # F = calculate_lift_force(x_upper, p_upper, x_lower, p_lower)
    # print("F", F)
    plt.show()
    #print(F)
    #plot_p_contour(x_contour, p_contour)


if __name__ == "__main__":
    import time
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))
