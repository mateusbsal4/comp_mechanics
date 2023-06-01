import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math

L = 3
h = 0.15
d = L/2
H = 2*L
lamb = 1.85
V = 100/3.6
R = L/2

def plot(delta, graph):
    fig = plt.figure()
    x = np.linspace(0, L + 2*d, int(2*L/delta))
    y = np.linspace(0, H, int(2*L/delta))
    print(x)
    x_1, y_1 = np.meshgrid(x, y)
    # print("len x", len(x))
    # print("len y", len(y))
    # print("x ", x_1)
    # print("y ", y_1)
    ax = fig.add_subplot(111, projection='3d')
    
    # Plot the surface
    ax.plot_surface(x_1, y_1, graph)

     #Set labels and title
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('Surface plot')
    #plt.pcolor(x_1, y_1, graph, cmap = 'jet')
    #plt.colorbar()
    plt.show()

def isInsideCar(x, y):
    return (y <= np.sqrt((L/2)**2 - (x-d-L/2)**2) + h) and y >= h
    

def mdf_psi(delta):
    n_lines = int(2*L/delta)        
    assert n_lines%2 != 0         #PRECISA ser ímpar
    n_columns = math.floor(L/delta)+1   
    psi = np.zeros((n_lines, n_columns))
    print("SHAPE PSI", psi.shape[0], psi.shape[1])
    for i in range(psi.shape[0]): # linhas, iteração de baixo para cima
        for j in range(psi.shape[1]): # colunas, iteração da esquerda para a direita
            psi_velho = psi[i][j]
            #divisão por regioes especificas
            if (i == psi.shape[0]-1 and j == 0): # canto superior
                psi[i][j] = 1
            elif (i == 0 and j == 0): # canto inferior
                psi[i][j] = 1
            elif (j == 0): # Parede (esquerda) do dominio
                psi[i][j] = 0.25*(2*psi[i][1]+psi[i+1][0]+psi[i-1][0])
            elif (i == psi.shape[0]-1): # Teto do dominio
                if (j == psi.shape[1]-1): # Caso esteja no MEIO do dominio (CUIDADO)
                    psi[i][j] = 0.25*(psi[i][j-1]+psi[i][j-1]+2*psi[i-1][j]+V*delta)
                else:
                    psi[i][j] = 0.25*(psi[i][j+1]+psi[i][j-1]+2*psi[i-1][j]+V*delta)
            elif (i == 0): # Chao do domínio
                psi[i][j] = 0
            elif isInsideCar(j*delta, i*delta):
                psi[i][j] = 0
            elif isInsideCar((j+1)*delta, i*delta): #Contorno está a direita
                n_center = n_columns - j -1                    #numero de deltas (inteiros) dentro do circulo na direção horizontal até o centro
                y_center = i*delta - h                  #altura do pto do contorno (e do pto a sua esquerda) ate o centro
                bij = delta - (np.sqrt(R**2-(y_center)**2)-n_center*delta)  #pitagoras
                psi[i][j] = (bij/(2*(bij+1)))*(psi[i+1][j]+psi[i-1][j]+(2/(bij+1))*psi[i][j-1])
            elif isInsideCar(j*delta, (i-1)*delta): #Contorno está abaixo
                n_center = n_columns - j -1 #numero de deltas (inteiros) dentro do circulo na direção horizontal até o centro
                y_center = (i-1)*delta - h   #altura do pto do contorno ate o centro
                aij = delta - (np.sqrt(R**2-(n_center*delta)**2)-y_center)  #pitagoras
                if (j == psi.shape[1]-1): # Caso esteja no MEIO do dominio (CUIDADO)
                    psi[i][j] = (aij/(2*(aij+1)))*(psi[i][j-1]+psi[i][j-1]+(2/aij)*psi[i+1][j])
                else:
                    psi[i][j] = (aij/(2*(aij+1)))*(psi[i][j+1]+psi[i][j-1]+(2/aij)*psi[i+1][j])
            elif isInsideCar((j+1)*delta, i*delta) and isInsideCar(j*delta, (i-1)*delta):       #contorno a direita e embaixo
                n_center = n_columns - j -1                    #numero de deltas (inteiros) dentro do circulo na direção horizontal até o centro
                y_center = i*delta - h                  #altura do pto do contorno (e do pto a sua esquerda) ate o centro
                bij = delta - (np.sqrt(R**2-(y_center)**2)-n_center*delta)
                n_center = n_columns - j -1 #numero de deltas (inteiros) dentro do circulo na direção horizontal até o centro
                y_center = (i-1)*delta - h   #altura do pto do contorno ate o centro
                aij = delta - (np.sqrt(R**2-(n_center*delta)**2)-y_center)  #pitagoras
                psi[i][j] = ((aij*bij)/(aij+bij))*((psi[i][j-1]/(bij+1))+(psi[i+1][j])/(aij+1))
            else: # Parte central do dominio
                if (j == psi.shape[1]-1): # Caso esteja no MEIO do dominio (CUIDADO)
                    psi[i][j] = 0.25*(psi[i+1][j]+psi[i-1][j]+psi[i][j-1]+psi[i][j-1])
                else:
                    psi[i][j] = 0.25*(psi[i+1][j]+psi[i-1][j]+psi[i][j+1]+psi[i][j-1])
            psi[i][j] = lamb*psi[i][j]+(1-lamb)*psi_velho
    return np.concatenate((psi, np.flip(psi, 1)[:, 1:]), axis=1)

def mdf_temp(delta):
    temp = np.zeros((int(2*L/delta) + 1,int(2*L/delta) + 1))
    for i in range(temp.shape[0]): # linhas, iteração de baixo para cima
        for j in range(temp.shape[1]): # colunas, iteração da esquerda para a direita
            #divisão por regioes especificas
            if (i == temp.shape[0]-1 and j == temp.shape[1]-1): # canto superior direita
                temp[i][j] = 21
            elif (i == temp.shape[0]-1 and j == 0): # canto superior esquerda
                temp[i][j] = 12
            elif (i == 0 and j == temp.shape[1]-1): # canto inferior direita
                temp[i][j] = 31
            elif (i == 0 and j == 0): # canto inferior esquerda
                temp[i][j] = 13
            elif (i == temp.shape[0]-1): # Teto do dominio
                temp[i][j] = 1
            elif (i == 0):
                temp[i][j] = 0 # Chao do domínio
            elif (j == temp.shape[1]-1): # Parede (direita) do dominio
                temp[i][j] = 2
            elif (j == 0): # Parede (esquerda) do dominio
                temp[i][j] = 9
            elif isInsideCar(j*delta, i*delta) and i*delta <= np.tan(np.pi/3) * j*delta + h - 3*np.tan(np.pi/3): # Motor do carro
                temp[i][j] = 80
                if not isInsideCar(j*delta, (i+1)*delta):
                    print("Contorno cima")
                if not isInsideCar(j*delta, (i-1)*delta):
                    print("Contorno baixo")
                if not isInsideCar((j+1)*delta, i*delta):
                    print("Contorno direita")
            elif isInsideCar(j*delta, i*delta): # Dentro do carro
                temp[i][j] = -1
                if not isInsideCar(j*delta, (i+1)*delta):
                    print("Contorno cima")
                if not isInsideCar(j*delta, (i-1)*delta):
                    print("Contorno baixo")
                if not isInsideCar((j+1)*delta, i*delta):
                    print("Contorno direita")
                if not isInsideCar((j-1)*delta, i*delta):
                    print("Contorno esquerda")
            else: # Parte central do dominio
                temp[i][j] = 3
    return temp

def main():
    print("Psi:")
    # print(mdf_psi(0.4))
    # print("T:")
    # print(mdf_temp(0.045))
    plot(0.045, mdf_psi(0.045))


if __name__ == "__main__":
    import time
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))
