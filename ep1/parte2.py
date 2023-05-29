import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

L = 3
h = 0.15
d = L/2
H = 2*L

def plot(delta, graph):
    fig = plt.figure()
    x = np.linspace(0, L + 2*d, int(2*L/delta)+1)
    y = np.linspace(0, H, int(2*L/delta) + 1)
    x_1, y_1 = np.meshgrid(x, y)

    plt.pcolor(x_1, y_1, graph, cmap = 'jet')
    plt.colorbar()
    plt.show()

def isInsideCar(x, y):
    return (y <= np.sqrt((L/2)**2 - (x-d-L/2)**2) + h) and y >= h

def mdf_psi(delta):
    psi = np.zeros((int(2*L/delta) + 1, int(2*L/(2*delta)) + 1))
    for i in range(psi.shape[0]): # linhas, iteração de baixo para cima
        for j in range(psi.shape[1]): # colunas, iteração da esquerda para a direita
            #divisão por regioes especificas
            if (i == psi.shape[0]-1 and j == 0): # canto superior
                psi[i][j] = 12
            elif (i == 0 and j == 0): # canto inferior
                psi[i][j] = 13
            elif (j == 0): # Parede (esquerda) do dominio
                psi[i][j] = 9
            elif (i == psi.shape[0]-1): # Teto do dominio
                if (j == psi.shape[1]-1): # Caso esteja no MEIO do dominio (CUIDADO)
                    pass
                psi[i][j] = 1
            elif (i == 0): # Chao do domínio
                psi[i][j] = 0
            elif isInsideCar(j*delta, i*delta):
                psi[i][j] = -1
                if not isInsideCar(j*delta, (i+1)*delta):
                    print("Contorno cima")
                if not isInsideCar(j*delta, (i-1)*delta):
                    print("Contorno baixo")
                if not isInsideCar((j-1)*delta, i*delta):
                    print("Contorno lateral")
            else: # Parte central do dominio
                if (j == psi.shape[1]-1): # Caso esteja no MEIO do dominio (CUIDADO)
                    pass
                psi[i][j] = 3
    return np.concatenate((psi, np.flip(psi,1)), axis=1)

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
            elif isInsideCar(j*delta, i*delta):
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
    print(mdf_psi(0.4))
    print("T:")
    print(mdf_temp(0.4))
    plot(0.01, mdf_temp(0.01))


if __name__ == "__main__":
    import time
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))