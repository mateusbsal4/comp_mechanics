import numpy as np
import matplotlib.pyplot as plt
from numba import njit


##Parameter definitions (SI)
M = 1783
a = 1.22
b = 1.5
Ic = 4000
e = 0.75
L = 0.5
A = 0.06
f = 0.35
f_e = 35 #2100/60
m_e = 20
w_e = 2*(np.pi)*f_e   
r = 0.045
Fn = m_e*(w_e**2)*r


@njit
def d1(t, V):
  if t>2:
    return 0
  w = 2*np.pi*(V/L)  
  return A*(1-np.cos(w*t))
@njit
def d2(t, V):
  if t>2:
    return 0
  w = 2*np.pi*(V/L)  
  return A*(1+np.cos(w*t))
@njit
def d1_dot(t, V):
  if t>2:
    return 0
  w = 2*np.pi*(V/L)  
  return A*w*np.sin(w*t)
@njit
def d2_dot(t, V): 
  if t>2:
    return 0
  w = 2*np.pi*(V/L)  
  return -A*w*np.sin(w*t)

@njit
def F(t, Y, k1, k2, c1, c2, V):
  K = np.zeros(4)
  x = Y[0]
  theta = Y[1]
  u = Y[2]
  v = Y[3]
  # Derivada primeira de x
  K[0] = u 
  # Derivada primeira de theta
  K[1] = v 
  # Derivada segunda de x
  K[2] = (-(k1+k2)*x+(k1*a-k2*b)*theta-(c1+c2)*u+(c1*a-c2*b)*v+k1*d1(t, V)+k2*d2(t, V)+c1*d1_dot(t, V)+c2*d2_dot(t, V)+Fn*np.sin(w_e*t))/M     
  # Derivada segunda de theta
  K[3] = ((k1*a-k2*b)*x-(k1*a**2+k2*b**2)*theta+(c1*a-c2*b)*u-(c1*a**2+c2*b**2)*v-k1*a*d1(t, V)+k2*b*d2(t, V)-c1*a*d1_dot(t, V)+c2*b*d2_dot(t, V)-Fn*(e*np.sin(w_e*t)+f*np.cos(w_e*t)))/Ic
  return K

@njit
def rk4_solver(Y0, t0, tf, h, k, c, V): 
    iterations  = int(np.floor((tf-t0)/h)) # Discretiza o tempo em intervalos de tamanho h
    Y = Y0
    for i in range(iterations): # Itera cada intervalo de tempo
        t_i = t0 + h*i
        # Aplicacao do método de Runge-Kutta de quarta ordem
        K1 = F(t_i, Y, k, k, c, c, V)
        K2 = F(t_i + 0.5*h, Y + 0.5*h*K1, k, k, c, c, V)
        K3 = F(t_i + 0.5*h, Y + 0.5*h*K2, k, k, c, c, V)
        K4 = F(t_i+ h, Y + h*K3, k, k, c, c, V)
        Y += (h/6)*(K1+2*K2+2*K3+K4)
    return Y


def rk4_plotter(h, V, k1, k2, c1, c2, scale_x, scale_theta, scale_xdot, scale_thetadot, scale_xddot, scale_thetaddot, Y0, t0 = 0, tf = 4, n =1000):
  t_sampled = np.linspace(t0, tf, n)
  Y = np.zeros((n,4))
  Y_dd = np.zeros((n, 2))
  for i in range(n):
    Y[i] = rk4_solver(Y0, t0, t_sampled[i], h, k1, c1, V) # Roda Runge-kutta para cada ponto 
    # Salva a solução obtida em cada uma das variáveis
    x = Y[i][0]
    theta = Y[i][1]
    u = Y[i][2]
    v = Y[i][3]
    # Obtem derivada segunda
    Y_dd[i] = np.array([((-(k1+k2)*x+(k1*a-k2*b)*theta-(c1+c2)*u+(c1*a-c2*b)*v+k1*d1(t_sampled[i], V)+k2*d2(t_sampled[i], V)+c1*d1_dot(t_sampled[i], V)+c2*d2_dot(t_sampled[i], V)+Fn*np.sin(w_e*t_sampled[i]))/M), (((k1*a-k2*b)*x-(k1*a**2+k2*b**2)*theta+(c1*a-c2*b)*u-(c1*a**2+c2*b**2)*v-k1*a*d1(t_sampled[i], V)+k2*b*d2(t_sampled[i], V)-c1*a*d1_dot(t_sampled[i], V)+c2*b*d2_dot(t_sampled[i], V)-Fn*(e*np.sin(w_e*t_sampled[i])+f*np.cos(w_e*t_sampled[i])))/Ic)])
  if (np.isnan(Y).any()):
     print("Overflow no vetor resultante")
  Z = np.transpose(Y)
  Y_dd = np.transpose(Y_dd)
  fig, ax = plt.subplots(figsize=(10,10))
  # Plota os pontos para essa iteração do runge-kutta
  ax.plot(t_sampled, scale_x*Z[0], label=r'$x(t) \times$' + "{}".format(scale_x))  
  ax.plot(t_sampled, scale_xdot*Z[2], label=r'$\dot{x}(t) \times$' + "{}".format(scale_xdot))  
  ax.plot(t_sampled, scale_xddot*Y_dd[0], label=r'$\ddot{x}(t) \times$' + "{}".format(scale_xddot))  
  ax.set_xlabel('t')  
  ax.set_ylabel('f(t)')  
  ax.set_title("Evolução temporal em x (h = {}, V = {}, k1=k2={}, c1=c2={})".format(h, V*3.6, k1, c1))  
  ax.legend() 
  
  fig, ax = plt.subplots(figsize=(10,10))
  ax.plot(t_sampled, scale_theta*Z[1], label=r'$\theta(t) \times$' + "{}".format(scale_theta))  
  ax.plot(t_sampled, scale_thetadot*Z[3], label=r'$\dot{\theta}(t) \times$' + "{}".format(scale_thetadot))
  ax.plot(t_sampled, scale_thetaddot*Y_dd[1], label=r'$\ddot{\theta}(t) \times$' + "{}".format(scale_thetaddot))   
  ax.set_xlabel('t')  
  ax.set_ylabel('f(t)')  
  ax.set_title("Evolução temporal em theta (h = {}, V = {}, k1=k2={}, c1=c2={})".format(h, V*3.6, k1, c1))  
  ax.legend() 

  plt.show()

def main():
    # Variáveis para cada alternativa
    h_s = 0.002    #small
    h_m = 0.004    #medium
    h_l = 0.02     #large
    V1 = 50/3.6
    V2 = 30/3.6
    V3 = 70/3.6
    k_a = 2.8 * 10**7
    c_a = 3 * 10**4
    k_b1  = 2*10**4
    c_b1  = 10**3
    k_b2  = 5*10**6
    c_b2  = 5 * 10**4
    k_b3  = 7*10**8
    c_b3  = 2.5*10**5


    print("Item a")
    print("Mostrando resultados para o passo pequeno h = 0.002s")
    rk4_plotter(h_s, V1, k_a, k_a, c_a, c_a, 10, 10, 0.1, 0.1, 0.001, 0.001, np.array([0, 0.09, 0, 0]))
    print("Mostrando resultados para o passo medio h = 0.004s")
    rk4_plotter(h_m, V1, k_a, k_a, c_a, c_a, 10, 10, 0.1, 0.1, 0.001, 0.001, np.array([0, 0.09, 0, 0]))
    print("Mostrando resultados para o passo grande h = 0.02s")
    rk4_plotter(h_l, V1, k_a, k_a, c_a, c_a, 10, 10, 0.1, 0.1, 0.001, 0.001, np.array([0, 0.09, 0, 0]))

    print("\n")

    print("Item b")
    print("i)")
    print("Mostrando resultados para V = 30km/h") 
    rk4_plotter(h_s, V2, k_a, k_a, c_a, c_a, 100, 100, 1, 1, 0.001, 0.001, np.array([0, 0.09, 0, 0]))
    print("Mostrando resultados para V = 70km/h") 
    rk4_plotter(h_s, V3, k_a, k_a, c_a, c_a, 10, 10, 1, 0.1, 0.001, 0.001, np.array([0, 0.09, 0, 0]))  
    
    print("ii)")
    print("Mostrando resultados para c = 1000 kg/s")
    rk4_plotter(h_s, V1, k_a, k_a, c_b1, c_b1, 10, 10, 0.1, 0.1, 0.001, 0.001, np.array([0, 0.09, 0, 0]))
    print("Mostrando resultados para c = 50000 kg/s")
    rk4_plotter(h_s, V1, k_a, k_a, c_b2, c_b2, 100, 100, 1, 0.1, 0.001, 0.001, np.array([0, 0.09, 0, 0]))
    print("Mostrando resultados para c = 250000 kg/s") 
    rk4_plotter(h_s, V1, k_a, k_a, c_b3, c_b3, 10, 10, 0.1, 0.1, 0.001, 0.001, np.array([0, 0.09, 0, 0]))  
    
    print("iii)")
    print("Mostrando resultados para k = 2000 N/m")
    rk4_plotter(h_s, V1, k_b1, k_b1, c_a, c_a, 10, 10, 0.1, 0.1, 0.001, 0.001, np.array([0, 0.09, 0, 0]))
    print("Mostrando resultados para k = 5000000 N/m")
    rk4_plotter(h_s, V1, k_b2, k_b2, c_a, c_a, 100, 10, 1, 1, 0.01, 0.001, np.array([0, 0.09, 0, 0]))
    print("Mostrando resultados para k = 700000000 N/m")
    rk4_plotter(h_s, V1, k_b3, k_b3, c_a, c_a, 10, 10, 0.1, 0.1, 0.001, 0.001, np.array([0, 0.09, 0, 0]))


if __name__ == "__main__":
    main()
