clear; clc; close all;  

delta_x = 0.005;
e = 0.03;
l = 0.7;
h = 0.4;
d = h/2;
a = 0.55;
b = 0.15;
c = l/2;
theta1 = atan2(h, c)*(180/pi)
theta2 = atan2(h, a-c)*(180/pi)
x = (delta_x*cosd(theta1))
y = (delta_x*cosd(theta1))
z =  e+c- delta_x*cosd(theta1)

w = delta_x*sind(theta1)
u = delta_x*sind(theta1)
v = h -delta_x*sind(theta1)

o = round(1000*(sqrt(h^2+c^2)))/1000;
r= round(1000*(sqrt(h^2+(a-c)^2)))/1000;

p = linspace(delta_x, o-delta_x, round((o-2*delta_x)/delta_x) +1); 
q = linspace(delta_x, r-delta_x, round((r-2*delta_x)/delta_x) +1);
 
#q = (h-sind(theta2)*q)';


coor = [(0:delta_x:e)' zeros(e/delta_x+1, 1); #coord 1 e 2
        (e + delta_x:delta_x:e+a)' zeros(a/delta_x, 1); # 2 e 3
        (a+e+delta_x:delta_x:l+e)' zeros(b/delta_x, 1); #3 e 4
        (l+e+delta_x: delta_x: l+2*e)' zeros(e/delta_x, 1); #4 e 5
        (l+e)*ones(d/delta_x, 1) (delta_x: delta_x: d)'; #coord entre nos 4 e 8
        (l+e+delta_x:delta_x: l+2*e)' d*ones(e/delta_x, 1); #coord 8 e 9
         (l+e)*ones(d/delta_x, 1) (d+delta_x: delta_x: h)'; # coord 8 e 13
         (l+e)*ones(e/delta_x, 1) (h+delta_x: delta_x: h+e)'; # coord 13 e 16
          (l+e-delta_x: -delta_x: a+e)' h*ones(b/delta_x, 1); # coord 13 e 11
          (a+e-delta_x: -delta_x: c+e)' h*ones(((0.2)/delta_x), 1); # coord 12 e 11
        (c+e)*ones(e/delta_x, 1) (h+delta_x: delta_x: h+e)'; #entre 11 e
        (c+e-delta_x: -delta_x: e)' h*ones(c/delta_x, 1); #entre 11 e 10
        e*ones(e/delta_x, 1) (h+delta_x: delta_x: h+e)'; #entre 10 e 14
        e*ones(d/delta_x, 1) (h-delta_x: -delta_x: d)'; #entre 10 e 7
        (e-delta_x: -delta_x: 0)' d*ones(e/delta_x, 1); #7 e 6
        e*ones(d/delta_x -1, 1) (d-delta_x: -delta_x: delta_x)'; #7 e 2
        (e+cosd(theta1)*p)' (sind(theta1)*p)'; # entre 2 e 11
        (c+e+cosd(theta2)*q)' (h-sind(theta2)*q)'; #entre 11 e 3
        #(c+ e+ (round(delta_x*cosd(theta2)*1000)/1000): (round(delta_x*cosd(theta2)*1000)/1000): e+a - (round(delta_x*cosd(theta2)*1000)/1000))' (h-(round(delta_x*sind(theta2)*1000)/1000): -(round(delta_x*sind(theta2)*1000)/1000): (round(delta_x*sind(theta2)*1000)/1000))'; 
        (a+e)*ones(h/delta_x-1, 1) (delta_x: delta_x: h-delta_x)';  #entre 3 e 12
        0 0.5155; #nó 55
        0.150 0.5155; #no 56
        0.580 0.5155; #no 57
        0.760 0.5155]; #No 58
        

        
                  
con = [(1: 1: e/delta_x)' (2: 1: e/delta_x+1)';      #1-4            
        (e/delta_x+a/delta_x+b/delta_x+1: 1: (2*e)/delta_x+a/delta_x+b/delta_x)' (e/delta_x+a/delta_x+b/delta_x+2: 1: (2*e)/delta_x+a/delta_x+b/delta_x+1)';  #10-13                
        ((2*e)/delta_x+a/delta_x+b/delta_x+d/delta_x+1: 1: (3*e)/delta_x+a/delta_x+b/delta_x+d/delta_x)' ((2*e)/delta_x+a/delta_x+b/delta_x+d/delta_x+2: 1: (3*e)/delta_x+a/delta_x+b/delta_x+d/delta_x+1)';    #16-19      
        ((3*e)/delta_x+a/delta_x+b/delta_x+(2*d)/delta_x+1: 1: (4*e)/delta_x+a/delta_x+b/delta_x+(2*d)/delta_x)' ((3*e)/delta_x+a/delta_x+b/delta_x+(2*d)/delta_x+2: 1: (4*e)/delta_x+a/delta_x+b/delta_x+(2*d)/delta_x + 1)'; #nós 22-25
       ((4*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x+1: 1: (5*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x)' ((4*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x+2: 1: (5*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x+1)'; #nó 31 - 34         
       ((5*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x + c/delta_x + 1: 1: (6*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x + c/delta_x)' ((5*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x + c/delta_x + 2: 1: (6*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x + c/delta_x + 1)'; #nós 37-40          
       ((6*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x + c/delta_x + d/delta_x+1: 1: (7*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x + c/delta_x + d/delta_x)' ((6*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x + c/delta_x + d/delta_x+2: 1: (7*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x + c/delta_x + d/delta_x+1)'; #Nós 43-46  
        #acabou azul   
        (e/delta_x+1: 1: e/delta_x+a/delta_x)' (e/delta_x+2: 1: e/delta_x+a/delta_x+1)'; #el 8 - discretizado do no 2 ate 3 ( 4- 7, 4 5, 5 6, 6 7) 
       (e/delta_x+a/delta_x+1: 1: e/delta_x+a/delta_x+b/delta_x)' (e/delta_x+a/delta_x+2: 1: e/delta_x+a/delta_x+b/delta_x+1)'; #el 9 - discretizado do no 3 ate 4  (7 -10) 
       e/delta_x+a/delta_x+b/delta_x+1 (2*e)/delta_x+a/delta_x+b/delta_x+2; #el 11 - nó 10 ao 14 (10, 14)
       ((2*e)/delta_x+a/delta_x+b/delta_x+2: 1: (2*e)/delta_x+a/delta_x+b/delta_x+d/delta_x)'   ((2*e)/delta_x+a/delta_x+b/delta_x+3: 1: (2*e)/delta_x+a/delta_x+b/delta_x+d/delta_x+1)';#el 11 - nó 14 ao 16 
       (2*e)/delta_x+a/delta_x+b/delta_x+d/delta_x+1 (3*e)/delta_x+a/delta_x+b/delta_x+d/delta_x+2; # nó 16 ao 20
       ((3*e)/delta_x+a/delta_x+b/delta_x+d/delta_x+2: 1: (3*e)/delta_x+a/delta_x+b/delta_x+(2*d)/delta_x)' ((3*e)/delta_x+a/delta_x+b/delta_x+d/delta_x+3: 1: (3*e)/delta_x+a/delta_x+b/delta_x+(2*d)/delta_x+1)';#nó 20 ao 22 (20 21, 21 22)
        (3*e)/delta_x+a/delta_x+b/delta_x+(2*d)/delta_x+1 (4*e)/delta_x+a/delta_x+b/delta_x+(2*d)/delta_x + 2; #22 a 26
        ((4*e)/delta_x+a/delta_x+b/delta_x+(2*d)/delta_x + 2:1:(4*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x)' ((4*e)/delta_x+a/delta_x+b/delta_x+(2*d)/delta_x + 3: 1: (4*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + 1)'; #26 ao 28
        ((4*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + 1: 1: (4*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x)' ((4*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + 2: 1: (4*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x+1)'; #nó 28 ao 31      
       (4*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x+1 (5*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x+2; #nós 31 e 35
       ((5*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x+2: 1: (5*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x + c/delta_x)' ((5*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x+3:1: (5*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x + c/delta_x + 1)'; #nós 35-37  
      (5*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x + c/delta_x + 1 (6*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x + c/delta_x + 2; #nós 37 e 41
       ((6*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x + c/delta_x + 2: 1: (6*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x + c/delta_x + d/delta_x)' ((6*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x + c/delta_x + 3: 1: (6*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x + c/delta_x + d/delta_x+1)';      #nóa 41 - 43 
       (6*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x + c/delta_x + d/delta_x+1 (7*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x + c/delta_x + d/delta_x+2; #nós 43 e 47
       ((7*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x + c/delta_x + d/delta_x+2: 1: (7*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x + c/delta_x + (2*d)/delta_x-1)' ((7*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x + c/delta_x + d/delta_x+3: 1: (7*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x + c/delta_x + (2*d)/delta_x)'; #nos  47 48 
       (7*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x + c/delta_x + (2*d)/delta_x e/delta_x+1; #nós 48 e 4   
          #acabou laranja
       e/delta_x+1 (7*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x + c/delta_x + (2*d)/delta_x+1; #nós 4 e 49
       ((7*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x + c/delta_x + (2*d)/delta_x+1: 1: (7*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x + c/delta_x + (2*d)/delta_x+1 + round(o/delta_x) -3)' ((7*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x + c/delta_x + (2*d)/delta_x +2: 1: (7*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x + c/delta_x + (2*d)/delta_x+1 + round(o/delta_x) -3 + 1)';#nós 49 e 50 
       ((7*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x + c/delta_x + (2*d)/delta_x+1 + round(o/delta_x) -3 + 1) ((4*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x+1); #nós 50 a 31   
        ((4*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x+1) ((7*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x + c/delta_x + (2*d)/delta_x+1 + round(o/delta_x) -3 + 2); #nós 31 a 51
        (((7*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x + c/delta_x + (2*d)/delta_x+1 + round(o/delta_x) -3+2):1: ((7*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x + c/delta_x + (2*d)/delta_x+1 + round(o/delta_x) -3+2 + round(r/delta_x) - 3))' (((7*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x + c/delta_x + (2*d)/delta_x+1 + round(o/delta_x) -3+3): 1: ((7*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x + c/delta_x + (2*d)/delta_x+1 + round(o/delta_x) -3+2 + round(r/delta_x) - 2))' #nós 51 52 
       ((7*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x + c/delta_x + (2*d)/delta_x+1 + round(o/delta_x) -3+2 + round(r/delta_x) - 2) (e/delta_x+a/delta_x+1); #nós 52 7
        (e/delta_x+a/delta_x+1) ((7*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x + c/delta_x + (2*d)/delta_x+1 + round(o/delta_x) -3+2 + round(r/delta_x) - 1); #nós 7 e 53
        ((7*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x + c/delta_x + (2*d)/delta_x+1 + round(o/delta_x) -3+2 + round(r/delta_x) - 1: 1: (7*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x + c/delta_x + (2*d)/delta_x+1 + round(o/delta_x) -3+2 + round(r/delta_x) - 1 + h/delta_x -3)' ((7*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x + c/delta_x + (2*d)/delta_x+1 + round(o/delta_x) -3+2 + round(r/delta_x): 1: (7*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x + c/delta_x + (2*d)/delta_x+1 + round(o/delta_x) -3+2 + round(r/delta_x) - 1 + h/delta_x -2)'; #53 a 54
        ((7*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x + c/delta_x + (2*d)/delta_x+1 + round(o/delta_x) -3+2 + round(r/delta_x) - 1 + h/delta_x -2) ((4*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + 1);# nós 54 e 28 
        #acabou verde
        ((7*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x + c/delta_x + (2*d)/delta_x+1 + round(o/delta_x) -3+2 + round(r/delta_x) - 1 + h/delta_x -1) ((6*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x + c/delta_x + 1); # 55 a 40
        ((6*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x + c/delta_x + 1) ((7*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x + c/delta_x + (2*d)/delta_x+1 + round(o/delta_x) -3+2 + round(r/delta_x) - 1 + h/delta_x);#40 a 56
        ((7*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x + c/delta_x + (2*d)/delta_x+1 + round(o/delta_x) -3+2 + round(r/delta_x) - 1 + h/delta_x -1) ((7*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x + c/delta_x + (2*d)/delta_x+1 + round(o/delta_x) -3+2 + round(r/delta_x) - 1 + h/delta_x); # 55 a 56
       ((7*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x + c/delta_x + (2*d)/delta_x+1 + round(o/delta_x) -3+2 + round(r/delta_x) - 1 + h/delta_x+1) ((7*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x + c/delta_x + (2*d)/delta_x+1 + round(o/delta_x) -3+2 + round(r/delta_x) - 1 + h/delta_x+2); #57 a 58
       ((7*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x + c/delta_x + (2*d)/delta_x+1 + round(o/delta_x) -3+2 + round(r/delta_x) - 1 + h/delta_x+2) ((4*e)/delta_x+a/delta_x+b/delta_x+(2*d)/delta_x + 1); #58 a 25
       ((7*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x + c/delta_x + (2*d)/delta_x+1 + round(o/delta_x) -3+2 + round(r/delta_x) - 1 + h/delta_x+1) ((4*e)/delta_x+a/delta_x+b/delta_x+(2*d)/delta_x + 1); #57 a 25       
        #acabou preto (A1)       
       ((7*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x + c/delta_x + (2*d)/delta_x+1 + round(o/delta_x) -3+2 + round(r/delta_x) - 1 + h/delta_x) ((5*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x+1); #56 a 34
        ((5*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x+1) ((7*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x + c/delta_x + (2*d)/delta_x+1 + round(o/delta_x) -3+2 + round(r/delta_x) - 1 + h/delta_x+1)]; #34 a 57
       #acabou molas
 
  
                  
         
plot_struct(coor, con);        
        
Nod = size(coor, 1)
Nel = size(con, 1)
Ngdl = 3* Nod

%%% Prop dos elementos

#contablização de quantos elementos de cada seçao dada
n_a3 = 7*(e/delta_x)
n_a2 = (2*l+2*h)/delta_x
n_a1 = round(o/delta_x) + round(r/delta_x) + h/delta_x
n_a4 = 6
n_molas = 2

sum = n_a1+n_a2+n_a3+n_a4+n_molas

#areas
data.A = [1e-6*(400-15^2)*ones(n_a3, 1);
          1e-6*(15*25-12.6*22.6)*ones(n_a2, 1);
          1e-6*(pi/4)*(25.4^2 - 22.4^2)*ones(n_a1, 1);
          1e-6*(150)*ones(n_a4, 1);
          -1*ones(n_molas, 1)];  %n se aplica para as molas
          
          
area = data.A
#modulo Young
data.E = 210e9*ones(Nel, 1);
#data.E
#densidade
data.rho = 7860*ones(Nel, 1);
#data.rho(n_a4+1:n_molas, 1) = -1;

#momentos de inercia

data.I = [1e-12*(20^4/12 - 15^4/12)*ones(n_a3, 1); #A3 - azuis
          1e-12*((15*25^3)/12-(12.6*22.6^3)/12)*ones(n_a2, 1); #A2 - laranja
          1e-12*(pi/4)*(25.4^4 - 22.4^4)*ones(n_a1, 1);#A1 - verde 
          1e-12*(5*30^3)/12*ones(n_a4, 1); #A4 preto - ja trelica
          -1*ones(n_molas, 1)];  %n se aplica para as molas
          
#inercia = data.I
#%para as molas
K = 200e3;
m = 0.5;


data.L = zeros(Nel, 1);
data.Q = zeros(Nel, 1);


  for e = 1:Nel
    x1 = coor(con(e, 1),1);
    x2 = coor(con(e, 2),1);
    
    y1 = coor(con(e, 1), 2);
    y2 = coor(con(e, 2), 2);
    
    data.L(e) = sqrt((x2-x1)^2 + (y2-y1)^2);
    data.Q(e) = atan2(y2-y1, x2-x1)*(180/pi);   
  endfor

#data.L
#data.Q
  
  
k0_trelica = [1 0 0 -1 0 0;
           0 0 0 0 0 0;
           0 0 0 0 0 0;
          -1 0 0 1 0 0;
           0 0 0 0 0 0;
           0 0 0 0 0 0];
m0_trelica = [2 0 0 1 0 0;
           0 0 0 0 0 0;
           0 0 0 0 0 0;
           1 0 0 2 0 0;
           0 0 0 0 0 0;
           0 0 0 0 0 0];
           

Kg = zeros(Ngdl, Ngdl);
Mg = zeros(Ngdl, Ngdl);

for e = 1:Nel 
  c = cosd(data.Q(e)); s = sind(data.Q(e));
  if e<=n_a1+n_a2+n_a3   #se for portico
    l = data.L(e);
    k0_viga = [0 0 0 0 0 0;
               0 12 6*l 0 -12 6*l;
               0 6*l 4*l^2 0 -6*l 2*l^2;
               0 0 0 0 0 0;
               0 -12 -6*l 0 12 -6*l;
               0 6*l 2*l^2 0 -6*l 4*l^2];                 
    m0_viga = [0 0 0 0 0 0;
               0 156 22*l 0 54 -13*l;
               0 22*l 4*l^2 0 13*l -3*l^2;
               0 0 0 0 0 0;
               0 54 13*l 0 156 -22*l;
               0 -13*l -3*l^2 0 -22*l 4*l^2];    
   ke = (data.E(e)*data.A(e)/data.L(e))*k0_trelica + (data.E(e)*data.I(e)/((data.L(e))^3))*k0_viga; 
   me = (data.rho(e)*data.A(e)*data.L(e)/6)*m0_trelica + (data.rho(e)*data.A(e)*data.L(e)/420)*m0_viga;
   T = [c -s 0 0 0 0; #u1
        s c 0 0 0 0; #w1 
        0 0 1 0 0 0; #phi1
        0 0 0 c -s 0; #u2         #matriz de rotacao p/elemento de portico - tirado da literatura
        0 0 0 -s c 0; #w2
        0 0 0 0 0 1]; #phi2
     # u1w1p1u2w2p2  
  else #se for trelica
    if e <= n_a4 #se for A1
      me = (data.rho(e)*data.A(e)*data.L(e)/6)*m0_trelica;
      ke = (data.E(e)*data.A(e)/data.L(e))*k0_trelica;      
    else #se for molas
      me = m*m0_trelica;
      ke = K*k0_trelica;
    endif
    T = [c s 0 0 0 0; #u1
      -  s c 0 0 0 0; #w1 
         0 0 1 0 0 0; #phi1
         0 0 0 c s 0; #u2
         0 0 0 -s c 0; #w2
         0 0 0 0 0 1]; #phi2
     # u1w1p1u2w2p2  
  endif

  %%Levando em conta inclinacoes
  ke = T'*ke*T;
  me = T'*me*T;
 
  nod1 = con(e, 1);
  nod2 = con(e, 2); 
  
  Kg(3*nod1-2:3*nod1, 3*nod1-2:3*nod1) = Kg(3*nod1-2:3*nod1, 3*nod1-2:3*nod1) + ke(1:3, 1:3);
  Kg(3*nod1-2:3*nod1, 3*nod2-2:3*nod2) = Kg(3*nod1-2:3*nod1, 3*nod2-2:3*nod2) + ke(1:3, 4:6);  
  Kg(3*nod2-2:3*nod2, 3*nod1-2:3*nod1) = Kg(3*nod2-2:3*nod2, 3*nod1-2:3*nod1) + ke(4:6, 1:3);
  Kg(3*nod2-2:3*nod2, 3*nod2-2:3*nod2) = Kg(3*nod2-2:3*nod2, 3*nod2-2:3*nod2) + ke(4:6, 4:6);  

  Mg(3*nod1-2:3*nod1, 3*nod1-2:3*nod1) = Mg(3*nod1-2:3*nod1, 3*nod1-2:3*nod1) + me(1:3, 1:3);
  Mg(3*nod1-2:3*nod1, 3*nod2-2:3*nod2) = Mg(3*nod1-2:3*nod1, 3*nod2-2:3*nod2) + me(1:3, 4:6);  
  Mg(3*nod2-2:3*nod2, 3*nod1-2:3*nod1) = Mg(3*nod2-2:3*nod2, 3*nod1-2:3*nod1) + me(4:6, 1:3);
  Mg(3*nod2-2:3*nod2, 3*nod2-2:3*nod2) = Mg(3*nod2-2:3*nod2, 3*nod2-2:3*nod2) + me(4:6, 4:6);   
endfor    

Kgm = Kg;

node1 = 1;  #v = 0
node2 = (2*e)/delta_x+a/delta_x+b/delta_x+1; #v = 0

node3 = (3*e)/delta_x+a/delta_x+b/delta_x+d/delta_x+1; #u = 0       #43 e 16
node4 = (7*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x + c/delta_x + d/delta_x+1; #u = 0

node5 = ((7*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x + c/delta_x + (2*d)/delta_x+1 + round(o/delta_x) -3+2 + round(r/delta_x) - 1 + h/delta_x -1);
node6 = ((7*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x + c/delta_x + (2*d)/delta_x+1 + round(o/delta_x) -3+2 + round(r/delta_x) - 1 + h/delta_x);
node7 = ((7*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x + c/delta_x + (2*d)/delta_x+1 + round(o/delta_x) -3+2 + round(r/delta_x) - 1 + h/delta_x +1);    #phi = 0 (55-58)
node8 = ((7*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x + c/delta_x + (2*d)/delta_x+1 + round(o/delta_x) -3+2 + round(r/delta_x) - 1 + h/delta_x +2);
node9 = (5*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x+1;              #Phi = 0
node10 = (6*e)/delta_x+a/delta_x+(2*b)/delta_x+(2*d)/delta_x + (a-c)/delta_x + c/delta_x + 1; 
node11 =  (4*e)/delta_x+a/delta_x+b/delta_x+(2*d)/delta_x + 1;


Kgm(:, 3*node3-2) = 0; Kgm(3*node3-2,:) = 0; Kgm(3*node3-2,3*node3-2) = 1;
Kgm(:, 3*node4-2) = 0; Kgm(3*node4-2,:) = 0; Kgm(3*node4-2,3*node4-2) = 1; #fixacoes em u

Kgm(:, 3*node1-1) = 0; Kgm(3*node1-1,:) = 0; Kgm(3*node1-1,3*node1-1) = 1;
Kgm(:, 3*node2-1) = 0; Kgm(3*node2-1,:) = 0; Kgm(3*node2-1,3*node2-1) = 1; #fixacoes em v


Kgm(:, 3*node5) = 0; Kgm(3*node5,:) = 0; Kgm(3*node5,3*node5) = 1;
Kgm(:, 3*node6) = 0; Kgm(3*node6,:) = 0; Kgm(3*node6,3*node6) = 1; #fixacoes em phi
Kgm(:, 3*node7) = 0; Kgm(3*node7,:) = 0; Kgm(3*node7,3*node7) = 1;
Kgm(:, 3*node8) = 0; Kgm(3*node8,:) = 0; Kgm(3*node8,3*node8) = 1; #fixacoes em phi
Kgm(:, 3*node9) = 0; Kgm(3*node9,:) = 0; Kgm(3*node9,3*node9) = 1;
Kgm(:, 3*node10) = 0; Kgm(3*node10,:) = 0; Kgm(3*node10,3*node10) = 1; #fixacoes em phi
Kgm(:, 3*node11) = 0; Kgm(3*node11,:) = 0; Kgm(3*node11,3*node11) = 1; #fixacoes em phi


Mgm(:, 3*node3-2) = 0; Mgm(3*node3-2,:) = 0; Mgm(3*node3-2,3*node3-2) = 1;
Mgm(:, 3*node4-2) = 0; Mgm(3*node4-2,:) = 0; Mgm(3*node4-2,3*node4-2) = 1; #fixacoes em u

Mgm(:, 3*node1-1) = 0; Mgm(3*node1-1,:) = 0; Mgm(3*node1-1,3*node1-1) = 1;
Mgm(:, 3*node2-1) = 0; Mgm(3*node2-1,:) = 0; Mgm(3*node2-1,3*node2-1) = 1; #fixacoes em v


Mgm(:, 3*node5) = 0; Mgm(3*node5,:) = 0; Mgm(3*node5,3*node5) = 1;
Mgm(:, 3*node6) = 0; Mgm(3*node6,:) = 0; Mgm(3*node6,3*node6) = 1; #fixacoes em phi
Mgm(:, 3*node7) = 0; Mgm(3*node7,:) = 0; Mgm(3*node7,3*node7) = 1;
Mgm(:, 3*node8) = 0; Mgm(3*node8,:) = 0; Mgm(3*node8,3*node8) = 1; #fixacoes em phi
Mgm(:, 3*node9) = 0; Mgm(3*node9,:) = 0; Mgm(3*node9,3*node9) = 1;
Mgm(:, 3*node10) = 0; Mgm(3*node10,:) = 0; Mgm(3*node10,3*node10) = 1; #fixacoes em phi
Mgm(:, 3*node11) = 0; Mgm(3*node11,:) = 0; Mgm(3*node11,3*node11) = 1; #fixacoes em phi

A = Mgm/Kgm;
[vec, val] = eig(A);

val = diag(val)(1:6, 1);
  ke = T'*ke*T;
  me = T'*me*T;
endfor

        
