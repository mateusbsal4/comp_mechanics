clear; clc; close all;  

coor = [ 0 0;
        30 0;
        580 0;
        730 0;
        760 0;
        0 200;
        30 200;
        730 200;
        760 200;
        30 400;
        380 400;
        580 400;
        730 400;
        30 430;
        380 430;
        730 430;
        0 515.5;
        180 515.5;
        580 515.5;
        760 515.5];
        
con  = [1 2;    %el 1
        4 5;    %el 2
        6 7;   %el 3 - A3
        8 9;     %...
        10 14;
        11 15;
        13 16;      % el 7 - ultimo A3
        2 3;
        3 4;
        2 7;
        4 8;
        7 10;
        8 13;
        10 11;
        11 12;
        12 13;
        2 11;
        3 11;
        3 12;
        14 17;
        14 18;
        17 18;
        16 19;
        16 20;
        19 20;
        15 18;
        15 19];
        
plot_struct(coor, con);        
        
Nod = size(coor, 1);
Nel = size(con, 1);
Ngdl = 2* Nod;

%%% Prop dos elementos
data.A = ones(Nel, 1);
data.A(1:7, 1) = 1e-6*(400-15^2);
data.A(8:16, 1) = 1e-6*(15*25-12.6*22.6);
data.A(17:19, 1) = 1e-6*(pi/4)*(25.4^2 - 22.4^2);
data.A(20:25, 1) = 1e-6*150;
data.A(26:27, 1) = -1;    %n se aplica para as molas


data.E = 210e9*ones(Nel, 1);
data.E(26:27, 1) = -1;

data.rho = 7860*ones(Nel, 1);
data.rho(26:27, 1) = -1;

data.I = ones(Nel, 1);
data.I(1:7, 1) = 1e-12*(20^4/12 - 15^4/12);
data.I(8:16, 1) = 1e-12*((15*25^3)/12-(12.6*22.6^3)/12);
data.I(17:19, 1) = 1e-12*(pi/4)*(25.4^4 - 22.4^4);
data.I(20:25, 1) = 1e-12*(5*30^3)/12;
data.I(26:27, 1) = -1;    %n se aplica para as molas


%para as molas
K = 200e3;
m = 0.5;


%comprimento e rotacao
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
           
%k0_vifa = [
for e = 20:Nel        #para trelicas
  if e<26
    me = (data.rho(e)*data.A(e)*data.L(e)/6)*m0_trelica;
    ke = (data.E(e)*data.A(e)/data.L(e))*k0_trelica;
  else 
    me = m*m0_trelica;
    ke = K*k0_trelica;
  endif
  c = cosd(data.Q(e)); s = sind(data.Q(e));
  T = [c s 0 0 0 0; #u1
      -s c 0 0 0 0; #w1 
       0 0 0 0 0 0; #phi1
       0 0 0 c s 0; #u2
       0 0 0 -s c 0; #w2
       0 0 0 0 0 0]; #phi2
     # u1w1p1u2w2p2
  ke = T'*ke*T;
  me = T'*me*T;
endfor

        