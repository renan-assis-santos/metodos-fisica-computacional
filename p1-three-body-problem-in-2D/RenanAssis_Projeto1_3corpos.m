## Definicao de constantes
G = (6.67408e-11)*((6.685e-12)^3)/((5.03e-31)*(3.17e-8)^2); # Constante gravitacional [UA^3/(Msolar*ano^2)]
M_sol = 1; # Massa do sol [Msolar]
M_terra = 3.003e-6; # Massa da terra [Msolar]
M_jupiter = 9.54791e-4; # Massa de Júpiter [Msolar]
#M_jupiter = 0.1*M_sol;
d_sol_terra = 1; # Distancia sol-terra [UA]
d_sol_jupiter = 5.2; # Distancia sol-jupiter [UA]
v_terra = 109040*(6.685e-9)/(1.14155e-4); # Velocidade de translacao da terra [UA/ano]
v_jupiter = 13.7*(6.685e-9)/(3.17e-8); # Velocidade de translacao de jupiter [UA/ano]
dt = 1/365.25; # Passo temporal [ano]
t = 0:dt:14; # Range de tempo [ano]
Nt = length(t);

## Inicializar vetores do Sol
xS = ones(1,Nt); # posicao X do Sol
yS = ones(1,Nt); # posicao Y do Sol
pxS = ones(1,Nt); # momento em X do Sol
pyS = ones(1,Nt); # momento em Y do Sol

## Inicializar vetores da Terra
xT = ones(1,Nt); # posicao X da Terra
yT = ones(1,Nt); # posicao Y da Terra
pxT = ones(1,Nt); # momento em X da Terra
pyT = ones(1,Nt); # momento em Y da Terra

## Inicializar vetores de Jupiter
xJ = ones(1,Nt); # posicao X de Jupiter
yJ = ones(1,Nt); # posicao Y de Jupiter
pxJ = ones(1,Nt); # momento em X de Jupiter
pyJ = ones(1,Nt); # momento em Y de Jupiter

## Condicoes iniciais do Sol
xS(1) = 0;
yS(1) = 0;
pxS(1) = 0;
pyS(1) = 0;

## Condicoes iniciais da Terra
xT(1) = d_sol_terra;
yT(1) = 0;
pxT(1) = 0;
pyT(1) = v_terra*M_terra;

## Condicoes iniciais de Jupiter
xJ(1) = 0;
yJ(1) = d_sol_jupiter;
pxJ(1) = -v_jupiter*M_jupiter;
pyJ(1) = 0;

## Funcao para calcular a forca do corpo 1 no corpo 2
function [res,res2] = F(G,M1,M2,x1,x2,y1,y2)
  res = (-G*M1*M2*(x2-x1))/(((x2-x1)^2+(y2-y1)^2)^(3/2));
  res2 = (-G*M1*M2*(y2-y1))/(((x2-x1)^2+(y2-y1)^2)^(3/2));
endfunction

## Funcao para calcular ks para a posicao
function [res,res2] = k_d(px,py,M,dt)
  res = px*dt/M;
  res2 = py*dt/M;
endfunction

## Funcao para calcular ks para o momento
function [res,res2] = k_p(F1x,F2x,F1y,F2y,dt)
  res = (F1x+F2x)*dt;
  res2 = (F1y+F2y)*dt;
endfunction

## Metodo Runge-Kutta 2a ordem
for i=1:Nt-1
  
  ## Definir forcas entre os corpos no tempo t
  [F_STx,F_STy] = F(G,M_terra,M_sol,xT(i),xS(i),yT(i),yS(i)); # Forca Terra>>Sol
  [F_TSx,F_TSy] = F(G,M_terra,M_sol,xS(i),xT(i),yS(i),yT(i)); # Forca Sol>>Terra
  [F_SJx,F_SJy] = F(G,M_jupiter,M_sol,xJ(i),xS(i),yJ(i),yS(i)); # Forca Jupiter>>Sol
  [F_JSx,F_JSy] = F(G,M_jupiter,M_sol,xS(i),xJ(i),yS(i),yJ(i)); # Forca Sol>>Jupiter
  [F_TJx,F_TJy] = F(G,M_terra,M_jupiter,xJ(i),xT(i),yJ(i),yT(i)); # Forca Jupiter>>Terra
  [F_JTx,F_JTy] = F(G,M_terra,M_jupiter,xT(i),xJ(i),yT(i),yJ(i)); # Forca Terra>>Jupiter
  
  ## Calc. k1s Sol
  [k1_xS,k1_yS] = k_d(pxS(i),pyS(i),M_sol,dt);
  [k1_pxS,k1_pyS] = k_p(F_STx,F_SJx,F_STy,F_SJy,dt);
  
  ## Calc. k1s Terra
  [k1_xT,k1_yT] = k_d(pxT(i),pyT(i),M_terra,dt);
  [k1_pxT,k1_pyT] = k_p(F_TSx,F_TJx,F_TSy,F_TJy,dt);
  
  ## Calc. k1s Jupiter
  [k1_xJ,k1_yJ] = k_d(pxJ(i),pyJ(i),M_jupiter,dt);
  [k1_pxJ,k1_pyJ] = k_p(F_JSx,F_JTx,F_JSy,F_JTy,dt);
  
  ## Calc. valores 1/2 passo Sol
  xS_meio = xS(i) + k1_xS/2;
  yS_meio = yS(i) + k1_yS/2;
  pxS_meio = pxS(i) + k1_pxS/2;
  pyS_meio = pyS(i) + k1_pyS/2;
  
  ## Calc. valores 1/2 passo Terra
  xT_meio = xT(i) + k1_xT/2;
  yT_meio = yT(i) + k1_yT/2;
  pxT_meio = pxT(i) + k1_pxT/2;
  pyT_meio = pyT(i) + k1_pyT/2;
  
  ## Calc. valores 1/2 passo Jupiter
  xJ_meio = xJ(i) + k1_xJ/2;
  yJ_meio = yJ(i) + k1_yJ/2;
  pxJ_meio = pxJ(i) + k1_pxJ/2;
  pyJ_meio = pyJ(i) + k1_pyJ/2;
  
  ## Definir forcas entre os corpos no tempo t+dt/2
  [F_STx,F_STy] = F(G,M_terra,M_sol,xT_meio,xS_meio,yT_meio,yS_meio); # Forca Terra>>Sol
  [F_TSx,F_TSy] = F(G,M_terra,M_sol,xS_meio,xT_meio,yS_meio,yT_meio); # Forca Sol>>Terra
  [F_SJx,F_SJy] = F(G,M_jupiter,M_sol,xJ_meio,xS_meio,yJ_meio,yS_meio); # Forca Jupiter>>Sol
  [F_JSx,F_JSy] = F(G,M_jupiter,M_sol,xS_meio,xJ_meio,yS_meio,yJ_meio); # Forca Sol>>Jupiter
  [F_TJx,F_TJy] = F(G,M_terra,M_jupiter,xJ_meio,xT_meio,yJ_meio,yT_meio); # Forca Jupiter>>Terra
  [F_JTx,F_JTy] = F(G,M_terra,M_jupiter,xT_meio,xJ_meio,yT_meio,yJ_meio); # Forca Terra>>Jupiter
  
  ## Calc. k2s Sol
  [k2_xS,k2_yS] = k_d(pxS_meio,pyS_meio,M_sol,dt);
  [k2_pxS,k2_pyS] = k_p(F_STx,F_SJx,F_STy,F_SJy,dt);
  
  ## Calc. k2s Terra
  [k2_xT,k2_yT] = k_d(pxT_meio,pyT_meio,M_terra,dt);
  [k2_pxT,k2_pyT] = k_p(F_TSx,F_TJx,F_TSy,F_TJy,dt);
  
  ## Calc. k2s Jupiter
  [k2_xJ,k2_yJ] = k_d(pxJ_meio,pyJ_meio,M_jupiter,dt);
  [k2_pxJ,k2_pyJ] = k_p(F_JSx,F_JTx,F_JSy,F_JTy,dt);
  
  ## Avancar no tempo Sol
  xS(i+1) = xS(i) + k2_xS;
  yS(i+1) = yS(i) + k2_yS;
  pxS(i+1) = pxS(i) + k2_pxS;
  pyS(i+1) = pyS(i) + k2_pyS;  
  
  ## Avancar no tempo Terra
  xT(i+1) = xT(i) + k2_xT;
  yT(i+1) = yT(i) + k2_yT;
  pxT(i+1) = pxT(i) + k2_pxT;
  pyT(i+1) = pyT(i) + k2_pyT;
  
  ## Avancar no tempo Jupiter
  xJ(i+1) = xJ(i) + k2_xJ;
  yJ(i+1) = yJ(i) + k2_yJ;
  pxJ(i+1) = pxJ(i) + k2_pxJ;
  pyJ(i+1) = pyJ(i) + k2_pyJ;
   
endfor

## Calc. energia cinetica e potencial dos astros 
K_Sol = (pxS.^2+pyS.^2)/(2*M_sol);
K_Terra = (pxT.^2+pyT.^2)/(2*M_terra);
K_Jupiter = (pxJ.^2+pyJ.^2)/(2*M_jupiter);
V_Sol = 0;
V_Terra = (-G*M_terra*M_sol)./(((xS-xT).^2+(yS-yT).^2).^(1/2));
V_Jupiter = (-G*M_jupiter*M_terra)./(((xT-xJ).^2+(yT-yJ).^2).^(1/2)) + ...
            (-G*M_jupiter*M_sol)./(((xS-xJ).^2+(yS-yJ).^2).^(1/2));

## Energias totais
K = K_Sol + K_Terra + K_Jupiter;
V = V_Sol + V_Terra + V_Jupiter;
E_tot = K + V;

## Vetor distancia dos astros em relacao ao sol
rST = sqrt((xT-xS).^2+(yT-yS).^2);
rSJ = sqrt((xJ-xS).^2+(yJ-yS).^2);

## Figura energias do sistema vs tempo
figure1 = figure(1);
plot(t,K,'LineWidth',1.5,t,V,'LineWidth',1.5,t,E_tot,'LineWidth',1.5);
set(gca,'FontSize',18);
xlabel('tempo(anos)','FontSize',18);
ylabel('Energia (Msolar*UA^{2}*Ano^{-2})','FontSize',18);
ylim([-8e-3 7e-3]);
title('Energia total do sistema','FontSize',18);
h = legend('Energia Cinética','Energia Potencial','Energia Total','location','north');
set(h,'FontSize',15);

## Figura distancia entre Jupiter e o Sol vs tempo
figure2 = figure(2);
plot(t,rSJ,'LineWidth',1.5); grid on;
set(gca,'FontSize',18);
xlabel('tempo(anos)','FontSize',18);
ylabel('|r(t)|(UA)','FontSize',18);
title('Distância de Júpiter em relação ao Sol ao longo de 14 anos','FontSize',18);

## Figura para a animacao
figure3 = figure(3);
plot(xS,yS,'r-','LineWidth',1.5,...
     xT,yT,'b-','LineWidth',1.5,...
     xJ,yJ,'k-','LineWidth',1.5);hold on; # Plota a orbita completa dos astros
count = 1; # contador para a animacao

## Loop para fazer a animação
for i = 1:Nt
  p1 = plot(xS(i),yS(i),'bo','LineWidth',2,...
            'MarkerEdgeColor','r',...
            'MarkerFaceColor',[1.0 1.0 0.0],...
            'MarkerSize',27,...
            xT(i),yT(i),'bo','LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor',[0.0 1.0 1.0],...
            'MarkerSize',10,...
            xJ(i),yJ(i),'bo','LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor',[0.5 0.5 0.5],...
            'MarkerSize',16); hold on; grid on;
  axis([-7 7 -7 7]);
  set(gca,'FontSize',18);
  xlabel('x(UA)','FontSize',18);
  ylabel('y(UA)','FontSize',18);
  title('Órbita da Terra e de Júpiter em torno do Sol','FontSize',18);
  pause(1e-6);
  M(count) = getframe(figure2);
  delete(p1);
  count = count + 1;
endfor
