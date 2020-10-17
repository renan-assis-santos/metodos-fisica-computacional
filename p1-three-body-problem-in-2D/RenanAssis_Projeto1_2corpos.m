## Definicao de constantes
G = (6.67408e-11)*((6.685e-12)^3)/((5.03e-31)*(3.17e-8)^2); # Constante gravitacional [UA^3/(Msolar*ano^2)]
M_sol = 1; # Massa do sol [Msolar]
M_terra = 3.003e-6; # Massa da terra [Msolar]
d_sol_terra = 1; # Distancia sol-terra [UA]
v_terra = 109040*(6.685e-9)/(1.14155e-4); # Velocidade de translacao da terra [UA/ano]
dt = 1/365.25; # Passo temporal [ano]
t = 0:dt:1; # Range de tempo [ano]
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

## Metodo Runge-Kutta 2a ordem
for i=2:Nt
  
  ## Calc. k1s Sol
  k1_xS = pxS(i-1)*dt/M_sol;
  k1_yS = pyS(i-1)*dt/M_sol;
  k1_pxS = ((-G*M_terra*M_sol*(xS(i-1)-xT(i-1)))/(((xT(i-1)-xS(i-1))^2+(yT(i-1)-yS(i-1))^2)^(3/2)))*dt;
  k1_pyS = ((-G*M_terra*M_sol*(yS(i-1)-yT(i-1)))/(((xT(i-1)-xS(i-1))^2+(yT(i-1)-yS(i-1))^2)^(3/2)))*dt;
  
  ## Calc. k1s Terra
  k1_xT = pxT(i-1)*dt/M_terra;
  k1_yT = pyT(i-1)*dt/M_terra;
  k1_pxT = ((-G*M_terra*M_sol*(xT(i-1)-xS(i-1)))/(((xS(i-1)-xT(i-1))^2+(yS(i-1)-yT(i-1))^2)^(3/2)))*dt;
  k1_pyT = ((-G*M_terra*M_sol*(yT(i-1)-yS(i-1)))/(((xS(i-1)-xT(i-1))^2+(yS(i-1)-yT(i-1))^2)^(3/2)))*dt;
  
  ## Calc. valores 1/2 passo Sol
  xS_meio = xS(i-1) + k1_xS/2;
  yS_meio = yS(i-1) + k1_yS/2;
  pxS_meio = pxS(i-1) + k1_pxS/2;
  pyS_meio = pyS(i-1) + k1_pyS/2;
  
  ## Calc. valores 1/2 passo Terra
  xT_meio = xT(i-1) + k1_xT/2;
  yT_meio = yT(i-1) + k1_yT/2;
  pxT_meio = pxT(i-1) + k1_pxT/2;
  pyT_meio = pyT(i-1) + k1_pyT/2;
  
  ## Calc. k2s Sol
  k2_xS = pxS_meio*dt/M_sol;
  k2_yS = pyS_meio*dt/M_sol;
  k2_pxS = ((-G*M_terra*M_sol*(xS_meio-xT_meio))/(((xT_meio-xS_meio)^2+(yT_meio-yS_meio)^2)^(3/2)))*dt;
  k2_pyS = ((-G*M_terra*M_sol*(yS_meio-yT_meio))/(((xT_meio-xS_meio)^2+(yT_meio-yS_meio)^2)^(3/2)))*dt;
  
  ## Calc. k2s Terra
  k2_xT = pxT_meio*dt/M_terra;
  k2_yT = pyT_meio*dt/M_terra;
  k2_pxT = ((-G*M_terra*M_sol*(xT_meio-xS_meio))/(((xS_meio-xT_meio)^2+(yS_meio-yT_meio)^2)^(3/2)))*dt;
  k2_pyT = ((-G*M_terra*M_sol*(yT_meio-yS_meio))/(((xS_meio-xT_meio)^2+(yS_meio-yT_meio)^2)^(3/2)))*dt;
  
  ## Avancar no tempo Sol
  xS(i) = xS(i-1) + k2_xS;
  yS(i) = yS(i-1) + k2_yS;
  pxS(i) = pxS(i-1) + k2_pxS;
  pyS(i) = pyS(i-1) + k2_pyS;

  ## Avancar no tempo Terra
  xT(i) = xT(i-1) + k2_xT;
  yT(i) = yT(i-1) + k2_yT;
  pxT(i) = pxT(i-1) + k2_pxT;
  pyT(i) = pyT(i-1) + k2_pyT;
  
endfor

## Calc. energia cinetica e potencial dos astros 
K_Sol = (pxS.^2+pyS.^2)/(2*M_sol);
K_Terra = (pxT.^2+pyT.^2)/(2*M_terra);
V_Sol = 0;
V_Terra = (-G*M_terra*M_sol)./(((xS-xT).^2+(yS-yT).^2).^(1/2));;

## Energias totais
K = K_Sol + K_Terra;
V = V_Sol + V_Terra;
E_tot = K + V;

## Vetor distancia da Terra em relacao ao Sol
rST = sqrt((xT-xS).^2+(yT-yS).^2);

## Figura energias do sistema vs tempo
figure1 = figure(1);
plot(t,K,'LineWidth',1.5,t,V,'LineWidth',1.5,t,E_tot,'LineWidth',1.5);
set(gca,'FontSize',18);
xlabel('tempo(anos)','FontSize',18);
ylabel('Energia (Msolar*UA^{2}*Ano^{-2})','FontSize',18);
ylim([-1.5e-4 1.5e-4]);
title('Energia total do sistema','FontSize',18);
h = legend('Energia Cinética','Energia Potencial','Energia Total','location','north');
set(h,'FontSize',15);

## Figura distancia entre a Terra e o Sol vs tempo
figure2 = figure(2);
plot(t,rST,'LineWidth',1.5); grid on;
set(gca,'FontSize',18);
xlabel('tempo(anos)','FontSize',18);
ylabel('|r(t)|(UA)','FontSize',18);
title('Distância entre a Terra e o Sol ao longo de 1 ano','FontSize',18);

## Figura animação
figure3  = figure(3);
count = 1; # contador para a animação
plot(xS,yS,'r-','LineWidth',1.5,...
     xT,yT,'b-','LineWidth',1.5);hold on; # Plota a órbita completa dos astros 

## Loop para fazer a animação
for i = 1:Nt
  p1 = plot(xS(i),yS(i),'bo','LineWidth',2,...
            'MarkerEdgeColor','r',...
            'MarkerFaceColor',[1.0 1.0 0.0],...
            'MarkerSize',36,...
            xT(i),yT(i),'bo','LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor',[0.0 1.0 1.0],...
            'MarkerSize',16); grid on;
  axis([-1.5 1.5 -1.5 1.5]);
  set(gca,'FontSize',18);
  xlabel('x(UA)','FontSize',18);
  ylabel('y(UA)','FontSize',18);
  title('Órbita da Terra','FontSize',18);
  pause(1e-6);
  M(count) = getframe(figure3);
  delete(p1);
  count = count + 1;
endfor
