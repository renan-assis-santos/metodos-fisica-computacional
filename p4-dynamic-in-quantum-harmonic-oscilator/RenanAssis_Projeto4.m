##############################################################
# Calcular a evolucao temporal do estado quantico Psi(x,t)
# sujeito a um potencial de oscilador harmonico unidimensional
##############################################################

clear;

function result = psi_n(n,m,w,h,x)
  
  % Solucoes psi_n(x) para a equacao de Schrodinger
  % independente do tempo em termos dos polinomios de Hermite

  cte = m*w/h;
  if n == 0
    result = (cte/pi)^(1/4)*exp(-(cte*(x).^2)/2); 
  elseif n == 1
    result = (cte/pi)^(1/4)*sqrt(2*cte)*x.*exp(-(cte*x.^2)/2); 
  else
    result = (cte/(4*pi))^(1/4)*(2*cte*x.^2-1).*exp(-(cte*x.^2)/2);
  endif
  
endfunction

function res = integracao_13simpson(f,dx,t,inicio,fim)
  
  % Integracao numerica pela regra 1/3 de simpson
  
  soma = f(inicio,t) + f(fim,t);
  for ii = inicio+1:fim-1
    if mod(ii,2) == 0 # Caso ii par
      soma = soma + 4*f(ii,t);
    else # Caso ii impar
      soma = soma + 2*f(ii,t);
    endif
  endfor
  res = (dx/3)*soma;
  
endfunction

## Definicoes iniciais
w = 50;              # Frequencia de oscilacao
T = 2*pi/w;          # Periodo
dt = 1e-4;           # Passo temporal
t = 0:dt:T;          # Vetor tempo
Nt = length(t);      # Numero de elementos do tempo
r = 0.075;           # Constante r
dx = sqrt(dt/(2*r)); # Passo da posicao
x = -1:dx:1;         # Vetor posicao
Nx = length(x);      # Numero de elementos da posicao
m = 1;               # Massa da particula
h = 1;               # Constante de Planck
k0 = 15;             # Numero de onda
Est = 0;             # Estado de Psi_n
onda_plana = true;   # Define se adiciona onda plana na condicao inicial

## Definicao da funcao de onda Psi(x,t)
Psi = zeros(Nx,Nt);
V = (m/2)*w^2*x.^2; # Potencial do oscilador harmonico
gradV = m*w^2*x;    # Gradiente do potencial

## Condicao inicial e de contorno
if onda_plana == false # Estado estacionario
  Psi(:,1) = psi_n(Est,m,w,h,x);
else # Adiciona onda plana no estado estacionario
  Psi(:,1) = psi_n(Est,m,w,h,x).*exp(i*k0*x);
endif
Psi(1,1) = Psi(Nx,1) = 0; # Condicao de contorno

## Partes real e imaginaria da funcao de onda
RePsi = real(Psi);               # Separar parte real de Psi
ImPsi = imag(Psi);               # Separar parte imaginaria de Psi
RePsi_12 = ImPsi_12 = zeros(Nx); # Definir vetores coluna para meio passo

## Calcular RePsi e ImPsi
for n = 1:Nt-1
  for ix = 2:Nx-1
    
    ## Calcular k1s
    k1R = -r*(ImPsi(ix+1,n)-2*ImPsi(ix,n)+ImPsi(ix-1,n))+V(ix)*ImPsi(ix,n)*dt;
    k1I = r*(RePsi(ix+1,n)-2*RePsi(ix,n)+RePsi(ix-1,n))-V(ix)*RePsi(ix,n)*dt;
    
    ## Atualizar meio passo
    RePsi_12(ix) = RePsi(ix,n) + k1R/2;
    ImPsi_12(ix) = ImPsi(ix,n) + k1I/2;
    
  endfor
  
  for ix = 2:Nx-1
    
    ## Calcular k2s
    k2R = -r*(ImPsi_12(ix+1)-2*ImPsi_12(ix)+ImPsi_12(ix-1))+V(ix)*ImPsi_12(ix)*dt;
    k2I = r*(RePsi_12(ix+1)-2*RePsi_12(ix)+RePsi_12(ix-1))-V(ix)*RePsi_12(ix)*dt;

    ## Atualizar passo
    RePsi(ix,n+1) = RePsi(ix,n) + k2R;
    ImPsi(ix,n+1) = ImPsi(ix,n) + k2I;
    
  endfor
endfor

## Redefinir funcao de onda complexa
Psi = complex(RePsi,ImPsi);
Psi2 = abs(Psi).^2; # |Psi|^2

## Calcular H*Psi e dPsi/dx
HRe = HIm = zeros(Nx,Nt);
dxPsi = zeros(Nx,Nt);
for n = 1:Nt
  for ix = 2:Nx-1
    HRe(ix,n) = -0.5*(RePsi(ix+1,n)-2*RePsi(ix,n)+RePsi(ix-1,n))/(dx^2) + V(ix)*RePsi(ix,n);
    HIm(ix,n) = -0.5*(ImPsi(ix+1,n)-2*ImPsi(ix,n)+ImPsi(ix-1,n))/(dx^2) + V(ix)*ImPsi(ix,n);
    dxPsi(ix,n) = (Psi(ix+1,n)-Psi(ix-1,n))/(2*dx);
  endfor
endfor
HPsi = complex(HRe,HIm);

## Calculo das integrais
CPsi_H_Psi = zeros(1,Nt);     # ConjPsi*H*Psi
CPsi_x_Psi = zeros(1,Nt);     # ConjPsi*x'*Psi
CPsi_p_Psi = zeros(1,Nt);     # ConjPsi*p*Psi
CPsi_gradV_Psi = zeros(1,Nt); # ConjPsi*V*Psi
N_t = zeros(1,Nt);            # Normalizacao
for it = 1:Nt
  CPsi_H_Psi(it) = integracao_13simpson(conj(Psi).*HPsi,dx,it,1,Nx);
  CPsi_x_Psi(it) = integracao_13simpson(conj(Psi).*x'.*Psi,dx,it,1,Nx);
  CPsi_p_Psi(it) = integracao_13simpson(conj(Psi)*(h/i).*dxPsi,dx,it,1,Nx);
  CPsi_gradV_Psi(it) = integracao_13simpson(conj(Psi).*gradV'.*Psi,dx,it,1,Nx);
  N_t(it) = integracao_13simpson(Psi2,dx,it,1,Nx);
endfor

## Calculo dos valores esperados
expv_E = CPsi_H_Psi./N_t;         # Valor esperado da energia
expv_x = CPsi_x_Psi./N_t;         # Valor esperado da posicao
expv_p = CPsi_p_Psi./N_t;         # Valor esperado do momento
expv_gradV = CPsi_gradV_Psi./N_t; # Valor esperado da derivada do potencial

## Derivadas temporais de <x> e <p>
dt_expv_x = dt_expv_p = zeros(1,Nt);
for n = 2:Nt-1
  dt_expv_x(n) = (expv_x(n+1) - expv_x(n-1))/(2*dt);
  dt_expv_p(n) = (expv_p(n+1) - expv_p(n-1))/(2*dt);
endfor

## Figura t=0
figure1 = figure(1,'position',[300 200 720 540]);
[Ax,h1,h2] = plotyy(x,Psi2(:,1),x,V);
xlabel('x','FontSize',18);
ylabel(Ax(1),'|\Psi(x,0)|^2','FontSize',18);
set(Ax(1),'ylim',[0 6]);
ylabel(Ax(2),'V(x)','FontSize',18);
set(Ax(2),'ylim',[0 300]);
set(Ax,'FontSize',18);
set(h2,'LineWidth',1.5);
set(h1,'LineWidth',1.5);
title('Estado fundamental','FontSize',18);

## Figura t=Nt
figure2 = figure(2,'position',[300 200 720 540]);
[Ax,h1,h2] = plotyy(x,Psi2(:,Nt),x,V);
xlabel('x','FontSize',18);
ylabel(Ax(1),'|\Psi(x,501)|^2','FontSize',18);
set(Ax(1),'ylim',[0 6]);
ylabel(Ax(2),'V(x)','FontSize',18);
set(Ax(2),'ylim',[0 300]);
set(Ax,'FontSize',18);
set(h2,'LineWidth',1.5);
set(h1,'LineWidth',1.5);
title('Estado fundamental','FontSize',18);

## Figura Normalizacao
figure3 = figure(3);
plot(t,N_t,'LineWidth',2);
xlim([0 t(Nt)]);
set(gca,'FontSize',16);
title('Normalizacao','FontSize',18);
xlabel('tempo(t)','FontSize',18);
ylabel('N(t)','FontSize',18);

## Figura <E>
figure4 = figure(4);
plot(t,real(expv_E),'LineWidth',1.5,...
     t,imag(expv_E),'LineWidth',1.5);
set(gca,'FontSize',16);
xlim([0 t(Nt)]);
xlabel('t','FontSize',18);
ylabel('<E(t)>','FontSize',18);
title('Valor esperado da energia ','FontSize',18);
h = legend('Parte real','Parte imaginaria','location','east');
set(h,'FontSize',16);

## Figura <x>
figure5 = figure(5);
plot(t,real(expv_x),'LineWidth',1.5,...
     t,imag(expv_x),'LineWidth',1.5);
set(gca,'FontSize',16);
xlim([0 t(Nt)]);
xlabel('t','FontSize',18);
ylabel('<x(t)>','FontSize',18);
title('Valor esperado da posicao','FontSize',18);
h = legend('Parte real','Parte imaginaria');
set(h,'FontSize',16);

## Figura d<x>/dt e <p>
figure6 = figure(6);
plot(t(2:Nt-1),dt_expv_x(2:Nt-1),'LineWidth',1.5,...
     t(2:Nt-1),expv_p(2:Nt-1),'LineWidth',1.5);
set(gca,'FontSize',16);
xlabel('t','FontSize',18);
ylabel('Valores','FontSize',18);
xlim([0 t(Nt)]);
title('Teorema de Ehrenfest Pt1','FontSize',18);
h = legend('d<x>/dt','<p>','location','north');
set(h,'FontSize',16);

## Figura d<p>/dt e -<gradV>
figure7 = figure(7);
plot(t,dt_expv_p,'LineWidth',1.5,...
     t,-expv_gradV,'LineWidth',1.5);
set(gca,'FontSize',16);
xlabel('t','FontSize',18);
ylabel('Valores','FontSize',18);
xlim([0 t(Nt)]);
title('Teorema de Ehrenfest Pt2','FontSize',18);
h = legend('d<p>/dt','-<gradV>','location','northwest');
set(h,'FontSize',16);

## Animacao
figure8 = figure(8);
count = 1;
for n = 1:(int16(Nt*0.016)):Nt
  [Ax,h1,h2] = plotyy(x,Psi2(:,n),x,V);
  xlabel('x','FontSize',16);
  ylabel(Ax(1),'|\Psi(x,t)|^2','FontSize',16);
  set(Ax(1),'ylim',[0 6]);
  #set(hAx(1),'YTick',[ymin:((ymax-ymin)/10):ymax]);
  ylabel(Ax(2),'V(x)','FontSize',16);
  set(Ax(2),'ylim',[0 300]);
  set(Ax(2),'YTick',[0:50:300]);
  title(['Tempo =' num2str(t(n),'%.4f')...
         ' Energia =' num2str(real(expv_E(n)),'%.1f')],...
         'FontSize',16);
  set(Ax,'FontSize',16);
  set(h2,'LineWidth',2);
  pause(1e-4);
  M(count) = getframe(figure8);
  delete(Ax); 
  delete(h1); 
  delete(h2);
  count = count + 1;
endfor