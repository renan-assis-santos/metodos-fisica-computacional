############################################################
# Resolucao numerica da equacao de Poisson em duas dimensões
# Para uma carga pontual entre duas placas FINITAS
# De um capacitor em um determinado potencial
############################################################

## Definicoes iniciais
L = 1;             # Lado do quadrado da região de dominio [m]
dx = dy = 0.05;    # Passo para as posicoes
x = 0:dx:L;        # Vetor posicao x
y = 0:dy:L;        # Vetor posicao y
Nx = length(x);    # Num elementos em x
Ny = length(y);    # Num elementos em y

## Caracteristicas da carga
Q_eps0 = 2;             # Valor de Q/eps_0 [N*m^2/C]
x_q = 0.5*L;            # Posicao x da carga
y_q = 0.5*L;            # Posicao y da carga
iq = round(x_q/dx)+1;   # Indice da posicao x da carga
jq = round(y_q/dx)+1;   # Indice da posicao y da carga
f = zeros(Nx,Ny);       # Grid para conter a carga (f = rho/eps_0)
f(iq,jq) = Q_eps0/dx^2; # Apenas no lugar da carga valor e diferente de 0

#### Opcao para mais cargas
##Q_eps0 = [2 -2];
##x_q = [0.9*L 0.1*L];
##y_q = [0.5*L 0.5*L];
##iq = round(x_q/dx)+1;
##jq = round(y_q/dx)+1;
##f = zeros(Nx,Ny);
##for num = 1:length(Q_eps0)
##  f(iq(num),jq(num)) = Q_eps0(num)/dx^2;
##endfor

## Inicializar vetor potencial
## Condicoes de contorno V = 0 ja aplicadas
## Chute inicial sendo V1(i,j) = 0
V = zeros(Nx,Ny);

## Condicoes de contorno para as placas
V(7,5:17) = -1;  # Placa da esquerda com -1V
V(15,5:17) = +1; # Placa da direita com 1V

## Criterios para convergencia
erro = 1e-3; # Erro maximo
DeltaV = 1;  # abs(V(i,j)-V_old(i,j))
n = 1;       # Numero de iteracoes,
nmax = 600;  # Numero maximo de iteracoes

## Loop principal
while DeltaV >= erro && n <= nmax # Criterios de parada
  
  V_old = V;           # Configuracao anterior do potencial
  DeltaV_old = DeltaV; # DeltaV anterior
  DeltaV = 0;          # Reiniciar DeltaV
  
  ## Preencher potencial
  for i = 2:Nx-1
    for j = 2:Ny-1
      V(i,j) = (1/4)*(V_old(i+1,j)+V_old(i-1,j)+V_old(i,j+1)+V_old(i,j-1)+f(i,j)*(dx^2));
    endfor
  endfor
  
  ## Manter condicoes de contorno para as placas a cada iteracao
  V(7,5:17) = -1;
  V(15,5:17) = +1;
  
  ## Atualizar DeltaV
  DeltaV = sum(sum(abs(V-V_old)));

  ## Atualizar iteracao
  n = n + 1;

endwhile

## Calculo campo eletrico
[Ex Ey] = gradient(-1*V',dx,dx);

## Plot curvas de nivel e campo eletrico
figure1 = figure(1);                 # Gera a primeira Figura
[c,w] = contour(x,y,V',50); hold on; # Gera as curvas de nivel do potencial
set(w,'Linewidth',1.5)               # Caracteristicas das curvas de nivel
barracor = colorbar;                 # Barra de cor para as curvas de nivel
set(barracor,'FontSize',16);         # Caracteristicas da barra de cor
if length(Q_eps0) == 1               # Titulo para uma carga
  title('Capacitor com uma carga pontual','FontSize',16);
else                                 # Titulo para mais cargas
  title('Capacitor com mais de uma carga pontual','FontSize',16);
endif
xlabel('x (m)','FontSize',16);       # Legenda eixo x
ylabel('y (m)','FontSize',16);       # Legenda eixo y
set(gca,'FontSize',16);              # Caracteristicas dos eixos
quiver(x,y,Ex,Ey,1.1); hold off;     # Campo vetorial do campo eletrico