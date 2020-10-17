#########################################################
# Calcular a magnetizacao media do modelo de Ising em 2D
# em funcao de h utilizando o algoritmo de metropolis
#########################################################

clear;

function config = gera_config(N)
  
  % Gera uma configuracao aleatoria NXN de spins -1 e 1
  % utilizando a funcao sinal em uma distribuicao normal de
  % N numeros aleatorios com media 0 e variancia 1
  
  config = sign(randn(N));
  
endfunction

function Eflip = calc_Eflip(config, n, m, J, h)
  
    % Calculo da diferenca de energia(Eflip) entre uma configuracao  
    % onde config(n,m) foi flipado e a configuracao original
    
    N = length(config);
    
    Eflip = 2*J*config(n,m)*(config(mod(n,N)+1,m) + ...
                             config(n,mod(m,N)+1) + ...
                             config(mod(n-2,N)+1,m) + ...
                             config(n,mod(m-2,N)+1)) + 2*h*config(n,m);
                             
endfunction

function result = somavizinhos(config)
  
  % Retorna uma matriz onde cada elemento (n,m) e a soma
  % dos primeiros vizinhos de config(n,m)
  
  result = circshift(config,[0 1]) ...
                   + circshift(config,[0 -1]) ...
                   + circshift(config,[1 0]) ...
                   + circshift(config,[-1 0]);
  
endfunction

function Em = calc_Em(config, J, h)
  
  % Calcula a energia media por sitio de uma configuracao de spins
  
  aux = config .* somavizinhos(config);
  soma_J = sum(aux(:)); # Soma do termo de J
  soma_h = sum(config(:)); # Soma do termo de h
  E_tot = 0.5*(- J * soma_J - h * soma_h); # Energia total da configuracao
  Em = E_tot/numel(config);
  
endfunction

function [mag_med, E_med, C, susc] = metropolis(config, J, h, T, Nvar)
  
  % Varre Nvar vezes uma configuracao de spins com o algoritmo de metropolis 
  % retornando a magnetizacao media por sitio, a energia media por sitio,
  % o calor específico e a suscetibilidade magnética, todas essas
  % grandezas em funcao da temperatura T
  
  config_nova = config;    # Cria nova configuracao baseada na original
  N = numel(config_nova);  # Numero de elementos da configuracao
  elementos_config = 1:N;  # Vetor com numero dos elementos
  perm_elementos = elementos_config(randperm(N)); # Permutar elementos
  m_alpha = zeros(1,Nvar); # Vetor com magnetizacao total para cada varredura
  E_alpha = zeros(1,Nvar); # Vetor com energia por sitio de cada varredura
  
  ## Fazer Nvar varreduras na configuracao
  for ivar = 1:Nvar
    
    ## Iterar aleatoriamente pelos spins
    for ind = 1:N
      [col,lin]  = ind2sub(size(config_nova),perm_elementos(ind));
      Eflip = calc_Eflip(config_nova,lin,col, J, h); # Calculo de Eflip
      if Eflip < 0
        config_nova(lin,col) *= -1; # Flipar spins com baixas energias
      else
        Pflip = exp(-Eflip/T); # Probabilidade de Boltzmann para Eflip
        if rand() < Pflip
          config_nova(lin,col) *= -1; # Flipar spins que aumentaram pouco a energia
        endif
      endif
    endfor
    
    m_alpha(ivar) = mean(config_nova(:));     # Armazenar densidade de magnetizacao
    E_alpha(ivar) = calc_Em(config_nova,J,h); # Armazenar energia media por sitio
    
  endfor
  
  mag_med = mean(m_alpha); # Magnetizacao media por sitio das configuracoes
  E_med = mean(E_alpha);   # Energia media por sitio das configuracoes
  C = var(E_alpha)/(T^2);  # Calor especifico
  susc = var(m_alpha)/T;   # Suscetibilidade 
  
endfunction

## Definicoes iniciais
J = 1;                   # Constante de troca 
h = 0;                   # Campo magnetico externo
Nspins = [10];           # Numero de spins em cada dimensao da configuracao    
Nvar = [1e3];            # Numero de varreduras da mesma configuracao
T = linspace(0.1,4,35);  # Valores de temperatura
Tc = 2*J/log(1+sqrt(2)); # Temperatura critica
NT = length(T);          # Dimensao do vetor temperatura

## Inicializar vetores
mag_med = zeros(1,NT);  # Vetor para as magnetizacoes medias por temperatura
E_med = zeros(1,NT);    # Vetor para as energias medias por temperatura
C_med = zeros(1,NT);    # Vetor para o calor especifico por temperatura
susc_med = zeros(1,NT); # Vetor para a suscetibilidade magnetica por temperatura

## Calculo principal
for Ns = 1:length(Nspins)
  for Nv = 1:length(Nvar)
    for iT = 1:NT # Iterar sobre diferentes temperaturas
      config = gera_config(Nspins(Ns)); # Gera a configuracao dos spins
      [mag_m,E_m,C_m,S_m] = metropolis(config,J,h,T(iT),Nvar(Nv));
      mag_med(iT) = mag_m;
      E_med(iT) = E_m;
      C_med(iT) = C_m;
      susc_med(iT) = S_m;
      fprintf('Ns=%i Nvar=%i T=%.2f\n',Nspins(Ns),Nvar(Nv),T(iT))
    endfor
    arq = strcat('N',num2str(Nspins(Ns)),'-NVAR',num2str(Nvar(Nv)));
    save(arq,'T','mag_med','E_med','C_med','susc_med')
  endfor
endfor

## Calculo dos expoentes criticos
#TTc = (T-Tc);
#p1 = polyfit(log(TTc),log(mag_med),1);

## Figura magnetizacao media x temperatura
figure1 = figure(1);
plot(T,abs(mag_med),'-o','linewidth',2);
xlabel('T','FontSize',18);
ylabel('<m>(T)','FontSize',18);
title('Magnetização média x Temperatura','FontSize',18);
set(gca,'FontSize',18);

## Figura energia media x temperatura
figure2 = figure(2);
plot(T,E_med,'-o','linewidth',2);
xlabel('T','FontSize',18);
ylabel('<E>(T)','FontSize',18);
title('Energia média x Temperatura','FontSize',18);
set(gca,'FontSize',18);

## Figura calor especifico medio x temperatura
figure3 = figure(3);
plot(T,C_med,'-o','linewidth',2);
xlabel('T','FontSize',18);
ylabel('C(T)','FontSize',18);
title('Calor específico por temperatura','FontSize',18);
set(gca,'FontSize',18);

## Figura suscetibilidade 
figure4 = figure(4);
plot(T,susc_med,'-o','linewidth',1.5);
xlabel('T','FontSize',18);
ylabel('\chi (T)','FontSize',18);
title('Suscetibilidade magnetica por temperatura','FontSize',18);
set(gca,'FontSize',18);