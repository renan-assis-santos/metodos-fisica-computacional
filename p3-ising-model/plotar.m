arq = {'N10-NVAR10000';'N15-NVAR10000';'N20-NVAR10000'}
#arq = {'N10-NVAR1000';'N10-NVAR5000';'N10-NVAR10000'}
for ind = 1:length(arq)
  load(arq{ind})
endfor


Tc = 2*J/log(1+sqrt(2));
ttc = (T-Tc)/Tc;
p1 = polyfit(log(ttc(20:45)),log(susc_med_20_1e4(20:45)),1)

## Figura magnetizacao media pela temperatura
figure1 = figure(1);
plot(T,abs(mag_med_10_1e4),'-o','linewidth',2,...
     T,abs(mag_med_15_1e4),'-o','linewidth',2,...
     T,abs(mag_med_20_1e4),'-o','linewidth',2);
xlabel('T','FontSize',18);
ylabel('<m>(T)','FontSize',18);
title('Magnetização média x Temperatura','FontSize',18);
set(gca,'FontSize',18);
h = legend('N = 10','N = 15','N = 20');
set(h,'FontSize',18);

## Log log
figure10 = figure(10)
loglog(ttc,abs(mag_med_10_1e4),'-o','linewidth',2,...
       ttc,abs(mag_med_15_1e4),'-o','linewidth',2,...
       ttc,abs(mag_med_20_1e4),'-o','linewidth',2);
xlabel('T','FontSize',18);
ylabel('<m>(T)','FontSize',18);
title('Magnetização média x Temperatura','FontSize',18);
set(gca,'FontSize',18);
h = legend('N = 10','N = 15','N = 20');
set(h,'FontSize',18);

## Figura energia media pela temperatura
figure2 = figure(2);
plot(T,E_med_10_1e4,'-o','linewidth',2,...
     T,E_med_15_1e4,'-o','linewidth',2,...
     T,E_med_20_1e4,'-o','linewidth',2);
xlabel('T','FontSize',18);
ylabel('<E>(T)','FontSize',18);
title('Energia média x Temperatura','FontSize',18);
set(gca,'FontSize',18);
h = legend('N = 10','N = 15','N = 20');
set(h,'FontSize',18);

## Figura calor especifico medio pela temperatura
figure3 = figure(3);
plot(T,C_med_10_1e4,'-o','linewidth',2,...
     T,C_med_15_1e4,'-o','linewidth',2,...
     T,C_med_20_1e4,'-o','linewidth',2);
xlabel('T','FontSize',18);
ylabel('C(T)','FontSize',18);
title('Calor específico por temperatura','FontSize',18);
xlim([1 4]);
##ylim([0 1]);
set(gca,'FontSize',18);
h = legend('N = 10','N = 15','N = 20');
set(h,'FontSize',18);

## Figura suscetibilidade magnetica medio pela temperatura
figure4 = figure(4);
plot(T,susc_med_10_1e4,'-o','linewidth',2,...
     T,susc_med_15_1e4,'-o','linewidth',2,...
     T,susc_med_20_1e4,'-o','linewidth',2);
xlabel('T','FontSize',18);
ylabel('\chi (T)','FontSize',18);
title('Suscetibilidade magnetica por temperatura','FontSize',18);
set(gca,'FontSize',18);
h = legend('N = 10','N = 15','N = 20');
set(h,'FontSize',18);