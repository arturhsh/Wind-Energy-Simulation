%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INEP                                          %%% 
%%% Semestre 2021.2 e 2022/1                      %%%
%%% Artur Henrique dos Santos Heerdt              %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Analise espectral e simulacao daS componentes %%%
%%% longa e rápida do vento.                      %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear 
clc 
format long

%% Tabela de Variáveis
%%% T: ; Table_50m: ; wind_data: ; Fs: ; Amp: ; Freq: ; 
%%% Phase: ; Fmax: ; Ts: ; Tsim: ; vel: ; wc: ; m1: ; m2: ;
%%% L: ; slope_regression: ; Ts: ; Tf: ; sd: ; Kf: ; V: ; M: ;
%%% N: ; r: ; k: ; m: ; Dw: ; Vel_med: ; P: ; H: ; sys: ; out: ; 
%%% intensidade_turb: ; ws: ; J: ; S: ; Kf2: .

%% Obtencao e transformacao da tabela de dados para array
T = readtable('C:\Users\artur\OneDrive\Documents\INEP\csv_vento_2019.csv');
Table_50m = T(:,6);
wind_data = table2array(Table_50m); %variÃ¡vel que contem vetor com os dados horÃ¡rios do vento

%% Análise espectral (componente longa)

Fs = 1/3600; %uma hora -> Frequencia de amostragem dos dados

[Amp Freq Phase] = fft_(wind_data, Fs); %Dados retirados da Transformada RÃ¡pida de Fourier

Fmax = Fs; %Estamos interessados em periodos acima de uma hora (10 minutos no exemplo do Matheus)

Ts = 1/Fs; %Periodo de amostragem

Tsim = 24*Ts*7;  


vel = Amp(1)*ones(1,floor(Tsim/Ts));

wcl = 2*pi*rand(1,length(Freq));


m1 = 0.4;
m2 = 0.25;

for i=1:length(vel)
    for j=2:length(Freq)
        vel(i) = vel(i)+ Amp(j)*cos(2*pi*Freq(j)*(i-1)*Ts+wcl(j)); %velocidade dos dados + espectral
    end
end


%% Simulacao turbulencia

L = 200;
slope_regression = 0.2; 
Ts = 1;

Tf = L./vel;

sd = slope_regression*vel;
Kf = ((2*pi*Tf/(Ts*beta(1/2,1/3)))).^(1/2); % B -> Beta Function

%% Equacoes 28 e 29 - Nichita
Ts = 1;
ws = 2*pi/Ts;
Dw = 0.002;

J = ws/(2*Dw); 

S = zeros(length(Tf),1);

for i=1:length(Tf)
    for k=1:J
        S(i) = S(i) + (((m1*Tf(i)*k*Dw)^2)+1)/(((Tf(i)*k*Dw)^2+1)*(((m2*Tf(i)*k*Dw)^2)+1));
    end
    Kf2(i) = (pi/(Ts*Dw*S(i)))^(1/2);
end

%% Continuacao Turbulencia

V = zeros(3600*length(vel),1);
M = 5000;
N = 100;
r = 1:M; %vector
k = 1:N+1; %vector 
m = 1:N;

Dw = 0.002;

for f=1:length(vel)
    wct =(rand(3600,1)-0.5)*2*sqrt(3);

    Vel_med((f-1)*3600+1:f*3600) = vel(f);
    
    V((f-1)*3600+1:f*3600) = vel(f); %velocidade espectral + turbulencia
    
    s = tf('s');
    sys(f) = Kf2(f)*(m1*Tf(f)*s+1)/((Tf(f)*s+1)*(m2*Tf(f)*s+1));

    b = c2d(sys(f), Ts, 'tustin');
    [n, d] = tfdata(b, 'v');


for k=1:length(wct)
      if(k==1)
         out(k) = n(1)*wct(k);
      elseif (k==2)
        out(k) =  n(1)*wct(k) + n(2)*wct(k-1)- d(2)*out(k-1) ; 
      else
  
      out(k) = n(1)*wct(k) + n(2)*wct(k-1) + n(3)*wct(k-2) - d(2)*out(k-1) - d(3)*out(k-2);
      end
end
out = out*sd(f);
   
    for m=1:3600   
        for i=1:1
            V((f-1)*3600+m) = V((f-1)*3600+m) + out(m); %trocar por H(i) out(m)
            if V((f-1)*3600+m) < 1e-2
                V((f-1)*3600+m) = 1e-2;
            end  
            if V((f-1)*3600+m) > 11
                V((f-1)*3600+m) = 11;
            end
        end

        
    end
    
    intensidade_turb(f) = std(V((f-1)*3600+1:f*3600))/mean(V((f-1)*3600+1:f*3600));
    if intensidade_turb(f) >=1, intensidade_turb(f) = 1; end

end

mean(intensidade_turb,'omitnan');

figure(1)
subplot(4,1,1)
plot(V,'b','Linewidth',2)
ylim([0 20])
grid on
 
% subplot(4,1,2)
% plot(vel,'b','Linewidth',2)
% grid on
% xlim([0 60])
% ylim([0 20])
%  
% 
% subplot(4,1,3)
% plot(wind_data, 'b', 'Linewidth', 2)
% 
% subplot(4,1,4)
% plot(intensidade_turb, 'b', 'Linewidth', 2)


%%
% plot(out)


max(wind_data)
max(V)

