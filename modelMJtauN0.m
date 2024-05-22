clear;  

global sd;
global N0;
global tau;
global rmax;
global rb;
global R0;

sd=0.02; %um; synaptic cleft
rb=.2; %um; radius of postsynaptic density
N0=83*6; %mol/um; vesicle glutamate content;
R0=83*2; %uM; initial concentration of receptor channel
tau=0.15*N0/498; %ms; decay time constant of glutamate release - modified to scale with N0 - 498 -> 6000
% D = 0.3 SEE LINE 170

% N0=1600/(6.022*10^23)/0.02*10^(15)*10^(6); 133mol/um; 1600 molecules/vesicle; 83mol/um gives 1000 molecules  
% 2500/(6.022*10^23)/0.02*10^(15)*10^(6); 208mol/um
% 2000/(6.022*10^23)/0.02*10^(15)*10^(6); 166mol/um
% R0=2000/(6.022*10^23*0.0max-rise-rate2)*(10^5)^3*10^6; 166uM; 2000 molecules/um^2
%----------------------------------------------

% a spatial mesh of (n_r+n_r2-1) points and a time mesh of n_t points
rmax=5; % um; larger than rb originally .4
n_r=500;
rmax2=50;  % originally 5
n_r2=100;   % these 4 numbers define the grid and variations indicated it gives accurate results

tmax=10; % ms
n_t=5000;

r1 = linspace(0,rmax,n_r);
r2 = linspace(rmax,rmax2,n_r2);
r = [r1 r2(2:end)];
t = linspace(0,tmax,n_t);
%----------------------------------------------


% Solve Equation
m = 0; %symmetry constant
sol = pdepe(m,@difpde,@pdex1ic,@pdex1bc,r,t);
% N1:A ; N2:R ; N3:AR ; N4:A2R ; N5:A2O ; N6:AD ; N7:A2D
N1 = sol(:,:,1);
N2 = sol(:,:,2);
N3 = sol(:,:,3);
N4 = sol(:,:,4);
N5 = sol(:,:,5);
N6 = sol(:,:,6);
N7 = sol(:,:,7);
%----------------------------------------------


% figure% 3D surface plot
% surf(r,t,N1)
% xlabel('r/\mum', 'FontSize',12)
% ylabel('t/ms', 'FontSize',12)
% zlabel('[A]/\muM', 'FontSize',12)
% title('A', 'FontSize',16)

% figure
% surf(r,t,N5)
% xlabel('r/\mum', 'FontSize',12)
% ylabel('t/ms', 'FontSize',12)
% zlabel('[A2O]/\muM', 'FontSize',12)
% title('A2O', 'FontSize',16)
%----------------------------------------------


% P0(t): average of A2O over postsynaptic area
P0=zeros(1,n_t);
dr=rmax/(n_r-1);
for i=1:n_t%time
    aver=0;
    for j=1:n_r%space
        aver=aver+2*pi*r(j)*dr*N5(i,j);
    end
    P0(i)=aver/(pi*rb^2);
end
%----------------------------------------------

% figure% A2O vs time
% plot(t,P0);
% xlabel('t/ms', 'FontSize',12)
% ylabel('Averaged A2O', 'FontSize',12)
%----------------------------------------------


% convert to postsynaptic current; holding voltage = 65 mV; single channel conductance = 7.6 pS
current=-65*7.6*10^(-3)*P0*10^(-6)*6.022*10^23*(10^(-5))^3*0.02*pi*rb^2; % postsynaptic current
[M,I]=max(P0);
Amp=-65*7.6*10^(-3)*M*10^(-6)*6.022*10^23*(10^(-5))^3*0.02*pi*rb^2; % amplitude
area=-trapz(t,current);
mEPSC=fittype(@(A1,a1,a2,x) -A1*(1-exp(-x/a1)).*exp(-x/a2) );
f=fit(t.',current.',mEPSC,'StartPoint',[10,1,1],'Lower',[0,0,0]);
coefficientValues = coeffvalues(f);
decaytime=coefficientValues(3);
b1=coefficientValues(2);
b2=coefficientValues(3);
tpeak=b1*log((b2/b1)+1);
Amp
Max_open_probability = M/R0;
[maxP0 , maxP0_i] = max(N5(:,1)/R0);
%maxP0
%maxP0_i
[maxP5 , maxP5_i] = max(N5(:,51)/R0);
%maxP5
%maxP5_i
%----------------------------------------------


% figure %postsynaptic current vs time
% plot(f,t.',current.')
hold on
plot(t.',current.')  % plot only current without fit
xlabel('t/ms', 'FontSize',12)
ylabel('I/pA', 'FontSize',12)
grid on
%----------------------------------------------

filecurrent = fopen('current.txt','w');% Write data to text file
fprintf(filecurrent,'%f\n',current);
fclose(filecurrent);

% rise time and normalized max rise rate
tmp1=0.1*M;
tmp2=0.9*M;
i=1;
while P0(i)<tmp1
    i=i+1;
end
j=i;
while P0(j)<tmp2
    j=j+1;
end

risetime10_90=t(j)-t(i);
risetime10_90
rise_rate=diff(P0)/(tmax/n_t);
max_rise_rate=max(rise_rate)/M;
%----------------------------------------------


%file1 = fopen('result-var-rb.txt','a');% Write data to text file
file1 = fopen('result.txt','a');% Write data to text file
fprintf(file1,'%f\t%f\t%f\t%f\t%f\t%d\t%d\t%d\t%d\t%d\n',-Amp,risetime10_90,area,max_rise_rate,decaytime,Max_open_probability,rb,N0,R0,tau);
fclose(file1);

%file1 = fopen('result2023.txt','a');% Write data to text file
%fprintf(file1,'%f\t%f\t%f\t%f\t%d\t%d\t%d\t%d\n',Amp,risetime10_90,area,decaytime, maxP0,maxP5,N0,R0);
%fclose(file1);
%----------------------------------------------


function fac=func(r)% Approximation to fac=1(r<rb); fac=0(r>rb)
    global rb;
    fac=1/(1+exp(500*(r-rb)));
end
%----------------------------------------------


function [c,f,s] = difpde(r,t,N,DNDr)% Equations to solve

% N1:A ; N2:R ; N3:AR ; N4:A2R ; N5:A2O ; N6:AD ; N7:A2D

global sd;
global N0;
global tau;

factor=.0485;  % scale factor to sloe rise-time - value of 1 means no change to GluR2 kinetics
D=0.3;    % um^2/ms; diffusion coefficient
% D1=D/0.8;   %  D=0.8 for bulk; so on rates a reduced proportionally to D; with bulk D=0.8 there is no reduction
k1=0.008; % from Krampfl et al., 2002, k1 = 0.008 and k2 = .004
k2=0.004; % on rates doubled to give lower EC50's ; Holm et al 2005
k3=0.004; %multiply by D1 to scale on-rate with D. Assume that Krampfl measurements were with D=.8
         % incorporating D into rates gave opposite result of dextran in Nielsen paper
            %Change on rates to 0.01; faster needed to produce some response
a=2.4;      % 5 in Jinbo's 
b=20*factor;    % b=20;
k_1=2;
k_2=4*factor;   % k_2=4;
k_3=0.0622;
d1=4.5*factor;  % d1=4.5;
d_1=0.007;
d2=0.6;
d_2=0.06;

c = [1; 1 ;1 ;1 ;1 ;1 ;1];
f = [D;0;0;0;0;0;0] .* DNDr; 
s = [D * DNDr(1)/r+k2*N(3)+N0/tau* exp(-t/tau)/(2*pi*sd^2)*exp(-(r).^2/(2*sd^2)) +  k_1*N(3) - k1*N(1)*N(2) + k_2*N(4) - k2*N(1)*N(3) + k_3*N(7) - k3*N(1)*N(6); ...
	k_1*N(3)-k1*N(1)*N(2);...
	k1*N(1)*N(2) + k_2*N(4) + d_2*N(6) - k_1*N(3) - k2*N(1)*N(3) - d2*N(3) ;...
	k2*N(1)*N(3) + d_1*N(7) + a*N(5) - k_2*N(4) - d1*N(4) - b*N(4) ;...
	b*N(4)-a*N(5) ;...
	d2*N(3) + k_3*N(7) - d_2*N(6) - k3*N(1)*N(6) ;...
	k3*N(1)*N(6) + d1*N(4) - k_3*N(7) - d_1*N(7) ];
end
%----------------------------------------------


function u0 = pdex1ic(r)% Initial Conditions
global R0;
u0 = [0;R0*func(r);0;0;0;0;0];
end
%----------------------------------------------


function [pl,ql,pr,qr] = pdex1bc(xl,ul,xr,ur,t)% Boundary Conditions
pl = [0;0;0;0;0;0;0];
ql = [1;1;1;1;1;1;1];
pr = [0; 0; 0; 0; 0; 0; 0];
qr = [1; 1; 1; 1; 1; 1; 1];
end
%----------------------------------------------