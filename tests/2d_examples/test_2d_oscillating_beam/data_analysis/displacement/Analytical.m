clc
clear 
syms c
time=1.0;
dt=0.01;
n=time/dt;
Ke=zeros(1,n);
Pe=zeros(1,n);
dis_=zeros(1,n);
t_=zeros(1,n);
num_material =1;
%% beam parameter
k = 9.375;
density_ = 1000;
thickness_ = 0.02;
nu_ = 0.4;
EA = 2e6;
EB = 1e6;
%% effective Young's Modulus
if num_material == 1
    E_= EA;
elseif num_material == 2
    E_ = (EA * EA + 14 * EA * EB + EB * EB)/ (8 * (EA + EB));
elseif num_material == 3
    E_ = (7.0 * EA  + EB )/ 8;
end

%% calculate beam frequency
omega_ = (k^4 * E_ * thickness_^2/(12*density_*(1-nu_^2)))^0.5;

%% calculate the amplitude
V_f = 0.01;
cs_ = ((E_/(3-6*nu_))/density_)^0.5;
v_max_ = V_f * cs_;
amp_ = v_max_/omega_;
%% calculate the beam displacement
 for i=1:n
     t_(i)=(i-1)*dt;
     %Ke(i)=(pi/8)*sigma(i)*sigma(i)*a(i)*b(i)*(a(i)*a(i)+b(i)*b(i))*1000;
     %Pe(i)=(pi/8)*1.0*1.0*a(i)*b(i)*(a(i)*a(i)+b(i)*b(i))*1000;
     dis_(i) = amp_ * cos(1.875 - omega_ * t_(i)- 0.304);
 end 
 t=t_';
 dis=dis_';
 %Ke1=Ke';
 %Pe1=Pe';
 %% plot the data
 plot(t,dis);
 filename = "Ana_num_material_" + num2str(num_material) + ".txt";
 fileID = fopen(filename,'w');
 fprintf(fileID,'%6s %6s\r\n','t1','displacement');
 for k=1:n
     fprintf(fileID,'%12.8f %12.8f\r\n',t(k),dis_(k));
 end
 fclose(fileID);
 
 
