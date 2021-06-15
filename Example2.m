 
 close all
clear all
global   A B K  Laplacian j_sim na   e x_broad ma  u_s gamma_bar

N=5; %number of agents




% agent model
% A=[0 1; 0 -0.4 ];
% B=[0.8 0.5]';
A=[0 -0.5; 0.5 0];
B=[0;1];
% A=[0 1; 0 0];
% B=[0;1];
% A=rand(2,2);
na=size(A,1);
% B=rand(2,1);
ma=size(B,2);

 
%% initial parameters
tau=0.01; % buffer time
h=0.02; % sampling period
mu=tau+h; 
zeta1=0.3; % exponential rate


Theta_bar=0.3;  
gamma_bar=0.02;
u_s=0.1;
psi_FDI=0.5;

Omega=0.3;

t_sim=0:0.01:20; 

alpha0=3;
alpha1=0.2;
alpha2=0.1;
sample_period=0.01;




%DOS ATTACKS
alpha=3.5;
for i =1:10
    d(i)=(i)+0.5*i*(i)-(1/alpha)*(i);
    d_plus_tau(i)=d(i)+(1/alpha)*(i);
end

 %% Laplacian 

Laplacian=[2 0 0 -1 -1; 0 2 -1 0 -1; 0 -1 2 -1 0; -1 0 -1 2 0; 0 0 0 0 0];

 Laplacian = [  2 -1 0 0 -1;
               -1 2 -1 0  0;
                0 -1 2 -1 0 ;
                0 0 -1 2 -1 ;
                -1 0 0 -1 2 ];
            
%             Laplacian=[2 -1 0 -1; -1 2 -1 0; 0 -1 2 -1; -1 0 -1 2];

L=Laplacian;

[aa,bb]=eig(L);
V_inverse=inv(aa);
J=bb(2:end,2:end); V=V_inverse(2:end,:);


temp=sort(eig(L));
[a,b]=eig(L);
W_inverse=inv(a);
J_tild=b(2:end,2:end); 
V=W_inverse(2:end,:);
J=b(1:end-1,1:end-1);

[aa,bb]=eig(L);
V_inverse=inv(aa);
J=bb(2:end,2:end); V=V_inverse(2:end,:);


%% decision variables
Q=sdpvar(na,na);

Y111=sdpvar(na,na); Y112=sdpvar(na,na,'full'); Y122=sdpvar(na,na);
Y211=sdpvar(na,na); Y212=sdpvar(na,na,'full'); Y222=sdpvar(na,na);
Y311=sdpvar(na,na); Y312=sdpvar(na,na,'full'); Y322=sdpvar(na,na);

Y1=[Y111 Y112; Y112' Y122];
Y2=[Y211 Y212; Y212' Y222];
Y3=[Y311 Y312; Y312' Y322];

P1=sdpvar(na,na);  P2=sdpvar(na,na); P3=sdpvar(na,na); P4=sdpvar(na,na); P5=sdpvar(na,na);

H=sdpvar(na,na,'full');

F1=sdpvar(na,na,'full'); F2=sdpvar(na,na,'full');
F=[F1' F2']';

E1=sdpvar(na,na,'full'); E2=sdpvar(na,na,'full');
E=[E1' E2']';

G1=sdpvar(na,na,'full'); G2=sdpvar(na,na,'full');
G=[G1' G2']';

S1=sdpvar(na,na,'full'); S2=sdpvar(na,na,'full');
S=[S1' S2']';

U1=sdpvar(na,na,'full'); U2=sdpvar(na,na,'full');
U=[U1' U2']';

D1=sdpvar(na,na,'full'); D2=sdpvar(na,na,'full');
D=[D1' D2']';

M=sdpvar(ma,na);


 %%%%%%%%% BLOCKS %%%%%%%%%%%%
W11=kron(eye(N-1),P2+P3+E1+E1'+S1'+S1+A*H'+H*A'+tau*Y111+mu*Y211+(mu-tau)*Y311+2*zeta1*P1);
W12=kron(eye(N-1),G1-E1+E2'+F1-S1+S2'+U1-D1+H*A'+tau*Y112+mu*Y212+(mu-tau)*Y312)+kron(J,B*M); 
W13= kron(eye(N-1),D1-F1);
W14=kron(eye(N-1),-G1-U1);
W15=kron(eye(N-1),P1-H'+H*A'); 
W16=kron(J*V,B*M);

T11=Theta_bar*kron(J*V,B*M);
T12=kron(V,B);
T13=-kron(V,B);
T14=-kron(V,B);
T15= zeros((N-1)*na,(N)*ma);
T16= zeros((N-1)*na,N*ma);

W21=W12';
W22=kron(eye(N-1),F2+F2'+G2+G2'+U2+U2'-S2-S2'-E2-E2'-D2-D2'+tau*Y122+mu*Y222+(mu-tau)*Y322)+kron(J,B*M)+kron(J,B*M)';   
W23=kron(eye(N-1),D2-F2);
W24=kron(eye(N-1),-G2-U2);
W25=-kron(eye(N-1),H')+kron(J,M'*B');  
W26=kron(J*V,B*M);  

T21=Theta_bar*kron(J*V,B*M);
T22=kron(V,B);
T23=-kron(V,B);
T24=-kron(V,B);
T25=gamma_bar*kron(J*V,M');
T26= zeros((N-1)*na,N*ma);
 
W31=W13'; W32=W23';
W33=kron(eye(N-1),-exp(-2*tau*zeta1)*P2);
W34=zeros((N-1)*na,(N-1)*na);  W35=zeros((N-1)*na,(N-1)*na); 
W36=zeros((N-1)*na,(N)*na);    
 
T31=zeros((N-1)*na,(N)*na);   
T32=zeros((N-1)*na,(N)*ma);   
T33=zeros((N-1)*na,(N)*ma);
T34=zeros((N-1)*na,(N)*ma);   
T35=zeros((N-1)*na,(N)*ma);   
T36=zeros((N-1)*na,(N)*ma);

W41=W14'; W42=W24';  W43=W34';
W44=kron(eye(N-1),-exp(-2*mu*zeta1)*P3);
W45=zeros((N-1)*na,(N-1)*na);  W46=zeros((N-1)*na,(N)*na);    

T41=zeros((N-1)*na,(N)*na);   
T42=zeros((N-1)*na,(N)*ma);   
T43=zeros((N-1)*na,(N)*ma);
T44=zeros((N-1)*na,(N)*ma);   
T45=zeros((N-1)*na,(N)*ma);   
T46=zeros((N-1)*na,(N)*ma);


W51=W15'; W52=W25'; W53=W35';   W54=W45';
W55=kron(eye(N-1),tau*P4+mu*P4+(mu-tau)*P5-H'-H);  
W56=kron(J*V,B*M);   

T51=Theta_bar*kron(J*V,B*M);
T52=kron(V,B);
T53=-kron(V,B);
T54=-kron(V,B);
T55=zeros((N-1)*na,(N)*ma);   
T56= zeros((N-1)*na,(N)*ma);  

W61=W16'; W62=W26'; W63=W36'; W64=W46'; W65=W56'; W66=-kron(eye(N),Q);  
 
T61=zeros((N)*na,(N)*na);   
T62=zeros((N)*na,(N)*ma);   
T63=zeros((N)*na,(N)*ma);
T64=zeros((N)*na,(N)*ma);   
T65=zeros((N)*na,(N)*ma);   
T66=gamma_bar*kron(V'*J*V,M');

 
 W=[W11 W12 W13 W14 W15 W16;
    W21 W22 W23 W24 W25 W26;
    W31 W32 W33 W34 W35 W36;
    W41 W42 W43 W44 W45 W46;
    W51 W52 W53 W54 W55 W56;
    W61 W62 W63 W64 W65 W66];

 T=[T11 T12 T13 T14 T15 T16;
    T21 T22 T23 T24 T25 T26;
    T31 T32 T33 T34 T35 T36;
    T41 T42 T43 T44 T45 T46;
    T51 T52 T53 T54 T55 T56;
    T61 T62 T63 T64 T65 T66];

R=blkdiag(-kron(eye(N),H+H')+eye(N*na),-eye(N*ma),-eye(N*ma),-eye(N*ma),-eye((N)*ma),-eye(N*ma));

C1=[ W  T; T'  R];
C2=[Y1 F;F' exp(-2*tau*zeta1)*P4];   C3=[Y1 E;E' exp(-2*tau*zeta1)*P4];  C4=[Y2 G;G' exp(-2*mu*zeta1)*P4];  
C5=[Y2 S;S' exp(-2*mu*zeta1)*P4];  C6=[Y3 U;U' exp(-2*mu*zeta1)*P5];  C7=[Y3 D;D' exp(-2*mu*zeta1)*P5];

zeta2=zeta1*(1-Omega)/Omega;
C8=[A*H'+H*A'-2*zeta2*P1+P2+P3    P1-H'+H*A';   (P1-H'+H*A')'    -H-H'+(tau+mu)*P4+(mu-tau)*P5];

constraints =[ C1<=0, C2>=0, C3>=0, C4>=0, C5>=0, C6>=0, C8<=0, C7>=0, P1>=0, P2>=0, P3>=0, P4>=0,  P5>=0, Q>=0 ];
 info=solvesdp(constraints,[], sdpsettings('solver','sdpt3'))


%%  Design parameters for leaders

K=double(M)*inv(double(H)')
Phi=inv(double(H))*double(Q)*inv(double(H)')
%
eig(kron(eye(N-1),A)+kron(J,B*K))
%%%%%%%%%%%


bar_P1=inv(double(H))*double(P1)*inv(double(H)');
bar_P2=inv(double(H))*double(P2)*inv(double(H)');
bar_P3=inv(double(H))*double(P3)*inv(double(H)');

b1=max(eig(bar_P1))+tau*max(eig(bar_P2))+mu*max(eig(bar_P3));
b2=min(eig(bar_P1));

b_ratio=sqrt(b1/b2);

c=sqrt(N*(u_s^2+psi_FDI^2+alpha0+alpha2)/(2*zeta1));
%%%%%%%%%%%

%% real-time consensus %%%%%%%%%%%%%%%%%%%%%% 
 
for i=1:N
       x_init=[6*i  i]'; %%% initial values

     X_init{i}= x_init;
    X_init_exp{i}=X_init{i}' ;
    x_stack{i}=[];
    u_stack{i}=[];
     event_time{i}=0;
end

%%% Initial values for X

  
  
X_init_save=X_init;

for i=1:N
     x_broad{i}=X_init{i}; x_broad_for_e{i}=X_init{i};
     x_old{i}=X_init{i};   xi5_old{i}=0;
end

counter=zeros(1,N);
Tk=zeros(1,N); 

z_store=[];
rate_upper_bound_store=[];
z_zero=norm(kron(V,eye(na))*reshape(cell2mat(X_init),[1 N*na])');

DoS=0;
DoS_history=0;


for iterations=1:length(t_sim)-2
     x_store=[];
        
for i=1:10
if t_sim(iterations)>d(i) && t_sim(iterations)<d_plus_tau(i)
    DoS=1;
end
end
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
 if DoS ==0  %%%%%%%% No DoS

    for j_sim=1:N
        
             if rand <=gamma_bar
               flag=1;
             else 
             flag=0;
             end
        
          K=double(M)*inv(double(H)');
             
     [t,xx] =ode45(@IS_dif,[t_sim(iterations) t_sim(iterations+1) t_sim(iterations+2) ],x_old{j_sim});
    
     x_stack{j_sim} = [ x_stack{j_sim} xx(2,:)'];

     current_x{j_sim}=xx(2,:)';
     x_store=[x_store; current_x{j_sim} ];
     e{j_sim}=(x_broad_for_e{j_sim})-current_x{j_sim};
     x_old{j_sim}=xx(2,:);
    
    end
   
if mod(iterations,h*100)==0 

for j_sim=1:N
%     z{j_sim}=norm((kron(Laplacian(j_sim,:),eye(na)))*cell2mat((x_broad)'));
      if   e{j_sim}'*Phi*e{j_sim}>alpha0*exp(-alpha1*t_sim(iterations))+alpha2 || DoS_history==1
     
     x_broad_for_e{j_sim}=current_x{j_sim};
     if rand <=Theta_bar
          x_broad{j_sim}=current_x{j_sim}+psi_FDI*tanh(0.1*j_sim*t_sim(iterations));
     else 
           x_broad{j_sim}=current_x{j_sim};
     end
          counter(j_sim)=counter(j_sim)+1; %%%%%%% event isnats
          event_time{j_sim}=[event_time{j_sim} t_sim(iterations)];
          Tk(j_sim)=t_sim(iterations);
      end
 
end

   DoS_history=0;

end

    end
     
    
    %%%%%%%%%%%%%%%%%%%
    
  if DoS ==1   %%%%%%%% in DoS
       
         K=zeros(ma,na);
      
   for j_sim=1:N
          
  [t,xx] =ode45(@IS_dif,[t_sim(iterations) t_sim(iterations+1) t_sim(iterations+2) ],x_old{j_sim});
    

     x_stack{j_sim} = [ x_stack{j_sim} xx(2,:)'];

     current_x{j_sim}=xx(2,:)';
     x_store=[x_store; current_x{j_sim} ];
     x_old{j_sim}=xx(2,:);      
                   
   end
              
          DoS_history=1;
  end   
     
     DoS=0;
    
    

for j_sim=1:N
   u_stack{j_sim} = [ u_stack{j_sim} K*kron(Laplacian(j_sim,:),eye(na))*cell2mat((x_broad)')];
end

z=norm(kron(V,eye(na))*reshape(cell2mat(current_x),[1 N*na])');
z_store=[z_store z];

rate_upper_bound=b_ratio*z_zero*exp(-zeta1*t_sim(iterations))+c;
rate_upper_bound_store=[rate_upper_bound_store rate_upper_bound];

% if <=thresh*norm(kron(V,eye(na))*reshape(cell2mat(X_init),[1 N*na])')
%     break
% end

end
%%

C = {'k','b','r','g','m',[.5 .6 .7],[.8 .2 .6]};



figure(1)
 plot(t_sim(1:iterations),x_stack{1}(1,1:end),'color','b','linewidth',2)
    hold on 
    grid on
     plot(t_sim(1:iterations),x_stack{1}(2,1:end),'color','r','linewidth',2)
      AX=legend('$r_i(t)$','$v_i(t)$','Interpreter','latex','FontSize',20,'fontweight','bold');

 for j_sim=1:N
   fig_h= plot(t_sim(1:iterations),x_stack{j_sim}(1,1:end),'color','b','linewidth',2);
set(fig_h.Annotation.LegendInformation, 'IconDisplayStyle', 'off')  
hold on
    fig_h=  plot(t_sim(1:iterations),x_stack{j_sim}(2,1:end),'color','r','linewidth',2);
    set(fig_h.Annotation.LegendInformation, 'IconDisplayStyle', 'off')  

 end

xlabel('Time (sec)', 'FontSize', 20) 
ylabel('States ', 'FontSize', 20)  
title('(a)', 'FontSize', 20)
 set(gca,'FontSize',20,'fontweight','bold')
  xlim([0, t_sim(iterations)])  
  
  height=ylim;
  height(1)=height(1)-1;
  height(2)=height(2)+1;
  

 for i=1:length(d)
 h1 = line([d(i) d(i)],[height(1) height(2)]);
h2 = line([d_plus_tau(i) d_plus_tau(i)],[height(1) height(2)]);
h1.Annotation.LegendInformation.IconDisplayStyle = 'off';
h2.Annotation.LegendInformation.IconDisplayStyle = 'off';
% Set properties of lines
set([h1 h2],'Color','k','LineWidth',1)
% Add a patch
pp=patch([d(i) d_plus_tau(i) d_plus_tau(i) d(i)],[height(1) height(1) height(2) height(2)],[192 192 192]/255);
pp.Annotation.LegendInformation.IconDisplayStyle = 'off';
% The order of the "children" of the plot determines which one appears on top.
% I need to flip it here.
set(gca,'children',flipud(get(gca,'children')))
 end
patch([d(i) d_plus_tau(i) d_plus_tau(i) d(i)],[height(1) height(1) height(2) height(2)],[192 192 192]/255);
 for j_sim=1:N
   fig_h= plot(t_sim(1:iterations),x_stack{j_sim}(1,1:end),'color','b','linewidth',2);
set(fig_h.Annotation.LegendInformation, 'IconDisplayStyle', 'off')  
hold on
    fig_h=  plot(t_sim(1:iterations),x_stack{j_sim}(2,1:end),'color','r','linewidth',2);
    set(fig_h.Annotation.LegendInformation, 'IconDisplayStyle', 'off')  

 end
ylim([height(1) height(2)])



figure(2)
%
for j_sim=1:N
for i=1:length(event_time{j_sim})
 plot([event_time{j_sim}(i) event_time{j_sim}(i) ],[ j_sim j_sim+0.2], 'color',C{j_sim},'linewidth',0.5)
 hold on
end
 plot([0 iterations/100],[j_sim j_sim], 'color', 'k')
end
%
xlabel('Events (sec)', 'FontSize', 20) 
ylabel('Agent', 'FontSize', 20)  
 set(gca,'FontSize',20,'fontweight','bold')
 xlim([0, t_sim(iterations)])
 set(gca,'FontSize',20,'fontweight','bold')
% end 
hold on
height=ylim;

  height(1)=0.5;
  height(2)=N+0.5;

 for i=1:length(d)
 h1 = line([d(i) d(i)],[height(1) height(2)]);
h2 = line([d_plus_tau(i) d_plus_tau(i)],[height(1) height(2)]);
h1.Annotation.LegendInformation.IconDisplayStyle = 'off';
h2.Annotation.LegendInformation.IconDisplayStyle = 'off';
% Set properties of lines
set([h1 h2],'Color','k','LineWidth',1)
% Add a patch
pp=patch([d(i) d_plus_tau(i) d_plus_tau(i) d(i)],[height(1) height(1) height(2) height(2)],[192 192 192]/255);
pp.Annotation.LegendInformation.IconDisplayStyle = 'off';
% The order of the "children" of the plot determines which one appears on top.
% I need to flip it here.
set(gca,'children',flipud(get(gca,'children')))
 end
patch([d(i) d_plus_tau(i) d_plus_tau(i) d(i)],[height(1) height(1) height(2) height(2)],[192 192 192]/255);
xlabel('Events (sec)', 'FontSize', 20) 
 set(gca,'FontSize',20,'fontweight','bold')
ylim([height(1) height(2)])
xlim([0 t_sim(iterations)])

for j_sim=1:N
for i=1:length(event_time{j_sim})
 plot([event_time{j_sim}(i) event_time{j_sim}(i) ],[ j_sim j_sim+0.2], 'color',C{j_sim},'linewidth',0.5)
 hold on
end
 plot([0 iterations/100],[j_sim j_sim], 'color', 'k')
end




figure(3)
 for j_sim=1:N
   fig_h= plot(t_sim(1:iterations),u_stack{j_sim}(1,1:end),'color',C{j_sim},'linewidth',2);
% set(fig_h.Annotation.LegendInformation, 'IconDisplayStyle', 'off')  
hold on
 end
    grid on
 AX=legend('$u_1(t)$','$u_2(t)$','$u_3(t)$','$u_4(t)$','$u_5(t)$','Interpreter','latex','FontSize',20,'fontweight','bold');
xlabel('Time (sec)', 'FontSize', 20) 
 set(gca,'FontSize',20,'fontweight','bold')
  xlim([0, t_sim(iterations)])  
  


tf=iterations/100
AE=mean(counter(1:N))
AIET=tf/AE
control_max=max(max(abs(cell2mat(u_stack'))))

% [norm(K) norm(Phi) iterations/100  AT  (iterations/100)/AT ]
