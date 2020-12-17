function  dx=IS_dif(t,x)

  global  A B K x_broad Laplacian j_sim  na  sint  ma u_s gamma_bar
  
dx=zeros(na,1) ;
sint=sin(j_sim*10*t);
cost=cos(j_sim*5*t);
 dx=A*x+(1-gamma_bar*sint)*B*K*kron(Laplacian(j_sim,:),eye(na))*cell2mat((x_broad)')+ B*u_s*cost; 
