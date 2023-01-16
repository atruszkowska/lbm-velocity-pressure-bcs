function f=mixer_BCs(rho,f,uLS_in,uLN_in,rho_out,Nr,Mc,geom)

% Mixer BCs

% compute the corner unknowns
% east boundary
% lower corner
f(geom.yE,Mc,5)=0.5*(rho(geom.yE+1,Mc)-(f(geom.yE,Mc,1)+f(geom.yE,Mc,2)+f(geom.yE,Mc,3)+f(geom.yE,Mc,4)+f(geom.yE,Mc,6)+f(geom.yE,Mc,8)+f(geom.yE,Mc,9)));
f(geom.yE,Mc,7)=f(geom.yE,Mc,5);
% upper corner
f(geom.yE+geom.W_E,Mc,6)=0.5*(rho(geom.yE+geom.W_E-1,Mc)-(f(geom.yE+geom.W_E,Mc,1)+f(geom.yE+geom.W_E,Mc,2)+f(geom.yE+geom.W_E,Mc,3)+f(geom.yE+geom.W_E,Mc,4)+f(geom.yE+geom.W_E,Mc,5)+f(geom.yE+geom.W_E,Mc,7)+f(geom.yE+geom.W_E,Mc,9)));
f(geom.yE+geom.W_E,Mc,8)=f(geom.yE+geom.W_E,Mc,6);

% new BCs
% East - outlet
fE=f(geom.yE+1:geom.yE+geom.W_E-1,Mc,:);
uE=-1+(fE(:,:,9)+fE(:,:,2)+fE(:,:,4)+2*(fE(:,:,1)+fE(:,:,5)+fE(:,:,8)))/rho_out;
ruE=uE*rho_out;

f(geom.yE+1:geom.yE+geom.W_E-1,Mc,3)=fE(:,:,1)-2/3*ruE;
f(geom.yE+1:geom.yE+geom.W_E-1,Mc,6)=fE(:,:,8)-1/6*ruE-0.5*(fE(:,:,4)-fE(:,:,2));
f(geom.yE+1:geom.yE+geom.W_E-1,Mc,7)=fE(:,:,5)-1/6*ruE-0.5*(fE(:,:,2)-fE(:,:,4));

% lower corner
f(geom.yE,Mc,3)=f(geom.yE,Mc,1)-2/3*ruE(1);
f(geom.yE,Mc,6)=f(geom.yE,Mc,8)-1/6*ruE(1)-0.5*(f(geom.yE,Mc,4)-f(geom.yE,Mc,2));
% upper corner
f(geom.yE+geom.W_E,Mc,3)=f(geom.yE+geom.W_E,Mc,1)-2/3*ruE(end);
f(geom.yE+geom.W_E,Mc,7)=f(geom.yE+geom.W_E,Mc,5)-1/6*ruE(end)-0.5*(f(geom.yE+geom.W_E,Mc,2)-f(geom.yE+geom.W_E,Mc,4));

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% compute the corner unknowns
% south left
f(1,geom.xN,6)=0.5*(rho(1,geom.xN)-(f(1,geom.xN,9)+f(1,geom.xN,1)+f(1,geom.xN,2)+f(1,geom.xN,3)+f(1,geom.xN,4)+f(1,geom.xN,5)+f(1,geom.xN,7)));
f(1,geom.xN,8)=f(1,geom.xN,6);
% north right 
f(Nr,geom.xS+geom.W_S,6)=0.5*(rho(Nr,geom.xS+geom.W_S)-(f(Nr,geom.xS+geom.W_S,9)+f(Nr,geom.xS+geom.W_S,1)+f(Nr,geom.xS+geom.W_S,2)+f(Nr,geom.xS+geom.W_S,3)+f(Nr,geom.xS+geom.W_S,4)+f(Nr,geom.xS+geom.W_S,5)+f(Nr,geom.xS+geom.W_S,7)));
f(Nr,geom.xS+geom.W_S,8)=f(Nr,geom.xS+geom.W_S,6);        
% south right 
f(1,geom.xN+geom.W_N,5)=0.5*(rho(1,geom.xN+geom.W_N)-(f(1,geom.xN+geom.W_N,9)+f(1,geom.xN+geom.W_N,1)+f(1,geom.xN+geom.W_N,2)+f(1,geom.xN+geom.W_N,3)+f(1,geom.xN+geom.W_N,4)+f(1,geom.xN+geom.W_N,6)+f(1,geom.xN+geom.W_N,8)));
f(1,geom.xN+geom.W_N,7)=f(1,geom.xN+geom.W_N,5);
% north left
f(Nr,geom.xN,5)=0.5*(rho(Nr,geom.xN)-(f(Nr,geom.xN,9)+f(Nr,geom.xN,1)+f(Nr,geom.xN,2)+f(Nr,geom.xN,3)+f(Nr,geom.xN,4)+f(Nr,geom.xN,6)+f(Nr,geom.xN,8)));
f(Nr,geom.xN,7)=f(Nr,geom.xN,5);

% new BCs
% % % % % % % % %     
rho_in=(f(Nr,geom.xS+1:geom.xS+geom.W_S-1,9)+f(Nr,geom.xS+1:geom.xS+geom.W_S-1,1)+f(Nr,geom.xS+1:geom.xS+geom.W_S-1,3)+2*(f(Nr,geom.xS+1:geom.xS+geom.W_S-1,2)+f(Nr,geom.xS+1:geom.xS+geom.W_S-1,5)+f(Nr,geom.xS+1:geom.xS+geom.W_S-1,6)))/(1+uLN_in);
f(Nr,geom.xS+1:geom.xS+geom.W_S-1,4)=f(Nr,geom.xS+1:geom.xS+geom.W_S-1,2)-2/3*rho_in*uLN_in;
f(Nr,geom.xS+1:geom.xS+geom.W_S-1,7)=f(Nr,geom.xS+1:geom.xS+geom.W_S-1,5)+0.5*(f(Nr,geom.xS+1:geom.xS+geom.W_S-1,1)-f(Nr,geom.xS+1:geom.xS+geom.W_S-1,3))-1/6*rho_in*uLN_in;
f(Nr,geom.xS+1:geom.xS+geom.W_S-1,8)=f(Nr,geom.xS+1:geom.xS+geom.W_S-1,6)-0.5*(f(Nr,geom.xS+1:geom.xS+geom.W_S-1,1)-f(Nr,geom.xS+1:geom.xS+geom.W_S-1,3))-1/6*rho_in*uLN_in;

f(Nr,geom.xS,4)=f(Nr,geom.xS,2)-2/3*rho_in(1)*uLN_in;
f(Nr,geom.xS,8)=f(Nr,geom.xS,6)-0.5*(f(Nr,geom.xS,1)-f(Nr,geom.xS,3))-1/6*rho_in(1)*uLN_in;

f(Nr,geom.xS+geom.W_S,4)=f(Nr,geom.xS+geom.W_S,2)-2/3*rho_in(end)*uLN_in;
f(Nr,geom.xS+geom.W_S,7)=f(Nr,geom.xS+geom.W_S,5)+0.5*(f(Nr,geom.xS+geom.W_S,1)-f(Nr,geom.xS+geom.W_S,3))-1/6*rho_in(end)*uLN_in;

% % % % % % % % %     
rho_out=(f(1,geom.xN+1:geom.xN+geom.W_N-1,9)+f(1,geom.xN+1:geom.xN+geom.W_N-1,1)+f(1,geom.xN+1:geom.xN+geom.W_N-1,3)+2*(f(1,geom.xN+1:geom.xN+geom.W_N-1,4)+f(1,geom.xN+1:geom.xN+geom.W_N-1,7)+f(1,geom.xN+1:geom.xN+geom.W_N-1,8)))/(1-uLS_in);
f(1,geom.xN+1:geom.xN+geom.W_N-1,2)=f(1,geom.xN+1:geom.xN+geom.W_N-1,4)+2/3*rho_out*uLS_in;
f(1,geom.xN+1:geom.xN+geom.W_N-1,5)=f(1,geom.xN+1:geom.xN+geom.W_N-1,7)-0.5*(f(1,geom.xN+1:geom.xN+geom.W_N-1,1)-f(1,geom.xN+1:geom.xN+geom.W_N-1,3))+1/6*rho_out*uLS_in;
f(1,geom.xN+1:geom.xN+geom.W_N-1,6)=f(1,geom.xN+1:geom.xN+geom.W_N-1,8)-0.5*(-f(1,geom.xN+1:geom.xN+geom.W_N-1,1)+f(1,geom.xN+1:geom.xN+geom.W_N-1,3))+1/6*rho_out*uLS_in;

f(1,geom.xN,2)=f(1,geom.xN,4)+2/3*rho_out(1)*uLS_in;
f(1,geom.xN,5)=f(1,geom.xN,7)-0.5*(f(1,geom.xN,1)-f(1,geom.xN,3))+1/6*rho_out(1)*uLS_in;

f(1,geom.xN+geom.W_N,2)=f(1,geom.xN+geom.W_N,4)+2/3*rho_out(end)*uLS_in;
f(1,geom.xN+geom.W_N,6)=f(1,geom.xN+geom.W_N,8)-0.5*(-f(1,geom.xN+geom.W_N,1)+f(1,geom.xN+geom.W_N,3))+1/6*rho_out(end)*uLS_in;
end
