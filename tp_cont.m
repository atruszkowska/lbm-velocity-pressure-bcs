% While - main time evolution loop
load('second_kerosene_water_ST_117.mat')
Max_Iter=2.7e5;
Cur_Iter=Cur_Iter;
StopFlag=false;
NUMBER=118;
tic
while (~StopFlag)
    Cur_Iter=Cur_Iter+1;
       
    % Composite density distribution function and macroscopic density
    rho=rho_1+rho_2; 
    rho_out_1=mean(rho_1(yE:yE+W_E,Mc)); 
    rho_out_2=mean(rho_2(yE:yE+W_E,Mc));
   
    % Species densities (based on species distribution functions)
    rho_1=density_comp(f_1,uLN_in,uLS_in,rho_out_1,Nr,Mc,mix_geom);
    rho_2=density_comp(f_2,uLN_in,uLS_in,rho_out_2,Nr,Mc,mix_geom);

    % Interaction potential (equal to density)  
    psi_1=rho_1; psi_2=rho_2;
    % Interaction forces summation computations, separate function file
    [Fxtemp_1 Fytemp_1]=Fsum(Nr,Mc,iabw1,jabw1,lena,psi_1);
    [Fxtemp_2 Fytemp_2]=Fsum(Nr,Mc,iabw1,jabw1,lena,psi_2);
    % Interaction forces, final form
    Fx_1=-G*psi_1.*Fxtemp_2; Fx_2=-G*psi_2.*Fxtemp_1;
    Fy_1=-G*psi_1.*Fytemp_2; Fy_2=-G*psi_2.*Fytemp_1;
    
        if Cur_Iter>1
             ux_1=zeros(Nr,Mc);uy_1=zeros(Nr,Mc); ux_2=zeros(Nr,Mc);uy_2=zeros(Nr,Mc);
            for ic=1:N_c-1
 
                ux_1(:,:)=ux_1(:,:)+C_x(ic).*f_1(:,:,ic);
                uy_1(:,:)=uy_1(:,:)+C_y(ic).*f_1(:,:,ic);
                ux_2(:,:)=ux_2(:,:)+C_x(ic).*f_2(:,:,ic);
                uy_2(:,:)=uy_2(:,:)+C_y(ic).*f_2(:,:,ic);

            end
        end
    
    % Velocity - boundaries correction
    % north
    ux_1(Nr,xN:xN+W_N)=0; uy_1(Nr,xN:xN+W_N)=uLN_in;
    ux_2(Nr,xN:xN+W_N)=0; uy_2(Nr,xN:xN+W_N)=uLN_in;
    % south 
    ux_1(1,xS:xS+W_S)=0; uy_1(1,xS:xS+W_S)=uLS_in;
    ux_2(1,xS:xS+W_S)=0; uy_2(1,xS:xS+W_S)=uLS_in;
    % east (pressure BC)
    fE=f_1(yE:yE+W_E,Mc,:);
    ux_1(yE:yE+W_E,Mc)=-1+(fE(:,:,9)+fE(:,:,2)+fE(:,:,4)+2*(fE(:,:,1)+fE(:,:,5)+fE(:,:,8)))/rho_out_1;
    uy_1(:,Mc)=0;
    
    fE=f_2(yE:yE+W_E,Mc,:);
    ux_2(yE:yE+W_E,Mc)=-1+(fE(:,:,9)+fE(:,:,2)+fE(:,:,4)+2*(fE(:,:,1)+fE(:,:,5)+fE(:,:,8)))/rho_out_2;
    uy_2(:,Mc)=0;
            
    % Composite velocity
    ux=(ux_1*omega_1+ux_2*omega_2)./(rho_1*omega_1+rho_2*omega_2);
    uy=(uy_1*omega_1+uy_2*omega_2)./(rho_1*omega_1+rho_2*omega_2);
    
	% Equilibrium velocities, if conditions to avoid division by 0
    if rho_1(ija)>0
        ux_1(ija)=ux(ija)+Fx_1(ija)./rho_1(ija)/omega_1+Fxs_1(ija)/omega_1; 
        uy_1(ija)=uy(ija)+Fy_1(ija)./rho_1(ija)/omega_1+Fys_1(ija)/omega_1;
    end
    
    if rho_2(ija)>0
        ux_2(ija)=ux(ija)+Fx_2(ija)./rho_2(ija)/omega_2+Fxs_2(ija)/omega_2; 
        uy_2(ija)=uy(ija)+Fy_2(ija)./rho_2(ija)/omega_2+Fys_2(ija)/omega_2;
    end
    
    % Function file for computing equilibrium distributions
    feq_1=f_eqilibrium(ux_1,uy_1,rho_1,ija,NxM,Nr,Mc,N_c);
    feq_2=f_eqilibrium(ux_2,uy_2,rho_2,ija,NxM,Nr,Mc,N_c);
    
	% Retrieving the macroscopic momentum (velocity after division with density, 
	% but this is not needed here since that�s the form used for computation of 
	% composite velocity and equilibrium velocity. Correct this (never mind here), 
	% i t doesn�t compute u at t=1 (which is fine, but still).      
	
    if rho_1(ija)>0
        ux_1(ija)=ux_1(ija)./rho_1(ija);
        uy_1(ija)=uy_1(ija)./rho_1(ija);
    elseif rho_2(ija)>0
        ux_2(ija)=ux_2(ija)./rho_2(ija);
        uy_2(ija)=uy_2(ija)./rho_2(ija);
    end
    
	% Collisions per species
	
    f_1=(1.-omega_1).*f_1+omega_1.*feq_1;
    f_2=(1.-omega_2).*f_2+omega_2.*feq_2;
    
	% Streaming in a separate function file, each species separately	
    f_1=stream_wWalls(uLN_in,uLS_in,rho_1,f_1,Nr,Mc,lenwint,iawint,jawint,lenobs,iobs,jobs,Channel2D,mix_geom);
    f_2=stream_wWalls(uLN_in,uLS_in,rho_2,f_2,Nr,Mc,lenwint,iawint,jawint,lenobs,iobs,jobs,Channel2D, mix_geom);
    
    % Boundary conditions
    f_1=mixer_BCs(rho_1,f_1,uLS_in,uLN_in,rho_out_1,Nr,Mc,mix_geom);
    f_2=mixer_BCs(rho_2,f_2,uLS_in,uLN_in,rho_out_2,Nr,Mc,mix_geom);
    
    %Checking for NaNs
    if isnan(mean(mean(rho_2)))
        disp('NaNs');
        break
    end
    
     if (mod(Cur_Iter,500)==0)
        save(['second_kerosene_water_ST_',num2str(NUMBER),'.mat'])
        NUMBER=NUMBER+1;
    end
    %Breaking condition check:
        if  (Cur_Iter > Max_Iter)
            StopFlag=true;
            break 
        end  
%       surf(rho_1), daspect([1 1 1]), view([0 90])
%       pause(1)
end
toc
