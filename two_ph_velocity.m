% T-mixer, two-phase flow

% Fluid 1 - dispersed
% Fluid 2 - continuous  

% saving tag
NUMBER=1;

% Channel setup
Lc=150;
Len_Channel_2D=Lc; % y direction
Channel_2D_half_width=300; Width=Channel_2D_half_width*2; % x direction
Channel2D=ones(Len_Channel_2D,Width); %Wet area
[Nr Mc]=size(Channel2D);

% Mixer geometry
% widths of inlets (north, west, south), actual width +1 (halfway
% bounceback)
W_N=60;
W_S=60;
% positions of inlets
% upper left corner
xN=2; yN=Nr;
% lower left corner
xS=2; yS=1;
% width of outlet (actual is +1)
W_E=60;
% lower left corner
xE=Mc; yE=Lc/2-W_S/2; 
% wrap it into structure
mix_geom.W_N=W_N;
mix_geom.W_S=W_S;
mix_geom.W_E=W_E;

mix_geom.xN=xN;
mix_geom.xS=xS;
mix_geom.xE=xE;

mix_geom.yN=yN;
mix_geom.yS=yS;
mix_geom.yE=yE;

% wall setup
Channel2D(:,1)=0; Channel2D(1:Lc/2-W_N/2-1,xN+W_N+1:Mc)=0;
Channel2D(Lc/2+W_N/2+1:Nr,xS+W_S+1:Mc)=0;

%Fluid and flow properties
cs2=1/3; % LBM constant
nu=1/6;
tau1=1;%3*nu+1/2;
tau2=1;%3*nu+1/2;
nu_1=(tau1-0.5)/3;
nu_2=(tau2-0.5)/3;

omega_1=1/tau1;
omega_2=1/tau2;
G=0.9;%0.625; %Interaction strength, these are repulsive interactions hence G>0
Gs_1=0.15; %Interactions with solids, "-" is attractive, "+", repulsive
Gs_2=-Gs_1;
% Normal velocity components
% Nope: 1:1
uLN_in=-1.3284e-03; % continuous fluid; 
uLS_in=3.8411e-04; % dispersed phase ;

% Indexing
Channel2D=logical(Channel2D);
% obstacles for Bounce Back 
Obstacles=bwperim(Channel2D,8); % perimeter of the grains for bounce back Bound.Cond.
% Wet locations etc.
[iabw1 jabw1]=find(Channel2D==1); % indices i,j, of active lattice locations 
lena=length(iabw1); % number of active location 
ija= (jabw1-1)*Nr+iabw1; % equivalent single index (i,j)->> ija for active locations
% absolute (single index) position of the obstacles in for bounce back in Channel2D
% Obstacles 
[iobs jobs]=find(Obstacles);lenobs=length(iobs); ijobs= (jobs-1)*Nr+iobs; % as above
% Internal wet locations : wet & ~obstables
% (i.e. internal wet lattice location non in contact with dray locations)
[iawint jawint]=find(( Channel2D==1 & ~Obstacles)); % indices i,j, of active lattice locations
lenwint=length(iawint); % number of internal (i.e. not border) wet locations
ijaint= (jawint-1)*Nr+iawint; % equivalent single index
NxM=Nr*Mc;

% Density - initial
density=2.0; %fluid density
density_1=zeros(Nr,Mc); density_2=zeros(Nr,Mc);
density_1(1:Lc/2-W_S/2,xS:xS+W_S)=2.0; density_2(find(density_1==0&Channel2D~=0))=2.0;
density_temp=density_1;
density_1(find(density_2~=0&Channel2D~=0))=0.06;
density_2(find(density_temp~=0&Channel2D~=0))=0.06;

% Directions
East=1; North=2; West=3; South=4; NE=5; Nw=6; SW=7; SE=8; RP=9;
N_c=9;
C_x=[1 0 -1 0 1 -1 -1 1 0];
C_y=[0 1 0 -1 1 1 -1 -1 0];
C=[C_x;C_y];
% Bounce back
ic_op=[3 4 1 2 7 8 5 6];
% Density weights
w0=16/36.;w1=4/36.;w2=1/36.;
W=[w1 w1 w1 w1 w2 w2 w2 w2 w0];
w3=1/9; w4=1/36;
Wpsi=[w3 w3 w3 w3 w4 w4 w4 w4 0];
cs2=1/3;cs2x2=2*cs2;cs4x2=2*cs2.^2;
f1=3.; f2=4.5; f3=1.5;
% Collisions constants
cs2=1/3;cs2x2=2*cs2;cs4x2=2*cs2.^2;
f1=3.; f2=4.5; f3=1.5;

% Declarative statements, DZ courtesy
f_1=zeros(Nr,Mc,N_c); 
f_2=zeros(Nr,Mc,N_c); 
temp=zeros(Nr,Mc,N_c);
feq=zeros(Nr,Mc,N_c);
feq_1=zeros(Nr,Mc,N_c);
feq_2=zeros(Nr,Mc,N_c);
Fxtemp_1=zeros(Nr,Mc);
Fxtemp_2=zeros(Nr,Mc);
Fytemp_1=zeros(Nr,Mc);
Fytemp_2=zeros(Nr,Mc);

ux=zeros(Nr,Mc);
uy=zeros(Nr,Mc);
rho_1=zeros(Nr,Mc);
rho_2=zeros(Nr,Mc);
Fx_1=zeros(Nr,Mc);
Fx_2=zeros(Nr,Mc);
Fy_1=zeros(Nr,Mc);
Fy_2=zeros(Nr,Mc);
Fxs_1=zeros(Nr,Mc); 
Fxs_2=zeros(Nr,Mc);
Fys_1=zeros(Nr,Mc); 
Fys_2=zeros(Nr,Mc);
ux_1=zeros(Nr,Mc);uy_1=zeros(Nr,Mc);
ux_2=zeros(Nr,Mc);uy_2=zeros(Nr,Mc);

% Initialization - start values in the wet area
for ia=1:lena % all the active locations
    i=iabw1(ia); j=jabw1(ia); % matrix coordinates of each active location
    f_1(i,j,:)=density_1(i,j)/9; 
    f_2(i,j,:)=density_2(i,j)/9;
end
rho=density_1+density_2; 

% The same for both species
[Fxs Fys]=solid_sum_extended(lena,iabw1,jabw1, Channel2D,Nr,Mc);
% This is due to the fact that psi=rho, so density term in the forcing term will cancel
% Different psi will cause it to be included in the loop
% So this will be incorporated into the computation of equilibrium velocity
% as tau_1*Fxs_1, no division by density, same for species 2
Fxs_1=-Gs_1*Fxs; Fys_1=-Gs_1*Fys; 
Fxs_2=-Gs_2*Fxs; Fys_2=-Gs_2*Fys;

% While - main time evolution loop
Max_Iter=2.7e6;
Cur_Iter=0;
StopFlag=false;
% initial outlet
rho_1=density_1;
rho_2=density_2;
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
        save(['slug_droplet_kerosene_water_',num2str(NUMBER),'.mat'])
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
% end