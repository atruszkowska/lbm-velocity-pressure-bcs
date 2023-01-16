function [Fx Fy]=Fsum(Nr,Mc,iabw1,jabw1,lena,psi)

%multicomponent interactions - computation of sums;
N_c=9;
F_x=zeros(Nr,Mc,N_c);
F_y=zeros(Nr,Mc,N_c);

C_x=[1 0 -1 0 1 -1 -1 1 0];
C_y=[0 1 0 -1 1 1 -1 -1 0];
w3=1/9; w4=1/36;
Wpsi=[w3 w3 w3 w3 w4 w4 w4 w4 0];
    
    for ic=1:1:N_c-1

         for ia=1:1:lena 
            i=iabw1(ia);  j=jabw1(ia); 
            if i~=1&&j~=1&&i~=Nr&&j~=Mc
            i2 = i+C_y(ic); j2 = j+C_x(ic); 
            F_x(i,j,ic)=Wpsi(ic)*psi(i2,j2)*C_x(ic);
            F_y(i,j,ic)=Wpsi(ic)*psi(i2,j2)*C_y(ic);
            else
            i2 = i; j2 = j; 
            F_x(i,j,ic)=Wpsi(ic)*psi(i2,j2)*C_x(ic);
            F_y(i,j,ic)=Wpsi(ic)*psi(i2,j2)*C_y(ic);  
            end
                      
                     
         end
         
           
    end   

        Fx=sum(F_x(:,:,1:8),3);
        Fy=sum(F_y(:,:,1:8),3);

end