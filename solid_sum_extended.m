function [Fxs Fys]=solid_sum(lena,iabw1,jabw1, Channel2D,Nr,Mc)
s=zeros(Nr,Mc); 
Fxs=zeros(Nr,Mc);
Fys=zeros(Nr,Mc);
N_c=9;
F_xs=zeros(Nr,Mc,N_c);
F_ys=zeros(Nr,Mc,N_c);

C_x=[1 0 -1 0 1 -1 -1 1 0];
C_y=[0 1 0 -1 1 1 -1 -1 0];
w3=1/9; w4=1/36;
Wpsi=[w3 w3 w3 w3 w4 w4 w4 w4 0];

for ic=1:1:N_c-1
	
		for ia=1:1:lena % obstacle
            
            i=iabw1(ia);  j=jabw1(ia); 
            if (i>=2&&i<=(Nr-1)&&j>=2&&j<=(Mc-1))
                i3 = i+C_y(ic); j3 = j+C_x(ic);                            

                if Channel2D(i3,j3)==0

                    s(i3,j3)=1;

                else

                    s(i3,j3)=0;

                end    

                F_xs(i,j,ic)=Wpsi(ic)*C_x(ic)*s(i3,j3);
                F_ys(i,j,ic)=Wpsi(ic)*C_y(ic)*s(i3,j3);
            end  
            
            if j==Mc||j==1
                i3 = i+C_y(ic); j3=j;                            

                if Channel2D(i3,j3)==0

                    s(i3,j3)=1;

                else

                    s(i3,j3)=0;

                end    

                F_xs(i,j,ic)=Wpsi(ic)*C_x(ic)*s(i3,j3);
                F_ys(i,j,ic)=Wpsi(ic)*C_y(ic)*s(i3,j3);
            end
            

        end
            
end 

	    Fxs=sum(F_xs(:,:,1:8),3);
        Fys=sum(F_ys(:,:,1:8),3);

end
