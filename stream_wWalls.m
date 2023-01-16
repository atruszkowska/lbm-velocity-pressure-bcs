function f=stream_wWalls(uLN_in,uLS_in,rho,f,Nr,Mc,lenwint,iawint,jawint,lenobs,iobs,jobs,Channel2D,geom)
N_c=9; feq=zeros(Nr,Mc,N_c);
%Bounce back
ic_op=[3 4 1 2 7 8 5 6];
C_x=[1 0 -1 0 1 -1 -1 1 0];
C_y=[0 1 0 -1 1 1 -1 -1 0];
 % Streaming
    feq=f;
    for ic=1:1:N_c-1
        ic2=ic_op(ic); % bounce-back BC
        temp=feq(:,:,ic);
        
     % All nodes except boundaries (covers the diagonal directions on
     % the boundaries that are suppose to be known)
     % non-bordering fluid locations
        for ia=1:1:lenwint
            i=iawint(ia); j=jawint(ia);
            if (j~=1&&j~=Mc)&&(i~=1&&i~=Nr)
                i2=i+C_y(ic); j2=j+C_x(ic);
                f(i2,j2,ic)=temp(i,j);
            end
        end
        
        % fluid locations on the borders of obstacles
        for ia=1:1:lenobs 
            i=iobs(ia);  j=jobs(ia); 
            i2 = i+C_y(ic); j2 = j+C_x(ic); 
            
            if (j~=Mc)&&(i2>1&&i2<Nr&&j2<Mc)
               if (Channel2D(i2,j2)==0) 
                    f(i,j,ic2) =temp(i,j); 
               else 
                    f(i2,j2,ic)=temp(i,j); 
               end 
            end
            
        % East boundary
            if (j==Mc)&&(ic==2||ic==4)
               if (Channel2D(i2,j2)==0) 
                    f(i,j,ic2) =temp(i,j); 
               else 
                    f(i2,j2,ic)=temp(i,j); 
               end 
            end
            
            if (j==Mc&&i>geom.yE&&i<geom.yE+geom.W_E)&&(ic==3||ic==6||ic==7)
                    f(i2,j2,ic)=temp(i,j); 
            end
           
            if (j==Mc&&i==geom.yE)&&(ic==3||ic==6)
                    f(i2,j2,ic)=temp(i,j); 
            end
           
            if (j==Mc&&i==geom.yE+geom.W_E)&&(ic==3||ic==7)
                    f(i2,j2,ic)=temp(i,j); 
            end
         % North/south boundaries
           if (i==1||i==Nr)&&(ic==1||ic==3)
               if (Channel2D(i2,j2)==0) 
                    f(i,j,ic2) =temp(i,j); 
               else 
                    f(i2,j2,ic)=temp(i,j); 
               end 
           end
           
           if (i==1&&j>geom.xN&&j<geom.xN+geom.W_N)&&(ic==6||ic==2||ic==5)
                    f(i2,j2,ic)=temp(i,j); 
           end
           
           if (i==1&&j==geom.xN)&&(ic==2||ic==5)
                    f(i2,j2,ic)=temp(i,j); 
           end
           
           if (i==1&&j==geom.xN+geom.W_N)&&(ic==2||ic==6)
                    f(i2,j2,ic)=temp(i,j); 
           end
           
           if (i==Nr&&j>geom.xS&&j<geom.xS+geom.W_S)&&(ic==7||ic==4||ic==8)
               f(i2,j2,ic)=temp(i,j); 
           end
           
           if (i==Nr&&j==geom.xS)&&(ic==4||ic==8)
               f(i2,j2,ic)=temp(i,j); 
           end
           
           if (i==Nr&&j==geom.xS+geom.W_S)&&(ic==7||ic==4)
               f(i2,j2,ic)=temp(i,j); 
           end
           
         end
   end
end

