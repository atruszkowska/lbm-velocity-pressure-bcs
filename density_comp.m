function rho=density_comp(f,uLN_in,uLS_in,rho_out,Nr,Mc,geom)

    % density computation for micromixer
    rho=sum(f,3);
    
    % density on the south boundary
    fS=f(1,geom.xS+1:geom.xS+geom.W_S-1,:);
    rho(1,geom.xS+1:geom.xS+geom.W_S-1)=(fS(:,:,9)+fS(:,:,1)+fS(:,:,3)+2*(fS(:,:,4)+fS(:,:,7)+fS(:,:,8)))/(1-uLS_in);
    % corners
    rho(1,geom.xS)=rho(1,geom.xS+1); rho(1,geom.xS+geom.W_S)=rho(1,geom.xS+geom.W_S-1);

    % north density - bulk
    fN=f(Nr,geom.xN+1:geom.xN+geom.W_N-1,:);
    rho(Nr,geom.xN+1:geom.xN+geom.W_N-1)=(fN(:,:,9)+fN(:,:,1)+fN(:,:,3)+2*(fN(:,:,2)+fN(:,:,6)+fN(:,:,5)))/(1+uLN_in);
    % corners
    rho(Nr,geom.xN)=rho(Nr,geom.xN+1); rho(Nr,geom.xN+geom.W_N)=rho(Nr,geom.xN+geom.W_N-1);

    % east density (pressure boundary)
    rho(geom.yE:geom.yE+geom.W_E,Mc)=rho_out;
    
end