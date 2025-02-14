function [A] = MatrixA(ssB,ssF,IsB,IsF)
    global I J daaF daaB;
    X = - IsB.*ssB ./ daaB;
    Y = - IsF.*ssF ./ daaF + IsB.*ssB ./ daaB;
    Z = IsF.*ssF ./ daaF;
    
    X(I,:) = -min(ssB(I,:),0)./daaB(I,:);
    Y(I,:) = -max(ssF(I,:),0)./daaF(I,:) + min(ssB(I,:),0)./daaB(I,:);
    Z(I,:) =  max(ssF(I,:),0)./daaF(I,:);
        ZZ = [0]; %This is needed because of the peculiarity of spdiags.
        for j=1:J
            ZZ = [ZZ ; Z(1:I-1,j);0];
        end
        
        YY = reshape(Y,I*J,1);
        
        XX = X(2:I,1);
        for j=2:J
            XX = [XX;0;X(2:I,j)];
        end
        A = spdiags(YY,0,I*J,I*J)+spdiags([ZZ;0],1,I*J,I*J)+spdiags([XX;0],-1,I*J,I*J);

end

