function [C] = MatrixC(mu)
global I J dz;
% 2.1.1. Constructing the Sparse Matrix for the Employed
chi =  - min(mu,0)/dz;
ups =  - max(mu,0)/dz + min(mu,0)/dz;
zet =  + max(mu,0)/dz;
% 2.1.1.1. Upperdiagonal
ud = zeros(I,1);
for j = 1:J-1
    ud = [ud; repmat(zet(j), I, 1)];
end
% 2.1.1.2. Centerdiagonal
cd = repmat(chi(1) + ups(1),I,1);
for j = 2:J-1
    cd = [cd; repmat(ups(j),I,1)];
end
cd = [cd; repmat(ups(J)+zet(J), I, 1)];
% 2.1.1.3. Lowerdiagonal
ld = repmat(chi(2), I, 1);
for j = 3:J
    ld=[ld; repmat(chi(j), I, 1)];
end
ld = [ld; zeros(I,1)];
C = spdiags(cd,0,I*J,I*J) + spdiags(ld,-I,I*J,I*J) + spdiags(ud,I,I*J,I*J);
end

