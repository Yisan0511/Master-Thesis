function [zefrKFE] = KFE4z(mue, muu, lamu, lame)

global zgrid J dz;

chi =  - min(mue,0)/dz;
ups =  - max(mue,0)/dz + min(mue,0)/dz;
zet =  + max(mue,0)/dz;
% 2.1.1.1. Upperdiagonal
ud = zeros(1,1);
for j = 1:J-1
    ud = [ud; zet(j)];
end
% 2.1.1.2. Centerdiagonal
cd = chi(1) + ups(1);
for j = 2:J-1
    cd = [cd; ups(j)];
end
cd = [cd; ups(J)+zet(J)];
% 2.1.1.3. Lowerdiagonal
ld = chi(2);
for j = 3:J
    ld=[ld; chi(j)];
end
ld = [ld; zeros(1,1)];
Ze = spdiags(cd,0,J,J) + spdiags(ld,-1,J,J) + spdiags(ud,1,J,J);

%--------------------------------------------------------------------------
chi =  - min(muu,0)/dz;
ups =  - max(muu,0)/dz + min(muu,0)/dz;
zet =  + max(muu,0)/dz;
% 2.1.1.1. Upperdiagonal
ud = zeros(1,1);
for j = 1:J-1
    ud = [ud; zet(j)];
end
% 2.1.1.2. Centerdiagonal
cd = chi(1) + ups(1);
for j = 2:J-1
    cd = [cd; ups(j)];
end
cd = [cd; ups(J)+zet(J)];
% 2.1.1.3. Lowerdiagonal
ld = chi(2);
for j = 3:J
    ld=[ld; chi(j)];
end
ld = [ld; zeros(1,1)];
Zu = spdiags(cd,0,J,J) + spdiags(ld,-1,J,J) + spdiags(ud,1,J,J);

ZZZ = blkdiag(sparse(Ze), sparse(Zu));
%--------------------------------------------------------------------------
Lam_ul = spdiags(-lame*ones(J,1),0,J,J);
Lam_ur = spdiags( lame*ones(J,1),0,J,J);
Lam_dl = spdiags( lamu*ones(J,1),0,J,J);
Lam_dr = spdiags(-lamu*ones(J,1),0,J,J);
Lam = [Lam_ul, Lam_ur; Lam_dl, Lam_dr];
Z = ZZZ + Lam;

M = J*2;
fix = 1;
ZT = Z';
null = zeros(M,1);
null(fix) = 1;
ZT(fix,:) = [zeros(1,fix-1),1,zeros(1, M-fix)];

% 5.2.2. Solve linear system
gVEC = ZT\null;
gsum = (sum(gVEC))*dz;
g = gVEC./ gsum;
ge = g(1:J);
zefrKFE = (zgrid' * ge)/(sum(ge));
end

