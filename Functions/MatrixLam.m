function [Lam] = MatrixLam(lamu)
    global I J lame;
    Lam_ul = spdiags(-lame*ones(I*J,1),0,I*J,I*J);
    Lam_ur = spdiags( lame*ones(I*J,1),0,I*J,I*J);
    Lam_dl = spdiags( lamu*ones(I*J,1),0,I*J,I*J);
    Lam_dr = spdiags(-lamu*ones(I*J,1),0,I*J,I*J);
    Lam = [Lam_ul, Lam_ur; Lam_dl, Lam_dr];
end

