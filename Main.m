%% ------------------------------------------------------------------------
%% Part 1: Preparation                                                    |
%% ------------------------------------------------------------------------

%% 0. Housekeeping
clear;
addpath('Figures');
addpath('Functions');

%% 1. Structural Parameter Settings

    % 1.1. Structural Parameters
    global I J daaF daaB dz lame amin amax gam aagrid zgrid;
    alf  = 0.3;     % Production Parameter          fr. Bardoczy
    bet  = 0.72;    % Worker's bargaining power     fr. Bardoczy    
    gam  = 1;       % Utility function              fr. Ben Moll Aiyagari Diffusion
    del  = 0.021;   % Depreciation Rate             fr. Ben Moll Aiyagari Diffusion
    eta  = 0.72;    % LM Tightness Elasticity       fr. Bardoczy
    phi  = 0.395;   % Flow Cost of Vacancy          fr. Bardoczy (It's \xi in Bardoczy)
    rho  = 0.01;    % Discount rate                 fr. Ben Moll Aiyagari Diffusion
    chi  = 1.7935;  % Matching Efficiency           fr. Bardoczy                
    lame = 0.1038;  % Separation Intensity          fr. Bardoczy and several other liter0ature
    hp   = 0.0001;  % Home Production

    % 1.2. Asset Grid
    amin = 0; 
    amax = 2500;
    I = 1500;
    J = 15;
    x = linspace(0,1,I)';
    coeff = 1.3; power = 2.8;
    xx  = x + coeff*x.^power;
    xmax = max(xx); 
    xmin = min(xx);
    agrid = (amax-amin)/(xmax - xmin)*xx + amin;
    aagrid = agrid*ones(1,J);
    daF = ones(I,1);
    daB = ones(I,1);
    daF(1:I-1) = agrid(2:I)-agrid(1:I-1);
    daB(2:I  ) = agrid(2:I)-agrid(1:I-1);
    daF(I) = daF(I-1); 
    daB(1) = daB(2);
    daaF = repmat(daF,1,J);
    daaB = repmat(daB,1,J);

    % 1.3. Productivity Grid
    ze =  1; 
    zu = 0.65;
    zmin = zu; 
    zmax = ze;
    zgrid = linspace(zmin, zmax, J)';
    zzgrid = ones(I,1) * zgrid';
    zgridREP= repelem(zgrid,I);
    dz = (zmax - zmin)/(J-1);
    the  = .25;
    theu = .50;
    mue = the*(ze - zgrid);
    muu = the*(zu - zgrid);

%% 2. The Transition Matrix for the Diffusion Process  
Ce = MatrixC(mue);
Cu = MatrixC(muu);
C = blkdiag(sparse(Ce), sparse(Cu));

%% 3. Solution Parameter Settings
tolHJB = 1e-10;  % tolerance level for HJB equations 12?
tolENT = 1e-3;  % 5 tolerance level for free entry condition
tolMKT = 1e-3;  % 5 tolerance level for asset market clearing
n = 300; % maximum number of iterations when solving HJB equation
N = 1000;  % maximum number of iterations in outer loop
Delta = 500;  % step size in HJB equation
dTheta = 1e-1;  % step size in updating theta
relax = .999;
kgap = 500;

%% ------------------------------------------------------------------------
%% Part 2: Solution                                                       | 
%% ------------------------------------------------------------------------
Tau = [0, 0.01, 0.02];

tic;
% (The Loop for Different Tax Rates)
for Taxation = 1:length(Tau) % First loop
    tau = Tau(Taxation);
    % Reset
    Delta = 500;
    relax = .999;
    kgap = 500;
    % 1. Initial Guess

        % 1.1. Updating Parameters

            % 1.1.1. Labor Market
            Theta = 1; 
            lamu = chi * (Theta ^ (1 - eta));
            q = chi * (Theta^(-eta));
            u = lame / (lame + lamu);
            v = Theta * u;

            % 1.1.2. Capital Market and Prices
            z = KFE4z(mue,muu,lamu,lame);
            k = (((rho+del)/(z*alf))^(1/(alf-1)))*1.01;
            r = z*alf*k^(alf-1);
            wj = bet * ((k^alf)* zgrid' - r * k);
            wg = repmat(wj,I,1);
            wgVEC = wg(:);

        % 1.2. Value functions
        We  = zeros(I,  J);
        Wu  = zeros(I,  J);

        WeF = zeros(I, J);
        WeB = zeros(I, J);
        WuF = zeros(I, J);
        WuB = zeros(I, J);

        Jf  = (repelem(zgrid,I) *k^alf - r*k - wgVEC) / (lame + r - del);
        
        % 1.3. Unemployment benefits
        upop = (lame)/(lamu+lame);
        [epop,upop] = (zDistribution(mue,muu,lamu,lame));
        govrev = (tau*wj)*epop;
        h = govrev/(sum(upop));
        h = h*ones(I,J);

% 2. Solution
% (The Outer Loop)
    for outerloop = 1:N

        % 2.1. Solving Consumer's HJB
        % (The Inner Loop)
        for innerloop = 1:n
            
            % 2.1.1. The Employed
            WeB(2:I,   :) = (We(2:I,   :)-We(1:I-1, :)) ./ (aagrid(2:I,:) - aagrid(1:I-1,:));
            WeF(1:I-1, :) = (We(2:I,   :)-We(1:I-1, :)) ./ (aagrid(2:I,:) - aagrid(1:I-1,:));
            WeB(1,:) = (wg(1,:)*(1-tau) + (r-del)*amin).^(-gam);
            WeF(I,:) = (wg(I,:)*(1-tau) + (r-del)*amax).^(-gam);
            ceB = max(WeB, 1e-6).^(-1/gam);
            ceF = max(WeF, 1e-6).^(-1/gam);
            ce0 = wg*(1-tau) + (r-del)*aagrid;
            Complex_Check(ceB);
            Complex_Check(ceF);
            seB = wg*(1-tau) + (r-del)*aagrid - ceB;
            seF = wg*(1-tau) + (r-del)*aagrid - ceF;
            HeF = log(ceF) + WeF.*seF;
            HeB = log(ceB) + WeB.*seB;
            IeEither = (1-(seF>0)) .* (1-(seB<0)); % using central 
            IeUnique = (seB<0).*(1-(seF>0)) + (1-(seB<0)).*(seF>0); % Ideal
            IeBoth   = (seB<0).*(seF>0); % Problematic
            IeB = IeUnique.*(seB<0) + IeBoth.*(HeB>=HeF);
            IeF = IeUnique.*(seF>0) + IeBoth.*(HeF>=HeB);
            Ie0 = IeEither;
            ce  = IeF.*ceF + IeB.*ceB + Ie0.*ce0;
            xie = wg*(1-tau) + (r-del).*aagrid - ce;
            Ae  = MatrixA(seB,seF,IeB,IeF);

            % 2.1.2. The Unemployed
            WuB(2:I,   :) = (Wu(2:I,  :)-Wu(1:I-1, :)) ./ (aagrid(2:I,:) - aagrid(1:I-1,:));
            WuF(1:I-1, :) = (Wu(2:I,  :)-Wu(1:I-1, :)) ./ (aagrid(2:I,:) - aagrid(1:I-1,:));
            WuB(1, :) = (hp + h(1,:) + (r-del)*amin).^(-gam); 
            WuF(I, :) = (hp + h(I,:) + (r-del)*amax).^(-gam); 
            cuB = max(WuB, 1e-6).^(-1/gam);
            cuF = max(WuF, 1e-6).^(-1/gam);
            cu0 = hp + h + (r-del).*aagrid;
            suB = hp + h + (r-del)*aagrid - cuB;
            suF = hp + h + (r-del)*aagrid - cuF;
            HuF = log(cuF) + WuF .* suF;
            HuB = log(cuB) + WuB .* suB;
            IuEither = (1-(suF>0)) .* (1-(suB<0));
            IuUnique = (suB<0).*(1-(suF>0)) + (1-(suB<0)).*(suF>0);
            IuBoth   = (suB<0).*(suF>0);
            IuB = IuUnique.*(suB<0) + IuBoth.*(HuB>=HuF);
            IuF = IuUnique.*(suF>0) + IuBoth.*(HuF>=HuB);
            Iu0 = IuEither;
            cu = IuF.*cuF + IuB.*cuB + Iu0.*cu0;
            xiu = hp + h + (r-del).*aagrid - cu;
            Au = MatrixA(suB,suF,IuB,IuF);
            
            c = [ce(:); cu(:)];
            if gam == 1
                utility = log(c);
            else
                utility = (c.^(1-gam)) / (1-gam);
            end

            % 2.1.4. Creating the new matrix Lambda
            Lam = MatrixLam(lamu);
            A = blkdiag(sparse(Ae), sparse(Au)) + C + Lam;

            % 2.1.5. Solving the Linear Algebra
            Dn = (1/Delta + rho)*speye(I*J*2) - A;
            Waz = [We(:);Wu(:)];
            dn = utility + Waz / Delta;
            WazNEW = Dn\dn;
            WeNEW = reshape(WazNEW(1:(I*J)), I, J);
            WuNEW = reshape(WazNEW((I*J)+1:end), I, J);
        
            distHJB = max(abs(WazNEW-Waz));
        
            if distHJB < tolHJB
                disp(['Value Function "Consumer HJB" Converged, Iteration = ', num2str(innerloop)]);
                break
            end
            % Update
            We = WeNEW; 
            Wu = WuNEW;
        end % Consumption innerloop ends 

        % 2.2. Solving Consumer's Fokker-Planck Equation    

            % 2.2.1. Setting up (fixing one value)
            da_tilde = 0.5*(daB+daF);
            da_tilde(1,:) = 0.5*daF(1,:);
            da_tilde(I,:) = 0.5*daB(I,:);
            da_stacked = repmat(da_tilde,2*J,1);
            grid_diag = dz*spdiags(da_stacked,0,I*J*2,I*J*2);
            M = I*J*2;
            fix = 1;
            AT = A';
            null = zeros(M,1);
            null(fix) = 1;
            AT(fix,:) = [zeros(1,fix-1),1,zeros(1, M-fix)];
            gVEC = AT\null;
            gsum = gVEC' * repmat(da_tilde,2*J,1)*dz;
            g = gVEC./ gsum;
            ge = g(1:I*J);
            gu = g(I*J+1:end);

        % 2.3. Solving Firm's HJB Equation
        wgVEC = wg(:);
        pi = repelem(zgrid,I)*(k^alf) - r*k - wgVEC;
        En = (1/Delta + lame + r - del)*speye(I*J) - Ae - Ce;
        for i = 1:n
            en = pi + Jf/Delta;
            Jfnew = En\en;
            distfirmHJB = max(abs(Jfnew - Jf));
            if distfirmHJB < tolHJB
                break
            end
            Jf = Jfnew;
        end

        % 2.5. Wage Bargaining
        JfIJ = reshape(Jf,I,J);
        dW = ce.^(-gam);
        ue = log(ce);
        dJf = [diff(JfIJ,1,1)./ diff(aagrid,1,1); zeros(1,J)];
        term1 = (zzgrid*(k^alf)  - r*k + dJf.*((r-del)*aagrid - ce)) ./ (1-dJf);
        term2 = (ue + dW.* ((r-del)*aagrid - ce) - rho*Wu) ./ dW;
        wgNEW = bet*term1 - (1-bet)*term2;
        Y = reshape(repelem(zgrid,I)*(k^alf) - r*k,I,J);
        wgNEW = min(wgNEW,Y);
        wgNEW = max((h+hp)/(1-tau),wgNEW);
        if any(wgNEW(:) < 0)
            disp('Non-positive number found in wgNEW.');
        end
        clear dW ue dJf term1 term2 Y;

        % 2.5. Market Clearing

            % 2.5.1. Free Entry Condition
            VV = -phi + (lamu/Theta) * Jf' * ((gu/u) .* repmat(da_tilde,J,1)*dz); 

            % 2.5.2. Reporting the Equilibirum
            revenu = sum(sum(tau*wg.*reshape(ge,I,J).*repmat(da_tilde,1,J)*dz));
            expend = sum(sum(     h.*reshape(gu,I,J).*repmat(da_tilde,1,J)*dz));

            % 2.5.3. Asset Market Clearing
                % 2.5.3.1. Asset Demand
                azazgrid = repmat(agrid,J,1);
                aggregatesavings = azazgrid' * ((ge+gu).*repmat(da_tilde,J,1)*dz);
                kas = aggregatesavings/(1-u);
                kgap = k - kas;

        if abs(VV) < tolENT && abs(kgap) < tolMKT
            disp(['Found the equilibrium! Theta = ', num2str(Theta), 'r = ', num2str(r)])
            break
        else
            disp(['OL:',num2str(outerloop),'VV:',num2str(VV),'kgap:',num2str(kgap),"tau:",num2str(tau),"maxrrate:",num2str(max(max(h./wg)))])   
        end 

        % 2.7. Updating Parameters
            % 2.7.1. Labor Market
            Theta = Theta + dTheta * VV;
            lamu = chi * (Theta ^ (1 - eta));
            q = chi * (Theta^(-eta));
            u = lame / (lame + lamu);
            v = Theta * u;
            
            % 2.7.2. Capital Market and Prices
            k = relax*k +(1-relax)*kas;
            z = KFE4z(mue,muu,lamu,lame);
            r = z*alf*k^(alf-1);
            wg = wgNEW; 
            clear wgNEW;
            
            % 2.7.3. Unemployment Insurance 
            govrev = sum(sum((tau*wg).*reshape(ge,I,J).*repmat(da_tilde,1,J)*dz));
            h = govrev/(sum(sum(reshape(gu,I,J).*repmat(da_tilde,1,J)*dz)));
            h = h*ones(I,J);
            disp(['maxnrr:',num2str(max(max((hp+h)./wg*(1-tau))))]);
    end

    Ge = reshape(ge,I,J) .* repmat(da_tilde,1,J) * dz; clear ge;
    Gu = reshape(gu,I,J) .* repmat(da_tilde,1,J) * dz; clear gu;
    
    avrwg = sum(sum(wg.*Ge))/(sum(sum(Ge)));
    maxnrr = max(max((h./(wg*(1-tau)))));
    minnrr = min(min((h./(wg*(1-tau)))));
    anrr = mean(mean(h))/(avrwg*(1-tau));
    avrast = sum(sum(aagrid.*(Ge+Gu)));
    avrcon = sum(sum(ce.*Ge + cu.*Gu));
    avrcone = sum(sum(ce.*Ge))/(sum(sum(Ge)));
    avrconu = sum(sum(cu.*Gu))/(sum(sum(Gu)));
    
    VA = sum(sum(((aagrid-avrast).^2).*(Ge+Gu))); 
    VC = sum(sum(((ce-avrcon).^2).*Ge + ((cu-avrcon).^2).*Gu));
    VCe = sum(sum(((ce-avrcone).^2).*(Ge/(sum(sum(Ge))))));
    VCu = sum(sum(((cu-avrconu).^2).*(Gu/(sum(sum(Gu))))));

    EJf = sum(sum(JfIJ .* Gu))/(sum(sum(Gu)));
    ur = lame/(lame+lamu);
    zemean = KFE4z(mue, muu, lamu, lame);
    zumean = sum(sum(zzgrid.*Gu))/(sum(sum(Gu)));
    
    taupct = num2str(round(100*tau,2));
    
% Plot 1: Distribution - Employed
    figure('Position', [100, 100, 800, 600]);
    surf(aagrid(1:100,:),zzgrid(1:100,:),Ge(1:100,:),'Edgealpha',.25);
    xlim([agrid(1), agrid(100)]);
    xlabel('Asset level (a)','FontSize',16);
    ylabel('Productivity (z)','FontSize',16);
    zlabel('Density','FontSize',16);
    zlim([0 0.05]);
    zticks(0:0.005:0.035);
    colormap(cool);
    filenamepng = sprintf('Figures/Surfplot - Distribution - Employed - %s.png', taupct);
    exportgraphics(gca, filenamepng, 'Resolution', 300);

% Plot 2: Distribution - Unemployed
    figure('Position', [100, 100, 800, 600]);
    surf(aagrid(1:100,:),zzgrid(1:100,:),Gu(1:100,:),'Edgealpha',.25);
    xlim([agrid(1), agrid(100)]);
    xlabel('Asset level (a)','FontSize',16);
    ylabel('Productivity (z)','FontSize',16);
    zlabel('Density','FontSize',16);
    colormap(cool);
    zlim([0 0.00075]);
    zticks(0:0.00025:0.00075);
    filenamepng = sprintf('Figures/Surfplot - Distribution - Unemployed - %s.png', taupct);
    exportgraphics(gca, filenamepng, 'Resolution', 300);

% Plot 3: Productivity Distribution
    figure('Position', [100, 100, 800, 600]);
    plot(zgrid,sum(Ge,1)/(sum(sum(Ge))),'b-', 'LineWidth', 1);
    hold on;
    plot(zgrid,sum(Gu,1)/(sum(sum(Gu))),'r-', 'LineWidth', 1);
    xlabel("Productivity",'FontSize',16);
    ylabel("Density (in each group)",'FontSize',16);
    legend("Employed", "Unemployed");
    filenamepng = sprintf('Figures/Plot - zDistribution - %s.png', taupct);
    exportgraphics(gca, filenamepng, 'Resolution', 300);

% Plot 4: The Savings Difference - Productivity (15)
    figure('Position', [100, 100, 800, 600]);
    plot(agrid(1:100),xie(1:100, 15), 'b-', 'LineWidth', 1.5);
    hold on;
    plot(agrid(1:100),xiu(1:100, 15), 'r-', 'LineWidth', 1.5);
    yline(0, '--k', 'LineWidth', 1);
    xlabel("Asset",'FontSize',16);
    ylabel("Savings",'FontSize',16);
    legend("Employed", "Unemployed");
    xlim([agrid(1), agrid(100)]);
    ylim([-2.0 2.0]);             
    yticks(-2.0:0.25:2.0);         
    filenamepng = sprintf('Figures/Plot - Difference in Savings - 15 - %s.png', taupct);
    exportgraphics(gca, filenamepng, 'Resolution', 300);

% Plot 5: The Savings Difference - Productivity (12)
    figure('Position', [100, 100, 800, 600]);
    plot(agrid(1:100),xie(1:100,12), 'b-', 'LineWidth', 1.5);
    hold on;
    plot(agrid(1:100),xiu(1:100,12), 'r-', 'LineWidth', 1.5);
    yline(0, '--k', 'LineWidth', 1);
    xlabel("Asset",'FontSize',16);
    ylabel("Savings",'FontSize',16);
    legend("Employed", "Unemployed");
    xlim([agrid(1), agrid(100)]);
    ylim([-2.0 2.0]);             
    yticks(-2.0:0.25:2.0);         
    filenamepng = sprintf('Figures/Plot - Difference in Savings - 12 - %s.png', taupct);
    exportgraphics(gca, filenamepng, 'Resolution', 300);

% Plot 6: The Savings Difference - Productivity (9)
    figure('Position', [100, 100, 800, 600]);
    plot(agrid(1:100),xie(1:100,9), 'b-', 'LineWidth', 1.5);
    hold on;
    plot(agrid(1:100),xiu(1:100,9), 'r-', 'LineWidth', 1.5);
    yline(0, '--k', 'LineWidth', 1);
    xlabel("Asset",'FontSize',16);
    ylabel("Savings",'FontSize',16);
    legend("Employed", "Unemployed");
    xlim([agrid(1), agrid(100)]);
    ylim([-2.0 2.0]);             
    yticks(-2.0:0.25:2.0);         
    filenamepng = sprintf('Figures/Plot - Difference in Savings - 09 - %s.png', taupct);
    exportgraphics(gca, filenamepng, 'Resolution', 300);

% Plot 7: The Savings Difference - Productivity (6)
    figure('Position', [100, 100, 800, 600]);
    plot(agrid(1:100),xie(1:100,6), 'b-', 'LineWidth', 1.5);
    hold on;
    plot(agrid(1:100),xiu(1:100,6), 'r-', 'LineWidth', 1.5);
    yline(0, '--k', 'LineWidth', 1);
    xlabel("Asset",'FontSize',16);
    ylabel("Savings",'FontSize',16);
    legend("Employed", "Unemployed");
    xlim([agrid(1), agrid(100)]);
    ylim([-2.0 2.0]);             
    yticks(-2.0:0.25:2.0);         
    filenamepng = sprintf('Figures/Plot - Difference in Savings - 06 - %s.png', taupct);
    exportgraphics(gca, filenamepng, 'Resolution', 300);

    % Plot 8: The Savings Difference - Productivity (3)
    figure('Position', [100, 100, 800, 600]);
    plot(agrid(1:100),xie(1:100,3), 'b-', 'LineWidth', 1.5);
    hold on;
    plot(agrid(1:100),xiu(1:100,3), 'r-', 'LineWidth', 1.5);
    yline(0, '--k', 'LineWidth', 1);
    xlabel("Asset",'FontSize',16);
    ylabel("Savings",'FontSize',16);
    legend("Employed", "Unemployed");
    xlim([agrid(1), agrid(100)]);
    ylim([-2.0 2.0]);             
    yticks(-2.0:0.25:2.0);
    filenamepng = sprintf('Figures/Plot - Difference in Savings - 03 - %s.png', taupct);
    exportgraphics(gca, filenamepng, 'Resolution', 300);

    % Plot 9: Job Values
    figure('Position', [100, 100, 800, 600]);
    surf(aagrid(1:100,:),zzgrid(1:100,:),JfIJ(1:100,:),'Edgealpha',.5,'FaceAlpha',.25);
    zlim([0 1.4]);   % Set z-axis limits
    zticks(0:0.2:1.4); % Set z-axis intervals
    xlabel('Asset level (a)','FontSize',16);
    ylabel('Productivity (z)','FontSize',16);
    zlabel('Job Value','FontSize',16);
    zlim([0 1.4]);
    zticks(0:0.35:1.4);    
    colormap(winter);
    filenamepng = sprintf('Figures/Plot - Job Value - %s.png', taupct);
    exportgraphics(gca, filenamepng, 'Resolution', 300);

    % Plot 10 Asset Distribution
    adistribution = sum(Ge+Gu,2);
    cumadistribution = cumsum(adistribution);
    fq = find(cumadistribution >= 0.25, 1, 'first');
    median = find(cumadistribution >= 0.5, 1, 'first');
    tq = find(cumadistribution >= 0.75, 1, 'first');
    figure('Position', [100, 100, 800, 600]);
    plot(agrid(1:100),adistribution(1:100), 'b-', 'LineWidth', 1.25);
    hold on;
    xline(agrid(median), '--k', 'LineWidth', 1.5);
    xline(agrid(fq), '--', 'Color','#4D4D4D', 'LineWidth', 1.5);
    xline(agrid(tq), '--', 'Color','#4D4D4D', 'LineWidth', 1.5);
    yline(0, '--k', 'LineWidth', 1);
    xlabel("Asset",'FontSize',16);
    ylabel("Density",'FontSize',16);
    xlim([agrid(1), agrid(100)]);
    ylim([0 0.035]);
    yticks(0:0.005:0.035);
    text(agrid(median), 0.0335, sprintf('%.2f',agrid(median)),'HorizontalAlignment','center','VerticalAlignment','top','FontSize',12,'Color','k','BackgroundColor', 'w');
    text(agrid(fq), 0.0335, sprintf('%.2f',agrid(fq)),'HorizontalAlignment','center','VerticalAlignment','top','FontSize',12,'Color','#4D4D4D','BackgroundColor', 'w');
    text(agrid(tq), 0.0335, sprintf('%.2f',agrid(tq)),'HorizontalAlignment','center','VerticalAlignment','top','FontSize',12,'Color','#4D4D4D','BackgroundColor', 'w');
    filenamepng = sprintf('Figures/Plot - Asset Distribution - %s.png', taupct);
    exportgraphics(gca, filenamepng, 'Resolution', 300);

% Plot 14 Consumption Difference
    figure('Position', [100, 100, 800, 600]);
    plot(agrid(1:100),100*((cu(1:100,15)-ce(1:100,15))./ce(1:100,15)), 'b-', 'LineWidth', 1);
    hold on;
    plot(agrid(1:100),100*((cu(1:100,3)-ce(1:100, 3))./ce(1:100, 3)), 'r-', 'LineWidth', 1);
    yline(0, '--k', 'LineWidth', 1);
    xlabel("Asset",'FontSize',16);
    ylabel("Percent (%)",'FontSize',16);
    legend("z = 1.0", "z = .73");
    xlim([agrid(1), agrid(100)]);
    ylim([-100, 25]);
    yticks(-100:25:25);
    filenamepng = sprintf('Figures/Plot - Difference in Consumption - %s.png', taupct);
    exportgraphics(gca, filenamepng, 'Resolution', 300);
end

toc;


