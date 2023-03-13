classdef Hutchinson
    properties
        dz              % Zenner pinning grain size
        dc              % Critical grain size
        T               % Current temperature
        t               % Current time
        strain          % 
        dhb             % Hot band grain size
        dgammagb        % Diameter of gamma fibre grains nucleated at grain boundaries
        dgammain        % Diameter of gamma fibre grains nuclearted intra-grain
        dhkl            % Diameter of grandom grains
        alld            % All grain diameters collate three types of grain together
        dmin
        dmax
        dmean           % Average grain diameter of all grains
        Ahkl            % Constant for nucleation rate for random orientation grains
        Agammain        % Constant for nucleation rate for gamma fibre grains
        B               % constant for nucleation rate for gamma grain on gb grate
        C               % Growth constant
        Q=370e3           % activation energy kJ/mol
        x               % Fraction recrystalised
        xhistory        % History of x as a function of time
        timehistory     % time steps 
        ctimehistory    % Coiling time history
        X0              % Incubation period, dimensionless
        TcoilCoolingStart      % Coilding temperature
        TcoilCoolingEnd 
        CoilCoolingRate     % Coiling colling rate K/s
        heatingRate     % Heating rate K/s
        Tstart          % 
        Thold           % Isothermal holding temperature for annealing
        Thistory
        timehold        % Isothermal holding time
        tpstart
        tend            % time to end
        wt_Ti           % Titanium content
        wt_S            % Sulphur content
        wt_N            % N content
        wt_C            % Carbon content wt%
        r_tic           % 
        r0              % r0 when cooled down to 500 degC after coiling
        xf = 0.96       % Fraction of recrystallisation considered to be complete
        dt=0.1         % time step
        stepcounter     % counter of time step
        volume          % volume of the materials, set to be volume of a hot band grain
        esurfTic        % Surface energy of TiC precipitate
        dzhistory       % Zenner pinning size for each time step
        k               % Constant for Hillert growth equation
        tptimed = false
        tcoil
        dz0history      %
        Vm = 12.1481e-6 % % Molar volume of TiC [m^3/mol]
        
    end
    
    properties(Constant)
         R= 8.31446261815324            % Avagadron gas constant (J⋅K^{-1}⋅mol^{-1})
%          kb = 1.38064852e-23             % Boltzman constant (m2 kg s-2 K-1)
%          k = 1e-5                        % m^2/s
    end

    
    methods
        function f = Hutchinson
        end
        
        function value = get.volume(hm)
            value = 1;
%             value = (pi*hm.dhb^3)/6;
        end
        
        function value = get.tend(hm)
            value = (hm.Thold-hm.Tstart)/hm.heatingRate + hm.timehold;
        end
        
        function f = temperature(hm,t)
            % Obtain temperature given a time, t
            % t must be between 0 and the end of isothermal holding
            if nargin < 2
                t = hm.t;
            end
            t_heating = (hm.Thold - hm.Tstart)/ hm.heatingRate;
            if t <= t_heating
                f = hm.Tstart + t * hm.heatingRate;
            elseif t > t_heating && t <= t_heating + hm.timehold
                f = hm.Thold;
            elseif t > t_heating + hm.timehold
                warning('Time is out of range. I do not record temperature after the isothermal holding.');
                f = [];
            end
        end
        
        function hm = coil(hm)
            %
            hm.tcoil= (hm.TcoilCoolingStart - hm.TcoilCoolingEnd)/hm.CoilCoolingRate;
            hm.r_tic = 0;
            for time = 0:hm.dt:hm.tcoil
                currentT = hm.TcoilCoolingStart - time*hm.CoilCoolingRate;
                hm = hm.precipitate(hm.r_tic, hm.dt, currentT);
            end
        end
        
%         function value = get.timehistory(hm)
%             totaltime = hm.Thold/hm.heatingRate + hm.timehold;
%             value = transpose(0:hm.dt:totaltime);
%         end

        
        function hm = anneal(hm)
            % 
            hm.x = 0;
            hm.t = 0;
            hm.dzhistory = zeros(1,1);
            hm.dgammagb = zeros(1,2);
            hm.dgammain = zeros(1,2);
            hm.dhkl = zeros(1,2);
            hm.stepcounter = 0;
            ph = waitbar(0, 'Annealing...');
            while hm.x < hm.xf && hm.t < hm.tend
                progress = hm.x / hm.xf;
                msg = ['Time: ', num2str(hm.t,'%05.2f'), ' s; ', 'x = ', num2str(hm.x, '%6.3e')];
                waitbar(progress, ph, msg);
                hm.stepcounter = hm.stepcounter + 1;
                hm.T = hm.temperature(hm.t);
                hm.timehistory(hm.stepcounter) = hm.t;
                hm.Thistory(hm.stepcounter) = hm.T;
                if hm.x >= hm.X0
                    hm = hm.nucleate('matrix');
                end
                hm = hm.nucleate('gb');
                if hm.T >=773 
                    if ~hm.tptimed
                        hm.tpstart = hm.t;
                        hm.tptimed = true;
                        r0anneal = hm.r_tic;
                    end
                    hm = hm.precipitate(r0anneal, hm.t-hm.tpstart, hm.T);
                    hm.dzhistory(hm.stepcounter) = hm.dz;
                end

                hm.dmean = hm.calcdmean;
                if hm.dmean <= hm.dz
                    hm = hm.grow;
                end
                hm.xhistory(hm.stepcounter) = hm.x;
                hm.t = hm.t + hm.dt;
            end
            delete(ph)
        end
        
        function value = get.dmin(hm)
            value = min([min(hm.dgammagb(:,2)), min(hm.dgammain(:,2)), min(hm.dhkl(:,2))]);
        end
        
        function value = get.dmax(hm)
            value = max([max(hm.dgammagb(:,2)), max(hm.dgammain(:,2)), max(hm.dhkl(:,2))]);
        end
        
        function value = get.alld(hm)
            value = vertcat(hm.dgammagb, hm.dgammain, hm.dhkl);
        end
        
        function hm = grow(hm)
            % Calculate grain growth
            
            growthRate = hm.C * hm.strain * (exp(-hm.Q / (hm.R * hm.T))) * (1 - hm.x);
            hm.dgammagb(:,2) = min(hm.dgammagb(:,2) + growthRate * hm.dt, hm.dz);
            hm.dgammain(:,2) = min(hm.dgammain(:,2) + growthRate * hm.dt, hm.dz);
            hm.dhkl(:,2) = min(hm.dhkl(:,2) + growthRate * hm.dt, hm.dz);
            doHillert = true;
            if doHillert
                v0 = sum(((pi*(hm.alld(:,2).^3))/6) .* hm.alld(:,1));
                options = optimset('TolX', 1e-18, 'MaxFunEvals', 50000, 'TolFun', 1e-25);
                dcfound = fminbnd(@dcObjFun, 0,hm.dmax, options);

    %             dcfound = fminsearch(@dcObjFun, 1e-9);
                hm.dc = dcfound;
                
                hm.dgammagb(hm.dgammagb(:,1)>0,2) = hm.dgammagb(hm.dgammagb(:,1)>0,2) + hm.hillertGrow(hm.dgammagb(hm.dgammagb(:,1)>0,2)) * hm.dt;
                hm.dgammain(hm.dgammain(:,1)>0,2) = hm.dgammain(hm.dgammain(:,1)>0,2)+ hm.hillertGrow(hm.dgammain(hm.dgammain(:,1)>0,2)) * hm.dt;
                hm.dhkl(hm.dhkl(:,1)>0,2) = hm.dhkl(hm.dhkl(:,1)>0,2) + hm.hillertGrow(hm.dhkl(hm.dhkl(:,1)>0,2)) * hm.dt;
                hm.dgammagb(hm.dgammagb(:,2)<=0,:) = [];
                hm.dgammain(hm.dgammain(:,2)<=0,:) = [];
                hm.dhkl(hm.dhkl(:,2)<=0,:) = [];
            end
            hm.x = sum(hm.alld(:,1).* ((4*pi*(hm.alld(:,2)).^3 )/3))/hm.volume;
            hm.x = min(hm.x, 1);
            
            function f = dcObjFun(x)
                hm.dc = x;
                hillertrate = hm.hillertGrow(hm.alld(:,2));
                d = hm.alld(:,2) + (hillertrate*hm.dt);
                d(d<=0) = 0;
                v = sum(((pi*d.^3)/6).*hm.alld(:,1));
%                 if ~isscalar(v)
%                     warning('v is not a scalar, check it out ...')
%                 end
                f = (v-v0).^2;
            end
        end
        
        function f = hillertGrow(hm, dgrains)                
                f = hm.k * (hm.dc^(-1) - dgrains.^(-1));
        end
        
        function f = calcdmean(hm)
            % Calculate the average grain size of all grain
            numgrains = sum(hm.dgammagb(:,1))+ sum(hm.dgammain(:,1)) + sum(hm.dhkl(:,1));
            if numgrains
                f = (sum(hm.dgammagb(:,2) .* hm.dgammagb(:,1)) + sum(hm.dgammain(:,2) .* hm.dgammain(:,1)) + sum(hm.dhkl(:,1) .* hm.dhkl(:,2))) /numgrains;
            else
                f = 0;
            end
        end
        
        function hm = precipitate(hm, r0, tp, temperature)
            % Calculate precipitate size and the Zenner pinning limit grain
            % size
            
            Tis = (10^(4.45-(10800/temperature)))/hm.wt_C;
            D = (1.42e-4)*exp(-2.32e5/(hm.R*temperature)); % Diffusivity
%             hm.esurfTic = 9.08; % [J/m^2]
            hm.r_tic = (((8/9)*hm.esurfTic * D * Tis * hm.Vm)*tp/(hm.R*temperature) + r0^3)^(1/3);
            atomNum_Ti  = 22;
            atomNum_S = 16; 
            atomNum_N = 7; 
            atomNum_C = 6;
            rho_tic = 4.93e3; % [kg/m^3]
            rho_steel = 7.8e3; 
            wt_Ti_TiC = hm.wt_Ti - Tis- hm.wt_N*atomNum_Ti/atomNum_N - hm.wt_S * atomNum_Ti / atomNum_S;
            vf = 0.01*wt_Ti_TiC * (atomNum_Ti + atomNum_C)*rho_steel/(rho_tic*atomNum_Ti);
            hm.dz = 1/((0.5*vf /hm.r_tic)+1/30);
        end
        
        function hm = nucleate(hm,where)
            %
            switch where
                case {'matrix', 'intragrain', 'ingrain', '1'}
                    Ndothkl = hm.Ahkl * hm.strain.^4 * exp(-hm.Q./(hm.R*hm.T))*(1-hm.x);
                    hm.dhkl(hm.stepcounter,1) = Ndothkl * hm.dt * hm.volume;
                    Ndotgammam = hm.Agammain * hm.strain.^4 * exp(-hm.Q./(hm.R*hm.T))*(1-hm.x);
                    hm.dgammain(hm.stepcounter,1) = Ndotgammam * hm.dt * hm.volume;
                case {'gb', 'GB', 'grain boundaries', 'GBs', 'gbs', '2'}
                    Ndotgammagb = hm.B * (hm.strain^4) * exp(-hm.Q./(hm.R*hm.T))*(1-hm.x)*...
                        ((2/hm.dhb)*(1+0.13*hm.strain+0.235*hm.strain^2 + 0.2*hm.strain.^3));
                    hm.dgammagb(hm.stepcounter,1) = Ndotgammagb * hm.dt * hm.volume;
            end
        end
        
        function plotx(hm)
        end
        
    end       
end