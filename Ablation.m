% Implementation of the atmospheric dispersion model in:

% Beech, Martin. "Towards an Understanding of the Fall 
% Circumstances of the Hoba Meteorite." Earth, Moon, 
% and Planets 111.1-2 (2013): 15-30.

% (c) 2015, Benjamin Black and Tushar Mittal. 
% Please email with questions or comments: 
% bblack@berkeley.edu,
% tmittal2@berkeley.edu


% Parameters
R       = 3390e3;                               % m, Mars radius
g       = 3.71;                                 % Mars g
Clift   = 0.0;                                  % Lift coefficient, Beech 2013

visc    = 0.75e-5;                              % m2/s, kinematic viscosity, Beech 2013
spec_heat_ratio=1.289;                          % Specific heat ratio for CO2, Beech 2013
rhoPh   = 1873;                                 % kg/m3, density of Phobos, Murchie and Erard 1996
drag    = 1.2 ;                                 % drag coefficient, Beech 2013
Vdark   = 2e3;                                  % m/s, Beech 2013
Meltenthalpy= 5e6;                              % J/kg, enthalpy of melting, Beech 2013
heatcoeff=0.01;                                 % heat coefficient, Beech 2013

r_mars = 3390e3 ;                               % meters, radius of Mars
r_phobos = 11.267e3 ;                           % meters, radius of Phobos
G=6.673e-11;                                    % Gravitational constant
M_mars=6.4185e23;                               % kg, mass of Mars
M_phobos=1.0659e16;                             % kg, mass of Phobos
 
% Initialization

tRange = [0,1];              % 1 second timestep
yZero = zeros(5,1);

figure
numsizes=9;                  % set range of bolides under consideration
xmin=0;

% Perform calculations for each bolide size

for ii=2:numsizes
    
    L=2^(ii+2);               % m, diameter of bolide
    A= 2*3.1415*L^2;          % m^2, area of projectile face
    
  
    yZero(2) = rhoPh*4/3*3.1415*(L/2)^3;        % kg, mass of bolide
    yZero(3) = pi/72;         % theta (radians, pi/72 is 2.5 degrees). This is the entry angle.   
    yZero(4) = 160e3;         % h
    yZero(5) = 0;             % X
    yZero(1) = (G*(M_mars+yZero(2))/(r_mars+yZero(4)))^.5;    % V, m/s 
    
    iii=2;
    t=0;
    h_t=yZero(4);
    x_t=yZero(5);
    m0=yZero(2);
    m_t=m0;
    inair=1;
    
    while (inair==1)
        
        % At each timestep, solve coupled equations with MATLAB ODE solver
        % to maintain stability.
        [myT,myY]=ode23s(@(myT,myY)impactRHS(myT,myY,R,g,Clift,A,drag,Vdark,Meltenthalpy,heatcoeff,L,visc,spec_heat_ratio),tRange,yZero);
        yZero=myY(end,:);
        
        % If mass or altitude drops to zero, terminate calculations
        if yZero(:,4)<=0
            inair=0;
        end
        if yZero(:,2)<=0
            inair=0;
        end
        
        % Save values
        t(iii)=t(iii-1)+myT(end);
        h_t(iii)=yZero(end,4);
        x_t(iii)=yZero(end,5);
        m_t(iii)=yZero(end,2);
        iii=iii+1;
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot dispersion, crater sizes, and mass loss due to ablation.
    
    subplot(3, 1, 1);
    plot(x_t/1e3,h_t/1e3,'Color',[(0.4*ii/numsizes) (0.4*ii/numsizes) (0.75*ii/numsizes)],'LineWidth',1)
    hold on
    crat_r= 1e-3*0.5*1.161*(g^-.22)*(yZero(end,1)^0.44)*(L^0.78)*(sin(yZero(end,3)))^(1/3);% Collins 2011 oblique impacts A.7
    ang=0:0.01:2*pi; 
    %xp=crat_r*13.2*cos(ang);
    xp=crat_r*cos(ang);
    yp=crat_r*sin(ang);
    x_crat=interp1(h_t,x_t,0);
    if xmin>0
        xmin=min([min(x_crat)/1e3 xmin]);
    else
        xmin=min(x_crat)/1e3;
    end
    
    xmax=2580;
    ylim([0 5])
    if (xmax>0)
        xmin=2500;
        xlim([xmin xmax])
    else
        xlim([xmin max(x_crat/1e3+xp)])
    end
    subplot(3, 1, 2);
    plot(x_crat/1e3+xp,0+yp,'Color',[(0.4*ii/numsizes) (0.4*ii/numsizes) (0.75*ii/numsizes)],'LineWidth',1);
    ylim([0 2.5])
    if (xmax>0)
        xmin=2500;
        xlim([xmin xmax])
    else
        xlim([xmin max(x_crat/1e3+xp)])
    end
    hold on
    subplot(3, 1, 3);
    plot(x_t/1e3,m_t/m0,'Color',[(0.4*ii/numsizes) (0.4*ii/numsizes) (0.75*ii/numsizes)],'LineWidth',1);
    ylim([.99 1])
    if (xmax>0)
        xmin=2500;
        xlim([xmin xmax])
    else
        xlim([xmin max(x_crat/1e3+xp)])
    end
    hold on
end


