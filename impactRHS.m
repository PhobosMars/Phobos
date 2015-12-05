function ret=impactRHS(t,y,R,g,Clift,A,drag,Vdark,Meltenthalpy,heatcoeff,L,visc,spec_heat_ratio);

    % set up return column vector
    ret = zeros(5,1);
    % get current values 
    V = y(1);
    m = y(2);
    theta = y(3);
    h = y(4);
    X = y(5);
    
    % Exponentially decaying atmosphere, fit manually to Schofield et al., 1997
    rhoatm=(1.76e-2)*10^(-1*h/1000*9/160);

    Tg=150; % Gas temperature in free field (K); see Clancy et al., 2000, JGR-Planets: "An intercomparison of ground?based millimeter, MGS TES, and Viking atmospheric temperature measurements"
    Tp=150; % Temperature of entering particle (K).
    
    csound=(spec_heat_ratio*8.314*Tg/0.044)^0.5;
    Re=L*V/visc;
    Mach=V/csound;
    Cinc=24/Re*(1+0.15*Re^.687);
    GRe=10^((2.5*(Re/312)^.6688)/(1+(Re/312)^0.6688)); % G(Re) Melosh and Goldin, LPSC Abstract 2457 2008
    HofM=(4.6/(1+Mach))+1.7*(Tp/Tg)^0.5; % H(M) Melosh and Goldin, LPSC Abstract 2457 2008
    
    % Drag coefficient: Melosh and Goldin, LPSC Abstract 2457 2008
    Cdrag=2+(Cinc-2)*exp(-3.07*((spec_heat_ratio)^.5)*Mach*GRe/Re)+(HofM*exp(-1*Re/2*Mach)/(Mach*(spec_heat_ratio)^.5));

    % return updates -- See Appendix 1, Beech 2013 Earth, Moon, Planets
    ret(1) = -1*drag*rhoatm*(A/m)*V^2 + g*sin(theta);
    if (V>=Vdark) % Must turn off ablation once V drops below dark-flight velocity
        ret(2) = -1*(heatcoeff/(2*Meltenthalpy))*rhoatm*A*V^3*(V^2-Vdark^2)/(V^2);
    else
        ret(2) = 0;
    end
    ret(3) = (1/(m*V))*(m*g*cos(theta)-0.5*Clift*rhoatm*A*V^2)-(V*cos(theta)/(R+h));
    ret(4) = -1*V*sin(theta);
    ret(5) = V*cos(theta)/(1+h/R);