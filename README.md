# DAVS-block-
DAVS block MATLAB script
%% Constants
F = 96485; % Faraday's constant (C/mol)
k = 1; % unit conversion constant (Vm3/C)
dt = 1; % time step for simulation
I_app_pe = 15;
%% Evaluated Active Material Particle Parameters
dR_pe = params.R_pe/10; % width of each layer
S_pe = 4*pi*(params.R_pe*(1:10)/10).^2; % Surface area of each layer (m2)
dV_pe = (4/3)*pi*((params.R_pe*(1:10)/10).^3 - (params.R_pe*(0:10-1)/10).^3); % Volume of each layer (m3)
Rd_pe = (k*10)./(4*pi*params.R_pe*params.D_pe*(1:10).^2); % diffusion resistance between the layers (Ohms)
N_a_pe = (3*params.Epsilon_pe*params.V_ed_pe)/(4*pi*params.R_pe.^3); % Number of active particles
Ibar0 = 0; % initial current input (A)
Ibar_10_pe = -I_app_pe/N_a_pe; % Input current to each particle (A) 
%% Simulation Parameters
t_pe = [Ibar_10_pe*ones(1,1740), Ibar_10_pe*zeros(1,3600), -Ibar_10_pe*ones(1,1020), -Ibar_10_pe*zeros(1,3624)]; % discharge and rest, charge and rest
timestep = length(t_pe);
cse_pe(1) = params.c0_pe;
c_pe(1, :) = params.c0_pe*ones(1,10);


% Calculating surface concentration due to charge/discharge


for j = 1:timestep-1

    % Evaluate Vn
    for n = 1:10
    V_pe(j,n) = (k * F) .* c_pe(j, n); % Voltage vector at each time step
    end

    % Evaluate Ibar
    for n = 1:10
    if n == 10
        Ibar_pe(j, n) = t_pe(j); % Ibar for end shell
    else    
        Ibar_pe(j, n) = (V_pe(j,n)-V_pe(j,n+1))./Rd_pe(n); % Ibar vector at each time step
    end
    end

    % Evaluate I
    for n = 1:10
    if n == 1
        I_pe(j,n) = -Ibar_pe(j,n); % I current for first layer
    else
        I_pe(j,n) = Ibar_pe(j,n-1) - Ibar_pe(j,n); % I current for intermediate layers
    end
    end

    % Evaluate concentration values for all layers
    for n = 1:10
    c_pe(j+1,n) = c_pe(j,n) + (I_pe(j,n)*(dt./(F*dV_pe(n)))); % concentration values for all layers
    end

    cse_pe(j+1) = 1.5*c_pe(j+1,end) - 0.5*c_pe(j+1,end-1);
end
sto_pe = cse_pe./params.c_max_pe;

OCP_pe = -(0.8090*sto_pe) + 4.4875 - (0.0428*tanh(18.5138*(sto_pe - 0.5542))) - (17.7326*tanh(15.789*(sto_pe - 0.3117))) + (17.5842*tanh(15.9308*(sto_pe - 0.3120)));

figure(1);
plot((1:length(t_pe))/60,t_pe);
grid on;
title('Current input for positive electrode particle')
xlabel('Time (min)'); ylabel('Curent (A)')

figure(2); 
plot((1:length(t_pe))/60,cse_pe/1000);
grid on;
title('Surface concentration of positive electrode particle')
xlabel('Time (min)'); ylabel('Concentration (kmol m^{-3})')

figure(3);
plot((1:length(t_pe))/60,OCP_pe);
grid on;
title('Open circuit potential for positive electrode')
xlabel('Time (min)'); ylabel('Voltage (V)')

%% Evaluated Active Material Particle Parameters
I_app_ne = 10;
dR_ne = params.R_ne/10; % width of each layer
S_ne = 4*pi*(params.R_ne*(1:10)/10).^2; % Surface area of each layer (m2)
dV_ne = (4/3)*pi*((params.R_ne*(1:10)/10).^3 - (params.R_ne*(0:10-1)/10).^3); % Volume of each layer (m3)
Rd_ne = (k*10)./(4*pi*params.R_ne*params.D_ne*(1:10).^2); % diffusion resistance between the layers (Ohms)
N_a_ne = (3*params.Epsilon_ne*params.V_ed_ne)/(4*pi*params.R_ne.^3); % Number of active particles
Ibar0 = 0; % initial current input (A)
Ibar_10_ne = I_app_ne/N_a_ne; % Input current to each particle (A)

%% Simulation Parameters
t_ne = [Ibar_10_ne*ones(1,1740), Ibar_10_ne*zeros(1,3600), -Ibar_10_ne*ones(1,1020), -Ibar_10_ne*zeros(1,3624)]; % discharge and rest, charge and rest
timestep = length(t_ne);
cse_ne(1) = params.c0_ne;
c_ne(1, :) = params.c0_ne*ones(1,10);


for j = 1:timestep-1

    % Evaluate Vn
    for n = 1:10
    V_ne(j,n) = (k * F) .* c_ne(j, n); % Voltage vector at each time step
    end

    % Evaluate Ibar
    for n = 1:10
    if n == 10
        Ibar_ne(j, n) = t_ne(j); % Ibar for end shell
    else    
        Ibar_ne(j, n) = (V_ne(j,n)-V_ne(j,n+1))./Rd_ne(n); % Ibar vector at each time step
    end
    end

    % Evaluate I
    for n = 1:10
    if n == 1
        I_ne(j,n) = -Ibar_ne(j,n); % I current for first layer
    else
        I_ne(j,n) = Ibar_ne(j,n-1) - Ibar_ne(j,n); % I current for intermediate layers
    end
    end

    % Evaluate concentration values for all layers
    for n = 1:10
    c_ne(j+1,n) = c_ne(j,n) + (I_ne(j,n)*(dt./(F*dV_ne(n)))); % concentration values for all layers
    end

    cse_ne(j+1) = 1.5*c_ne(j+1,end) - 0.5*c_ne(j+1,end-1);
end
sto_ne = cse_ne./params.c_max_ne;

OCP_ne = (1.9793*exp(-39.3631*sto_ne)) + 0.2482 - (0.0909*tanh(29.8538*(sto_ne - 0.1234))) - (0.04478*tanh(14.9159*(sto_ne - 0.2769))) - (0.0205*tanh(30.4444*(sto_ne - 0.6103)));


figure(4);
plot((1:length(t_ne))/60,t_ne);
grid on;
title('Current input for negative electrode particle')
xlabel('Time (min)'); ylabel('Curent (A)')

figure(5); 
plot((1:length(t_ne))/60,cse_ne/1000);
grid on;
title('Surface concentration of negative electrode particle')
xlabel('Time (min)'); ylabel('Concentration (kmol m^{-3})')

figure(6);
plot((1:length(t_ne))/60,OCP_ne);
grid on;
title('Open circuit potential for negative electrode')
xlabel('Time (min)'); ylabel('Voltage (V)')
%% OCV evaluation of LG M50 cell
OCV = OCP_pe - OCP_ne;

figure(7);
plot((1:length(t_ne))/60,OCV);
grid on;
title('Open circuit voltage for cell')
xlabel('Time (min)'); ylabel('Voltage (V)')
