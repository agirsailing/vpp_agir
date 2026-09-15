function env = get_env_params()

    %% Constants
    env.g = 9.81;                           % m/s^2
    env.kn2ms = 0.514444;
    env.ms2kn = 1/env.kn2ms;

    %% Water (fresh, 22 degC, Lake Garda)
    env.water.rho = 997.8;                  % kg/m^3  (sea water: 1025)
    env.water.nu = 0.955e-6;                % m^2/s   (sea water 15C: 1.19e-6)

    %% Air (25 degC, 65 m elevation -> 100.55 kPa)
    env.air.rho = 1.175;                    % kg/m^3  (ISA SL: 1.225)
    env.air.nu = 1.56e-5;                   % m^2/s
    
    %% Wind
    env.wind.TWS_kn = 4:2:18;
    env.wind.TWA_deg = 30:10:180;
end