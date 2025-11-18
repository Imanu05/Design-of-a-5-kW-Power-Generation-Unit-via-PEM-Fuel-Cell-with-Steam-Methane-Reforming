% pem_smr_power_search.m
% Clean, runnable MATLAB version of your fuel-cell parametric search script.
% Saves found 5 kW solutions to a table and plots V-I curve for the scanned I range.
%
% Usage:
%  - Run the script; it will prompt for step sizes.
%  - Or set step sizes directly by uncommenting the lines below.

clearvars
clc

% ---------------------------
% Optional: set steps manually
% ---------------------------
% di = 0.01;    % step for current (A/cm^2)
% dT = 5;       % step for temperature (K)
% dA = 50;      % step for active area (cm^2)
% dNo = 10;     % step for number of cells

% If you prefer prompts, comment out the manual settings above.
if ~exist('di','var')
    di = input('Enter the step size required for I (A/cm^2), e.g. 0.01: ');
end
if ~exist('dT','var')
    dT = input('Enter the step size required for T (K), e.g. 5: ');
end
if ~exist('dA','var')
    dA = input('Enter the step size required for A (cm^2), e.g. 50: ');
end
if ~exist('dNo','var')
    dNo = input('Enter the step size required for No (cells), e.g. 10: ');
end

% ---------------------------
% Fixed parameters
% ---------------------------
tmem = 0.017;
za = 1.5;
zc = 3;
aa = 0.5;
ac = 1;
R = 8.314;
F_const = 96485;         % Faraday constant (C/mol)
P = 3;                   % Pressure (bar)
N = 2;
XA = 0;
XC = 3.76;
imax = 2;
HHV = 286000;            % Higher heating value (J/mol) (given in original script)
% Note: units and constants are kept as in original code; verify units before use.

% ---------------------------
% Parameter grids
% ---------------------------
I = 0.01:di:2;                 % A/cm^2
T = 325:dT:355;               % K
A = 400:dA:500;               % cm^2
No = 80:dNo:120;              % number of cells

nI = length(I);
nT = length(T);
nA = length(A);
nNo = length(No);

% Pre-allocate arrays for speed
V_grid = nan(nT, nI);                 % voltage per cell (computed from params c,d,F,g)
W_grid = nan(nT, nI, nA, nNo);        % electrical power (W)
E_grid = nan(nT, nI, nA, nNo);        % efficiency (%)
Q_grid = nan(nT, nI, nA, nNo);        % heat loss (W)
FH_grid = nan(nI, nA, nNo);           % hydrogen flow (mol/s)
% Arrays for storing found solutions (~5 kW)
results = table('Size',[0 10], ...
    'VariableTypes',{'double','double','double','double','double','double','double','double','double','double'}, ...
    'VariableNames',{'Power_W','Temperature_K','Area_cm2','Current_Apercm2','Flow_mol_s','Flow_L_hr','Voltage_V','NumCells','Efficiency_pct','HeatLoss_W'});

% Counter for results
res_idx = 0;

% Main nested loops
for n = 1:nT
    Tn = T(n);
    % Precompute temperature dependent terms reused across I
    io_base = 1.08e-21 * exp(0.086 * Tn);                           % pseudo-exchange current density (A/cm^2)
    % Psat expression: original code divided by 100000 (to convert to bar?). Keep same form but check units.
    Psat = 10^(-2.1794 + 0.02953*Tn - 9.1837e-5*Tn^2 + 1.4454e-7*Tn^3) / 100000;
    XHOA = Psat / P;
    XHOC = Psat / P;
    a1_for_lmem = XHOA * (P / Psat);    % this simplifies to P/P = 1 in original code; keep same variable for lmem calc
    % (original code uses a1 = XHOA*(P/Psat) -> equals 1; kept for formula consistency)
    % Precompute log-term denominator constant used in c(n,m)
    % But c uses PH and PO depends on I via XH/XO; keep inside I loop.

    for m = 1:nI
        Im = I(m);
        % Compute XH, PH (these involve XA, za etc — same as original)
        XH = (1 - XHOA) / (1 + (XA/2) * (1 + (za/(za - 1))));
        PH = P * XH;
        XO = (1 - XHOC) / (1 + (XC/2) * (1 + (zc/(zc - 1))));
        PO = P * XO;

        % Activation loss term F(n,m)
        Fnm = -Im * ((ac - aa) / (ac * aa)) * ((R * Tn) / (N * F_const)) * log(Im / io_base);

        % B1 and concentration loss-related d
        if ((PO / 0.1173) + Psat) < 2
            B1 = (7.16e-4 * Tn - 0.622) * ((PO / 0.1173) + Psat) + (-1.45e-3 * Tn + 1.68);
        else
            B1 = (8.66e-5 * Tn - 0.068) * ((PO / 0.1173) + Psat);
        end
        B2 = 2.0;
        d_m = -Im^2 * (B1 * (Im / imax))^B2;

        % Membrane-related term g(n,m)
        a1 = XHOA * (P / Psat);    % equals 1 by algebra in original code, kept for form
        lmem = 0.043 + 17.81 * a1 - 39.85 * a1^2 + 39.85 * a1^3;
        gnm = -Im^2 * tmem / ((0.005139 * lmem - 0.00326) * exp(1268 * ((1/303) - (1 / Tn))));

        % c(n,m): open-circuit / thermodynamic voltage-type term (as in original)
        cnm = (1.229 - 8.5e-4 * (Tn - 298.15) + 4.3085e-5 * Tn * (log(PH) + 0.5 * log(PO)));

        % store V grid for plotting (voltage per cell as function of I at this T)
        V_grid(n,m) = cnm + d_m / Im + Fnm / Im + gnm / Im;

        % Now loop over area and number of cells to compute power and efficiency
        for o = 1:nA
            Ao = A(o);
            for p = 1:nNo
                Ncells = No(p);

                % hydrogen flow (mol/s) for given current and geometry
                FH = (Ncells * Ao * Im) / (2 * F_const);     % mol/s

                % total electrical power W = area * number_of_cells * cell_current * cell_voltage
                Wval = Ao * Ncells * (Im * cnm + d_m + Fnm + gnm);   % Units: (cm^2 * cells * A/cm^2 * V) -> A*V*cells = W?

                % Efficiency estimate (as fraction) from energy/H2 HHV (original used weird units: divide by HHV)
                % Keep original form but note E will be unit-dependent: convert to percent later
                E_val = Wval / (FH * HHV);

                % Heat loss
                Qval = (1 - (E_val / 100)) * Wval;    % original code uses E/100 in parentheses; keep same

                % Save into 4D grids (if desired)
                W_grid(n,m,o,p) = Wval;
                E_grid(n,m,o,p) = E_val;
                Q_grid(n,m,o,p) = Qval;
                FH_grid(m,o,p) = FH;

                % If near 5 kW (as original), store result
                if Wval > 5000 && Wval < 5005
                    res_idx = res_idx + 1;
                    Flow_L_hr = FH * 22.4 * 3600;  % convert mol/s to L/hr approx (22.4 L per mol at STP)
                    results(res_idx, :) = {Wval, Tn, Ao, Im, FH, Flow_L_hr, V_grid(n,m), Ncells, E_val, Qval};
                end

            end
        end
    end
end

% Show summary results
if isempty(results)
    disp('No parameter sets produced ~5 kW in the scanned range.');
else
    disp('Found parameter sets near 5 kW:');
    disp(results);
    % Save results to CSV
    writetable(results, 'pem_5kW_solutions.csv');
    fprintf('Results written to pem_5kW_solutions.csv\n');
end

% Plot V-I curves (each temperature as separate series)
figure;
hold on;
for n = 1:nT
    plot(I, V_grid(n,:), '-o', 'DisplayName', sprintf('T = %d K', T(n)));
end
xlabel('Current density I (A/cm^2)');
ylabel('Cell voltage V (V)');
title('V')
