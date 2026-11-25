%% WORKSPACE INITIALIZATION
close all
clear all
clc

%% LOAD DATA
connectivity = load("connectivity.mat");
FRF_H1 = load("FRF_H1");
punti_laser = load("punti_laser.mat");    
modal_output = load("modal_output.mat");

% Extract frequency and FRF data
f = FRF_H1.f;
G = FRF_H1.FRF_z;

% Define the frequency range of interest (OBC method)
idx_min = 163;
idx_max = 656;
f = f(idx_min:idx_max);
G = G(idx_min:idx_max, :);
[rows, columns] = size(G);

%% PLOT MAGNITUDE AND PHASE

% Magnitude (individual sensors)
figure;
subplot(2, 1, 1);
semilogy(f, abs(G), "color", [0.5 0.5 0.5]);
xlabel('Frequency (Hz)', 'Interpreter', 'Latex');
ylabel('Module (m/N)', 'Interpreter', 'Latex');
xlim tight
ylim tight
grid on;

% Phase (individual sensors)
subplot(2, 1, 2);
plot(f, wrapTo2Pi(angle(G)), "color", [0.5 0.5 0.5]);
xlabel('Frequency (Hz)', 'Interpreter', 'Latex');
ylabel('Phase (rad)', 'Interpreter', 'Latex');
xlim tight;
ylim padded;
yticks(0 : pi/2 : 2*pi);
grid on;

% Compute spatial mean across sensors
Gabs_mean = mean(abs(G), 2);

% Magnitude with spatial mean
figure;
semilogy(f, abs(G), "color", [0.5 0.5 0.5]);
hold on
semilogy(f, Gabs_mean, 'b');
xlabel('Frequency (Hz)', 'Interpreter', 'Latex');
ylabel('Module (m/N)', 'Interpreter', 'Latex');
xlim tight
ylim tight
grid on;

%% EIGENFREQUENCY IDENTIFICATION

% Partition frequency spectrum to locate modal peaks
G_mean_int1 = Gabs_mean(1:109);
G_mean_int2 = Gabs_mean(110:133);
G_mean_int3 = Gabs_mean(134:166);
G_mean_int4 = Gabs_mean(167:248);

max1 = max(G_mean_int1);
max2 = max(G_mean_int2);
max3 = max(G_mean_int3);
max4 = max(G_mean_int4);

maxx_val = [max1 max2 max3 max4];
idx_peak = zeros(1, 4);
for j = 1:length(maxx_val)
    for i = 1:length(Gabs_mean)
        if Gabs_mean(i) == maxx_val(j)
            idx_peak(j) = i;
        end
    end
end

n_modes = length(idx_peak);

figure;
semilogy(f, abs(G), "color", [0.5 0.5 0.5]);
hold on
semilogy(f, Gabs_mean, 'b');
hold on
plot(f(idx_peak), maxx_val, 'or' , 'MarkerFaceColor', 'r');
xlabel('Frequency (Hz)', 'Interpreter', 'Latex');
ylabel('Module (m/N)', 'Interpreter', 'Latex');
xlim tight
ylim tight
grid on;

% Extract experimental eigenfrequencies
f_exp = f(idx_peak);
w_exp = 2*pi*f(idx_peak);

disp('WT natural frequencies BEFORE scaling:');
for i = 1:n_modes
    fprintf('Mode %d: %f Hz\n', i, f_exp(i));
end

% Frequency scaling from wind tunnel model (1:70) to full-scale structure
lambda = 1/70;
norm = sqrt(lambda);
f_WT_scaled = (f_exp*norm)';

disp('WT natural frequencies AFTER scaling:');
for i = 1:n_modes
    fprintf('Mode %d: %f Hz\n', i, f_WT_scaled(i));
end

% Full-scale FEM natural frequencies
f_fem = modal_output.modpar.freq(1:4);

disp('FEM natural frequencies:');
for i = 1:n_modes
    fprintf('Mode %d: %f Hz\n', i, f_fem(i));
end

%% PEAK PICKING AND ANTI-NODAL SENSOR SELECTION
vector_j_max = zeros(4, 1);
vector_frf_max = zeros(4, 1);

for i = 1:length(w_exp)
    vet_sup = abs(G(idx_peak(i), :));
    j_max = 1;
    frf_max = vet_sup(j_max);

    for j = 2:length(vet_sup)
        if vet_sup(j) > frf_max
            j_max = j;
            frf_max = vet_sup(j_max);
        end
    end

    vector_j_max(i) = j_max;
    vector_frf_max(i) = frf_max;
end

vector_threshold = [vector_frf_max(1) * 0.8, vector_frf_max(2) * 0.7, vector_frf_max(3) * 0.8, vector_frf_max(4) * 0.7];

% Mode 1 sensors
mode = 1;
[Gw1_cut, vector_idx_sensors_w1] = find_sensors(columns, G, idx_peak, mode, vector_threshold);

% Mode 2 sensors
mode = 2;
[Gw2_cut, vector_idx_sensors_w2] = find_sensors(columns, G, idx_peak, mode, vector_threshold);

% Mode 3 sensors
mode = 3;
[Gw3_cut, vector_idx_sensors_w3] = find_sensors(columns, G, idx_peak, mode, vector_threshold);

% Mode 4 sensors
mode = 4;
[Gw4_cut, vector_idx_sensors_w4] = find_sensors(columns, G, idx_peak, mode, vector_threshold);

%% PLOT SELECTED SENSORS PER MODE

figure; 
subplot(2, 2, 1);
semilogy(f, abs(G), "color", [0.5 0.5 0.5]);
hold on
semilogy(f, abs(Gw1_cut), 'b');
hold on
xline(f(idx_peak(1)), '--');
title('MODE SHAPE 1: FRF in anti-nodal positions');
xlabel('Frequency (Hz)', 'Interpreter', 'Latex');
ylabel('Module (m/N)', 'Interpreter', 'Latex');
xlim tight
ylim tight
grid on;

subplot(2, 2, 2);
semilogy(f, abs(G), "color", [0.5 0.5 0.5]);
hold on
semilogy(f, abs(Gw2_cut), 'b');
hold on
xline(f(idx_peak(2)), '--');
title('MODE SHAPE 2: FRF in anti-nodal positions');
xlabel('Frequency (Hz)', 'Interpreter', 'Latex');
ylabel('Module (m/N)', 'Interpreter', 'Latex');
xlim tight
ylim tight
grid on;

subplot(2, 2, 3);
semilogy(f, abs(G), "color", [0.5 0.5 0.5]);
hold on
semilogy(f, abs(Gw3_cut), 'b');
hold on
xline(f(idx_peak(3)), '--');
title('MODE SHAPE 3: FRF in anti-nodal positions');
xlabel('Frequency (Hz)', 'Interpreter', 'Latex');
ylabel('Module (m/N)', 'Interpreter', 'Latex');
xlim tight
ylim tight
grid on;

subplot(2, 2, 4);
semilogy(f, abs(G), "color", [0.5 0.5 0.5]);
hold on
semilogy(f, abs(Gw4_cut), 'b');
hold on
xline(f(idx_peak(4)), '--');
title('MODE SHAPE 4: FRF in anti-nodal positions');
xlabel('Frequency (Hz)', 'Interpreter', 'Latex');
ylabel('Module (m/N)', 'Interpreter', 'Latex');
xlim tight
ylim tight
grid on;

% Complex-domain average of selected FRFs per mode
Gw1_av = mean(Gw1_cut, 2);
Gw2_av = mean(Gw2_cut, 2);
Gw3_av = mean(Gw3_cut, 2);
Gw4_av = mean(Gw4_cut, 2);

% Verify complex mean FRFs
figure; 
subplot(2, 2, 1);
semilogy(f, abs(Gw1_cut), 'b');
hold on
semilogy(f, abs(Gw1_av), 'r');
hold on
xline(f(idx_peak(1)), '--');
title('MODE SHAPE 1: complex mean');
xlabel('Frequency (Hz)', 'Interpreter', 'Latex');
ylabel('Module (m/N)', 'Interpreter', 'Latex');
xlim tight
ylim tight
grid on;

subplot(2, 2, 2);
semilogy(f, abs(Gw2_cut), 'b');
hold on
semilogy(f, abs(Gw2_av), 'r');
hold on
xline(f(idx_peak(2)), '--');
title('MODE SHAPE 2: complex mean');
xlabel('Frequency (Hz)', 'Interpreter', 'Latex');
ylabel('Module (m/N)', 'Interpreter', 'Latex');
xlim tight
ylim tight
grid on;

subplot(2, 2, 3);
semilogy(f, abs(Gw3_cut), 'b');
hold on
semilogy(f, abs(Gw3_av), 'r');
hold on
xline(f(idx_peak(3)), '--');
title('MODE SHAPE 3: complex mean');
xlabel('Frequency (Hz)', 'Interpreter', 'Latex');
ylabel('Module (m/N)', 'Interpreter', 'Latex');
xlim tight
ylim tight
grid on;

subplot(2, 2, 4);
semilogy(f, abs(Gw4_cut), 'b');
hold on
semilogy(f, abs(Gw4_av), 'r');
hold on
xline(f(idx_peak(4)), '--');
title('MODE SHAPE 4: complex mean');
xlabel('Frequency (Hz)', 'Interpreter', 'Latex');
ylabel('Module (m/N)', 'Interpreter', 'Latex');
xlim tight
ylim tight
grid on;

%% DAMPING RATIO ESTIMATION (HALF-POWER BANDWIDTH METHOD)
csi_exp = zeros(1, length(w_exp));

% Mode 1
mode = 1;
csi_exp = peak_picking(Gw1_av, idx_peak, mode, csi_exp, w_exp, f);

% Mode 2
mode = 2;
csi_exp = peak_picking(Gw2_av, idx_peak, mode, csi_exp, w_exp, f);

% Mode 3
mode = 3;
csi_exp = peak_picking(Gw3_av, idx_peak, mode, csi_exp, w_exp, f);

% Mode 4
mode = 4;
csi_exp = peak_picking(Gw4_av, idx_peak, mode, csi_exp, w_exp, f);

disp('Adimensional Damping Ratios (Half-Power Method):');
for i = 1:n_modes
    fprintf('Mode %d: %.2f Hz -> csi = %.2f%%\n', i, f_WT_scaled(i), csi_exp(i)*100);
end

%% EXPERIMENTAL MODE SHAPE EXTRACTION
% Extract mode shapes from the imaginary part of the FRF at resonance: phi_ki(x) = Imag(G_ij(omega))
phi_sens = zeros(n_modes, columns);
for i = 1:n_modes
    if i==3
        phi_sens(i,:) = imag(G(idx_peak(i), :));
    else
        phi_sens(i,:) = -imag(G(idx_peak(i), :));
    end
end

phi_sens_max = zeros(n_modes, 1);
for i = 1:n_modes
    massimo = abs(phi_sens(i, 1));
    for j = 2:columns
        if abs(phi_sens(i, j)) > massimo
            massimo = abs(phi_sens(i, j));
        end
    end
    phi_sens_max(i) = massimo;
end

% Normalize mode shapes to unit maximum displacement
phi_sens_norm = phi_sens./phi_sens_max;

%% SENSOR LAYOUT VISUALIZATION ON FEM MESH
load("FRF_H1.mat")
load("punti_laser.mat")
load("modal_output.mat")
load("connectivity.mat")

% Extract canopy surface connectivity from FE model
[~,Loc1] = ismember(surface_canopy.Joint1,nodes.ID);
[~,Loc2] = ismember(surface_canopy.Joint2,nodes.ID);
[~,Loc3] = ismember(surface_canopy.Joint3,nodes.ID);
nodi123 = [Loc1,Loc2,Loc3];

% Scale factors for mode shape visualization
scalaFEM = 10;
scalaLASER = 10;

% Plot FEM geometry overlaid with laser measurement points
figure;
cc=ones(length(nodes.X),1);
patch('Faces',nodi123,'Vertices',[nodes.X,nodes.Y,nodes.Z],'CData',cc,'FaceColor','interp','EdgeColor','none')
axis tight
hold on;
plot3(x,y,z+1,'ro','MarkerSize',4,'MarkerFaceColor','r')

%% FEM VS EMA MODE SHAPE COMPARISON
phi_fem = [];
for ii = 1:4
    % FEM mode shape extraction
    mode = ii;
    mode_sel = modeshapes(modeshapes.No == mode,:);
    [~,Locb] = ismember(mode_sel.ID,nodes.ID);
    modedef  = mode_sel{Locb,[{'uX'},{'uY'},{'uZ'}]};

    fff = max(abs(modedef(:,3)));

    % EMA spatial interpolation onto FEM nodes
    F = scatteredInterpolant(x, y, z, phi_sens_norm(ii,:)', 'natural', 'nearest');
    phi_exp_FEM(ii,:) = F(nodes.X, nodes.Y, nodes.Z);

    % Side-by-side comparison plot
    figure;
    subplot(1, 2, 1);
    axis equal
    grid on
    hold on
    view(2)

    cc=(modedef(:,3))/fff;
    patch('Faces',nodi123,'Vertices',[nodes.X,nodes.Y,nodes.Z]+modedef/fff*scalaFEM,'CData',cc,'FaceColor','interp','EdgeColor','none')
    title(sprintf('FEM -- Mode %d: f=%5.3f Hz', mode,modpar.freq(mode)) )
    nmap = 10;
    map     = jet(nmap);
    colormap(gca,map)
    clim([-1 1])
    cb = colorbar('eastoutside');
    axis tight

    phi_fem = [phi_fem cc];

    subplot(1, 2, 2);
    axis equal
    grid on
    hold on
    view(2)
    patch('Faces', nodi123, 'Vertices', [nodes.X, nodes.Y, nodes.Z], 'FaceVertexCData', phi_exp_FEM(ii,:)', 'FaceColor', 'interp', 'EdgeColor', 'none');
    title(sprintf('EMA -- Mode %d: f=%5.3f Hz csi=%.2f%%', ii, f_WT_scaled(ii), csi_exp(ii)*100))
    nmap = 10;
    map     = jet(nmap);
    colormap(gca,map)
    clim([-1 1])
    cb = colorbar('eastoutside');
    axis tight
end
phi_fem = phi_fem';

%% MODAL ASSURANCE CRITERION (MAC) AND FREQUENCY ERROR ANALYSIS

% Coordinates of experimental laser points and FEM nodes
laser_position = [punti_laser.x(:), punti_laser.y(:), punti_laser.z(:)];
nodes_position = [connectivity.nodes.X(:), connectivity.nodes.Y(:), connectivity.nodes.Z(:)];

nodes_number = length(nodes_position(:,1));

% Euclidean distance matrix between measurement points and FEM nodes
distance = zeros(columns, nodes_number);
for i = 1:columns
    for j = 1:nodes_number
       distance(i,j) = sqrt((laser_position(i,1) - nodes_position(j,1))^2 + (laser_position(i,2) - nodes_position(j,2))^2 + (laser_position(i,3) - nodes_position(j,3))^2);
    end
end

% Map each laser measurement point to its nearest FEM node
opt_nodes = [];
for i = 1:columns
    [~, R] = find(distance==min(distance(i,:)));
    opt_nodes = [opt_nodes R];
end

% Extract FEM mode shapes at matched sensor locations
phi_fem_reduced = phi_fem(:, opt_nodes);
nodes_values = nodes_position(opt_nodes, :);

% Compute MAC components
num_tot = [];
den1_tot = [];
den2_tot = [];

for i = 1:n_modes
    num = sum(conj(phi_fem_reduced(i,:)).*phi_sens_norm(i,:));
    den1 = sum(conj(phi_fem_reduced(i,:)).*phi_fem_reduced(i,:));
    den2 = sum(conj(1i*phi_sens_norm(i,:)).*phi_sens_norm(i,:));
    den1_tot = [den1_tot den1];
    den2_tot = [den2_tot den2];
    num_tot = [num_tot num];
end

% Compute diagonal MAC values
for i = 1:n_modes
    mac(i) = (abs(num_tot(i))^2) / ( (abs(den1_tot(i))) * (abs(den2_tot(i))) );
end

average_mac = mean(mac);
disp(['Average MAC: ', num2str(average_mac)]);

% Eigenfrequency comparison (rad/s)
w_fem = 2 * pi * modal_output.modpar.freq;
w_exp_scaled = 2 * pi * f_WT_scaled;

% Percentage error
for i = 1:n_modes
    error_perc(i) = (abs(w_fem(i) - w_exp_scaled(i))*100)/(w_fem(i));
end

% Absolute error
for i = 1:n_modes
    error_perc(i) = abs(w_fem(i) - w_exp_scaled(i));
end

table_of_w = [w_exp_scaled(:), w_fem(1:4)];

%% COORDINATE MODAL ASSURANCE CRITERION (CoMAC)
n_nodes = size(phi_fem_reduced, 2);
n_modi = size(phi_fem_reduced, 1);

num_tot_comac = [];
den1_tot_comac = [];
den2_tot_comac = [];

for j = 1:n_nodes
    num = sum(abs(phi_fem_reduced(:, j).*phi_sens_norm(:,j)))^2;
    den1 = sum(phi_fem_reduced(:, j).^2);
    den2 = sum(phi_sens_norm(:, j).^2);
    den1_tot_comac = [den1_tot_comac den1];
    den2_tot_comac = [den2_tot_comac den2];
    num_tot_comac = [num_tot_comac num];
end

comac = zeros(1, n_nodes);
for j = 1:n_nodes
    comac(j) = num_tot_comac(j) / (den1_tot_comac(j) * den2_tot_comac(j));
end

%% PLOT CoMAC
comac = comac';
[X, Y] = meshgrid(x, y);
comac_grid = diag(comac);
mask = not(comac_grid);
comac_grid(mask) = NaN;

figure;
cc=ones(length(nodes.X), 1);
patch('Faces', nodi123, 'Vertices', [nodes.X, nodes.Y, nodes.Z], 'CData', cc, 'FaceColor', 'none', 'EdgeColor', 'none');
axis tight;
hold on;
scatter(X(:), Y(:), 50, comac_grid(:), 'filled');

% Custom red-to-green colormap
c1 = (linspace(1, 0, 20))';
c2 = (linspace(0, 1, 20))';
c3 = zeros(20, 1);
map = [c1 c2 c3];
colorbar;
colormap(map);
clim([0 1]);