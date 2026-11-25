%% CLEANING THE FILE
close all
clear all
clc


%% LOAD THE DATA

connectivity = load("connectivity.mat");
FRF_H1 = load("FRF_H1");
punti_laser = load("punti_laser.mat");    
modal_output = load("modal_output.mat");


% Exploit with a different name the contents of the data
f = FRF_H1.f;
G = FRF_H1.FRF_z;

% Define the range of interest (OBC method)
idx_min = 163;
idx_max = 656;
f = f(idx_min:idx_max);
G = G(idx_min:idx_max, :);
[rows, columns] = size(G);


%% PLOT THE MAGNITUDE AND THE PHASE

% Magnitude *without* the mean
figure;
subplot(2, 1, 1);
semilogy(f, abs(G), "color", [0.5 0.5 0.5]);
xlabel('Frequency (Hz)', 'Interpreter', 'Latex');
ylabel('Module (m/N)', 'Interpreter', 'Latex');
xlim tight
ylim tight
grid on;

% Phase *without* the mean
subplot(2, 1, 2);
plot(f, wrapTo2Pi(angle(G)), "color", [0.5 0.5 0.5]);
xlabel('Frequency (Hz)', 'Interpreter', 'Latex');
ylabel('Phase (rad)', 'Interpreter', 'Latex');
xlim tight;
ylim padded;
yticks(0 : pi/2 : 2*pi);
grid on;

% Calculation of the mean
Gabs_mean = mean(abs(G), 2); % mean across the columns and we get a column vector

% Magnitude *with* the mean
figure;
semilogy(f, abs(G), "color", [0.5 0.5 0.5]);
hold on
semilogy(f, Gabs_mean, 'b');
xlabel('Frequency (Hz)', 'Interpreter', 'Latex');
ylabel('Module (m/N)', 'Interpreter', 'Latex');
xlim tight
ylim tight
grid on;


%% DETERMINE THE EIGENFREQUENCIES

% To determine the eigenfrequencies, I use findpeaks
G_mean_int1 = Gabs_mean(1:109);
G_mean_int2 = Gabs_mean(110:133);
G_mean_int3= Gabs_mean(134:166);
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

% Calculate the eigenfrequencies
f_exp = f(idx_peak);
w_exp = 2*pi*f(idx_peak);

% Print on the command window the values
disp('WT natural frequencies BEFORE scaling:');
for i = 1:n_modes
    fprintf('Mode %d: %f Hz\n', i, f_exp(i));
end

% SCALING THE FREQUENCY (necessary for our analysis): Pass from the wind tunnel's frequency to building's frequency
lambda = 1/70;
norm = sqrt(lambda);
f_WT_scaled = (f_exp*norm)';

% Print on the command window the values
disp('WT natural frequencies AFTER scaling:');
for i = 1:n_modes
    fprintf('Mode %d: %f Hz\n', i, f_WT_scaled(i));
end

% Frequencies of the building in scale 1:1
f_fem = modal_output.modpar.freq(1:4);

% Print on the command window the values
disp('FEM natural frequencies:');
for i = 1:n_modes
    fprintf('Mode %d: %f Hz\n', i, f_fem(i));
end

%% FIND THE PEAKS AND HIGHLIGHT THE SENSORS OF INTEREST FOR EACH MODE
% Apply the Peak Picking Method
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

% mode shape 1
mode = 1;
[Gw1_cut, vector_idx_sensors_w1] = find_sensors(columns, G, idx_peak, mode, vector_threshold);

% mode shape 2
mode = 2;
[Gw2_cut, vector_idx_sensors_w2] = find_sensors(columns, G, idx_peak, mode, vector_threshold);

% mode shape 3
mode = 3;
[Gw3_cut, vector_idx_sensors_w3] = find_sensors(columns, G, idx_peak, mode, vector_threshold);

% mode shape 4
mode = 4;
[Gw4_cut, vector_idx_sensors_w4] = find_sensors(columns, G, idx_peak, mode, vector_threshold);


%% PLOTS OF THE SENSORS ASSOCIATED TO EACH MODE

% Verify that for each natural frequency we've chosen the FRF with the highest magnitude
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

% Average in the complex domain of the FRF for each mode
Gw1_av = mean(Gw1_cut, 2); % mean across the columns (sensors) and the result is a column vector
Gw2_av = mean(Gw2_cut, 2);
Gw3_av = mean(Gw3_cut, 2);
Gw4_av = mean(Gw4_cut, 2);

% Verify that for each natural frequency we've done correctly the mean
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


%% APPLY THE PEAK PICKING METHOD FOR EACH MODE
% Calculation of csi_exp
csi_exp = zeros(1, length(w_exp));

% mode shape 1
mode = 1;
csi_exp = peak_picking(Gw1_av, idx_peak, mode, csi_exp, w_exp, f);

% mode shape 2
mode = 2;
csi_exp = peak_picking(Gw2_av, idx_peak, mode, csi_exp, w_exp, f);

% mode shape 3
mode = 3;
csi_exp = peak_picking(Gw3_av, idx_peak, mode, csi_exp, w_exp, f);

% mode shape 4
mode = 4;
csi_exp = peak_picking(Gw4_av, idx_peak, mode, csi_exp, w_exp, f);

% Plot on the command window the results of all the peaks
disp('Adimensional Damping Ratios (Half-Power Method):');
for i = 1:n_modes
    fprintf('Mode %d: %.2f Hz -> csi = %.2f%%\n', i, f_WT_scaled(i), csi_exp(i)*100);
end

%% FIND THE MODAL SHAPES
% To find the mode shapes, we will apply the following formula: phi_ki(x)=Imag(G_ij(omega)
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

% Normalize the quantity
phi_sens_norm = phi_sens./phi_sens_max;


%% REPRESENTATION OF THE POSITION OF SENSORS
% Load again the data to acquire the correct name
load("FRF_H1.mat")
load("punti_laser.mat")
load("modal_output.mat")
load("connectivity.mat")

% nodes in the Canopy from FE
[~,Loc1] = ismember(surface_canopy.Joint1,nodes.ID);
[~,Loc2] = ismember(surface_canopy.Joint2,nodes.ID);
[~,Loc3] = ismember(surface_canopy.Joint3,nodes.ID);
nodi123 = [Loc1,Loc2,Loc3];

% scale factor to display mode shapes
scalaFEM = 10;
scalaLASER = 10;

% plot FEM model + experimental measuring point
figure;
cc=ones(length(nodes.X),1);
patch('Faces',nodi123,'Vertices',[nodes.X,nodes.Y,nodes.Z],'CData',cc,'FaceColor','interp','EdgeColor','none')
axis tight
hold on;
plot3(x,y,z+1,'ro','MarkerSize',4,'MarkerFaceColor','r')


%% FEM vs EMA
phi_fem = [];
for ii = 1:4
    
    %FEM
    mode = ii;
    mode_sel = modeshapes(modeshapes.No == mode,:);
    [~,Locb] = ismember(mode_sel.ID,nodes.ID);
    modedef  = mode_sel{Locb,[{'uX'},{'uY'},{'uZ'}]};

    fff = max(abs(modedef(:,3)));

    % EMA - Doing the interpolation
    F = scatteredInterpolant(x, y, z, phi_sens_norm(ii,:)', 'natural', 'nearest');
    phi_exp_FEM(ii,:) = F(nodes.X, nodes.Y, nodes.Z);

    % PLOT
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


%% COMPARE THE MODE SHAPES AND THE EIGENFREQUENCIES 

% It has to be performed between:
% FOR FEM => nodes
% FOR EXPERIMENTAL => phi_exp_FEM; punti_laser
% We must consider that our node is our sensor of the FEM. 
% We must highlight only the sensors which are over the sensor

% Define the position of the lasers and the nodes
laser_position = [punti_laser.x(:), punti_laser.y(:), punti_laser.z(:)];
nodes_position = [connectivity.nodes.X(:), connectivity.nodes.Y(:), connectivity.nodes.Z(:)];

% Define the total number of nodes
nodes_number = length(nodes_position(:,1));

% Create the distance matrix
distance = zeros(columns, nodes_number);
for i = 1:columns
    for j = 1:nodes_number
       distance(i,j) = sqrt((laser_position(i,1) - nodes_position(j,1))^2 + (laser_position(i,2) - nodes_position(j,2))^2 + (laser_position(i,3) - nodes_position(j,3))^2);
    end
end

% Create the vector of the neareast distance between each sensor
opt_nodes = [];
for i = 1:columns
    [~, R] = find(distance==min(distance(i,:)));
    opt_nodes = [opt_nodes R];
end

% Define the mode shapes associated to the nodes of our interest
phi_fem_reduced = phi_fem(:, opt_nodes);
nodes_values = nodes_position(opt_nodes, :);

% Do the MAC (now that we have both the phi)
% phi_sens and phi_fem_reduced must be considered
num_tot = [];
den1_tot = [];
den2_tot = [];

for i = 1:n_modes  % non cambia niente dal coniugato del complesso, è normale, no?
    num = sum(conj(phi_fem_reduced(i,:)).*phi_sens_norm(i,:));
    den1 = sum(conj(phi_fem_reduced(i,:)).*phi_fem_reduced(i,:));
    den2 = sum(conj(1i*phi_sens_norm(i,:)).*phi_sens_norm(i,:));
    den1_tot = [den1_tot den1];
    den2_tot = [den2_tot den2];
    num_tot = [num_tot num];
end

% Calculate the MAC
for i = 1:n_modes
    mac(i) = (abs(num_tot(i))^2) / ( (abs(den1_tot(i))) * (abs(den2_tot(i))) );
end

% Calculate the average MAC value for all modes (forse inutile)
average_mac = mean(mac);
disp(['Average MAC: ', num2str(average_mac)]);


% Compare the eigenfrequencies and its percentage error
w_fem = 2 * pi * modal_output.modpar.freq;
w_exp_scaled = 2 * pi * f_WT_scaled;

% Calculate percentage error
for i = 1:n_modes
    error_perc(i) = (abs(w_fem(i) - w_exp_scaled(i))*100)/(w_fem(i));
end

% Calculate absolute error
for i = 1:n_modes
    error_perc(i) = abs(w_fem(i) - w_exp_scaled(i));
end

table_of_w = [w_exp_scaled(:), w_fem(1:4)];

%% CoMAC
% Do the COMAC (now that we have both the phi)
n_nodes = size(phi_fem_reduced, 2);
n_modi = size(phi_fem_reduced, 1);

num_tot_comac = [];
den1_tot_comac = [];
den2_tot_comac = [];

for j = 1:n_nodes  % non cambia niente dal coniugato del complesso, è normale, no?
    num = sum(abs(phi_fem_reduced(:, j).*phi_sens_norm(:,j)))^2;
    den1 = sum(phi_fem_reduced(:, j).^2);
    den2 = sum(phi_sens_norm(:, j).^2);
    den1_tot_comac = [den1_tot_comac den1];
    den2_tot_comac = [den2_tot_comac den2];
    num_tot_comac = [num_tot_comac num];
end

% Calculate the CoMAC
comac = zeros(1, n_nodes);
for j = 1:n_nodes
    comac(j) = num_tot_comac(j) / (den1_tot_comac(j) * den2_tot_comac(j));
end

%% Plot CoMAC
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

% Map
c1 = (linspace(1, 0, 20))';
c2 = (linspace(0, 1, 20))';
c3 = zeros(20, 1);
map = [c1 c2 c3];
colorbar;
colormap(map);
clim([0 1]);

%% CODICE DELLA PROF

% %% PURE FEM RESULTS (FROM POMARANZI)
% phi_fem = [];
% 
% for ii = 1:4
% 
%     %FEM
%     mode = ii;
%     mode_sel = modeshapes(modeshapes.No == mode,:);
%     [~,Locb] = ismember(mode_sel.ID,nodes.ID);
%     modedef  = mode_sel{Locb,[{'uX'},{'uY'},{'uZ'}]};
% 
%     fff = max(abs(modedef(:,3)));
% 
%     figure;
%     axis equal
%     grid on
%     hold on
%     view(2)
% 
%     cc=(modedef(:,3))/fff;
%     patch('Faces',nodi123,'Vertices',[nodes.X,nodes.Y,nodes.Z]+modedef/fff*scalaFEM,'CData',cc,'FaceColor','interp','EdgeColor','none')
%     title(sprintf('FEM -- Mode %i: f=%5.3f Hz', mode,modpar.freq(mode)) )
%     nmap = 10;
%     map     = jet(nmap);
%     colormap(gca,map)
%     clim([-1 1])
%     cb = colorbar('eastoutside');
%     axis tight
% 
%     phi_fem = [phi_fem cc];
% end
% phi_fem = phi_fem';
% 
% 
% 
% %% DO THE INTERPOLATION "NODES-SENSORS"
% 
% % Do the scattering and then project/determine the values of scattering at
% % the points of the nodes
% for i = 1:n_modes
%     F = scatteredInterpolant(x, y, z, phi_sens_norm(i,:)', 'natural', 'nearest');
%     phi_exp_FEM(i,:) = F(nodes.X, nodes.Y, nodes.Z);
% end
% 
% % Plot the FEM exploited through the experimental data
% for i = 1:n_modes
%     figure;
%     axis equal
%     grid on
%     hold on
%     view(2)
%     patch('Faces', nodi123, 'Vertices', [nodes.X, nodes.Y, nodes.Z], 'FaceVertexCData', phi_exp_FEM(i,:)', 'FaceColor', 'interp', 'EdgeColor', 'none');
%     title(sprintf('FEM -- Mode %i: f=%5.3f Hz csi=%.2f%%', i, f_WT_scaled(i), csi_exp(i)*100))
%     nmap = 10;
%     map     = jet(nmap);
%     colormap(gca,map)
%     clim([-1 1])
%     cb = colorbar('eastoutside');
%     axis tight
% end

%% COMPARISON BETWEEN FEM AND EXPERIMENTAL RESULTS  (ERA SBAGLIATO!!!)

% cc_abs = modedef(:,3) / max(abs(modedef(:,3)));
% Apply the new frequency for the
% for i = 1:n_modes
%     error(i,:) = (phi_exp_FEM(i,:) - phi_fem(i,:));
% end
% 
% % Plot the FEM exploited through the experimental data
% for i = 1:n_modes
%     figure;
%     patch('Faces', nodi123, 'Vertices', [nodes.X, nodes.Y, nodes.Z], 'FaceVertexCData', error(i,:)', 'FaceColor', 'interp', 'EdgeColor', 'none');
%     title(sprintf('FEM -- Mode %i: f=%5.3f Hz', i, f(idx_peak(i))))
%     nmap = 10;
%     map = jet(nmap);
%     colormap(gca,map)
%     clim([0 1])
%     clim([-1 1])
%     cb = colorbar('eastoutside');
%     axis tight
% end