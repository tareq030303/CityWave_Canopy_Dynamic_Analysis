function [Gw_cut, vector_idx_sensors] = find_sensors(columns, G, idx_peak, mode, vector_threshold)
arguments (Input)
    columns
    G
    idx_peak
    mode
    vector_threshold
end

arguments (Output)
    Gw_cut
    vector_idx_sensors
end

vector_idx_sensors = [];

for j = 1:columns

    if abs(G(idx_peak(mode), j)) > vector_threshold(mode)

         vector_idx_sensors = [vector_idx_sensors j];
    end
end

Gw_cut = G(:, vector_idx_sensors);

end