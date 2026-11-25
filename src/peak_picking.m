function csi_exp = peak_picking(Gw_av, idx_peak, mode, csi_exp, w_exp, f)
arguments (Input)
    Gw_av
    idx_peak
    mode
    csi_exp
    w_exp
    f
end

arguments (Output)
    csi_exp
end

thres = abs(Gw_av(idx_peak(mode))) / sqrt(2);

% avanti
esc = false;
idx = idx_peak(mode) + 1;

while esc == false

    if abs(Gw_av(idx)) < thres

        esc = true;
        idx_w2 = idx - 1;
    else

        idx = idx + 1;
    end
end

% indietro
esc = false;
idx = idx_peak(mode) - 1;

while esc == false

    if abs(Gw_av(idx)) < thres

        esc = true;
        idx_w1 = idx + 1;
    else

        idx = idx - 1;
    end
end

csi_exp(mode) = ((f(idx_w2)*2*pi)^2 - (f(idx_w1)*2*pi)^2) / (4*(w_exp(mode)^2));

end