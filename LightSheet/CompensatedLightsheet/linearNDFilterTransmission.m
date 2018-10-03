function linearNDFilterTransmission(ND_start,ND_end,ND_length,ND_samples,beam_pos,beam_length)

    if nargin < 1
        ND_start = 0; %underestimated
    end
    if nargin < 2
        ND_end = 4;
    end
    if nargin < 3
        ND_length = 50 * 1e-3; %[m] overestimated
    end
    if nargin < 4
        ND_samples = 1000;
    end
    if nargin < 5
        beam_pos = 10 * 1e-3; %[m]
    end
    if nargin < 6
        beam_length = 2.4*4 * 1e-3; %[m]
    end
    
    xRange = [1:ND_samples] / ND_samples * ND_length;
    xStep = xRange(2) - xRange(1);
    ND_range = [1:ND_samples] / ND_samples * (ND_end - ND_start);
    ND_transmission = 10.^(-1 .* ND_range);
    
    figure(1);
    subplot(2,1,1);
    plot(xRange * 1e3,ND_range);
    xlabel('x [mm]');
    ylabel('ND strength');
    title('ND filter properties');
    subplot(2,1,2);
    plot(xRange * 1e3,ND_transmission);
    xlabel('x [mm]');
    ylabel('ND transmisison');

    beam_pos_start_pixels = find(abs(xRange - beam_pos) <= (xStep / 2),1,'first');
    beam_pos_end_pixels = find(abs(xRange - (beam_pos + beam_length)) <= (xStep / 2),1,'first');

    figure(2);
    plot(xRange(beam_pos_start_pixels:beam_pos_end_pixels) * 1e3...
        ,ND_transmission(beam_pos_start_pixels:beam_pos_end_pixels) / min(ND_transmission(beam_pos_start_pixels:beam_pos_end_pixels)));
    xlabel('x [mm]');
    ylabel('ND transmisison');
    title('Beam profile');
end