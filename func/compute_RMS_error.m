function [Mic_pos_err, Source_pos_err] = compute_RMS_error(g)
Mic_pos_err = 0;
Source_pos_err = 0;

    for i = 2:g.M
        mic_idx = (i - 1) * 5;

        Mic_pos_err = Mic_pos_err + ((g.x_gt(mic_idx + 1) - g.x(mic_idx + 1))^2 + ...
        (g.x_gt(mic_idx + 2) - g.x(mic_idx + 2))^2 + (g.x_gt(mic_idx + 3) - g.x(mic_idx + 3))^2);
          
    end

    % Loop over all edges
    for eid = 1:length(g.edges)
        edge = g.edges(eid);
        % pose-pose constraint
        if (strcmp(edge.type, 'P') ~= 0)
            x = g.x(edge.fromIdx:edge.fromIdx+2);
            x_gt = g.x_gt(edge.fromIdx:edge.fromIdx+2);

            Source_pos_err = Source_pos_err + ((x_gt(1) - x(1))^2 + (x_gt(2) - x(2))^2 + (x_gt(3) - x(3))^2);
        end           
    end
    Source_pos_err = (Source_pos_err / ((length(g.x)-40) / 3))^0.5;
    Mic_pos_err = (Mic_pos_err / (g.M - 1))^0.5;
end