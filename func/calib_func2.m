function calib_func2(input)
% This file performs least square SLAM
% 
%% parameters

% gonna estimate clock drift?
est_drift_on = 0;
% gonna estimate starting time delay?
est_delay_on = 0;
% display starting time delay estimation error?
display_delay_error_on = 0;
% display norm(dx) for each iteration?
display_norm_dx_on = 0;

% the maximum number of iterations
numIterations = 50;

% maximum allowed dx
EPSILON = input.eps;%1.5/1e-2;

% Error
err = 0;

% load the graph into the variable "g"
load(input.graph_file);

% if est_drift_on is not enabled, assign the ground truth values
if est_drift_on<1
    for n = 2:g.M
        g.x(5*(n-1)+5) = g.x_gt(5*(n-1)+5);
    end
end

% if est_delay_on is not enabled, assign the ground truth values
if est_delay_on<1
    for n = 2:g.M
        g.x(5*(n-1)+4) = g.x_gt(5*(n-1)+4);
    end
end

sample_num = 100;%
threshold = 1; % 菱形1.2， 螺旋线3.5,random 1.5
decrease_threshold = -1;
total_sample_num = (length(g.edges)+1)/2;

for i = 1 : (length(g.edges))
    g.edges(i).toIdx_original = g.edges(i).toIdx;
end
%用g_online去迭代，计算eigenvalue与信息增益
g.M_x = 2;
g.M_y = 2;
g.M_z = 2;
g_online = struct();
g_online.M_x = g.M_x;
g_online.M_y = g.M_y;
g_online.M_z = g.M_z;
g_online.M = g.M;
% g_online.mic_dis = g.mic_dis;
g_online.x_gt = g.x_gt(1:(5 * g.M + 3 * sample_num),:); 
g_online.x = g.x(1:(5 * g.M + 3 * sample_num),:);
g_online.edges = g.edges(1:(2 * sample_num - 1)); 
g_online.idLookup = g.idLookup(1:(g.M + sample_num)); 

batch = 2;%一批处理两个点
%% start slSLAM



%     dx_norms = zeros(numIterations, 1);
    % carry out the iterations
    for i = 1:numIterations
    %   disp(['Performing iteration ', num2str(i)]);
      [Mic_pos_err, Source_pos_err] = compute_RMS_error(g_online);
      disp([Mic_pos_err, Source_pos_err]);
      % solve the dx 
      % H_mic是替代H进行信息增益计算的，只包含了对于麦克风阵列的信息
      [dx,H,H_mic] = linearize_and_solve_with_H(g_online,est_delay_on,est_drift_on,i);

      % TODO: apply the solution to the state vector g.x
      g_online.x = g_online.x + dx;
%       dx_norms(i) = norm(dx);
      
      % compute the rotation matrix
      rot_yaw = -atan2(g_online.x((g_online.M_x-1)*5+2),g_online.x((g_online.M_x-1)*5+1));
      rot_pitch = atan2(g_online.x((g_online.M_x-1)*5+3),sqrt(g_online.x((g_online.M_x-1)*5+1)^2+g_online.x((g_online.M_x-1)*5+2)^2));
      M_half = transform_matrix_from_trans_ypr(0,0,0,rot_yaw,rot_pitch,0);
      M_y_p_hom = M_half*[g_online.x((g_online.M_y-1)*(g_online.M_x)*5+1:(g_online.M_y-1)*(g_online.M_x)*5+3);1];
      rot_roll = -atan2(M_y_p_hom(3),M_y_p_hom(2));
      M_transform = transform_matrix_from_trans_ypr(0,0,0,rot_yaw,rot_pitch,rot_roll);
      % rotate the mic positions
      for n=2:g_online.M
          g_online.x(5*(n-1)+1:5*(n-1)+3) = [eye(3) zeros(3,1)]*M_transform*[g_online.x(5*(n-1)+1:5*(n-1)+3);1];
      end
      % rotate the sound src positions
      for n=1:(size(g_online.x,1)-5*g_online.M)/3
          g_online.x(5*g_online.M+3*(n-1)+1:5*g_online.M+3*(n-1)+3) = [eye(3) zeros(3,1)]*M_transform*[g_online.x(5*g_online.M+3*(n-1)+1:5*g_online.M+3*(n-1)+3);1];
      end
          
      % display estimation error of mic delay if asked
      if display_delay_error_on > 0    
          x_3_error = (g.x(9:5:g.M*5-1) - g.x_gt(9:5:g.M*5-1));
          disp('estimation error of starting time delay: ');
          x_3_error'
      end

      % TODO: implement termination criterion as suggested on the sheet
      if display_norm_dx_on>0
        disp(['norm(dx) = ' num2str(norm(dx))]);
      end
      
      if (norm(dx)<EPSILON)
        break;
      end


    
    end

% plot the current state of the graph
figure;
plot_graph_with_cov(g_online, i, H);
title(input.fig.title);
view(input.fig.view_a, input.fig.view_e);
end




