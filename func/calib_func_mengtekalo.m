function [sample_point,final_Mic_pos_err,final_Source_pos_err] = calib_func_mengtekalo(input)
% This file performs least square SLAM
% 
%% parameters

% gonna estimate clock drift?
est_drift_on = 1;
% gonna estimate starting time delay?
est_delay_on = 1;
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

sample_num = 100 ;%菱形50,螺旋线50
decrease_threshold = -1;


for i = 1 : (length(g.edges))
    g.edges(i).toIdx_original = g.edges(i).toIdx;
end
%用g_online去迭代，计算eigenvalue与信息增益
g_online = struct();
g_online.M_x = g.M_x;
g_online.M_y = g.M_y;
g_online.M_z = g.M_z;
g_online.M = g.M;
g_online.mic_dis = g.mic_dis;
g_online.x_gt = g.x_gt(1:(5 * g.M + 3 * sample_num),:); 
g_online.x = g.x(1:(5 * g.M + 3 * sample_num),:);
g_online.edges = g.edges(1:(2 * sample_num - 1)); 
g_online.idLookup = g.idLookup(1:(g.M + sample_num)); 

batch = 2;%一批处理两个点
%% start slSLAM

%% 从头开始向后迭代


initial_sample_num = 2;

current_x_index = 5 * g.M + 3 * initial_sample_num; 
current_edge_index = 2 * initial_sample_num;
current_idLook_index = g.M + initial_sample_num;

for a = 0:(  floor(sample_num / batch) - 2) %这里从0开始是为了初始化

    if a > 0 
        disp(a);
        x_gt_temp1 = g_online.x_gt;
        x_temp1 = g_online.x;
        edges_temp1 = g_online.edges;
        idLookup_temp1 = g_online.idLookup;
%         x_gt_temp = g_online.x_gt( (current_x_index + 1):(current_x_index+3*batch),: );
%         x_temp = g_online.x( (current_x_index + 1):(current_x_index+3*batch),: );
%         edges_temp = g_online.edges( (current_edge_index+1):(current_edge_index+2*batch));
%         idLookup_temp = g_online.idLookup( (current_idLook_index+1):(current_idLook_index+batch));
        g_online.x_gt( (current_x_index + 1):(current_x_index+3*batch),: ) = [];
        g_online.x( (current_x_index + 1):(current_x_index+3*batch),: ) = [];
        g_online.edges( (current_edge_index+1):(current_edge_index+2*batch)) = [];
        g_online.idLookup( (current_idLook_index+1):(current_idLook_index+batch)) = [];
        %删去以后需要改 idLookup 和 edge的p测量值，ptoIdx,pfromIxd，和l的toIdx
        for i = (current_idLook_index + 1):length(g_online.idLookup) 
            g_online.idLookup(i).offset = g_online.idLookup(i).offset - 3 * batch;
        end
%         g_online.idLookup( (current_idLook_index+1):end) = g.idLookup( (current_idLook_index+1):(end - batch));

        delta_POSE = g_online.x(current_x_index + 1:current_x_index + 3) - g_online.x(current_x_index - 2:current_x_index);
        g_online.edges(current_edge_index).measurement = delta_POSE;

        for i = (current_edge_index+1):length(g_online.edges)
            g_online.edges(i).toIdx = g_online.edges(i).toIdx - 3 * batch;
        end    

        for i = (current_edge_index+2):2:(length(g_online.edges)-1)
            g_online.edges(i).fromIdx = g_online.edges(i).fromIdx - 3 * batch;
        end
 

        %不管是不是删掉，都要往下把gx给赋值


        disp(current_x_index);
    end

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

    %% 这里判断是否提供了足够的信息增益
    [Mic_pos_err, Source_pos_err] = compute_RMS_error(g_online);
    if (isnan(Mic_pos_err)) || (Mic_pos_err > 1)
    break;
    end
    information_gain = 0;
    % a=0用来初始化
    if a == 0
        last_eigValues = eig(full(H_mic));
    else
        now_eigValues = eig(full(H_mic));

        for j = 1:length(now_eigValues)-5 -14 % 始终有五个特征值为0，对应原点麦克风
        information_gain = information_gain + log(now_eigValues(j)) ;
        end
        for j = 1:length(last_eigValues)-5-14 %没有对两个参数进行估计，每个mic只有三个参数估计，需要再减2*7
          information_gain = information_gain  - log(last_eigValues(j));
        end
        % 判断信息增益，如果删去以后降低了很多，就保留该批
        if information_gain < decrease_threshold

            g_online.x_gt = x_gt_temp1;%多加2个点(一批)的值
            g_online.x = x_temp1;
            g_online.edges = edges_temp1;
            g_online.idLookup = idLookup_temp1;

            current_x_index = current_x_index + 3 * batch;
            current_edge_index = current_edge_index + 2 * batch;
            current_idLook_index = current_idLook_index + batch;            

        else
             last_eigValues = now_eigValues;

        end
        disp(information_gain);

    end

if a == (  floor(sample_num / batch) - 2)
    num = size(g_online.x,1);
    sample_point = (num - 40)/3;
    disp(sample_point);%一共250个点
    final_Mic_pos_err = Mic_pos_err;
    final_Source_pos_err = Source_pos_err;
end
% plot the current state of the graph
% figure;
% plot_graph_with_cov(g_online, i, H);
% title(input.fig.title);
% view(input.fig.view_a, input.fig.view_e);
end
num = size(g_online.x,1);
sample_point = (num - 40)/3;
final_Mic_pos_err = Mic_pos_err;
final_Source_pos_err = Source_pos_err;
end

