function [H,H_mic] = compute_GaussNewton(g_online,est_delay_on,est_drift_on,display_delay_error_on,display_norm_dx_on,EPSILON,numIterations)


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
    
      if display_norm_dx_on>0
        disp(['norm(dx) = ' num2str(norm(dx))]);
      end
      
      if (norm(dx)<EPSILON)
        break;
      end

    end

end