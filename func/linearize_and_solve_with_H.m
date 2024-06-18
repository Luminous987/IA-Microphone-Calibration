% performs one iteration of the Gauss-Newton algorithm
% each constraint is linearized and added to the Hessian

function [dx,H,H_mic] = linearize_and_solve_with_H(g,est_delay_on,est_drift_on,i)

nnz = nnz_of_graph(g);

% allocate the sparse H and the vector b
H = spalloc(length(g.x), length(g.x), nnz);
b = zeros(length(g.x), 1);
H_mic = zeros(5*g.M,5*g.M);
global Cov_e_inv_l_total Cov_e_inv_p_total;
if i == 1
Cov_e_inv_l_total = [];
Cov_e_inv_p_total = [];
end
current_p_index = 1;
current_l_index = 1;
needToAddPrior = true;
abnormal = 0;
% compute the addend term to H and b for each of our constraints
for eid = 1:length(g.edges)
  edge = g.edges(eid);

  % pose-pose constraint
  if (strcmp(edge.type, 'P') ~= 0)
    % edge.fromIdx and edge.toIdx describe the location of
    % the first element of the pose in the state vector
    % You should use also this index when updating the elements
    % of the H matrix and the vector b.
    % edge.measurement is the measurement
    % edge.information is the information matrix
    x1 = g.x(edge.fromIdx:edge.fromIdx+2);  % the first robot pose
    x2 = g.x(edge.toIdx:edge.toIdx+2);      % the second robot pose

    % Computing the error and the Jacobians
    % e the error vector
    % A Jacobian wrt x1
    % B Jacobian wrt x2
    [e, A, B] = linearize_pose_pose_constraint(x1, x2, edge.measurement);
    
    
%     if i == 1
%         Cov_e_p = (e * e');
%         Cov_e_inv_p = inv(Cov_e_p); 
%         Cov_e_inv_p_total = [Cov_e_inv_p_total;Cov_e_inv_p];
%     end
%     
%     % 替换原有的信息矩阵
%     edge.information = Cov_e_inv_p_total(current_p_index:current_p_index+2,:);
%     current_p_index = current_p_index + 3;

        Cov_e_p = zeros(3,3);
        lambda_p = 1 * 1e-7;
        adjusted_p_lambda = 0;
        max_p_error = 0.5;

        for i = 1 : 3
            if abs(e(i,1)) > max_p_error
                adjusted_p_lambda = lambda_p * 500;
            else
                adjusted_p_lambda = lambda_p;
            end

            Cov_e_p(i,i) = e(i,1) * e(i,1); 
        end
        Cov_e_inv_p = inv(Cov_e_p + adjusted_p_lambda * eye(size(Cov_e_p))); 
        edge.information = 2 * Cov_e_inv_p / 10^-1;
          

    % TODO: compute and add the term to H and b
    b(edge.fromIdx:edge.fromIdx+2) = (b(edge.fromIdx:edge.fromIdx+2)' + (e')*edge.information*A)';
    b(edge.toIdx:edge.toIdx+2) = (b(edge.toIdx:edge.toIdx+2)' + (e')*edge.information*B)';

    H(edge.fromIdx:edge.fromIdx+2,edge.fromIdx:edge.fromIdx+2) = H(edge.fromIdx:edge.fromIdx+2,edge.fromIdx:edge.fromIdx+2) + A'*edge.information*A;
    H(edge.fromIdx:edge.fromIdx+2,edge.toIdx:edge.toIdx+2) = H(edge.fromIdx:edge.fromIdx+2,edge.toIdx:edge.toIdx+2) + A'*edge.information*B;
    H(edge.toIdx:edge.toIdx+2,edge.fromIdx:edge.fromIdx+2) = H(edge.toIdx:edge.toIdx+2,edge.fromIdx:edge.fromIdx+2) + B'*edge.information*A;
    H(edge.toIdx:edge.toIdx+2,edge.toIdx:edge.toIdx+2) = H(edge.toIdx:edge.toIdx+2,edge.toIdx:edge.toIdx+2) + B'*edge.information*B;

    

  % pose-landmark constraint
  elseif (strcmp(edge.type, 'L') ~= 0)
    % edge.fromIdx and edge.toIdx describe the location of
    % the first element of the pose and the landmark in the state vector
    % You should use also this index when updating the elements
    % of the H matrix and the vector b.
    % edge.measurement is the measurement
    % edge.information is the information matrix
    x1 = g.x(edge.toIdx:edge.toIdx+2);   % the robot pose
    x2 = g.x(edge.fromIdx:edge.fromIdx+(5*g.M-1));     % the landmark

    % Computing the error and the Jacobians
    % e the error vector
    % A Jacobian wrt x1
    % B Jacobian wrt x2
    [e, A, B] = linearize_pose_landmark_constraint(x1, x2, edge.measurement,edge.toIdx,est_delay_on,est_drift_on,g);

%     if i == 1
%         Cov_e_l = (e * e');
%         Cov_e_inv_l = inv(Cov_e_l); 
%         Cov_e_inv_l_total = [Cov_e_inv_l_total; Cov_e_inv_l];
%     end
%   
%    edge.information = Cov_e_inv_l_total(current_l_index:current_l_index+6,:);
%    current_l_index = current_l_index + 7;

        lambda_l = 1 * 1e-9;
        max_error = 0.5;
        Cov_e_l = zeros(7,7);
        for i =1 : 7
%             if abs(e(i,1)) > max_error
%                 e(i,1) = sign(e(i,1)) * max_error * 0;
%             end
            Cov_e_l(i,i) = e(i,1) * e(i,1);
        end
%         Cov_e_l = (e * e');
        Cov_e_inv_l = inv(Cov_e_l + lambda_l * eye(size(Cov_e_l))) ;
        edge.information = 0.5 * Cov_e_inv_l / 10^-2;




%     % 调整信息矩阵，这里设置阈值为1，可根据具体情况调整
%     threshold = 1*10^-4;
%     edge.information = adjustInformationMatrix(e, edge.information, threshold);
%         

%     % 应用Huber核调整
%     W = zeros(length(e), length(e));  % 初始化权重矩阵
%         for i = 1:length(e)
%             [~, drho] = huber(e(i), k);
%             W(i,i) = drho / e(i);  % 当e(i)非零时安全
%             e(i) = drho;           % 更新误差向量
%         end
%         
%     % 更新雅可比矩阵
%     A = W * A;
%     B = W * B;
    % 这里的J是论文中17式，写作分块矩阵J = [A B],H=J*J^T
    % compute and add the term to H and b
    b(edge.fromIdx:edge.fromIdx+(5*g.M-1)) = (b(edge.fromIdx:edge.fromIdx+(5*g.M-1))' + (e')*edge.information*A)'; % b += JΩe,Ω是information
    b(edge.toIdx:edge.toIdx+2) = (b(edge.toIdx:edge.toIdx+2)' + (e')*edge.information*B)';

    H(edge.fromIdx:edge.fromIdx+(5*g.M-1),edge.fromIdx:edge.fromIdx+(5*g.M-1)) = H(edge.fromIdx:edge.fromIdx+(5*g.M-1),edge.fromIdx:edge.fromIdx+(5*g.M-1)) + A'*edge.information*A;
    H(edge.fromIdx:edge.fromIdx+(5*g.M-1),edge.toIdx:edge.toIdx+2) = H(edge.fromIdx:edge.fromIdx+(5*g.M-1),edge.toIdx:edge.toIdx+2) + A'*edge.information*B;
    H(edge.toIdx:edge.toIdx+2,edge.fromIdx:edge.fromIdx+(5*g.M-1)) = H(edge.toIdx:edge.toIdx+2,edge.fromIdx:edge.fromIdx+(5*g.M-1)) + B'*edge.information*A;
    H(edge.toIdx:edge.toIdx+2,edge.toIdx:edge.toIdx+2) = H(edge.toIdx:edge.toIdx+2,edge.toIdx:edge.toIdx+2) + B'*edge.information*B;
    % H += J^TΩJ 用分块矩阵计算
    H_mic(edge.fromIdx:edge.fromIdx+(5*g.M-1),edge.fromIdx:edge.fromIdx+(5*g.M-1)) = H_mic(edge.fromIdx:edge.fromIdx+(5*g.M-1),edge.fromIdx:edge.fromIdx+(5*g.M-1)) + A'*edge.information*A;
  end

end

if (needToAddPrior)
  % add the prior for one pose of this edge
  % This fixes one node to remain at its current location

  % 1st mic
  H(1:5,1:5) = eye(5);
  
  % 2nd mic
  H(5*(g.M_x-1)+2,5*(g.M_x-1)+2) = 1;
  H(1:5*(g.M_x-1)+2-1,5*(g.M_x-1)+2) = zeros(5*(g.M_x-1)+2-1,1);
  H(5*(g.M_x-1)+2,1:5*(g.M_x-1)+2-1) = zeros(1,5*(g.M_x-1)+2-1);
  H(5*(g.M_x-1)+3,5*(g.M_x-1)+3) = 1;
  H(1:5*(g.M_x-1)+3-1,5*(g.M_x-1)+3) = zeros(5*(g.M_x-1)+3-1,1);
  H(5*(g.M_x-1)+3,1:5*(g.M_x-1)+3-1) = zeros(1,5*(g.M_x-1)+3-1);
  
  % 3rd mic
  H(5*(g.M_y-1)*g.M_x+3,5*(g.M_y-1)*g.M_x+3) = 1;
  H(1:5*(g.M_y-1)*g.M_x+3-1,5*(g.M_y-1)*g.M_x+3) = zeros(5*(g.M_y-1)*g.M_x+3-1,1);
  H(5*(g.M_y-1)*g.M_x+3,1:5*(g.M_y-1)*g.M_x+3-1) = zeros(1,5*(g.M_y-1)*g.M_x+3-1);

  % if don't need to estimate drift, adjust the information matrix to
  % have no correlation with other state variables
  if est_drift_on<1
      for n = 2:g.M
          H(5*(n-1)+5,1:(5*(n-1)+5-1)) = zeros(1,(5*(n-1)+5-1));
          H(1:(5*(n-1)+5-1),5*(n-1)+5) = zeros((5*(n-1)+5-1),1);
          H(5*(n-1)+5,5*(n-1)+5)=1;
      end
  end

  % if don't need to estimate drift, adjust the information matrix to
  % have no correlation with other state variables
  if est_delay_on<1
      for n = 2:g.M
          H(5*(n-1)+4,1:(5*(n-1)+4-1)) = zeros(1,(5*(n-1)+4-1));
          H(1:(5*(n-1)+4-1),5*(n-1)+4) = zeros((5*(n-1)+4-1),1);
          H(5*(n-1)+4,5*(n-1)+4)=1;
      end
  end
end

% solve the linear system, whereas the solution should be stored in dx
% Remember to use the backslash operator instead of inverting H

%   disp(H);
%   [eigenVectors, eigenValues] = eigs(H,10,'smallestabs');
%   disp('eigenvalues');
%   disp(eigenValues);
    
%   % eigenValues 包含最小的特征值
%   min_eigenvalue = eigenValues(1,1); 
%   disp(min_eigenvalue); % 显示最小特征值

% disp(rank(H));
dx = H\(-b);

end
% function [rho, drho] = huber(r, k)
%     % Huber核函数和它的导数
%     abs_r = abs(r);
%     if abs_r <= k
%         rho = 0.5 * r^2;
%         drho = r;
%     else
%         rho = k * (abs_r - 0.5 * k);
%         drho = k * sign(r);
%     end
% end
function Omega_adj = adjustInformationMatrix(e, Omega, threshold)
    n = length(e);  % 误差向量的长度
    Omega_adj = Omega;  % 初始化调整后的信息矩阵为原信息矩阵
    
    sigma = std(e);  % 计算误差的标准差
    for i = 1:n
        if abs(e(i)) > threshold * sigma  % 如果误差大于阈值乘以标准差
            Omega_adj(i, i) = Omega(i, i) * 0.001;  % 减少该误差项的权重
        else
            Omega_adj(i, i) = Omega(i, i);  % 保持原权重
        end
    end
end
