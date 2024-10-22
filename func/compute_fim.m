function [J_G, J1, J2, T, J1_only_offset, J1_only_drift, J1_no_offset_drift, ...
    J_G_only_offset, J_G_only_drift, J_G_no_offset_drift, ...
    FIM_only_offset, FIM_only_drift, FIM_no_offset_drift,  ...
    J_G_col_num, J_G_rank,FIM,FIM_eigs,min_FIM_eig,rank_deficiency,min_eigen] = compute_fim(g)

J_G = [];
W_inv = [];

rank_deficiency = [];
min_eigen = [];

J_G_rank1 =[];%自己加的
J_G_row_num1 = [];
J_G_col_num1 = [];

for eid = 1:length(g.edges)
  % 把所有的麦克风和声源(robot pose)看成图中的节点，对每个节点进行遍历，判断边的约束类型并进行分别处理
  % disp([num2str(eid),'/',num2str(length(g.edges))])
  edge = g.edges(eid);

  % pose-pose constraint
  if (strcmp(edge.type, 'P') ~= 0) % edge.type是P，P时返回1，不等于0
    % edge.fromIdx and edge.toIdx describe the location of
    % the first element of the pose in the state vector
    % You should use also this index when updating the elements
    % of the H matrix and the vector b.
    % edge.measurement is the measurement
    % edge.information is the information matrix
    x1 = g.x_gt(edge.fromIdx:edge.fromIdx+2);  % the first robot pose 这里的x是x,y,z(3x1),代表了sk和sk-1，在下面的if中会有不同的含义
    x2 = g.x_gt(edge.toIdx:edge.toIdx+2);      % the second robot pose 看完这个函数去看图是怎么编写的

    % Computing the error and the Jacobians
    % e the error vector
    % A Jacobian wrt x1
    % B Jacobian wrt x2
    [e, A, B] = linearize_pose_pose_constraint(x1, x2, edge.measurement);% J1对pose的偏导为0，J2对pose的偏导为正负I，与pose-landmark一起可以计算出总的J
    
    % compute and add the term to H and b
    J_G = [J_G, zeros(size(J_G,1),3);...
           zeros(3,size(J_G,2)+3)]; % 在右侧添加3列全0列，底部添加3行全0行，保证增加约束后维度正确
    J_G(end-2:end, end-5:end) = [A,B]; %这一步对应这17式，一行一行地对J进行添加，这步是添加倒数第二行，下一个if是对最后一行添加Lk和Tk，p-p约束和p-l约束交叉排列，使代码实现非常优雅
    
    W_inv = [W_inv, zeros(size(W_inv,1),3);...
             zeros(3,size(W_inv,2)), edge.information];

  % pose-landmark constraint
  elseif (strcmp(edge.type, 'L') ~= 0)
    % edge.fromIdx and edge.toIdx describe the location of
    % the first element of the pose and the landmark in the state vector
    % You should use also this index when updating the elements
    % of the H matrix and the vector b.
    % edge.measurement is the measurement
    % edge.information is the information matrix
    x1 = g.x_gt(edge.toIdx:edge.toIdx+2);   % the robot pose,相当于声源的移动
    x2 = g.x_gt(edge.fromIdx:edge.fromIdx+(5*g.M-1));     % the landmark,麦克风阵列，8个阵列，x2取1-40

    % Computing the error and the Jacobians
    % e the error vector
    % A Jacobian wrt x1
    % B Jacobian wrt x2
    % 这里的edge.measurement是对于每个时刻的声源的TDOA，g.M有8个，故有7个TDOA；写一个声源的运动，用TDOA的公式加噪声生成一下；
    [e, A, B] = linearize_pose_landmark_constraint(x1, x2, edge.measurement,edge.toIdx,1,1,g); % 这个函数是用来计算A：H/Lk；B：Tk的
    
    A = A(:,6:end);
    
    % compute and add the term to H and b
    if isempty(J_G)
        J_G = zeros(g.M-1,size(A,2)+size(B,2));
    else
        J_G = [J_G;...
               zeros(g.M-1,size(J_G,2))]; % 这里添加一行全零行，不用加列
    end
    J_G(end-(g.M-1)+1:end,1:5*(g.M-1)) = A; % Lk 对于Fig4 797x275, 797 = (N-1)K+3(K-1)=10k-3,K=80,N=8 
    J_G(end-(g.M-1)+1:end,end-3+1:end) = B; % Tk 
    
    W_inv = [W_inv, zeros(size(W_inv,1),g.M-1);...
             zeros(g.M-1,size(W_inv,2)), edge.information];
         
%     if isempty(rank_deficiency)
      J_G_col_num = size(J_G,2);
      J_G_rank = rank(J_G);

      J_G_rank1 = [J_G_rank1,J_G_rank];%自己加的,判断什么时候rank_deficiency不受到行数的影响；从K=5到K=21时，列不满秩
      J_G_row_num1 = [J_G_row_num1,size(J_G,1)];
      J_G_col_num1 = [J_G_col_num1,size(J_G,2)];

      rank_deficiency = [rank_deficiency,J_G_col_num - J_G_rank];  %Fig4时K=21时为0
      
      FIM = J_G'*W_inv*J_G;
      FIM_eigs = eig(FIM);
      min_eigen = [min_eigen, norm(min(FIM_eigs))];
%     end

  end
  
  
  
end

J_G_col_num = size(J_G,2);
J_G_rank = rank(J_G);

FIM = J_G'*W_inv*J_G;
FIM_eigs = eig(FIM);

min_FIM_eig = min(FIM_eigs);

% situation in which only time offset / clock drift or nothing is present
J1 = J_G(:,1:5*(g.M-1));
J2 = J_G(:,5*(g.M-1)+1:end);

J1_only_offset = zeros(size(J1,1),size(J1,2)-(g.M-1));
J1_only_drift = zeros(size(J1,1),size(J1,2)-(g.M-1));
J1_no_offset_drift = zeros(size(J1,1),size(J1,2)-2*(g.M-1));
for n=1:g.M-1
    J1_only_offset(:,4*(n-1)+1:4*n) = J1(:,5*(n-1)+1:5*(n-1)+4);
    J1_only_drift(:,4*(n-1)+1:4*n) = [J1(:,5*(n-1)+1:5*(n-1)+3),J1(:,5*(n-1)+5)];
    J1_no_offset_drift(:,3*(n-1)+1:3*n) = J1(:,5*(n-1)+1:5*(n-1)+3);
end

J_G_only_offset = [J1_only_offset, J2];
J_G_only_drift = [J1_only_drift, J2];
J_G_no_offset_drift = [J1_no_offset_drift, J2];

% T matrix
T = zeros(size(J2,1),3);
K = size(J2,2)/3;
for k=1:K
    T = T+J2(:,3*(k-1)+1:3*(k-1)+3);
end


FIM_only_offset = J_G_only_offset'*W_inv*J_G_only_offset;
FIM_only_drift = J_G_only_drift'*W_inv*J_G_only_drift;
FIM_no_offset_drift = J_G_no_offset_drift'*W_inv*J_G_no_offset_drift;

end