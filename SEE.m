% Single edge embedding on a network of n nodes and h hops
%function []= SEE()

% Functions used: node, path array, edge embed, OMP, OMP modified, 
%valid path, packet

    clc;clear;
    n = 7; % no. of nodes
    % Define path to be travelled
    %path = [2 4 1 3 n];
    hop_len = [3 4 5 6 7]; % Required hop length
    h = hop_len(2); % choose hop length according to path
    %h = length(path)- 1;
    
    
    % Form Topology
    % Node distribution
    
    %{
    l = 100; b = 100; %length breadth
    points = [];
    for i=1:n
        points = [points; [randi(l,1,1) randi(b,1,1)]];
    end
    
    % Distance between coordinates: Adjacency Matrix
    Adj = squareform(pdist(points,'euclidean'));
    G = graph(Adj);

    [T,pred] = minspantree(G); % Form MST using graph G
    maxweight = max(T.Edges.Weight); % Max edge weight of MST
    reachable_adj_mat = Adj; 

    % Form connected graph using max edge weight
    %reachable_adj_mat(reachable_adj_mat <= maxweight+10 & reachable_adj_mat >0) = 1;
    %reachable_adj_mat(reachable_adj_mat > maxweight+10) = 0;
    %graph_adj_mat = reachable_adj_mat;
    %}
 graph_adj_mat =  [0     0     0     1     1     1     1;
     0     0     0     1     0     1     0;
     0     0     0     1     0     0     1;
     1     1     1     0     0     1     0;
     1     0     0     0     0     0     1;
     1     1     0     1     0     0     0;
     1     0     1     0     1     0     0];
    connected_G = graph(graph_adj_mat);

    %graph_adj_mat = ones(n) - diag(ones([1,n])); %adjacency matrix for complete graph
    %connected_G = graph(graph_adj_mat~=0);  %form graph using adj mat
    figure(1);
    plot(connected_G);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Simulation
    % Simulation Parameters
    
    E = (n-1)^2; % Total no of directed edges
    EE = n^2;    % Modified no of directed edges for additional constraint
    
    mu=0;sigma=1; % gaussian parameters
    error_threshold = 200; % Error Rate = error_threshold/no_of_pkts
    no_of_pkts = 10^3; % Total no of pkts txed
    
    
    temp_m = 2*h; 
    m=[]; % form array for various values of no of rows of matrix Ar(m,n) 
    for ii=1:5
        m =[m temp_m+(ii-2)*temp_m/2];
    end
    %m= [9 12 15 18 21 24];
    m=[8 12 16];
    % Define n nodes 
    nodes_1 =[]; nodes_2=[];
    for i= 1:n
        nodes_1 =[nodes_1 node(i,[],graph_adj_mat(i,:))];
        nodes_2 = [nodes_2 node(i,[],graph_adj_mat(i,:))]; 
    end 
    
    % Destination
    dest = nodes_1(n);
    %fprintf("Destination node: %d",dest.Node_id);
    
    
    % Finding optimal path with required hop length
    temp_matrix = graph_adj_mat;
    temp_matrix(temp_matrix == 0) = inf;
    path=[1     6     4     3     7];
    %[path,src] = find_path(temp_matrix,nodes_1,h);
    if isempty(path)
        fprintf('No path available\n\n');
        %CS1();
    else
        fprintf('Path choosen:');disp(path) 
    end
    
    
    src = nodes_1(path(1));
    fprintf('Path choosen:');disp(path) 
    
    % Path array similar to x (y=Ax) for verify
    Path_arr = path_array(path,E,n); %change
    %Path_arrE = path_array1(path,EE,n); % path func for n^2
    %fprintf("\npath array:");
    %disp(Path_arr)
    
    fprintf("\n Single edge embedding\n");
    
    error_rate = zeros(1,length(m));    % SE CVX
    error_rate_OMP = zeros(1,length(m));    % SE OMP
    error_rate_OMP_mod_adj = zeros(1,length(m));    % SE modified OMP with adjacency matrix
    error_rate_OMP_mod = zeros(1,length(m));    % SE modified OMP with brute force
    error_rate_OMP_topo = zeros(1,length(m));   % SE OMP with topology knowledge
    error_rate_gOMP = zeros(1,length(m));       % SE generalised OMP
    error_rate_gOMP_mod_adj = zeros(1,length(m)); % SE gOMP with list and adjacency matrix path constraint
    %error_rate_cosamp=zeros(1,length(m));
    %error_rate_stOMP=zeros(1,length(m));
    
    for m_index =1: length(m)
        error_count = 0;    
        error_count_OMP=0;  
        error_count_OMP_mod_adj=0; 
        error_count_OMP_mod=0; 
        error_count_OMP_topo =0;
        error_count_gOMP=0;
        error_count_gOMP_mod_adj=0;
       % error_count_stOMP=0;
       % error_count_cosamp=0;
        error_rate(m_index) = 0;
        error_rate_OMP(m_index) = 0;
        error_rate_OMP_mod_adj(m_index) = 0;
        error_rate_OMP_mod(m_index) = 0;
        error_rate_OMP_topo(m_index)=0;
        error_rate_gOMP(m_index)=0;
        error_rate_gOMP_mod_adj(m_index)=0;
      %  error_rate_stOMP(m_index)=0;
      %  error_rate_cosamp(m_index)=0;
    
        pkt_count = 0;

        sigma=1/m(m_index);
        for pkt_i=1:no_of_pkts

            pkt =  Packet;
            pkt_count = pkt_count +1;
            Ar = normrnd(mu,sigma,[m(m_index),E]); %change: for n^2

            % Assigning columns to edges
            Edge_id = reshape(Ar,[m(m_index),E/(n-1),n-1]);% change /(n-1),n-1]);

            % Distributing edge ids to nodes
            for ni= 1:n
                if ni==n
                    edge=Ar; %change
                    %edge = Edge_id(:,:,ni);
                else
                    edge=Edge_id(:,:,ni);
                end
                nodes_1(ni).Edge_id = edge;
            end

            % Embed edge ids
            pkt = edge_embed(m(m_index),path,pkt,nodes_1); %change
            y = pkt.provenance;
            %fprintf("Final provenance\n");disp(y)

            % Recovery using OMP
            x_OMP = OMP(h,y,Ar);
            %{
            for k=1:length(x_OMP)
                if abs(x_OMP(k))<=0.001
                    rec_x_OMP(k)=0;
                else
                    rec_x_OMP(k)=1;
                end
            end
            %}
            
            x_OMP(abs(x_OMP)<=0.001)=0;
            x_OMP(abs(x_OMP)>0.001)=1;
            
            %{
            if rec_x_OMP' ~= x_OMP
                fprintf("bug")
                break
            end
            %}
            
            %if rec_x_OMP' == Path_arr
            if ~isequal(x_OMP, Path_arr)
               error_count_OMP = error_count_OMP +1; % Increment count
            end

            
            % Recovery using gOMP
            N = floor(min(h,m(m_index)/h));
            x_gOMP = gOMP(h,y,Ar,N);
            x_gOMP(abs(x_gOMP)<=0.001)=0;
            x_gOMP(abs(x_gOMP)>0.001)=1;
            if ~isequal(x_gOMP, Path_arr)
               error_count_gOMP = error_count_gOMP +1; % Increment count
            end
            
            %{
            % Recovery using StOMP
            delta = m(m_index)/E;
            rho = h/m(m_index);
            S = 10;
            alpha_0 = delta*(1-rho)/S;
            q = min((m(m_index)-h)/h,0.5);

            [x_stOMP, iters] = SolveStOMP(Ar, y, E, 'FAR', alpha_0, S, 1);
            
            
            
            %x_stOMP = StOMP(h,y,Ar);
            x_stOMP(abs(x_stOMP)<=0.001)=0;
            x_stOMP(abs(x_stOMP)>0.001)=1;
            if ~isequal(x_stOMP, Path_arr)
               error_count_stOMP = error_count_stOMP +1; % Increment count
            else
           %     fprintf("bug")
            end
            
            
            % Recovery using CoSaMP
            x_cosamp = CoSaMP(h,y,Ar);
            x_cosamp(abs(x_cosamp)<=0.001)=0;
            x_cosamp(abs(x_cosamp)>0.001)=1;
            if ~isequal(x_cosamp, Path_arr)
               error_count_cosamp = error_count_cosamp +1; % Increment count
            end
            
            
            %}
            
            
            
            % Recovery using modified OMP
            x_OMP_mod = OMP_modified(h,y,Ar);
            x_OMP_mod(abs(x_OMP_mod)<=0.001)=0;
            x_OMP_mod(abs(x_OMP_mod)>0.001)=1;
            x_OMP_mod_adj = x_OMP_mod;
            orig_x_OMP_mod = x_OMP_mod;
           
            % verifying using adjacency (n-1)^2 -> n^2
            B=[];
            for ii=1:n
                count=0;
                for j=1:EE
                    k=n*(ii-1)+1;
                    count=k+n;
                    if j>=k && j<count 
                        B(ii,j)=1;
                    else 
                        B(ii,j)=0;
                    end
                end
            end
            C=[];
            for ii=1:n
                C=[C;eye(n)];
            end
            
            
            % Extra indices 
            temp_e=zeros(1,n);
            temp_f=zeros(1,n-1);
            temp_g=zeros(1,n-1);
            
            for i=1:n
                temp_e(i)=(i-1)*n + i;
            end
            
            for i=1:n-1
                temp_f(i)=i*n;
                temp_g(i)= n*(n-1)+i;
            end
            
            temp_dest=zeros(1,n-1); % indices mapped to i->dest
            for i=1:(n-1)
                temp_dest(i)=i*(n-1);
            end
            
            ts_adj = tic;       % Time elapse check
            
            jj=[];              % Invalid path indices from set of solutions
            for j=1:size(x_OMP_mod_adj,2)
                %index = find(x_OMP_mod_adj(:,j)==1);
                %if size(index,1)~=h
                 %     jj=[jj j];
                %else      
                    
                    x_nminus1_sq = x_OMP_mod_adj(:,j);
                    %{
                    x_i2_dest = x_nminus1_sq(temp_dest);
                    x_nminus1_sq(temp_dest)=[];
                    x_n_sq = zeros(n^2,1);
                    u=1;v=1;
                    for k=1:n^2
                        if (ismember(k,temp_e)) || (ismember(k,temp_f))
                             x_n_sq(k)=0;    
                        elseif (ismember(k,temp_g))
                             x_n_sq(k)= x_i2_dest(u);
                             u=u+1;
                        else
                             %fprintf("k=%d",k);
                             %fprintf(" v=%d\n",v);
                             x_n_sq(k)= x_nminus1_sq(v);
                             %v = mod(v+1,length(x_nminus1_sq)+2);
                             v=v+1;
                        end
                        %fprintf("k=%d",k);
                    end
                    %}
                    x_n_sq = embed_extra_zeros(n,x_nminus1_sq,temp_e,temp_f,temp_g,temp_dest);
                    %{
                    A_mat = diag(x_n_sq);
                    Z = B*A_mat*C;
                    nn=zeros(1,h); % array containing i,j element of Z^i,i=1 to h
                    for ii=1:h
                        Za = Z^ii;
                        nn(ii)=Za(end,src.Node_id);
                    end

                    c = zeros(1,h);
                    c(h)=1;
                    %}
                    opt = Adjacency_constraint(x_n_sq,B,C,h,src.Node_id);
                    if opt==1
                       %fprintf("Valid h len path");
                    else
                        %fprintf("Invalid");
                        jj=[jj j];
                    end
                %end
                
            end
            x_OMP_mod_adj(:,jj)=[];
            te_adj(pkt_i) = toc(ts_adj) ;
    %{        
            
            % verifying using valid path func
            ji=[];
            ts_brute = tic;
            for j=1:size(x_OMP_mod,2)
                opt = valid_path(x_OMP_mod(:,j),n,h,src.Node_id);
                if opt==1
                else
                    ji=[ji j];
                end
           end
           x_OMP_mod(:,ji)=[];
           te_brute(pkt_i) = toc(ts_brute) ;
           %display(ji);
           if x_OMP_mod_adj ~= x_OMP_mod
               display(x_OMP_mod_adj);
               display(x_OMP_mod)
           end
      %}     
           if isempty(x_OMP_mod_adj)
                rec_x_OMP_adj = orig_x_OMP_mod(:,1);
           else
                rec_x_OMP_adj = x_OMP_mod_adj(:,1);
           end
        %{    
            if isempty(x_OMP_mod)
                rec_x_OMP_valid = orig_x_OMP_mod(:,1);
            else
                rec_x_OMP_valid = x_OMP_mod(:,1);
            end 
          %} 
            if ~isequal(rec_x_OMP_adj, Path_arr)
               error_count_OMP_mod_adj = error_count_OMP_mod_adj +1; % Increment count
            end
            %{
            if ~isequal(rec_x_OMP_valid , Path_arr)
               error_count_OMP_mod = error_count_OMP_mod +1; % Increment count
            end
            %}
            
            % Recovery using modified gOMP
            lastwarn('','');
            N = floor(min(h,m(m_index)/h));
            x_gOMP_mod = gOMP_modified(h,y,Ar,N);
            x_gOMP_mod(abs(x_gOMP_mod)<=0.001)=0;
            x_gOMP_mod(abs(x_gOMP_mod)>0.001)=1;
            x_gOMP_mod_adj = x_gOMP_mod;
            orig_x_gOMP_mod = x_gOMP_mod;
            jj=[];              % Invalid path indices from set of solutions
            for j=1:size(x_gOMP_mod_adj,2)
                    x_nminus1_sq = x_gOMP_mod_adj(:,j);
                    x_n_sq = embed_extra_zeros(n,x_nminus1_sq,temp_e,temp_f,temp_g,temp_dest);
                    
                    opt = Adjacency_constraint(x_n_sq,B,C,h,src.Node_id);
                    if opt==0
                        jj=[jj j];
                    end
            end
            x_gOMP_mod_adj(:,jj)=[];
            if isempty(x_gOMP_mod_adj)
                rec_x_gOMP_adj = orig_x_gOMP_mod(:,1);
            else
                rec_x_gOMP_adj = x_gOMP_mod_adj(:,1);
            end
            if ~isequal(rec_x_gOMP_adj, Path_arr)
               error_count_gOMP_mod_adj = error_count_gOMP_mod_adj +1; % Increment count
            end
            
            
            
            % Recovery using OMP when topology is known
            % Selecting columns of Ar using topology adjacency matrix
            temp_k=[];
            kk=1;
            for i=1:(n-1)
                for j=1:n
                    if i==j
                    elseif graph_adj_mat(i,j)==1
                        temp_k=[temp_k kk];
                        kk=kk+1;
                    else
                        kk=kk+1;
                    end
                end
            end
            A_topo = Ar(:,temp_k);
            x_topo = OMP_modified(h,y,A_topo);
            x_topo(abs(x_topo)<=0.001)=0;
            x_topo(abs(x_topo)>0.001)=1;
            
            % Embed zeros for x_topo to x_(n-1)^2
            
            x_OMP_top=zeros(E,size(x_topo,2));
            x_check = x_OMP_top;
            temp_check = temp_k;
            for uu=1:size(x_topo,2)
                x_check(temp_check,uu) = x_topo(:,uu); 
                kk=1;ki=1;  % Count variables
                for i=1:(n-1)
                    for j=1:n
                        if i==j
                        elseif graph_adj_mat(i,j)==1
                            x_OMP_top(kk,uu)=x_topo(ki,uu);
                            kk=kk+1;
                            ki=ki+1;
                        else
                             x_OMP_top(kk,uu)=0;
                             kk=kk+1;
                        end
                    end
                end
            end
            if ~isequal(x_check, x_OMP_top)
                fprtintf("Bug")
                break
            end
            orig_x_OMP_top = x_OMP_top;
            ji=[];
            for j=1:size(x_OMP_top,2)
                x_nminus1_sq = x_OMP_top(:,j);
                x_n_sq = embed_extra_zeros(n,x_nminus1_sq,temp_e,temp_f,temp_g,temp_dest);
                opt1 = Adjacency_constraint(x_n_sq,B,C,h,src.Node_id);
              %  opt2=valid_path(x_OMP_top(:,j),n,h,src.Node_id);
                    if opt1==1
                       %fprintf("Valid h len path");
                       
                    %elseif opt2~=opt1
                    %   fprintf("BUg");
                    %   break;
                    else
                        %fprintf("Invalid");
                         ji=[ji j];
                    end    
                
            end
            x_OMP_top(:,ji)=[];
           
           
            if isempty(x_OMP_top)
                rec_x_OMP_top = orig_x_OMP_top(:,1);
            else
                rec_x_OMP_top = x_OMP_top(:,1);
            end
            
            if ~isequal(rec_x_OMP_top, Path_arr)
               error_count_OMP_topo = error_count_OMP_topo +1; % Increment count
            end
            
            
            
            
            
            
            
            
            
            
            
            
            % Recovery using CVX
            %x = cvx_solver(E,y,Ar);
           %{ 
            cvx_begin quiet
                cvx_precision default
                variable x(E) %binary  %change
                minimize( norm(x, 1 ) )
                subject to
                    Ar*x == y;
                    
            cvx_end
            %{
            for k=1:length(x)
                if abs(x(k))<=0.001
                    rec_x(k)=0;
                else
                    rec_x(k)=1;
                end
            end
            %}
            
            x(abs(x)<=0.001)=0;
            x(abs(x)>0.001)=1;
           
            % Comparing recovered path wth original path
            if isequal(x, Path_arr)
                %fprintf("Path matched\n")
            else
                error_count = error_count+1; %Increment if path recovered is different from path travelled
            end
             %}
            if error_count_OMP_topo == error_threshold % when threshold reached
                break
            end
        end
        fprintf("size gomp list:");display(size(orig_x_gOMP_mod));
        fprintf("size omp list:") ;display(size(orig_x_OMP_mod));
        t_Adj(m_index) = mean(te_adj);
%        t_Brute(m_index) = mean(te_brute);
        error_rate(m_index) = error_count/pkt_count;
        error_rate_OMP(m_index) = error_count_OMP/pkt_count;
        error_rate_OMP_topo(m_index) = error_count_OMP_topo/pkt_count;
        error_rate_OMP_mod(m_index) =  error_count_OMP_mod/pkt_count;
        error_rate_OMP_mod_adj(m_index) =  error_count_OMP_mod_adj/pkt_count;
        error_rate_gOMP(m_index) =  error_count_gOMP/pkt_count;
        error_rate_gOMP_mod_adj(m_index)=error_count_gOMP_mod_adj/pkt_count;
       % error_rate_stOMP(m_index) = error_count_stOMP/pkt_count;
       % error_rate_cosamp(m_index) = error_count_cosamp/pkt_count;
        %fprintf("Error Rate:%f for column size %d\n",error_rate(m_index),m(m_index));
        fprintf("Error Rate:%f for column size OMP %d\n",error_rate_OMP(m_index),m(m_index));
        fprintf("Error Rate:%f for column size OMP topo %d\n",error_rate_OMP_topo(m_index),m(m_index));
        %fprintf("Error Rate:%f for column size mod OMP %d\n",error_rate_OMP_mod(m_index),m(m_index));
        fprintf("Error Rate:%f for column size mod OMP adj %d\n",error_rate_OMP_mod_adj(m_index),m(m_index));
        fprintf("Error Rate:%f for column size gOMP  %d\n",error_rate_gOMP(m_index),m(m_index));
        fprintf("Error Rate:%f for column size mod gOMP  %d\n",error_rate_gOMP_mod_adj(m_index),m(m_index));
        %fprintf("Error Rate:%f for column size stomp %d\n",error_rate_stOMP(m_index),m(m_index));
        %fprintf("Error Rate:%f for column size cosamp  %d\n",error_rate_cosamp(m_index),m(m_index));
    end
    
      
    
    
%% PLOTS    
    
    
    figure(2)
    %semilogy(m,error_rate, 'mo-', 'LineWidth', 2);
    
    semilogy(m,error_rate_OMP, 'bs-', 'LineWidth', 2);
    hold on
    semilogy(m,error_rate_gOMP, 'ms-', 'LineWidth', 2);
    %semilogy(m,error_rate_OMP_mod, 'go-', 'LineWidth', 2);
    semilogy(m,error_rate_OMP_mod_adj, 'gs-', 'LineWidth', 2);
    semilogy(m,error_rate_OMP_topo, 'rs-', 'LineWidth', 2);
    semilogy(m,error_rate_gOMP_mod_adj, 'ks-', 'LineWidth', 2);
    %semilogy(m,error_rate_stOMP, 'yo-', 'LineWidth', 2);
    %semilogy(m,error_rate_cosamp, 'ko-', 'LineWidth', 2);
    
    axis([0 120 0 1]);
    
    grid on
    legend('SE OMP','SE gOMP' ,'SE modified OMP adjacency','SE modified OMP adjacency & topology knowledge','SE mod gomp')
    %legend('SE OMP','SE modified OMP using path func','SE modified OMP adjacency');
    %legend('L=4','L=8','L=12');
    title('Error rate vs number of rows');
    xlabel('Column size');
    ylabel('Error rate')
    %hold off
%end
%DEE(path,n)
infos ={'Results for SEE: (II)1. OMP 2. Modified OMP with path constraints (complete graph) in form of:a) adjacency matrix: DE to SE b)brute force 3. (2.a) with topology knowledge';
    'n= no. of nodes';'h = path hop length/sparsity';'E= no. of possible edges';'m = column size of Ar for OMP';'graph_adj_mat = topology adjacency matrix';'Error rate for different cases'};
%filename = 'SEE_omp.mat';
%save SEE_omp_gomp.mat infos n E h graph_adj_mat path m error_rate_OMP error_rate_OMP_mod_adj error_rate_OMP_topo error_rate_gOMP error_rate_gOMP_mod_adj