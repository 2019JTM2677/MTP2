% Double edge embedding
%function []= DEE(path,n,h)
% Functions used: node, path array de, double edge embed, OMP, OMP modified, 
%DE to SE, DE valid path, packet
    n = 5; % no. of nodes
    %path = [2 4 1 3 n];
    hop_len = [3 4 5 6 7]; % Required hop length
    h = hop_len(2);
    %h = length(path)- 1; %  hop length
    
    
    % Form Topology
    % Node distribution
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
    reachable_adj_mat(reachable_adj_mat <= maxweight+10 & reachable_adj_mat >0) = 1;
    reachable_adj_mat(reachable_adj_mat > maxweight+10) = 0;
    graph_adj_mat = reachable_adj_mat;
    connected_G = graph(graph_adj_mat);

    
    
    %graph_adj_mat = ones(n) - diag(ones([1,n]));
    %connected_G = graph(graph_adj_mat~=0);
    
    figure(1);
    plot(connected_G);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Simulation Parameters
    DE = ((n-2)^2)*(n-1); % No. of Double edges
    EE = n^2; 
    
    mu = 0; sigma = 5;  % gaussian parameters

    error_threshold = 200; % Error rate= error_threshold/ no_of_pkts
    no_of_pkts = 10^4; % Total no of pkts txed
    m=[4 8 12 16 20 24];
    %m=[1 3 6 9 11];
    
    % Define n nodes
    nodes_de=[];  % All varE used for additional constraint
    for i= 1:n
        nodes_de = [nodes_de node(i,[],graph_adj_mat(i,:))];
    end

    % Destination
    dest = nodes_de(n);
    %fprintf("Destination node: %d",dest.Node_id);

    % Finding optimal path with required hop length
    temp_matrix = graph_adj_mat;
    temp_matrix(temp_matrix == 0) = inf;
    path=[];
    [path,src] = find_path(temp_matrix,nodes_de,h);
    if isempty(path)
        fprintf('No path available\n\n');
        %CS1();
    else
        fprintf('Path choosen:');disp(path) 
    end
    
    
    %src = nodes_1(path(1));
    
    % Penultimate node
    penultimate = nodes_de(path(h));
    
    % Path array similar to x (y=Ax) for verify
    Path_arr_de = path_array_de(path,DE,n);
    %fprintf("path array for double edge:");
    %disp(Path_arr_de')
    error_rate_de=zeros(1,length(m));   %DE CVX
    error_rate_Ode=zeros(1,length(m));  % DE OMP
    error_rate_OMP_mod=zeros(1,length(m)); % DE modified OMP with brute force
    error_rate_OMP_mod_adj=zeros(1,length(m)); % DE modified OMP with adjacency matrix
    error_rate_OMP_topo = zeros(1,length(m)); % SE OMP with topology knowledge
    
    fprintf("\n Double edge embedding\n");
    for m_index =1: length(m)
        error_count_de=0;
        error_count_Ode=0;
        error_count_OMP_mod=0;
        error_count_OMP_mod_adj=0;
        error_count_OMP_topo =0;
        error_rate_de(m_index) = 0;
        error_rate_Ode(m_index) = 0;
        error_rate_OMP_mod(m_index) = 0;
        error_rate_OMP_mod_adj(m_index) = 0;
        error_rate_OMP_topo(m_index)=0;
        
        pkt_count = 0;

        for pkt_i=1:no_of_pkts
            pkt_de =  Packet;
            Ar_de = normrnd(mu,sigma,[m(m_index),DE]);
            
            % Assigning columns to edges
            Edge_id_de = reshape(Ar_de,[m(m_index),(n-2)^2,n-1]);
            pkt_count = pkt_count +1;

            for ni= 1:n
                if ni==n
                    dbl_edge = Ar_de;
                else
                    dbl_edge = Edge_id_de(:,:,ni);
                end
                nodes_de(ni).Edge_id = dbl_edge;
            end

            % Embed double edge ids
            pkt_de = double_edge_embed(m(m_index),path,pkt_de,nodes_de);
            b_de = pkt_de.provenance;
            %fprintf("Final provenance\n");disp(b_de)

            % Recovery using OMP double edge
            if rem(length(path),2)==0    
                h_omp = length(path)/2 -1;
            else
                h_omp = floor(length(path)/2);
            end
           
            x_Ode = OMP(h_omp,b_de,Ar_de);
            %{
            for k=1:length(x_Ode)
                if abs(x_Ode(k))<=0.001
                    rec_x_Ode(k)=0;
                else
                   rec_x_Ode(k)=1;
                end
            end
            %}
            x_Ode(abs(x_Ode)<=0.001)=0;
            x_Ode(abs(x_Ode)>0.001)=1;
            
            if isequal(x_Ode, Path_arr_de)
               %fprintf("Path matched")
            else
               error_count_Ode = error_count_Ode +1; % Increment count
            end

            % Recovery using modified OMP
            x_OMP_mod = OMP_modified(h_omp,b_de,Ar_de);
            x_OMP_mod(abs(x_OMP_mod)<=0.001)=0;
            x_OMP_mod(abs(x_OMP_mod)>0.001)=1;
            x_OMP_mod_adj = x_OMP_mod;
            orig_x_OMP_mod = x_OMP_mod;
            
            % Verifying using adajacency matrix
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
            
            ts_adj = tic;
            temp_a=zeros(1,n^2);temp_b=zeros(1,n^2);
            temp_c=zeros(1,n^2);temp_d=zeros(1,n^2);
            ii=0;
            for j=1:n
                for i=1:n
                    ii=ii+1;
                    temp_a(ii)= n^2*(j-1)+ n*(j-1)+i;
                    temp_b(ii)= n^2*j -i + 1;
                    temp_c(ii)= n^2*(j-1)+ j + n*(i-1);
                    temp_d(ii)= n^2*(j-1)+ (n+1)*(i-1)+1;

                end
            end
            temp_e = n^2*(n-1)+(1:n^2);
            A_1=zeros(n^2,n^3); %a->b->i
            A_2=zeros(n^2,n^3); %i_>a->b

            for i=1:n^2
                j=(1:n)+n*(i-1);
                A_1(i,j)=1;
            end

            for i=1:n^2
                g=ceil(i/n);
                k= mod(i,n);
                if k==0
                    k=n;
                end
                j=g+(0:n-1)*n+(n^2*(k-1));
                A_2(i,j)=1;
            end

            jj=[];
            for j=1:size(x_OMP_mod_adj,2) 
                x_n_sq = DE_to_SE(x_OMP_mod_adj(:,j),n,h,penultimate.Node_id,temp_a,temp_b,temp_c,temp_d,temp_e,A_1,A_2);
                A_mat = diag(x_n_sq);
                Z = B*A_mat*C;
                nn=zeros(1,h); 
                for hi=1:h
                    Za = Z^hi;
                    nn(hi)=Za(end,src.Node_id);
                end
                c = zeros(1,h);
                c(h)=1;
                if nn==c
                   %fprintf("Valid h len path");
                else
                    %fprintf("Invalid");
                    jj=[jj j];
                end
            end
            x_OMP_mod_adj(:,jj)=[];
            te_adj(pkt_i) = toc(ts_adj) ;
            
            % verifying using valid path func
            ts_brute = tic;
            ji=[];
            for j=1:size(x_OMP_mod,2)
                opt = DE_valid_path(x_OMP_mod(:,j),n,h,src.Node_id,dest.Node_id,penultimate.Node_id);
                if opt==1
                else
                    ji=[ji j];
                end
            end
            x_OMP_mod(:,ji)=[];
            te_brute(pkt_i) = toc(ts_brute) ;

            if isempty(x_OMP_mod_adj)
                rec_x_OMP_adj = orig_x_OMP_mod(:,1);
            else
                rec_x_OMP_adj = x_OMP_mod_adj(:,1);
            end

            if isempty(x_OMP_mod)
                rec_x_OMP_valid = orig_x_OMP_mod(:,1);
            else
                rec_x_OMP_valid = x_OMP_mod(:,1);
            end

            if isequal(rec_x_OMP_adj, Path_arr_de)
            %    fprintf("Path matched\n")
            else
               error_count_OMP_mod_adj = error_count_OMP_mod_adj +1; % Increment count
            end
            
            if isequal(rec_x_OMP_valid, Path_arr_de)
            %    fprintf("Path matched\n")
            else
               error_count_OMP_mod = error_count_OMP_mod +1; % Increment count
            end
            
            if ~isequal(rec_x_OMP_adj,Path_arr_de) && isequal(rec_x_OMP_valid,Path_arr_de)
                fprintf("bug");
                break
                
            end
            
            
            % Recovery using OMP when topology is known
            % Selecting columns of Ar using topology adjacency matrix
            temp_k=[];
            kk=1;
            for i=1:(n-1)
                for j=1:(n-1)
                    if i==j
                    elseif graph_adj_mat(i,j)==1
                        for b=1:n
                            if (b==i) || (b==j)
                                %fprintf("Not embedding:%d %d %d %d",i,j,b,graph_adj_mat(b,i));
                            elseif graph_adj_mat(b,i)==1
                                temp_k=[temp_k kk];
                                kk=kk+1;
                                %fprintf("Embed: %d %d %d %d",i,j,b,graph_adj_mat(b,i));
                            else
                                kk=kk+1;
                            end
                        end
                    else 
                        for l=1:n-2
                            kk=kk+1;
                        end
                    end
                end
            end
            A_topo = Ar_de(:,temp_k);
            x_topo = OMP_modified(h_omp,b_de,A_topo);
            x_topo(abs(x_topo)<=0.001)=0;
            x_topo(abs(x_topo)>0.001)=1;
            
            x_OMP_top=zeros(DE,size(x_topo,2));
            for uu=1:size(x_topo,2)
                x_OMP_top(temp_k,uu) = x_topo(:,uu);
            end
            
            orig_x_OMP_top = x_OMP_top;
            ji=[];
            for j=1:size(x_OMP_top,2) 
                x_n_sq = DE_to_SE(x_OMP_top(:,j),n,h,penultimate.Node_id,temp_a,temp_b,temp_c,temp_d,temp_e,A_1,A_2);
                A_mat = diag(x_n_sq);
                Z = B*A_mat*C;
                nn=zeros(1,h); 
                for hi=1:h
                    Za = Z^hi;
                    nn(hi)=Za(end,src.Node_id);
                end
                c = zeros(1,h);
                c(h)=1;
                if nn==c
                   %fprintf("Valid h len path");
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
            
            if isequal(rec_x_OMP_top, Path_arr_de)
            %    fprintf("Path matched\n")
            else
               error_count_OMP_topo = error_count_OMP_topo +1; % Increment count
            end
            
            % Recovery using CVX double edge
            
            %x_de = cvx_solver(DE,b_de,Ar_de);
            cvx_begin quiet
                cvx_precision default
                variable x_de(DE) %binary
                minimize( norm(x_de, 1 ) )
                subject to
                    Ar_de*x_de == b_de;
            cvx_end
            for k=1:length(x_de)
                if abs(x_de(k))<=0.001
                    rec_x_de(k)=0;
                else
                   rec_x_de(k)=1;
                end
            end

            if rec_x_de' == Path_arr_de
                %fprintf("Path matched\n")
            else
                error_count_de = error_count_de +1; %Increment if path recovered is different from path travelled
            end
            %}
            if error_count_OMP_topo == error_threshold% when threshold reached
                break
            end

        end
        
        t_Adj(m_index) = mean(te_adj);
        t_Brute(m_index) = mean(te_brute);
        error_rate_de(m_index) = error_count_de/pkt_count;
        error_rate_Ode(m_index) = error_count_Ode/pkt_count;
        error_rate_OMP_mod(m_index) =  error_count_OMP_mod/pkt_count;
        error_rate_OMP_mod_adj(m_index) =  error_count_OMP_mod_adj/pkt_count;
        error_rate_OMP_topo(m_index) = error_count_OMP_topo/pkt_count;
        
        fprintf("Error Rate:%f for column size %d\n",error_rate_de(m_index),m(m_index));
        fprintf("Error Rate:%f for column size %d\n",error_rate_Ode(m_index),m(m_index));
        fprintf("Error Rate:%f for column size mod OMP %d\n",error_rate_OMP_mod(m_index),m(m_index));
        fprintf("Error Rate:%f for column size mod OMP adj %d\n",error_rate_OMP_mod_adj(m_index),m(m_index));
        fprintf("Error Rate:%f for column size OMP topo %d\n",error_rate_OMP_topo(m_index),m(m_index));
    end
    % end

    %% PLOTS
    figure(2)
    %a1=semilogy(m,error_rate, 'm+-', 'LineWidth', 2);
    %a2=semilogy(m,error_rate_OMP, 'bd-', 'LineWidth', 2);
    %a3=semilogy(m,error_rate_de, 'g+-', 'LineWidth', 2);
    semilogy(m,error_rate_de, 'm+-', 'LineWidth', 2);
    hold on
    semilogy(m,error_rate_Ode, 'b+-', 'LineWidth', 2);
    
    %semilogy(m,error_rate_OMP_mod, 'g+-', 'LineWidth', 2);
    semilogy(m,error_rate_OMP_mod_adj, 'g+-', 'LineWidth', 2);
    semilogy(m,error_rate_OMP_topo, 'r+-', 'LineWidth', 2);
    %axis([0 120 0 1]);
    %hold off
    grid on
    %M1="n=5,h=4"; M2="n=9,h=4";
    %legend('SE CVX','SE OMP','DE CVX','DE OMP');
    %legend('DE OMP','DE OMP modified','DE modified OMP adjacency')
    legend('DE CVX','DE OMP','DE OMP modified','DE modified OMP+topology knowledge')
    %legend('h=3','h=4','h=5')
    title('Error rate vs provenance size DE OMP n=5');
    xlabel('Column size, m');
    ylabel('Reconstruction Error rate');
    hold off
%end

infos ={'Results for DEE: (I) CVX (II)1. OMP 2. Modified OMP with path constraints (complete graph) in form of:a) adjacency matrix: DE to SE b)brute force 3. (2.a) with topology knowledge';
    'n= no. of nodes';'h = path hop length';'DE= no. of possible double edges';'m = column size of Ar for OMP';'graph_adj_mat = topology adjacency matrix';'Error rate for different cases'};
filename = 'DEE_omp.mat';
save DEE_omp.mat infos n DE h graph_adj_mat path m error_rate_Ode error_rate_OMP_mod_adj error_rate_OMP_mod error_rate_OMP_topo