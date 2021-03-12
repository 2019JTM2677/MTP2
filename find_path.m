function [path_selected,source]= find_path(adj_mat,nodes,hop_len) 
    
    n = length(nodes);
    dest = nodes(n);
    dest.Node_id
    flag=0;
    for j=1:n
        src = nodes(j);
        %src.Node_id
        [shortestPaths,totalCosts] =  kShortestPath(adj_mat, src.Node_id,dest.Node_id, 12); 
        for i = 1: length(shortestPaths)
            fprintf('Path # %d:\n',i);
            disp(shortestPaths{i})
            fprintf('Cost of path %d is %5.2f\n\n',i,totalCosts(i));
            if totalCosts(i)== hop_len
                path_selected = shortestPaths{i};
                flag=1;
                break
            end
        end
        if flag==1
            source = src;
            break
        end
        
    end
    fprintf("Source node: %d",src.Node_id);
end
