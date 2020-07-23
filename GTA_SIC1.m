function [BER,SER] = GTA_SIC1(n_symbols,n_iterations,bits_all,n_all,H_all,par)

N = 2*par.MT;
L = 2^(log2(par.M)/2); % Converting M-QAM in a PAM
% Alphabet
A=-(L-1):2:(L-1);
energy= (L*L-1)/3; % PAM energy
n_errors = zeros(length(par.SNRdB_list),1); %Numero de Erros para cada SNR -> Simbolos
n_errors_bits = zeros(length(par.SNRdB_list),1); %Numero de Erros para cada SNR -> Bits

for loop1=1:n_iterations
    for ind_db=1:length(par.SNRdB_list)
        sigma=sqrt(N*energy*10^(-par.SNRdB_list(ind_db)/10));
        %Inicialização Erro
        err = 0;
        err_bits = 0;

        
        for loop2=1:n_symbols            
            Hc = H_all(:,:,loop1); % normal random matrix with zero mean and unit variance
            H = [real(Hc) -imag(Hc) ; imag(Hc) real(Hc)]; %Transformei Matriz complexa a Matriz real
            % generate transmit symbol
            idx = bi2de(bits_all(:,:,loop2),'left-msb')+1;
            x = par.symbols(idx).';             
            % apply noisy linear
            yc = Hc*x+sigma*n_all(:,loop2,loop1);
            y = [real(yc);imag(yc)];
            
            % SIC (Successive Interference Cancelation)
            x_decod=zeros(N,1);
            original_ind=(1:N)'; % Keep original order of variables 
            
            for k_sic=1:N                    
                % ////////////////////////
                % /////// Tree BP ////////
                % ////////////////////////                
                % Compute weights of total graph 
                C=sigma^2*pinv(H'*H+(sigma^2/energy)*eye(N-k_sic+1)); % Covariance matrix                 
                % MMSE estimator
                %z=inv(H'*H+(sigma^2/energy)*eye(N-k_sic+1))*H'*y;
                z=sigma^(-2)*C*H'*y; % This way avoids to do an inverse again
                
                if k_sic==N % Last iteration is special because there is only one node (tree is not necessary)
                    xi=A;
                    belief_x=f_xi(xi,z(1),C(1,1));
                    % Decoder
                    [max_belief,ind]=max(belief_x,[],2);
                    x_decod(original_ind(1))=A(ind);
                    break
                end
                
                % We look for the smallest diagonal element of C. Root node
                % will be the index of this element
                [minC,root]=min(diag(C));                
                w=C.^2./(diag(C)*diag(C)');
                
                % As we look for maximum weights, rewrite diagonals of w as -1
                % (we are not interested in weights of nodes with theirselves)
                
                w_aux=w-2*diag(diag(w));
                
                % Compute maximum spanning                
                parent=zeros(N-k_sic+1,1); % Keep father of each node, 
                % i.e. parent(i) is the parent node of i   
                
                % Assume that root node is node root
                parent(root)=0;
                inTree=root; % Keep nodes added to the tree 
                Tree=zeros(N-k_sic+1); % Edges of tree are keep as '1'
                
                for k=1:N-k_sic
                    % We look for maximal weights of nodes inside the tree
                    [max_ws,inds]=max(w_aux(inTree,:),[],2);
                    % We only save the maximum weigth of previous weights
                    [max_w,ind]=max(max_ws);
                    % Write as -1 the edges weigth inside the tree (to not consider
                    % them in future iterations) 
                    new_node=inds(ind);
                    w_aux(inTree,new_node)=-1;
                    w_aux(new_node,inTree)=-1;
                    % Add new node to the tree 
                    parent(new_node)=inTree(ind);
                    inTree=[inTree new_node];
                    Tree(new_node,parent(new_node))=1;
                    Tree(parent(new_node),new_node)=1;
                end
                
                
                % f(xi; z, C)
                f_xi=@(xi,zi,Cii)(exp(-0.5*((xi-zi).^2)/Cii));
                
                % f(xi|xj; z, C) (where xj=x_p(i), i.e, the parent node of
                % i)
                f_xi_xj=@(xi,zi,xj,zj,Cij,Cii,Cjj)(exp(-0.5*(((xi-zi)-(Cij/Cjj)*(xj-zj)).^2)/(Cii-(Cij^2/Cjj))));
                
                % Look for the tree leaves and its edges 
                n_edges=sum(Tree,1)'; % Number of edges of each node 
                n_children=n_edges-1; % Number of chilgren of each node 
                n_children(root)=n_children(root)+1;
                leaf_node=n_children==0; % Leaves are kept as '1' 
                
                % ////////////////////
                % /////// BP /////////
                % ////////////////////
                
                % //////////////////////////////////////////////
                % ///// DOWNWARD BP MESSAGES (m_i-->p(i)) //////
                % //////////////////////////////////////////////                
                messages_down=ones(N-k_sic+1,L); % Each column is associated with an alphabet symbol
                % Initialized as '1' to merge (25) and (26). This means we
                % consider messages from leaves to children as '1'                
                i_nodes=find(leaf_node==1); % We start from leaves                 
                n_messages_to_father=zeros(N-k_sic+1,1); % Keep the number os messages sended to each parent node
                
                while (isempty(i_nodes)==0)
                    i_nodes_sig=[];
                    for k=1:length(i_nodes)
                        i=i_nodes(k);
                        pi=parent(i);
                        xi=A;
                        children_i=find(parent==i);
                        messages_children=prod(messages_down(children_i,:),1); % multiply messages from children to parent node i
                        
                        for kk=1:L
                            xpi=A(kk);
                            messages_down(i,kk)=sum(f_xi_xj(xi,z(i),xpi,z(pi),C(i,pi),C(i,i),C(pi,pi)).*messages_children);
                        end
                        
                        n_messages_to_father(pi)=n_messages_to_father(pi)+1;
                        
                        if (n_edges(pi)-1==n_messages_to_father(pi) && pi~=root) % If parent pi is the root node, it can not added as next node because downward messages end in this node
                            i_nodes_sig=[i_nodes_sig; pi]; % pi is added as next node because it has already received all messages from its children and it can send a message to its parent 
                        end
                    end
                    i_nodes=i_nodes_sig;
                end
                % ////////////////////////////
                % ///// BELIEF VARIABLE //////
                % ////////////////////////////                
                % After message passing, we compute the belief root
                % variable                 
                xi=A;
                children_i=find(parent==root);
                belief_x=f_xi(xi,z(root),C(root,root)).*prod(messages_down(children_i,:),1);                
                % ////////////////////
                % ///// DECODER //////
                % ////////////////////                
                % Decoding                
                [max_belief,ind]=max(belief_x,[],2);
                x_decod(original_ind(root))=A(ind);                
                % Interference Cancelation
                y=y-H(:,root)*x_decod(original_ind(root));
                % We copy last column of H to root column and delete the
                % last column
                H(:,root)=H(:,end);
                H=H(:,1:N-k_sic);                
                original_ind(root)=original_ind(N-k_sic+1); 
            end            
            x_decod_c = x_decod(1:par.MT)+1i*x_decod(par.MT+1:2*par.MT);
            [none,idxhat] = min(abs(x_decod_c*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
            x_hat = par.symbols(idxhat).'; %Decodificação do x_zf - Quantização
            bits_hat = par.bits(idxhat,:); %Simbolos a Bits
            %Conteo Erros
            err = err + sum(x~=x_hat); %Erro De Simbolos
            err_bits = err_bits + sum(sum(bits_all(:,:,loop2)~=bits_hat)); %Erro De Bits
        end
        n_errors(ind_db)=err;
        n_errors_bits(ind_db)=err_bits;
    end
end
SER = n_errors/(par.MT*n_iterations*n_symbols); %SER         
BER = n_errors_bits/(par.MT*par.Q*n_iterations*n_symbols); %BER         