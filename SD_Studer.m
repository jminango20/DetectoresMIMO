function [BER,SER] = SD_Studer(n_symbols , n_iterations , bits_all , n_all , H_all , par) %Este vai chamar à função sphdec

n_errors=zeros(length(par.SNRdB_list),1); %Numero de Erros para cada SNR -> Simbolos
n_errors_bits=zeros(length(par.SNRdB_list),1); %Numero de Erros para cada SNR -> Bits

for loop1 = 1:n_iterations   
    H = H_all(:,:,loop1); %Matriz Canal com Distribuição Normal com media zero e varianza unitaria
    
    for ind_db = 1:length(par.SNRdB_list)
        sigma = sqrt(par.MT*par.Es*10^(-par.SNRdB_list(ind_db)/10));        
        %Inicialização Erro
        err = 0;
        err_bits = 0;
        
        %parfor loop2=1:n_symbols
        for loop2 = 1:n_symbols
            % generate transmit symbol
            idx = bi2de(bits_all(:,:,loop2),'left-msb')+1;
            x = par.symbols(idx).'; 
            % Sinal Recebido = Canal*SinalEnviado + Ruido
            y = H*x + sigma*n_all(:,loop2,loop1);
            
            %% Começa Detector Esferico
            %Pre-Processamento com MMSE
            conteo = 0;
            G = pinv( H'*H + (sigma^2/par.Es)*eye(par.MT) )*H';
            x_mmse = G*y; %MMSE
            [none,idxhat] = min(abs(x_mmse*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
            x_mmse = par.symbols(idxhat).'; %Decodificação do x_mmse - Quantização            
            % -- initialization  
            %Radius = inf;
            Radius = norm(y-H*x_mmse)^2;  %Raio com Valor Inicial
            PA = zeros(par.MT,1); % path
            ST = zeros(par.MT,length(par.symbols)); % stack  
            % -- preprocessing
            [Q,R] = qr(H,0);  
            y_hat = Q'*y;    
            % -- add root node to stack
            Level = par.MT; 
            ST(Level,:) = abs(y_hat(Level)-R(Level,Level)*par.symbols.').^2;  
            % -- begin sphere decoder
            while ( Level<=par.MT )          
                % -- find smallest PED in boundary    
                [minPED,idx] = min( ST(Level,:) );
                % -- only proceed if list is not empty
                if minPED<inf
                    ST(Level,idx) = inf; % mark child as tested        
                    NewPath = [ idx ; PA(Level+1:end,1) ]; % new best path           
                    % -- search child
                    if ( minPED<Radius )
                        % -- valid candidate found
                        if ( Level>1 )                  
                            % -- expand this best node
                            PA(Level:end,1) = NewPath;
                            Level = Level-1; % downstep
                            DF = R(Level,Level+1:end) * par.symbols(PA(Level+1:end,1)).';
                            ST(Level,:) = minPED + abs(y_hat(Level)-R(Level,Level)*par.symbols.'-DF).^2;
                        else                            
                            % -- valid leaf found
                            idxML = NewPath;
                            bits_hat = par.bits(idxML',:);
                            % -- update radius (radius reduction)
                            Radius = minPED;
                            conteo = 1;
                        end                        
                    end                   
                else                    
                    % -- no more childs to be checked
                    Level=Level+1;      
                end                
            end %Fim Esferico    
            if conteo==1
                x_sd = par.symbols(idxML).'; %Decodificação do x_SD - Quantização
            else
                x_sd = x_mmse;
                [none,idxhat] = min(abs(x_sd*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
                bits_hat = par.bits(idxhat,:); %Simbolos a Bits
            end
            %Conteo Erros
            err = err + sum(x~=x_sd); %Erro De Simbolos
            err_bits = err_bits + sum(sum(bits_all(:,:,loop2)~=bits_hat)); %Erro De Bits
        end %fim Loop2-> n_symbols
        n_errors(ind_db)=err;
        n_errors_bits(ind_db)=err_bits;
    end
end
            
SER = n_errors/(par.MT*n_iterations*n_symbols); %SER         
BER = n_errors_bits/(par.MT*par.Q*n_iterations*n_symbols); %BER         