function [BER,SER] = MLD(n_symbols , n_iterations , bits_all , n_all , H_all , par)

n_errors=zeros(length(par.SNRdB_list),1); %Numero de Erros para cada SNR -> Simbolos
n_errors_bits=zeros(length(par.SNRdB_list),1); %Numero de Erros para cada SNR -> Bits

for loop1 = 1:n_iterations   
    H = H_all(:,:,loop1); %Matriz Canal com Distribuição Normal com media zero e varianza unitaria
    
    for ind_db = 1:length(par.SNRdB_list)
        sigma = sqrt(par.MT*par.Es*10^(-par.SNRdB_list(ind_db)/10));
        %N0 = par.MT*par.Es*10^(-par.SNRdB_list(ind_db)/10);         
        
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
            
            % Começo Detecção MLD
            x_ml = par.symbols(1)*ones(par.MT,1);            
            min_norm = norm(H*x_ml-y);
            x_decod = x_ml;
            
            for k = 1:(par.M)^par.MT %Numero de Possiveis Combinações
                for kk = par.MT:-1:1
                    if (mod(k-1,(par.M)^(par.MT-kk))==0 && k~=1)
                        ind = find(par.symbols == x_ml(kk)); %Indice do Elemento do Alfabeto
                        if ind==par.M
                            ind = 0;
                        end
                        x_ml(kk) = par.symbols(ind+1);
                    end
                end
                
                % Solução MLD-> Distancia Euclideana
                norm_act = norm(H*x_ml-y);
                if norm_act<min_norm
                    x_decod = x_ml;
                    min_norm = norm_act;
                end
            end
            [none,idxhat] = min(abs(x_decod*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
            bits_hat = par.bits(idxhat,:); %Simbolos a Bits
            %Conteo Erros
            err = err + sum(x~=x_decod); %Erro De Simbolos
            err_bits = err_bits + sum(sum(bits_all(:,:,loop2)~=bits_hat)); %Erro De Bits
        end %fim Loop2-> n_symbols
        n_errors(ind_db)=err;
        n_errors_bits(ind_db)=err_bits;
    end
end
            
SER = n_errors/(par.MT*n_iterations*n_symbols); %SER         
BER = n_errors_bits/(par.MT*par.Q*n_iterations*n_symbols); %BER         