function [BER , SER] = EP1(n_symbols,n_iterations,L_EP,bits_all,n_all,H_all,par,beta)

% Function: EP
% number of EP iterations (L_EP>0).  

N = 2*par.MT;
L = 2^(log2(par.M)/2); % Converting M-QAM in a PAM
% Alphabet
A=-(L-1):2:(L-1);
energy= (L*L-1)/3; % PAM energy
n_errors = zeros(length(par.SNRdB_list),1); %Numero de Erros para cada SNR -> Simbolos
n_errors_bits = zeros(length(par.SNRdB_list),1); %Numero de Erros para cada SNR -> Bits

for loop1=1:n_iterations    
    Hc = H_all(:,:,loop1); % normal random matrix with zero mean and unit variance
    H = [real(Hc) -imag(Hc) ; imag(Hc) real(Hc)]; %Transformei Matriz complexa a Matriz real
    
    for ind_db=1:length(par.SNRdB_list)
        sigma=sqrt(N*energy*10^(-par.SNRdB_list(ind_db)/10));
        %Inicialização Erro
        err = 0;
        err_bits = 0;        
        
        %parfor loop2=1:n_symbols
        for loop2=1:n_symbols            
            % generate transmit symbol
            idx = bi2de(bits_all(:,:,loop2),'left-msb')+1;
            x = par.symbols(idx).';             
            % apply noisy linear
            yc = Hc*x+sigma*n_all(:,loop2,loop1);
            y = [real(yc);imag(yc)];            
            
            % Initialization
            gamma = zeros(N,1);
            GAMMA = ones(N,1)/energy;
            cavity_i = zeros(N,L);
            
            % Constants
            %beta = 0.2;
            epsilon = 5e-7;
            
            % meand and variance vectors of all marginal qi(x) at iteration
            % l=0
            C_q = inv(sigma^(-2)*H'*H+diag(GAMMA));
            var_qi = diag(C_q);
            nu_qi = C_q*(sigma^(-2)*H'*y+gamma);          
            
            for l=1:L_EP                
                % Mean (ti) and variance (hi_2) of each cavity marginal
                hi_2 = var_qi./(1-var_qi.*GAMMA);
                ti = hi_2.*(nu_qi./var_qi-gamma);                    
                % Compute cavity marginal for each xi
                % Note: Both 'for' do the same thing, but the first one by
                % columns and the other by rows. Depend on the number of
                % columns (N) or rows (L), we will interest in one or
                % another for. 
                if L<=N
                    for k=1:L
                        cavity_i(:,k)=(1./sqrt(2*pi*hi_2)).*exp(-0.5*(A(k)-ti).^2./hi_2);
                    end
                else
                    for k=1:N
                        cavity_i(k,:)=(1/sqrt(2*pi*hi_2(k)))*exp(-0.5*(A-ti(k)).^2/hi_2(k));
                    end
                end
                
                % Normalize distribution (33)                
                K = sum(cavity_i,2);
                cavity_i = diag(1./K)*cavity_i;
                
                % Mean (nu_pi) and variance (var_pi) of distribution (33)
                nu_pi = sum((diag(A)*cavity_i')',2);
                var_pi = sum(((diag(A)*ones(L,N))'-diag(nu_pi)*ones(N,L)).^2.*cavity_i,2);
                var_pi = max([var_pi,epsilon*ones(N,1)],[],2);
                
                % Compute new values for gamma and GAMMA                
                GAMMA_aux = beta*((1./var_pi)-(1./hi_2))+(1-beta)*GAMMA;
                gamma_aux = beta*((nu_pi./var_pi)-(ti./hi_2))+(1-beta)*gamma;
                
                for k=1:N
                    if GAMMA_aux(k)>=0
                        GAMMA(k)=GAMMA_aux(k);
                        gamma(k)=gamma_aux(k);
                    end
                end
                
                % meand and variance vectors of all marginal qi(x) at next
                % iteration l
                C_q = inv(sigma^(-2)*H'*H+diag(GAMMA));
                var_qi = diag(C_q);
                nu_qi = C_q*(sigma^(-2)*H'*y+gamma);
            end
            
            % Decoder EP
            belief_x = zeros(N,L);
            for i=1:L
                belief_x(:,i) = abs(nu_qi-A(i));
            end
            
            [min_belief,ind] = min(belief_x,[],2);
            x_decod = A(ind)';
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