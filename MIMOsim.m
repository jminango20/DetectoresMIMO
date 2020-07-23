function MIMOsim(varargin)
  % -- Parametros por default/Parametros do Usuario
  if isempty(varargin)
    disp('Usando Parametros Simulação Por Default...')        
    % Parametros por default
    par.simName = 'ERR_4x4_16QAM'; % Nome da Simulação (usado para salvar resultados)
    par.runId = 0; % simulação ID (usado para reproduzir resultados)
    par.MR = 32; % Antenas Receptoras
    par.MT = 32; % Antenas Transmissoras (não maiores que MR!) 
    par.mod = '16QAM'; % Tipo Modulação: 'BPSK','QPSK','16QAM','64QAM'
    par.SNRdB_list = 10:3:25; % Valores da SNR [dB] a ser simulados
    par.SNRdB_fixo = 5; % Valores da SNR [dB] a ser simulado Fixo
    par.beta = 0:0.1:1; %Parametro smooth beta.    
    %par.detector = {'MLD','ZF','MMSE','ZF_SIC1','MMSE_SIC','SD_Studer'}; % Definição detectore(s) a serem simulados
    %par.detector = {'ZF','GTA','GTA_SIC','EP','SD_Studer'}; % Definição detectore(s) a serem simulados
    %par.detector = {'EP','EP1','ZF'}; % Definição detectore(s) a serem simulados
    par.detector = {'EP-Beta'}; % Definição detectore(s) a serem simulados    

  else      
    disp('usando parametros já estabelecidos...')    
    par = varargin{1}; % só argumento de estrutura
  end

  % -- Inizialização
  % use runId random seed (enables reproducibility)
  %rng(par.runId); 

  % Alfabeto da Constelação com Mapeamento Gray-mapped 
  switch (par.mod)
    case 'BPSK',
      par.symbols = [ -1 1 ];
    case 'QPSK', 
      par.symbols = [ -1-1i,-1+1i, ...
                      +1-1i,+1+1i ];
    case '16QAM',
      par.symbols = [ -3-3i,-3-1i,-3+3i,-3+1i, ...
                      -1-3i,-1-1i,-1+3i,-1+1i, ...
                      +3-3i,+3-1i,+3+3i,+3+1i, ...
                      +1-3i,+1-1i,+1+3i,+1+1i ];
    case '64QAM',
      par.symbols = [ -7-7i,-7-5i,-7-1i,-7-3i,-7+7i,-7+5i,-7+1i,-7+3i, ...
                      -5-7i,-5-5i,-5-1i,-5-3i,-5+7i,-5+5i,-5+1i,-5+3i, ...
                      -1-7i,-1-5i,-1-1i,-1-3i,-1+7i,-1+5i,-1+1i,-1+3i, ...
                      -3-7i,-3-5i,-3-1i,-3-3i,-3+7i,-3+5i,-3+1i,-3+3i, ...
                      +7-7i,+7-5i,+7-1i,+7-3i,+7+7i,+7+5i,+7+1i,+7+3i, ...
                      +5-7i,+5-5i,+5-1i,+5-3i,+5+7i,+5+5i,+5+1i,+5+3i, ...
                      +1-7i,+1-5i,+1-1i,+1-3i,+1+7i,+1+5i,+1+1i,+1+3i, ...
                      +3-7i,+3-5i,+3-1i,+3-3i,+3+7i,+3+5i,+3+1i,+3+3i ];                         
  end

  % Obtenção average de Energia de Simbolo
  par.Es = mean(abs(par.symbols).^2);   
  % precomputo bits
  par.Q = log2(length(par.symbols)); % numero de bits por simbolo
  par.M = length(par.symbols); % number of simbolos da modulação
  par.bits = de2bi(0:length(par.symbols)-1,par.Q,'left-msb');
  par.n_iterations = 5;     % numero de diferentes canais de H
  par.n_symbols = 500;    % numero de vetores Transmitidos por cada canal de H


  %% Parametros da Simulação
  n = sqrt(0.5)*(randn(par.MR,par.n_symbols,par.n_iterations)+1i*randn(par.MR,par.n_symbols,par.n_iterations)); %Vetor Ruidos
  H = sqrt(0.5)*(randn(par.MR,par.MT,par.n_iterations)+1i*randn(par.MR,par.MT,par.n_iterations)); %Vetor Canal Rayleigh
  bits = randi([0 1],par.MT,par.Q,par.n_symbols,par.n_iterations);  %Bits Gerados para a Transmissão

  %% Parametros Dos Resultados
  if length(par.detector)~=1
      res.SER = zeros(length(par.SNRdB_list),length(par.detector)); % symbol error rate
      res.BER = zeros(length(par.SNRdB_list),length(par.detector)); % bit error rate
      res.time = zeros(length(par.detector),1); % Tempo Ejecução de Cada Algoritmo
  else %para EP em 3D, em funcion de diferentes betas
      res.SER = zeros(length(par.SNRdB_list),length(par.beta)); % symbol error rate
      res.BER = zeros(length(par.SNRdB_list),length(par.beta)); % bit error rate
      res.time = zeros(length(par.detector),length(par.beta)); % Tempo Ejecução de Cada Algoritmo      
  end
  
  
  
  %% -- COMEÇA A SIMULAÇÃO  
  %parfor i=1:n_iterations      
  for i=1:par.n_iterations
      i
      
      for d=1:length(par.detector)         
        switch (par.detector{d}) % Escolha Algoritmos de Detecção

%% DETECTOR MLD            
            case 'MLD', % MLD Detecção - Tradicional            
                tic;
                [BER , SER] = MLD(par.n_symbols , 1 , bits(:,:,:,i) , n(:,:,i) , H(:,:,i) , par);
                res.SER(:,d) = res.SER(:,d) + SER/par.n_iterations; 
                res.BER(:,d) = res.BER(:,d) + BER/par.n_iterations; 
                res.time(d,1) = res.time(d,1) + toc;
                
%% DETECTOR ESFERICO            
            case 'SD', % MLD Detecção - Esferico            
                tic;
                [BER , SER] = SD(par.n_symbols , 1 , bits(:,:,:,i) , n(:,:,i) , H(:,:,i) , par);
                res.SER(:,d) = res.SER(:,d) + SER/par.n_iterations; 
                res.BER(:,d) = res.BER(:,d) + BER/par.n_iterations; 
                res.time(d,1) = res.time(d,1) + toc;
                
            case 'SD_Studer', % MLD Detecção - Esferico Studer: Mais Rapido que do Acima            
                tic;
                [BER , SER] = SD_Studer(par.n_symbols , 1 , bits(:,:,:,i) , n(:,:,i) , H(:,:,i) , par);
                res.SER(:,d) = res.SER(:,d) + SER/par.n_iterations; 
                res.BER(:,d) = res.BER(:,d) + BER/par.n_iterations; 
                res.time(d,1) = res.time(d,1) + toc;

%% DETECTOR ZERO FORCING            
            case 'ZF', % ZF Detecção - Tradicional            
                tic;
                [BER , SER] = ZF(par.n_symbols , 1 , bits(:,:,:,i) , n(:,:,i) , H(:,:,i) , par);
                res.SER(:,d) = res.SER(:,d) + SER/par.n_iterations; 
                res.BER(:,d) = res.BER(:,d) + BER/par.n_iterations; 
                res.time(d,1) = res.time(d,1) + toc;                

            case 'ZF_SIC1', % ZF-SIC Detecção - Tradicional            
                tic;
                [BER , SER] = ZF_SIC1(par.n_symbols , 1 , bits(:,:,:,i) , n(:,:,i) , H(:,:,i) , par);
                res.SER(:,d) = res.SER(:,d) + SER/par.n_iterations; 
                res.BER(:,d) = res.BER(:,d) + BER/par.n_iterations; 
                res.time(d,1) = res.time(d,1) + toc;                

%% DETECTOR MMMSE                            
            case 'MMSE', % MMSE Detecção - Tradicional            
                tic;
                [BER , SER] = MMSE(par.n_symbols , 1 , bits(:,:,:,i) , n(:,:,i) , H(:,:,i) , par);
                res.SER(:,d) = res.SER(:,d) + SER/par.n_iterations; 
                res.BER(:,d) = res.BER(:,d) + BER/par.n_iterations; 
                res.time(d,1) = res.time(d,1) + toc;                
                
            case 'MMSE_SIC1', % MMSE-SIC Detecção - Tradicional            
                tic;
                [BER , SER] = MMSE_SIC1(par.n_symbols , 1 , bits(:,:,:,i) , n(:,:,i) , H(:,:,i) , par);
                res.SER(:,d) = res.SER(:,d) + SER/par.n_iterations; 
                res.BER(:,d) = res.BER(:,d) + BER/par.n_iterations; 
                res.time(d,1) = res.time(d,1) + toc;

%% DETECTOR GTA                            
            case 'GTA', % GTA Detecção - Tradicional            
                tic;
                [BER , SER] = GTA1(par.n_symbols , 1 , bits(:,:,:,i) , n(:,:,i) , H(:,:,i) , par);
                res.SER(:,d) = res.SER(:,d) + SER/par.n_iterations; 
                res.BER(:,d) = res.BER(:,d) + BER/par.n_iterations; 
                res.time(d,1) = res.time(d,1) + toc;          

%% DETECTOR GTA-SIC                            
            case 'GTA_SIC', % GTA Detecção - Tradicional            
                tic;
                [BER , SER] = GTA_SIC1(par.n_symbols , 1 , bits(:,:,:,i) , n(:,:,i) , H(:,:,i) , par);
                res.SER(:,d) = res.SER(:,d) + SER/par.n_iterations; 
                res.BER(:,d) = res.BER(:,d) + BER/par.n_iterations; 
                res.time(d,1) = res.time(d,1) + toc;          
                
%% DETECTOR EP-> Expected Propagation                            
            case 'EP', % EP-> Expected Propagation
                tic;
                L_EP = 10;
                beta = 0.2;
                [BER , SER] = EP1(par.n_symbols , 1 , L_EP, bits(:,:,:,i) , n(:,:,i) , H(:,:,i) , par, beta);
                res.SER(:,d) = res.SER(:,d) + SER/par.n_iterations; 
                res.BER(:,d) = res.BER(:,d) + BER/par.n_iterations; 
                res.time(d,1) = res.time(d,1) + toc;                          

%% DETECTOR EP-> Expected Propagation                            
            case 'EP1', % EPProva-> Vou Trabalhar sem o Parametro beta, ou seja, beta=1
                tic;
                L_EP = 10;
                beta = 1;
                [BER , SER] = EPProva(par.n_symbols , 1 , L_EP, bits(:,:,:,i) , n(:,:,i) , H(:,:,i) , par, beta);
                res.SER(:,d) = res.SER(:,d) + SER/par.n_iterations; 
                res.BER(:,d) = res.BER(:,d) + BER/par.n_iterations; 
                res.time(d,1) = res.time(d,1) + toc;                          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%EP com Diferentes Betas
%% DETECTOR EP-> Expected Propagation com Diferentes Betas                           
            case 'EP-Beta', % EPProva-> Vou Trabalhar sem o Parametro beta, ou seja, beta=1
                tic;
                L_EP = 10;
                [BER , SER] = EP_Beta(par.n_symbols , 1 , L_EP, bits(:,:,:,i) , n(:,:,i) , H(:,:,i) , par);
                res.SER = res.SER + SER/par.n_iterations; 
                res.BER = res.BER + BER/par.n_iterations; 
                res.time = res.time + toc;                          
                
                
                
            otherwise,            
                error('tipo de par.detector não definido.')      
        end
      end %fim d->detector
      
  end %fim i->n_iteration
  
  %% MOSTRA DOS RESULTADOS
  save([ par.simName '_' num2str(par.runId) ],'par','res');    
  
%% -- Mostra Resultados (Generação Matlab plot)  -> BER
  marker_style = {'bo-','rs--','mv-.','kp:','g*-','c>--','yx:'};
  figure(1)
  mesh(par.SNRdB_list,par.beta,res.BER');
  set(gca, 'ZScale', 'log')
  
  for d=1:length(par.detector)
      if d==1
          semilogy(par.SNRdB_list,res.BER(:,d),marker_style{d},'LineWidth',2)
          hold on      
      else          
          semilogy(par.SNRdB_list,res.BER(:,d),marker_style{d},'LineWidth',2)
      end    
  end
  %hold off
  grid on
  xlabel('SNR [dB]','FontSize',12)
  ylabel('bit error rate (BER)','FontSize',12)
  axis([min(par.SNRdB_list) max(par.SNRdB_list) 1e-5 1])
  legend(par.detector,'FontSize',12)
  title(['Comparação de Varios Detectores em ',num2str(par.MT),'x',num2str(par.MR),' sistema, ', par.mod, ' simbolos']);
  set(gca,'FontSize',12)

  %% -- Mostra Resultados (Generação Matlab plot)  -> SER
  figure(2)
  for d=1:length(par.detector)
      if d==1
          semilogy(par.SNRdB_list,res.SER(:,d),marker_style{d},'LineWidth',2)
          hold on      
      else          
          semilogy(par.SNRdB_list,res.SER(:,d),marker_style{d},'LineWidth',2)
      end    
  end
  %hold off
  grid on
  xlabel('SNR [dB]','FontSize',12)
  ylabel('Symbol Error Rate (SER)','FontSize',12)
  axis([min(par.SNRdB_list) max(par.SNRdB_list) 1e-5 1])
  legend(par.detector,'FontSize',12)
  title(['Comparação de Varios Detectores em ',num2str(par.MT),'x',num2str(par.MR),' sistema, ', par.mod, ' simbolos']);
  set(gca,'FontSize',12)
  
  %% -- Mostra Resultados (Generação Matlab plot)  -> Tempo Simulação
  figure(3)
  bar(res.time/60)
  set(gca,'XTickLabel',par.detector);
  xlabel('Detectores'), ylabel('Tempo Computacional (min)')

  
     
