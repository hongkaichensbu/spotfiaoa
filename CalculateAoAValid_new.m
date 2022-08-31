% generate plot for each csi measurements in a folder for different number
% of paths


clear;

addpath('matlab');
addpath('spotfi');

% choose the folder
%dataDir = '../../data/19-6-11/20/';
%dataDir = '../../data/19-6-17/csi/';
%dataDir = 'data/2020-9-3/';
%dataDir = 'data/2020-9-5/';
%dataDir = 'data/2020-9-10/';
%dataDir = 'data/2020-10-19/';
dataDir = '';

num_paths = 4;%8;

%To change to offset 00, 11, 12, 21, 22, use '-off-00', '-off-11','-off-12', '-off-21', '-off-22' 
% off_cfg = '-off-00' 
% off_cfg = '-off-11' 
%off_cfg = '-off-12'
%off_cfg = '-off-21' 
%off_cfg = '-off-22' 
% change format of the csi matrix
% off_cfg = '-format-2-off-00' 
% off_cfg = '-format-2-off-11' 
% off_cfg = '-format-2-off-12'
% off_cfg = '-format-2-off-21' 
% off_cfg = '-format-2-off-22' 
% load('calibration/offset_12.mat');
% load('calibration/offset_13.mat');

off_cfg = '-417off-11' 
% off_cfg = '-format-2-417off-11'
load('calibration/2021-4-17-offset_12.mat');
load('calibration/2021-4-17-offset_13.mat');

% off_cfg = ""

%%%%%%%%%%%
%  setup for spotifi
%%%%%%%%%%
fc = 5.75e9; % center frequency
% fc = 5.2e9;
M = 3;    % number of rx antennas
fs = 40e6; % channel bandwidth
% fs = 20e6;
c = 3e8;  % speed of light
d = 2.6e-2;  % distance between adjacent antennas in the linear antenna array
% dTx = 2.6e-2;
SubCarrInd = [-58,-54,-50,-46,-42,-38,-34,-30,-26,-22,-18,-14,-10,-6,-2,2,6,10,14,18,22,26,30,34,38,42,46,50,54,58]; % WiFi subcarrier indices at which CSI is available
N = length(SubCarrInd); % number of subcarriers
% subCarrSize = 128;  % total number fo
fgap = 312.5e3; % frequency gap in Hz between successive subcarriers in WiFi
lambda = c/fc;  % wavelength
T = 1; % number of transmitter antennas

% MUSIC algorithm requires estimating MUSIC spectrum in a grid. paramRange captures parameters for this grid
% For the following example, MUSIC spectrum is caluclated for 101 ToF (Time of flight) values spaced equally between -25 ns and 25 ns. MUSIC spectrum is calculated for for 101 AoA (Angle of Arrival) values between -90 and 90 degrees.
paramRange = struct;
paramRange.GridPts = [101 181 1]; % number of grid points in the format [number of grid points for ToF (Time of flight), number of grid points for angle of arrival (AoA), 1]
paramRange.delayRange = [-100 100]*1e-9; % lowest and highest values to consider for ToF grid.
paramRange.angleRange = 90*[-1 1]; % lowest and values to consider for AoA grid.
do_second_iter = 0;
% paramRange.seconditerGridPts = [1 51 21 21];
paramRange.K = floor(M/2)+1; % parameter related to smoothing.
paramRange.L = floor(N/2); % parameter related to smoothing.
paramRange.T = 1;
paramRange.deltaRange = [0 0];

maxRapIters = Inf;
useNoise = 0;
paramRange.generateAtot = 2;
%paramRange.circularTx = 1;

%nPaths=4;
nSnaps=800;
totalAOA=[];%zeros(nPaths*nSnaps,1);
totalTOF=[];%zeros(nPaths*nSnaps,1);
index=1;
S_i=1;


% get all the files in the folder
files = dir(fullfile(dataDir,'*.mat'));
files = files(1);
% files = files(2);

%count = 0;
%poolobj = parpool(4);
%addAttachedFiles(poolobj,{'spotfi/GetQnBackscatter.m'});

% parpool(num_paths)
length(files)

colorstring = 'kbgrymc';

% go through each file
for k = 1:length(files)
    %fprintf('file count: %d\n', count); 
    filename = files(k).name
    
    filedata = strcat(dataDir,filename);
    
    load(filedata);
    
    data_size = length(csi_real(:,1));
    
    fill_v = zeros(1,5);
    
    csi_all = complex(csi_real,csi_imag);
    
    % each file generate plot for path #1-5
%     parfor i = 1:num_paths
    for i = num_paths

        nPaths=i;
        
        %totalAOA=[];%zeros(nPaths*nSnaps,1);
        %totalTOF=[];%zeros(nPaths*nSnaps,1);
        
        totalAOA = zeros(data_size,nPaths);
        totalTOF = zeros(data_size,nPaths);
        
        index=1;
 %       fig1=figure(1);
        
        %aoaSave = zeros(nPaths*data_size,1);
        figname = strrep(filename,'.mat','');
        
        %for ii=S_i:50:S_i+nSnaps
        for ii=1:data_size %1:20
            fprintf("File: %s nPath: %d #: %d\n",filename,nPaths,ii);
            
            
            csi = reshape(csi_all(ii,:),[3,30]);
%             csi = reshape(csi_all(ii,:), 30, [])';
            
            %plot(db(abs(squeeze(csi(1,:,:)).')))
            %db(get_eff_SNRs(csi), 'pow')
            
            %legend('RX Antenna A', 'RX Antenna B', 'RX Antenna C', 'Location', 'SouthEast' );
            %xlabel('Subcarrier index');
            %ylabel('SNR [dB]');
            %drawnow
            %pause(0.1)
            
            %csi_new = phaseOffset(squeeze(csi),phase_diff_12,phase_diff_13);
			
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% Change this for different offset combinations
             %To change to offset 11, 12, 21, 22, use '-off-11','-off-12', '-off-21', '-off-22' 
            csi_new = phaseOffset(squeeze(csi),offset_12(1,:),offset_13(1,:));
%             csi_new = phaseOffset(squeeze(csi),offset_12(1,:),offset_13(2,:));
%             csi_new = phaseOffset(squeeze(csi),offset_12(2,:),offset_13(1,:));
%             csi_new = phaseOffset(squeeze(csi),offset_12(2,:),offset_13(2,:));
            
            % Use this for no offset '-off-00'
%             csi_new = squeeze(csi);  
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			
            %sample_csi_trace=[reshape(csi_new(1,1,:),[30 1]);reshape(csi_new(1,2,:),[30 1]);reshape(csi_new(1,3,:),[30 1])];
            sample_csi_trace=[reshape(csi_new(1,:),[30,1]);reshape(csi_new(2,:),[30,1]);reshape(csi_new(3,:),[30,1])];
            %sample_csi_trace=[reshape(csi_new(1,:),[30,1]);reshape(csi_new(2,:),[30,1])];
            
            % sample CSI trace is a 90x1 vector where first 30 elements correspond to subcarriers for first rx antenna, second 30 correspond to CSI from next 30 subcarriers and so on.
            % replace sample_csi_trace with CSI from Intel 5300 converted to 90x1 vector
            %sample_csi_traceTmp = load('sample_csi_trace');
            %sample_csi_traceTmp = load('lcor_T1_e');
            %sample_csi_traceTmp = load('lab_e');
            %sample_csi_traceTmp = load('lab_s_r');
            %sample_csi_trace = sample_csi_traceTmp.sample_csi_trace;
            %sample_csi_trace = sample_csi_traceTmp.csi_traced;
            %sample_csi_trace = sample_csi_traceTmp.y;
            
            
            
            % ToF sanitization code (Algorithm 1 in SpotFi paper)
            csi_plot = reshape(sample_csi_trace, N, M);
            [PhsSlope, PhsCons] = removePhsSlope(csi_plot,M,SubCarrInd,N);
            ToMult = exp(1i* (-PhsSlope*repmat(SubCarrInd(:),1,M) - PhsCons*ones(N,M) ));
            csi_plot = csi_plot.*ToMult;
            relChannel_noSlope = reshape(csi_plot, N, M, T);
            sample_csi_trace_sanitized = relChannel_noSlope(:);
            
            % MUSIC algorithm for estimating angle of arrival
            % aoaEstimateMatrix is (nComps x 5) matrix where nComps is the number of paths in the environment. First column is ToF in ns and second column is AoA in degrees as defined in SpotFi paper
            aoaEstimateMatrix = backscatterEstimationMusic(sample_csi_trace_sanitized, M, N, c, fc,...
                T, fgap, SubCarrInd, d, paramRange, maxRapIters, useNoise, do_second_iter, ones(nPaths));
            tofEstimate = aoaEstimateMatrix(:,1); % ToF in nanoseconds
            aoaEstimate = aoaEstimateMatrix(:,2); % AoA in degrees
            
            %totalAOA((index-1)*nPaths+1:index*nPaths) = aoaEstomate;
            %totalTOF((index-1)*nPaths+1:index*nPaths)= tofEstimate;
            %index = index +1;
            
            
            totalAOA(ii, 1:length(aoaEstimate) ) = aoaEstimate.';
            totalTOF(ii, 1:length(tofEstimate) ) = tofEstimate.';
            
            %%% Save the data every 1000 csi packets. incase power pc
            %%% shutdown..
            if ii==1000
                matname_aoa = strcat(dataDir,'aoadata/',figname,'-aoa-nPaths-',num2str(i),off_cfg,'-40mhz.mat');
                parsave(matname_aoa, totalAOA);
        
                matname_tof = strcat(dataDir,'aoadata/',figname,'-tof-nPaths-',num2str(i),off_cfg,'-40mhz.mat');
                parsave(matname_tof, totalTOF);
            
            elseif mod(ii,1000) == 0
                matname_aoa = strcat(dataDir,'aoadata/',figname,'-aoa-nPaths-',num2str(i),off_cfg,'-40mhz.mat');
                parsaveapp(matname_aoa, totalAOA);
        
                matname_tof = strcat(dataDir,'aoadata/',figname,'-tof-nPaths-',num2str(i),off_cfg,'-40mhz.mat');
                parsaveapp(matname_tof, totalTOF);
            end
            
            %n_path = 5;
            %totalAOA = [totalAOA, aoaEstomate(1:n_path)];
            %totalTOF = [totalTOF, tofEstimate(1:n_path)];
            
            %for j = 1:length(aoaEstimate)
            %    plot(tofEstimate(j),aoaEstimate(j),'*','Color', colorstring(j))
            %    hold on
            %end
            %axis([-110 110 -100 100])
            
            %if length(aoaEstimate) == 1
            %    legend('#1')
            %elseif length(aoaEstimate) == 2
            %    legend('#1','#2')
            %elseif length(aoaEstimate) == 3
            %    legend('#1','#2','#3')
            %elseif length(aoaEstimate) == 4
            %    legend('#1','#2','#3','#4')
            %elseif length(aoaEstimate) == 5
            %    legend('#1','#2','#3','#4','#5')
            %end
            %hold on
            %figname = strrep(filename,'.dat','');
            %%figname = strrep(figname,'_','-'); %%%%%%%%
%             title_name = strcat(figname, '-nPaths-',num2str(i),off_cfg)
%             title(title_name)
%             xlabel('ToF')
%             ylabel('AoA')
%             
%             
%             figname = strcat(dataDir,'plots/',figname,'-nPaths-',num2str(i),'-',num2str(ii),off_cfg,'.png');
%             saveas(fig1,figname);
            %clf(fig1);
            
%             if isempty(aoaSave)
%                 aoaSave = aoaEstimateMatrix;
%             else
%                 aoaSave = cat(1,aoaSave,aoaEstimateMatrix);
%                 
%                 if length(aoaEstimate) < i
%                     for j=1:(i-length(aoaEstimate))
%                         aoaSave = cat(1,aoaSave,fill_v)
%                     end
%                 end
%             end
            
        end
        
%        clf(fig1);
        
        % save fig as jpg
        figname = strrep(filename,'.mat','');
        %%figname = strrep(figname,'_','-'); %%%%%%%%%%%%%%%%%
        %matname = strcat(dataDir,'aoadata/',figname,'-nPaths-',num2str(i),off_cfg,'.mat');
        %figname = strcat('bosch/plots/',figname,'-nPaths-',num2str(i),'.png');
        %saveas(fig1,figname);
        %parsave(matname, aoaSave);
        %clear figure
        %clf(fig1);
        
        matname_aoa = strcat(dataDir,'aoadata/',figname,'-aoa-nPaths-',num2str(i),off_cfg,'-40mhz.mat');
        if data_size < 1000
            parsave(matname_aoa, totalAOA);
        else
            parsaveapp(matname_aoa, totalAOA);
        end
        
        matname_tof = strcat(dataDir,'aoadata/',figname,'-tof-nPaths-',num2str(i),off_cfg,'-40mhz.mat');
        if data_size < 1000
            parsave(matname_tof, totalTOF);
        else
            parsaveapp(matname_tof, totalTOF);
        end
        
        
    end
    
    %count = count + 1;
end