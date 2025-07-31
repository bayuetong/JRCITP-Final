
runTests=1;
sensor_test=true;
hearbeat_test=true;
speed_test=false;
real_time=true;
ts=0.06;
img_filename = '~/Desktop/matlab/images/jarcat.jpg';
%function unified_usrp_test_30k(runTests, img_filename)

p = gcp('nocreate');
if (isempty(p)); parpool(24); end

addpath("~/Desktop/matlab/Sensor_code"); addpath("~/Desktop/matlab/test_data");
addpath('/home/xzhang/Documents/MATLAB/Examples/R2023b/5g/NewRadioPUSCHThroughputExample'); addpath('~/Desktop/matlab/images');
addAttachedFiles(gcp,["read_complex_binary.m" "ProcessBuffer.m"])

%% Tests Variables %%
% runTests = 1; % Boolean for running tests
strCOM = ['/dev/ttyACM1' ...
    ''];
% img_filename = '~/Desktop/matlab/images/jarcat.jpg';

hf=load ('hf_vs_RC_03_024_513_16_666HZ.dat');
%hf_hb=load('hf_vs_RC_14_15_513_16_666HZ.dat');

hf_hb=load('hf_vs_RC_13_07_257_16_666HZ.dat');

hf_sensor=load('hf_vs_RC_03_024_257_10HZ.dat');

hf_25HZ=load('hf_vs_BR_RC_03_024_257_25HZ.dat');
hf_33HZ=load('hf_vs_BR_RC_03_024_257_33_333HZ.dat');
hf_hb_25HZ=load('hf_vs_HB_RC_13_07_257_25HZ.dat');
hf_hb_33HZ=load('hf_vs_HB_RC_13_07_257_33_333HZ.dat');

if ts==0.03
    hf=hf_33HZ;
    hf_hb=hf_hb_33HZ;

elseif ts==0.04
            hf=hf_25HZ;
    hf_hb=hf_hb_25HZ;
end

%% USRP signal analysis variables initialisation %%
load ('U1WaveQPSK100M_upsample_200Msps')
load ('waveforminfo_30k.mat')
load ('5G_30k_signal.mat')
% hf=load ('hf_vs_RC_03_024_513_16_666HZ.dat');
% hf_hb=load('hf_vs_RC_14_15_513_16_666HZ.dat');
% hf_sensor=load('hf_vs_RC_03_024_257_10HZ.dat');

N_antenna = 1;

%% Parameters for simulation

INUM=1000;  %% number of Monte-Carlo tests
%----- boilerplate image transmission params
cj=sqrt(-1);
doa_enabled=false;
use_decoded_source_signal=true;
channel_interpolation_used=false;
self_interfernce_cancellation_enabled=false;

simParameters = struct();       % Clear simParameters variable to contain all key simulation parameters
simParameters.NFrames = 1;      % Number of 10 ms frames
simParameters.SNRIn = [-5 0 5]; % SNR range (dB)
simParameters.PerfectChannelEstimator = false;
% simParameters.DisplaySimulationInformation = true;
% simParameters.DisplayDiagnostics = false;

% Set waveform type and PDSCH numerology (SCS and CP type)
simParameters.Carrier = nrCarrierConfig;         % Carrier resource grid configuration
simParameters.Carrier.NSizeGrid = 270;            % Bandwidth in number of resource blocks (51 RBs at 30 kHz SCS for 20 MHz BW)
simParameters.Carrier.SubcarrierSpacing = 30;    % 15, 30, 60, 120 (kHz)
simParameters.Carrier.CyclicPrefix = 'Normal';   % 'Normal' or 'Extended' (Extended CP is relevant for 60 kHz SCS only)
simParameters.Carrier.NCellID = 1;               % Cell identity

% PDSCH/DL-SCH parameters
simParameters.PDSCH = nrPDSCHConfig;      % This PDSCH definition is the basis for all PDSCH transmissions in the BLER simulation
simParameters.PDSCHExtension = struct();  % This structure is to hold additional simulation parameters for the DL-SCH and PDSCH

% Define PDSCH time-frequency resource allocation per slot to be full grid (single full grid BWP)
simParameters.PDSCH.PRBSet = 0:simParameters.Carrier.NSizeGrid-1;                 % PDSCH PRB allocation
simParameters.PDSCH.SymbolAllocation = [0,simParameters.Carrier.SymbolsPerSlot];  % Starting symbol and number of symbols of each PDSCH allocation
simParameters.PDSCH.MappingType = 'A';     % PDSCH mapping type ('A'(slot-wise),'B'(non slot-wise))

% Scrambling identifiers
simParameters.PDSCH.NID = simParameters.Carrier.NCellID;
simParameters.PDSCH.RNTI = 1;

% PDSCH resource block mapping (TS 38.211 Section 7.3.1.6)
simParameters.PDSCH.VRBToPRBInterleaving = 0; % Disable interleaved resource mapping
simParameters.PDSCH.VRBBundleSize = 4;

% Define the number of transmission layers to be used
simParameters.PDSCH.NumLayers = 1;            % Number of PDSCH transmission layers

% Define codeword modulation and target coding rate
% The number of codewords is directly dependent on the number of layers so ensure that
% layers are set first before getting the codeword number
if simParameters.PDSCH.NumCodewords > 1                             % Multicodeword transmission (when number of layers being > 4)
    simParameters.PDSCH.Modulation = {'16QAM','16QAM'};             % 'QPSK', '16QAM', '64QAM', '256QAM'
    simParameters.PDSCHExtension.TargetCodeRate = [490 490]/1024;   % Code rate used to calculate transport block sizes
else
    simParameters.PDSCH.Modulation = 'QPSK';                       % 'QPSK', '16QAM', '64QAM', '256QAM'
    simParameters.PDSCHExtension.TargetCodeRate = 490/1024;         % Code rate used to calculate transport block sizes
end

% DM-RS and antenna port configuration (TS 38.211 Section 7.4.1.1)
simParameters.PDSCH.DMRS.DMRSPortSet = 0:simParameters.PDSCH.NumLayers-1; % DM-RS ports to use for the layers
simParameters.PDSCH.DMRS.DMRSTypeAPosition = 2;      % Mapping type A only. First DM-RS symbol position (2,3)
simParameters.PDSCH.DMRS.DMRSLength = 1;             % Number of front-loaded DM-RS symbols (1(single symbol),2(double symbol))
simParameters.PDSCH.DMRS.DMRSAdditionalPosition = 2; % Additional DM-RS symbol positions (max range 0...3)
simParameters.PDSCH.DMRS.DMRSConfigurationType = 2;  % DM-RS configuration type (1,2)
simParameters.PDSCH.DMRS.NumCDMGroupsWithoutData = 1;% Number of CDM groups without data
simParameters.PDSCH.DMRS.NIDNSCID = 1;               % Scrambling identity (0...65535)
simParameters.PDSCH.DMRS.NSCID = 0;                  % Scrambling initialization (0,1)

% PT-RS configuration (TS 38.211 Section 7.4.1.2)
simParameters.PDSCH.EnablePTRS = 0;                  % Enable or disable PT-RS (1 or 0)
simParameters.PDSCH.PTRS.TimeDensity = 1;            % PT-RS time density (L_PT-RS) (1, 2, 4)
simParameters.PDSCH.PTRS.FrequencyDensity = 2;       % PT-RS frequency density (K_PT-RS) (2 or 4)
simParameters.PDSCH.PTRS.REOffset = '00';            % PT-RS resource element offset ('00', '01', '10', '11')
simParameters.PDSCH.PTRS.PTRSPortSet = [];           % PT-RS antenna port, subset of DM-RS port set. Empty corresponds to lower DM-RS port number

% Reserved PRB patterns, if required (for CORESETs, forward compatibility etc)
simParameters.PDSCH.ReservedPRB{1}.SymbolSet = [];   % Reserved PDSCH symbols
simParameters.PDSCH.ReservedPRB{1}.PRBSet = [];      % Reserved PDSCH PRBs
simParameters.PDSCH.ReservedPRB{1}.Period = [];      % Periodicity of reserved resources

% Additional simulation and DL-SCH related parameters
% PDSCH PRB bundling (TS 38.214 Section 5.1.2.3)
simParameters.PDSCHExtension.PRGBundleSize = [];     % 2, 4, or [] to signify "wideband"

% HARQ process and rate matching/TBS parameters
simParameters.PDSCHExtension.XOverhead = 6*simParameters.PDSCH.EnablePTRS; % Set PDSCH rate matching overhead for TBS (Xoh) to 6 when PT-RS is enabled, otherwise 0
simParameters.PDSCHExtension.NHARQProcesses = 16;    % Number of parallel HARQ processes to use
simParameters.PDSCHExtension.EnableHARQ = true;      % Enable retransmissions for each process, using RV sequence [0,2,3,1]

% LDPC decoder parameters
% Available algorithms: 'Belief propagation', 'Layered belief propagation', 'Normalized min-sum', 'Offset min-sum'
simParameters.PDSCHExtension.LDPCDecodingAlgorithm = 'Normalized min-sum';
simParameters.PDSCHExtension.MaximumLDPCIterationCount = 6;

% Define the overall transmission antenna geometry at end-points
% If using a CDL propagation channel then the integer number of antenna elements is
% turned into an antenna panel configured when the channel model object is created
simParameters.NTxAnts = 1;                        % Number of PDSCH transmission antennas (1,2,4,8,16,32,64,128,256,512,1024) >= NumLayers
if simParameters.PDSCH.NumCodewords > 1           % Multi-codeword transmission
    simParameters.NRxAnts = 8;                    % Number of UE receive antennas (even number >= NumLayers)
else
    simParameters.NRxAnts = 1;                    % Number of UE receive antennas (1 or even number >= NumLayers)
end

% % Define the general CDL/TDL propagation channel parameters
% simParameters.DelayProfile = 'CDL-C';   % Use CDL-C model (Urban macrocell model)
% simParameters.DelaySpread = 300e-9;
% simParameters.MaximumDopplerShift = 5;

% waveformInfo = nrOFDMInfo(simParameters.Carrier); % Get information about the baseband waveform after OFDM modulation step

% Set up redundancy version (RV) sequence for all HARQ processes
if simParameters.PDSCHExtension.EnableHARQ
    rvSeq = [0 2 3 1];
else
    % HARQ disabled - single transmission with RV=0, no retransmissions
    rvSeq = 0;
end

% Create DL-SCH encoder system object to perform transport channel encoding
encodeDLSCH = nrDLSCH;
encodeDLSCH.MultipleHARQProcesses = true;
encodeDLSCH.TargetCodeRate = simParameters.PDSCHExtension.TargetCodeRate;

% Create DL-SCH decoder system object to perform transport channel decoding
% Use layered belief propagation for LDPC decoding, with half the number of
% iterations as compared to the default for belief propagation decoding
decodeDLSCH = nrDLSCHDecoder;
decodeDLSCH.MultipleHARQProcesses = true;
decodeDLSCH.TargetCodeRate = simParameters.PDSCHExtension.TargetCodeRate;
decodeDLSCH.LDPCDecodingAlgorithm = simParameters.PDSCHExtension.LDPCDecodingAlgorithm;
decodeDLSCH.MaximumLDPCIterationCount = simParameters.PDSCHExtension.MaximumLDPCIterationCount;

% Take full copies of the simulation-level parameter structures so that they are not
% PCT broadcast variables when using parfor
simLocal = simParameters;

% Take copies of channel-level parameters to simplify subsequent parameter referencing
carrier = simLocal.Carrier;
pdsch = simLocal.PDSCH;
pdschextra = simLocal.PDSCHExtension;
decodeDLSCHLocal = decodeDLSCH;  % Copy of the decoder handle to help PCT classification of variable
decodeDLSCHLocal.reset();        % Reset decoder at the start of each SNR point

% Specify the fixed order to cycle through the HARQ process IDs
harqSequence = 0:pdschextra.NHARQProcesses-1;

% Initialize the state of all HARQ processes
harqEntity = HARQEntity(harqSequence,rvSeq,pdsch.NumCodewords);

% Total number of slots in the simulation period
NSlots = simLocal.NFrames * carrier.SlotsPerFrame;


N=waveformInfo.Nfft;
CP_len=waveformInfo.CyclicPrefixLengths;
NP=min(CP_len);
NNP=N+NP; %% block length (OFDM symbol plus CP)
K=3240;
M=280;

%% Set system parameters
fc=40000; %% center frequency (MHz) (802.11ad at channel 2: 60.48GHz)
w_bw=100; %% bandwidth (MHz) (11ad at 60GHz)
sample_rate=waveformInfo.SampleRate;
samp_period=1/sample_rate; %% sampling period in seconds (overs-ampling factor 8/7)


%% Some parameters for radar algorithms and radar resolutions
K1=1; %% the IFFT size for range estimation is K1*N (K1=1)
K2=8; %% the FFT size for Doppler estimation is K2*M (larger K2 is better for single target. K2=16)
NK=N*K1; %% IFFT size for range estimation
MK=M*K2; %% FFT size for Doppler estimation

MKNNP_samp_fc = MK*NNP*samp_period/150*fc;

%% Maximum allowable range and speed (speed can be positive or negative)
f_allow=0.5/(NNP*samp_period); %% Maximum allowable Doppler frequency (Hz)
v_allow=f_allow*540/fc; %% Maximum allowable speed of the target (km/h)
v_allow_mps=v_allow/3.6; %% Maximum allowable speed of the target (m/s)
t_allow=(NP-8)*samp_period; %% Maximum allowable round-trip time (CP minus the filter) (second)
r_allow=t_allow*(3e8); %% Maximum allowable round-trip range (meter)

%% Set your own limits for round-trip delay/range and Doppler/speed
v_max=min(v_allow,400); %% Maximum speed of the target (km/h)
f_max=2*fc*v_max*10/(3*3600); %% Maximum Doppler frequency (Hz)
t_max=min(t_allow,1e-6); %% Maximum round-trip delay time (second)
r_max=t_max*(3e8); %% Maximum round-trip range (meter)

fd_min=1/(MK*NNP*samp_period); %% Minimum Doppler frequency (Hz)

t_min=1e-9; %% Minimum round-trip delay time (second)

%% empirical half width of the mainlobe is 4 in sampling periods
v_min=fd_min*540/fc; %% Minimum speed of the target (km/h)
r_min=t_min*(3e8); %% Minimum round-trip range (meter)

tic

%% Compute the resolutions (some of resolutions here means resolvable accuracy for single target)
%% For multiple targets, resolution may mean resolvable distance between different targets
Doppler_eff=NNP*samp_period*f_max; %% the effect of Doppler within one OFDM block (normalized)
Doppler_res=1/(2*M*NNP*samp_period);  %% the Doppler resolution by FFT method (Hz)
Doppler_res_largeFFT=1/(2*MK*NNP*samp_period);  %% the Doppler resolution by large FFT method (Hz)
Speed_res=Doppler_res*(150/fc);  %% the speed resolution by FFT method (m/s)
Speed_res_largeFFT=Doppler_res_largeFFT*(150/fc);  %% the speed resolution by large FFT method (m/s)
Range_res=1e6*samp_period*300/2;  %% the round-trip range resolution by FFT method (m)
Range_res_largeFFT=1e6*samp_period*300/(2*K1);  %% the round-trip range resolution by FFT method (m)

w_speed=ones(1,M); %% window for speed (rectangular-no window, better in general)
%% ambiguity function for speed is ideal (K2=1 case) if no window (rectangular window) is used


%% Cancel vector for self-interference cancellation
%% For range
fw_range=abs(ifft(ones(1,N),NK)); %% ifft of the window vector for range
%fw_range=abs(ifft(w_range,NK)); %% ifft of the window vector for range
fw_range=fw_range/fw_range(1); %% normalized
%fw_range=fw_range/max(fw_range); %% normalized
LL2=1; %% filter length in symbols (short filter. need to check the filter)

d_cali=K1*(LL2-1)/2; %% calibrated delay samples (verified by simulations)
%r_cancel=[fw_range(NK-d_cali+1:NK),fw_range(1:NK/2)].^2; %% for range and power
%r_cancel=r_cancel.';

%% For speed
fw_speed=fft(w_speed,MK); %% fft of the window vector for speed
fw_speed=abs(fw_speed/fw_speed(1)); %% normalized
v_cancel=(fw_speed.').^2; %% for speed and power (better)

firstSC = (N/2) - (K/2) + 1;
B_null = [1:(firstSC-1) (firstSC+K):N];
B_data=[firstSC:(firstSC+K)-1];

range_result=zeros(1,INUM);
speed_result=zeros(1,INUM);
U=B_data;
s=ResourceGrids;
angle_result = zeros(1, INUM);

r_array=((0:NK-1).'/K1)*(10^6*samp_period)*150; %% range
v_array=(-MK/2:MK/2-1).'*150/(MK*NNP*fc*samp_period);

Orig1 = resample(Wave1, 1228800, 2000000);

%% Estimate the range and speed from 2D vector
%% Find the maximum sample number for Doppler frequency and range
M_start=floor(fd_min*(MK*NNP*samp_period)); %% start point of Doppler frequency (excluding the self-interference)
N_start=floor(t_min*K1/samp_period); %% excluding the self-interference
f_max=2*fc*v_max*10/(3*3600); %% the considered maximum Doppler frequency
M_max=max(ceil(f_max*MK*NNP*samp_period)+3,M_start+1); %% the maximum sample number for Doppler
N_max=ceil(K1*t_max/samp_period)+3; %% the maximum sample number for range

N_sqrt = sqrt(N);


%% Main Tests
if runTests
    [~, respiratorTemp] = RunTest(strCOM,sensor_test);
    image_encoder(img_filename);


    f = fopen("~/Desktop/matlab/images/imgbitstream2", 'r');
    bitstream = fscanf(f, '%c'); %read bitstream as chars
    img_orig_length = length(bitstream);
    img_padded_length = ceil(length(bitstream)/NSlots) * NSlots;

    % converts char array to doubles array (padded with 0's)
    a = zeros(1, img_padded_length);
    for i=1:img_orig_length
        a(i) = str2double(bitstream(i));
    end
    img_input=reshape(a, img_padded_length/NSlots, NSlots); %divide the bitstream by the number of slots

    wtx = 1;
    [carrier, WaveImg1, dmrsAntSymbols, dmrsAntIndices, dmrsSymbols_all, dmrsIndices_all, trBlk_tx, trBlkSizes, pdsch_all, pdschIndices, pdschIndices_all, encodeDLSCH] = img_propagation(img_padded_length, img_input, NSlots, carrier, simLocal, pdsch, encodeDLSCH, decodeDLSCHLocal, pdschextra, harqEntity);
end



%% Calculate Syncronisation Index
v=read_complex_binary(['/media/ramdisk/usrp_replay_capture' num2str(0) '.dat']);
testMod=resample(v, 1228800, 2000000);
a=fft(testMod(1:1+length(Orig1)-1));
c=zeros(1228800,1);
b=fft(Orig1);
for i=1:1228800
    c(i)= a(i)*b(i)';
end
CorrVal=ifft(c);
index=find(CorrVal == max(CorrVal));

% Range-Doppler Map calcuation (based on one data snippet)
v=read_complex_binary(['/media/ramdisk/usrp_replay_capture' num2str(0) '.dat']);
testMod=resample(v, 1228800, 2000000);

ch_vec_fft_FCCR_all = zeros(NK,MK,N_antenna);
ch_vec_fft_orig=zeros(NK,MK);
ch_vector_orig=zeros(N,M);

s_full=zeros(N,M);
s_full(B_data,:)=s;
s1 = zeros(4096, 280);

a=fft(testMod(1:1+length(Orig1)-1));
c=zeros(1228800,1);
b=fft(Orig1);
for i=1:1228800
    c(i)= a(i)*b(i)';
end
CorrVal=ifft(c);
index=find(CorrVal == max(CorrVal));
s01=testMod(index:index+length(Orig1)-1);

e_l_LOS=0;
v_l_LOS=0;

current_pos=0;
for m=1:M
    cp_len_pos = mod(m,14);
    if (cp_len_pos==0)
        cp_len_pos = 14;
    end
    s1(:,m)=fft(s01(current_pos+CP_len(cp_len_pos)+1:current_pos+CP_len(cp_len_pos)+4096))/sqrt(N);
    current_pos = current_pos + CP_len(cp_len_pos) + 4096;
    s1(:,m) = circshift(s1(:,m),floor(size(s1(:,m))/2));
end

for m=1:M
    tem=zeros(N,1);
    tem(U,:)=s1(U,m).*conj((s_full(U,m))); %% for cyclic convolution
    ch_ave_r=0;
    for n=U
        ch_ave_r=ch_ave_r+tem(n)*exp(cj*2*pi*e_l_LOS*n/N)*exp(-cj*2*pi*((m-1)*NNP+NP+n)*v_l_LOS);
    end
    tem_cancel=zeros(N,1 );
    for n=U
        tem_cancel(n)=tem(n)-ch_ave_r/length(U)*exp(-cj*2*pi*e_l_LOS*n/N)*exp(cj*2*pi*((m-1)*NNP+NP+n)*v_l_LOS);
    end
    ch_vector_orig(:,m)=ifft(tem_cancel); %% time domain matched filtering (cyclic correlation)
end


for mm=1:NK
    ch_vec_fft_orig(mm,:)=sqrt(1/MK)*fft(ch_vector_orig(mm,:).',MK).'; %% FFT on the differential correlation
end

ch_vec_fft=ch_vec_fft_orig;
% normalises the graph?
ch_vec_fft=abs(ch_vec_fft);
ch_vec_fft=ch_vec_fft/max(max(ch_vec_fft));

        [row_max1,row_index1]=max(abs(ch_vec_fft(N_start+1:N_max,M_start+1:M_max)));
        %% find the maximum value and index of each column (for positive speed)
        [row_max2,row_index2]=max(abs(ch_vec_fft(N_start+1:N_max,MK-M_max+1:MK-M_start)));
        %% find the maximum value and index of each column (for negative speed)
        [total_max1,total_index1]=max(row_max1); %% find the maximum value in the row
        [total_max2,total_index2]=max(row_max2); %% find the maximum value in the row
        if total_max1>=total_max2
            range_index=row_index1(total_index1); %% range index of the maximum point
        else
            range_index=row_index2(total_index2); %% range index of the maximum point
        end

        %% Range and speed estimation by using the maximum point of the 2D-FFT
        %% This gives the best performance (ML estimations)
        %% Speed estimation
        if total_max1>=total_max2 %% speed is positive
            v_est=(total_index1-1+M_start-1)/(MK*NNP*samp_period)*150/fc; %% the estimated speed (m/s)
        else %% speed is negative
            v_est=(-M_max-1+total_index2)/(MK*NNP*samp_period)*150/fc; %% the estimated speed (m/s) (verified)
        end

        %% Range estimation
        %% The calibrated results (calibrated with average range error almost zero for all K1)
        %% Due to filter and oversize IFFT, calibrating is necessary
        peak_points=((range_index-1)/K1)*(10^6*samp_period)...
            +N_start*(10^6*samp_period)/K1; %% estimated round-trip delays (us)
        %    r_error4(k_snr)=r_error4(k_snr)+abs(peak_points-ch_delay(2)*10^6)*300; %% absolute range error (meter)

        %% the estimation for the second largest target
        range=peak_points*150;
x_ar=ch_vec_fft; %% the range-Doppler map
for n=1:NK
    x_ar(n,1:MK/2)=ch_vec_fft(n,MK/2+1:MK);
    x_ar(n,MK/2+1:MK)=ch_vec_fft(n,1:MK/2);
end
% figure;
% mesh(v_array, r_array(1:20),abs(x_ar(1:20,:)).^2)
% axis('tight');
% grid on;
% xlabel('Speed (m/s)','fontsize',12);
% ylabel('Range (m)','fontsize',12);
% title('JRC RDM');
% textStr1 = sprintf('range = %f m', range); 
% textStr2 = sprintf('speed = %f m/s', v_est); 
% z=max(max(x_ar(1:20,:)));
% txt = {[textStr1], [textStr2]};
%     text(v_est,range, z/2,txt, 'FontSize', 12);

% Perform signal analysis
tic
disp('Running signal analysis...');
parfor NUM=1:INUM
%    for NUM=1:INUM
    [angle,angle_1,angle_2,speed,range]=signalAnalysis(NUM, Orig1, index, NK, MK, NNP, fc, N_antenna, M, s, B_data, CP_len, N_sqrt, N, U, M_start, M_max, N_start, N_max, MKNNP_samp_fc, K1, samp_period, cj, NP);
    angle_result(NUM)=angle;
            angle_result_1(NUM)=angle_1;
        angle_result_2(NUM)=angle_2;
        range_result(NUM)=range;
        speed_result(NUM)=speed;
end
toc




%%%% Graph Plots %%%%
%%  Antenna results  %%
x1=angle_result;
x1=unwrap(x1);
x1=x1-mean(x1);
%offset=20; x1=x1(offset:end); % removes the first few signals
%x1=x1-mean(x1);
x2=conv(x1,hf,'same');
x3=conv(x1,hf_hb,'same');
% Reordering the speed dimension to cater for positive and negative speed





% screen_size = get(0, 'ScreenSize');
% f1 = tiledlayout(2,4);
% set(f1, 'Position', [0 0 screen_size(3) screen_size(4) ] );
if (hearbeat_test==true)
    tiledlayout(2,5);
elseif (sensor_test==true)
tiledlayout(2,4);
else
    tiledlayout(2,3);
end
nexttile;
%fitiledlayout(2,3);gure;
mesh(v_array, r_array(1:20),abs(x_ar(1:20,:)).^2)
axis('tight');
grid on;
xlabel('Speed (m/s)','fontsize',12);
ylabel('Range (m)','fontsize',12);
title('JRC RDM');
textStr1 = sprintf('range = %f m', range); 
textStr2 = sprintf('speed = %f m/s', v_est); 
z=max(max(x_ar(1:20,:)));
txt = {[textStr1], [textStr2]};
    text(v_est,range, z/2,txt, 'FontSize', 12);

%figure;
nexttile;
%plot(0.06:0.06:30,x1);
plot(ts:ts:INUM*ts,x1);
xlabel('Time (s)'); ylabel('Phase Change');
title('JRC raw breath signal');

%figure;
if speed_test==true
nexttile;
%plot(0.06:0.06:30,x2);
plot(ts:ts:INUM*ts,speed_result);
xlabel('Time (s)'); ylabel('Speed m/s');
title('JRC Speed Result');
else
nexttile;
%plot(0.06:0.06:30,x2);
plot(ts:ts:INUM*ts,x2);
xlabel('Time (s)'); ylabel('Phase Change');
title('JRC breath signal after a filter');
end


%% Plot the spectrum of the filtered signal
N_fft=8192; %% FFT size
ND=1000;
spec_x1=fftshift(fft(x2,N_fft))/sqrt(ND);
%fs=16.6667;
fs=1/ts;
f_rs=fs; %% the sampling rate (Hz)
desLen=N_fft;

%figure;
nexttile;
x_axis=(-desLen/2:desLen/2-1)*f_rs/desLen;
% x_axis=(0:desLen-1)*f_rs/desLen;
semilogy(x_axis, abs(spec_x1).^2)
% plot(x_axis, abs(spec_x1).^2)
xlabel('Frequency (Hz)');
ylabel('PSD');
axis tight;
title('JRC: PSD of the breath signal ');



toc

%% Find the peaks of the signal
f_start=0.1; %% the start frequency point for searching (Hz)
MS=floor(f_start/f_rs*desLen); %% the number of frequency points to be ignored
y_pf=abs(spec_x1(N_fft/2+MS:N_fft)).^2; %% the signal PSD at positive frequency
[pks,locs] = max(y_pf);

%% Calculate the frequency of the peak point and rate per minute
disp('JRC results: ')
p_freq_breath=(locs+MS-2)*f_rs/N_fft;
p_rate_breath_JRC=p_freq_breath*60;  %% rate per minute


% nexttile;
% plot(0.06:0.06:18,x3);
% xlabel('Time (s)'); ylabel('Phase Change');
% title('JRC heart signal after a filter');

%% Plot the spectrum of the filtered signal
if (hearbeat_test==true)
N_fft=8192; %% FFT size
%ND=1000;
spec_x3=fftshift(fft(x3,N_fft))/sqrt(ND);
%fs=16.6667;
fs=1/ts;
f_rs=fs; %% the sampling rate (Hz)
desLen=N_fft;

%figure;
nexttile;
x_axis=(-desLen/2:desLen/2-1)*f_rs/desLen;
% x_axis=(0:desLen-1)*f_rs/desLen;
semilogy(x_axis, abs(spec_x3).^2)
% plot(x_axis, abs(spec_x1).^2)
xlabel('Frequency (Hz)');
ylabel('PSD');
axis tight;
title('JRC: PSD of the heart beat signal ');

%% Find the peaks of the signal
f_start=0.1; %% the start frequency point for searching (Hz)
MS=floor(f_start/f_rs*desLen); %% the number of frequency points to be ignored
y_pf=abs(spec_x3(N_fft/2+MS:N_fft)).^2; %% the signal PSD at positive frequency
[pks,locs] = max(y_pf);

%% Calculate the frequency of the peak point and rate per minute
disp('JRC results: ')
p_freq_heart=(locs+MS-2)*f_rs/N_fft;
p_rate_heart_JRC=p_freq_heart*60;  %% rate per minute
end
% % 
%  fh = figure(); 
% % % add axis to bottom, left corner
%  ah = axes(fh,'position',[0,0,1,1]); 
% % add text
% % x.bytes = 9000; 
% % filesizeini = s.bytes;
% %figure
% p_rate_heart_JRC=1;
% p_rate_heart_Sensor=2;
% textStr1 = sprintf('JRC breath rate = %f per min', p_rate_heart_JRC); 
% textStr2 = sprintf('Sensor breath rate = %f per min', p_rate_heart_Sensor); 
% txt = {[textStr1], [textStr2]};
% %    text(0.5,0.5,txt, 'FontSize', 12);
% %text(ah,.2,.2,txt,'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle')
% text(ah,.3,.5,txt)
% % Turn off axis
% ah.Visible = 'off';


if sensor_test==true
    %%  Respirator results  %%
    x2_sensor=respiratorTemp;
x2_sensor=conv(x2_sensor,hf_sensor,'same');
    
    ND=INUM;
    spec_x1_sensor=fftshift(fft(x2_sensor,N_fft))/sqrt(ND);
    fs=10;
    f_rs=fs; %% the sampling rate (Hz)
    x_axis=(-desLen/2:desLen/2-1)*f_rs/desLen;
    % x_axis=(0:desLen-1)*f_rs/desLen;
    % plot(x_axis, abs(spec_x1).^2)\
   % figure;
   nexttile;
    semilogy(x_axis, abs(spec_x1_sensor).^2)
    xlabel('Frequency (Hz)');
    ylabel('PSD');
    axis tight;
    title('Sensor: PSD of the breath signal ');

 %   figure;
 nexttile;
    plot(0.1:0.1:60,respiratorTemp);
    xlabel('Time (s)')
    ylabel('Power')
    title('Respiration Monitor');


    %% Find the peaks of the signal
    f_start=0.1; %% the start frequency point for searching (Hz)
    MS=floor(f_start/f_rs*desLen); %% the number of frequency points to be ignored
    y_pf=abs(spec_x1_sensor(N_fft/2+MS:N_fft)).^2; %% the signal PSD at positive frequency
    [pks,locs] = max(y_pf);

    %% Calculate the frequency of the peak point and rate per minute
    disp("Respirator values: ");
    p_freq_breath=(locs+MS-2)*f_rs/N_fft;
    p_rate_breath_Sensor=p_freq_breath*60 ; %% rate per minute


textStr1 = sprintf('JRC breath rate = %f times per min', p_rate_breath_JRC); 
textStr2 = sprintf('Sensor breath rate = %f times per min', p_rate_breath_Sensor); 
textStr3 = sprintf('range = %f m', range); 
textStr4 = sprintf('speed = %f m/s', v_est); 
textStr5 = sprintf('JRC heart rate = %f times per min', p_rate_heart_JRC); 
%txt = {[textStr3], [textStr4],[textStr1], [textStr2]};
%txt = {[textStr1], [textStr2],[textStr5]};
txt = {[textStr1], [textStr2]};
if (hearbeat_test==true)
    txt = {[textStr1], [textStr2],[textStr5]};


 nexttile;
    plot(p_rate_breath_Sensor,'*'); hold on
    plot(p_rate_breath_JRC,'-ro'); hold on
    plot(p_rate_heart_JRC,'-gd'); hold off
    xlabel('Test Number')
    ylabel('Breath Rate Per Minute')
    legend('Respirator Sensor','JRC breath rate Detection','JRC heartbeat rate Detection');
 text(0.2,p_rate_breath_Sensor,txt, 'FontSize', 8);

end

 %   figure;
 nexttile;
 %    plot(p_rate_breath_Sensor,'*'); hold on
 %    plot(p_rate_breath_JRC,'-ro'); hold off
 %    xlabel('Test Number')
 %    ylabel('Breath Rate Per Minute')
 %    legend('Respirator Sensor','JRC Detection');
 % text(0.2,p_rate_breath_Sensor,txt, 'FontSize', 8);



h1 = animatedline;
h2 = animatedline;
h3 = animatedline;
h1.Color = 'r';
h2.Color = 'g';
h3.Color = 'm';

% x = linspace(0,4*pi,10000);
% y = sin(x);

xlabel('seconds(s)')
ylabel('breath rate per min')
legend('JRC','Sensor')
%legend('JRC','Sensor','heartbeat')
%v=VideoWriter('breathrate.avi');
%open(v);

   NUM2=300;
   len=300;
   axis([len*ts INUM*ts 0 100])
for NUM=len+1:INUM
tic

   %    [angle,speed,range]=signalAnalysis(NUM, Orig1, index, NK, MK, NNP, fc, N_antenna, M, s, B_data, CP_len, N_sqrt, N, U, M_start, M_max, N_start, N_max, MKNNP_samp_fc, K1, samp_period, cj, NP);
   % angle_result(NUM)=angle;

x1=angle_result(NUM-len:NUM);
x1=unwrap(x1);
x1=x1-mean(x1);
%x1=x1-mean(x1);
x2=conv(x1,hf,'same');

%% Plot the spectrum of the filtered signal
N_fft=8192; %% FFT size
%ND=1000;
spec_x1=fftshift(fft(x2,N_fft))/sqrt(ND);
%fs=16.6667;
%fs=1/0.06;
fs=1/ts;
f_rs=fs; %% the sampling rate (Hz)
desLen=N_fft;

%% Find the peaks of the signal
f_start=0.1; %% the start frequency point for searching (Hz)
MS=floor(f_start/f_rs*desLen); %% the number of frequency points to be ignored
y_pf=abs(spec_x1(N_fft/2+MS:N_fft)).^2; %% the signal PSD at positive frequency
[pks,locs] = max(y_pf);
p_freq_breath=(locs+MS-2)*f_rs/N_fft;
p_rate_breath_JRC=p_freq_breath*60 ; %% rate per minute
    addpoints(h1,NUM*ts,p_rate_breath_JRC);

if (hearbeat_test==true)
legend('JRC','Sensor','heartbeat');
x3=conv(x1,hf_hb,'same');
spec_x3=fftshift(fft(x3,N_fft))/sqrt(ND);
%% Find the peaks of the signal
f_start=0.1; %% the start frequency point for searching (Hz)
MS=floor(f_start/f_rs*desLen); %% the number of frequency points to be ignored
y_pf=abs(spec_x3(N_fft/2+MS:N_fft)).^2; %% the signal PSD at positive frequency
[pks,locs] = max(y_pf);

%[x1,y1]=sort(x);
i=0;
%p_freq_heart=(y(y1(end-i))+MS-2)*f_rs/N_fft;
% while ((abs(p_freq_heart/p_freq_breath)<3.1 && abs(p_freq_heart/p_freq_breath)>3) || (abs(p_freq_heart/p_freq_breath)<4.1 && abs(p_freq_heart/p_freq_breath)>4))
% i=i+1;
% p_freq_heart=(y(y1(end-i))+MS-2)*f_rs/N_fft;
% end


p_freq_heart=(locs+MS-2)*f_rs/N_fft;



p_rate_heart_JRC=p_freq_heart*60;  %% rate per minute

    addpoints(h3,NUM*ts,p_rate_heart_JRC);

end




    %toc

%if (NUM>=500)
 NUM2=floor(NUM*ts/0.1);
if NUM2>1
   % x2_sensor=respiratorTemp(NUM2-300+1:NUM2);;
     x2_sensor=respiratorTemp(1:NUM2);;
x2_sensor=conv(x2_sensor,hf_sensor,'same');
    
    ND=INUM;
    spec_x1_sensor=fftshift(fft(x2_sensor,N_fft))/sqrt(ND);
    fs=10;
    f_rs=fs; %% the sampling rate (Hz)
        %% Find the peaks of the signal
    f_start=0.1; %% the start frequency point for searching (Hz)
    MS=floor(f_start/f_rs*desLen); %% the number of frequency points to be ignored
    y_pf=abs(spec_x1_sensor(N_fft/2+MS:N_fft)).^2; %% the signal PSD at positive frequency
    [pks,locs] = max(y_pf);

    %% Calculate the frequency of the peak point and rate per minute
%    disp("Respirator values: ");
    p_freq_breath=(locs+MS-2)*f_rs/N_fft;
    p_rate_breath_Sensor=p_freq_breath*60;  %% rate per minute

    addpoints(h2,NUM2*0.1,p_rate_breath_Sensor);
  %  drawnow limitrate
end

    drawnow limitrate
    % frame=getframe(gcf);
    % writeVideo(v,frame)
end
drawnow




    % 
%  fh = figure(); 
% 
%  ah = axes(fh,'position',[0,0,1,1]); 
% 
% txt = {[textStr3], [textStr4],[textStr1], [textStr2]};
% %    text(0.5,0.5,txt, 'FontSize', 12);
% %text(ah,.2,.2,txt,'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle')
% text(ah,.3,.5,txt)
% % Turn off axis
% ah.Visible = 'off';

else


%  fh = figure(); 
% 
%  ah = axes(fh,'position',[0,0,1,1]); 
% % p_rate_heart_JRC=1;
% % p_rate_heart_Sensor=2;
% textStr1 = sprintf('JRC breath rate = %f per min', p_rate_heart_JRC); 
% txt = {[textStr3], [textStr4],[textStr1]};
% text(ah,.3,.5,textStr1)
% % Turn off axis
% ah.Visible = 'off';
textStr1 = sprintf('JRC breath rate = %f per min', p_rate_breath_JRC); 
%textStr2 = sprintf('Sensor breath rate = %f per min', p_rate_heart_Sensor); 
textStr3 = sprintf('range = %f m', range); 
textStr4 = sprintf('speed = %f m/s', v_est); 
txt = {[textStr3], [textStr4],[textStr1]};
 nexttile;
   % plot(p_rate_heart_Sensor,'*'); hold on
    plot(p_rate_heart_JRC,'-ro'); hold off
    xlabel('Test Number')
    ylabel('Breath Rate Per Minute')
    legend('JRC Detection');
 text(0.2,p_rate_heart_JRC-2,txt, 'FontSize', 12);


end
%% Image Reconstruction %%
disp("Reconstructing image...");
pathFilters = [];
% Timing offset, updated in every slot for perfect synchronization and
% when the correlation is strong for practical synchronization

v=read_complex_binary('/media/ramdisk/usrp_replay_capture9999.dat');
Orig1 = resample(WaveImg1, 1228800, 2000000);
testMod=resample(v, 1228800, 2000000);

a=fft(testMod(1:1+length(Orig1)-1));
c=zeros(1228800,1);
b=fft(Orig1);
for i=1:1228800
    c(i)= a(i)*b(i)';
end
CorrVal=ifft(c);
index=find(CorrVal == max(CorrVal));
% index=39;
y=testMod(index:index+length(Orig1)-1);

s01= y;
current_pos=0;
for m=1:M
    cp_len_pos = mod(m,14);
    if (cp_len_pos==0)
        cp_len_pos = 14;
    end
    s1(:,m)=fft(s01(current_pos+CP_len(cp_len_pos)+1:current_pos+CP_len(cp_len_pos)+4096))/sqrt(N);
    current_pos = current_pos + CP_len(cp_len_pos) + 4096;
    s1(:,m) = circshift(s1(:,m),floor(size(s1(:,m))/2));
end

rxWaveform_all=s01;
estChannel=[];
rx=[];
trBlk_rx=[];
dmrsSymbols_rx=[];

for slot=1:waveformInfo.SlotsPerFrame
    resetSoftBuffer(decodeDLSCHLocal,0,harqEntity.HARQProcessID);
    pdsch=pdsch_all{slot};
    pdschIndices=pdschIndices_all(:,slot);
    dmrsSymbols = dmrsSymbols_all(:,slot);
    dmrsIndices = dmrsIndices_all(:,slot);
    % dmrsAntSymbols=dmrsAntSymbols_all(:,slot);
    rxWaveform_perslot = rxWaveform_all((slot-1)*61440+1:(slot-1)*61440+61440);

    offset=0;
    rxWaveform = rxWaveform_perslot(1+offset:end,:);
    rxGrid = nrOFDMDemodulate(carrier,rxWaveform);
    [K,L,R] = size(rxGrid);
    if (L < carrier.SymbolsPerSlot)
        rxGrid = cat(2,rxGrid,zeros(K,carrier.SymbolsPerSlot-L,R));
    end

    if (simLocal.PerfectChannelEstimator)
        estChannelGrid = nrPerfectChannelEstimate(carrier,pathGains,pathFilters,offset,sampleTimes);

        % Get perfect noise estimate (from the noise realization)
        noiseGrid = nrOFDMDemodulate(carrier,noise(1+offset:end ,:));
        noiseEst = var(noiseGrid(:));

        % Get PDSCH resource elements from the received grid and channel estimate
        [pdschRx,pdschHest,~,pdschHestIndices] = nrExtractResources(pdschIndices,rxGrid,estChannelGrid);

        % Apply precoding to channel estimate
        pdschHest = hPRGPrecode(size(estChannelGrid),carrier.NStartGrid,pdschHest,pdschHestIndices,permute(wtx,[2 1 3]));
    else
        % Practical channel estimation between the received grid and each transmission layer, using the PDSCH DM-RS for each
        % layer. This channel estimate includes the effect of transmitter precoding
        [estChannelGrid,noiseEst] = nrChannelEstimate(carrier,rxGrid,dmrsIndices,dmrsSymbols,'CDMLengths',pdsch.DMRS.CDMLengths);

        % [estChannelGrid,noiseEst] = nrChannelEstimate(carrier,rxGrid,dmrsIndices,dmrsAntSymbols,'CDMLengths',pdsch.DMRS.CDMLengths);

        % Get PDSCH resource elements from the received grid and
        % channel estimate
        [pdschRx,pdschHest] = nrExtractResources(pdschIndices,rxGrid,estChannelGrid);
        estChannel=[estChannel estChannelGrid];
    end

    % Equalization
    [pdschEq,csi] = nrEqualizeMMSE(pdschRx,pdschHest,noiseEst);

    rx_tmp=zeros(3240,14);
    rx_tmp(pdschIndices)=pdschEq;
    rx_eq=sign(real(rx_tmp))*0.7071+sign(imag(rx_tmp))*0.7071*cj;
    rx_eq(dmrsAntIndices) = dmrsAntSymbols;
    dmrsSymbols_rx=[dmrsSymbols_rx dmrsAntSymbols];
    rx=[rx rx_eq];

    % Decode PDSCH physical channel
    [dlschLLRs,rxSymbols] = nrPDSCHDecode(carrier,pdsch,pdschEq,noiseEst);

    % Scale LLRs by CSI
    csi = nrLayerDemap(csi); % CSI layer demapping
    for cwIdx = 1:pdsch.NumCodewords
        Qm = length(dlschLLRs{cwIdx})/length(rxSymbols{cwIdx}); % bits per symbol
        csi{cwIdx} = repmat(csi{cwIdx}.',Qm,1);                 % expand by each bit per symbol
        dlschLLRs{cwIdx} = dlschLLRs{cwIdx} .* csi{cwIdx}(:);   % scale by CSI
    end

    % Decode the DL-SCH transport channel
    decodeDLSCHLocal.TransportBlockLength = trBlkSizes;
    [decbits,blkerr] = decodeDLSCHLocal(dlschLLRs,pdsch.Modulation,pdsch.NumLayers,harqEntity.RedundancyVersion,harqEntity.HARQProcessID);

    %      trBlk_rx=[trBlk_rx decbits];
    % HARQ processing
    for cwIdx = 1:pdsch.NumCodewords
        % If new data for current process and codeword then create a new DL-SCH transport block
        if harqEntity.NewData(cwIdx)
            trBlk = decbits;
            setTransportBlock(encodeDLSCH,trBlk,cwIdx-1,harqEntity.HARQProcessID);
        end
    end
    trBlk_rx=[trBlk_rx trBlk];
end

tmp=int8(trBlk_tx)-trBlk_rx;
ber=sum(sum(abs(tmp)))/39936/10;

% Convert the received bitstream back to an image
bitstream = uint8(reshape(trBlk_rx(1:img_padded_length/NSlots, :), 1, []));
bitstream2 = char(bitstream + '0'); % converts unit8 array into a char array

byteArray = zeros(1, img_orig_length / 8);
for i = 1:length(byteArray)
    byteArray(i) = bin2dec(bitstream2((i-1)*8+1:i*8));
end

img_extension = img_filename(regexp(img_filename, '\.'): end);

fid = fopen(['~/Desktop/matlab/test_data/temp' img_extension], 'wb'); % hacky way to convert the bytearray into an image
fwrite(fid, byteArray, 'uint8');
fclose(fid);

img = imread(['~/Desktop/matlab/test_data/temp' img_extension]);
nexttile
imshow(img);
movegui([2100 1000]);

%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% User-defined Functions' Space %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [M, N] = RunTest(strCOM,sensor_test)
funcs = {@rfnoc_test, @sensorReading};
params = {1, strCOM};
vals = cell(1,2);
% Runs both the Antenna and Respirator tests simultaneously
if (sensor_test==true)
    parfor i=1:2
        vals{i} = funcs{i}(params{i});
    end
    M = vals{1}; N = vals{2}; % N contains the values of the respiration monitor
else
    parfor i=1:1
        vals{i} = funcs{i}(params{i});
    end
    M = vals{1}; N = vals{2}; % N contains the values of the respiration monitor
end
disp("Tests completed!");
end




function tmp = sensorReading(strCOM)
% sensor params
if (exist('sensor','var'))
    fclose(sensor); delete(sensor); clear sensor;
end
sensor = serial(strCOM, 'BaudRate', 115200);
fopen(sensor);
% fprintf(sensor,'s'); % stop the sensor
% fprintf(sensor,hex2dec('2F'));
Pt_start_Buf = 1;
Pt_stop_Buf = 15;

% --- sensor calibration for each round ---
% Read the first 100 points to derive max,min,mean values
% Sets the baseline for the remaining data collection
disp('Start getting baseline data for calibration....');
tic
first10 = zeros(100,1);
sensor_raw_all=[];
for i=1:100
    [sensor_raw,count] = fread(sensor,15,'uint8');
    sensor_raw_all=[sensor_raw_all;sensor_raw ];
    if(count<10) break;end;
    [sensor_reading,Pt_start_Buf,Pt_stop_Buf] = ProcessBuffer(sensor_raw,Pt_start_Buf,Pt_stop_Buf);
    first10(i)=sensor_reading;
end
baseline = mean(first10);
range = max(first10) - min(first10);
toc
disp('Finished collecting baseline data.');

% --- start figure ---
% figure('Name','Respiration Monitor','NumberTitle','off');
% DispWindowSize = 30*10;
% plt_wndw = 1:DispWindowSize;
% hCurve = plot(plt_wndw, plt_wndw.*0,'b:');hold on;
% hCurve_Filtered = plot(plt_wndw, plt_wndw.*0,'r');hold off;
% set(gca,'XLimMode','manual','Xlim',[1 DispWindowSize],'YGrid','on');

% Reading data from sensor and detect peaks
tmp=zeros(600,1);
for i=1:600
    [sensor_raw,count] = fread(sensor,100,'uint8');
    if sensor.BytesAvailable >= 100; [sensor_raw_dump,count_dump] = fread(sensor,100,'uint8'); end % clear residual data
    if(count<10); break; end
    [sensor_reading,Pt_start_Buf,Pt_stop_Buf] = ProcessBuffer(sensor_raw,Pt_start_Buf,Pt_stop_Buf); % extract sensor reading from buffer
    tmp(i) = (sensor_reading-baseline) / range;
end

fclose(sensor);
delete(sensor);
clear sensor;
end


function [angle_result,angle_result_1,angle_result_2,speed_result,range_result]= signalAnalysis(NUM, Orig1, index, NK, MK, NNP, fc, N_antenna, M, s, B_data, CP_len, N_sqrt, N, U, M_start, M_max, N_start, N_max, MKNNP_samp_fc, K1, samp_period, cj, NP)

v=read_complex_binary(['/media/ramdisk/usrp_replay_capture' num2str(NUM-1) '.dat']);
testMod=resample(v, 1228800, 2000000);
%CorrVal = zeros(1, length(testMod) - length(Orig1));
% a=fft(testMod(1:1+length(Orig1)-1));
% c=zeros(1228800,1);
% b=fft(Orig1);
% for i=1:1228800
%    c(i)= a(i)*b(i)';
% end
% CorrVal=ifft(c);
% index=find(CorrVal == max(CorrVal));
% index=39;
s01=testMod(index:index+length(Orig1)-1);
s1 = zeros(4096, 280);

ch_vec_fft_FCCR_all = zeros(NK,MK,N_antenna);
%ch_vec_fft_FCCR_all_no_coding = zeros(NK,MK,N_antenna);
%ch_cancel_est_all = zeros(N,M,N_antenna);

current_pos=0;
for m=1:M
    cp_len_pos = mod(m,28);
    if (cp_len_pos==0)
        cp_len_pos = 28;
    end
    s1(:,m)=fft(s01(current_pos+CP_len(cp_len_pos)+1:current_pos+CP_len(cp_len_pos)+4096))/N_sqrt;
    current_pos = current_pos + CP_len(cp_len_pos) + 4096;
    s1(:,m) = circshift(s1(:,m),floor(size(s1(:,m))/2));
end

%% Radar range and speed estimations based on the FCCR method
% ch_vector=zeros(N,M);  %% correlation at all blocks
s_full=zeros(N,M);
s_full(B_data,:)=s;
ch_vector_orig=zeros(N,M);
% ch_vector_ch=zeros(N,M);

e_l_LOS=0;
v_l_LOS=0;

for m=1:M
    tem=zeros(N,1);
    tem(U,:)=s1(U,m).*conj((s_full(U,m))); %% for cyclic convolution
    % ch_vector(:,m)=ifft(tem); %% time domain matched filtering (cyclic correlation)
    %  tem(U,:)=tmp(y)*s_full(U,m).*conj((s_full(U,m))); %% for cyclic convolution
    % tem(U,:)=(s1(U,m)-tmp(y)*s_full(U,m)).*conj((s_full(U,m)));

    % old equation vv
    % ch_vector_orig(:,m)=ifft(tem-sum(tem)/3240); %% time domain matched filtering (cyclic correlation)

    % tem(U,:)=s1(U,m)./((s_full(U,m)));
    % ch_vector_ch(B_data,m)=tem(B_data);

    ch_ave_r=0;
    for n=U
        ch_ave_r=ch_ave_r+tem(n)*exp(cj*2*pi*e_l_LOS*n/N)*exp(-cj*2*pi*((m-1)*NNP+NP+n)*v_l_LOS);
    end
    %            if abs(ch_delay(1)/samp_period-ch_delay(2)/samp_period)<=4
    %                 ch_ave_r=0;
    %            end
    tem_cancel=zeros(N,1 );
    for n=U
        tem_cancel(n)=tem(n)-ch_ave_r/length(U)*exp(-cj*2*pi*e_l_LOS*n/N)*exp(cj*2*pi*((m-1)*NNP+NP+n)*v_l_LOS);
    end

    ch_vector_orig(:,m)=ifft(tem_cancel); %% time domain matched filtering (cyclic correlation)

end


%%  2D search based on the correlations/(differential correlations)
%%  FFT on the speed dimension
% ch_vec_fft=zeros(NK,MK);
% for mm=1:NK
%     ch_vec_fft(mm,:)=sqrt(1/MK)*fft((ch_cancel(mm,:).*w_speed).',MK).'; %% FFT on the differential correlation
% end

ch_vec_fft_orig=zeros(NK,MK);
for mm=1:NK
    %ch_vec_fft_orig(mm,:)=sqrt_of_mk*fft((ch_vector_orig(mm,:).*w_speed).',MK).'; %% FFT on the differential correlation
    ch_vec_fft_orig(mm,:)=sqrt(1/MK)*fft(ch_vector_orig(mm,:).',MK).'; %% FFT on the differential correlation

end
%ch_vec_fft=ch_vec_fft_ch;
ch_vec_fft=ch_vec_fft_orig;
%% Plot the 2D radar profile of the last test
% r_array=((0:NK-1).'/K1)*(10^6*samp_period)*150; %% range
% v_array=(-MK/2:MK/2-1).'*150/(MK*NNP*fc*samp_period);
ch_vec_fft_FCCR=ch_vec_fft;
ch_vec_fft_FCCR_all(:,:)=ch_vec_fft_FCCR;
% ch_vec_fft=abs(ch_vec_fft);
%ch_vec_fft=ch_vec_fft/max(max(ch_vec_fft));

for A_i=1:N_antenna
    ch_vec_fft=ch_vec_fft_FCCR_all(:,:,A_i);
    %  ch_vec_fft=ch_vec_fft_FCCR_all_sum;


    %% Find the maximum point of the 2D vector
    %% The range is confined in [N_start+1,N_max]
    %% The speed is confined in [M_start+1,M_max] (positive) and [MK-M_max+1:MK-M_start] (negative)
    [row_max1,row_index1]=max(abs(ch_vec_fft(N_start+1:N_max,M_start+1:M_max)));

    %% find the maximum value and index of each column (for positive speed)
    [row_max2,row_index2]=max(abs(ch_vec_fft(N_start+1:N_max,MK-M_max+1:MK-M_start)));

    %% find the maximum value and index of each column (for negative speed)
    [total_max1,total_index1]=max(row_max1); %% find the maximum value in the row
    [total_max2,total_index2]=max(row_max2); %% find the maximum value in the row
    if total_max1>=total_max2
        range_index=row_index1(total_index1); %% range index of the maximum point
    else
        range_index=row_index2(total_index2); %% range index of the maximum point
    end

    %% Range and speed estimation by using the maximum point of the 2D-FFT
    %% This gives the best performance (ML estimations)
    %% Speed estimation
    if total_max1>= total_max2%% speed is positive
        %replaced the multiplication and division with a fixed variable
        v_est=(total_index1-1+M_start)/(MK*NNP*samp_period)*150/fc; %% the estimated speed (m/s)
    else %% speed is negative
        v_est=(-M_max-1+total_index2)/(MK*NNP*samp_period)*150/fc; %% the estimated speed (m/s) (verified)
    end
    %   v_error4(k_snr)=v_error4(k_snr)+abs(v_av(2)*10/36-v_est); %% speed error (second target) (m/s)
    A_tmp=ch_vec_fft(N_start+1:N_max,M_start+1:M_max);

    angle_result= angle(A_tmp(range_index,total_index1));

                angle_result_1=real(A_tmp(range_index,total_index1));
            angle_result_2=imag(A_tmp(range_index,total_index1));

    speed_result= v_est;
    %% Range estimation
    %% The calibrated results (calibrated with average range error almost zero for all K1)
    %% Due to filter and oversize IFFT, calibrating is necessary
    peak_points=((range_index-1)/K1)*(10^6*samp_period)+N_start*(10^6*samp_period)/K1; %% estimated round-trip delays (us)
    %    r_error4(k_snr)=r_error4(k_snr)+abs(peak_points-ch_delay(2)*10^6)*300; %% absolute range error (meter)

    %% the estimation for the second largest target
    range_result=peak_points*150;

end
end


function [carrier, WaveImg1, dmrsAntSymbols, dmrsAntIndices, dmrsSymbols_all, dmrsIndices_all, trBlk_tx, trBlkSizes, pdsch_all, pdschIndices, pdschIndices_all, encodeDLSCH] = img_propagation(img_padded_length, img_input, NSlots, carrier, simLocal, pdsch, encodeDLSCH, decodeDLSCHLocal, pdschextra, harqEntity)

% Loop over the entire waveform length
waveform=[];

dmrsSymbols_all=[];
dmrsIndices_all=[];
trBlk_tx=[];
pdsch_all={};
pdschIndices_all=[];


for nslot = 0:NSlots-1
    % Update the carrier slot numbers for new slot
    carrier.NSlot = nslot;

    % Calculate the transport block sizes for the transmission in the slot
    [pdschIndices,pdschIndicesInfo] = nrPDSCHIndices(carrier, pdsch);
    pdschIndices_all=[pdschIndices_all pdschIndices];
    trBlkSizes = nrTBS(pdsch.Modulation,pdsch.NumLayers,numel(pdsch.PRBSet),pdschIndicesInfo.NREPerPRB,pdschextra.TargetCodeRate,pdschextra.XOverhead);

    % HARQ processing
    for cwIdx = 1:pdsch.NumCodewords
        % If new data for current process and codeword then create a new DL-SCH transport block
        if harqEntity.NewData(cwIdx)
            trBlk = randi([0 1], trBlkSizes(cwIdx), 1);
            trBlk(1: img_padded_length/NSlots) = img_input(:, nslot+1);
            trBlk_tx=[trBlk_tx trBlk];
            setTransportBlock(encodeDLSCH, trBlk, cwIdx-1, harqEntity.HARQProcessID);
            % harqEntity_all{nslot+1}=harqEntity;
            % If new data because of previous RV sequence time out then flush decoder soft buffer explicitly
            if harqEntity.SequenceTimeout(cwIdx)
                resetSoftBuffer(decodeDLSCHLocal,cwIdx-1,harqEntity.HARQProcessID);
            end
        end
    end

    % Encode the DL-SCH transport blocks
    codedTrBlocks = encodeDLSCH(pdsch.Modulation, pdsch.NumLayers, pdschIndicesInfo.G, harqEntity.RedundancyVersion, harqEntity.HARQProcessID);

    % Create resource grid for a slot
    pdschGrid = nrResourceGrid(carrier,simLocal.NTxAnts);
    % codedTrBlocks_all=[codedTrBlocks_all codedTrBlocks];

    % PDSCH modulation and precoding
    pdschSymbols = nrPDSCH(carrier,pdsch,codedTrBlocks);
    pdschAntSymbols=pdschSymbols;
    pdschAntIndices=pdschIndices;

    % PDSCH mapping in grid associated with PDSCH transmission period
    pdschGrid(pdschAntIndices) = pdschAntSymbols;

    % PDSCH DM-RS precoding and mapping
    dmrsSymbols = nrPDSCHDMRS(carrier,pdsch);
    dmrsIndices = nrPDSCHDMRSIndices(carrier,pdsch);

    % [dmrsAntSymbols,dmrsAntIndices] = hPRGPrecode(size(pdschGrid),carrier.NStartGrid,dmrsSymbols,dmrsIndices,wtx);
    dmrsAntSymbols=dmrsSymbols;
    dmrsAntIndices=dmrsIndices;
    pdschGrid(dmrsAntIndices) = dmrsAntSymbols;
    dmrsSymbols_all=[dmrsSymbols_all dmrsSymbols];
    dmrsIndices_all=[dmrsIndices_all dmrsIndices];

    % PDSCH PT-RS precoding and mapping
    ptrsSymbols = nrPDSCHPTRS(carrier,pdsch);
    ptrsIndices = nrPDSCHPTRSIndices(carrier,pdsch);
    ptrsAntSymbols=ptrsSymbols;
    ptrsAntIndices=ptrsIndices;
    pdschGrid(ptrsAntIndices) = ptrsAntSymbols;

    % OFDM modulation
    txWaveform = nrOFDMModulate(carrier,pdschGrid);
    waveform=[waveform; txWaveform];
    %ResourceGrids=[ResourceGrids pdschGrid];
    pdsch_all{nslot+1}=pdsch;
end

WaveImg1=resample(waveform, 2000000, 1228800);
save ~/Desktop/rfnoc_test/U1Wave_Image_QPSK50M_200Msps.mat WaveImg1;

disp('Propagating image bitstream...');
system('python3 ~/Desktop/rfnoc_test/image_rfnoc_replay.py');
end

function a = rfnoc_test(~)
disp("Running rfnoc_test.py...");
system('python3 ~/Desktop/rfnoc_test/rfnoc_replay_30k.py');
a = -1;
end

function image_encoder(img_filename)
disp('Encoding image...');
system(['python3 ~/Desktop/matlab/images/img_encoder.py ' img_filename]); % todo: nativize this script in matlab?
end
