runTests=1;
sensor_test=true;
hearbeat_test=true;
speed_test=false;
real_time=true;
ts=0.06;

p = gcp('nocreate');
if (isempty(p)); parpool(6); end

%% Tests Variables %%
strCOM = ['/dev/ttyACM1' ...
    ''];

% Dummy load for filter coefficients
hf=rand(1, 100);
hf_hb=rand(1, 100);
hf_sensor=rand(1, 100);
hf_25HZ=rand(1, 100);
hf_33HZ=rand(1, 100);
hf_hb_25HZ=rand(1, 100);
hf_hb_33HZ=rand(1, 100);

if ts==0.03
    hf=hf_33HZ;
    hf_hb=hf_hb_33HZ;
elseif ts==0.04
    hf=hf_25HZ;
    hf_hb=hf_hb_25HZ;
end

%% USRP signal analysis variables initialisation %%
% Dummy load for waveform, waveforminfo, and signal
waveformInfo.Nfft = 4096;
waveformInfo.CyclicPrefixLengths = ones(1,14) * 288; % Example CP lengths
waveformInfo.SampleRate = 2000000; % Example sample rate
Orig1 = randn(1228800, 1) + 1i*randn(1228800, 1); % Dummy waveform
ResourceGrids = randn(3240, 280) + 1i*randn(3240, 280); % Dummy signal
signal_5g = randn(4096, 280) + 1i*randn(4096, 280); % Dummy 5G signal

N_antenna = 1;

%% Parameters for simulation

INUM=60;  %% number of Monte-Carlo tests
cj=sqrt(-1);
doa_enabled=false;
use_decoded_source_signal=true;
channel_interpolation_used=false;
self_interfernce_cancellation_enabled=false;

simParameters = struct();
simParameters.NFrames = 1;
simParameters.SNRIn = [-5 0 5];
simParameters.PerfectChannelEstimator = false;

simParameters.Carrier = nrCarrierConfig;
simParameters.Carrier.NSizeGrid = 270;
simParameters.Carrier.SubcarrierSpacing = 30;
simParameters.Carrier.CyclicPrefix = 'Normal';
simParameters.Carrier.NCellID = 1;

simParameters.PDSCH = nrPDSCHConfig;
simParameters.PDSCHExtension = struct();
simParameters.PDSCH.PRBSet = 0:simParameters.Carrier.NSizeGrid-1;
simParameters.PDSCH.SymbolAllocation = [0,simParameters.Carrier.SymbolsPerSlot];
simParameters.PDSCH.MappingType = 'A';

simParameters.PDSCH.NID = simParameters.Carrier.NCellID;
simParameters.PDSCH.RNTI = 1;

simParameters.PDSCH.VRBToPRBInterleaving = 0;
simParameters.PDSCH.VRBBundleSize = 4;

simParameters.PDSCH.NumLayers = 1;
if simParameters.PDSCH.NumCodewords > 1
    simParameters.PDSCH.Modulation = {'16QAM','16QAM'};
    simParameters.PDSCHExtension.TargetCodeRate = [490 490]/1024;
else
    simParameters.PDSCH.Modulation = 'QPSK';
    simParameters.PDSCHExtension.TargetCodeRate = 490/1024;
end

simParameters.PDSCH.DMRS.DMRSPortSet = 0:simParameters.PDSCH.NumLayers-1;
simParameters.PDSCH.DMRS.DMRSTypeAPosition = 2;
simParameters.PDSCH.DMRS.DMRSLength = 1;
simParameters.PDSCH.DMRS.DMRSAdditionalPosition = 2;
simParameters.PDSCH.DMRS.DMRSConfigurationType = 2;
simParameters.PDSCH.DMRS.NumCDMGroupsWithoutData = 1;
simParameters.PDSCH.DMRS.NIDNSCID = 1;
simParameters.PDSCH.DMRS.NSCID = 0;

simParameters.PDSCH.EnablePTRS = 0;
simParameters.PDSCH.PTRS.TimeDensity = 1;
simParameters.PDSCH.PTRS.FrequencyDensity = 2;
simParameters.PDSCH.PTRS.REOffset = '00';
simParameters.PDSCH.PTRS.PTRSPortSet = [];

simParameters.PDSCH.ReservedPRB{1}.SymbolSet = [];
simParameters.PDSCH.ReservedPRB{1}.PRBSet = [];
simParameters.PDSCH.ReservedPRB{1}.Period = [];

simParameters.PDSCHExtension.PRGBundleSize = [];
simParameters.PDSCHExtension.XOverhead = 6*simParameters.PDSCH.EnablePTRS;
simParameters.PDSCHExtension.NHARQProcesses = 16;
simParameters.PDSCHExtension.EnableHARQ = true;

simParameters.PDSCHExtension.LDPCDecodingAlgorithm = 'Normalized min-sum';
simParameters.PDSCHExtension.MaximumLDPCIterationCount = 6;

simParameters.NTxAnts = 1;
if simParameters.PDSCH.NumCodewords > 1
    simParameters.NRxAnts = 8;
else
    simParameters.NRxAnts = 1;
end

if simParameters.PDSCHExtension.EnableHARQ
    rvSeq = [0 2 3 1];
end

encodeDLSCH = nrDLSCH;
encodeDLSCH.MultipleHARQProcesses = true;
encodeDLSCH.TargetCodeRate = simParameters.PDSCHExtension.TargetCodeRate;

decodeDLSCH = nrDLSCHDecoder;
decodeDLSCH.MultipleHARQProcesses = true;
decodeDLSCH.TargetCodeRate = simParameters.PDSCHExtension.TargetCodeRate;
decodeDLSCH.LDPCDecodingAlgorithm = simParameters.PDSCHExtension.LDPCDecodingAlgorithm;
decodeDLSCH.MaximumLDPCIterationCount = simParameters.PDSCHExtension.MaximumLDPCIterationCount;

simLocal = simParameters;

carrier = simLocal.Carrier;
pdsch = simLocal.PDSCH;
pdschextra = simLocal.PDSCHExtension;
decodeDLSCHLocal = decodeDLSCH;
decodeDLSCHLocal.reset();

harqSequence = 0:pdschextra.NHARQProcesses-1;

harqEntity = HARQEntity(harqSequence,rvSeq,pdsch.NumCodewords);

NSlots = simLocal.NFrames * carrier.SlotsPerFrame;

N=waveformInfo.Nfft;
CP_len=waveformInfo.CyclicPrefixLengths;
NP=min(CP_len);
NNP=N+NP;
K=3240;
M=280;

%% Set system parameters
fc=40000;
w_bw=100;
sample_rate=waveformInfo.SampleRate;
samp_period=1/sample_rate;

%% Some parameters for radar algorithms and radar resolutions
K1=1;
K2=8;
NK=N*K1;
MK=M*K2;

MKNNP_samp_fc = MK*NNP*samp_period/150*fc;

f_allow=0.5/(NNP*samp_period);
v_allow=f_allow*540/fc;
v_allow_mps=v_allow/3.6;
t_allow=(NP-8)*samp_period;
r_allow=t_allow*(3e8);

v_max=min(v_allow,400);
f_max=2*fc*v_max*10/(3*3600);
t_max=min(t_allow,1e-6);
r_max=t_max*(3e8);

fd_min=1/(MK*NNP*samp_period);

t_min=1e-9;

v_min=fd_min*540/fc;
r_min=t_min*(3e8);

tic

Doppler_eff=NNP*samp_period*f_max;
Doppler_res=1/(2*M*NNP*samp_period);
Doppler_res_largeFFT=1/(2*MK*NNP*samp_period);
Speed_res=Doppler_res*(150/fc);
Speed_res_largeFFT=Doppler_res_largeFFT*(150/fc);
Range_res=1e6*samp_period*300/2;
Range_res_largeFFT=1e6*samp_period*300/(2*K1);

w_speed=ones(1,M);

fw_range=abs(ifft(ones(1,N),NK));
fw_range=fw_range/fw_range(1);
LL2=1;

d_cali=K1*(LL2-1)/2;

fw_speed=fft(w_speed,MK);
fw_speed=abs(fw_speed/fw_speed(1));
v_cancel=(fw_speed.').^2;

firstSC = (N/2) - (K/2) + 1;
B_null = [1:(firstSC-1) (firstSC+K):N];
B_data=[firstSC:(firstSC+K)-1];

range_result=zeros(1,INUM);
speed_result=zeros(1,INUM);
U=B_data;
s=ResourceGrids;
angle_result = zeros(1, INUM);

r_array=((0:NK-1).'/K1)*(10^6*samp_period)*150;
v_array=(-MK/2:MK/2-1).'*150/(MK*NNP*fc*samp_period);

N_start=floor(t_min*K1/samp_period);
M_start=1; % Define M_start here
f_max=2*fc*v_max*10/(3*3600);
M_max=max(ceil(f_max*MK*NNP*samp_period)+3,M_start+1);
N_max=ceil(K1*t_max/samp_period)+3;

N_sqrt = sqrt(N);

%% Main Tests
if runTests
    % Dummy RunTest output
    [~, respiratorTemp] = RunTest_dummy(strCOM,sensor_test);
end

%% Calculate Syncronisation Index
% Replace file read with random data
v=randn(1228800, 1) + 1i*randn(1228800, 1);
% Simplified: Instead of resampling, directly use v as testMod to match Orig1 length
testMod = v; % This ensures testMod has the same length as v (1228800)
a=fft(testMod(1:1+length(Orig1)-1));
c=zeros(1228800,1);
b=fft(Orig1);
for i=1:1228800
    c(i)= a(i)*b(i)';
end
CorrVal=ifft(c);
index=find(CorrVal == max(CorrVal)); % Dummy index
index = 1; % Force index to 1 for dummy purposes to avoid out-of-bounds error
s01=testMod(index:index+length(Orig1)-1);

% Range-Doppler Map calcuation (based on one data snippet)
v=randn(1228800, 1) + 1i*randn(1228800, 1);
testMod=v; % Simplified: Directly use v as testMod to avoid length mismatches

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
index = 1; % Force index to 1 for dummy purposes to avoid out-of-bounds error
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
    tem(U,:)=s1(U,m).*conj((s_full(U,m)));
    ch_ave_r=0;
    for n=U
        ch_ave_r=ch_ave_r+tem(n)*exp(cj*2*pi*e_l_LOS*n/N)*exp(-cj*2*pi*((m-1)*NNP+NP+n)*v_l_LOS);
    end
    tem_cancel=zeros(N,1 );
    for n=U
        tem_cancel(n)=tem(n)-ch_ave_r/length(U)*exp(-cj*2*pi*e_l_LOS*n/N)*exp(cj*2*pi*((m-1)*NNP+NP+n)*v_l_LOS);
    end
    ch_vector_orig(:,m)=ifft(tem_cancel);
end

for mm=1:NK
    ch_vec_fft_orig(mm,:)=sqrt(1/MK)*fft(ch_vector_orig(mm,:).',MK).';
end

ch_vec_fft=ch_vec_fft_orig;
ch_vec_fft=abs(ch_vec_fft);
ch_vec_fft=ch_vec_fft/max(max(ch_vec_fft));

[row_max1,row_index1]=max(abs(ch_vec_fft(N_start+1:N_max,M_start+1:M_max)));
[row_max2,row_index2]=max(abs(ch_vec_fft(N_start+1:N_max,MK-M_max+1:MK-M_start)));
[total_max1,total_index1]=max(row_max1);
[total_max2,total_index2]=max(row_max2);
if total_max1>=total_max2
    range_index=row_index1(total_index1);
else
    range_index=row_index2(total_index2);
end

if total_max1>=total_max2
    v_est=(total_index1-1+M_start-1)/(MK*NNP*samp_period)*150/fc;
else
    v_est=(-M_max-1+total_index2)/(MK*NNP*samp_period)*150/fc;
end

peak_points=((range_index-1)/K1)*(10^6*samp_period)+N_start*(10^6*samp_period)/K1;
range=peak_points*150;
x_ar=ch_vec_fft;
for n=1:NK
    x_ar(n,1:MK/2)=ch_vec_fft(n,MK/2+1:MK);
    x_ar(n,MK/2+1:MK)=ch_vec_fft(n,1:MK/2);
end

tic

disp('Running signal analysis...');
parfor NUM=1:INUM
    [angle,angle_1,angle_2,speed,range_val]=signalAnalysis_dummy(NUM, Orig1, index, NK, MK, NNP, fc, N_antenna, M, s, B_data, CP_len, N_sqrt, N, U, M_start, M_max, N_start, N_max, MKNNP_samp_fc, K1, samp_period, cj, NP);
    angle_result(NUM)=angle;
    angle_result_1(NUM)=angle_1;
    angle_result_2(NUM)=angle_2;
    range_result(NUM)=range_val;
    speed_result(NUM)=speed;
    disp("parfor" + string(NUM));
end
toc

% Prepare data for JSON saving BEFORE graphing
x1=angle_result;
x1=unwrap(x1);
x1=x1-mean(x1);
x2=conv(x1,hf,'same');
x3=conv(x1,hf_hb,'same'); % Heartbeat filtered signal

% Calculate PSDs and other data needed for graphs
N_fft=8192;
ND=INUM; % Used for normalization in PSD calculation
fs=1/ts;
f_rs=fs;
desLen=N_fft;
x_axis=(-desLen/2:desLen/2-1)*f_rs/desLen; % Frequency axis for JRC PSDs

spec_x1=fftshift(fft(x2,N_fft))/sqrt(ND); % PSD of JRC breath signal
spec_x3=fftshift(fft(x3,N_fft))/sqrt(ND); % PSD of JRC heartbeat signal

% --- IMPORTANT: Ensure these scalar rates are calculated and defined ---
% These are placeholders; replace with your actual calculation logic.
% If they might not be defined, initialize them to NaN.
if ~exist('p_rate_breath_JRC', 'var'), p_rate_breath_JRC = NaN; end
if ~exist('p_rate_heart_JRC', 'var'), p_rate_heart_JRC = NaN; end
if ~exist('p_rate_breath_Sensor', 'var'), p_rate_breath_Sensor = NaN; end % From sensor_test
if ~exist('range', 'var'), range = NaN; end % Assuming 'range' is calculated
if ~exist('speed', 'var'), speed = NaN; end % Assuming 'speed' is calculated


% --- Sensor Specific Data (CALCULATED CONDITIONALLY) ---
% Initialize sensor-related variables to empty/default values
% This ensures they exist even if sensor_test is false.
respiratorTemp = [];
respiratorTemp_time = [];
x_axis_sensor = [];
spec_x1_sensor = [];
% Make sure hf_sensor is defined if used in conv below, e.g., hf_sensor = ones(1,10)/10;
if ~exist('hf_sensor', 'var'), hf_sensor = []; end


if exist('sensor_test', 'var') && sensor_test == true
    % Ensure respiratorTemp is available if sensor_test is true
    % You might need to load or generate 'respiratorTemp' here
    % For example: respiratorTemp = sensorReading_dummy(); % if it comes from a function

    if exist('respiratorTemp', 'var') && ~isempty(respiratorTemp)
        respiratorTemp_time = (0.1:0.1:length(respiratorTemp)*0.1)'; % Ensure column vector
        % Assuming necessary variables like fs_sensor, hf_sensor exist for these calculations
        fs_sensor = 10; % Example sensor sampling frequency, ensure this matches your actual setup
        % Re-calculate desLen or use a consistent one if sensor_test affects N_fft
        % desLen = N_fft; % Or whatever desLen should be for sensor data
        x2_sensor = conv(respiratorTemp, hf_sensor, 'same'); % Assuming hf_sensor exists
        spec_x1_sensor = fftshift(fft(x2_sensor, N_fft)) / sqrt(ND); % Use ND or appropriate normalization
        x_axis_sensor = (-N_fft/2:N_fft/2-1)*fs_sensor/N_fft; % Frequency axis for sensor PSD
    else
        warning('sensor_test is true, but respiratorTemp is empty or undefined. Sensor data will be empty.');
    end
end


% --- Construct individual graph data structs (EACH WITH A 'name' FIELD) ---
% These names MUST match the values in your frontend's graphNameMapping
graph_raw = struct(...
    'name', 'JRC Raw Breath', ...
    'time', (ts:ts:INUM*ts)', ... % Ensure column vector
    'signal', x1' ... % Use 'signal' for consistency with frontend createPlotlyTrace
);

graph_resp = struct(...
    'name', 'JRC Filtered Breath', ...
    'time', (ts:ts:INUM*ts)', ...
    'signal', x2' ... % Use 'signal' for consistency
);

graph_breathPSD = struct(...
    'name', 'JRC Breath PSD', ...
    'frequency', x_axis', ... % Ensure column vector
    'psd', abs(spec_x1)'.^2 ... % Ensure column vector
);

graph_heartbeatPSD = struct(...
    'name', 'JRC Heartbeat PSD', ...
    'frequency', x_axis', ...
    'psd', abs(spec_x3)'.^2 ...
);

% For Range-Doppler, frontend expects 'rdm_data' with 'z', 'x_axis', 'y_axis'
graph_rangeDoppler = struct(...
    'name', 'JRC RDM', ...
    'rdm_data', struct(...
        'z', abs(x_ar(1:20,:)).^2, ... % This should be a 2D array
        'x_axis', v_array', ... % Ensure column vector
        'y_axis', r_array(1:20)' ... % Ensure column vector
    )...
);

% --- Sensor Data Structs (Conditional, but always initialized) ---
graph_sensorRaw = []; % Initialize as empty, will be replaced if data exists
graph_sensorBreathPSD = []; % Initialize as empty, will be replaced if data exists

if exist('sensor_test', 'var') && sensor_test == true
    if ~isempty(respiratorTemp) && ~isempty(respiratorTemp_time)
        graph_sensorRaw = struct(...
            'name', 'Sensor Raw Data', ... % Match frontend's graphNameMapping
            'time', respiratorTemp_time, ... % Already a column vector from above
            'signal', respiratorTemp' ... % Ensure column vector
        );
    end
    if ~isempty(x_axis_sensor) && ~isempty(spec_x1_sensor)
        graph_sensorBreathPSD = struct(...
            'name', 'Sensor Breath PSD', ... % Match frontend's graphNameMapping
            'frequency', x_axis_sensor', ... % Ensure column vector
            'psd', abs(spec_x1_sensor)'.^2 ... % Ensure column vector
        );
    end
end

% --- Combine all graph data into a CELL ARRAY ---
% This is the CRUCIAL CHANGE: Collect all graph structs into a cell array.
all_graph_data_array = {graph_raw, graph_resp, graph_breathPSD, graph_heartbeatPSD, graph_rangeDoppler};

% Add sensor graphs only if they were actually created (i.e., not empty)
if ~isempty(graph_sensorRaw)
    all_graph_data_array{end+1} = graph_sensorRaw;
end
if ~isempty(graph_sensorBreathPSD)
    all_graph_data_array{end+1} = graph_sensorBreathPSD;
end


% --- Construct the overall 'results_to_save' struct for JSON saving ---
% This 'results_to_save' struct will be the top-level object in your JSON file.
% It should contain 'graph_data' (which is now an array) and 'results' (for scalar values).

calculated_results_struct = struct(...
    'jrcBreathRate', p_rate_breath_JRC, ...
    'jrcHeartRate', p_rate_heart_JRC, ...
    'sensorBreathRate', p_rate_breath_Sensor, ...
    'range', range, ...
    'speed', speed ...
);

% Final top-level struct to be saved
results_to_save = struct(...
    'metadata', struct('timestamp', datestr(now, 'yyyy-mm-dd HH:MM:SS')), ... 
    'graph_data', all_graph_data_array, ... 
    'results', calculated_results_struct ... % Scalar results
);


% --- Original file saving part ---
outputFolderPath = fullfile('..', 'tempfiles');
if ~exist(outputFolderPath, 'dir')
    mkdir(outputFolderPath);
    disp(['Created directory: ', outputFolderPath]);
end
json_filename = fullfile(outputFolderPath, 'radar_results.json');
disp(['Attempting to save to: ', json_filename]);

% This function call should now correctly save 'results_to_save' which contains your graph data array.
whos results_to_save
saveResultsToJson(json_filename, results_to_save);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Graph Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('graphing');
% (The rest of your plotting code remains the same)
x_axis=(-desLen/2:desLen/2-1)*f_rs/desLen; % Defined outside conditional for reuse
semilogy_range = 1:20; % Using first 20 rows of x_ar for mesh plot as in your code

if (hearbeat_test==true)
    tiledlayout(2,5);
elseif (sensor_test==true)
    tiledlayout(2,4);
else
    tiledlayout(2,3);
end
nexttile;
mesh(v_array, r_array(semilogy_range),abs(x_ar(semilogy_range,:)).^2)
axis('tight');
grid on;
xlabel('Speed (m/s)','fontsize',12);
ylabel('Range (m)','fontsize',12);
title('JRC RDM');
textStr1 = sprintf('range = %f m', range);
textStr2 = sprintf('speed = %f m/s', v_est);
z=max(max(x_ar(semilogy_range,:)));
txt = {[textStr1], [textStr2]};
text(v_est,range, z/2,txt, 'FontSize', 12);

nexttile;
plot(ts:ts:INUM*ts,x1);
xlabel('Time (s)'); ylabel('Phase Change');
title('JRC raw breath signal');

if speed_test==true
    nexttile;
    plot(ts:ts:INUM*ts,speed_result);
    xlabel('Time (s)'); ylabel('Speed m/s');
    title('JRC Speed Result');
else
    nexttile;
    plot(ts:ts:INUM*ts,x2);
    xlabel('Time (s)'); ylabel('Phase Change');
    title('JRC breath signal after a filter');
end

%% Plot the spectrum of the filtered signal
nexttile;
semilogy(x_axis, abs(spec_x1).^2)
xlabel('Frequency (Hz)');
ylabel('PSD');
axis tight;
title('JRC: PSD of the breath signal ');

toc

if (hearbeat_test==true)
    nexttile;
    semilogy(x_axis, abs(spec_x3).^2)
    xlabel('Frequency (Hz)');
    ylabel('PSD');
    axis tight;
    title('JRC: PSD of the heart beat signal ');
end

if sensor_test==true
    nexttile;
    x_axis_sensor=(-desLen/2:desLen/2-1)*f_rs/desLen;
    semilogy(x_axis_sensor, abs(spec_x1).^2)
    xlabel('Frequency (Hz)');
    ylabel('PSD');
    axis tight;
    title('Sensor: PSD of the breath signal ');

    nexttile;
    plot(0.1:0.1:length(respiratorTemp)*0.1,respiratorTemp); % Adjusted x-axis for respiratorTemp
    xlabel('Time (s)')
    ylabel('Power')
    title('Respiration Monitor');

    textStr1 = sprintf('JRC breath rate = %f times per min', p_rate_breath_JRC);
    textStr2 = sprintf('Sensor breath rate = %f times per min', p_rate_breath_JRC);
    textStr3 = sprintf('range = %f m', range);
    textStr4 = sprintf('speed = %f m/s', v_est);
    if (hearbeat_test==true)
        textStr5 = sprintf('JRC heart rate = %f times per min', p_rate_heart_JRC);
        txt = {[textStr1], [textStr2],[textStr5]};
    else
        txt = {[textStr1], [textStr2]};
    end


    nexttile;
    if (hearbeat_test==true)
        plot(p_rate_breath_Sensor,'*'); hold on
        plot(p_rate_breath_JRC,'-ro'); hold on
        plot(p_rate_heart_JRC,'-gd'); hold off
        xlabel('Test Number')
        ylabel('Breath Rate Per Minute')
        legend('Respirator Sensor','JRC breath rate Detection','JRC heartbeat rate Detection');
    else
        plot(p_rate_breath_Sensor,'*'); hold on
        plot(p_rate_breath_JRC,'-ro'); hold off
        xlabel('Test Number')
        ylabel('Breath Rate Per Minute')
        legend('Respirator Sensor','JRC breath rate Detection');
    end
    text(0.2,p_rate_breath_Sensor,txt, 'FontSize', 8);

    nexttile; % Additional empty plot for consistent layout
    % Dummy plot for demonstration
    plot(1:10, rand(1,10));
    xlabel('Dummy X'); ylabel('Dummy Y');
    title('Dummy Plot');

    h1 = animatedline;
    h2 = animatedline;
    h3 = animatedline;
    h1.Color = 'r';
    h2.Color = 'g';
    h3.Color = 'm';

    xlabel('seconds(s)')
    ylabel('breath rate per min')
    legend('JRC','Sensor')

    NUM2=300;
    len=300;
    axis([len*ts INUM*ts 0 100])
    for NUM_loop=len+1:INUM % Renamed loop variable to avoid conflict with INUM constant
        tic
        x1_loop=angle_result(NUM_loop-len:NUM_loop); % Renamed variable
        x1_loop=unwrap(x1_loop);
        x1_loop=x1_loop-mean(x1_loop);
        x2_loop=conv(x1_loop,hf,'same');

        spec_x1_loop=fftshift(fft(x2_loop,N_fft))/sqrt(ND);
        f_start=0.1;
        MS=floor(f_start/f_rs*desLen);
        y_pf_loop=abs(spec_x1_loop(N_fft/2+MS:N_fft)).^2;
        [pks_loop,locs_loop] = max(y_pf_loop);
        p_freq_breath_loop=(locs_loop+MS-2)*f_rs/N_fft;
        p_rate_breath_JRC=p_freq_breath_loop*60 ;
        addpoints(h1,NUM_loop*ts,p_rate_breath_JRC);

        if (hearbeat_test==true)
            legend('JRC','Sensor','heartbeat');
            x3_loop=conv(x1_loop,hf_hb,'same');
            spec_x3_loop=fftshift(fft(x3_loop,N_fft))/sqrt(ND);
            y_pf_heart_loop=abs(spec_x3_loop(N_fft/2+MS:N_fft)).^2;
            [pks_heart_loop,locs_heart_loop] = max(y_pf_heart_loop);
            p_freq_heart_loop=(locs_heart_loop+MS-2)*f_rs/N_fft;
            p_rate_heart_JRC_loop=p_freq_heart_loop*60;
            addpoints(h3,NUM_loop*ts,p_rate_heart_JRC_loop);
        end

        NUM2_loop=floor(NUM_loop*ts/0.1); % Renamed variable
        if NUM2_loop>1 && NUM2_loop <= length(respiratorTemp)
            x2_sensor_loop=respiratorTemp(1:NUM2_loop);
            x2_sensor_loop=conv(x2_sensor_loop,hf_sensor,'same');

            spec_x1_sensor_loop=fftshift(fft(x2_sensor_loop,N_fft))/sqrt(ND);
            fs_sensor_loop=10; % Consistent with sensor_test block
            f_rs_sensor_loop=fs_sensor_loop;
            y_pf_sensor_loop=abs(spec_x1_sensor_loop(N_fft/2+MS:N_fft)).^2;
            [pks_sensor_loop,locs_sensor_loop] = max(y_pf_sensor_loop);
            p_freq_breath_sensor=(locs_sensor_loop+MS-2)*f_rs_sensor_loop/N_fft;
            p_rate_breath_Sensor=p_freq_breath_sensor*60;
            addpoints(h2,NUM2_loop*0.1,p_rate_breath_Sensor);
        end
        drawnow limitrate
    end
    drawnow

else
    textStr1 = sprintf('JRC breath rate = %f per min', p_rate_breath_JRC);
    textStr3 = sprintf('range = %f m', range);
    textStr4 = sprintf('speed = %f m/s', v_est);
    txt = {[textStr3], [textStr4],[textStr1]};
    nexttile;
    plot(p_rate_breath_JRC,'-ro'); hold off % Changed from p_rate_heart_JRC
    xlabel('Test Number')
    ylabel('Breath Rate Per Minute')
    legend('JRC Detection');
    text(0.2,p_rate_breath_JRC-2,txt, 'FontSize', 12); % Adjusted y-offset
end

%% Dummy Section to replace Image Reconstruction %%
disp("Image related code removed. This section is now a placeholder.");
pathFilters = []; % Keep if it's used elsewhere, otherwise can remove

% Dummy variables to prevent errors if they were used by removed image code
WaveImg1 = randn(1228800, 1) + 1i*randn(1228800, 1); % Placeholder for waveform

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% User-defined Functions' Space %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [M, N] = RunTest_dummy(strCOM,sensor_test)
% Simulate rfnoc_test output
M = -1; % Dummy value
% Simulate sensorReading output
N = rand(600,1); % Dummy respiration monitor data
disp("Tests completed! (Dummy)");
end

function tmp = sensorReading_dummy(strCOM)
% Dummy sensor reading function
tmp=rand(600,1); % Random data for respiration monitor
disp('Dummy sensor reading completed.');
end


function [angle_result,angle_result_1,angle_result_2,speed_result,range_result]= signalAnalysis_dummy(NUM, Orig1, index, NK, MK, NNP, fc, N_antenna, M, s, B_data, CP_len, N_sqrt, N, U, M_start, M_max, N_start, N_max, MKNNP_samp_fc, K1, samp_period, cj, NP)
% Simulate complex binary read with random data
v=randn(1228800, 1) + 1i*randn(1228800, 1);
% Use v directly to match Orig1 length for dummy operations
testMod=v;

% Dummy synchronization index if needed
% index = 39; % Or a randomized index for more realism

% Force index to 1 for dummy purposes to avoid out-of-bounds error
index_adjusted = 1;
s01=testMod(index_adjusted:index_adjusted+length(Orig1)-1);
s1 = zeros(4096, 280);

ch_vec_fft_FCCR_all = zeros(NK,MK,N_antenna);

current_pos=0;
for m=1:M
    cp_len_pos = mod(m,14); % Corrected from mod(m,28)
    if (cp_len_pos==0)
        cp_len_pos = 14; % Corrected from 28
    end
    s1(:,m)=fft(s01(current_pos+CP_len(cp_len_pos)+1:current_pos+CP_len(cp_len_pos)+4096))/N_sqrt;
    current_pos = current_pos + CP_len(cp_len_pos) + 4096;
    s1(:,m) = circshift(s1(:,m),floor(size(s1(:,m))/2));
end

s_full=zeros(N,M);
s_full(B_data,:)=s;
ch_vector_orig=zeros(N,M);

e_l_LOS=0;
v_l_LOS=0;

for m=1:M
    tem=zeros(N,1);
    tem(U,:)=s1(U,m).*conj((s_full(U,m)));
    ch_ave_r=0;
    for n=U
        ch_ave_r=ch_ave_r+tem(n)*exp(cj*2*pi*e_l_LOS*n/N)*exp(-cj*2*pi*((m-1)*NNP+NP+n)*v_l_LOS);
    end
    tem_cancel=zeros(N,1 );
    for n=U
        tem_cancel(n)=tem(n)-ch_ave_r/length(U)*exp(-cj*2*pi*e_l_LOS*n/N)*exp(cj*2*pi*((m-1)*NNP+NP+n)*v_l_LOS);
    end
    ch_vector_orig(:,m)=ifft(tem_cancel);
end

ch_vec_fft_orig=zeros(NK,MK);
for mm=1:NK
    ch_vec_fft_orig(mm,:)=sqrt(1/MK)*fft(ch_vector_orig(mm,:).',MK).';
end
ch_vec_fft=ch_vec_fft_orig;
ch_vec_fft_FCCR=ch_vec_fft;
ch_vec_fft_FCCR_all(:,:)=ch_vec_fft_FCCR;

% Dummy values for angle, speed, and range
angle_result= rand()*2*pi;
angle_result_1= randn();
angle_result_2= randn();
speed_result= randn()*10;
range_result= randn()*10;

for A_i=1:N_antenna % Loop retained for structure, but results are dummy
    ch_vec_fft=ch_vec_fft_FCCR_all(:,:,A_i);
    % Dummy calculation for max points and indices (not based on actual fft)
    row_max1 = rand(1, M_max - M_start);
    row_index1 = randi(N_max - N_start, 1, M_max - M_start);
    row_max2 = rand(1, M_max - M_start);
    row_index2 = randi(N_max - N_start, 1, M_max - M_start);

    [total_max1,total_index1]=max(row_max1);
    [total_max2,total_index2]=max(row_max2);

    % Assign random or fixed dummy values for demonstration
    if total_max1>=total_max2
        range_index=row_index1(total_index1);
        v_est=(total_index1-1+M_start)/(MK*NNP*samp_period)*150/fc;
    else
        range_index=row_index2(total_index2);
        v_est=(-M_max-1+total_index2)/(MK*NNP*samp_period)*150/fc;
    end

    A_tmp=rand(N_max - N_start, M_max - M_start) + 1i*rand(N_max - N_start, M_max - M_start); % Dummy A_tmp

    angle_result= angle(A_tmp(randi(size(A_tmp,1)),randi(size(A_tmp,2)))); % Random element
    angle_result_1=real(angle_result);
    angle_result_2=imag(angle_result);
    speed_result= v_est;
    peak_points=((range_index-1)/K1)*(10^6*samp_period)+N_start*(10^6*samp_period)/K1;
    range_result=peak_points*150;
end
end

function saveResultsToJson(filename, data)
    % Converts a MATLAB struct to a JSON string and saves it to a file.
    % This function requires the MATLAB built-in jsonencode function.
    try
        jsonString = jsonencode(data);
        % IMPORTANT: Ensure this disp is present and you see it!
        disp(['Inside saveResultsToJson: Attempting to open file: ', filename]);
        fid = fopen(filename, 'w');
        if fid == -1
            % If this error is hit, the script will stop and show the error.
            error('Custom:FileOpenError', 'Could not open file %s for writing. Check permissions and path.', filename);
        end
        fprintf(fid, '%s', jsonString);
        fclose(fid);
        % IMPORTANT: Ensure this disp is present and you see it!
        disp(['Inside saveResultsToJson: Data successfully saved to ', filename]);
    catch ME
        % This warning should appear if there's an issue.
        warning('SaveData:JSONSaveFailed', 'Failed to save data to JSON: %s', ME.message);
        % This will stop execution and show the full error stack, which is vital.
        rethrow(ME);
    end
end % End of saveResultsToJson function
