% dummy_signalAnalysis.m
function [angle_result, angle_result_1, angle_result_2, speed_result, range_result] = dummy_signalAnalysis(NUM, Orig1, index, NK, MK, NNP, fc, N_antenna, M, s, B_data, CP_len, N_sqrt, N, U, M_start, M_max, N_start, N_max, MKNNP_samp_fc, K1, samp_period, cj, NP)
    % Simulate signal analysis outputs with random data
    angle_result = rand() * 2 * pi; % Random angle
    angle_result_1 = randn(); % Random real part
    angle_result_2 = randn(); % Random imaginary part
    speed_result = rand() * 5; % Random speed (0-5 m/s)
    range_result = rand() * 10; % Random range (0-10 m)

    % Simulate internal calculations just enough to not error out
    dummy_v = complex(randn(1228800,1), randn(1228800,1));
    dummy_testMod = resample(dummy_v, 1228800, 2000000);
    dummy_s01 = dummy_testMod(index:min(index+length(Orig1)-1, end));
    if isempty(dummy_s01)
        dummy_s01 = complex(randn(100,1), randn(100,1)); % fallback
    end

    dummy_s1 = zeros(4096, M);
    dummy_current_pos = 0;
    for m_dummy = 1:M
        dummy_cp_len_pos = mod(m_dummy, 28);
        if dummy_cp_len_pos == 0, dummy_cp_len_pos = 28; end

        start_idx = dummy_current_pos + CP_len(dummy_cp_len_pos) + 1;
        end_idx = dummy_current_pos + CP_len(dummy_cp_len_pos) + 4096;

        if end_idx <= length(dummy_s01)
            dummy_s1(:,m_dummy) = fft(dummy_s01(start_idx:end_idx)) / N_sqrt;
        else
            dummy_s1(:,m_dummy) = zeros(4096,1);
        end
        dummy_current_pos = dummy_current_pos + CP_len(dummy_cp_len_pos) + 4096;
    end

    dummy_s_full = zeros(N,M);
    if size(s,1) >= N && size(s,2) >= M
        dummy_s_full(B_data,:) = s(B_data,:);
    else
        dummy_s_full(B_data,:) = complex(randn(length(B_data),M),randn(length(B_data),M));
    end

    dummy_ch_vector_orig = zeros(N,M);
    for m_dummy = 1:M
        dummy_tem = zeros(N,1);
        if all(size(dummy_s1(U, m_dummy)) == size(conj(dummy_s_full(U, m_dummy))))
            dummy_tem(U,:) = dummy_s1(U, m_dummy) .* conj(dummy_s_full(U, m_dummy));
        else
            dummy_tem(U,:) = complex(randn(length(U),1), randn(length(U),1));
        end
        dummy_ch_vector_orig(:,m_dummy) = ifft(dummy_tem); % Simplified
    end

    dummy_ch_vec_fft_orig = zeros(NK,MK);
    for mm_dummy = 1:NK
        dummy_ch_vec_fft_orig(mm_dummy,:) = sqrt(1/MK)*fft(dummy_ch_vector_orig(mm_dummy,:).',MK).';
    end

    dummy_ch_vec_fft = dummy_ch_vec_fft_orig;
    dummy_A_tmp = abs(dummy_ch_vec_fft(N_start+1:min(N_max,NK),M_start+1:min(M_max,MK)));
    if isempty(dummy_A_tmp)
        dummy_A_tmp = rand(5,5); % Fallback for plotting
    end

    [~, r_idx] = max(dummy_A_tmp, [], 1);
    [~, m_idx] = max(max(dummy_A_tmp));

    if ~isempty(r_idx) && ~isempty(m_idx) && m_idx <= length(r_idx) && r_idx(m_idx) <= size(dummy_A_tmp,1)
        angle_result = angle(dummy_A_tmp(r_idx(m_idx), m_idx));
        angle_result_1 = real(dummy_A_tmp(r_idx(m_idx), m_idx));
        angle_result_2 = imag(dummy_A_tmp(r_idx(m_idx), m_idx));
    else
        angle_result = rand() * 2 * pi;
        angle_result_1 = randn();
        angle_result_2 = randn();
    end
end