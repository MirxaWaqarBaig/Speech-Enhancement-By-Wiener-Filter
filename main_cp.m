clear all 
clc
%% Read the data from the excel file. 
%From third Column according to instruction
% I've changed the readmatrix to xlsread
% Only used the "C" column data
vowel_u = xlsread('vowel_u_all_trial_subs.xls','C:C'); 
vowel_i = xlsread('vowel_i_all_trial_subs.xls','C:C');
vowel_e = xlsread('vowel_e_all_trial_subs.xls','C:C');
vowel_ae = xlsread('vowel_ae_all_trial_subs.xls','C:C');
vowel_a = xlsread('vowel_a_all_trial_subs .xls','C:C');
%% Plotting the Original Signals from the excel file column three
figure(1)
subplot(2,3,1)
plot(vowel_u) %this is your signal
title('Vowel U Signal')
subplot(2,3,2)
plot(vowel_i) %this is your signal
title('Vowel I Signal')
subplot(2,3,3)
plot(vowel_e) %this is your signal
title('Vowel E Signal')
subplot(2,3,4)
plot(vowel_ae) %this is your signal
title('Vowel AE Signal')
subplot(2,3,5)
plot(vowel_a) %this is your signal
title('Vowel A Signal')
%% Saving data in WAV format. I re-wrote this portion of saving Audio files
%% Vowel U
%%%% Uncomment it and run If you want to save some other data 
%%%% I already saved it so don't need to run it again
% xlim([0 200704])
% sound(vowel_u,40000) %play the signal
% audiowrite('Vowel_U_Column3.wav',vowel_u,40000,'BitsPerSample',16);
%% Vowel I
% xlim([0 200704])
% % sound(vowel_i,40000) %play the signal
% audiowrite('Vowel_I_Column3.wav',vowel_i,40000,'BitsPerSample',16);
 %% Vowel E
% xlim([0 200704])
% % sound(vowel_e,40000) %play the signal
% audiowrite('Vowel_E_Column3.wav',vowel_e,40000,'BitsPerSample',16);
%% Vowel AE
% xlim([0 200704])
% % sound(vowel_ae,40000) %play the signal
% audiowrite('Vowel_AE_Column3.wav',vowel_ae,40000,'BitsPerSample',16);
%% Vowel A
% xlim([0 200704])
% % sound(vowel_a,40000) %play the signal
% audiowrite('Vowel_A_Column3.wav',vowel_a,40000,'BitsPerSample',16);
%% Plotting the time domain graphs of Vowels present in the third column
%% Vowel Graph in Time Domain
%%%%%%%% Vowel U
[vowel_u_aud,fs_u] = audioread('Vowel_U_Column3.wav');
t_u = linspace(0,length(vowel_u_aud)/fs_u,length(vowel_u_aud));
figure(2)
subplot(2,3,1)
plot(t_u,vowel_u_aud) %Time domain plot of Signal
title('Time Domain graph of Vowel U Signal')
xlabel('Time');
ylabel('Amplitude');
%%%%%%%%% Vowel I
[vowel_i_aud,fs_i] = audioread('Vowel_I_Column3.wav');
t_i= linspace(0,length(vowel_i_aud)/fs_i,length(vowel_i_aud));
% figure('WindowState','maximized')
subplot(2,3,2)
plot(t_i,vowel_i_aud) %Time domain plot of Signal
title('Time Domain graph of Vowel I Signal')
xlabel('Time');
ylabel('Amplitude');
%%%%%%%%% Vowel E
[vowel_e_aud,fs_e] = audioread('Vowel_E_Column3.wav');
t_e= linspace(0,length(vowel_e_aud)/fs_e,length(vowel_e_aud));
% figure('WindowState','maximized')
subplot(2,3,3)
plot(t_e,vowel_e_aud) %Time domain plot of Signal
title('Time Domain graph of Vowel E Signal')
xlabel('Time');
ylabel('Amplitude');
%%%%%%%%% Vowel AE
[vowel_ae_aud,fs_ae] = audioread('Vowel_AE_Column3.wav');
t_ae= linspace(0,length(vowel_ae_aud)/fs_ae,length(vowel_ae_aud));
% figure('WindowState','maximized')
subplot(2,3,4)
plot(t_ae,vowel_ae_aud) %Time domain plot of Signal
title('Time Domain graph of Vowel AE Signal')
xlabel('Time');
ylabel('Amplitude');
%%%%%%%%% Vowel A
[vowel_a_aud,fs_a] = audioread('Vowel_A_Column3.wav');
t_a= linspace(0,length(vowel_a_aud)/fs_a,length(vowel_a_aud));
% figure('WindowState','maximized')
subplot(2,3,5)
plot(t_a,vowel_a_aud) %Time domain plot of Signal
title('Time Domain graph of Vowel A Signal')
xlabel('Time');
ylabel('Amplitude');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Computing SNR values
snr_u = snr(vowel_u_aud);
snr_i = snr(vowel_i_aud);
snr_e = snr(vowel_e_aud);
snr_ae = snr(vowel_e_aud);
snr_a = snr(vowel_a_aud);

%% Plotting the Frequency domain graphs of Vowels present in the third column
%%%%%%%%%% Vowel U Frequency Domain
nfft = 1024; % Samples taken
f_u = linspace(0,fs_u,nfft);
Y_U = abs(fft(vowel_u_aud,nfft));
figure(3);
subplot(2,3,1)
plot(f_u(1:nfft/2),Y_U(1:nfft/2));
title('Freuqncy Domain graph of Vowel U Signal')
xlabel('Frequency');
ylabel('Abs');
%%%%%%%%%%%% Vowel I Frequency Domain
f_i = linspace(0,fs_i,nfft);
Y_I = abs(fft(vowel_i_aud,nfft));
% figure();
subplot(2,3,2)
plot(f_u(1:nfft/2),Y_I(1:nfft/2));
title('Freuqncy Domain graph of Vowel I Signal')
xlabel('Frequency');
ylabel('Abs');
%%%%%%%%%%%% Vowel E Frequency Domain
f_e = linspace(0,fs_e,nfft);
Y_E = abs(fft(vowel_e_aud,nfft));
% figure();
subplot(2,3,3)
plot(f_u(1:nfft/2),Y_E(1:nfft/2));
title('Freuqncy Domain graph of Vowel E Signal')
xlabel('Frequency');
ylabel('Abs');
%%%%%%%%%%%% Vowel AE Frequency Domain
f_ae = linspace(0,fs_ae,nfft);
Y_AE = abs(fft(vowel_ae_aud,nfft));
% figure();
subplot(2,3,4)
plot(f_u(1:nfft/2),Y_AE(1:nfft/2));
title('Freuqncy Domain graph of Vowel AE Signal')
xlabel('Frequency');
ylabel('Abs');
%%%%%%%%%%%% Vowel A Frequency Domain
f_a = linspace(0,fs_a,nfft);
Y_A = abs(fft(vowel_a_aud,nfft));
% figure();
subplot(2,3,5)
plot(f_u(1:nfft/2),Y_A(1:nfft/2));
title('Freuqncy Domain graph of Vowel A Signal')
xlabel('Frequency');
ylabel('Abs');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Adding Noise to the Signals 
%%%% Vowel U  
vowel_u_noisy = vowel_u_aud+(randn(size(vowel_u_aud)));
audiowrite("vowel_U_Noise_c_2.wav",vowel_u_noisy,40000,'BitsPerSample',16)
% sound(vowel_u_noisy);
figure(4)
subplot(2,3,1)
plot(vowel_u_noisy);
title('Noise Added Vowel U Signal')
%%%% Vowel I
vowel_i_noisy = vowel_i_aud+(randn(size(vowel_i_aud)));
audiowrite("vowel_I_Noise_c_2.wav",vowel_i_noisy,40000,'BitsPerSample',16)
% sound(vowel_i_noisy);
subplot(2,3,2)
plot(vowel_i_noisy);
title('Noise Added Vowel I Signal')
%%%%% Vowel E
vowel_e_noisy = vowel_e_aud+(randn(size(vowel_e_aud)));
audiowrite("vowel_E_Noise_c_2.wav",vowel_e_noisy,40000,'BitsPerSample',16)
% sound(vowel_e_noisy);
subplot(2,3,3)
plot(vowel_e_noisy);
title('Noise Added Vowel E Signal')
%%%%% Vowel AE
vowel_ae_noisy = vowel_ae_aud+(randn(size(vowel_ae_aud)));
audiowrite("vowel_AE_Noise_c_2.wav",vowel_ae_noisy,40000,'BitsPerSample',16)
% sound(vowel_i_noisy);
subplot(2,3,4)
plot(vowel_ae_noisy);
title('Noise Added Vowel AE Signal')
%%%%% Vowel A
vowel_a_noisy = vowel_a_aud+(randn(size(vowel_a_aud)));
audiowrite("vowel_A_Noise_c_2.wav",vowel_a_noisy,40000,'BitsPerSample',16)
% sound(vowel_i_noisy);
subplot(2,3,5)
plot(vowel_i_noisy);
title('Noise Added Vowel A Signal')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SNR of Noise Signal 
snr_u_noise = snr(vowel_u_noisy);
snr_i_noise = snr(vowel_i_noisy);
snr_e_noise = snr(vowel_e_noisy);
snr_ae_noise = snr(vowel_e_noisy);
snr_a_noise = snr(vowel_a_noisy);
%% Now implementing the de-noising method of Wiener and ploting graphs
%% voice enhacement of the vowels with wiener
%%%% Vowel U
u0 = vowel_u_noisy;
[u1 , u2] = WienerNoiseReduction(vowel_u_noisy,40000,1);
%writing to a file
%audiowrite("vowel_u_TSNR.wav",u1,40000,'BitsPerSample',16)
%audiowrite("vowel_u_HRNR.wav",u2,40000,'BitsPerSample',16)
%%
%%%% Vowel I
I0 = vowel_i_noisy;
[I1 , I2] = WienerNoiseReduction(vowel_i_noisy,40000,1);
%writing to a file
%audiowrite("vowel_i_TSNR.wav",I1,40000,'BitsPerSample',16)
%audiowrite("vowel_i_HRNR.wav",I2,40000,'BitsPerSample',16)
%%
%%%% Vowel E
E0 = vowel_e_noisy;
[E1 , E2] = WienerNoiseReduction(vowel_e_noisy,40000,1);
%writing to a file
%audiowrite("vowel_e_TSNR.wav",E1,40000,'BitsPerSample',16)
%audiowrite("vowel_e_HRNR.wav",E2,40000,'BitsPerSample',16)
%%
%%%% Vowel AE
AE0 = vowel_ae_noisy;
[AE1 , AE2] = WienerNoiseReduction(vowel_ae_noisy,40000,1);
%writing to a file
%audiowrite("vowel_ae_TSNR.wav",AE1,40000,'BitsPerSample',16)
%audiowrite("vowel_ae_HRNR.wav",AE2,40000,'BitsPerSample',16)
%%
%%%% Vowel A
A0 = vowel_a_noisy;
[A1 , A2] = WienerNoiseReduction(vowel_a_noisy,40000,1);
%writing to a file
%audiowrite("vowel_a_TSNR.wav",u1,40000,'BitsPerSample',16)
%audiowrite("vowel_a_HRNR.wav",u2,40000,'BitsPerSample',16)
%% SNR of TSNR and HRNR 
snr_u_ts = snr(u1);
snr_i_ts = snr(I1);
snr_e_ts = snr(E1);
snr_ae_ts = snr(AE1);
snr_a_ts = snr(A1);
%%
snr_u_hr = snr(u2);
snr_i_hr = snr(I2);
snr_e_hr = snr(E2);
snr_ae_hr = snr(AE2);
snr_a_hr = snr(A2);
%%
%%% Now I need to just plot time and frequency domain graph of TSNR and
%%% HRNR
%% Plotting the time domain graphs of TSNR & HRNR Vowels present in the third column
%% TSNR Graph in Time Domain
%%%%%%%% Vowel U
[vowel_u_aud_tsnr,fs_u_tsnr] = audioread('vowel_u_TSNR.wav');
t_u_tsnr = linspace(0,length(vowel_u_aud_tsnr)/fs_u_tsnr,length(vowel_u_aud_tsnr));
figure(5)
subplot(2,3,1)
plot(t_u_tsnr,vowel_u_aud_tsnr) %Time domain plot of Signal
title('Time Domain graph of TSNR Vowel U Signal')
xlabel('Time');
ylabel('Amplitude');
%%%%%%%%% Vowel I
[vowel_i_aud_tsnr,fs_i_tsnr] = audioread('vowel_i_TSNR.wav');
t_i_tsnr= linspace(0,length(vowel_i_aud_tsnr)/fs_i_tsnr,length(vowel_i_aud_tsnr));
% figure('WindowState','maximized')
subplot(2,3,2)
plot(t_i_tsnr,vowel_i_aud_tsnr) %Time domain plot of Signal
title('Time Domain graph of TSNR Vowel I Signal')
xlabel('Time');
ylabel('Amplitude');
%%%%%%%%% Vowel E
[vowel_e_aud_tsnr,fs_e_tsnr] = audioread('vowel_e_TSNR.wav');
t_e_tsnr= linspace(0,length(vowel_e_aud_tsnr)/fs_e_tsnr,length(vowel_e_aud_tsnr));
% figure('WindowState','maximized')
subplot(2,3,3)
plot(t_e_tsnr,vowel_e_aud_tsnr) %Time domain plot of Signal
title('Time Domain graph of TSNR Vowel E Signal')
xlabel('Time');
ylabel('Amplitude');
%%%%%%%%% Vowel AE
[vowel_ae_aud_tsnr,fs_ae_tsnr] = audioread('vowel_ae_TSNR.wav');
t_ae_tsnr= linspace(0,length(vowel_ae_aud_tsnr)/fs_ae_tsnr,length(vowel_ae_aud_tsnr));
% figure('WindowState','maximized')
subplot(2,3,4)
plot(t_ae_tsnr,vowel_ae_aud_tsnr) %Time domain plot of Signal
title('Time Domain graph of TSNR Vowel AE Signal')
xlabel('Time');
ylabel('Amplitude');
%%%%%%%%% Vowel A
[vowel_a_aud_tsnr,fs_a_tsnr] = audioread('vowel_a_TSNR.wav');
t_a_tsnr= linspace(0,length(vowel_a_aud_tsnr)/fs_a_tsnr,length(vowel_a_aud_tsnr));
% figure('WindowState','maximized')
subplot(2,3,5)
plot(t_a_tsnr,vowel_a_aud_tsnr) %Time domain plot of Signal
title('Time Domain graph of TSNR Vowel A Signal')
xlabel('Time');
ylabel('Amplitude');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HRNR Graph in Time Domain
%%%%%%%% Vowel U
[vowel_u_aud_hrnr,fs_u_hrnr] = audioread('vowel_u_HRNR.wav');
t_u_hrnr = linspace(0,length(vowel_u_aud_hrnr)/fs_u_hrnr,length(vowel_u_aud_hrnr));
figure(6)
subplot(2,3,1)
plot(t_u_hrnr,vowel_u_aud_hrnr) %Time domain plot of Signal
title('Time Domain graph of HRNR Vowel U Signal')
xlabel('Time');
ylabel('Amplitude');
%%%%%%%%% Vowel I
[vowel_i_aud_hrnr,fs_i_hrnr] = audioread('vowel_i_HRNR.wav');
t_i_hrnr= linspace(0,length(vowel_i_aud_hrnr)/fs_i_hrnr,length(vowel_i_aud_hrnr));
% figure('WindowState','maximized')
subplot(2,3,2)
plot(t_i_hrnr,vowel_i_aud_hrnr) %Time domain plot of Signal
title('Time Domain graph of HRNR Vowel I Signal')
xlabel('Time');
ylabel('Amplitude');
%%%%%%%%% Vowel E
[vowel_e_aud_hrnr,fs_e_hrnr] = audioread('vowel_e_HRNR.wav');
t_e_hrnr= linspace(0,length(vowel_e_aud_hrnr)/fs_e_hrnr,length(vowel_e_aud_hrnr));
% figure('WindowState','maximized')
subplot(2,3,3)
plot(t_e_hrnr,vowel_e_aud_hrnr) %Time domain plot of Signal
title('Time Domain graph of HRNR Vowel E Signal')
xlabel('Time');
ylabel('Amplitude');
%%%%%%%%% Vowel AE
[vowel_ae_aud_hrnr,fs_ae_hrnr] = audioread('vowel_ae_HRNR.wav');
t_ae_hrnr= linspace(0,length(vowel_ae_aud_hrnr)/fs_ae_hrnr,length(vowel_ae_aud_hrnr));
% figure('WindowState','maximized')
subplot(2,3,4)
plot(t_ae_hrnr,vowel_ae_aud_hrnr) %Time domain plot of Signal
title('Time Domain graph of HRNR Vowel AE Signal')
xlabel('Time');
ylabel('Amplitude');
%%%%%%%%% Vowel A
[vowel_a_aud_hrnr,fs_a_hrnr] = audioread('vowel_a_HRNR.wav');
t_a_hrnr= linspace(0,length(vowel_a_aud_hrnr)/fs_a_hrnr,length(vowel_a_aud_hrnr));
% figure('WindowState','maximized')
subplot(2,3,5)
plot(t_a_hrnr,vowel_a_aud_hrnr) %Time domain plot of Signal
title('Time Domain graph of HRNR Vowel A Signal')
xlabel('Time');
ylabel('Amplitude');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%