%% ECE 160 - CHALLENGE 3: "Playback Rate Modification of Audio Signals" %%

%NAME 1: ANTONIO JAVIER SAMANIEGO JURADO, PERM: 6473680
%NAME 2: PABLO MARTIN GARCIA, PERM: 6473607


function [output_audio] = audiomod(input_audio,rate,fs) 

%% Set Parameters: %%

length_input_audio=length(input_audio);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Set Window Length below: %%%%%%%%%%
windowLength=1000; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fftLength=windowLength;
hopLength=round(windowLength*0.2); 

%Create window:
w=hamming(windowLength);


%% Perform STFT: %%

stft_input_audio=zeros(ceil((fftLength+1)/2),fix((length_input_audio-windowLength)/hopLength)+1);

shifted_point=0;
for i=1:size(stft_input_audio,2)
    %Window position:
    i_start=shifted_point+1; 
    i_stop=shifted_point+windowLength;
    shifted_point=shifted_point+hopLength;
    
    %Take part of the signal corresponding to window position and multiply.
    input_after_window=w.*input_audio(i_start:i_stop);
    
    %Perform FFT to that part of the signal.
    fft_after_window=fft(input_after_window,fftLength); 
    stft_input_audio(:,i)=fft_after_window(1:size(stft_input_audio,1)); %Each FFT corresponds to a certain group of samples in time, which results in the spectrogram.
end




%% Interpolate STFT for Every Group of Samples Depending on Rate: %%

%Set which time positions of the spectrogram to perform interpolation based on rate:
time_jumps=linspace(1,size(stft_input_audio,2)-1,ceil((size(stft_input_audio,2))/rate));

stft_phase = angle(stft_input_audio(:,1));
for i=1:length(time_jumps)

    j=ceil(time_jumps(i)-1); 
    
    %Take neighbor columns, to interpolate their frequencies.
    stft_instant1=stft_input_audio(:,j+1); 
    stft_instant2=stft_input_audio(:,j+2);
    
    %Interpolate magnitude and use interpolated phase to update the column.
    magnitude_interp=abs(stft_instant1)+(abs(stft_instant2)); 
    stft_interp(:,i)=magnitude_interp.*exp(1j*stft_phase); 
    stft_phase=stft_phase+(angle(stft_instant2)-angle(stft_instant1)); %Interpolate/Update phase.

end
    
    
    
%% Inverse STFT: Signal Reconstruction %%

length_output=length(time_jumps)*hopLength+windowLength;
output_audio=zeros(1,length_output);

w=w';

%Set positions that will adjust output signal based on hops and window overlap:
hops=linspace(0,(hopLength*(length(time_jumps)-1)),ceil(((hopLength*(length(time_jumps))))/hopLength));

for k=1:length(time_jumps)
    
    m=hops(k);
    
    inverse_window=stft_interp(1:end-1,k)';
    
    %Calculate and concatenate the N/2 complex conjugate points of the FFT of each
    %window (each spectrogram column), as we need them to reconstruct the initial signal:
    conj_inverse_window=zeros(1,(windowLength/2));
    for n=(windowLength/2):-1:2
        conj_inverse_window(n)=conj(inverse_window(n));
    end
    conc_inverse_window=[inverse_window,conj_inverse_window];
    
    %Set window position for output signal:
    window_pos=(m+1):(m+windowLength);
    
    %Perform IFFT for each window:
    inverse_fft_window=real(ifft(conc_inverse_window));
    
    %Take part of the output signal corresponding to window position and multiply.
    output_audio(window_pos)=output_audio(window_pos)+w.*inverse_fft_window;
end
   
%% Plots: %%

%Spectogram of input audio:
t_init=windowLength/2;
t_final=windowLength/2+(size(stft_input_audio,2)-1)*hopLength;
t_stft=(t_init:hopLength:t_final)/fs;

freq=(0:size(stft_input_audio,1)-1)*fs/fftLength;

figure(1)
imagesc(t_stft,freq,10*log10(abs(stft_input_audio)))
view(270,90)
title('Spectrogram - Input Audio')
xlabel('Time (s)')
ylabel('Frequency (Hz)')

%Input/Ouput Audio - Time Domain
t_input=(0:length(input_audio)-1)/fs;
t_output=(0:length(output_audio)-1)/fs;

figure(2)
subplot(2,2,1)
plot(t_input,input_audio);
title('Input Audio - Time Domain')
xlabel('Time (s)')
ylabel('Amplitude')
subplot(2,2,2)
plot(t_output,output_audio);
title('Output Audio - Time Domain')
xlabel('Time (s)')
ylabel('Amplitude')

%Input/Ouput Audio - Frequency Domain
NFFT1=2^nextpow2(length(input_audio));
S1=fft(input_audio,NFFT1)/length(input_audio);
NFFT2=2^nextpow2(length(output_audio));
S2=fft(output_audio,NFFT2)/length(output_audio);

freq1=fs/2*linspace(0,1,NFFT1/2+1);
freq2=fs/2*linspace(0,1,NFFT2/2+1);

subplot(2,2,3)
plot(freq1,2*abs(S1(1:NFFT1/2+1)),'red') 
title('Input Audio - Frequency Domain')
xlabel('Frequency (Hz)')
ylabel('|S(f)|')
subplot(2,2,4)
plot(freq2,2*abs(S2(1:NFFT2/2+1)),'red') 
title('Output Audio - Frequency Domain')
xlabel('Frequency (Hz)')
ylabel('|S(f)|')


end