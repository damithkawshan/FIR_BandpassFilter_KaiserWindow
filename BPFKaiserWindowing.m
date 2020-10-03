%%
close all;
clear all;
clc

% for this some parameters are calculated using the numbers in the array indx
indx=[1 6 0 1 3 4] % input your index no as a matrix

% getting A B C 
A=indx(4);
B=indx(5);
C=indx(6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%parameters for the BPF
A_p=0.05+0.01*A;
A_a=40+B;
Wp1=C*100+400;
Wp2=C*100+700;
Wa1=C*100+250;
Wa2=C*100+800;
Ws=2*(C*100+1200);
fs=Ws/(2*pi);%sampling freq
L=2048;% no of samples for plot input signals

%%%%%%%%%%%%%%%%%%%%%%%%%%Kaiser Window Generation%%%%%%%%%%%%%%%%%%%%%%%
%%
% defining parameters as needed
% get critical transition width
Bt=min((Wp1-Wa1),(Wa2-Wp2));

Wc1=Wp1-Bt/2;%ideal lower cutoff freq
Wc2= Wp2+Bt/2;%ideal upper cutoff freq

%get actual passband ripple
delt_P=(10^(0.05*A_p)-1)/(10^(0.05*A_p)+1);
delt_A=10^(-0.05*A_a);

delta=min(delt_P,delt_A);

%actual stopband attenuation & passband attenuation calculation
Aa= - 20*log10(delta);

%get parameter alpha
if Aa<=21
    alpha=0;
elseif Aa<=50
    alpha=0.5842*(Aa-21)^0.4+0.07886*(Aa-21);
else
    alpha=0.1102*(Aa-8.7);
end

%get parameter D
if Aa<=21
    D=0.9222;
else
    D=(Aa-7.95)/14.36;
end

%select N- Order of the filter
temp_var=(Ws*D/Bt)+1;
if (temp_var-fix(temp_var)) ==0
    if mod(temp_var,2)==1
        N=temp_var;
    else
        N=temp_var+1;
    end
else
    if mod(fix(temp_var),2)==0
        N=fix(temp_var)+1;
    else
        N=fix(temp_var)+2;
    end
end

%get M
M=(N-1)/2;

%set non zero domain for window
n=linspace(-M,M,2*M+1);

%get the beta values
beta=alpha*sqrt(1-(n/M).^2);

num=getBessel(beta);
denom=getBessel(alpha);

%kiaser window generation
w_n=num/denom;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%ideal BPF generation
h_n=sin(Wc2*n/fs)./(pi*n)-sin(Wc1*n/fs)./(pi*n);
h_n(M+1)=(Wc2-Wc1)*2/Ws;

%resulting FIR BP Filter
hw_n=h_n.*w_n;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Analysis%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
[h,f]=freqz(hw_n);%frequency response of filter in normalized freq range

f=f*(Ws/(2*pi));%freq domain correction
H_db=20*log10(abs(h));%get the magnitude in db scale

%causal impulse response of the filter
figure;
stem(n(M+1:2*M+1),hw_n(M+1:2*M+1));
xlabel('n (Samples)');ylabel('Amplitude');title('impulse response of the filter')

minPassRip=-A_p/2*ones(1,length(f));
maxStopRip=-A_a*ones(1,length(f));
maxPassRip=A_p/2*ones(1,length(f));

%frequency response of the filter
figure;
plot(f,H_db);
hold on;
plot(f,minPassRip);
hold on;
plot(f,maxPassRip);
hold on;
plot(f,maxStopRip);
xlabel('Frequency (rad/s)');ylabel('Magnitude (dB)');title('Frequency response of the bandpass FIR filter (Kaiser Windowed)')

figure;
plot(f,abs(h));
xlabel('Frequency (rad/s)');ylabel('Magnitude');title('BP filter');
%%%%%%%%%%%%%%%%%%%%%%validating results%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
n=(0:L-1);%integer domain
f=(-L/2:L/2-1)*(Ws/L);%frequency domain

%input signal
x_n=sin(n*(Wc1+Wc2)/(2*fs))+sin(n*(0+Wa1)/(2*fs))+sin(n*(Wa2+Ws/2)/(2*fs));
t_n=sin(n*(Wc1+Wc2)/(2*fs));%ideal output

%%%%%%%%%%%%%%%%%%%%%%%frequency domain %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get DFT of signals
X_n=fftshift(fft(x_n,L));
H_n=fftshift(fft(hw_n,L));
T_n=fftshift(fft(t_n,L));

%convolve two signals in time domain (freq domain multiplication)
Y_n=(X_n.*H_n);

%plot results
figure;
subplot(3,1,1);plot(f,abs(X_n/L));
xlabel('Frequency (rad/s)');ylabel('Magnitude');title('Input signal');
subplot(3,1,2);plot(f,abs(H_n));
xlabel('Frequency (rad/s)');ylabel('Magnitude');title('BP filter');
subplot(3,1,3);plot(f,abs(Y_n/L));
xlabel('Frequency (rad/s)');ylabel('Magnitude');title('Output signal');

figure;
subplot(2,1,1);plot(f,abs(Y_n/L));
xlabel('Frequency (rad/s)');ylabel('Magnitude');title('Output of the FIR BPF in Frequency Domain');
subplot(2,1,2);plot(f,abs(T_n/L));
xlabel('Frequency (rad/s)');ylabel('Magnitude');title('Output of the Ideal BPF in Frequency Domain');

%%%%%%%%%%%%%%%%%%%%%%%%%get time domain results%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%arrange the filter
hw_n1=zeros(1,L);
hw_n1(floor(L/2)-M:floor(L/2)+M)=hw_n;

%filter the signal
y_n=conv(x_n,hw_n,'full');

%plot results using subplots
figure;
subplot(3,1,1);plot(n,x_n);
ylim([min(x_n) max(x_n)]);xlim([0 L]);xlabel('Samples (n)');ylabel('Amplitude');title('Input signal');
subplot(3,1,2);plot(n,hw_n1);
ylim([min(hw_n1) max(hw_n1)]);xlim([0 L]);xlabel('Samples (n)');ylabel('Amplitude');title('FIR BPF in time domain ');
subplot(3,1,3);plot(n,y_n(M:L+M-1));
ylim([-1 1]);xlim([0 L]);xlabel('Samples (n)');ylabel('Amplitude');title('Output of the FIR BPF');

figure;
subplot(2,1,1);plot(n,y_n(M:L+M-1));
xlabel('Samples (n)');ylabel('Amplitude');title('Output of the FIR BPF');ylim([-1 1]);xlim([0 L]);
subplot(2,1,2);
plot(n,t_n);
ylim([-1 1]);xlim([0 L]);xlabel('Samples (n)');ylabel('Amplitude');title('Output of the ideal BPF');