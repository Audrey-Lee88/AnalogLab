% run sdrsetup
%% Exercise 1
fs = 300e3; % this is the sample rate
fc = 90.8e6; % this is the center frequency

t = linspace(1,fs*5,fs*5+160);
t = t.';
x = zeros(fs*5,1); % empty vector to store data

% create object for RTL-SDR receiver
rx = comm.SDRRTLReceiver('CenterFrequency',fc, 'EnableTunerAGC', false, 'TunerGain', 35,  'SampleRate', fs);

counter = 1; % initialize a counter
while(counter < length(x)) % while the buffer for data is not full
    rxdata = rx();   % read from the RTL-SDR
    x(counter:counter + length(rxdata)-1) = rxdata; % save the samples returned
    counter = counter + length(rxdata); % increment counter
end
% the data are returned as complex numbers
% separate real and imaginary part, and remove any DC offset
y_I = real(x)-mean(real(x));
y_Q = imag(x)-mean(imag(x));
y_I2 = y_I;

%% part c
plot_FT(y_I2,fs); % Plot FFT of y_I
hold on;
xlabel("Frequency (Hz)")
ylabel("Magnitude")
hold off;
%% part d
% plot y_I in time domain
t2 = t./fs;
plot(t2,y_I2)
hold on;
xlabel("Time (Sec)")
ylabel("Magnitude")
hold off;
%% part e
% multiply by cos to put blob at center
y_I2 = y_I2 .* cos(2.*pi.*100000.*t);
plot_FT(y_I2,fs); % Plot FFT of y_I after cos
hold on;
xlabel("Frequency (Hz)")
ylabel("Magnitude")
hold off;
%%
% take derivative of signal
dydx = diff(y_I2(:))./diff(t(:));
% zero out negative components
dydx(dydx < 0) = 0;
% normalize it
dydx = dydx./max(dydx);
% plot it in time domain
plot(t(1:end-1,:)/fs,dydx)
hold on;
xlabel("Time (Sec)")
ylabel("Magnitude")
hold off;
%% part f
% dydx = dydx .* cos(2.*pi.*100000.*t);

% create a LPF
vals = linspace(-5,5,101);
h = sinc(vals);
% plot_FT(h,fs)

output_I = conv(dydx, h);
output_I_norm = output_I/max(output_I);
% plot_FT(output_I,fs);

% plot in time domain both derivative and LPF
plot(t(1:end-1,:)/fs,dydx)
hold on;
plot(t(1:end-1,:)/fs,output_I_norm(1:length(t)-1,:))
xlabel("Time (Sec)")
ylabel("Magnitude")
legend("Derivative","LPF")
hold off;

%% part g
% get the output and play the sound
output = output_I - mean(output_I);
output_norm = output ./ (max(output).*10);

out = decimate(output_norm,4);
sound(out,fs/4)

%%
t_by4 = linspace(t(1),t(end),round(length(output_I_norm)/4)).';
%%

% before decim
plot(t(1:end-1,:)/fs,output_norm(1:length(t)-1,:), "ro-")
hold on;
% after decim
plot(t_by4/fs,out, "bo-")
xlabel("Time (Sec)")
ylabel("Magnitude")
legend("Before Decimation", "After Decimation")
hold off;

%% save audio file
filename = 'exercise1.wav';
audiowrite(filename,out,fs/4);

%% replay saved audio file
[y,Fs] = audioread(filename);
sound(y,Fs);
%% Exercise 2
fs = 300e3; % this is the sample rate
fc = 107.9e6; % this is the center frequency

t = linspace(1,fs*5,fs*5+160);
t = t.';
x = zeros(fs*5,1); % empty vector to store data

% create object for RTL-SDR receiver
rx2 = comm.SDRRTLReceiver('CenterFrequency',fc, 'EnableTunerAGC', false, 'TunerGain', 35,  'SampleRate', fs);

counter = 1; % initialize a counter
while(counter < length(x)) % while the buffer for data is not full
    rxdata = rx2();   % read from the RTL-SDR
    x(counter:counter + length(rxdata)-1) = rxdata; % save the samples returned
    counter = counter + length(rxdata); % increment counter
end
% the data are returned as complex numbers
% separate real and imaginary part, and remove any DC offset
y_I_2 = real(x)-mean(real(x));
y_Q = imag(x)-mean(imag(x));
%% part b
% take derivative of y_I and y_Q
dyIdx =  diff(y_I_2(:))./diff(t(:));
dyQdx = diff(y_Q(:))./diff(t(:));

% plug into m_hat equation
m_hat = dyQdx .* y_I_2(1:length(dyQdx)) - (dyIdx .* y_Q(1:length(dyQdx)));
%%
% plot m_hat in time domain
plot(t(1:end-1,:)/fs,m_hat)
hold on
xlabel("Time (Sec)")
ylabel("Magnitude")
hold off
%%
% play m_hat sound
out_m = 0.4 * decimate(m_hat,4) / max(m_hat);
sound(out_m,fs/4)

%% save audio file
filename = 'exercise2.wav';
audiowrite(filename,out_m,fs/4);

%% replay saved audio file
[y,Fs] = audioread(filename);
sound(y,Fs);