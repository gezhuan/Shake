
%% Read the raw signals
clear all 
close all
clc
beads_num = [1,2,3,4,5,6,7,8,9,10,12,15,18,21,23,28,33,38,45,55,65,75]; %% Num of glass beads in box
fs = 48000; %% sampling frequency of the apparatus
dt = 1/fs;

% Define the frequency range of interest
freq_range = [0, 15000]; % in Hz

% Initialize storage for 3D data
num_freq_points = 5000; % Number of points in the frequency domain for plotting
fspan_resampled = linspace(freq_range(1), freq_range(2), num_freq_points);
signal_strengths = zeros(length(beads_num), num_freq_points);

for index = 1:length(beads_num)

    filename = "r_ball_" + num2str(beads_num(index)) + "_0.m4a";

    [signal_rawA, ~] = audioread(filename);
    nums=0;
    Ts=18;
    Tf=20;
    for j=1:length(signal_rawA)
        if j>=Ts/dt&j<=Tf/dt
            nums=nums+1;
            signal_raw(nums)=signal_rawA(j);
        end
    end
    Length = length(signal_raw);

    % Compute FFT and magnitude
    signal_raw_dft = fft(signal_raw);
    signal_raw_dft_abs = abs(signal_raw_dft) / Length; % Normalized magnitude
    fspan = fs / Length * (0:Length - 1); % Frequency domain in Hz

    % Interpolate to the resampled frequency span
    signal_raw_dft_abs_resampled = interp1(fspan, signal_raw_dft_abs, fspan_resampled, 'linear', 0);

    % Store the interpolated signal strength
    signal_strengths(index, :) = signal_raw_dft_abs_resampled;
end

% Create a 3D plot
[X, Y] = meshgrid(fspan_resampled / 1000, beads_num); % Convert frequency to kHz
Z = signal_strengths;

figure(1)
surf(X, Y, Z, 'EdgeColor', 'none');
xlabel("Freq (kHz)");xlim([3,10]);
ylabel("Num of beads")
zlabel("signal magnitude (a.u.)")
colorbar()
view(3)
grid on;
