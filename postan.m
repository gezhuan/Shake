close all
clear all
clc


df=dir('r_ball_*m4a')
Ts=15;
Tf=20;
fs=48000;
dt=1/fs;
freq_range = [0, 15000]; % in Hz
beads_num=[1,2,3,4,5,6,7,8,9,10,12,15,18,21,23,28,33,38,45,55,65,75];
% Initialize storage for 3D data
num_freq_points = 5000; % Number of points in the frequency domain for plotting
fspan_resampled = linspace(freq_range(1), freq_range(2), num_freq_points);
signal_strengths = zeros(length(beads_num), num_freq_points);

for i=1:length(df)
    clear signal Fs
    [signalA, Fs] = audioread(df(i).name);
    nums=0;

    for j=1:length(signalA)
        if j>=Ts/dt&j<=Tf/dt
            nums=nums+1;
            signal(nums)=signalA(j);
        end
    end
    
    for j=1:length(signal)
        data(i,j)=abs(signal(j));
        time(i,j)=j*dt;
    end
    [p{i},x{i}] = hist(abs(signal),80);
    Sig{i}=signal;

    Length = length(signal);

    % Compute FFT and magnitude
    signal_raw_dft = fft(signal);
    signal_raw_dft_abs = abs(signal_raw_dft) / Length; % Normalized magnitude
    fspan = fs / Length * (0:Length - 1); % Frequency domain in Hz

    % Interpolate to the resampled frequency span
    signal_raw_dft_abs_resampled = interp1(fspan, signal_raw_dft_abs, fspan_resampled, 'linear', 0);

    % Store the interpolated signal strength
    signal_strengths(i, :) = signal_raw_dft_abs_resampled;

end


NPLOT=10;
cog=jet(NPLOT);
for j=1:NPLOT
    col{j}=cog(j,:);
end
legendLabels = cell(1,NPLOT);  

% 循环绘图并构建图例标签  

  
% 循环结束后，添加图例  

h=figure(1); 
subplot(1,2,1)
for j=1:NPLOT
    plot(x{j},p{j},'color',col{j});
    legendLabels{j} = num2str(beads_num(j));  
    hold on
end
legend(legendLabels); 
set(gca, 'YScale', 'log') % 将Y轴设置为对数刻度
set(gca, 'XScale', 'log') % 将Y轴设置为对数刻度

xlabel('强度');  
ylabel('数量');  

subplot(1,2,2)
for j=1:NPLOT
    plot([time(j,:)],Sig{j},'color',col{j})
    hold on;
end
set(gcf, 'unit', 'centimeters', 'position', [4 4 32 12]);


figure(2)
[X, Y] = meshgrid(fspan_resampled / 1000, beads_num); % Convert frequency to kHz
Z = signal_strengths;

surf(X, Y, Z, 'EdgeColor', 'none');
xlabel("Freq (kHz)");xlim([0,10]);
ylabel("Num of beads")
zlabel("signal magnitude (a.u.)")
colorbar()
view(3)
grid on;

