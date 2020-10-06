clear all
clc;

%% Radar Specifications 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency of operation = 77GHz
% Max Range = 200m
% Range Resolution = 1 m
% Max Velocity = 100 m/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%speed of light = 3e8
%% User Defined Range and Velocity of target
% *%TODO* :
% define the target's initial position and velocity. Note : Velocity
% remains contant
ir = 110; %Initial range cannot exceed max value of 200 m 100
iv = -20; %Initial velocity can be any value in range -70 to 70 m/s 50

%% FMCW Waveform Generation
dres = 1; %Range Resolution
c = 3e8; %Speed of Light
Rmax = 200; %Max Range

% *%TODO* :
%Design the FMCW waveform by giving the specs of each of its parameters.
% Calculate the Bandwidth (B), Chirp Time (Tchirp) and Slope (slope) of the FMCW
% chirp using the requirements above.
B_sweep = c/(2*dres); 
Tchirp = (5.5*2*Rmax)/c;
slope = B_sweep/Tchirp;

%Operating carrier frequency of Radar 
fc= 77e9;             %carrier freq

                                                          
%The number of chirps in one sequence. Its ideal to have 2^ value for the ease of running the FFT
%for Doppler Estimation. 
Nd=128;                   % #of doppler cells OR #of sent periods % number of chirps

%The number of samples on each chirp. 
Nr=1024;                  %for length of time OR # of range cells

% Timestamp for running the displacement scenario for every sample on each
% chirp
t=linspace(0,Nd*Tchirp,Nr*Nd); %total time for samples


%Creating the vectors for Tx, Rx and Mix based on the total samples input.
Tx=zeros(1,length(t)); %transmitted signal
Rx=zeros(1,length(t)); %received signal
Mix = zeros(1,length(t)); %beat signal

%Similar vectors for range_covered and time delay.
r_t=zeros(1,length(t));
td=zeros(1,length(t));


%% Signal generation and Moving Target simulation
% Running the radar scenario over the time. 

for i=1:length(t)         
    
    
    % *%TODO* :
    %For each time stamp update the Range of the Target for constant velocity. 
    rt= ir + t(i)* iv; %range at particular instance of time 
    td = 2*rt/c; %time taken for travel
    % *%TODO* :
    %For each time sample we need update the transmitted and
    %received signal. 
    Tx(i) = cos(2*pi*(fc*t(i) + ((slope*t(i)^2)/2)));
    Rx (i) = cos(2*pi*(fc*(t(i)-td) + ((slope*(t(i)-td)^2)/2)));
    
    % *%TODO* :
    %Now by mixing the Transmit and Receive generate the beat signal
    %This is done by element wise matrix multiplication of Transmit and
    %Receiver Signal
    Mix(i) = Tx(i).*Rx(i);
    
end

%% RANGE MEASUREMENT


 % *%TODO* :
%reshape the vector into Nr*Nd array. Nr and Nd here would also define the size of
%Range and Doppler FFT respectively.
Mix_Matrix = reshape(Mix, [Nr,Nd]);
 % *%TODO* :
%run the FFT on the beat signal along the range bins dimension (Nr) and
%normalize.
s_fft = fft(Mix_Matrix,Nr)./Nr;
 % *%TODO* :
% Take the absolute value of FFT output
s_fft = abs(s_fft);
 % *%TODO* :
% Output of FFT is double sided signal, but we are interested in only one side of the spectrum.
% Hence we throw out half of the samples.
s_fft  = s_fft(1:length(t)/2+1);

%plotting the range
figure ('Name','Range from First FFT')
subplot(2,1,1)

 % *%TODO* :
 % plot FFT output 

plot(s_fft) 
axis([0 200 0 1]);





%% RANGE DOPPLER RESPONSE
% The 2D FFT implementation is already provided here. This will run a 2DFFT
% on the mixed signal (beat signal) output and generate a range doppler
% map.You will implement CFAR on the generated RDM


% Range Doppler Map Generation.

% The output of the 2D FFT is an image that has reponse in the range and
% doppler FFT bins. So, it is important to convert the axis from bin sizes
% to range and doppler based on their Max values.

Mix=reshape(Mix,[Nr,Nd]);

% 2D FFT using the FFT size for both dimensions.
sig_fft2 = fft2(Mix,Nr,Nd);

% Taking just one side of signal from Range dimension.
sig_fft2 = sig_fft2(1:Nr/2,1:Nd);
sig_fft2 = fftshift (sig_fft2);
RDM = abs(sig_fft2);
RDM = 10*log10(RDM) ;


%use the surf function to plot the output of 2DFFT and to show axis in both
%dimensions
doppler_axis = linspace(-100,100,Nd);
range_axis = linspace(-200,200,Nr/2)*((Nr/2)/400);
figure,surf(doppler_axis,range_axis,RDM);
title('Range Doppler Response Map')


%% CFAR implementation

%Slide Window through the complete Range Doppler Map

% *%TODO* :
%Select the number of Training Cells in both the dimensions.
Tr = 16; % 12 16 10
Td = 8; % 6 8 4
% *%TODO* :
%Select the number of Guard Cells in both dimensions around the Cell under 
%test (CUT) for accurate estimation
Gr = 5; %5 8
Gd = 2; %2 4

% *%TODO* :
% offset the threshold by SNR value in dB
offset = 7; %1.4
% *%TODO* :
%Create a vector to store noise_level for each iteration on training cells
noise_level = zeros(1,1);


% *%TODO* :
%design a loop such that it slides the CUT across range doppler map by
%giving margins at the edges for Training and Guard Cells.
%For every iteration sum the signal level within all the training
%cells. To sum convert the value from logarithmic to linear using db2pow
%function. Average the summed values for all of the training
%cells used. After averaging convert it back to logarithimic using pow2db.
%Further add the offset to it to determine the threshold. Next, compare the
%signal under CUT with this threshold. If the CUT level > threshold assign
%it a value of 1, else equate it to 0.


   % Use RDM[x,y] as the matrix from the output of 2D FFT for implementing
   % CFAR
GridSize = (2*Tr+2*Gr+1)*(2*Td*Gd*1);

%The training cell size without guard cells
TrainingCells = (2*Tr+2*Gr+1)*(2*Td*Gd*1)-(2*Gr+1)*(2*Gd+1);

%Getting the size of RDM to declare threshold_cfar matrix to hold
%thresholds and signal_cfar to hold signal value after thresholding
[r,c] = size(RDM);

% Vector to hold threshold values 
threshold_cfar = zeros(r,c);

%Vector to hold final signal after thresholding
signal_cfar = zeros(r,c);

%Cell Under Test

for i = (Tr+Gr+1):(r -2*Tr-2*Gr)
    for j = (Td+Gd+1) : (c - 2*Td-2*Gd)
        %Summing the signal level across the complete grid.
        threshold_cfar(i,j) = sum(sum(db2pow(RDM(i-(Tr+Gr) : i+(Tr+Gr),j-(Td+Gd) : j+(Td+Gd)))));           
        %Now subtracting the Guard Cells noise from the calculated noise sum of whole grid.
        threshold_cfar(i,j) = threshold_cfar(i,j) - sum(sum(db2pow(RDM((i-Gr):(i+Gr),(j-Gd):(j+Gd)))));
        %taking the average of noise
        threshold_cfar(i,j) = threshold_cfar(i,j)/TrainingCells;
        %Adding offset and converting back to log scale.
        threshold_cfar(i,j) = offset + pow2db(threshold_cfar(i,j));
        %Compare CUT(RDM(i,j) against threshold_cfar
        if(RDM(i,j) > threshold_cfar(i,j))
            signal_cfar(i,j) = 1;
        %else
            %signal_cfar(i,j) = 0;
        end
    end
    
end






% *%TODO* :
% The process above will generate a thresholded block, which is smaller 
%than the Range Doppler Map as the CUT cannot be located at the edges of
%matrix. Hence,few cells will not be thresholded. To keep the map size same
% set those values to 0. 
 
%As I am using the signal_cfar which I have declared, I am not using the
%RDM variable, hence it is not necessary.




% *%TODO* :
%display the CFAR output using the Surf function like we did for Range
%Doppler Response output.
%figure,surf(doppler_axis,range_axis,'replace this with output');
%colorbar;
figure('Name', 'CFAR output')
surf(doppler_axis, range_axis, signal_cfar);
colorbar;
title('CFAR Range Doppler Map')


