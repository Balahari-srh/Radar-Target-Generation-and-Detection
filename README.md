# Radar-Target-Generation-and-Detection

###Range from First FFT
![Image of Range from First FFT](https://github.com/Balahari-srh/Radar-Target-Generation-and-Detection/blob/main/images/Range%20from%20First%20FFT.jpg)

###Range Doppler Response
![Image of Range Doppler Response Map](https://github.com/Balahari-srh/Radar-Target-Generation-and-Detection/blob/main/images/Range%20Dopple%20Response%20Map.jpg)

Write a README outlining the following:
1.Implementation steps for the 2D CFAR process.



2.Selection of Training, Guard cells and offset.
  -The training cells, Guard Cells and Offset are chosen after different trials based on the scenario.
  -When offset is too small, it is unable to seperate the signals from noise
  -When Training Cells is too small, we keep missing the actual target signal.
  -When Training Cells is too large, it is unable to find the target signal and we get a plane.
  -When Guard Cells are too smal, ther is very little leakage of target signal.
  -When Guard Cells are too high, there is leakage to other bins.
```matlab
  %Select the number of Training Cells in both the dimensions.
  Tr = 16; % 12 16 10
  Td = 8; % 6 8 4
  %Select the number of Guard Cells in both dimensions around the Cell under
  %test (CUT) for accurate estimation
  Gr = 5; %5 8
  Gd = 2; %2 4
  % offset the threshold by SNR value in dB
  offset = 7;
  %The training cell size without guard cells
  TrainingCells = (2*Tr+2*Gr+1)*(2*Td*Gd*1)-(2*Gr+1)*(2*Gd+1);

```

3.Steps taken to suppress the non-thresholded cells at the edges.
  -Two matrix are declared of same size as RDM, to hold threshold values(threshold_cfar) and the final output signal(Signal_cfar).
  ```matlab
    % Vector to hold threshold values
    threshold_cfar = zeros(r,c);

    %Vector to hold final signal after thresholding
    signal_cfar = zeros(r,c);
```
  -The threshold is calculated by traversing along training cells.
  -The total threshold is calculated and the values of guard cells are later subtracted.
  -Once the noise level is calculated for training cells, the average is calculated by dividing the values by number of training cells.
  -Offset is added to the noise, also converting back it to logarithmic scale using pow2db.
  -The signals under CUT are compared against threshold.
  -If CUT level>threshold assigned a value 1.(I am not assigning 0 as they are already a zeros matrix same of RDM size.)
```matlab

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

```
###Final CFAR Output
![Image of CFAR Range Doppler Map](https://github.com/Balahari-srh/Radar-Target-Generation-and-Detection/blob/main/images/CFAR%20Range%20Doppler%20Map.jpg)
