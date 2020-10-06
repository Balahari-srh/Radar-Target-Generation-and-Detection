# Radar-Target-Generation-and-Detection



###Range from First FFT
![Image of Range from First FFT](https://github.com/Balahari-srh/Radar-Target-Generation-and-Detection/blob/main/images/Range%20from%20First%20FFT.jpg)

###Range Doppler Response
![Image of Range Doppler Response Map](https://github.com/Balahari-srh/Radar-Target-Generation-and-Detection/blob/main/images/Range%20Dopple%20Response%20Map.jpg)
```matlab
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

```
###Final CFAR Output
![Image of CFAR Range Doppler Map](https://github.com/Balahari-srh/Radar-Target-Generation-and-Detection/blob/main/images/CFAR%20Range%20Doppler%20Map.jpg)
