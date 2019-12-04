% A.Crimi 2019
% 
% Code related to Kara et al. 2020
%

%%%% Red vs Green staining %%%%

function [n, n_norm, mean_int, n_pixel, n_green, intensities ] =  red_vs_green( red_file, green_file)

       % Common threshold for binarization
       bin_th = 0.1; 

       % Load files
       red_c = imread(red_file);
       green_c = imread(green_file);
     %  figure; imagesc(green_c)
    %   figure; imagesc(red_c)
   
       % Subtract background
       background_red = imopen(red_c,strel('disk',15));
       background_green = imopen(green_c,strel('disk',15));
       red_c = red_c ./ background_red;
       green_c = green_c ./ background_green;

       % Remove noise
       red_c = medfilt2(red_c);
       green_c = medfilt2(green_c);  
       
       % Automatic Contrast adjustment 
       red_adj =imadjust(red_c);       
       green_adj =imadjust(green_c);       
    %   figure; imagesc(green_adj)
    %   figure; imagesc(red_adj)
       
       % Binarize images and compute difference
       BW_red = imbinarize(red_adj,bin_th );
       BW_green = imbinarize(green_adj,bin_th);
       
       % In in Orange images there is an area larger  than a threshold it
       % means there is an artefact
       [L_red,n_red] = bwlabel(BW_red);
 
       maxvalue = 0;
       for mm = 1 : n_red
           [r,c] = find(L_red==mm);
           if(length(r) > maxvalue)
              maxvalue = length(r);
           end
       end 
       
       threshold_artFX = 150000;
     %  if (maxvalue < threshold_artFX)
        if (1)
           difference = BW_red - BW_green;
           dif = difference>0; %Discart negative values as these represent the difference green - red.
 %          figure; imagesc(dif)

           % Perform morphological operators
           % The dilation is a bit bigger to avoid over partitioning of the
           % same cells
           SE_ero = strel('rectangle',[1 1]); 
           BW_erode = imerode(dif,SE_ero);  
           SE_dil = strel('rectangle',[1 1]); 
           BW_final = imdilate(BW_erode,SE_dil);
       %    figure;imagesc(BW_final);

           % Label cells and count them
           [L,n] = bwlabel(BW_final);
   %       figure;imagesc(L);

           % Compute mean intensity within cells
           [r,c] = find(BW_final>0);

           % If no pixels left set every features to zero
           if r >0 
                    selected_pixels_int = [] ;
                    n_pixel = length(r);
                    % Normalize the total count of pixels 
                    n_pixel =  n_pixel / n;

                    for kk = 1 : n_pixel
                        selected_pixels_int(end+1) = red_adj(r(kk),c(kk));
                    end
                    mean_int = mean(selected_pixels_int);

                    % Intracells intensities
                    intensities = 0;
           %{
           intensities = zeros(n,1);
           for ll = 1 : n
               [loc_r,loc_c] = find(L==ll);
               for jj = 1 : length(loc_r)
                   temp_stack = red_adj(loc_r(jj),loc_c(jj));
               end
               intensities(ll) = mean(temp_stack);
           end
           %}

                   % Normalize the number of detected cells over the number of cells in
                   % the Green Channel
                   [L_green,n_green] = bwlabel(BW_green);
                   n_norm = n / n_green;
           else
                n = 0;
                mean_int = 0;
                n_pixel = 0;
                n_norm = 0;
                n_green = 0;
                intensities = 0;
           end
       else
           disp('Warning, artefact present')
           n  = NaN ;
           n_norm  = NaN ;
           mean_int  = NaN ;
           n_pixel = NaN ;
           n_green = NaN ;
           intensities = NaN ;
       end
       
end
