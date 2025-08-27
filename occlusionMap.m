%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Author: Patricia Amado Caballero
% Email: patricia.amado@uva.es
% Date: 2025-08-27

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function heatmap = occlusionMap(net, image, predictedLabel, occ_size, occ_stride, occ_pixel)
       % occ_size, occ_stride, occ_pixel: Occlusion Parameters

    [height, width] = size(image); 

    % Compute the dimensions of the heatmap
    output_height = ceil((height - occ_size) / occ_stride);
    output_width = ceil((width - occ_size) / occ_stride);

     % Initialize the heatmap
    heatmap= zeros(output_height, output_width);

    % Occlusion Algorithm
    for h = 0:(output_height-1)
        for w = 0:(output_width-1)
            h_start = h * occ_stride + 1;
            w_start = w * occ_stride + 1;
            h_end = h_start + occ_size - 1;
            w_end = w_start + occ_size - 1;

            if h_end > height || w_end > width
                continue;
            end

            input_image_occluded = image; 
            input_image_occluded(h_start:h_end, w_start:w_end) = occ_pixel;

         % Ensure the occluded image has the correct format for the network
         % For example, if the network expects HxWx1, use reshape:
         % input_image_occluded_reshaped = reshape(input_image_occluded, [height, width, 1]);
            [~, scores] = classify(net, input_image_occluded); 
            
          % Extract the probability of the opposite class
          % If predictedLabel is 0, get probability of class 1
          % If predictedLabel is 1, get probability of class 0
          % Assuming 'scores' is [prob_class0, prob_class1]

            if double(predictedLabel) == 0
                prob_opposite_class = scores(2); 
            else
                prob_opposite_class = scores(1); 
            end

            % Asigne heatmap
            heatmap(h + 1, w + 1) = prob_opposite_class; 
        end
    end

    % Resize the heatmap to the original image size
    heatmap = imresize(heatmap, [height, width], 'bicubic');

end
