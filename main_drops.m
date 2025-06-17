clear; clc; close all

%% Cargar carpetas
folderList = {'C1-3_FTest20','C1-2_Dtest6','C1-4_Dtest21'};  % Lista de carpetas

%C1-2_Dtest6
%C1-3_FTest20
%C1-4_Dtest21
%CNB_040425_02_bf_20250510_2154
%CNB_040425_03_20250511_1148
%CNB_040425_05_bf_20250511_1330
%CNB_040425_04_20250511_1307
%CNB_040425_06_20250511_1518

% Inicializar estructura para almacenar datos
dropData = struct();

% Ruta base donde están las carpetas (ajusta esto si necesario)
baseFolder = '/Users/anabellehassan/Documents/Biomedical eng/TFG/Segment/';

for i = 1:length(folderList)
    dropFolder = folderList{i};
    dataFolder = fullfile(baseFolder, dropFolder);
    
    if ~isfolder(dataFolder)
        error('La carpeta "%s" no existe.', dataFolder);
    end
    
    fprintf('Procesando carpeta: %s\n', dataFolder);
    
    % Upload the files
    data = load(fullfile(dataFolder, 'combined_segmentation.mat'));
    bboxes = load(fullfile(dataFolder, 'bboxes.mat'));
    areas = load(fullfile(dataFolder, 'areas.mat'));
    drop = load(fullfile(dataFolder, 'first_mask_drop.mat'));
    
    % making sure that we are not using duplicate values
    [unique_bboxes, unique_idx] = unique(bboxes.bboxes, 'rows');
    if length(unique_idx) < size(bboxes.bboxes, 1)
        fprintf('Drop %s: Se eliminaron %d bounding boxes duplicados.\n', ...
            dropFolder, size(bboxes.bboxes, 1) - length(unique_idx));
    end
    % Actualizar estructuras asociadas
    bboxes.bboxes = unique_bboxes;
    areas.areas = areas.areas(unique_idx);
    
    % Cargar la imagen .tif (basado en el nombre de la carpeta)
    [~, folderName] = fileparts(dataFolder);
    parts = split(folderName, '_');
    lastIsTime = ~isempty(regexp(parts{end}, '^\d{4}$', 'once'));
    secondLastIsDate = ~isempty(regexp(parts{end-1}, '^\d{8}$', 'once'));
    if length(parts) > 2 && lastIsTime && secondLastIsDate
        baseName = strjoin(parts(1:end-2), '_');
    else
        baseName = folderName;
    end
    tifFileName = [baseName, '.tif'];
    tifFullPath = fullfile(dataFolder, tifFileName);
    original_image = imread(tifFullPath);
    
    % Calcular información de la primera máscara
    first_mask_info.Area = double(drop.area);
    first_mask_info.Radius = sqrt(first_mask_info.Area / pi);
    [y, x] = find(drop.segmentation);
    first_mask_info.Center = [mean(x), mean(y)];
    
    % Calcular elipses concéntricas
    % Calculate region properties of the original segmentation
    props = regionprops(drop.segmentation, 'MajorAxisLength', 'MinorAxisLength', 'Orientation', 'Centroid');

    % If multiple regions found, merge into one binary mask
    if length(props) > 1
        warning('Multiple regions found in drop.segmentation. Merging into a single region.');
        drop.segmentation = drop.segmentation > 0;  % Merge all non-zero areas into one mask
        props = regionprops(drop.segmentation, 'MajorAxisLength', 'MinorAxisLength', 'Orientation', 'Centroid');  % Recalculate props
    end

    % Now safe to use props(1) because only one region exists
    majorAxis = props(1).MajorAxisLength;
    minorAxis = props(1).MinorAxisLength;
    orientation = -deg2rad(props(1).Orientation);
    centroid = props(1).Centroid;
    num_ellipses = 10;
    
    ellipse_info = struct('ScaleFactor', {}, 'MajorAxis', {}, 'MinorAxis', {}, 'Centroid', {}, 'Orientation', {}, 'Points', {});
    for j = 1:num_ellipses
        scale_factor = (num_ellipses - j + 1) / num_ellipses;
        scaled_majorAxis = scale_factor * majorAxis;
        scaled_minorAxis = scale_factor * minorAxis;
        theta = linspace(0, 2 * pi, 100);
        x_standard = (scaled_majorAxis / 2) * cos(theta);
        y_standard = (scaled_minorAxis / 2) * sin(theta);
        R = [cos(orientation), -sin(orientation); sin(orientation), cos(orientation)];
        rotated_points = R * [x_standard; y_standard];
        x_rotated = rotated_points(1, :) + centroid(1);
        y_rotated = rotated_points(2, :) + centroid(2);
        
        ellipse_info(j).ScaleFactor = scale_factor;
        ellipse_info(j).MajorAxis = scaled_majorAxis;
        ellipse_info(j).MinorAxis = scaled_minorAxis;
        ellipse_info(j).Centroid = centroid;
        ellipse_info(j).Orientation = orientation;
        ellipse_info(j).Points = [x_rotated; y_rotated];
    end
    
    % Calcular centroides de los bounding boxes
    num_bboxes = size(bboxes.bboxes, 1);
    centroids = zeros(num_bboxes, 2);
    for idx = 1:num_bboxes
        bbox = bboxes.bboxes(idx, :);
        x_start = max(1, round(bbox(1)));
        y_start = max(1, round(bbox(2)));
        x_end = min(size(data.combined_segmentation, 2), round(bbox(1) + bbox(3) - 1));
        y_end = min(size(data.combined_segmentation, 1), round(bbox(2) + bbox(4) - 1));
        mask_img = zeros(size(data.combined_segmentation));
        mask_img(y_start:y_end, x_start:x_end) = 1;
        
        I_sum = sum(mask_img(:));
        if I_sum > 0
            [x, y] = meshgrid(1:size(mask_img, 2), 1:size(mask_img, 1));
            x_center = sum(mask_img(:) .* x(:)) / I_sum;
            y_center = sum(mask_img(:) .* y(:)) / I_sum;
            centroids(idx, :) = [x_center, y_center];
        else
            centroids(idx, :) = [NaN, NaN];
        end
    end
    
    % Guardar en la estructura
    dropData(i).name = dropFolder;
    dropData(i).combined_segmentation = single(data.combined_segmentation);
    dropData(i).bboxes = bboxes.bboxes;
    dropData(i).areas = areas.areas;
    dropData(i).drop = drop;
    dropData(i).original_image = original_image;
    dropData(i).first_mask_info = first_mask_info;
    dropData(i).ellipse_info = ellipse_info;
    dropData(i).centroids = centroids;
end

fprintf('Carga de datos completada para todas las gotas.\n');

%% Display masks and ellipses
figure('Units', 'normalized', 'Position', [0, 0.2, 1, 0.8]);

for i = 1:length(dropData)
    % Subplot 1: Display the combined segmentation mask
    subplot(length(dropData), 2, (i-1)*2 + 1);
    imagesc(dropData(i).combined_segmentation);
    colormap('jet');
    colorbar;
    title(['Combined Segmentation Mask - ', dropData(i).name], 'FontSize', 12);
    axis image;
    set(gca, 'FontSize', 10);
    
    % Subplot 2: Display individual masks
    subplot(length(dropData), 2, (i-1)*2 + 2);
    label_image = bwlabel(dropData(i).combined_segmentation > 0);
    colored_masks_rgb = label2rgb(label_image, 'jet', 'k', 'shuffle');
    imshow(colored_masks_rgb);
    title(['Individual Masks - ', dropData(i).name], 'FontSize', 12);
    axis image;
    set(gca, 'FontSize', 10);
end

%% Display Ellipses

figure('Name', 'Concentric Ellipses for Each Drop', 'NumberTitle', 'off', 'Units', 'normalized', 'Position', [0, 0.2, 1, 0.8]);
for i = 1:length(dropData)
    subplot(1, length(dropData), i);
    label_image = bwlabel(dropData(i).combined_segmentation > 0);
    colored_masks_rgb = label2rgb(label_image, 'jet', 'k', 'shuffle');
    imshow(colored_masks_rgb);
    hold on;
    for j = 1:10
        [x_ellipse, y_ellipse] = pol2cart(linspace(0, 2*pi, 100), j * dropData(i).first_mask_info.Radius / 10);
        plot(x_ellipse + dropData(i).first_mask_info.Center(1), y_ellipse + dropData(i).first_mask_info.Center(2), 'w-', 'LineWidth', 1.5);
    end
    hold off;
    title(['Ellipses - ', dropData(i).name], 'FontSize', 12);
    axis image;
end

%% Crystal Density Comparison

num_ellipses = 10;  % Asegúrate que este valor es coherente
R_over_Rc = linspace(0.1, 1.0, num_ellipses);  % centro → borde (creciente)

figure('Name', 'Crystal Density Comparison', 'NumberTitle', 'off', 'Units', 'normalized', 'Position', [0, 0.2, 1, 0.8]);

% Subplot 1: Individual Crystal Density for Each Drop
subplot(1, 2, 1);
hold on;
all_densities = [];
for i = 1:length(dropData)
    crystal_density = rand(1, num_ellipses);  % Reemplaza esto con tus datos reales
    crystal_density = flip(crystal_density);  % ← Asegura que va centro → borde
    plot(R_over_Rc, crystal_density, 'LineWidth', 1.5, 'DisplayName', dropData(i).name);
    all_densities = [all_densities; crystal_density];
end
hold off;
legend;
xlabel('R / R_c');
ylabel('Crystal Density');
title('Crystal Density per Drop');
xticks(R_over_Rc);
xticklabels(arrayfun(@(r) sprintf('%.1f', r), R_over_Rc, 'UniformOutput', false));
grid on;

% Subplot 2: Mean, Max, and Min Crystal Density
subplot(1, 2, 2);
mean_density = mean(all_densities, 1);
max_density = max(all_densities, [], 1);
min_density = min(all_densities, [], 1);
hold on;
plot(R_over_Rc, mean_density, 'k-', 'LineWidth', 2, 'DisplayName', 'Mean Density');
plot(R_over_Rc, max_density, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Max Density');
plot(R_over_Rc, min_density, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Min Density');
hold off;
legend;
xlabel('R / R_c');
ylabel('Crystal Density');
title('Mean, Max, and Min Crystal Density');
xticks(R_over_Rc);
xticklabels(arrayfun(@(r) sprintf('%.1f', r), R_over_Rc, 'UniformOutput', false));
grid on;


%% Calcular el tensor de estructura para cada gota

for i = 1:length(dropData)
    num_bboxes = size(dropData(i).bboxes, 1);
    tensor_results = struct('BBox', {}, 'Tensor', {}, 'Lambda', {}, 'Orientation', {}, 'Elongation', {}, 'Distance', {}, 'Misalignment', {});
    
    for idx = 1:num_bboxes
        bbox = dropData(i).bboxes(idx, :);
        x_start = max(1, round(bbox(1)));
        y_start = max(1, round(bbox(2)));
        x_end = min(size(dropData(i).original_image, 2), round(bbox(1) + bbox(3) - 1));
        y_end = min(size(dropData(i).original_image, 1), round(bbox(2) + bbox(4) - 1));
        mask_img = zeros(size(dropData(i).original_image));
        mask_img(y_start:y_end, x_start:x_end) = 1;
        
        masked_img = double(dropData(i).original_image) .* mask_img;
        
        I_sum = sum(masked_img(:));
        if I_sum == 0
            continue;
        end
        
        [X, Y] = meshgrid(1:size(mask_img, 2), 1:size(mask_img, 1));
        x_center = sum(masked_img(:) .* X(:)) / I_sum;
        y_center = sum(masked_img(:) .* Y(:)) / I_sum;
        
        T_xx = sum((X(:) - x_center).^2 .* masked_img(:)) / I_sum;
        T_yy = sum((Y(:) - y_center).^2 .* masked_img(:)) / I_sum;
        T_xy = sum((X(:) - x_center) .* (Y(:) - y_center) .* masked_img(:)) / I_sum;
        T = [T_xx, T_xy; T_xy, T_yy];
        
        [eig_vecs, eig_vals] = eig(T);
        lambda1 = eig_vals(1, 1);
        lambda2 = eig_vals(2, 2);
        [lambda_max, max_idx] = max([lambda1, lambda2]);
        n_max = eig_vecs(:, max_idx);
        if n_max(2) < 0
            n_max = -n_max;
        end
        
        elongation = lambda_max / min(lambda1, lambda2);
        orientation = atan2(n_max(2), n_max(1));
        
        center_vector = [dropData(i).first_mask_info.Center(1) - x_center; dropData(i).first_mask_info.Center(2) - y_center];
        distance = norm(center_vector);
        center_vector_norm = center_vector / distance;
        
        misalignment_angle = abs((n_max(1) * center_vector(2)) - (n_max(2) * center_vector(1))) / (norm(n_max) * norm(center_vector));
        
        if isnan(misalignment_angle)
            misalignment_angle = 0;
        end
        
        tensor_results(idx).BBox = bbox;
        tensor_results(idx).Tensor = T;
        tensor_results(idx).Lambda = [lambda1, lambda2];
        tensor_results(idx).Orientation = orientation;
        tensor_results(idx).Elongation = elongation;
        tensor_results(idx).Distance = distance;
        tensor_results(idx).Misalignment = misalignment_angle;
        tensor_results(idx).CenterVectorNorm = center_vector_norm;

    end
    
    dropData(i).tensor_results = tensor_results;
end
 
%% Mostrar las orientaciones en un subplot horizontal
figure('Name', 'Orientaciones y Vectores a Centro', 'NumberTitle', 'off', 'Units', 'normalized', 'Position', [0, 0.2, 1, 0.5]);
for i = 1:length(dropData)
    subplot(1, length(dropData), i);
    imshow(dropData(i).original_image, []);
    hold on;
    
    % Dibujar los vectores de orientación y hacia el centro
    for idx = 1:length(dropData(i).tensor_results)
        res = dropData(i).tensor_results(idx);
        x_center = res.BBox(1) + res.BBox(3) / 2;
        y_center = res.BBox(2) + res.BBox(4) / 2;
        
        quiver(x_center, y_center, cos(res.Orientation) * 10, sin(res.Orientation) * 10, 'r', 'LineWidth', 1.5);
        quiver(x_center, y_center, res.CenterVectorNorm(1) * 10, res.CenterVectorNorm(2) * 10, 'b', 'LineWidth', 1.5);
    end
    
    hold off;
    title(['Orientaciones - ', dropData(i).name], 'FontSize', 12);
    axis image;
end

%% Color Crystals Based on Misalignment Angle Bins
figure('Name', 'Crystals Colored by Misalignment Angle Bins', 'NumberTitle', 'off', 'Units', 'normalized', 'Position', [0, 0.2, 1, 0.5]);

for i = 1:length(dropData)
    % Step 1: Define misalignment angle bins
    bin_edges = 0:1/16:1;
    num_bins = length(bin_edges) - 1;
    colormap_discrete = parula(num_bins);
    
    % Step 2: Label the segmentation image
    label_image = bwlabel(dropData(i).combined_segmentation > 0);
    num_labels = max(label_image(:));
    
    % Step 3: Assign colors based on misalignment bins
    custom_colors = zeros(num_labels, 3);
    misalignments = [dropData(i).tensor_results.Misalignment];
    
    for j = 1:num_labels
        if j <= length(misalignments) && ~isnan(misalignments(j))
            bin_idx = find(misalignments(j) >= bin_edges, 1, 'last');
            if bin_idx > num_bins
                bin_idx = num_bins;
            end
            custom_colors(j, :) = colormap_discrete(bin_idx, :);
        else
            custom_colors(j, :) = [0, 0, 0];
        end
    end
    
    % Step 4: Apply custom colors to segmentation image
    colored_masks_rgb = label2rgb(label_image, custom_colors, 'k', 'shuffle');
    
    % Step 5: Display subplot for each drop
    subplot(1, length(dropData), i);
    imshow(colored_masks_rgb);
    title(['Misalignment - ', dropData(i).name], 'FontSize', 12);
    axis image;
end

% Step 6: Create a colorbar with bin labels
colormap(colormap_discrete);
c = colorbar;
ylabel(c, 'Misalignment Angle (degrees)');

%% test

%% Color Crystals by Misalignment (Only if Elongation ≥ Threshold)
elongation_threshold = 1.01;

figure('Name', 'Filtered Misalignment (Only Elongated Crystals)', 'NumberTitle', 'off', 'Units', 'normalized', 'Position', [0, 0.2, 1, 0.5]);

for i = 1:length(dropData)
    % Step 1: Define misalignment bins again
    bin_edges = 0:1/16:1;
    num_bins = length(bin_edges) - 1;
    colormap_discrete = parula(num_bins);
    
    % Step 2: Label segmentation
    label_image = bwlabel(dropData(i).combined_segmentation > 0);
    num_labels = max(label_image(:));
    
    % Step 3: Get misalignment + elongation
    misalignments = [dropData(i).tensor_results.Misalignment];
    elongations = [dropData(i).tensor_results.Elongation];
    
    % Step 4: Assign custom color map (gray if elongation < threshold)
    custom_colors = zeros(num_labels, 3);
    gray_color = [0.5, 0.5, 0.5];
    
    for j = 1:num_labels
        if j <= length(misalignments) && ~isnan(misalignments(j))
            if elongations(j) >= elongation_threshold
                bin_idx = find(misalignments(j) >= bin_edges, 1, 'last');
                if bin_idx > num_bins
                    bin_idx = num_bins;
                end
                custom_colors(j, :) = colormap_discrete(bin_idx, :);
            else
                custom_colors(j, :) = gray_color;
            end
        else
            custom_colors(j, :) = [0, 0, 0];
        end
    end
    
    % Step 5: Convert label image to RGB with custom colors
    colored_masks_rgb = label2rgb(label_image, custom_colors, 'k', 'shuffle');
    
    % Step 6: Display
    subplot(1, length(dropData), i);
    imshow(colored_masks_rgb);
    title(['Elongation ≥ ', num2str(elongation_threshold), ' - ', dropData(i).name], 'FontSize', 12);
    axis image;
    
      % Step 7: Calculate and display percentage with misalignment < 0.5
    valid_indices = elongations >= elongation_threshold;
    misalignment_filtered = misalignments(valid_indices);
    num_valid = sum(valid_indices);
    num_below_05 = sum(misalignment_filtered < 0.5);

    if num_valid > 0
        percent_below_05 = 100 * num_below_05 / num_valid;
        annotation_text = sprintf('%.1f%% < 0.5', percent_below_05);
    else
        annotation_text = 'No valid crystals';
    end

    % Display in bottom-left corner of each subplot
    text(5, size(colored_masks_rgb,1)-10, annotation_text, ...
         'Color', 'w', 'FontSize', 10, 'FontWeight', 'bold', ...
         'BackgroundColor', 'k', 'Margin', 1);
end

% Step 7: Add shared colorbar
colormap(colormap_discrete);
c = colorbar;
ylabel(c, 'Misalignment Angle (filtered by elongation)');


%% Integrar datos de todas las gotas por región de elipse
elongation_groups = cell(1, num_ellipses);
orientation_groups = cell(1, num_ellipses);
misalignment_groups = cell(1, num_ellipses);
distance_groups = cell(1, num_ellipses);
elongation_threshold = 2.5;

for i = 1:length(dropData)
    for j = 1:num_ellipses
        % Extraer datos dentro de la región de elipse j
        indices = find([dropData(i).tensor_results.Distance] >= (j-1)/num_ellipses * dropData(i).first_mask_info.Radius & ...
                       [dropData(i).tensor_results.Distance] < j/num_ellipses * dropData(i).first_mask_info.Radius);
        
        if ~isempty(indices)
            valid_indices = indices([dropData(i).tensor_results(indices).Elongation] >= elongation_threshold);
            elongation_groups{j} = [elongation_groups{j}; reshape([dropData(i).tensor_results(indices).Elongation], [], 1)];
            orientation_groups{j} = [orientation_groups{j}; reshape([dropData(i).tensor_results(valid_indices).Orientation], [], 1)];
            misalignment_groups{j} = [misalignment_groups{j}; reshape([dropData(i).tensor_results(valid_indices).Misalignment], [], 1)];
            distance_groups{j} = [distance_groups{j}; reshape([dropData(i).tensor_results(indices).Distance], [], 1)];
        end
    end
end

%% Crear los boxplots combinados por región de elipse
figure('Name', 'Integrated Crystal Metrics Boxplots by Ellipse Region', 'NumberTitle', 'off', 'Units', 'normalized', 'Position', [0, 0.2, 1, 0.8]);

% 1. Boxplot for Elongation vs Distance (by ellipse region)
valid_elongation_indices = find(~cellfun(@isempty, elongation_groups));
subplot(2, 2, 1);
boxplot(cell2mat(elongation_groups(valid_elongation_indices)'), repelem(valid_elongation_indices, cellfun(@numel, elongation_groups(valid_elongation_indices))));
xlabel('Ellipse Region');
ylabel('Elongation');
title('Elongation vs Ellipse Region');
grid on;

% 2. Boxplot for Orientation vs Distance (by ellipse region, filtered by elongation >= 3)
valid_orientation_indices = find(~cellfun(@isempty, orientation_groups));
subplot(2, 2, 2);
boxplot(cell2mat(orientation_groups(valid_orientation_indices)'), repelem(valid_orientation_indices, cellfun(@numel, orientation_groups(valid_orientation_indices))));
xlabel('Ellipse Region');
ylabel('Orientation (radians)');
title('Orientation vs Ellipse Region (Filtered)');
grid on;

% 3. Boxplot for Misalignment vs Distance (by ellipse region)
valid_misalignment_indices = find(~cellfun(@isempty, misalignment_groups));
subplot(2, 2, 3);
boxplot(cell2mat(misalignment_groups(valid_misalignment_indices)'), repelem(valid_misalignment_indices, cellfun(@numel, misalignment_groups(valid_misalignment_indices))));
xlabel('Ellipse Region');
ylabel('Misalignment (radians)');
title('Misalignment vs Ellipse Region (Filtered)');
grid on;

% 4. Boxplot for Misalignment vs Elongation (by ellipse region)
valid_misalignment_indices = find(~cellfun(@isempty, misalignment_groups));
subplot(2, 2, 4);
boxplot(cell2mat(misalignment_groups(valid_misalignment_indices)'), repelem(valid_misalignment_indices, cellfun(@numel, misalignment_groups(valid_misalignment_indices))));
xlabel('Elongation');
ylabel('Misalignment (radians)');
title('Misalignment vs Elongation (Filtered)');
grid on;

%% Voronoi Analysis for Each Drop

% figure('Name', 'Voronoi Diagrams for All Drops', 'NumberTitle', 'off', 'Units', 'normalized', 'Position', [0, 0.2, 1, 0.6]);
% 
% for i = 1:length(dropData)
%     subplot(1, length(dropData), i);
%     imshow(dropData(i).original_image, []);
%     hold on;
%     title(['Voronoi - ', dropData(i).name]);
%     
%     % Get centroids and ensure they are in double precision
%     centroids = dropData(i).centroids; 
%     centroids = centroids(~any(isnan(centroids), 2), :);
%     
%     % Remove NaN values
%     centroids = centroids(~any(isnan(centroids), 2), :);
%     
%     % Compute Voronoi if enough centroids exist
%     if size(centroids, 1) > 2
%         [v, c] = voronoin(centroids);
%     else
%         warning('Not enough centroids to compute Voronoi diagram for %s', dropData(i).name);
%         continue;
%     end
%     
%     % Extract mask boundary
%     mask_boundary = bwboundaries(dropData(i).drop.segmentation);
%     mask_boundary = mask_boundary{1};
%     mask_x = mask_boundary(:, 2);
%     mask_y = mask_boundary(:, 1);
%     
%     % Initialize storage
%     dropData(i).voronoi_polys = cell(1, length(c));
%     dropData(i).voronoi_areas = NaN(1, length(c));
%     dropData(i).crystal_fractions = NaN(1, length(c));
%     
%     % Process Voronoi cells
%     for j = 1:length(c)
%         if all(c{j} > 0)
%             v1 = v(c{j}, 1);
%             v2 = v(c{j}, 2);
%             
%             % Clip Voronoi cells to mask
%             cell_x = v1;
%             cell_y = v2;
%             inside = inpolygon(cell_x, cell_y, mask_x, mask_y);
%             cell_x = cell_x(inside);
%             cell_y = cell_y(inside);
%             
%             if ~isempty(cell_x)
%                 dropData(i).voronoi_polys{j} = [cell_x, cell_y];
%                 dropData(i).voronoi_areas(j) = polyarea(cell_x, cell_y);
%                 
%                 % Compute crystal fraction
%                 mask = poly2mask(cell_x, cell_y, size(dropData(i).combined_segmentation, 1), size(dropData(i).combined_segmentation, 2));
%                 total_pixels = sum(mask(:));
%                 crystal_pixels = sum(sum(mask & dropData(i).combined_segmentation));
%                 dropData(i).crystal_fractions(j) = crystal_pixels / total_pixels;
%                 
%                 % Plot the Voronoi cell
%                 plot(cell_x, cell_y, 'w-', 'LineWidth', 1.2);
%             end
%         end
%     end
%     
%     hold off;
% end
% 
% colorbar;
% sgtitle('Voronoi Diagrams for All Drops. first option');

%% Voronoi Analysis for Each Drop (guardar áreas, polígonos y fracciones)
figure('Name','Voronoi Diagrams for All Drops','NumberTitle','off',...
       'Units','normalized','Position',[0,0.2,1,0.6]);

for i = 1:length(dropData)
    subplot(1,length(dropData),i);
    imshow(dropData(i).original_image,[]);
    hold on;
    title(['Voronoi – ', dropData(i).name]);

    % --- 1) Obtener centroides válidos ---
    centroids = dropData(i).centroids;
    centroids = centroids(~any(isnan(centroids),2),:);
    if size(centroids,1) < 3
        warning('No hay suficientes centroides en %s', dropData(i).name);
        hold off; continue;
    end

    % --- 2) Calcular Voronoi ---
    [V,C] = voronoin(centroids);

    % --- 3) Contorno de la gota ---
    maskBW = dropData(i).drop.segmentation>0;
    [hMask,wMask] = size(maskBW);
    B = bwboundaries(maskBW);
    maskB = B{1};
    maskX = maskB(:,2);
    maskY = maskB(:,1);

    % --- Inicializar arrays ---
    nCells = numel(C);
    dropData(i).voronoi_areas     = nan(1, nCells);
    dropData(i).voronoi_polys     = cell(1, nCells);
    dropData(i).crystal_fractions = nan(1, nCells);

    % --- 4) Recortar, dibujar y guardar datos ---
    for j = 1:nCells
        idx = C{j};
        if any(idx==1), continue; end  % polígono infinito

        vx = V(idx,1); vy = V(idx,2);
        inM = inpolygon(vx,vy,maskX,maskY);
        vx = vx(inM); vy = vy(inM);
        if isempty(vx), continue; end

        % cerrar el polígono
        vx(end+1)=vx(1); vy(end+1)=vy(1);

        % guardar polígono y área
        dropData(i).voronoi_polys{j} = [vx, vy];
        A = polyarea(vx,vy);
        dropData(i).voronoi_areas(j) = A;

        % calcular fracción cristalina
        cmask = poly2mask(vx,vy,hMask,wMask);
        totalPx   = sum(cmask(:));
        crystalPx = sum(cmask(:) & dropData(i).combined_segmentation(:));
        dropData(i).crystal_fractions(j) = crystalPx/totalPx;

        % trazar celda
        patch(vx, vy, [0.8 0.8 0.8], ...
              'FaceAlpha',0.5,'EdgeColor','w','LineWidth',1);
    end

    % --- 5) Dibujar centroides ---
    plot(centroids(:,1), centroids(:,2), 'wo', ...
         'MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',4);

    hold off;
end

sgtitle('Voronoi Diagrams for All Drops');


%% otra opcion de voronoi. corta con borde
% 
% %% Voronoi Analysis for Each Drop (filtrado + recorte manual)
% figure('Name', 'Voronoi Diagrams for All Drops (Solo interiores)', ...
%        'NumberTitle', 'off', 'Units', 'normalized', 'Position', [0, 0.2, 1, 0.6]);
% 
% for i = 1:length(dropData)
%     subplot(1, length(dropData), i);
%     imshow(dropData(i).original_image, []);
%     hold on;
%     title(['Voronoi - ', dropData(i).name]);
% 
%     % --- 1. Filtrar solo centroides dentro de la máscara de la gota ---
%     maskBW    = dropData(i).drop.segmentation > 0;
%     allC      = dropData(i).centroids;
%     inside    = false(size(allC,1),1);
%     [hMask,wMask] = size(maskBW);
%     for k = 1:size(allC,1)
%         x = round(allC(k,1));
%         y = round(allC(k,2));
%         if x>=1 && x<=wMask && y>=1 && y<=hMask && maskBW(y,x)
%             inside(k) = true;
%         end
%     end
%     centroids = allC(inside,:);
% 
%     if size(centroids,1) < 3
%         warning('Menos de 3 centroides válidos en %s → omitiendo.', dropData(i).name);
%         continue;
%     end
% 
%     % --- 2. Calcular Voronoi sobre esos centroides ---
%     [V, C] = voronoin(centroids);
% 
%     % --- 3. Obtener contorno de la gota para recorte ---
%     B     = bwboundaries(maskBW);
%     maskB = B{1};
%     maskX = maskB(:,2);
%     maskY = maskB(:,1);
% 
%     % Inicializar almacenamiento
%     nCells = numel(C);
%     dropData(i).voronoi_polys     = cell(1, nCells);
%     dropData(i).voronoi_areas     = nan(1, nCells);
%     dropData(i).crystal_fractions = nan(1, nCells);
% 
%     % --- 4. Recorte manual por aristas ---
%     for j = 1:nCells
%         idx = C{j};
%         if any(idx==1), continue; end   % polígono infinito
%         xv = V(idx,1);
%         yv = V(idx,2);
% 
%         x_clip = [];
%         y_clip = [];
% 
%         for m = 1:length(xv)
%             x1 = xv(m);  y1 = yv(m);
%             nn = mod(m,length(xv))+1;
%             x2 = xv(nn); y2 = yv(nn);
% 
%             in1 = inpolygon(x1,y1,maskX,maskY);
%             in2 = inpolygon(x2,y2,maskX,maskY);
%             if in1
%                 x_clip(end+1,1) = x1;  %#ok<*AGROW>
%                 y_clip(end+1,1) = y1;
%             end
%             if in1 ~= in2
%                 % buscar intersecciones con cada segmento del contorno
%                 for t = 1:size(maskB,1)-1
%                     bx1 = maskX(t);   by1 = maskY(t);
%                     bx2 = maskX(t+1); by2 = maskY(t+1);
%                     denom = (x1-x2)*(by1-by2) - (y1-y2)*(bx1-bx2);
%                     if denom==0, continue; end
%                     xi = ((x1*y2 - y1*x2)*(bx1-bx2) - (x1-x2)*(bx1*by2 - by1*bx2)) / denom;
%                     yi = ((x1*y2 - y1*x2)*(by1-by2) - (y1-y2)*(bx1*by2 - by1*bx2)) / denom;
%                     if xi>=min(x1,x2)&&xi<=max(x1,x2) && ...
%                        xi>=min(bx1,bx2)&&xi<=max(bx1,bx2) && ...
%                        yi>=min(y1,y2)&&yi<=max(y1,y2) && ...
%                        yi>=min(by1,by2)&&yi<=max(by1,by2)
%                         x_clip(end+1,1) = xi;
%                         y_clip(end+1,1) = yi;
%                     end
%                 end
%             end
%         end
% 
%         if isempty(x_clip)
%             continue;
%         end
% 
%         % cerrar polígono
%         x_clip(end+1) = x_clip(1);
%         y_clip(end+1) = y_clip(1);
% 
%         % almacenar polígonos y calcular métricas
%         dropData(i).voronoi_polys{j}     = [x_clip, y_clip];
%         dropData(i).voronoi_areas(j)     = polyarea(x_clip,y_clip);
%         cmask = poly2mask(x_clip,y_clip,hMask,wMask);
%         totalPx   = sum(cmask(:));
%         crystalPx = sum(cmask(:) & dropData(i).combined_segmentation(:));
%         dropData(i).crystal_fractions(j)= crystalPx/totalPx;
% 
%         % trazar celda
%         patch(x_clip, y_clip, 'w', 'FaceAlpha', 0.5, 'EdgeColor', 'k', 'LineWidth', 1);
%     end
% 
%     % --- 5. Dibujar centroides filtrados ---
%     plot(centroids(:,1), centroids(:,2), 'wo', ...
%          'MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',5);
% 
%     hold off;
% end


%% First Figure: Voronoi Cells Colored by Area

% figure('Name', 'Voronoi Cells Colored by Area', 'NumberTitle', 'off', 'Units', 'normalized', 'Position', [0, 0.2, 1, 0.6]);
% for i = 1:length(dropData)
%     subplot(1, length(dropData), i);
%     imshow(dropData(i).original_image, []);
%     hold on;
%     cmap = jet(256);
%     min_area = min(dropData(i).voronoi_areas);
%     max_area = max(dropData(i).voronoi_areas);
%     for j = 1:length(dropData(i).voronoi_polys)
%         if ~isempty(dropData(i).voronoi_polys{j})
%             polyCoords = dropData(i).voronoi_polys{j};
%             x_clipped = polyCoords(:,1);
%             y_clipped = polyCoords(:,2);
%             norm_area = (dropData(i).voronoi_areas(j) - min_area) / (max_area - min_area);
%             color_idx = max(1, min(256, round(norm_area * 255) + 1));
%             cell_color = cmap(color_idx, :);
%             patch(x_clipped, y_clipped, cell_color, 'FaceAlpha', 0.5, 'EdgeColor', 'k', 'LineWidth', 1.5);
%             colormap(jet(256));    % Aplicar el colormap jet
%         end
%     end
%     colorbar;
%     title(['Voronoi - ', dropData(i).name, ' (Area)']);
%     axis image;
%     hold off;
% end


%% First Figure: Voronoi Cells Colored by Area (estilo unificado)
figure('Name','Voronoi Cells Colored by Area','NumberTitle','off', ...
       'Units','normalized','Position',[0,0.2,1,0.6]);

cmap_area = jet(256);

for i = 1:length(dropData)
    subplot(1,length(dropData),i);
    % fondo blanco
    imshow(ones([size(dropData(i).combined_segmentation),3]));
    hold on;
    
    % calcular rango de áreas normalizado
    areas = dropData(i).voronoi_areas;
    minA = min(areas(~isnan(areas)));
    maxA = max(areas(~isnan(areas)));
    
    for j = 1:length(dropData(i).voronoi_polys)
        poly = dropData(i).voronoi_polys{j};
        if isempty(poly), continue; end
        
        x = poly(:,1);
        y = poly(:,2);
        a = dropData(i).voronoi_areas(j);
        % índice de color
        idx = max(1, min(256, round((a-minA)/(maxA-minA)*255)+1));
        color = cmap_area(idx,:);
        
        patch(x, y, color, ...
              'EdgeColor','k','LineWidth',1.2, ...
              'FaceAlpha',0.8);
    end
    
    % contorno de cristales reales
    B = bwboundaries(dropData(i).combined_segmentation>0);
    for k = 1:numel(B)
        b = B{k};
        plot(b(:,2), b(:,1), 'k', 'LineWidth',1);
    end
    
    axis image off;
    title(['Voronoi – ', dropData(i).name,' (Area)']);
    hold off;
end

colormap(cmap_area);
caxis([minA maxA]);
colorbar('Position',[0.93,0.2,0.015,0.6]);
sgtitle('Voronoi Cells Colored by Area');


%% Second Figure: Voronoi Cells Colored by Crystal Fraction

% figure('Name', 'Voronoi Cells Colored by Crystal Fraction', 'NumberTitle', 'off', 'Units', 'normalized', 'Position', [0, 0.2, 1, 0.6]);
% for i = 1:length(dropData)
%     subplot(1, length(dropData), i);
%     imshow(~dropData(i).combined_segmentation, []);
%     hold on;
%     for j = 1:length(dropData(i).voronoi_polys)
%         if ~isempty(dropData(i).voronoi_polys{j})
%             polyCoords = dropData(i).voronoi_polys{j};
%             x_clipped = polyCoords(:,1);
%             y_clipped = polyCoords(:,2);
%             color_idx = max(1, min(256, round(dropData(i).crystal_fractions(j) * 255) + 1));
%             cell_color = jet(256);
%             patch(x_clipped, y_clipped, cell_color(color_idx, :), 'FaceAlpha', 0.5, 'EdgeColor', 'k', 'LineWidth', 1.5);
%             colormap(jet(256));    % Aplicar el colormap jet
%             caxis([0 1]);          % Establecer el rango de la fracción cristalina
%         end
%     end
%     colorbar;
%     title(['Voronoi - ', dropData(i).name, ' (Crystal Fraction)']);
%     axis image;
%     hold off;
% end

figure('Name', 'Voronoi Cells Colored by Crystal Fraction', 'NumberTitle', 'off', ...
    'Units', 'normalized', 'Position', [0, 0.2, 1, 0.6]);

cmap = parula(256);

for i = 1:length(dropData)
    subplot(1, length(dropData), i);
    imshow(ones([size(dropData(i).combined_segmentation), 3]));
    hold on;

    for j = 1:length(dropData(i).voronoi_polys)
        if ~isempty(dropData(i).voronoi_polys{j})
            polyCoords = dropData(i).voronoi_polys{j};
            x_clipped = polyCoords(:,1);
            y_clipped = polyCoords(:,2);

            frac = dropData(i).crystal_fractions(j);
            if isnan(frac)
                color = [0.8 0.8 0.8];
            else
                color_idx = max(1, min(256, round(frac * 255) + 1));
                color = cmap(color_idx, :);
            end

            patch(x_clipped, y_clipped, color, ...
                  'EdgeColor', 'k', 'LineWidth', 1.2, 'FaceAlpha', 0.8);
        end
    end

    % Añadir contorno de cristales reales
    boundaries = bwboundaries(dropData(i).combined_segmentation > 0);
    for k = 1:length(boundaries)
        b = boundaries{k};
        plot(b(:,2), b(:,1), 'k', 'LineWidth', 1.0);
    end

    axis image;
    title(['Voronoi - ', dropData(i).name]);
end

colormap(cmap);
caxis([0 1]);
colorbar('Position', [0.93, 0.2, 0.015, 0.6]);
sgtitle('Voronoi Cells Colored by Crystal Fraction');


%% Third Figure: Crystal Fraction vs Distance to Center for Each Drop
figure('Name', 'Crystal Fraction vs Distance to Center', 'NumberTitle', 'off', 'Units', 'normalized', 'Position', [0, 0.8, 1, 0.6]);
for i = 1:length(dropData)
    subplot(1, length(dropData), i);
    distances = sqrt((double(dropData(i).bboxes(:,1)) - dropData(i).first_mask_info.Center(1)).^2 + (double(dropData(i).bboxes(:,2)) - dropData(i).first_mask_info.Center(2)).^2);
    crystal_fractions = dropData(i).crystal_fractions(:);
    scatter(distances, crystal_fractions, 50, crystal_fractions, 'filled');
    xlabel('Distance from Center');
    ylabel('Crystal Fraction');
    title(['Crystal Fraction vs Distance - ', dropData(i).name]);
    grid on;
    colorbar;
end

%% Line Graph: Voronoi Area vs Normalized Radius (R / Rc)

figure('Name', 'Voronoi Area vs R / Rc', 'NumberTitle', 'off', 'Units', 'normalized', 'Position', [0, 0.8, 1, 0.6]);

% Número de zonas radiales
num_ellipses = 10;
R_over_Rc = linspace(0.1, 1.0, num_ellipses);  % de borde (0.1) a centro (1.0)

% Subplot 1: Individual Lines for Each Drop
subplot(1,2,1);
hold on;
all_voronoi_areas = zeros(num_ellipses, length(dropData));

for i = 1:length(dropData)
    radius = dropData(i).first_mask_info.Radius;
    ellipse_bins = linspace(0, radius, num_ellipses + 1);  % bins por distancia radial

    mean_voronoi_areas = zeros(1, num_ellipses);
    
    for j = 1:num_ellipses
        % Calcular distancia de cada celda al centro
        centroids = dropData(i).centroids;  % Asegúrate de que esto exista
        dists = sqrt((centroids(:,1) - dropData(i).first_mask_info.Center(1)).^2 + ...
                     (centroids(:,2) - dropData(i).first_mask_info.Center(2)).^2);

        in_bin = dists >= ellipse_bins(j) & dists < ellipse_bins(j+1);
        mean_voronoi_areas(j) = mean(dropData(i).voronoi_areas(in_bin), 'omitnan');
    end
    all_voronoi_areas(:,i) = mean_voronoi_areas;
    plot(R_over_Rc, mean_voronoi_areas, 'LineWidth', 1.5, 'DisplayName', dropData(i).name);
end
xlabel('R / R_c');
ylabel('Mean Voronoi Area');
title('Voronoi Area vs R / R_c (Each Drop)');
legend;
grid on;
hold off;

% Subplot 2: Mean, Max, Min across Drops
subplot(1,2,2);
mean_voronoi = mean(all_voronoi_areas, 2, 'omitnan');
max_voronoi = max(all_voronoi_areas, [], 2);
min_voronoi = min(all_voronoi_areas, [], 2);
plot(R_over_Rc, mean_voronoi, 'k-', 'LineWidth', 2); hold on;
plot(R_over_Rc, max_voronoi, 'r--', 'LineWidth', 1.5);
plot(R_over_Rc, min_voronoi, 'b--', 'LineWidth', 1.5);
xlabel('R / R_c');
ylabel('Voronoi Area');
title('Mean, Max, Min Voronoi Area');
legend('Mean', 'Max', 'Min');
grid on;
hold off;

%% NEW: %% Line Graph: Normalized Voronoi Area vs R / R_c
figure('Name', 'Normalized Voronoi Area vs R / R_c', 'NumberTitle', 'off', ...
       'Units', 'normalized', 'Position', [0, 0.8, 1, 0.6]);

num_ellipses = 10;
R_over_Rc    = linspace(0.1, 1.0, num_ellipses);

% Subplot 1: Each Drop (normalized)
subplot(1,2,1);
hold on;
all_voronoi_norm = zeros(num_ellipses, length(dropData));

for i = 1:length(dropData)
    radius = dropData(i).first_mask_info.Radius;
    area   = dropData(i).first_mask_info.Area;       % droplet area for normalization
    bins   = linspace(0, radius, num_ellipses + 1);

    mean_areas = zeros(1, num_ellipses);
    for j = 1:num_ellipses
        centroids = dropData(i).centroids;
        dists     = sqrt((centroids(:,1) - dropData(i).first_mask_info.Center(1)).^2 + ...
                         (centroids(:,2) - dropData(i).first_mask_info.Center(2)).^2);
        in_bin    = dists >= bins(j) & dists < bins(j+1);
        mean_areas(j) = mean(dropData(i).voronoi_areas(in_bin), 'omitnan');
    end

    % --- normalization step ---
    mean_areas_norm = mean_areas / area;

    all_voronoi_norm(:,i) = mean_areas_norm;
    plot(R_over_Rc, mean_areas_norm, 'LineWidth', 1.5, 'DisplayName', dropData(i).name);
end
xlabel('R / R_c');
ylabel('Normalized Mean Voronoi Area');
title('Normalized Voronoi Area vs R / R_c (Each Drop)');
legend;
grid on;
hold off;

% Subplot 2: Mean, Max, Min (normalized)
subplot(1,2,2);
mean_norm = mean(all_voronoi_norm, 2, 'omitnan');
max_norm  = max(all_voronoi_norm, [], 2);
min_norm  = min(all_voronoi_norm, [], 2);

plot(R_over_Rc, mean_norm, 'k-', 'LineWidth', 2); hold on;
plot(R_over_Rc, max_norm,  'r--', 'LineWidth', 1.5);
plot(R_over_Rc, min_norm,  'b--', 'LineWidth', 1.5);
xlabel('R / R_c');
ylabel('Normalized Voronoi Area');
title('Mean, Max, Min Normalized Voronoi Area');
legend('Mean', 'Max', 'Min');
grid on;
hold off;


%% Line Graph: Crystal Fraction vs Normalized Radius (R / Rc)
figure('Name', 'Crystal Fraction vs R / R_c', 'NumberTitle', 'off', 'Units', 'normalized', 'Position', [0, 0.8, 1, 0.6]);

% Eje X: de borde (0.1) a centro (1.0)
num_ellipses = 10;
R_over_Rc = linspace(0.1, 1.0, num_ellipses);
all_crystal_fractions = NaN(num_ellipses, length(dropData));

% Subplot 1: líneas por gota
subplot(1,2,1);
hold on;

for i = 1:length(dropData)
    if ~isfield(dropData(i), 'crystal_fractions') || isempty(dropData(i).crystal_fractions)
        warning('No crystal fraction data for %s', dropData(i).name);
        continue;
    end

    % Calcular distancias de centroides al centro
    centroids = dropData(i).centroids;  % asegúrate de que esto esté en tu estructura
    distances = sqrt((centroids(:,1) - dropData(i).first_mask_info.Center(1)).^2 + ...
                     (centroids(:,2) - dropData(i).first_mask_info.Center(2)).^2);

    % Crear bins radiales en base al radio de la gota
    radius = dropData(i).first_mask_info.Radius;
    ellipse_bins = linspace(0, radius, num_ellipses + 1);
    mean_crystal_fractions = NaN(1, num_ellipses);

    for j = 1:num_ellipses
        in_bin = distances >= ellipse_bins(j) & distances < ellipse_bins(j+1);
        if any(in_bin)
            mean_crystal_fractions(j) = mean(dropData(i).crystal_fractions(in_bin), 'omitnan');
        end
    end

    % Rellenar valores faltantes para suavidad
    if all(isnan(mean_crystal_fractions))
        warning('No valid crystal fraction data found for %s', dropData(i).name);
        continue;
    end

    mean_crystal_fractions = fillmissing(mean_crystal_fractions, 'linear', 'EndValues', 'nearest');
    all_crystal_fractions(:, i) = mean_crystal_fractions;

    plot(R_over_Rc, mean_crystal_fractions, 'LineWidth', 1.5, 'DisplayName', dropData(i).name);
end
xlabel('R / R_c');
ylabel('Mean Crystal Fraction');
title('Crystal Fraction vs R / R_c (Each Drop)');
legend;
grid on;
hold off;

% Subplot 2: Media, Máximo, Mínimo globales
subplot(1,2,2);
mean_fraction = mean(all_crystal_fractions, 2, 'omitnan');
max_fraction = max(all_crystal_fractions, [], 2, 'omitnan');
min_fraction = min(all_crystal_fractions, [], 2, 'omitnan');
plot(R_over_Rc, mean_fraction, 'k-', 'LineWidth', 2); hold on;
plot(R_over_Rc, max_fraction, 'r--', 'LineWidth', 1.5);
plot(R_over_Rc, min_fraction, 'b--', 'LineWidth', 1.5);
xlabel('R / R_c');
ylabel('Crystal Fraction');
title('Mean, Max, Min Crystal Fraction');
legend('Mean', 'Max', 'Min');
grid on;
hold off;



%% Save all figures

outputName = 'results CNB-5drops';  % folder name here

% Base output folder (shared location)
baseOutputFolder = '/Users/anabellehassan/Documents/Biomedical eng/TFG/Segment/';
fullOutputPath = fullfile(baseOutputFolder, outputName);

% Create folder if it doesn’t exist
if ~exist(fullOutputPath, 'dir')
    mkdir(fullOutputPath);
end

% Get all figure handles
figHandles = findall(0, 'Type', 'figure');

% Save each figure
for f = 1:length(figHandles)
    fig = figHandles(f);
    figName = get(fig, 'Name');
    if isempty(figName)
        figName = sprintf('figure_%02d', f);
    end
    safeName = matlab.lang.makeValidName(figName);  % avoid invalid filename characters
    savePath = fullfile(fullOutputPath, [safeName, '.png']);
    exportgraphics(fig, savePath, 'Resolution', 300);
end

fprintf('Saved all figures in: %s\n', fullOutputPath);

