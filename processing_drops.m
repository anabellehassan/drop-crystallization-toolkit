%% processing_drops.m - Load, clean, analyze, and save structured data for each droplet.
clear; clc; close all;

%% Load folders, Configuration
folderList = {'C1-2_Dtest6', 'C1-3_FTest20', 'C1-4_Dtest21'}; 

%C1-2_Dtest6
%C1-3_FTest20
%C1-4_Dtest21
%CNB_040425_02_bf_20250510_2154
%CNB_040425_03_20250511_1148
%CNB_040425_05_bf_20250511_1330
%CNB_040425_04_20250511_1307
%CNB_040425_06_20250511_1518

% Base folder where all droplet folders are stored
baseFolder = '/Users/anabellehassan/Documents/Biomedical eng/TFG/Segment/';

% Output folder for saving dropData and tables, with timestamp
timestamp = datestr(now, 'ddmmyy_HHMM');  % Format: 06040625_1027
outputFolder = fullfile(baseFolder, ['processed_output_', timestamp]);
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

%% General loop over droplets

for i = 1:length(folderList)
  % Load files
    dropName = folderList{i};
    dropPath = fullfile(baseFolder, dropName);
    fprintf('Loading: %s\n', dropName);
    
    if ~isfolder(dropPath)
        warning('Folder not found: %s — skipping.', dropName);
        continue;
    end
    
    try
        data = load(fullfile(dropPath, 'combined_segmentation.mat'));
        bboxes = load(fullfile(dropPath, 'bboxes.mat'));
        areas = load(fullfile(dropPath, 'areas.mat'));
        drop = load(fullfile(dropPath, 'first_mask_drop.mat'));
    catch
        warning('Missing required files in %s — skipping.', dropName);
        continue;
    end
    
    % Load .tif image
    tifName = [dropName, '.tif'];
    tifPath = fullfile(dropPath, tifName);
    if exist(tifPath, 'file')
        original_image = imread(tifPath);
    else
        warning('Missing .tif for %s — using zero image.', dropName);
        original_image = zeros(512); 
    end

  % Process information at drop level
  
    % Clean bboxes and areas (in case there are duplicate bboxes)
    [unique_bboxes, unique_idx] = unique(bboxes.bboxes, 'rows');
    clean_areas = areas.areas(unique_idx);

    if length(unique_idx) < size(bboxes.bboxes, 1)
        fprintf('Removed %d duplicate bounding boxes\n', ...
                size(bboxes.bboxes,1) - length(unique_idx));
    end

    % Compute drop general geometry from first mask (dropMask)
    dropMask = drop.segmentation;
    dropArea = double(drop.area);
    dropRadius = sqrt(dropArea / pi);
    [Y, X] = find(dropMask);
    dropCenter = [mean(X), mean(Y)];

    % Save all
    dropInfo.name = dropName;
    dropInfo.folder = dropPath;
    dropInfo.combined_segmentation = single(data.combined_segmentation);
    dropInfo.original_image = original_image;
    dropInfo.drop = drop;
    dropInfo.Area = dropArea;
    dropInfo.Radius = dropRadius;
    dropInfo.Center = dropCenter;
    dropInfo.numCrystals = size(unique_bboxes, 1);
    dropInfo.crystal_bboxes_list = unique_bboxes;
    dropInfo.crystal_areas_list = clean_areas;
     
    dropData(i) = dropInfo;
end

fprintf('\nFinished loading and preprocessing all droplets.\n');

% To add later:
% dropInfo.meanArea = NaN;
% dropInfo.meanPerimeter = NaN;

%% Ellipse Region Generation Using Drop Mask Geometry

numEllipses = 10; 
theta = linspace(0, 2*pi, 100);  % ellipse perimeter resolution

for i = 1:length(dropData)
    % full mask to find ellipse
    mask = dropData(i).drop.segmentation > 0;
    props = regionprops(mask, 'MajorAxisLength', 'MinorAxisLength', ...
                               'Orientation', 'Centroid');

    if isempty(props)
        warning('Could not compute ellipse for drop %s — skipping.', dropData(i).name);
        dropData(i).ellipseInfo = [];
        continue;
    end

    % Extract shape properties
    MA  = props(1).MajorAxisLength;
    ma  = props(1).MinorAxisLength;
    ori = -deg2rad(props(1).Orientation);  % negative for image coord correction
    cen = props(1).Centroid;

    % Build ellipses from outermost to innermost
    ellipseStruct = struct('Points', []);
    for k = 1:numEllipses
        scale = (numEllipses - k + 1) / numEllipses;
        a = scale * MA / 2;
        b = scale * ma / 2;

        % Parametric ellipse points before rotation
        pts = [a * cos(theta); b * sin(theta)];

        % Rotate and translate
        R = [cos(ori), -sin(ori); sin(ori), cos(ori)];
        pts_rotated = R * pts + cen';

        % Store
        ellipseStruct(k).Points = pts_rotated;
    end

    % Save into dropData
    dropData(i).ellipseInfo = ellipseStruct;
end

fprintf('\nElliptical region masks computed and stored per drop.\n');

%% Visual Check — Ellipses on Original Images

% for i = 1:length(dropData)
%     img = dropData(i).original_image;
%     ellipses = dropData(i).ellipseInfo;
%     figure('Name', ['Ellipses - ', dropData(i).name], 'NumberTitle', 'off');
%     imshow(img, []); hold on;
%     for k = 1:length(ellipses)
%         P = ellipses(k).Points;
%         plot(P(1,:), P(2,:), '-', 'LineWidth', 1.5);
%     end
%     title(['Elliptical Regions: ', dropData(i).name], 'Interpreter', 'none');
%     hold off;
% end

%% Area/Size Calculations Normalization — with debug

for i = 1:length(dropData)

    areas = dropData(i).crystal_areas_list;
    dropArea = sum(dropMask(:));

    % Normalize areas (correct: by area of the drop)
    normalizedAreas = double(areas) / dropArea;

    % Store in drop level
    dropData(i).meanArea = mean(areas); % pixels
    dropData(i).meanArea_norm = mean(normalizedAreas); % adimensional

    if ~isfield(dropData(i), 'crystals') || isempty(dropData(i).crystals)
        dropData(i).crystals = repmat(struct('normalized_area', []), length(areas), 1);
    end

    for j = 1:length(areas)
        dropData(i).crystals(j).normalized_area = normalizedAreas(j);
    end
end

fprintf('\n✅ Normalized areas computed and stored.\n');
disp([dropData.meanArea_norm])


%% Crystal Area Density per Elliptical Region

for i = 1:length(dropData)
    drop = dropData(i);

    nRegions = length(drop.ellipseInfo) - 1;
    seg = drop.combined_segmentation > 0;
    [H, W] = size(seg);

    crystal_area = zeros(1, nRegions);
    region_area  = zeros(1, nRegions);

    for k = 1:nRegions
        % Define ring mask between ellipses
        P_outer = drop.ellipseInfo(k).Points;
        P_inner = drop.ellipseInfo(k+1).Points;

        m_outer = poly2mask(P_outer(1,:), P_outer(2,:), H, W);
        m_inner = poly2mask(P_inner(1,:), P_inner(2,:), H, W);

        ring_mask = m_outer & ~m_inner;

        % Count pixels
        region_area(k)  = sum(ring_mask(:));
        crystal_area(k) = sum(seg(ring_mask));
    end

    % Final center region (innermost ellipse)
    P_last = drop.ellipseInfo(end).Points;
    m_last = poly2mask(P_last(1,:), P_last(2,:), H, W);
    region_area(end+1)  = sum(m_last(:));
    crystal_area(end+1) = sum(seg(m_last));

    % Area density = fraction of pixels in region that belong to crystals
    area_density = crystal_area ./ region_area;

    % Store
    dropData(i).areaDensityProfile       = area_density; % vector [1 × N], one value per ring
    % dropData(i).areaDensityRegionAreas   = region_area; % total area (pixels) per ring
    % dropData(i).areaDensityCrystalAreas  = crystal_area; % crystal area (pixels) per ring
end

fprintf('\nArea density per elliptical region computed and stored.\n');

%% Tensor Metrics + Shape Descriptors (per crystal)

for i = 1:length(dropData)
    drop = dropData(i);
    segMask = drop.combined_segmentation;
    center = drop.Center(:);  % column vector
    R0 = drop.Radius;
    img = double(drop.original_image);

    nCrystals = size(drop.crystal_bboxes_list, 1);
    updatedCrystals = drop.crystals;

    for j = 1:nCrystals
        % Bounding box [x, y, w, h]
        bb = round(drop.crystal_bboxes_list(j, :));
        x1 = max(1, bb(1));
        y1 = max(1, bb(2));
        x2 = min(size(segMask,2), x1 + bb(3) - 1);
        y2 = min(size(segMask,1), y1 + bb(4) - 1);

        % Extract patch and mask, convert everything to double
        rawImgPatch = double(img(y1:y2, x1:x2));
        %localMask = (segMask(y1:y2, x1:x2) == j);
        
        % Get cropped segmentation patch
        segCrop = segMask(y1:y2, x1:x2);
        % Find the most common non-zero label (assumes one dominant crystal in bbox)
        labels = segCrop(segCrop > 0);
        if isempty(labels)
            warning('Crystal %d in drop %s has no valid label — skipping.', j, drop.name);
            continue;
        end
        crystalLabel = mode(labels);  % dominant label in the region

        % Create local mask for this crystal
        localMask = double(segCrop == crystalLabel);
        % Multiply and store final image patch (only crystal j)
        localImg = rawImgPatch .* localMask;
        
        % Skip if mask is empty or total intensity is zero
        if all(localImg(:) == 0) || sum(localMask(:)) == 0
            warning('Crystal %d in drop %s has empty mask — skipping.', j, drop.name);
            continue;
        end

        % Compute intensity-weighted centroid
        [W, H] = meshgrid(x1:x2, y1:y2);
        W = double(W);
        H = double(H);
        S  = sum(localImg(:));
        xc = sum(localImg(:) .* W(:)) / S;
        yc = sum(localImg(:) .* H(:)) / S;

        dx = W(:) - xc;
        dy = H(:) - yc;

        % Inertia tensor
        Txx = sum((dx.^2) .* localImg(:)) / S;
        Tyy = sum((dy.^2) .* localImg(:)) / S;
        Txy = sum((dx .* dy) .* localImg(:)) / S;
        T = [Txx, Txy; Txy, Tyy];

        % Eigen decomposition
        [V,D] = eig(T);
        lambdas = diag(D); % eigenvalues
        [lambda_max, mi] = max(lambdas); % max eigenvalue (major axis)
        lambda_min = min(lambdas); % min eigenvalue (minor axis)
        v_max = V(:,mi); % eigenvector of major axis
        if v_max(2)<0, v_max = -v_max; end % make y-positive for consistency

        % Metrics
        elongation = lambda_max / max(lambda_min, 1e-6);  % avoid /0, how stretched
        orientation = atan2(v_max(2), v_max(1));

        % Radial vector
        vec = center - [xc; yc];
        dist = norm(vec); % distance in pixels
        r_norm = dist / R0; % normalized to drop radius
        
        % Misalignment = sin(θ) between crystal and radial direction
        misalignment = abs(v_max(1)*vec(2) - v_max(2)*vec(1)) / (norm(v_max)*norm(vec)); 

        % Store in crystal struct
        updatedCrystals(j).elongation = elongation;
        updatedCrystals(j).orientationAngle = orientation;
        updatedCrystals(j).misalignment = misalignment;
        updatedCrystals(j).distanceToCenter = dist;
        updatedCrystals(j).r_normalized = r_norm;
        updatedCrystals(j).centroid = [xc, yc];

    end

    % Save back to drop
    dropData(i).crystals = updatedCrystals;
    
    % Extract all valid elongation and misalignment values
    elongVals = [updatedCrystals.elongation];
    misalVals = [updatedCrystals.misalignment];

    % Remove NaNs (just in case)
    elongVals = elongVals(~isnan(elongVals));
    misalVals = misalVals(~isnan(misalVals));

    % Store mean metrics at drop level
    dropData(i).meanElongation = mean(elongVals);
    dropData(i).meanMisalignment = mean(misalVals);
end

fprintf('\nTensor and shape metrics computed and stored per crystal.\n');

%% Plot tensor metrics
% % Pick one drop to test (replace 1 with any valid index)
% drop = dropData(3);
% 
% % Open figure
% figure('Name',['Orientation Arrows — ', drop.name], 'NumberTitle','off');
% imshow(drop.original_image, []); hold on;
% title(['Crystal Orientation Arrows — ', drop.name]);
% 
% % Loop through crystals
% for j = 1:length(drop.crystals)
%     c = drop.crystals(j);
% 
%     % Skip if orientation was not computed
%     if ~isfield(c, 'orientationAngle') || isnan(c.orientationAngle)
%         continue;
%     end
% 
%     % Use centroid if available
%     if isfield(c, 'centroid') && all(~isnan(c.centroid))
%         x0 = c.centroid(1);
%         y0 = c.centroid(2);
%     else
%         % Fallback: use bounding box center
%         bb = drop.crystal_bboxes_list(j,:);
%         x0 = bb(1) + bb(3)/2;
%         y0 = bb(2) + bb(4)/2;
%     end
% 
%     % Orientation vector (arrow)
%     angle = c.orientationAngle;
%     len = 15;  % arrow length for visibility
%     dx = len * cos(angle);
%     dy = len * sin(angle);
% 
%     % Plot arrow
%     quiver(x0, y0, dx, dy, 0, 'r', 'LineWidth', 1);
% end
% 
% hold off;


%% Perimeter per crystal (from segmentation)

for i = 1:length(dropData)
    drop = dropData(i);
    segMask = drop.combined_segmentation;
    nCrystals = size(drop.crystal_bboxes_list, 1);
    img = double(drop.original_image);  % in case needed later

    updatedCrystals = drop.crystals;

    for j = 1:nCrystals
        % Extract bounding box
        bb = round(drop.crystal_bboxes_list(j, :));
        x1 = max(1, bb(1));
        y1 = max(1, bb(2));
        x2 = min(size(segMask,2), x1 + bb(3) - 1);
        y2 = min(size(segMask,1), y1 + bb(4) - 1);

        % Get cropped segmentation patch
        segCrop = segMask(y1:y2, x1:x2);
        labels = segCrop(segCrop > 0);
        if isempty(labels)
            updatedCrystals(j).perimeter = NaN;
            continue;
        end
        crystalLabel = mode(labels);

        % Create local mask
        localMask = (segCrop == crystalLabel);

        % Compute perimeter
        stats = regionprops(localMask, 'Perimeter');
        if ~isempty(stats)
            updatedCrystals(j).perimeter = stats(1).Perimeter;
        else
            updatedCrystals(j).perimeter = NaN;
        end
    end

    dropData(i).crystals = updatedCrystals;
end

fprintf('\n✅ Crystal-level perimeter computed and stored.\n');

%% Mean Perimeter per Drop

for i = 1:length(dropData)
    perims = [dropData(i).crystals.perimeter];
    perims = perims(~isnan(perims));  % Remove NaNs
    if isempty(perims)
        dropData(i).meanPerimeter = NaN;
        dropData(i).meanPerimeter_norm = NaN;
    else
        dropData(i).meanPerimeter = mean(perims);
        dropData(i).meanPerimeter_norm = dropData(i).meanPerimeter / dropArea;
    end
end

fprintf('\n✅ Mean perimeter per drop computed and stored.\n');


%% Voronoi Analysis with Validity Filtering

fprintf('\n Starting Voronoi analysis...\n');

for i = 1:length(dropData)
    drop = dropData(i);

    % Extract valid centroids
    centroids = arrayfun(@(c) c.centroid, drop.crystals, 'UniformOutput', false);
    centroids = cat(1, centroids{:});
    valid = ~any(isnan(centroids), 2);
    centroids = centroids(valid, :);

    % Skip if not enough centroids
    if size(centroids, 1) < 4
        dropData(i).voronoi_valid = false;
        warning('Drop %s skipped: not enough centroids.', drop.name);
        continue;
    end

    % Compute Voronoi diagram
    [V, C] = voronoin(centroids);
    dropMask = drop.drop.segmentation > 0;
    [H, W] = size(dropMask);
    B = bwboundaries(dropMask);
    bd = B{1}; bx = bd(:,2); by = bd(:,1);

    % Initialize storage
    numCells = numel(C);
    polys = cell(1, numCells);
    areas = nan(1, numCells);
    fracs = nan(1, numCells);
    validCell = false(1, numCells);
    clippedCount = 0;

    for j = 1:numCells
        idx = C{j};

        % Skip infinite or empty cells
        if isempty(idx) || any(idx == 1)
            continue;
        end

        xj = V(idx,1);
        yj = V(idx,2);

        % Clip polygon to droplet mask
        in = inpolygon(xj, yj, bx, by);
        xj = xj(in);
        yj = yj(in);

        if numel(xj) < 3
            clippedCount = clippedCount + 1;
            continue;
        end

        % Close polygon
        xj(end+1) = xj(1);
        yj(end+1) = yj(1);
        polys{j} = [xj, yj];

        % Compute area (normalized to drop)
        dropArea = sum(dropMask(:));
        areas(j) = polyarea(xj, yj) / dropArea;

        % Compute crystal pixel fraction in the region
        mcell = poly2mask(xj, yj, H, W);
        crystalMask = drop.combined_segmentation > 0;
        fracs(j) = sum(mcell(:) & crystalMask(:)) / sum(mcell(:));

        validCell(j) = true;
    end

    % Evaluate clipping ratio
    totalFinite = sum(~cellfun(@(c) isempty(c) || any(c==1), C));
    clipRatio = clippedCount / max(totalFinite, 1);

    if clipRatio > 0.3
        dropData(i).voronoi_valid = false;
        warning('Drop %s skipped: %.0f%% of Voronoi cells were clipped.', ...
            drop.name, 100*clipRatio);
        continue;
    end

    % Save crystal-level metrics
    dropData(i).voronoi_valid = true;
    dropData(i).voronoi_polys = polys;
    dropData(i).voronoi_areas = areas;
    dropData(i).crystal_fractions = fracs;
    dropData(i).voronoi_mean_area = mean(areas(validCell), 'omitnan');
    dropData(i).voronoi_mean_fraction = mean(fracs(validCell), 'omitnan');

    % Attach to individual crystals if count matches
    count = sum(valid);
    for j = 1:count
        if j <= numel(drop.crystals)
            dropData(i).crystals(j).voronoi_area = areas(j);
            dropData(i).crystals(j).crystal_fraction = fracs(j);
        end
    end
    
    % Plot diagnostic figure
    figure('Name',['Voronoi Diagnostic — ', drop.name], 'NumberTitle','off');
    imshow(drop.original_image, []); hold on;
    title(['Voronoi — ', drop.name]);

    % Draw Voronoi cells
    for j = 1:numCells
        if ~validCell(j), continue; end
        poly = polys{j};
        plot(poly(:,1), poly(:,2), 'g-', 'LineWidth', 1.5);
    end

    % Draw centroids
    plot(centroids(:,1), centroids(:,2), 'wo', 'MarkerSize', 4, 'LineWidth', 1);

    % Highlight centroids with no valid cell
    for j = 1:size(centroids,1)
        if j > numel(validCell) || ~validCell(j)
            plot(centroids(j,1), centroids(j,2), 'rx', 'MarkerSize', 8, 'LineWidth', 1.5);
        end
    end

    hold off;
end

fprintf('\n✅ Voronoi analysis complete. Valid cases saved.\n');

%% Conditions

% TEMPORARY MANUAL INPUT — to be automated later
dropData(1).composition = 'D';  % or {'mucin', 'NaCl'} values
dropData(1).temperature = 20;   % in °C
dropData(1).humidity    = 40;   % in %
dropData(1).imaging = 'brightfield';

dropData(2).composition = 'F';
dropData(2).temperature = 20;
dropData(2).humidity    = 40;
dropData(2).imaging = 'brightfield';

dropData(3).composition = 'D';
dropData(3).temperature = 22;
dropData(3).humidity = 40;
dropData(3).imaging = 'brightfield';


%%
% Collect drop-level metrics including areaDensityProfile
dropTable = table;
% Conversion factors
radiusFactor = 1.3004;
areaFactor = 1.691;

for i = 1:length(dropData)
    d = dropData(i);
    nRegions = length(d.areaDensityProfile);
    row = table;

    % Basic info
    row.DropID = string(d.name);
    row.Area = d.Area;
    row.Area_um2 = d.Area * areaFactor; % in um2
    row.Radius = d.Radius;
    row.Radius_um = d.Radius * radiusFactor; % in um
    row.numCrystals = d.numCrystals;
    row.meanArea = d.meanArea;
    row.meanArea_norm = d.meanArea_norm;
    row.meanPerimeter = d.meanPerimeter;
    row.meanPerimeter_norm = d.meanPerimeter / d.Radius;
    row.meanElongation = d.meanElongation;
    row.meanMisalignment = d.meanMisalignment;
    row.voronoi_mean_area = d.voronoi_mean_area;
    row.voronoi_mean_fraction = d.voronoi_mean_fraction;
    row.temperature = d.temperature;
    row.humidity = d.humidity;
    row.composition = string(d.composition);

    % Add area density values: one column per region
    for k = 1:nRegions
        colName = sprintf('areaDensity_R%d', k);
        row.(colName) = d.areaDensityProfile(k);
    end

    dropTable = [dropTable; row];
end

% Save table
writetable(dropTable, fullfile(outputFolder, 'drop_level_table.csv'));
fprintf('\nDrop-level table with area density per region saved.\n');

%% Initialize crystal table fields
crystalRows = [];

for i = 1:length(dropData)
    dropID = dropData(i).name;
    radius = dropData(i).Radius;
    crystals = dropData(i).crystals;
    bboxes = dropData(i).crystal_bboxes_list;
    areas = dropData(i).crystal_areas_list;
    dropArea = sum(dropMask(:));

    for j = 1:length(crystals)
        % Crystal name
        crystalID = sprintf('%s-%03d', dropID, j);

        % Bounding box
        bb = bboxes(j, :);  % [x, y, w, h]

        % Area normalization
        rawArea = double(areas(j));
        normArea = rawArea / dropArea;

        % Perimeter normalization
        rawPerim = crystals(j).perimeter;
        normPerim = rawPerim / dropArea;

        % Row struct
        row = struct(...
            'CrystalID', crystalID, ...
            'DropID', dropID, ...
            'bbox_x', bb(1), ...
            'bbox_y', bb(2), ...
            'bbox_w', bb(3), ...
            'bbox_h', bb(4), ...
            'area_pixels', rawArea, ...
            'area_norm', normArea, ...
            'perimeter_pixels', rawPerim, ...
            'perimeter_norm', normPerim, ...
            'elongation', crystals(j).elongation, ...
            'misalignment', crystals(j).misalignment, ...
            'r_normalized', crystals(j).r_normalized ...
        );

        % Optionally add Voronoi data if present
        if isfield(crystals(j), 'voronoi_area')
            row.voronoi_area = crystals(j).voronoi_area;
        end
        if isfield(crystals(j), 'crystal_fraction')
            row.crystal_fraction = crystals(j).crystal_fraction;
        end

        crystalRows = [crystalRows; row];
    end
end

% Convert to table
crystalTable = struct2table(crystalRows);

% Export
writetable(crystalTable, fullfile(outputFolder, 'crystalTable.csv'));
