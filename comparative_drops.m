%% comparative_drops.m - Compare two or more experimental groups using clean data
clear; clc; close all;

%% 1. Setup and Initialization
baseFolder = '/Users/anabellehassan/Documents/Biomedical eng/TFG/Segment/';
folderList = {'C1-2_Dtest6','C1-4_Dtest21','C1-3_FTest20','CNB_040425_02_bf_20250510_2154',...
    'CNB_040425_03_20250511_1148','CNB_040425_05_bf_20250511_1330'};  % modify as needed

%CNB_040425_02_bf_20250510_2154
%CNB_040425_03_20250511_1148
%CNB_040425_05_bf_20250511_1330
%CNB_040425_04_20250511_1307 no está bien
%CNB_040425_06_20250511_1518 no es el mejor
numFolders  = numel(folderList);
numEllipses = 10;

% Preallocate dropData
dropData(numFolders) = struct();

%% 2. Load Data and Compute Basic Geometry
for i = 1:numFolders
    name = folderList{i};
    path = fullfile(baseFolder, name);
    fprintf('Loading: %s\n', name);
    
    % Load segmentation and bounding boxes
    segData = load(fullfile(path, 'combined_segmentation.mat'));
    bbData  = load(fullfile(path, 'bboxes.mat'));
    areaData= load(fullfile(path, 'areas.mat'));
    maskData= load(fullfile(path, 'first_mask_drop.mat'));

    % Remove duplicate boxes
    [bbUniq, ia] = unique(bbData.bboxes, 'rows');
    bbAreas = areaData.areas(ia);

    % Read original image
    tifs = dir(fullfile(path, '*.tif'));
    img = imread(fullfile(path, tifs(1).name));

    % Droplet mask properties
    mask = maskData.segmentation > 0;
    dropData(i).dropMask = maskData.segmentation > 0;
    areaVal  = double(maskData.area);
    radius   = sqrt(areaVal/pi);
    [ry, rx] = find(mask);
    center   = [mean(rx), mean(ry)];

    % Compute concentric ellipses
    
    % Extracts and stores key parameters of the Drop Mask for ellipse
    % fitting
    props = regionprops(mask, 'MajorAxisLength','MinorAxisLength','Orientation','Centroid');
    if numel(props)>1
        mask = mask>0;
        props = regionprops(mask,'MajorAxisLength','MinorAxisLength','Orientation','Centroid');
    end
    MA = props(1).MajorAxisLength;
    ma = props(1).MinorAxisLength;
    ori = -deg2rad(props(1).Orientation);
    cen = props(1).Centroid;
    
    % Begins a loop to compute the parameters of each concentric ellipse
    for k = 1:numEllipses
        sf = (numEllipses-k+1)/numEllipses; % scaling factor
        a  = sf*MA/2; % a and b are the semi-axes of the scaled ellipse
        b  = sf*ma/2;
        theta = linspace(0,2*pi,100); % generates 100 points along the perimeter of the ellipse
        pts = [a*cos(theta); b*sin(theta)]; % parametric equations
        R   = [cos(ori), -sin(ori); sin(ori), cos(ori)]; %2D rotation matrix to rotate ellipse by the droplet orientation angle
        pts = R*pts + cen'; % translates it to be centered at the droplet centroid
        ellipseInfo(k) = struct('Points',pts); % stores the coordinates of the ellipse in a structure
    end
    % These ellipses get smaller as k increases, moving inward from edge to center
    
    % Compute centroids of each crystal using bbox
    numBB = size(bbUniq,1);
    bbCentroids = nan(numBB,2);
    for k = 1:numBB
        bb = round(bbUniq(k,:));  % [x, y, w, h]
        % Crop the segmentation mask to the bounding box region
        x1 = bb(1); y1 = bb(2); w = bb(3); h = bb(4);
        x2 = x1 + w - 1;
        y2 = y1 + h - 1;
        % Make sure the bounding box stays within image bounds
        x1 = max(1, x1); y1 = max(1, y1);
        x2 = min(size(segData.combined_segmentation,2), x2);
        y2 = min(size(segData.combined_segmentation,1), y2);
        % Extract the crystal mask within the bounding box
        localMask = segData.combined_segmentation(y1:y2, x1:x2) > 0;
        % Compute the centroid from the actual mask
        props = regionprops(localMask, 'Centroid');
        if ~isempty(props)
            % Centroid is relative to the cropped region—convert to full image coords
            cx = props(1).Centroid(1) + x1 - 1;
            cy = props(1).Centroid(2) + y1 - 1;
            bbCentroids(k,:) = [cx, cy];
        end
    end


    % Assign to structure
    dropData(i).name         = name;
    dropData(i).segmentation = single(segData.combined_segmentation);
    dropData(i).bboxes       = bbUniq;
    dropData(i).areas        = bbAreas;
    dropData(i).firstMask    = struct('Radius',radius,'Center',center);
    dropData(i).ellipseInfo  = ellipseInfo;
    dropData(i).bbCentroids  = bbCentroids;
    dropData(i).originalImg  = img;

    % Assign condition group
    if contains(name,'est','IgnoreCase',true)
        dropData(i).group = 'Conditions1';
    else
        dropData(i).group = 'Conditions2';
    end
end

%% 3. Split into Groups
group1 = dropData(strcmp({dropData.group}, 'Conditions1'));
group2 = dropData(strcmp({dropData.group}, 'Conditions2'));
groups = {'Conditions1','Conditions2'};

%% DEBUG BBs & Centroids (usa x,y,w,h como en tu estructura)
for i = 1:numFolders
    name   = folderList{i};
    path   = fullfile(baseFolder, name);
    
    % 1) Carga de bboxes.mat
    bbData = load(fullfile(path,'bboxes.mat'));
    bboxes = bbData.bboxes;          % [nBB x 4]: [x y w h]
    nBB    = size(bboxes,1);
    
    % 2) Leo la imagen original
    tifs = dir(fullfile(path,'*.tif'));
    img  = imread(fullfile(path,tifs(1).name));
    
    % 3) Ploteo
    figure('Name',sprintf('DEBUG BB+Centroids — %s',name),...
           'NumberTitle','off');
    imshow(img,[]); hold on;
    
    for k = 1:nBB
        % extraigo tu x,y,w,h
        x = bboxes(k,1);
        y = bboxes(k,2);
        w = bboxes(k,3);
        h = bboxes(k,4);
        
        % a) dibujo la caja
        rectangle('Position',[x y w h],...
                  'EdgeColor','r','LineWidth',1.5);
        
        % b) calculo y dibujo el centro
        cx = x + w/2;
        cy = y + h/2;
        plot(cx,cy,'bo','MarkerFaceColor','b','MarkerSize',6);
        
        % c) marco el índice
        text(cx+3, cy+3, sprintf('%d',k),...
             'Color','y','FontWeight','bold','FontSize',10);
    end
    
    title(sprintf('DEBUG: %d boxes cargadas en %s',nBB,name),...
          'Interpreter','none');
    axis image off; hold off;
end

%% 4. Display Masks & Ellipses
for gi = 1:2
    grp = eval(sprintf('group%d',gi));
    lbl = groups{gi};
    % Combined & Individual Masks
    figure('Name',['Masks - ', lbl]);
    for j=1:numel(grp)
        subplot(numel(grp),2,(j-1)*2+1);
        imagesc(grp(j).segmentation); axis image; colormap jet; colorbar;
        title(['Combined Mask - ',grp(j).name]);
        subplot(numel(grp),2,(j-1)*2+2);
        L = bwlabel(grp(j).segmentation>0);
        imshow(label2rgb(L,'jet','k','shuffle'));
        title(['Individual Masks - ',grp(j).name]);
    end
    % Concentric Ellipses
    figure('Name',['Ellipses - ', lbl]);
    for j=1:numel(grp)
        subplot(1,numel(grp),j);
        imshow(grp(j).segmentation,[]); hold on;
        for k=1:numEllipses
            pts=grp(j).ellipseInfo(k).Points;
            plot(pts(1,:),pts(2,:),'w-','LineWidth',1.5);
        end
        axis image; title(grp(j).name);
    end
end


%% --------------------------- AREA-GRAL ---------------------------------

%% 5. Compute Ring-Wise Crystal Density for Each Group

% Preallocate container
densities = cell(1,2);
Rnorm = linspace(0.1, 1, numEllipses);

for g = 1:2
    % select drops in this group
    grp = dropData(strcmp({dropData.group}, groups{g}));
    n   = numel(grp);
    D   = nan(n, numEllipses);
    
    for i = 1:n
        seg = grp(i).segmentation > 0;
        [M,N] = size(seg);
        
        % outer rings
        for k = 1:numEllipses-1
            P1 = grp(i).ellipseInfo(k).Points;
            P2 = grp(i).ellipseInfo(k+1).Points;
            m1 = poly2mask(P1(1,:),P1(2,:),M,N);
            m2 = poly2mask(P2(1,:),P2(2,:),M,N);
            ring   = m1 & ~m2;
            D(i,k) = sum(seg(ring)) / sum(ring(:));
        end
        
        % innermost
        PL = grp(i).ellipseInfo(end).Points;
        mL = poly2mask(PL(1,:),PL(2,:),M,N);
        D(i,end) = sum(seg(mL)) / sum(mL(:));
    end
    
    densities{g} = D;
end

%% 5a. Plot Per-Drop Curves + Black Mean for Each Group

for g = 1:2
    grp = dropData(strcmp({dropData.group}, groups{g}));
    D   = densities{g};
    n   = size(D,1);
    
    % pick distinct colors for each drop
    cmap = lines(n);
    
    figure('Name',['Density — ', groups{g}],'NumberTitle','off');
    
    % per-drop subplot
    subplot(1,2,1); hold on;
      for i = 1:n
        plot(Rnorm, fliplr(D(i,:)), 'Color', cmap(i,:), 'LineWidth',1.5, ...
             'DisplayName', grp(i).name);
      end
      if n>1
        legend('Interpreter','none','Location','best');
      end
      xlabel('r / R_c'); ylabel('Normalized Area');
      title(['Per Drop — ', groups{g}]);
      grid on;
    hold off;
    
    % mean±std subplot
    subplot(1,2,2); hold on;
      m = mean(fliplr(D),1,'omitnan');
      s = std(fliplr(D),0,1,'omitnan');
      
      % shade ±σ only if >1 drop
      if n>1
        x_fill = [Rnorm, fliplr(Rnorm)];
        y_fill = [m+s,   fliplr(m-s)];
        fill(x_fill, y_fill, [0.5 0.5 0.5], ...
             'FaceAlpha',0.3,'EdgeColor','none');
      end
      
      plot(Rnorm, m, 'k-', 'LineWidth',2, 'DisplayName','Mean');
      legend('Location','best');
      xlabel('r / R_c'); ylabel('Normalized Area');
      title('Mean ± Std');
      grid on;
    hold off;
end

%% 6. Comparative Mean ± Std Between the Two Groups

m1 = mean(densities{1},1,'omitnan');  s1 = std(densities{1},0,1,'omitnan');
m2 = mean(densities{2},1,'omitnan');  s2 = std(densities{2},0,1,'omitnan');

figure('Name','Comparison — Mean ± Std','NumberTitle','off'); hold on;
  x_fill = [Rnorm, fliplr(Rnorm)];
  
  % group1 (blue)
  if size(densities{1},1)>1
    fill(x_fill, [m1+s1, fliplr(m1-s1)], [0,0.4470,0.7410], ...
         'FaceAlpha',0.15,'EdgeColor','none');
  end
  p1 = plot(Rnorm, m1, '-', 'Color',[0,0.4470,0.7410], 'LineWidth',2);
  
  % group2 (orange)
  if size(densities{2},1)>1
    fill(x_fill, [m2+s2, fliplr(m2-s2)], [0.8500,0.3250,0.0980], ...
         'FaceAlpha',0.15,'EdgeColor','none');
  end
  p2 = plot(Rnorm, m2, '-', 'Color',[0.8500,0.3250,0.0980], 'LineWidth',2);
  
  legend([p1,p2], groups, 'Location','best');
  xlabel('r / R_c'); ylabel('Normalized Area');
  grid on;
hold off;

%% 7. Number of Crystals per Drop

% Split by conditions
dataLow  = dropData(strcmp({dropData.group}, 'Conditions1'));
dataHigh = dropData(strcmp({dropData.group}, 'Conditions2'));
groups   = {'Conditions1','Conditions2'};

% Count crystals = number of bounding boxes per drop
countsLow  = cellfun(@(b) size(b,1), {dataLow.bboxes})';
countsHigh = cellfun(@(b) size(b,1), {dataHigh.bboxes})';

% Boxplot of per-drop counts
figure('Name','Crystal Count Distribution','NumberTitle','off');
allCounts  = [countsLow; countsHigh];
groupLabel = [ repmat(groups(1), numel(countsLow),  1);
               repmat(groups(2), numel(countsHigh), 1) ];
boxplot(allCounts, groupLabel, 'LabelOrientation','inline');
ylabel('Number of Crystals (bboxes)');
title('Distribution of Crystal Counts per Drop');
grid on;

% Bar plot of mean ± SEM
meanLow  = mean(countsLow);
semLow   = std(countsLow)/sqrt(numel(countsLow));
meanHigh = mean(countsHigh);
semHigh  = std(countsHigh)/sqrt(numel(countsHigh));

figure('Name','Mean Crystal Count ± SEM','NumberTitle','off');
bar(1:2, [meanLow, meanHigh], 'FaceColor',[0.7 0.7 0.7]); hold on;
errorbar(1:2, [meanLow, meanHigh], [semLow, semHigh], 'k.', 'LineWidth',1.5);
hold off;
set(gca, 'XTick',1:2, 'XTickLabel',groups);
ylabel('Mean Number of Crystals');
title('Average Crystal Count by Group');
grid on;

%% 8. Mean Crystal Size per Drop

% Compute mean crystal area per drop (using pre-loaded ‘areas’)
nLow  = numel(dataLow);
nHigh = numel(dataHigh);
%meanSizeLow  = cellfun(@mean, {dataLow.areas})';
%meanSizeHigh = cellfun(@mean, {dataHigh.areas})';

% Normalize crystal size by actual droplet area in pixels
meanSizeLow = arrayfun(@(d) mean(d.areas) / sum(d.dropMask(:)), dataLow)';
meanSizeHigh = arrayfun(@(d) mean(d.areas) / sum(d.dropMask(:)), dataHigh)';

% Bar plot of mean ± SEM
meanLow  = mean(meanSizeLow);
semLow   = std(meanSizeLow)/sqrt(nLow);
meanHigh = mean(meanSizeHigh);
semHigh  = std(meanSizeHigh)/sqrt(nHigh);

figure('Name','Mean Crystal Size ± SEM','NumberTitle','off');
bar(1:2, [meanLow, meanHigh], 'FaceColor',[0.7 0.7 0.9]); hold on;
errorbar(1:2, [meanLow, meanHigh], [semLow, semHigh], 'k.', 'LineWidth',1.5);
hold off;
set(gca, 'XTick',1:2, 'XTickLabel',groups);
ylabel('Mean Crystal Area (norm)');
title('Average Crystal Size by Group');
grid on;

%% 9. Mean Crystal Perimeter per Drop

% Compute mean crystal perimeter per drop
nLow  = numel(dataLow);
nHigh = numel(dataHigh);
meanPerimLow  = zeros(nLow,1);
meanPerimHigh = zeros(nHigh,1);

for i = 1:nLow
    mask = dataLow(i).segmentation > 0;
    cc   = bwconncomp(mask);
    stats = regionprops(cc, 'Perimeter');
    perims = [stats.Perimeter];
    dropArea = sum(dataLow(i).dropMask(:));
    meanPerimLow(i) = mean(perims) / dropArea;
end

for i = 1:nHigh
    mask = dataHigh(i).segmentation > 0;
    cc   = bwconncomp(mask);
    stats = regionprops(cc, 'Perimeter');
    perims = [stats.Perimeter];
    dropArea = sum(dataHigh(i).dropMask(:));
    meanPerimHigh(i) = mean(perims) / dropArea;
end


% Bar plot of mean ± SEM
meanLow  = mean(meanPerimLow);
semLow   = std(meanPerimLow)/sqrt(nLow);
meanHigh = mean(meanPerimHigh);
semHigh  = std(meanPerimHigh)/sqrt(nHigh);

figure('Name','Mean Crystal Perimeter ± SEM','NumberTitle','off');
bar(1:2, [meanLow, meanHigh], 'FaceColor',[0.9 0.7 0.7]); hold on;
errorbar(1:2, [meanLow, meanHigh], [semLow, semHigh], 'k.', 'LineWidth',1.5);
hold off;
set(gca, 'XTick',1:2, 'XTickLabel',groups);
ylabel('Mean Perimeter (norm)');
title('Average Crystal Perimeter by Group');
grid on;

%% ---------------------- TENSOR STRUCTURE -----------------------------

%% 10. Crystal Shape Tensor Metrics

% Parameters
numR       = numEllipses;    % from Section 1
nDrops     = numel(dropData);

% Prepare grouping cell-arrays: 2 conditions × numR regions
tensorGroups = struct();
for f = {'BBox','Tensor','Lambda','Elongation','Orientation','Distance','Misalignment','CenterVectorNorm'}
    fld = f{1};
    tensorGroups.(fld) = cell(2, numR);
end

% Loop over drops
for i = 1:nDrops
    % Determine group index (1 or 2)
    gi = find(strcmp(dropData(i).group, groups), 1);
    if isempty(gi), continue; end
    
    bb    = dropData(i).bboxes;          % [nBB × 4]
    img   = double(dropData(i).originalImg);
    center= dropData(i).firstMask.Center(:);
    R0    = dropData(i).firstMask.Radius;
    nBB   = size(bb,1);
    
    % Preallocate struct array
    TR(nBB) = struct( ...
        'BBox',[], 'Tensor',[], 'Lambda',[], ...
        'Elongation',[], 'Orientation',[], ...
        'Distance',[], 'Misalignment',[], ...
        'CenterVectorNorm',[] );
    
    % Compute metrics per bounding box
    for k = 1:nBB
        % Get bounding box coordinates
        x1 = max(1, round(bb(k,1)));
        y1 = max(1, round(bb(k,2)));
        x2 = min(size(img,2), round(bb(k,1) + bb(k,3) - 1));
        y2 = min(size(img,1), round(bb(k,2) + bb(k,4) - 1));

        % Crop segmentation and image to bbox
        seg_crop = dropData(i).segmentation(y1:y2, x1:x2) > 0;  % crystal mask
        img_crop = double(img(y1:y2, x1:x2));                   % grayscale

        % Apply mask
        Im = img_crop .* seg_crop;
        S = sum(Im(:));
        if S == 0, continue; end

        % Compute intensity-weighted centroid
        [W, H] = meshgrid(x1:x2, y1:y2);
        W = double(W);
        H = double(H);
        xc = sum(Im(:) .* W(:)) / S;
        yc = sum(Im(:) .* H(:)) / S;

        % Compute moment arms
        dx = W(:) - xc;
        dy = H(:) - yc;

        % Compute inertia tensor
        Txx = sum((dx.^2) .* Im(:)) / S;
        Tyy = sum((dy.^2) .* Im(:)) / S;
        Txy = sum((dx .* dy) .* Im(:)) / S;
        T   = [Txx, Txy; Txy, Tyy];

        
        % 11.4 Eigen decomposition
        [V,D] = eig(T);
        lambdas = diag(D);
        [lambda_max,mi] = max(lambdas);
        lambda_min      = min(lambdas);
        v_max = V(:,mi);
        if v_max(2)<0, v_max = -v_max; end
        
        % 11.5 Metrics
        elong   = lambda_max / lambda_min;
        ori     = atan2(v_max(2), v_max(1));
        vec     = center - [xc; yc];
        dist    = norm(vec);
        vec_n   = vec / dist;
        misal   = abs(v_max(1)*vec(2) - v_max(2)*vec(1)) / (norm(v_max)*norm(vec));
        
        % 11.6 Store in TR
        TR(k).BBox             = bb(k,:);
        TR(k).Tensor           = T;
        TR(k).Lambda           = lambdas;
        TR(k).Elongation       = elong;
        TR(k).Orientation      = ori;
        TR(k).Distance         = dist;
        TR(k).Misalignment     = misal;
        TR(k).CenterVectorNorm = vec_n;
        
        % 11.7 Determine radial region j
        j = min( ceil(dist/(R0/numR)), numR );
        
        % 11.8 Append to tensorGroups
        for fld = fieldnames(tensorGroups)'
            f = fld{1};
            % 11.8 Append to tensorGroups, field by field
            fields = fieldnames(tensorGroups);
            for ff = 1:numel(fields)
                fld = fields{ff};
                val = TR(k).(fld);
                % Retrieve existing cell
                C = tensorGroups.(fld);
                cellval = C{gi,j};
                if isempty(cellval)
                    % First entry: make it the same type
                    C{gi,j} = val;
                else
                    % Subsequent entries: vertically concatenate
                    % (compatible shapes assumed: row vectors stack into matrix, structs into struct array, etc.)
                    C{gi,j} = [cellval; val];
                end
                tensorGroups.(fld) = C;
            end
        end
    end
    
    % Save back to dropData if desired
    dropData(i).tensor_results = TR;  
end


%% 11. Orientations & Radial Vectors per Composition Group

% 11.0 – Extraer los dos grupos correctamente
dataLow  = dropData(strcmp({dropData.group}, 'Conditions1'));
dataHigh = dropData(strcmp({dropData.group}, 'Conditions2'));
groupData = {dataLow, dataHigh};
groups    = {'Conditions1','Conditions2'};

% 11.1 – Plot orientation (rojo) vs. radial vector (azul) para cada cristal
for g = 1:2
    dataG = groupData{g};
    figure('Name', ['Orientations & Radial Vectors — ' groups{g}], ...
           'NumberTitle','off','Units','normalized','Position',[0,0.2,1,0.5]);
    
    for i = 1:numel(dataG)
        subplot(1, numel(dataG), i);
        imshow(dataG(i).originalImg, []);  
        hold on;
        
        % 11.2 – Cargo bboxes y tensor_results
        bboxes = dataG(i).bboxes;          % [nBB x 4] = [x y w h]
        TR      = dataG(i).tensor_results; % vector struct
        
        nBB = size(bboxes,1);
        
        for k = 1:nBB
            % 11.3 – Extraer x,y,w,h y calcular el centro
            x = bboxes(k,1);
            y = bboxes(k,2);
            w = bboxes(k,3);
            h = bboxes(k,4);
            cx = x + w/2;
            cy = y + h/2;
            
            % 11.4 – Escala para la longitud de flecha
            scale = double(w) / 4;
            
            % 11.5 – Sacar orientación y vector radial
            ori_k    = TR(k).Orientation;        % en radianes
            vecRad_n = TR(k).CenterVectorNorm;   % [2×1] norm
            
            % 11.6 – Flecha de orientación (rojo)
            u_ori = [cos(ori_k), sin(ori_k)] * scale;
            quiver(cx, cy, u_ori(1), u_ori(2), 0, ...
                   'Color','r','LineWidth',1.2,'MaxHeadSize',0.5);
            
            % 11.7 – Flecha radial hacia el centro (azul)
            u_rad = vecRad_n' * scale;
            quiver(cx, cy, u_rad(1), u_rad(2), 0, ...
                   'Color','b','LineWidth',1.2,'MaxHeadSize',0.5);
        end
        
        hold off;
        axis image off;
        title(sprintf('Drop %s', dataG(i).name), 'Interpreter','none');
    end
end


%% 12. Color Crystals by Misalignment Angle Bins

% Define bin edges and colormap
binEdges = linspace(0, 1, 17); % 16 equal-width bins
cmap     = parula(16);

% Loop over composition groups
for g = 1:2
    dataG = groupData{g};  % from Section 12
    figure('Name',['Misalignment Bins — ' groups{g}], ...
           'NumberTitle','off','Units','normalized','Position',[0,0.2,1,0.5]);
    
    % Loop over each drop in group
    for i = 1:numel(dataG)
        subplot(1, numel(dataG), i);
        
        % Label each connected crystal region from the segmentation
        segMask = dataG(i).segmentation > 0;  
        L       = bwlabel(segMask);
        numLab  = max(L(:));
        
        % Assign a color to each label based on its misalignment
        misVals = [dataG(i).tensor_results.Misalignment];
        colors  = zeros(numLab, 3);
        for lbl = 1:numLab
            if lbl <= numel(misVals) && ~isnan(misVals(lbl))
                binIdx = find(misVals(lbl) >= binEdges, 1, 'last');
                binIdx = min(binIdx,16);
                colors(lbl,:) = cmap(binIdx,:);
            else
                colors(lbl,:) = [0,0,0];
            end
        end
        
        % Convert label matrix to RGB image with custom palette
        rgb = label2rgb(L, colors, 'k', 'shuffle');
        imshow(rgb);
        title(sprintf('Drop %s', dataG(i).name), 'Interpreter','none');
        axis image off;
    end
    
    % Shared colorbar
    colormap(cmap);
    cb = colorbar('southoutside');
    cb.Ticks      = linspace(0,1,16);
    cb.TickLabels = arrayfun(@(b)sprintf('%.2f',b), binEdges(1:end-1),'UniformOutput',false);
    cb.Label.String = 'Misalignment';
end

%% 13. new eliminating non elongated ones (arreglar luego)
%% 13. Boxcharts of Tensor Metrics by Region & Composition (optional filtering)
% --- Config ---
excludeLowElongation = true;    % Set to false to include all crystals
elongationThreshold  = 1.5;     % Minimum elongation value to include

% Metrics and labels
metrics  = {'Elongation','Misalignment','Orientation'};
ylabels  = {'Elongation (\lambda_{max}/\lambda_{min})', ...
            'Misalignment (sin \theta)', ...
            'Orientation (rad)'};
offset   = 0.10;  % horizontal shift between groups

% Colors
col1 = [0, 0.4470, 0.7410];      % MATLAB blue
col2 = [0.8500, 0.3250, 0.0980]; % MATLAB orange

for m = 1:numel(metrics)
    metric    = metrics{m};
    ylabelTxt = ylabels{m};

    figure('Name',[metric ' by Region & Composition'], 'NumberTitle','off');
    hold on;

    % Store invisible handles for legend
    hLegend = gobjects(1,2);
    hLegend(1) = plot(nan, nan, 's', 'MarkerFaceColor', col1, 'MarkerEdgeColor', 'none');
    hLegend(2) = plot(nan, nan, 's', 'MarkerFaceColor', col2, 'MarkerEdgeColor', 'none');
    
    % Set Y-axis limits per metric
    switch metric
        case 'Elongation'
            ylim([0 30]);  % adjust depending on your expected max
        case 'Misalignment'
            ylim([0 1]);   % misalignment is sin(θ), always [0,1]
        case 'Orientation'
            ylim([0 pi]);  % orientation in radians
    end

    for g = 1:2
        thisCol = col1; if g==2, thisCol = col2; end

        for j = 1:numR
            vals = tensorGroups.(metric){g,j};

            % Apply elongation filter if enabled and metric is not Elongation
            if excludeLowElongation && ~strcmp(metric,'Elongation')
                elongVals = tensorGroups.Elongation{g,j};
                validIdx  = elongVals >= elongationThreshold;
                vals = vals(validIdx);
            end

            if isempty(vals), continue; end

            x = j + (g - 1.5) * offset; % fixed position per region
            boxchart(x * ones(size(vals)), vals, ...
                     'BoxWidth', offset * 0.8, ...
                     'BoxFaceColor', thisCol, ...
                     'BoxFaceAlpha', 0.6, ...
                     'LineWidth', 1.2,...
                    'MarkerStyle', 'none');
        end
    end

    % Set x-axis ticks and labels as R / R0
    xticks(1:numR);
    xticklabels(arrayfun(@(i) sprintf('%.2f', i / numR), 1:numR, 'UniformOutput', false));
    xlabel('r / R_0');
    ylabel(ylabelTxt);
    title([metric ' vs. Radial Region by Composition']);
    legend(hLegend, groups, 'Location', 'northwest');
    set(gca, 'TickLabelInterpreter', 'none');  % optional: avoid LaTeX
    xtickangle(30);  % tilt labels to prevent overlap
    grid on;
    hold off;
end

%% 14. Misalignment vs. Elongation — Comparative View

% 1) Gather data
E_all = []; M_all = []; G_all = [];
for g = 1:2
    for r = 1:numR
        Es = tensorGroups.Elongation{g,r};
        Ms = tensorGroups.Misalignment{g,r};
        E_all = [E_all; Es(:)];
        M_all = [M_all; Ms(:)];
        G_all = [G_all; repmat(g, numel(Es),1)];
    end
end

% --- before Section 15c, after gathering E_all ---

% choose number of bins
nQ = 10;

% compute quantile edges
q = quantile(E_all, linspace(0,1,nQ+1));

% enforce strictly increasing (in case of duplicates at extremes)
binEdges = unique(q);
numBins  = numel(binEdges)-1;

% generate human‐readable labels
binLabels = arrayfun(@(i) sprintf('%.2f–%.2f', ...
                    binEdges(i), binEdges(i+1)), ...
                    1:numBins, 'UniformOutput',false);

% Digitize elongations
[~,~,binIdx] = histcounts(E_all, binEdges);

% 4) Plot
figure('Name','Misalignment by Elongation Bin & Composition','NumberTitle','off');
hold on;

col1 = [0,0.4470,0.7410];    % blue
col2 = [0.8500,0.3250,0.0980];% orange
offset = 0.15;

% invisible points for legend
hBlue   = plot(nan, nan, 's', 'MarkerFaceColor',col1, 'MarkerEdgeColor','none');
hOrange = plot(nan, nan, 's', 'MarkerFaceColor',col2, 'MarkerEdgeColor','none');

% 5) Plot the paired boxcharts
for g = 1:2
    if g == 1
        thisCol = col1;
    else
        thisCol = col2;
    end
    
    for b = 1:numBins
        idx = (G_all==g) & (binIdx==b);
        if ~any(idx), continue; end
        
        x = b + (g-1.5)*offset;  % left/right shift
        h = boxchart(x*ones(sum(idx),1), M_all(idx), 'BoxWidth', offset*0.8);
        h.BoxFaceColor = thisCol;
        h.BoxFaceAlpha = 0.6;
        h.LineWidth    = 1.2;
    end
end

% 6) Final formatting
xticks(1:numBins);
xticklabels(binLabels);
xlabel('Elongation Bin (\lambda_{max}/\lambda_{min})');
ylabel('Misalignment (sin θ)');
title('Misalignment by Elongation Bin & Composition');
legend([hBlue, hOrange], groups, 'Location','northwest');
grid on;
hold off;


%% --------------------------- VORONOI ---------------------------------

%% 15. Voronoi Analysis for Each Drop

% (Make sure you added in your load‐loop:
%   dropData(i).dropMask = maskData.segmentation > 0;
% )

% Preallocate storage
for i = 1:numel(dropData)
    dropData(i).voronoi_polys     = {};
    dropData(i).voronoi_areas     = [];
    dropData(i).crystal_fractions = [];
end

figure('Name','Voronoi Diagrams for Each Drop','NumberTitle','off', ...
       'Units','normalized','Position',[0,0.2,1,0.6]);

for i = 1:numel(dropData)
    % 1) Extract and clean crystal centroids
    ctr = dropData(i).bbCentroids;               % your crystal centroids
    ctr = ctr(~any(isnan(ctr),2),:);
    if size(ctr,1) < 3
        warning('Drop %s has fewer than 3 centroids → skipping.', dropData(i).name);
        continue;
    end
    
    % 2) Compute Voronoi diagram
    [V, C] = voronoin(ctr);

    % 3) Get droplet boundary for clipping
    maskBW = dropData(i).dropMask;              % <-- use dropMask here
    [hMask, wMask] = size(maskBW);
    B      = bwboundaries(maskBW);
    bd     = B{1};
    bx     = bd(:,2);
    by     = bd(:,1);

    % Prepare arrays
    nCells = numel(C);
    polys  = cell(1,nCells);
    areas  = nan(1,nCells);
    fracs  = nan(1,nCells);

    % 4) Clip, store, and plot each finite cell
    subplot(1,numel(dropData),i);
    imshow(dropData(i).originalImg,[]);        % background image
    hold on
    title(['Voronoi – ', dropData(i).name]);

    crystalMask = dropData(i).segmentation > 0;

    for j = 1:nCells
        idx = C{j};
        if any(idx==1), continue; end        % skip infinite cells
        
        xj = V(idx,1); yj = V(idx,2);
        in = inpolygon(xj,yj,bx,by);
        xj = xj(in); yj = yj(in);
        if numel(xj) < 3, continue; end

        % close polygon
        xj(end+1) = xj(1);
        yj(end+1) = yj(1);

        % store
        polys{j} = [xj,yj];
        % Normalize by droplet area
        dropArea = sum(maskBW(:));
        areas(j) = polyarea(xj,yj) / dropArea;

        % areas(j) = polyarea(xj,yj);

        % crystal fraction
        mcell     = poly2mask(xj,yj,hMask,wMask);
        totalPx   = sum(mcell(:));
        crystalPx = sum(mcell(:) & crystalMask(:));
        fracs(j)  = crystalPx/totalPx;

        % draw cell
        patch(xj, yj, [0.7 0.7 0.7], 'FaceAlpha',0.4, 'EdgeColor','w');
    end

    % 5) Overlay centroids
    plot(ctr(:,1), ctr(:,2), 'wo','MarkerFaceColor','w','MarkerSize',4);
    hold off

    % Save back into dropData
    dropData(i).voronoi_polys     = polys;
    dropData(i).voronoi_areas     = areas;
    dropData(i).crystal_fractions = fracs;
end

sgtitle('Voronoi Diagrams Clipped to Droplet Mask');

%% 16. Voronoi Cells Colored by Area

figure('Name','Voronoi Cells Colored by Area','NumberTitle','off', ...
       'Units','normalized','Position',[0,0.2,1,0.6]);

cmap_area = jet(256);

for i = 1:numel(dropData)
    subplot(1,numel(dropData),i);
    
    % white background
    imshow(ones([size(dropData(i).segmentation),3]));
    hold on;
    
    % normalized area range for this drop
    areas = dropData(i).voronoi_areas;
    valid = ~isnan(areas);
    minA  = min(areas(valid));
    maxA  = max(areas(valid));
    
    % draw each Voronoi cell colored by its area
    for j = 1:numel(dropData(i).voronoi_polys)
        poly = dropData(i).voronoi_polys{j};
        if isempty(poly), continue; end
        
        x = poly(:,1);
        y = poly(:,2);
        a = areas(j);
        
        % map area to colormap index
        idx = round( (a - minA)/(maxA - minA)*255 ) + 1;
        idx = max(1, min(256, idx));
        cellColor = cmap_area(idx,:);
        
        patch(x, y, cellColor, ...
              'EdgeColor','k','LineWidth',1.2, ...
              'FaceAlpha',0.8);
    end
    
    % overlay actual crystal boundaries
    B = bwboundaries(dropData(i).segmentation > 0);
    for k = 1:numel(B)
        b = B{k};
        plot(b(:,2), b(:,1), 'k', 'LineWidth',1);
    end
    
    axis image off;
    % title(sprintf('Voronoi %s', dropData(i).name), 'Interpreter','none');
    hold off;
end

% colorbar (using last computed minA/maxA)
colormap(cmap_area);
caxis([minA maxA]);
colorbar('Position',[0.93,0.2,0.015,0.6]);
% sgtitle('Voronoi Cells Colored by Area');

%% 17. Voronoi Cells Colored by Crystal Fraction in Cell

figure('Name','Voronoi Cells Colored by Crystal Fraction','NumberTitle','off', ...
       'Units','normalized','Position',[0,0.2,1,0.6]);

cmap_frac = jet(256);    % ← use the same jet colormap

for i = 1:numel(dropData)
    subplot(1,numel(dropData),i);
    imshow(ones([size(dropData(i).segmentation),3]));
    hold on;
    
    % get valid fractions and their range
    fracs = dropData(i).crystal_fractions;
    valid = ~isnan(fracs);
    minF  = min(fracs(valid));
    maxF  = max(fracs(valid));
    
    % draw each cell colored by its crystal fraction
    for j = 1:numel(dropData(i).voronoi_polys)
        poly = dropData(i).voronoi_polys{j};
        if isempty(poly) || ~valid(j), continue; end
        
        x = poly(:,1);
        y = poly(:,2);
        f = fracs(j);
        
        % map frac → [1..256]
        idx = round((f - minF)/(maxF - minF)*255) + 1;
        idx = max(1, min(256, idx));
        cellColor = cmap_frac(idx,:);
        
        patch(x, y, cellColor, ...
              'EdgeColor','none',...
              'FaceAlpha',0.75);
    end
    
    % overlay true crystal boundaries
    B = bwboundaries(dropData(i).segmentation > 0);
    for k = 1:numel(B)
        b = B{k};
        plot(b(:,2), b(:,1), 'k', 'LineWidth',1);
    end
    
    axis image off;
    title(sprintf('Voronoi – %s (Frac.)', dropData(i).name), 'Interpreter','none');
    hold off;
end

% use same jet colormap and scaling
colormap(cmap_frac);
caxis([minF maxF]);
hcb = colorbar('Position',[0.93,0.2,0.015,0.6]);
ylabel(hcb,'Crystal Fraction');
sgtitle('Voronoi Cells Colored by Crystal Fraction');

%% 18. Crystal Fraction vs Normalized Distance to Center

figure('Name','Crystal Fraction vs r/R_c','NumberTitle','off', ...
       'Units','normalized','Position',[0,0.8,1,0.6]);

% build a global range for the color scale
allFracs = [];
for i = 1:numel(dropData)
    allFracs = [allFracs; dropData(i).crystal_fractions(:)];
end
vmin = min(allFracs(~isnan(allFracs)));
vmax = max(allFracs(~isnan(allFracs)));

DropName = {'Drop1','Drop2','Drop3'};

for i = 1:numel(dropData)
    subplot(1,numel(dropData),i);
    
    % droplet center & radius
    center = dropData(i).firstMask.Center;  
    R0     = dropData(i).firstMask.Radius;
    
    % crystal centroids and corresponding fractions
    C = dropData(i).bbCentroids;         % [nBoxes×2]
    F = dropData(i).crystal_fractions(:);
    
    % filter to valid entries
    valid = ~any(isnan(C),2) & ~isnan(F);
    C = C(valid,:);
    F = F(valid);
    
    % normalized distance
    d = sqrt((C(:,1)-center(1)).^2 + (C(:,2)-center(2)).^2);
    Rn = d / R0;
    
    % scatter, colored by fraction
    sc = scatter(Rn, F, 30, F, 'filled');
    sc.MarkerEdgeColor = 'k';
    sc.LineWidth       = 0.5;
    colormap(jet(256));
    caxis([vmin vmax]);
    colorbar;
    hold on;
    
    % lowess smoothing
    if numel(Rn) > 10
        [Rs, idx] = sort(Rn);
        Fs = F(idx);
        smoothF = smooth(Rs, Fs, 0.2, 'lowess');
        plot(Rs, smoothF, 'k-', 'LineWidth',1.5);
    end
    
    hold off;
    xlabel('r / R_c');
    ylabel('Crystal Fraction');
    %title(dropData(i).name, 'Interpreter','none');
    title(DropName(i));
    xlim([0 1]);
    ylim([vmin vmax]);
    grid on;
end

%sgtitle('Crystal Fraction vs Normalized Distance');

%% 19. Voronoi Cell Area vs Normalized Radius — By Composition

figure('Name','Voronoi Area vs R/R_c by Composition','NumberTitle','off', ...
       'Units','normalized','Position',[0,0.8,1,0.8]);

numR   = 10;
Rnorm  = linspace(0.1,1.0,numR);
groups = {'NaCl22-Mucin1','NaCl26-Mucin0.3'};
groupNames = {'Conditions1','Conditions2'};        % display names


for g = 1:2
    % filter drops for this composition
    grpDrops = dropData(strcmp({dropData.group}, groups{g}));
    nDrops   = numel(grpDrops);
    
    % preallocate
    allAreas = nan(numR, nDrops);
    
    % compute per-drop radial profiles
    for i = 1:nDrops
        center = grpDrops(i).firstMask.Center;
        R0     = grpDrops(i).firstMask.Radius;
        bins   = linspace(0, R0, numR+1);
        
        C = grpDrops(i).bbCentroids;
        A = grpDrops(i).voronoi_areas(:);
        
        for j = 1:numR
            d     = sqrt((C(:,1)-center(1)).^2 + (C(:,2)-center(2)).^2);
            inBin = d >= bins(j) & d < bins(j+1);
            allAreas(j,i) = mean(A(inBin), 'omitnan');
        end
    end
    
    % Subplot A: individual curves
    subplot(2,2,(g-1)*2+1);
    hold on;
    for i = 1:nDrops
        plot(Rnorm, allAreas(:,i), 'LineWidth',1.5, 'DisplayName', grpDrops(i).name);
    end
    hold off;
    xlabel('r / R_c');
    ylabel('Mean Voronoi Cell Area');
    %title([groups{g} ' — Each Drop']);
    title([groupNames{g} ' — Each Drop']);
    legend('Interpreter','none','Location','best');
    grid on;
    
    % Subplot B: summary (mean ± range)
    subplot(2,2,(g-1)*2+2);
    mA = mean(allAreas,2,'omitnan');
    xA = max(allAreas,[],2,'omitnan');
    nA = min(allAreas,[],2,'omitnan');
    plot(Rnorm, mA, 'k-','LineWidth',2); hold on;
    plot(Rnorm, xA, 'r--','LineWidth',1.5);
    plot(Rnorm, nA, 'b--','LineWidth',1.5);
    hold off;
    xlabel('r / R_c');
    ylabel('Voronoi Cell Area');
    %title([groups{g} ' — Mean, Max, Min']);
    title([groupNames{g} ' — Mean, Max, Min']);
    legend({'Mean','Max','Min'},'Location','best');
    grid on;
end

sgtitle('Voronoi Cell Area vs Normalized Radius by Composition');

%% 20. Comparative Mean ± 1σ: Voronoi Cell Area vs Normalized Radius

% Parameters
numR   = 10;
Rnorm  = linspace(0.1,1.0,numR);
groups = {'NaCl22-Mucin1','NaCl26-Mucin0.3'};
col1   = [0,0.4470,0.7410];
col2   = [0.8500,0.3250,0.0980];

% Preallocate storage for each group
meanAll  = cell(1,2);
stdAll   = cell(1,2);

for g = 1:2
    grpDrops = dropData(strcmp({dropData.group}, groups{g}));
    nDrops   = numel(grpDrops);
    allAreas = nan(numR, nDrops);
    
    % build radial profiles
    for i = 1:nDrops
        center = grpDrops(i).firstMask.Center;
        R0     = grpDrops(i).firstMask.Radius;
        bins   = linspace(0, R0, numR+1);
        
        C = grpDrops(i).bbCentroids;
        A = grpDrops(i).voronoi_areas(:);
        
        for j = 1:numR
            d     = sqrt((C(:,1)-center(1)).^2 + (C(:,2)-center(2)).^2);
            inBin = d >= bins(j) & d < bins(j+1);
            allAreas(j,i) = mean(A(inBin), 'omitnan');
        end
    end
    
    meanAll{g} = mean(allAreas, 2, 'omitnan');
    stdAll{g}  = std(allAreas, 0, 2, 'omitnan');
end

% Plot comparison
figure('Name','Voronoi Area Mean \pm 1\sigma by Composition','NumberTitle','off');
hold on;
% Group 1
fill([Rnorm, fliplr(Rnorm)], ...
     [meanAll{1} + stdAll{1}; flipud(meanAll{1} - stdAll{1})], ...
     col1, 'FaceAlpha',0.2, 'EdgeColor','none');
p1 = plot(Rnorm, meanAll{1}, '-', 'Color',col1, 'LineWidth',2, 'DisplayName',groupNames{1});

% Group 2
fill([Rnorm, fliplr(Rnorm)], ...
     [meanAll{2} + stdAll{2}; flipud(meanAll{2} - stdAll{2})], ...
     col2, 'FaceAlpha',0.2, 'EdgeColor','none');
p2 = plot(Rnorm, meanAll{2}, '-', 'Color',col2, 'LineWidth',2, 'DisplayName',groupNames{2});

hold off;
xlabel('r / R_c');
ylabel('Mean Voronoi Cell Area (px^2)');
title('Comparison of Voronoi Area vs Normalized Radius');
legend([p1, p2],'Location','best');
grid on;

%% 21. Crystal Fraction vs Normalized Radius — By Composition
% 
% figure('Name','Crystal Fraction vs R/R_c by Composition','NumberTitle','off', ...
%        'Units','normalized','Position',[0,0.8,1,0.8]);
% 
% numR   = 10;
% Rnorm  = linspace(0.1,1.0,numR);
% groups = {'NaCl22-Mucin1','NaCl26-Mucin0.3'};
% col1   = [0,0.4470,0.7410];
% col2   = [0.8500,0.3250,0.0980];
% groupColor = {col1, col2};
% 
% for g = 1:2
%     % select drops of this composition
%     grpDrops = dropData(strcmp({dropData.group}, groups{g}));
%     nDrops   = numel(grpDrops);
%     
%     % prepare storage
%     allFracs = nan(numR, nDrops);
%     
%     % compute per‐drop radial profiles
%     for i = 1:nDrops
%         center = grpDrops(i).firstMask.Center;
%         R0     = grpDrops(i).firstMask.Radius;
%         bins   = linspace(0, R0, numR+1);
%         
%         C = grpDrops(i).bbCentroids;           % crystal centroids
%         F = grpDrops(i).crystal_fractions(:);  % cell fractions
%         
%         for j = 1:numR
%             d     = sqrt((C(:,1)-center(1)).^2 + (C(:,2)-center(2)).^2);
%             inBin = d >= bins(j) & d < bins(j+1);
%             allFracs(j,i) = mean(F(inBin), 'omitnan');
%         end
%         % interpolate missing
%         allFracs(:,i) = fillmissing(allFracs(:,i), 'linear', 'EndValues','nearest');
%     end
%     
%     % Subplot A: individual curves
%     subplot(2,2,(g-1)*2+1);
%     hold on;
%     for i = 1:nDrops
%         plot(Rnorm, allFracs(:,i), 'Color', groupColor{g}, ...
%              'LineWidth',1.5, 'DisplayName', grpDrops(i).name);
%     end
%     hold off;
%     xlabel('R / R_c');
%     ylabel('Mean Crystal Fraction');
%     % title([groups{g} ' — Each Drop']);
%     title([groupNames{g}]);
%     legend('Interpreter','none','Location','best');
%     grid on;
%     
%     % Subplot B: summary (mean ± range)
%     subplot(2,2,(g-1)*2+2);
%     mF = mean(allFracs,2,'omitnan');
%     xF = max(allFracs,[],2,'omitnan');
%     nF = min(allFracs,[],2,'omitnan');
%     plot(Rnorm, mF, 'k-','LineWidth',2); hold on;
%     plot(Rnorm, xF, 'r--','LineWidth',1.5);
%     plot(Rnorm, nF, 'b--','LineWidth',1.5);
%     hold off;
%     xlabel('R / R_c');
%     ylabel('Crystal Fraction');
%     % title([groups{g} ' — Mean, Max, Min']);
%     title([groupNames{g} ' — Mean, Max, Min']);
%     legend({'Mean','Max','Min'},'Location','best');
%     grid on;
% end
% 
% %sgtitle('Crystal Fraction vs Normalized Radius by Composition');

%%
figure('Name','Crystal Fraction vs R/R_c by Composition','NumberTitle','off', ...
       'Units','normalized','Position',[0,0.8,1,0.8]);

numR   = 10;
Rnorm  = linspace(0.1,1.0,numR);
groups = {'NaCl22-Mucin1','NaCl26-Mucin0.3'};
groupNames = {'Condition 1', 'Condition 2'};  % for plot titles only

for g = 1:2
    % Select droplets of this group
    grpDrops = dropData(strcmp({dropData.group}, groups{g}));
    nDrops   = numel(grpDrops);

    % Preallocate
    allFracs = nan(numR, nDrops);

    % Compute radial crystal fraction profiles
    for i = 1:nDrops
        center = grpDrops(i).firstMask.Center;
        R0     = grpDrops(i).firstMask.Radius;
        bins   = linspace(0, R0, numR+1);

        C = grpDrops(i).bbCentroids;
        F = grpDrops(i).crystal_fractions(:);

        for j = 1:numR
            d     = sqrt((C(:,1)-center(1)).^2 + (C(:,2)-center(2)).^2);
            inBin = d >= bins(j) & d < bins(j+1);
            allFracs(j,i) = mean(F(inBin), 'omitnan');
        end

        % Interpolate missing values
        allFracs(:,i) = fillmissing(allFracs(:,i), 'linear', 'EndValues','nearest');
    end

    % Subplot A: individual droplet curves with unique colors
    subplot(2,2,(g-1)*2+1);
    hold on;
    cmap = lines(nDrops);  % assign different color to each droplet
    for i = 1:nDrops
        plot(Rnorm, allFracs(:,i), ...
             'Color', cmap(i,:), ...
             'LineWidth', 1.5, ...
             'DisplayName', ['Drop ' num2str(i)]);
    end
    hold off;
    xlabel('R / R_c');
    ylabel('Mean Crystal Fraction');
    title([groupNames{g}]);
    legend('Interpreter','none','Location','best');
    grid on;

    % Subplot B: group-level summary (mean ± range)
    subplot(2,2,(g-1)*2+2);
    mF = mean(allFracs, 2, 'omitnan');
    xF = max(allFracs, [], 2, 'omitnan');
    nF = min(allFracs, [], 2, 'omitnan');
    plot(Rnorm, mF, 'k-',  'LineWidth', 2, 'DisplayName', 'Mean'); hold on;
    plot(Rnorm, xF, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Max');
    plot(Rnorm, nF, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Min');
    hold off;
    xlabel('R / R_c');
    ylabel('Crystal Fraction');
    title([groupNames{g} ' — Mean, Max, Min']);
    legend('Location','best');
    grid on;
end

%sgtitle('Crystal Fraction vs Normalized Radius by Composition');


%% 22. Comparative Mean ± 1σ: Crystal Fraction vs Normalized Radius

% Parameters
numR   = 10;
Rnorm  = linspace(0.1,1.0,numR);
groups = {'NaCl22-Mucin1','NaCl26-Mucin0.3'};
col1   = [0,0.4470,0.7410];
col2   = [0.8500,0.3250,0.0980];

% Preallocate storage for each group
meanFracAll = cell(1,2);
stdFracAll  = cell(1,2);

for g = 1:2
    grpDrops = dropData(strcmp({dropData.group}, groups{g}));
    nDrops   = numel(grpDrops);
    allFracs = nan(numR, nDrops);
    
    % build radial profiles
    for i = 1:nDrops
        center = grpDrops(i).firstMask.Center;
        R0     = grpDrops(i).firstMask.Radius;
        bins   = linspace(0, R0, numR+1);
        
        C = grpDrops(i).bbCentroids;
        F = grpDrops(i).crystal_fractions(:);
        
        for j = 1:numR
            d     = sqrt((C(:,1)-center(1)).^2 + (C(:,2)-center(2)).^2);
            inBin = d >= bins(j) & d < bins(j+1);
            allFracs(j,i) = mean(F(inBin), 'omitnan');
        end
        allFracs(:,i) = fillmissing(allFracs(:,i), 'linear', 'EndValues','nearest');
    end
    
    meanFracAll{g} = mean(allFracs, 2, 'omitnan');
    stdFracAll{g}  = std(allFracs, 0, 2, 'omitnan');
end

% Plot comparison
figure('Name','Crystal Fraction Mean \pm 1σ by Composition','NumberTitle','off');
hold on;

% Group 1 band + line
fill([Rnorm, fliplr(Rnorm)], ...
     [meanFracAll{1}+stdFracAll{1}; flipud(meanFracAll{1}-stdFracAll{1})], ...
     col1, 'FaceAlpha',0.2, 'EdgeColor','none');
p1 = plot(Rnorm, meanFracAll{1}, '-', 'Color',col1, 'LineWidth',2, 'DisplayName',groupNames{1});

% Group 2 band + line
fill([Rnorm, fliplr(Rnorm)], ...
     [meanFracAll{2}+stdFracAll{2}; flipud(meanFracAll{2}-stdFracAll{2})], ...
     col2, 'FaceAlpha',0.2, 'EdgeColor','none');
p2 = plot(Rnorm, meanFracAll{2}, '-', 'Color',col2, 'LineWidth',2, 'DisplayName',groupNames{2});

hold off;
xlabel('R / R_c');
ylabel('Mean Crystal Fraction');
title('Comparison of Crystal Fraction vs Normalized Radius');
legend([p1, p2],'Location','best');
grid on;

%% 23. Save All Figures to Custom Folder

% Specify folder name
outputName = 'results concentration yes';  

% Base output folder (shared location)
baseOutputFolder = '/Users/anabellehassan/Documents/Biomedical eng/TFG/Segment/';
fullOutputPath   = fullfile(baseOutputFolder, outputName);

% Create folder if it doesn’t exist
if ~exist(fullOutputPath, 'dir')
    mkdir(fullOutputPath);
end

% Get all figure handles
figHandles = findall(0, 'Type', 'figure');

% Save each figure as high-res PNG
for f = 1:numel(figHandles)
    fig     = figHandles(f);
    figName = get(fig, 'Name');
    if isempty(figName)
        figName = sprintf('figure_%02d', f);
    end
    safeName = matlab.lang.makeValidName(figName);  % sanitize
    savePath = fullfile(fullOutputPath, [safeName, '.png']);
    
    % export at 300 DPI
    exportgraphics(fig, savePath, 'Resolution', 300);
end

fprintf('Saved all figures in: %s\n', fullOutputPath);
