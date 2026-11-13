function [roi, data] = read_voi_HERMES(file, data)
%READ_VOI_HERMES  Read VOI time-activity curves exported from HERMES (.xlsx).
%
%   [roi, data] = READ_VOI_HERMES(file, data) reads all VOIs from a HERMES
%   Excel export. Each VOI occupies a block of rows; the block boundaries
%   are identified from the volume column (column 16), which is empty
%   outside the VOI blocks. Within a block, column 2 holds the frame start
%   time (s) and column 6 the mean activity concentration (Bq/ml).
%
%   Output struct roi (one cell per VOI)
%     name      VOI names (column 1 of the first row of each block)
%     tac       mean activity concentration per frame (kBq/ml)
%     size      VOI volume, as given in the volume column
%     midtimes  {1}: frame mid-times (s)
%     duration  {1}: frame durations (s)
%
%   Fields added to data
%     nr_frames  number of frames
%     t_start    frame start times (s)
%     t_end      frame end times (s)
%     t_tac      frame mid-times (s)
%     duration   frame durations (s); the last frame is assumed to be as
%                long as the second-to-last frame

[~, ~, raw] = xlsread(file);

% Locate the first and last row of each VOI block from the volume column
vol_col = cell2mat(raw(3:end, 16));
missing = [isnan(vol_col); true];   % treat rows beyond the end as missing
position_vois = [];
for i = 2:length(vol_col)
    if ~missing(i) && missing(i-1)          % first row of a block
        position_vois(end+1) = i + 2;       %#ok<AGROW>
    elseif ~missing(i) && missing(i+1)      % last row of a block
        position_vois(end+1) = i + 2;       %#ok<AGROW>
    end
end
block_start = position_vois(1:2:end);
block_end   = position_vois(2:2:end);
nvoi        = length(block_start);
numframes   = block_end(1) - block_start(1);

voi_names = raw(block_start);
fprintf('VOIs found in %s:\n', file);
fprintf('  %s\n', voi_names{:});

% Column 1: frame start times; columns 2..nvoi+1: mean activity per VOI
tac_data = cell(numframes, nvoi+1);
tac_data(:, 1) = raw(block_start(1)+1:block_end(1), 2);
for k = 1:nvoi
    tac_data(:, k+1) = raw(block_start(k)+1:block_end(k), 6);
end

% Numbers stored as text (possibly with spaces or decimal commas)
for idx = 1:numel(tac_data)
    if ischar(tac_data{idx})
        value = tac_data{idx};
        value = value(~isspace(value));
        value = strrep(value, ',', '.');
        tac_data{idx} = str2double(value);
    end
end
tac_data = cell2mat(tac_data);

for i = 1:nvoi
    roi.name{i} = char(voi_names(i));
    roi.tac{i}  = tac_data(:, i+1) / 1000;          % Bq/ml -> kBq/ml
    roi.size{i} = vol_col(position_vois(i*2-1));
end

% Frame timing (s), derived from the frame start times
data.t_start   = tac_data(:, 1);
data.nr_frames = numframes;
data.duration  = [diff(data.t_start); data.t_start(end) - data.t_start(end-1)];
data.t_end     = [data.t_start(2:end); data.t_start(end) + data.duration(end)];
data.t_tac     = (data.t_start + data.t_end) / 2;
roi.midtimes{1} = data.t_tac;
roi.duration{1} = data.duration;

end
