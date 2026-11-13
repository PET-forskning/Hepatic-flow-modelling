function roi = read_roi_voiager(file)
%READ_ROI_VOIAGER  Read time-activity curves from a VOIager *_stat.txt file.
%
%   roi = READ_ROI_VOIAGER(file) reads all VOIs from a VOIager statistics
%   file (exported from GE VOIager; file name ending in _stat.txt).
%
%   Output struct roi (one cell per VOI)
%     name      VOI names
%     tac       mean activity concentration per frame (kBq/ml; converted
%               from Bq/ml)
%     midtimes  frame mid-times (s)
%     duration  frame durations (s)
%     size      VOI volume (ml), NaN if not found in the file
%
%   Mark Lubberink, Uppsala University Hospital (2010, volumes added 2024).

fprintf('Reading %s\n', file);
lines = textread(file, '%s', 'delimiter', '\n', 'headerlines', 0);

% Number of frames, from the 'Number Of Frames' header entry
i = find_line(lines, 'Number Of Frames');
nfr = str2double(regexp(lines{i}, '(\d+)\s*$', 'tokens', 'once'));

% The column header of the data table is the first line containing 'Mean'
hl = find_line(lines, 'Mean') - 1;

% Re-read the file from the column header onwards
lines = textread(file, '%s', 'delimiter', '\n', 'headerlines', hl);
header = lines{1};

% VOI names are found between the unit labels in the column header, e.g.
% "... [s] Aorta [BQML] Liver [BQML]"
s = strfind(header, '[s]');
b = strfind(header, '[BQML]');
if isempty(b)
    b = strfind(header, '[Bq/ml]');
end
nroi = length(b);
roi.name{1} = header(max(s)+3:b(1)-1);
for i = 2:nroi
    roi.name{i} = header(b(i-1)+6:b(i)-1);
end

% Data table: frame number, duration, mid-time, one column per VOI
fmt = [repmat('%n ', 1, 3+nroi) ' '];
cols = cell(1, 3+nroi);
[cols{:}] = textread(file, fmt, nfr, 'headerlines', hl+1);
dur = cols{2};
mid = cols{3};

% VOI volumes (mm3) on the line following 'Volume [mm3]'
iv = find_line(lines, 'Volume [mm3]', false);
if ~isempty(iv)
    volumes = sscanf(lines{iv+1}, '%f')';
else
    volumes = nan(1, nroi);
end

for i = 1:nroi
    roi.name{i}     = sscanf(roi.name{i}, '%s');   % remove white space
    roi.tac{i}      = cols{3+i} / 1000;            % Bq/ml -> kBq/ml
    roi.midtimes{i} = mid;
    roi.duration{i} = dur;
    roi.size{i}     = volumes(i) / 1000;           % mm3 -> ml
end

end


function i = find_line(lines, pattern, required)
% Index of the first line containing PATTERN.
if nargin < 3
    required = true;
end
i = find(~cellfun(@isempty, strfind(lines, pattern)), 1);
if isempty(i) && required
    error('read_roi_voiager:missingEntry', 'No line containing "%s" found.', pattern);
end
end
