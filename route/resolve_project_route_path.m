function route_path = resolve_project_route_path(project_root, route_file, route_direction)
%RESOLVE_PROJECT_ROUTE_PATH Resolve a route MAT file in the project tree.

    if nargin < 1 || isempty(project_root)
        project_root = fileparts(fileparts(mfilename('fullpath')));
    end
    if nargin < 2 || isempty(route_file)
        error('route_file must not be empty.');
    end
    if nargin < 3
        route_direction = '';
    end

    if isstruct(route_file) && isfield(route_file, 'file')
        route_file = route_file.file;
    end
    route_file = char(string(route_file));

    if isfile(route_file)
        route_path = route_file;
        return;
    end

    route_root = fullfile(project_root, 'route');
    line7_root = fullfile(route_root, 'Guangxhou_line_7');
    candidate_paths = cell(0, 1);

    normalized_direction = normalize_optional_direction(route_direction);
    if ~isempty(normalized_direction)
        candidate_paths{end+1, 1} = fullfile(line7_root, normalized_direction, route_file); %#ok<AGROW>
    end

    candidate_paths{end+1, 1} = fullfile(line7_root, route_file); %#ok<AGROW>
    candidate_paths{end+1, 1} = fullfile(route_root, route_file); %#ok<AGROW>

    if isempty(normalized_direction)
        candidate_paths{end+1, 1} = fullfile(line7_root, 'up', route_file); %#ok<AGROW>
        candidate_paths{end+1, 1} = fullfile(line7_root, 'down', route_file); %#ok<AGROW>
    else
        other_direction = ternary(strcmp(normalized_direction, 'up'), 'down', 'up');
        candidate_paths{end+1, 1} = fullfile(line7_root, other_direction, route_file); %#ok<AGROW>
    end

    candidate_paths = unique(candidate_paths, 'stable');
    for idx = 1:numel(candidate_paths)
        if exist(candidate_paths{idx}, 'file') == 2
            route_path = candidate_paths{idx};
            return;
        end
    end

    route_path = which(route_file);
    if ~isempty(route_path)
        return;
    end

    error('Route file not found: %s', route_file);
end

function route_direction = normalize_optional_direction(route_direction)
    if nargin < 1 || isempty(route_direction)
        route_direction = '';
        return;
    end

    route_direction = lower(strtrim(char(string(route_direction))));
    switch route_direction
        case {'up', 'u'}
            route_direction = 'up';
        case {'down', 'd'}
            route_direction = 'down';
        otherwise
            error('route_direction must be ''up'' or ''down'', got: %s', route_direction);
    end
end

function out = ternary(cond, true_value, false_value)
    if cond
        out = true_value;
    else
        out = false_value;
    end
end