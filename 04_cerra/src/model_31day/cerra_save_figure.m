function cerra_save_figure(h, figureDir, baseName)
%CERRA_SAVE_FIGURE Save PNG and editable MATLAB figure.

    if nargin < 2
        error('cerra_save_figure requires a figure handle and output path.');
    end
    if nargin < 3 || isempty(baseName)
        [figureDir, baseName, extension] = fileparts(char(figureDir));
        if isempty(figureDir)
            figureDir = '.';
        end
        if isempty(baseName)
            error('The figure output path must include a file name.');
        end
        if ~isempty(extension) && ~strcmpi(extension, '.png')
            error('Two-argument calls must provide a PNG path.');
        end
    end

    ensure_directory(figureDir);
    pngFile = fullfile(figureDir, [char(baseName), '.png']);
    figFile = fullfile(figureDir, [char(baseName), '.fig']);

    % The original two-argument interface bug called ensure_directory with
    % the complete PNG path.  A failed legacy run may therefore have left
    % an empty directory named "*.png" where this file must be written.
    % rmdir without the recursive flag removes only an empty directory, so
    % unexpected contents are never deleted silently.
    if isfolder(pngFile)
        [removed, message] = rmdir(pngFile);
        if ~removed
            error('Cannot replace legacy PNG directory "%s": %s', ...
                pngFile, message);
        end
    end
    if isfolder(figFile)
        [removed, message] = rmdir(figFile);
        if ~removed
            error('Cannot replace legacy FIG directory "%s": %s', ...
                figFile, message);
        end
    end

    try
        exportgraphics(h, pngFile, 'Resolution', 180);
    catch
        print(h, pngFile, '-dpng', '-r180');
    end
    savefig(h, figFile);
    close(h);
end


