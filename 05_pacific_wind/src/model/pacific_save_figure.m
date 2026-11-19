function pacific_save_figure(h, figureDir, baseName)
%PACIFIC_SAVE_FIGURE Save PNG and editable MATLAB figure.

    ensure_directory(figureDir);
    pngFile = fullfile(figureDir, [char(baseName), '.png']);
    figFile = fullfile(figureDir, [char(baseName), '.fig']);
    try
        exportgraphics(h, pngFile, 'Resolution', 180);
    catch
        print(h, pngFile, '-dpng', '-r180');
    end
    savefig(h, figFile);
    close(h);
end
