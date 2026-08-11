% Write the current figure to docs/figures/<name>.pdf as vector graphics. Strips the
% hover toolbar first, or it is baked into the PDF. Used by the example drivers that
% produce the report's figures.
function export_figure(name)
    figdir = fullfile(fileparts(mfilename('fullpath')), 'docs', 'figures');
    if ~isfolder(figdir); mkdir(figdir); end
    set(findall(gcf, 'Type', 'axes'), 'Toolbar', []);
    exportgraphics(gcf, fullfile(figdir, [name '.pdf']), 'ContentType', 'vector');
end
