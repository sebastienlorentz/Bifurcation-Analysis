% Give a group of line handles a single legend entry: the first carries `name`,
% the rest are hidden from the legend. Collapses the repeated "DisplayName on the
% first, HandleVisibility off on the others" idiom in the plotters.
function legend_once(h, name)
    for i = 1:numel(h)
        if i == 1
            set(h(1), 'DisplayName', name);
        else
            set(h(i), 'HandleVisibility', 'off');
        end
    end
end
