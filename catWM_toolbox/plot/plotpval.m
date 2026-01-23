function h = plotpval(ax, xpos, ypos, pval, color, pexact,pname,minp)
% h = plotpval(ax, xpos, ypos, pval, color, pexact,pname,nperm)
% add pval next to significance bar in a figure
% adapted from mysigstar.m, written by Anne Urai 
% (in turn adapted from sigstar)

if ~exist('minp','var'); minp = 1/1000; end %assume 1/1000 as the lowest possible non-zero p-value
if ~exist('pname','var'); pname = ''; end %add name of test
if ~exist('pexact','var'); pexact = false; end %print exact pvalue instead of stars?
if ~exist('ax', 'var'); ax = gca; end
if ~exist('color', 'var'); color = 'k'; end

% use white background
txtBg = 'none';
fz = 16; fontweight = 'bold';
if ~pexact
    % get sig stars
if pval < 1e-3
    txt = '***';
elseif pval < 1e-2
    txt = '**';
elseif pval < 0.05
    txt = '*';
elseif ~isnan(pval),
    % this should be smaller
    txt = 'n.s.';
    %txt = '';
    fz = 11; fontweight = 'normal';
else
    return
end

else
    %get exact p value
    if pval == 0
        nzeros = ceil(log10(1/minp));%estimate leading zeros
        txt = sprintf(['p<%0.' num2str(nzeros) 'f'],minp);
        fz  = 11; fontweight = 'normal';
    elseif ~isnan(pval)
        txt = sprintf('p=%1.3f',pval);
        fz  = 11; fontweight = 'normal';
    end
end

if ~isempty(pname), txt = [txt ' (' pname ')']; end
% draw the stars in the bar
h = text(double(mean(xpos)), double(mean(ypos)), txt, ...
    'horizontalalignment', 'center', 'backgroundcolor', ...
    txtBg, 'margin', 1, 'fontsize', fz, 'fontweight', fontweight, 'color', color, 'Parent', ax);
end
