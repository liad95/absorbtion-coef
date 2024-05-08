function changecolorbar(src, ~)
h = gco;
xl = xlim; 
yl = ylim;
ix = find(h.XData>=xl(1) &  h.XData<=xl(2));
iy = find(h.YData>=yl(1) &  h.YData<=yl(2));
C = h.CData(ix,iy);
caxis([min(C(:)) max(C(:))]);
end