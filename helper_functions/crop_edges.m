function im = crop_edges(im,fwid)
%
%

im   = im(fwid+1:end-fwid,fwid+1:end-fwid,:); 