function show_depths(fname)

%% 
%fname = 'buda.depths.dat';

Z2 = read_depths(fname);

imagesc(Z2)
colormap(gray)

fclose(fp);
end
