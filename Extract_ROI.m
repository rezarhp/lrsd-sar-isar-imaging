function roi = Extract_ROI(scene,Nx0,Ny0)

% Crop the scene image to extract the ROI stated in the center.

nrows= size(scene,1);
ncols= size(scene,2);
rss= ceil(nrows/2)-Ny0/2+1:ceil(nrows/2)+Ny0/2;
css= ceil(ncols/2)-Nx0/2+1:ceil(ncols/2)+Nx0/2;
rss(rss<1)=1; rss(rss>nrows)=nrows;
css(css<1)=1; css(css>ncols)=ncols;
roi=scene(rss, css); % Crop the scene image.
 
end

