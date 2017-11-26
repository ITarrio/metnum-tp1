function show_normales()

x=dlmread('buda_x.csv'); N=x;
y=dlmread('buda_y.csv'); N(:,:,2)=y;
z=dlmread('buda_z.csv'); N(:,:,3)=z;

[height,width,~] = size(N);

[X,Y] = meshgrid(1:width,1:height);

%%% codigo para subsamplear
p = 2; % paso / modificar para ver más o menos datos

N = N(1:p:end, 1:p:end,:);
X = X(1:p:end, 1:p:end,:);
Y = Y(1:p:end, 1:p:end,:);

%%% fin codigo para subsamplear

Z = zeros(size(N,1),size(N,2));

figure;
quiver3(X, Y, Z, N(:,:,1),N(:,:,2),N(:,:,3))

figure
quiver(N(:,:,1),N(:,:,2))

end