function show_model()

filename = 'caballo.depth.csv'; % Apuntar al archivo de texto con sus Z
Z = dlmread(filename);

[height,width] = size(Z);

[X,Y] = meshgrid(1:width,1:height);

%%% codigo para subsamplear
p = 2; % paso / modificar para ver más o menos datos

X = X(1:p:end, 1:p:end,:);
Y = Y(1:p:end, 1:p:end,:);
Z = Z(1:p:end, 1:p:end,:);
%%% fin codigo para subsamplear

figure,surf(X,Y,Z);
figure,mesh(X,Y,Z);

end
