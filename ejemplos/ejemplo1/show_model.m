function show_model(fname)

%filename = strcat(fname, '_z.csv'); % Apuntar al archivo de texto con sus Z
%Z = dlmread(filename);
Z = dlmread(strcat(fname, '.depth.csv'));

[height,width] = size(Z);

[X,Y] = meshgrid(1:width,1:height);

%%% codigo para subsamplear
% p = 8; % paso / modificar para ver más o menos datos

% X = X(1:p:end, 1:p:end,:);
% Y = Y(1:p:end, 1:p:end,:);
% Z = Z(1:p:end, 1:p:end,:);
%%% fin codigo para subsamplear

figure,surf(X,Y,Z);
figure,mesh(X,Y,Z);

end