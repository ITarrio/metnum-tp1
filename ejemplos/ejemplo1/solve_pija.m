function solve_pija(M, v, zer)

Z = linsolve(M'*M, M'*v);
Z = Z';

% INSERTAR ZEROS
for i = 1:length(zer)
	puto = zer(i) + 1;
	if puto >= length(Z)
		Z = [Z 0];
	else
		Z = [Z(1:puto) 0 Z(puto+1:end)];
	end
end

sol = vec2mat(Z, 66)';

[height,width] = size(sol)
pause

[X,Y] = meshgrid(1:width,1:height);

%%% codigo para subsamplear
% p = 8; % paso / modificar para ver más o menos datos

% X = X(1:p:end, 1:p:end,:);
% Y = Y(1:p:end, 1:p:end,:);
% Z = Z(1:p:end, 1:p:end,:);
%%% fin codigo para subsamplear

figure,surf(X,Y,sol);
figure,mesh(X,Y,sol);