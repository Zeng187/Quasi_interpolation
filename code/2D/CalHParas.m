function [Bases,Hparas,Err] = CalHParas(Points,Normals,k,kdtree)
%CALHPARAS 预计算2D曲线重建所需二次曲线参数
%   Points为数据点n*2矩阵
%   Normals为各点单位法向量n*2矩阵
%   k为近邻数量
%   kdtree为数据点kd树
%   Bases为各点局部坐标系的基(n*4)
%   Hparas为各点建立的二次曲线参数(n*1)
%   Err为各点的二次曲线拟合误差(n*1)
n = size(Points,1);
Bases = zeros(n,4);
Hparas = zeros(n,1);
Err = zeros(n,1);
for i=1:n
    % 法向量顺时针旋转90度即为局部坐标系的u基，法向量本身为v基
    Bases(i,1) = Normals(i,2);
    Bases(i,2) = -Normals(i,1);
    Bases(i,3:4) = Normals(i,:);
    % 计算最佳曲线参数即min(sum((vi-a*ui^2).^2))
    indices = knnsearch(kdtree, Points(i,:), 'K', k(i));
    U = (Points(indices,:)-Points(i,:))*Bases(i,1:2)';%基向量单位正交
    V = (Points(indices,:)-Points(i,:))*Bases(i,3:4)';
    Hparas(i) = (U.^2\V);

end

end

