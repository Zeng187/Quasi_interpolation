function [H] = Calh(x,y,Points,Bases,Hparas)
%CALH 给定目标点和邻近数据点，计算目标点在邻近数据点处的h值(v-v(u))
%   x,y为目标点坐标
%   Points为需要处理的数据点(k*2)
%   Bases为各点局部坐标系的基(k*4)
%   Hparas为各点建立的二次曲线参数(k*1)
%   H为各点处h值(k*1)
U = sum(([x,y]-Points).*Bases(:,1:2), 2);
V = sum(([x,y]-Points).*Bases(:,3:4), 2);
H = V - Hparas.*U.^2;
end

