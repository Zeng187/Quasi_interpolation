function [H] = Calh(x,y,Points,Bases,Hparas)
%CALH ����Ŀ�����ڽ����ݵ㣬����Ŀ������ڽ����ݵ㴦��hֵ(v-v(u))
%   x,yΪĿ�������
%   PointsΪ��Ҫ��������ݵ�(k*2)
%   BasesΪ����ֲ�����ϵ�Ļ�(k*4)
%   HparasΪ���㽨���Ķ������߲���(k*1)
%   HΪ���㴦hֵ(k*1)
U = sum(([x,y]-Points).*Bases(:,1:2), 2);
V = sum(([x,y]-Points).*Bases(:,3:4), 2);
H = V - Hparas.*U.^2;
end

