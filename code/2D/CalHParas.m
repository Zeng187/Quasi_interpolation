function [Bases,Hparas,Err] = CalHParas(Points,Normals,k,kdtree)
%CALHPARAS Ԥ����2D�����ؽ�����������߲���
%   PointsΪ���ݵ�n*2����
%   NormalsΪ���㵥λ������n*2����
%   kΪ��������
%   kdtreeΪ���ݵ�kd��
%   BasesΪ����ֲ�����ϵ�Ļ�(n*4)
%   HparasΪ���㽨���Ķ������߲���(n*1)
%   ErrΪ����Ķ�������������(n*1)
n = size(Points,1);
Bases = zeros(n,4);
Hparas = zeros(n,1);
Err = zeros(n,1);
for i=1:n
    % ������˳ʱ����ת90�ȼ�Ϊ�ֲ�����ϵ��u��������������Ϊv��
    Bases(i,1) = Normals(i,2);
    Bases(i,2) = -Normals(i,1);
    Bases(i,3:4) = Normals(i,:);
    % ����������߲�����min(sum((vi-a*ui^2).^2))
    indices = knnsearch(kdtree, Points(i,:), 'K', k(i));
    U = (Points(indices,:)-Points(i,:))*Bases(i,1:2)';%��������λ����
    V = (Points(indices,:)-Points(i,:))*Bases(i,3:4)';
    Hparas(i) = (U.^2\V);

end

end

