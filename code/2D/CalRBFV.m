function [value] = CalRBFV(x,y,kdtree,Points,Bases,Hparas,maxSupportSize,SupportSizes,ApproxArgs,s,C,para)
%CALRBFV 计算给定点处的拟合值
[indices,dists] = rangesearch(kdtree, [x,y], maxSupportSize);
indices = indices{1};
dists = dists{1};
valsum = 0;
rbfsum = 0;
H = Calh(x,y,Points(indices,:),Bases(indices,:),Hparas(indices));
for i=1:length(indices)
    if dists(i)>=SupportSizes(indices(i))
        continue
    end
    if(para==1)
        if dists(i)<1e-7 && ApproxArgs(indices(i))<1e-7
            valsum = C(indices(i))+H(i);
            rbfsum = 1;
            break;
        end
        value = rbf(dists(i)/SupportSizes(indices(i)),ApproxArgs(indices(i)),s);
    elseif(para==2)   
        value = rbf_wend(dists(i)/SupportSizes(indices(i)));
    else
        value = rbf_shp(dists(i)/SupportSizes(indices(i)));
        if(value==inf)
            valsum=(C(indices(i))+H(i));
            rbfsum=1;
            break;
        end
    end
    valsum = valsum + (C(indices(i))+H(i))*value;
    rbfsum = rbfsum + value;
end
% 离群点处理
if rbfsum == 0
    indices = knnsearch(kdtree, [x,y], 'K', 16);
    H = Calh(x,y,Points(indices,:),Bases(indices,:),Hparas(indices));
    value = sum(H);
else
    value = valsum/rbfsum;
end
end

