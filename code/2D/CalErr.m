function [Err] = CalErr(x,y,Points,Bases,Hparas,kdtree,maxSupportSize,SupportSizes,rthis)
Err=0;
[indices,dists] = rangesearch(kdtree, [x,y], maxSupportSize);
indices = indices{1};
dists = dists{1};
r=dists(:)./SupportSizes(indices(:));
H = Calh(x,y,Points(indices,:),Bases(indices,:),Hparas(indices));
n=0;
for i=1:length(indices)
    if dists(i)>=SupportSizes(indices(i))
        continue
    end
    n=n+1;
    %value= rbf(dists(i)/SupportSizes(indices(i)),1,-0.5);
    if(r(i)>1e-7)
        %Err=Err+abs(H(i))*((1-r(i))^4)/(sum((1-r).^4));
%         Err=Err+abs(H(i))*((1-r(i)))/(sum((1-r)));
        Err=Err+abs(H(i))^2*((1-r(i))^4);
    end
end
end
