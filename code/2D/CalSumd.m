function [d,count] = CalSumd(Points,x,xx,y,yy)
%返回所有不超过8个点的叶子节点的对角线长度，以及叶子节点个数
d=0;count=0;
if(length(Points)>8 && length(Points)~=0)
    midx=(x+xx)/2;
    midy=(y+yy)/2;
    Indicesx1=find(Points(:,1)<=midx);Indicesx2=setdiff([1:length(Points)]',Indicesx1);
    Indicesy1=find(Points(:,2)<=midy);Indicesy2=setdiff([1:length(Points)]',Indicesy1);
    [d1,count1]=CalSumd(Points(intersect(Indicesx1,Indicesy1),:),x,midx,y,midy);
    [d2,count2]=CalSumd(Points(intersect(Indicesx1,Indicesy2),:),x,midx,midy,yy);
    [d3,count3]=CalSumd(Points(intersect(Indicesx2,Indicesy1),:),midx,xx,y,midy);
    [d4,count4]=CalSumd(Points(intersect(Indicesx2,Indicesy2),:),midx,xx,midy,yy);
    count=count1+count2+count3+count4;
    d=d1+d2+d3+d4;
elseif(length(Points)~=0)
    count=1;d=sqrt((xx-x)^2+(yy-y)^2);
end
end