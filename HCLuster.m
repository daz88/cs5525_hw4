function [  ] = HCLuster( op )
%HCLUSTER Summary of this function goes here
%   Detailed explanation goes here
data=[1 0 1 1 0;1 1 0 1 0;0 0 1 1 0;0 1 0 1 0;1 0 1 0 1;0 1 1 0 0]
%1:single,2:complete,3:group average
if op==1
    disp('single');
    p1=pdist(data,@RC);
    l1=linkage(p1,'single');
    dendrogram(l1);
elseif op==2
    disp('complete');
    p2=pdist(data,@SMC);
    l2=linkage(p2,'complete');
    dendrogram(l2);
elseif op==3
    disp('group');
    p3=pdist(data,@JC);
    l3=linkage(p3,'weighted');
    dendrogram(l3);
end

end

