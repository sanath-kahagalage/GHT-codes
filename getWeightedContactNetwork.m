clear all;  clc; close all;

p = importdata('dataSets/DTLZ7_8D.txt');
pu = importdata('dataSets/DTLZ7_8D.txt');

ids = 0:length(p(:,1))-1;
p = [ids' p];
m = length(p(:,1));
[cIdx,cD] = knnsearch(p,p,'K',2);
conall = [];
xt = mean(pu(:,1))/1;
yt = mean(pu(:,2))/1;
zt = mean(pu(:,3))/1;
at = mean(pu(:,4))/1;
bt = mean(pu(:,5))/1;
ct = mean(pu(:,6))/1;
dt = mean(pu(:,7))/1;
et = mean(pu(:,8))/1;
thresold = xt^2 + yt^2+zt^2+at^2+bt^2+ct^2+dt^2+et^2;
for i = 1:(length(p(:,1))-1)
    con = [];
    for j = i+1:length(p(:,1))
        x = abs(p(i,2)- p(j,2));
        y = abs(p(i,3)- p(j,3));
        z = abs(p(i,4)- p(j,4));
        a = abs(p(i,5)- p(j,5));
        b = abs(p(i,6)- p(j,6));
        c = abs(p(i,7)- p(j,7));
        d = abs(p(i,8)- p(j,8));
        e = abs(p(i,9)- p(j,9));
        dis = x^2 + y^2 + z^2 + a^2 + b^2 + c^2 + d^2 + e^2;        
        if x < xt && y < yt && z < zt && a < at && b < bt && c < ct && d < dt && e < et  && dis < thresold
            con1 = [p(i,1) p(j,1) sqrt(dis)];
            con = [con;con1];
        end
    end
    conall = [conall;con];
end
xd = (pu(cIdx(:,1),1)- pu(cIdx(:,2) ,1)).^2;
yd = (pu(cIdx(:,1),2)- pu(cIdx(:,2) ,2)).^2;
zd = (pu(cIdx(:,1),3)- pu(cIdx(:,2) ,3)).^2;
ad = (pu(cIdx(:,1),4)- pu(cIdx(:,2) ,4)).^2;
bd = (pu(cIdx(:,1),5)- pu(cIdx(:,2) ,5)).^2;
cd = (pu(cIdx(:,1),6)- pu(cIdx(:,2) ,6)).^2;
dd = (pu(cIdx(:,1),7)- pu(cIdx(:,2) ,7)).^2;
ed = (pu(cIdx(:,1),8)- pu(cIdx(:,2) ,8)).^2;
conNearest = [cIdx(:,1)-1 cIdx(:,2)-1 sqrt(xd + yd + zd + ad + bd+ cd + dd+ ed)];
ad = sparse(m,m);
cnn = [conall;conNearest];
cnn = unique(cnn,'rows');
for k = 1:length(cnn(:,1))
    i = cnn(k,1)+1;
    j = cnn(k,2)+1;
    ad(i,j) = 1;
    ad(j,i) = 1;
end
deg = sum(ad);
hist(deg)
min(deg)
mean(deg)
save('results\degreeDistribution_DTLZ7_8D.mat','deg')

A = cnn;
fid = fopen(sprintf('contactData/DTLZ%d_8D.txt',7),'wt');
for ii = 1:size(A,1)
    fprintf(fid,'%g\t',A(ii,:));
    fprintf(fid,'\n');
end
fclose(fid);


