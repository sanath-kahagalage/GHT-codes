clear all;
clc;
tic

tree = importdata('results/gomory-hu-tree_DTLZ7_3D.txt');
fronts = importdata('dataSets/DTLZ7_3D.txt');
capacity = importdata('results/capacity_DTLZ7_3D.txt');
aveDeg = importdata('results/degreeDistribution_DTLZ7_3D_new.mat');
ids = find(capacity==0);
AA = tree;

ratio = [];

for index = 1:length(AA)

        G = graph(AA(:,1)+1, AA(:,2)+1);
    G = rmedge(G,AA(index,1)+1, AA(index,2)+1);

   [S,C] = conncomp(G);
    ratio1 = min(C)/max(C);
    ratio = [ratio;ratio1];
end
ratio;
indAll = [];
for j = 1:10
    ind = find(ratio>=j*.1-0.1&ratio<=j*.1);
    cap  = [ind capacity(ind)];
    inM = find(min(cap(:,2))==cap(:,2));
    indAll = [indAll;cap(inM,1)];
end

silAll = [];
for j = 2:length(indAll)
    comb = nchoosek(indAll,j);
    sil = [];
    for i = 1:length(comb(:,1))
        G = graph(tree(:,1)+1, tree(:,2)+1);
        G = rmedge(G,tree(comb(i,:),1)+1, tree(comb(i,:),2)+1);
        [Sn,Cn] = conncomp(G);
        a = silhouette(fronts,Sn);
        sil = [sil; mean(a)];
    end

    silAll = [silAll; max(sil) find(sil==max(sil))];
end

numEdgesCorresMaxSil = 2 %find(silAll(:,1)==max(silAll(:,1))) % 5 for DTLZ7
G = graph(tree(:,1)+1, tree(:,2)+1);
comb = nchoosek(indAll,numEdgesCorresMaxSil+1)
G = rmedge(G,tree(comb(silAll(numEdgesCorresMaxSil,2),:),1)+1, tree(comb(silAll(numEdgesCorresMaxSil,2),:),2)+1);
[Sn,Cn] = conncomp(G);
a = silhouette(fronts,Sn)

f1 = [];
f2 = [];
f3 = [];
for i=1:length(Cn)
    ind = find(Sn==i);
    adClust = aveDeg(ind);
    inMdClust = find(adClust==max(adClust));
    f1MD = fronts(ind(inMdClust),1);
    f2MD = fronts(ind(inMdClust),2);
    f3MD = fronts(ind(inMdClust),3);
    f1 = [f1;f1MD];
    f2 = [f2;f2MD];
    f3 = [f3;f3MD];
%     ave = [ave mean(a(ind))];
%     colormap jet(8) %,'MarkerEdgeColor',col*i,'MarkerFaceColor',col*i
%     figure(1)
%     scatter3(fronts(ind,1), fronts(ind,2), fronts(ind,3),'o','filled')
%     hold on
end

%% find extreme solutions
allIds = 1:length(AA(:,1));
inCap = [allIds' capacity];
[Y,I] = sort(inCap(:,2),'ascend');
IC = inCap(I,:);
minCapInd = find(IC(:,2)<=prctile(IC(:,2),10));

exInd = [];
for k = minCapInd
ind = IC(k,1);
G = graph(tree(:,1)+1, tree(:,2)+1);
G = rmedge(G,tree(ind,1)+1, tree(ind,2)+1);
[Sn,Cn] = conncomp(G);

for i=1:length(Cn)
    ind = find(Sn==i);
    if length(ind)==min(Cn)
        exInd = [exInd;ind'];
    end

end
end
scatter3(fronts(:,1),fronts(:,2),fronts(:,3), 25, 'o', 'MarkerEdgeColor',[0.5 .5 .5],...
              'MarkerFaceColor',[0.7 .7 .7],...
              'LineWidth',1.5);
hold on
scatter3(f1(:,1),f2(:,1),  f3(:,1),25,'s','MarkerEdgeColor',[1 0 0],...
              'MarkerFaceColor',[1 0 0],...
              'LineWidth',1.5)
hold on
scatter3(fronts(exInd,1), fronts(exInd,2), fronts(exInd,3),25,'s','MarkerEdgeColor',[1 0 0],...
              'MarkerFaceColor',[1 0 0],...
              'LineWidth',1.5)
xlabel('$f_1$', 'Interpreter','latex')
ylabel('$f_2$', 'Interpreter','latex')
zlabel('$f_3$', 'Interpreter','latex')
view([32 35]);
h = legend('Pareto fronts{    }{    }','Representative solustions{    }{    }',...
    'Location','southoutside','Orientation','horizontal')
legend boxoff
set(h,'Interpreter','latex')
My_LGD = legend;
My_LGD.NumColumns = 2; 
set(gca,'FontSize',12,'linewidth',1,'TickLabelInterpreter', 'LaTex');
saveas(gcf,'results/ParetpFronts_GHT_rep_sols_DTLZ7.png')
% saveas(gcf,'results/ParetpFronts_kmeans_DTLZ1.png') 
% close all

%% save solutions
repSol = [f1(:,1) f2(:,1) f3(:,1); fronts(exInd,1) fronts(exInd,2) fronts(exInd,3)];
size(repSol)
save('results/repSol_DTLZ7_3D_GHT.mat', 'repSol')
toc
