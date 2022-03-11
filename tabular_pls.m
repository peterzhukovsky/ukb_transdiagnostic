addpath('D:\Canada_2020\UK_biobank\reports');
clear; dataimport
cognitive_ordered.x6350_2_0(cognitive_ordered.x6350_2_0> (nanmean(cognitive_ordered.x6350_2_0)+4*nanstd(cognitive_ordered.x6350_2_0)) | cognitive_ordered.x6350_2_0<100)=NaN; %Duration to complete alphanumeric path 
cognitive_ordered.x6348_2_0(cognitive_ordered.x6348_2_0<100 | cognitive_ordered.x6348_2_0> (nanmean(cognitive_ordered.x6350_2_0)+3*nanstd(cognitive_ordered.x6348_2_0)) )=NaN; %Duration to complete numeric/easy path 
cognitive_ordered.tmt_cor=(cognitive_ordered.x6350_2_0 + 5*cognitive_ordered.x6351_2_0);% - (cognitive_ordered.x6348_2_0 +5*cognitive_ordered.x6349_2_0) ;
%%% requires tabular_ukb to have ran
%ix=(cellfun('isempty', clinical)' | sum(isnan(ms_ordered)')'>0 ) ; a=age; a(ix)=[];s=sex;s(ix)=[];clin=clinical; clin(ix)=[];
ix=(cellfun('isempty', clinical)'| isnan(ica_partial(:,1))) ; a=age; a(ix)=[];s=sex;s(ix)=[];clin=clinical; clin(ix)=[];
cognitive=cognitive_ordered; cognitive(ix,:)=[]; minimal=minimal_ordered; minimal(ix,:)=[]; MDD_prs=MDD_prs_ordered; MDD_prs(ix)=[];

%X=ms_ordered; X(ix,:)=[]; 
X=ica_partial; X(ix,:)=[]; 
%%% what are the y variables: alphanumtrails, fluid intelligence,  word pairs (PAL), symbol digit, %%%% pairs match tower  cognitive.x21004_2_0,
Y=[cognitive.tmt_cor, cognitive.x20016_2_0, cognitive.x20197_2_0, cognitive.x23324_2_0];%   Y=[cognitive.tmt_cor, cognitive.x21004_2_0, cognitive.x20197_2_0]; % cognitive.x399_2_2+cognitive.x399_2_1, 
%%% no controls? remove them here
 index=~(  ismember(clin','ahc')   ); % %
 index=~ismember(clin','ahc') & ~ismember(clin', 'depanx') & ~ismember(clin', 'ahc') & ~ismember(clin', 'str');
Y=Y(index,:); X=X(index,:); minimal=minimal(index,:); c=clin(index); MDD_prs=MDD_prs(index);
%%% remove nans
naninx=sum(isnan(Y)')'; Y=Y(naninx==0,:);  X=X(naninx==0,:); minimal=minimal(naninx==0,:);c=c(naninx==0);MDD_prs=MDD_prs(naninx==0);
Y=zscore(Y); X=zscore(X); Y(:,1)=Y(:,1)*-1;
ncomp=3

%%% regress out sex and age and potentially other covariates
mdl = fitlm(table(minimal.x21003_2_0, minimal.x31_0_0, minimal.x54_2_0,Y(:,1)) ); Y(:,1)=mdl.Residuals.Raw;
mdl = fitlm(table(minimal.x21003_2_0, minimal.x31_0_0, minimal.x54_2_0,Y(:,2)) ); Y(:,2)=mdl.Residuals.Raw;
mdl = fitlm(table(minimal.x21003_2_0, minimal.x31_0_0, minimal.x54_2_0,Y(:,3)) ); Y(:,3)=mdl.Residuals.Raw;
mdl = fitlm(table(minimal.x21003_2_0, minimal.x31_0_0, minimal.x54_2_0,Y(:,4)) ); Y(:,4)=mdl.Residuals.Raw;
clear x; for i=1:length(X(1,:)); mdl = fitlm(table(minimal.x21003_2_0, minimal.x31_0_0, minimal.x25741_2_0, minimal.x54_2_0,X(:,i)) ); x(:,i)=mdl.Residuals.Raw;end

[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(x,Y,ncomp);PCTVAR
%%% explore the PLS components and their correlation with Y
%XS VS YS FIGURE:%               figure(31); scatter(XS(ismember(c','dep'),1), YS(ismember(c','dep'),1), 5,'b', 'filled'); hold on;scatter(XS(ismember(c','depanx'),1), YS(ismember(c','depanx'),1), 5, 'r', 'filled'); hold on;scatter(XS(ismember(c','anx'),1), YS(ismember(c','anx'),1), 5, 'k', 'filled'); mdl= fitlm(XS(:, 1), YS(:,1)); [ypred,yci] = predict(mdl,XS(:, 1), 'Alpha',0.001); hold on; plot(XS(:, 1), ypred, 'k', 'LineWidth', 2);plot(XS(:, 1), yci(:,1), 'k', 'LineWidth', 0.5);plot(XS(:, 1), yci(:,2), 'k', 'LineWidth', 0.5);
% and their distributions: %        figure; plot_histogram_shaded(XS,'Alpha',0.3,'color',[0.5 0.5 0.5], 'Normalization', 'pdf');figure; plot_histogram_shaded(YS,'Alpha',0.3,'color',[0.5 0.5 0.5], 'Normalization', 'pdf');
figure(1);imagesc(corr(XS(:,3),Y)); colormap bone ; corr(XS, Y) %1.TMT  2.Tower  3.PAL  4.DSST
combs=allcomb([1:3], [1:4]);figure(13); 
for i=1:length(combs)
hold off; ix1=combs(i,1);ix2=combs(i,2); subplot(3,4,i); 
scatter(XS(ismember(c','dep'),ix1), Y(ismember(c','dep'),ix2), 5,'b', 'filled'); 
hold on;scatter(XS(ismember(c','depanx'),ix1), Y(ismember(c','depanx'),ix2), 5, 'r', 'filled')
hold on;scatter(XS(ismember(c','anx'),ix1), Y(ismember(c','anx'),ix2), 5, 'k', 'filled');
mdl= fitlm(XS(:, ix1), Y(:,ix2)); [ypred,yci] = predict(mdl,XS(:, ix1), 'Alpha',0.001); hold on
plot(XS(:, ix1), ypred, 'k', 'LineWidth', 2);
plot(XS(:, ix1), yci(:,1), 'k', 'LineWidth', 0.5);
plot(XS(:, ix1), yci(:,2), 'k', 'LineWidth', 0.5);
end; clear mdl ypred yci mdl ix1 ix2 combs

corr(minimal.x21003_2_0, XS)
anova1(XS(:,1), c'); max(corr(X, Y))

%% permutation testing
permutations=5000;   
allobservations=Y; 
for ncomp=1:3
    
    parfor n = 1:permutations
    % selecting either next combination, or random permutation
    permutation_index = randperm(length(allobservations));
    % creating random sample based on permutation index
    randomSample = allobservations(permutation_index,:);
    % running the PLS for this permutation
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(x,randomSample,ncomp);
    Rsq(n) = sum(PCTVAR(2,:));
    Rsq1(n) = sum(PCTVAR(1,:));
    end
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(x,Y,ncomp);
    p(ncomp)=sum(sum(PCTVAR(2,:))<Rsq')/permutations
    p_1(ncomp)=sum(sum(PCTVAR(1,:))<Rsq1')/permutations
end
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(x,Y,ncomp);PCTVAR

% figure to show the permutation distribution vs actual distribution
figure(2); histogram(Rsq); xlim([0.02 0.06]); sum(PCTVAR')
%% bootstrapping to get the func connectivity weights for PLS1, 2 and 3
dim=3
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(x,Y,dim);PCTVAR
PLS1w=stats.W(:,1);
PLS2w=stats.W(:,2);
PLS3w=stats.W(:,2);

bootnum=5000;
PLS1weights=[];
PLS2weights=[];
PLS3weights=[];

parfor i=1:bootnum
    i;
    myresample = randsample(size(x,1),size(x,1),1);
    res(i,:)=myresample; %store resampling out of interest
    Xr=x(myresample,:); % define X for resampled subjects
    Yr=Y(myresample,:); % define X for resampled subjects
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(Xr,Yr,dim); %perform PLS for resampled data
      
    newW=stats.W(:,1);%extract PLS1 weights
    if corr(PLS1w,newW)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
        newW=-1*newW;
    end
    PLS1weights=[PLS1weights,newW];%store (ordered) weights from this bootstrap run
    
   
    newW=stats.W(:,2);%extract PLS2 weights
    if corr(PLS2w,newW)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
        newW=-1*newW;
    end
    PLS2weights=[PLS2weights,newW]; %store (ordered) weights from this bootstrap run    
    
    newW=stats.W(:,3);%extract PLS2 weights
    if corr(PLS3w,newW)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
        newW=-1*newW;
    end
    PLS3weights=[PLS3weights,newW]; %store (ordered) weights from this bootstrap run    
end

PLS1sw=std(PLS1weights');
PLS2sw=std(PLS2weights');
PLS3sw=std(PLS3weights');

plsweights1=PLS1w./PLS1sw';
plsweights2=PLS2w./PLS2sw'; 
plsweights3=PLS3w./PLS3sw';

sum(plsweights1>3)
sum(plsweights2<-3)


%% Visuals
D=21; ICs={'IC1';'IC7';'IC9';'IC13';'IC14';'IC20';'IC21';'IC5';'IC6';'IC16';'IC10';'IC11';'IC12';'IC17';'IC3';'IC2';'IC4';'IC8';'IC19';'IC15';'IC18'};
combs=allcomb(ICs, ICs);
for i=1:441; comb(i)=cellstr(strcat(combs{i,2}, combs{i,1})); end
ICs=reshape(comb, [21,21]); clear comb combs i
ICs_vector=ICs(triu(ones(D),1)==1)';

ic_mapping=readtable('D:\Canada_2020\OASIS\reports\fmri\ic_order_mapping.csv');
tmp=ic_mapping.original_order_num;

myLabel=ica2yeo7.Yeo7N;myLabel=myLabel([1 7 9 13 14 20 21 5 6 16 10 11 12 17 3 2 4 8 19 15 18]);
myLabel=repmat({''}, [21,1]); figure
upperlim=4; lowerlim=-4; %upperlim=nodes(205,component); lowerlim=nodes(5,component);
Weights=plsweights1(tmp); Weights(Weights<upperlim & Weights>lowerlim)=0;%Weights=stats.W(:,component); 
Weights(Weights<0)=0; Weights_square= zeros(21); Weights_square(triu(ones(21),1)>0) = abs(Weights);
hold on; myColorMap=zeros(21,3);myColorMap(:,3)=1; circularGraph(Weights_square, 'Colormap',myColorMap, 'Label',myLabel);
Weights=plsweights1(tmp); Weights(Weights<upperlim & Weights>lowerlim)=0;
Weights(Weights>0)=0; Weights_square= zeros(21); Weights_square(triu(ones(21),1)>0) = abs(Weights);
hold on; myColorMap=zeros(21,3);myColorMap(:,1)=1; circularGraph(Weights_square, 'Colormap',myColorMap, 'Label',myLabel);

%% Comparing PLSs
load('pls_dep.mat');load('ypca_dep.mat');
plsweights(1,:)=plsweights1;
ypca=ypca_dep;%[~,ypca,~] = pca(Y);ypca=ypca(:,1); corr(ypca,Y)
CorrDep=corr(x,ypca); CorrDep=r2z(CorrDep);
load('pls_anx.mat');load('ypca_anx.mat');
plsweights(2,:)=plsweights1;
ypca=ypca_anx;%[~,ypca,~] = pca(Y);ypca=ypca(:,1); corr(ypca,Y)
CorrAnx=corr(x,ypca); CorrAnx=r2z(CorrAnx);
DepVsAnx=(CorrDep-CorrAnx)/sqrt( 1/(426-3) + 1/(1895-3) ); %2.58

corr(plsweights')
dep_vs_anx=plsweights(1,:)-plsweights(2,:);
figure;scatter(dep_vs_anx', DepVsAnx);corr(dep_vs_anx', DepVsAnx) 

%dep_vs_anx(abs(plsweights(1,:))<3 & abs(plsweights(2,:))<3)=0;
d=dep_vs_anx;
d(abs(DepVsAnx)<1.96)=0;
d(abs(dep_vs_anx)<4)=0;

myLabel=ica2yeo7.Yeo7N;myLabel=myLabel([1 7 9 13 14 20 21 5 6 16 10 11 12 17 3 2 4 8 19 15 18]);
myLabel=repmat({''}, [21,1]); figure
upperlim=4; lowerlim=-4; %upperlim=nodes(205,component); lowerlim=nodes(5,component);
Weights=dep_vs_anx(tmp); Weights(Weights<upperlim & Weights>lowerlim)=0;%Weights=stats.W(:,component); 
Weights(Weights>0)=0; Weights_square= zeros(21); Weights_square(triu(ones(21),1)>0) = abs(Weights);
hold on; myColorMap=0.5*ones(21,3);myColorMap(:,3)=1; circularGraph(Weights_square, 'Colormap',myColorMap, 'Label',myLabel);
Weights=dep_vs_anx(tmp); Weights(Weights<upperlim & Weights>lowerlim)=0;
Weights(Weights<0)=0; Weights_square= zeros(21); Weights_square(triu(ones(21),1)>0) = abs(Weights);
hold on; myColorMap=0.5*ones(21,3);myColorMap(:,3)=1; circularGraph(Weights_square, 'Colormap',myColorMap, 'Label',myLabel);

myLabel=ica2yeo7.Yeo7N;myLabel=myLabel([1 7 9 13 14 20 21 5 6 16 10 11 12 17 3 2 4 8 19 15 18]);
myLabel=repmat({''}, [21,1]); hold on
upperlim=1.96; lowerlim=-1.96; %upperlim=nodes(205,component); lowerlim=nodes(5,component);
Weights=DepVsAnx(tmp); Weights(Weights<upperlim & Weights>lowerlim)=0;%Weights=stats.W(:,component); 
Weights(Weights<0)=0; Weights_square= zeros(21); Weights_square(triu(ones(21),1)>0) = abs(Weights);
hold on; myColorMap=0.7*ones(21,3);%myColorMap(:,3)=1; 
circularGraph(Weights_square, 'Colormap',myColorMap, 'Label',myLabel);
Weights=DepVsAnx(tmp); Weights(Weights<upperlim & Weights>lowerlim)=0;
Weights(Weights>0)=0; Weights_square= zeros(21); Weights_square(triu(ones(21),1)>0) = abs(Weights);
hold on; myColorMap=0.7*ones(21,3);%myColorMap(:,1)=1; 
circularGraph(Weights_square, 'Colormap',myColorMap, 'Label',myLabel);

myLabel=ica2yeo7.Yeo7N;myLabel=myLabel([1 7 9 13 14 20 21 5 6 16 10 11 12 17 3 2 4 8 19 15 18]);
myLabel=repmat({''}, [21,1]); hold on
upperlim=1.96; lowerlim=-1.96; %upperlim=nodes(205,component); lowerlim=nodes(5,component);
Weights=d(tmp); Weights(Weights<upperlim & Weights>lowerlim)=0;%Weights=stats.W(:,component); 
Weights(Weights<0)=0; Weights_square= zeros(21); Weights_square(triu(ones(21),1)>0) = abs(Weights);
hold on; myColorMap=zeros(21,3);%myColorMap(:,3)=1; 
circularGraph(Weights_square, 'Colormap',myColorMap, 'Label',myLabel);
Weights=d(tmp); Weights(Weights<upperlim & Weights>lowerlim)=0;
Weights(Weights>0)=0; Weights_square= zeros(21); Weights_square(triu(ones(21),1)>0) = abs(Weights);
hold on; myColorMap=zeros(21,3);%myColorMap(:,1)=1; 
circularGraph(Weights_square, 'Colormap',myColorMap, 'Label',myLabel);


%%%% example visuals
ICs={'IC1';'IC2';'IC3';'IC4';'IC5';'IC6';'IC7';'IC8';'IC9';'IC10';'IC11';'IC12';'IC13';'IC14';'IC15';'IC16';'IC17';'IC18';'IC19';'IC20';'IC21'};
combs=allcomb(ICs, ICs);D=21
for i=1:441; comb(i)=cellstr(strcat(combs{i,2}, combs{i,1})); end
ICs=reshape(comb, [21,21]); clear comb combs i
ICs_vector=ICs(triu(ones(D),1)==1)';

roi=9; ICs_vector(roi)
load('pls_dep.mat'); plsweights1(roi)
[~,ypca,~] = pca(Y);ypca=ypca(:,1); 
figure; scatter(X(:,roi), ypca, 2, 'k', 'filled'); mdl= fitlm(X(:,roi), ypca); [ypred,yci] = predict(mdl,X(:,roi), 'Alpha',0.001); hold on; plot(X(:,roi), ypred, 'k', 'LineWidth', 2);%plot(X(:,roi), yci(:,1), 'k', 'LineWidth', 0.5);plot(X(:,roi), yci(:,2), 'k', 'LineWidth', 0.5);
corr(X(:,roi), ypca)
load('pls_anx.mat');plsweights1(roi)
[~,ypca,~] = pca(Y);ypca=ypca(:,1); 
figure; scatter(X(:,roi), ypca, 2, 'k', 'filled'); mdl= fitlm(X(:,roi), ypca); [ypred,yci] = predict(mdl,X(:,roi), 'Alpha',0.001); hold on; plot(X(:,roi), ypred, 'k', 'LineWidth', 2);%plot(X(:,roi), yci(:,1), 'k', 'LineWidth', 0.5);plot(X(:,roi), yci(:,2), 'k', 'LineWidth', 0.5);
corr(X(:,roi), ypca)

%% train on MDD and predict ANX
load('pls_dep.mat');
beta_dep=BETA;
load('pls_anx.mat');
y_pred = [ones(size(X,1),1) X]*BETA;
r=corr(Y, y_pred)
y_pred = [ones(size(X,1),1) X]*beta_dep;
r=corr(Y, y_pred)
r(eye(4)==1)

%% In MDD clinical featurs vs XS
PHQ2=(minimal.x2050_2_0 + minimal.x2060_2_0)>=2 ; sum(PHQ2>0)
anova1(XS(:,1), PHQ2)
corr(XS(:,1), minimal.x20433_0_0, 'rows', 'pairwise')
corr(XS(:,1), cognitive.x20127_0_0, 'rows', 'pairwise')

%% Hold out data
%%
%% Hold out data in all cases and in all HC
l=length(X(:,1)); 
%randomsample=randperm(l); X=X(randomsample,:);Y=Y(randomsample,:);
l=round(0.25*(l),0); 

for h=1:4
xstart=(h-1)*l+1;
xend=h*l;
ixtest=zeros(1,length(X(:,1)));
ixtest(xstart:xend)=1;
X_train=x(ixtest==0,:); X_test=x(ixtest==1 ,:);
Y_train=Y(ixtest==0,:); Y_test=Y(ixtest==1 ,:);
%l=length(X(:,1)); l=round(0.25*l,0); X_test=x(1:l,:); X_train=x(l+1:length(X(:,1)) ,:);Y_test=Y(1:l,:); Y_train=Y(l+1:length(X(:,1)) ,:);
dim=3;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(X_train,Y_train,dim);PCTVAR;
r=corr(XS, Y_train); accuracytrain(h,:)=r(1,:);
y_pred = [ones(size(X_test ,1),1) X_test]*BETA;
r=corr(Y_test, y_pred);accuracyholdout(h,:)=r(eye(4)==1);
end
mean(accuracyholdout)

figure(11); hold off; for i=1:4; subplot(4,1,i); scatter(y_pred(:,i), Y_test(:,i), 'x', 'k'); 
    p = polyfit(y_pred(:,i),Y_test(:,i),1); pred = polyval(p,y_pred(:,i)); hold on; plot(y_pred(:,i),pred,'r','LineWidth',3); set(gca,'xtick',[]); set(gca,'ytick',[]); end





















%% visualize the weights using circularGraph(x);
component=1
nodes=dlmread('icad25_nodesordered.csv'); tmp=nanmean(ica_partial(id_healthy==0,:))';
nodes=[plsweights1,plsweights2, (1:210)',tmp, nodes]; nodes=sortrows(nodes, component);

myLabel=ica2yeo7.Yeo7N;figure
upperlim=3; lowerlim=-3; %upperlim=nodes(205,component); lowerlim=nodes(5,component);
Weights=plsweights1; Weights(Weights<upperlim & Weights>lowerlim)=0;%Weights=stats.W(:,component); 
Weights(Weights<0)=0; Weights_square= zeros(21); Weights_square(triu(ones(21),1)>0) = abs(Weights);
hold on; myColorMap=zeros(21,3);myColorMap(:,3)=1; circularGraph(Weights_square, 'Colormap',myColorMap, 'Label',myLabel);
Weights=plsweights1; Weights(Weights<upperlim & Weights>lowerlim)=0;
Weights(Weights>0)=0; Weights_square= zeros(21); Weights_square(triu(ones(21),1)>0) = abs(Weights);
hold on; myColorMap=zeros(21,3);myColorMap(:,1)=1; circularGraph(Weights_square, 'Colormap',myColorMap, 'Label',myLabel);

%%% average connectivity visual
Weights=tmp;
Weights(Weights>-0.4)=0; Weights_square= zeros(21); Weights_square(triu(ones(21),1)>0) = abs(Weights);
figure; myColorMap=zeros(21,3);myColorMap(:,3)=1; circularGraph(Weights_square, 'Colormap',myColorMap, 'Label',myLabel);
Weights=tmp;
Weights(Weights<0.4)=0; Weights_square= zeros(21); Weights_square(triu(ones(21),1)>0) = abs(Weights);
hold on; myColorMap=zeros(21,3);myColorMap(:,1)=1; circularGraph(Weights_square, 'Colormap',myColorMap, 'Label',myLabel);

%% alternative visuals
D=21 
ICs={'IC1';'IC7';'IC9';'IC13';'IC14';'IC20';'IC21';'IC5';'IC6';'IC16';'IC10';'IC11';'IC12';'IC17';'IC3';'IC2';'IC4';'IC8';'IC19';'IC15';'IC18'};
combs=allcomb(ICs, ICs);
for i=1:441; comb(i)=cellstr(strcat(combs{i,2}, combs{i,1})); end
ICs=reshape(comb, [21,21]); clear comb combs i
ICs_vector=ICs(triu(ones(D),1)==1)';

ic_mapping=readtable('D:\Canada_2020\OASIS\reports\fmri\ic_order_mapping.csv');
tmp=ic_mapping.original_order_num;

plsweights1=plsweights1(tmp);


Tdep_square = zeros(21); Tdep_square(triu(ones(21),1)>0) = plsweights1; %Tanx_square(Tanx_square==0)=NaN;
Panova_square=zeros(21);Panova_square(triu(ones(21),1)>0) = plsweights1; %= p_perm_mddprs; 
tmp=Tdep_square; tmp(abs(Panova_square)<3)=0;
figure(22); imagesc(tmp); colorbar; caxis([-5 5]); 
tmp2=ica2yeo7.Yeo7N([1 7 9 13 14 20 21 5 6 16 10 11 12 17 3 2 4 8 19 15 18]); set(gca, 'XTick', 1:21, 'XTickLabel', tmp2 , 'XTickLabelRotation',90); set(gca, 'YTick', 1:21, 'YTickLabel', tmp2);





%% PRS associations

corr(XS(:,1), YS(:,1))
corr(MDD_prs, Y, 'rows','complete')
scatter(clean(MDD_prs_r,3), X(:,1) )

[b,bint,r] = regress(MDD_prs,[ones(length(Y),1), minimal.x21003_2_0, minimal.x31_0_0, minimal.x54_2_0 ] ); MDD_prs_r=r;

T=table(MDD_prs_r, c', XS(:,2));
mdl = fitlm(T,'Var3~Var2+MDD_prs_r')
PHQ2=(minimal.x2050_2_0 + minimal.x2060_2_0)>=2 ; sum(PHQ2>0)
anova1(XS(:,1), PHQ2)









%% regress out sex and age
ix=(cellfun('isempty', clinical)'| isnan(ica_partial(:,1))) ; a=age; a(ix)=[];s=sex;s(ix)=[];clin=clinical; clin(ix)=[];
cognitive=cognitive_ordered; cognitive(ix,:)=[];

%X=thicknessFS; X(ix,:)=[]; 
X=ica_partial; X(ix,:)=[]; 
%%% what are the y variables: symbol digit, alphanumtrails, tower, pairs match, word pairs (PAL)
Y=horzcat(a, s,cognitive.x23324_2_0);% , cognitive.x6350_2_0, cognitive.x21004_2_0); %, cognitive.x399_2_2+cognitive.x399_2_1, cognitive.x20197_2_0
%%% no controls? remove them here
index=~(  ismember(clin','hc') | ismember(clin', 'str')  ); %index=~ismember(clin','hc'); %
Y=Y(index,:); X=X(index,:);
%%% remove nans
naninx=sum(isnan(Y)')'; Y=Y(naninx==0,:);  X=X(naninx==0,:);
%Y=zscore(Y); X=zscore(X);

clear x y
[b,bint,r] = regress(Y(:,3),[ones(length(Y(:,3)),1), Y(:,1:2)] ); y(:,1)=r;
for i=1:length(X(1,:))
   [b,bint,r] = regress(X(:,i),[ones(length(Y(:,3)),1), Y(:,1:2)] ); x(:,i)=r; 
end


ncomp=4
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(X,Y,ncomp);PCTVAR
corr(XS, Y)
figure; scatter(XS(:,2), y(:,2))

%% norming the scores
cognitive_hc=cognitive(  ismember(clin','hc'),: );
minimal=minimal_ordered; minimal(ix,:)=[];
minimal_hc=minimal(  ismember(clin','hc'),: );
% personalized scores
for i=1:length(Y(:,1))
ss=Y(i,2);
aa=Y(i,1);
idcontrol=(minimal_hc.x31_0_0==ss & minimal_hc.x21003_2_0<=(aa+2) & minimal_hc.x21003_2_0>=(aa-2) );
Y_control=cognitive_hc(idcontrol,:);
    Yz(i,1)=(Y(i,3)-nanmean( Y_control.x23324_2_0))/nanstd( Y_control.x23324_2_0 );
    Yz(i,2)=(Y(i,4)-nanmean( Y_control.x6350_2_0 ))/nanstd( Y_control.x6350_2_0 );
    Yz(i,3)=(Y(i,5)-nanmean( Y_control.x21004_2_0 ))/nanstd( Y_control.x21004_2_0 );
    %figure;histogram(Yz);
end


%% ALTERNATIVE COMPARING PLS
[~,ypca,~] = pca(Y);
ypca=[ypca(:,1),ypca(:,1)]; 
ypca(strcmp(c, 'dep'),1)=0;
ypca(strcmp(c, 'anx'),2)=0;
Y=ypca;
