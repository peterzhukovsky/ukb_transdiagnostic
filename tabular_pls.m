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
index=~ismember(clin','hc') & ~ismember(clin', 'dep') & ~ismember(clin', 'depanx') & ~ismember(clin', 'str'); %  index=~(  ismember(clin','hc')   ); %      | ismember(clin', 'str')
Y=Y(index,:); X=X(index,:); minimal=minimal(index,:); c=clin(index); MDD_prs=MDD_prs(index);
%%% remove nans
naninx=sum(isnan(Y)')'; Y=Y(naninx==0,:);  X=X(naninx==0,:); minimal=minimal(naninx==0,:);c=c(naninx==0);MDD_prs=MDD_prs(naninx==0);
Y=zscore(Y); X=zscore(X);
ncomp=3
%%% regress out sex and age and potentially other covariates
mdl = fitlm(table(minimal.x21003_2_0, minimal.x31_0_0, minimal.x54_2_0,Y(:,1)) ); Y(:,1)=mdl.Residuals.Raw;
mdl = fitlm(table(minimal.x21003_2_0, minimal.x31_0_0, minimal.x54_2_0,Y(:,2)) ); Y(:,2)=mdl.Residuals.Raw;
mdl = fitlm(table(minimal.x21003_2_0, minimal.x31_0_0, minimal.x54_2_0,Y(:,3)) ); Y(:,3)=mdl.Residuals.Raw;
mdl = fitlm(table(minimal.x21003_2_0, minimal.x31_0_0, minimal.x54_2_0,Y(:,4)) ); Y(:,4)=mdl.Residuals.Raw;
clear x; for i=1:length(X(1,:)); mdl = fitlm(table(minimal.x21003_2_0, minimal.x31_0_0, minimal.x25741_2_0, minimal.x54_2_0,X(:,i)) ); x(:,i)=mdl.Residuals.Raw;end

[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(x,Y,ncomp);PCTVAR
%%% explore the PLS components and their correlation with Y
figure(1);imagesc(corr(XS,Y)); corr(XS, Y) %1.TMT  2.Tower  3.PAL  4.DSST
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

%% alternative visuals

Tdep_square = zeros(21); Tdep_square(triu(ones(21),1)>0) = plsweights1; %Tanx_square(Tanx_square==0)=NaN;
Panova_square=zeros(21);Panova_square(triu(ones(21),1)>0) = plsweights1; %= p_perm_mddprs; 
tmp=Tdep_square; tmp(abs(Panova_square)<3)=0;
figure(22); imagesc(tmp); colorbar; caxis([-5 5]); set(gca, 'XTick', 1:21, 'XTickLabel', ica2yeo7.Yeo7N, 'XTickLabelRotation',90);
 set(gca, 'YTick', 1:21, 'YTickLabel', ica2yeo7.Yeo7N);

