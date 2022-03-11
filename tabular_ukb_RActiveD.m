addpath('D:\Canada_2020\UK_biobank\reports');
%% permutation testing GLMs
clear; dataimport
permutations=1000;
ica_partial=dlmread('D:\Canada_2020\UK_biobank\data\ica_d25_par.csv');
PHQ2=(minimal_ordered.x2050_2_0 + minimal_ordered.x2060_2_0)>=2 ;
ix=cellfun('isempty', clinical)'| isnan(ica_partial(:,1)) | (strcmp(clinical, 'dep')' & PHQ2==0);
ica_partial(ix,:)=[]; a=age; a(ix)=[];s=sex;s(ix)=[];clin=clinical; clin(ix)=[];
cognitive=cognitive_ordered; cognitive(ix,:)=[]; minimal=minimal_ordered; minimal(ix,:)=[]; 


%% subset controls
%%%%% matching cases to controls 
case_ids=find(~strcmp(clin, 'ahc')); to_keep=zeros(length(case_ids),1);
for i=1:sum(~strcmp(clin, 'ahc'))
   age_case=a(case_ids(i));sex_case=s(case_ids(i));
   ix_to_keep=find(strcmp(clin, 'ahc')' & a==age_case & s==sex_case);
   ix_to_keep=ix_to_keep( randperm(length(ix_to_keep)) ); 
   for j=1:length(ix_to_keep) % this ensures we dont have any duplicate matched controls, ie each person with same age and sex gets a new match
       tmp(j)=  sum(ix_to_keep(j)==to_keep)>0;
   end
   ix_to_keep(tmp>0)=[]; to_keep(i)=ix_to_keep(1); clear tmp
end; clear ix_to_keep case_ids age_case sex_case
ix=1:length(a); ix(  vertcat(to_keep, find(~strcmp(clin, 'ahc'))')  )=[];
length(a)-length(ix)

ica_partial(ix,:)=[]; a(ix)=[];s(ix)=[]; clin(ix)=[]; cognitive=cognitive_ordered; minimal(ix,:)=[];  cognitive(ix,:)=[]; % cutting data
1-sum(s(strcmp(clin, 'ahc')))/3149 % percent female
figure; histogram(a(strcmp(clin, 'ahc'))); hold on; histogram(a(~strcmp(clin, 'ahc'))) % age histograms

for i=1:length(ica_partial(1,:))
        allobservations=ica_partial(:,i);i
        parfor n = 1:permutations; 
        permutation_index = randperm(length(allobservations));
        randomSample = allobservations(permutation_index,:);
        T=mktbl(a,s,clin',  minimal.x25741_2_0,  minimal.x54_2_0, randomSample); T=sortrows(T, 3);
        mdl = fitlm(T,'y~x1+x1^2+x1*x2+x2+x3+x4+x5');Td(1,n)=mdl.Coefficients.tStat(5); Tda(1,n)=mdl.Coefficients.tStat(6); Ts(1,n)=mdl.Coefficients.tStat(7); Ta(1,n)=mdl.Coefficients.tStat(4);
        anovamdl=anova(mdl); F(1,n)=anovamdl.F(3);
        end; 
        F_perm(i,:)=F; Tperm_dep(i,:)=Td;Tperm_depan(i,:)=Tda;Tperm_str(i,:)=Ts;Tperm_anx(i,:)=Ta;
        clear Td Tda Ts Ta Tmddp F2 F3 F4
  clear T; T=mktbl(a,s,clin', minimal.x25741_2_0, minimal.x54_2_0,  allobservations);T=sortrows(T, 3);
  mdl = fitlm(T,'y~x1+x1^2+x1*x2+x2+x3+x4+x5'); anovamdl=anova(mdl);
  p_perm_anova(i)=1-sum(anovamdl.F(3)>F_perm(i,:))/permutations; 
  Tage(1,i)=mdl.Coefficients.tStat(2); 
  Tsex(1,i)=mdl.Coefficients.tStat(3);
  Tdep(1,i)=mdl.Coefficients.tStat(5); Pdep(1,i)=mdl.Coefficients.pValue(5);
  Tdepanx(1,i)=mdl.Coefficients.tStat(6); Pdepanx(1,i)=mdl.Coefficients.pValue(6);
  Tstr(1,i)=mdl.Coefficients.tStat(7); Pstr(1,i)=mdl.Coefficients.pValue(7);
  Tanx(1,i)=mdl.Coefficients.tStat(4); Panx(1,i)=mdl.Coefficients.pValue(4);
end
% effect sizes
2*5/sqrt(6298) 
%%% sign testing correlations - potential covariates: 1) simple 2) with three prs scores and 3) Years of education, income, urban vs rural environment, marital status, stressful life events and substance use family history variables for psychiatric disease
r=corr(Tperm_dep, Tperm_anx); r=r(eye(permutations)==1); P(1,1)=1-sum(corr(Tdep', Tanx')>r)/permutations; % column1 for dep vs anx
r=corr(Tperm_dep, Tperm_str); r=r(eye(permutations)==1); P(1,3)=1-sum(corr(Tdep', Tstr')>r)/permutations; % column1 for dep vs anx
r=corr(Tperm_str, Tperm_anx); r=r(eye(permutations)==1); P(2,3)=1-sum(corr(Tstr', Tanx')>r)/permutations; % column1 for dep vs anx
r=corr(Tperm_depan, Tperm_dep); r=r(eye(permutations)==1); P(1,2)=1-sum(corr(Tdepanx', Tdep')>r)/permutations; % column1 for dep vs anx
r=corr(Tperm_depan, Tperm_anx); r=r(eye(permutations)==1); P(2,2)=1-sum(corr(Tdepanx', Tanx')>r)/permutations; % column1 for dep vs anx
r=corr(Tperm_depan, Tperm_str); r=r(eye(permutations)==1); P(3,3)=1-sum(corr(Tdepanx', Tstr')>r)/permutations; % column1 for dep vs anx

corr([Tdep', Tanx',Tdepanx',Tstr']) %, -1*Tmci'
                        %% visuals
figure(3);subplot(2,3,1); scatter(Tdep', Tanx','.','k'); ylim([-6 6]); corr(Tdep', Tanx')
P = polyfit(Tdep', Tanx',1); x0 = min(Tdep') ; x1 = max(Tdep') ;xi = linspace(x0,x1) ; yi = P(1)*xi+P(2); hold on; plot(xi,yi,'k') ;
subplot(2,3,2);scatter(Tdep', Tstr','.','k');ylim([-6 6]);  corr(Tdep', Tstr')
P = polyfit(Tdep', Tstr',1); x0 = min(Tdep') ; x1 = max(Tdep') ;xi = linspace(x0,x1) ; yi = P(1)*xi+P(2); hold on; plot(xi,yi,'k') ;
subplot(2,3,3);scatter(Tanx', Tstr','.','k'); ylim([-6 6]); corr(Tanx', Tstr')
P = polyfit(Tanx', Tstr',1); x0 = min(Tanx') ; x1 = max(Tanx') ;xi = linspace(x0,x1) ; yi = P(1)*xi+P(2); hold on; plot(xi,yi,'k') ;
subplot(2,3,4); scatter(Tdep', Tdepanx','.','k'); ylim([-6 6]);corr(Tdep', Tdepanx')
P = polyfit(Tdep', Tdepanx',1); x0 = min(Tdep') ; x1 = max(Tdep') ;xi = linspace(x0,x1) ; yi = P(1)*xi+P(2); hold on; plot(xi,yi,'k') ;
subplot(2,3,5); scatter(Tanx', Tdepanx','.','k'); ylim([-6 6]);corr(Tanx', Tdepanx')
P = polyfit(Tanx', Tdepanx',1); x0 = min(Tanx') ; x1 = max(Tanx') ;xi = linspace(x0,x1) ; yi = P(1)*xi+P(2); hold on; plot(xi,yi,'k') ;
subplot(2,3,6); scatter(Tstr', Tdepanx','.','k'); ylim([-6 6]); corr(Tstr', Tdepanx')
P = polyfit(Tstr', Tdepanx',1); x0 = min(Tstr') ; x1 = max(Tstr') ;xi = linspace(x0,x1) ; yi = P(1)*xi+P(2); hold on; plot(xi,yi,'k') ;


%% %% %% %%%% %%%% %%%% %%%% %%
%% depression vs anxiety vs stress -  structural
%% %% %% %%%% %%%% %%%% %%%% %%
clear; dataimport
ct_ordered=dlmread('D:\Canada_2020\UK_biobank\reports\ordered\ct_ordered.csv');
permutations=1000;
PHQ2=(minimal_ordered.x2050_2_0 + minimal_ordered.x2060_2_0)>=2 ;
ix=cellfun('isempty', clinical)' | sum(isnan(ct_ordered)')'>0| (strcmp(clinical, 'dep')' & PHQ2==0);
ct_ordered(ix,:)=[]; a=age; a(ix)=[];s=sex;s(ix)=[];clin=clinical; clin(ix)=[];
cognitive=cognitive_ordered; cognitive(ix,:)=[]; minimal=minimal_ordered; minimal(ix,:)=[];  MDD_prs=MDD_prs_ordered; MDD_prs(ix)=[]; ANX_prs=ANX_prs_ordered; ANX_prs(ix)=[]; PTSD_prs=PTSD_prs_ordered; PTSD_prs(ix)=[];

%%%% subset controls
%%%%% reviewers suggested matching cases to controls 
case_ids=find(~strcmp(clin, 'ahc')); to_keep=zeros(length(case_ids),1);
for i=1:sum(~strcmp(clin, 'ahc'))
   age_case=a(case_ids(i));sex_case=s(case_ids(i));
   ix_to_keep=find(strcmp(clin, 'ahc')' & a==age_case & s==sex_case);
   ix_to_keep=ix_to_keep( randperm(length(ix_to_keep)) ); 
   for j=1:length(ix_to_keep) % this ensures we dont have any duplicate matched controls, ie each person with same age and sex gets a new match
       tmp(j)=  sum(ix_to_keep(j)==to_keep)>0;
   end
   ix_to_keep(tmp>0)=[]; to_keep(i)=ix_to_keep(1); clear tmp
end; clear ix_to_keep case_ids age_case sex_case
ix=1:length(a); ix(  vertcat(to_keep, find(~strcmp(clin, 'ahc'))')  )=[];
length(a)-length(ix)

ct_ordered(ix,:)=[]; a(ix)=[];s(ix)=[]; clin(ix)=[]; cognitive=cognitive_ordered; minimal(ix,:)=[];  cognitive(ix,:)=[]; 
1-sum(s(strcmp(clin, 'ahc')))/3570
figure; histogram(a(strcmp(clin, 'ahc'))); hold on; histogram(a(~strcmp(clin, 'ahc')))

for i=1:length(ct_ordered(1,:))
        allobservations=ct_ordered(:,i);i
        parfor n = 1:permutations; 
        permutation_index = randperm(length(allobservations));
        randomSample = allobservations(permutation_index,:);
        T=mktblct(a,s,clin',    minimal.x54_2_0, randomSample);T=sortrows(T, 3);
        mdl = fitlm(T,'y~x1+x1^2+x1*x2+x2+x3+x4 ');
        anovamdl=anova(mdl); F(1,n)=anovamdl.F(3);
        end; 
        F_perm(i,:)=F; 
  clear T; T=mktblct(a,s,clin', minimal.x54_2_0, allobservations);T=sortrows(T, 3);
  mdl = fitlm(T,'y~x1+x1^2+x1*x2+x2+x3+x4'); anovamdl=anova(mdl);
  p_perm_anova(i)=1-sum(anovamdl.F(3)>F_perm(i,:))/permutations; 
  Tage(1,i)=mdl.Coefficients.tStat(2); 
  Tsex(1,i)=mdl.Coefficients.tStat(3);
  Tdep(1,i)=mdl.Coefficients.tStat(5); Pdep(1,i)=mdl.Coefficients.pValue(5);
  Tdepanx(1,i)=mdl.Coefficients.tStat(6); Pdepanx(1,i)=mdl.Coefficients.pValue(6);
  Tstr(1,i)=mdl.Coefficients.tStat(7); Pstr(1,i)=mdl.Coefficients.pValue(7);
  Tanx(1,i)=mdl.Coefficients.tStat(4); Panx(1,i)=mdl.Coefficients.pValue(4);
end
%figure;imagesc(Tanx);colorbar; colormap default

%%%sign regions:
tmp=Pdepanx; tmp2=Tdepanx;
table(label_names_all(p_perm_anova<0.05 & tmp<0.0125),...
label_names_HOA(p_perm_anova<0.05 & tmp<0.0125), ...
tmp2(p_perm_anova<0.05 & tmp<0.0125)', tmp(p_perm_anova<0.05 & tmp<0.0125)',...
2*tmp2(p_perm_anova<0.05 & tmp<0.0125)'./sqrt(12214)    ); clear tmp tmp2; sortrows(ans, 2)


figure(3);subplot(2,3,1); scatter(Tdep', Tanx'); corr(Tdep', Tanx')
subplot(2,3,2);scatter(Tdep', Tstr'); corr(Tdep', Tstr')
subplot(2,3,3);scatter(Tanx', Tstr'); corr(Tanx', Tstr')
subplot(2,3,4); scatter(Tdep', Tdepanx'); corr(Tdep', Tdepanx')
subplot(2,3,5); scatter(Tanx', Tdepanx'); corr(Tanx', Tdepanx')
subplot(2,3,6); scatter(Tstr', Tdepanx'); corr(Tstr', Tdepanx')

corr([Tdep', Tanx',Tdepanx',Tstr']) %, -1*Tmci'

%% Remapping to a better ordering
D=21 
ICs={'IC1';'IC7';'IC9';'IC13';'IC14';'IC20';'IC21';'IC5';'IC6';'IC16';'IC10';'IC11';'IC12';'IC17';'IC3';'IC2';'IC4';'IC8';'IC19';'IC15';'IC18'};
combs=allcomb(ICs, ICs);
for i=1:441; comb(i)=cellstr(strcat(combs{i,2}, combs{i,1})); end
ICs=reshape(comb, [21,21]); clear comb combs i
ICs_vector=ICs(triu(ones(D),1)==1)';

ic_mapping=readtable('D:\Canada_2020\OASIS\reports\fmri\ic_order_mapping.csv');
tmp=ic_mapping.original_order_num;

Tanx=Tanx(tmp);Tdep=Tdep(tmp);Tstr=Tstr(tmp);Tdepanx=Tdepanx(tmp); p_perm_anova=p_perm_anova(tmp);
Panx=Panx(tmp);Pdep=Pdep(tmp);Pstr=Pstr(tmp);Pdepanx=Pdepanx(tmp); 


Tdep_square = zeros(21); Tdep_square(triu(ones(21),1)>0) = Tdep; %Tanx_square(Tanx_square==0)=NaN;
Panova_square=zeros(21);Panova_square(triu(ones(21),1)>0) = p_perm_anova; %Panx_square(Panx_square==0)=NaN;
Pdep_square=zeros(21);Pdep_square(triu(ones(21),1)>0) = Pdep; Pdep_square(Pdep_square==0)=NaN;
tmp=Tdep_square; tmp(Panova_square>0.05 | Pdep_square>0.0125)=0;
tmp2=ica2yeo7.Yeo7N([1 7 9 13 14 20 21 5 6 16 10 11 12 17 3 2 4 8 19 15 18]);
figure(2); imagesc(tmp); colorbar; caxis([-5 5]); set(gca, 'XTick', 1:21, 'XTickLabel', tmp2 , 'XTickLabelRotation',90);
set(gca, 'YTick', 1:21, 'YTickLabel', tmp2);
