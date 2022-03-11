%% permutation testing GLMs
clear; dataimport
permutations=1000;
ica_partial=dlmread('D:\Canada_2020\UK_biobank\data\ica_d25_par.csv');
ix=cellfun('isempty', clinical)'| isnan(ica_partial(:,1)) |isnan(whitebritish_ordered); %| ~(  ismember(minimal_ordered.x21000_0_0,'White')   ) ; 
a=age; a(ix)=[];s=sex;s(ix)=[];clin=clinical; clin(ix)=[];
ica_partial(ix,:)=[];cognitive=cognitive_ordered; cognitive(ix,:)=[]; minimal=minimal_ordered; minimal(ix,:)=[]; MDD_prs=MDD_prs_ordered; MDD_prs(ix)=[]; ANX_prs=ANX_prs_ordered; ANX_prs(ix)=[]; PTSD_prs=PTSD_prs_ordered; PTSD_prs(ix)=[];whitebritish=whitebritish_ordered;whitebritish(ix)=[];


%% subset controls
%%%%% matching cases to controls 
case_ids=find(~strcmp(clin, 'ahc')); to_keep=zeros(length(case_ids),1);
for i=1:sum(~strcmp(clin, 'ahc'))
   age_case=a(case_ids(i));sex_case=s(case_ids(i));
   ix_to_keep=find(strcmp(clin, 'ahc')' & a==age_case & s==sex_case);
   ix_to_keep=ix_to_keep( randperm(length(ix_to_keep)) ); 
   for j=1:length(ix_to_keep)
       tmp(j)=  sum(ix_to_keep(j)==to_keep)>0;
   end
   ix_to_keep(tmp>0)=[]; to_keep(i)=ix_to_keep(1); clear tmp
end; clear ix_to_keep case_ids age_case sex_case
ix=1:length(a); ix(  vertcat(to_keep, find(~strcmp(clin, 'ahc'))')  )=[];
length(a)-length(ix)

ica_partial(ix,:)=[]; a(ix)=[];s(ix)=[]; clin(ix)=[]; minimal(ix,:)=[]; % cutting data
MDD_prs(ix)=[]; ANX_prs(ix)=[]; PTSD_prs(ix)=[];whitebritish(ix)=[];
1-sum(s(strcmp(clin, 'ahc')))/sum(~strcmp(clin, 'ahc')) % percent female
figure; histogram(a(strcmp(clin, 'ahc'))); hold on; histogram(a(~strcmp(clin, 'ahc'))) % age histograms

for i=1:length(ica_partial(1,:))
        allobservations=ica_partial(:,i);i
        parfor n = 1:permutations; 
        permutation_index = randperm(length(allobservations));
        randomSample = allobservations(permutation_index,:);
        T=mktbl(a,s,clin',  minimal.x25741_2_0,  minimal.x54_2_0, MDD_prs,ANX_prs,PTSD_prs, randomSample); T=sortrows(T, 3);
        mdl = fitlm(T,'y~x1+x1^2+x1*x2+x2+x3+x4+x5+x6+x7+x8');Td(1,n)=mdl.Coefficients.tStat(5); Tda(1,n)=mdl.Coefficients.tStat(6); Ts(1,n)=mdl.Coefficients.tStat(7); Ta(1,n)=mdl.Coefficients.tStat(4);
        anovamdl=anova(mdl); F(1,n)=anovamdl.F(3);
        end; 
        F_perm(i,:)=F; Tperm_dep(i,:)=Td;Tperm_depan(i,:)=Tda;Tperm_str(i,:)=Ts;Tperm_anx(i,:)=Ta;
        clear Td Tda Ts Ta Tmddp F2 F3 F4
  clear T; T=mktbl(a,s,clin', minimal.x25741_2_0, minimal.x54_2_0, MDD_prs,ANX_prs,PTSD_prs, allobservations);T=sortrows(T, 3);
  mdl = fitlm(T,'y~x1+x1^2+x1*x2+x2+x3+x4+x5+x6+x7+x8'); anovamdl=anova(mdl);
  p_perm_anova(i)=1-sum(anovamdl.F(3)>F_perm(i,:))/permutations; 
  Tage(1,i)=mdl.Coefficients.tStat(2); 
  Tsex(1,i)=mdl.Coefficients.tStat(3);
  Tdep(1,i)=mdl.Coefficients.tStat(5); Pdep(1,i)=mdl.Coefficients.pValue(5);
  Tdepanx(1,i)=mdl.Coefficients.tStat(6); Pdepanx(1,i)=mdl.Coefficients.pValue(6);
  Tstr(1,i)=mdl.Coefficients.tStat(7); Pstr(1,i)=mdl.Coefficients.pValue(7);
  Tanx(1,i)=mdl.Coefficients.tStat(4); Panx(1,i)=mdl.Coefficients.pValue(4);
end
% effect sizes
2*5/sqrt(7602) 
%%% sign testing correlations - potential covariates: 1) simple 2) with three prs scores and 3) Years of education, income, urban vs rural environment, marital status, stressful life events and substance use family history variables for psychiatric disease
r=corr(Tperm_dep, Tperm_anx); r=r(eye(permutations)==1); P(1,1)=1-sum(corr(Tdep', Tanx')>r)/permutations; % column1 for dep vs anx
r=corr(Tperm_dep, Tperm_str); r=r(eye(permutations)==1); P(1,3)=1-sum(corr(Tdep', Tstr')>r)/permutations; % column1 for dep vs anx
r=corr(Tperm_str, Tperm_anx); r=r(eye(permutations)==1); P(2,3)=1-sum(corr(Tstr', Tanx')>r)/permutations; % column1 for dep vs anx
r=corr(Tperm_depan, Tperm_dep); r=r(eye(permutations)==1); P(1,2)=1-sum(corr(Tdepanx', Tdep')>r)/permutations; % column1 for dep vs anx
r=corr(Tperm_depan, Tperm_anx); r=r(eye(permutations)==1); P(2,2)=1-sum(corr(Tdepanx', Tanx')>r)/permutations; % column1 for dep vs anx
r=corr(Tperm_depan, Tperm_str); r=r(eye(permutations)==1); P(3,3)=1-sum(corr(Tdepanx', Tstr')>r)/permutations; % column1 for dep vs anx

corr([Tdep', Tanx',Tdepanx',Tstr']) %, -1*Tmci'
                        %% visuals
%index=ismember(clin,'hc'); figure(1);histogram(a(index)); hold on;yyaxis right ; %histogram(a(~index)); clear index
Tdep_square = zeros(21); Tdep_square(triu(ones(21),1)>0) = Tdep; %Tanx_square(Tanx_square==0)=NaN;
Panova_square=zeros(21);Panova_square(triu(ones(21),1)>0) = p_perm_anova; %Panx_square(Panx_square==0)=NaN;
Pdep_square=zeros(21);Pdep_square(triu(ones(21),1)>0) = Pdep; Pdep_square(Pdep_square==0)=NaN;
tmp=Tdep_square; tmp(Panova_square>0.05 | Pdep_square>0.0125)=0;
figure; imagesc(tmp); colorbar; caxis([-5 5]); set(gca, 'XTick', 1:21, 'XTickLabel', ica2yeo7.Yeo7N, 'XTickLabelRotation',90);
set(gca, 'YTick', 1:21, 'YTickLabel', ica2yeo7.Yeo7N);

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

nodes=dlmread('icad25_nodesordered.csv'); tmp=nanmean(ica_partial(id_healthy==0,:))';
nodes=[Tdep', Tdepanx', Tanx', (1:210)',tmp, nodes]; nodes=sortrows(nodes, 1);

%%%% distributions plotting
roi=find(Tdepanx==min(Tdepanx)); Tdep(roi)
allobservations=ica_partial(:,roi); Tdepanx(roi)
figure('Color','w'); Tanx(roi)
hold on;Tstr(roi)
subplot(4,1,1);plot_histogram_shaded(allobservations(strcmp(clin, 'ahc')),'Alpha',0.3, 'Normalization', 'pdf');
subplot(4,1,1);plot_histogram_shaded(allobservations(strcmp(clin, 'dep')),'Alpha',0.3, 'Normalization', 'pdf'); xlim([-1.5 3]);
subplot(4,1,2);plot_histogram_shaded(allobservations(strcmp(clin, 'ahc')),'Alpha',0.3, 'Normalization', 'pdf');
subplot(4,1,2);plot_histogram_shaded(allobservations(strcmp(clin, 'depanx')),'Alpha',0.3, 'Normalization', 'pdf');xlim([-1.5 3]);
subplot(4,1,3);plot_histogram_shaded(allobservations(strcmp(clin, 'ahc')),'Alpha',0.3, 'Normalization', 'pdf');
subplot(4,1,3);plot_histogram_shaded(allobservations(strcmp(clin, 'anx')),'Alpha',0.3, 'Normalization', 'pdf');xlim([-1.5 3]);
subplot(4,1,4);plot_histogram_shaded(allobservations(strcmp(clin, 'ahc')),'Alpha',0.3, 'Normalization', 'pdf');
subplot(4,1,4);plot_histogram_shaded(allobservations(strcmp(clin, 'str')),'Alpha',0.3, 'Normalization', 'pdf');xlim([-1.5 3]);
hold off

%ALTERNATIVE:%figure('Color','w');hold on; plot_histogram_shaded(allobservations(strcmp(clin, 'dep')),'Alpha',0.3, 'Normalization', 'pdf');plot_histogram_shaded(allobservations(strcmp(clin, 'depanx')),'Alpha',0.3, 'Normalization', 'pdf');plot_histogram_shaded(allobservations(strcmp(clin, 'anx')),'Alpha',0.3, 'Normalization', 'pdf');plot_histogram_shaded(allobservations(strcmp(clin, 'str')),'Alpha',0.3, 'Normalization', 'pdf');plot_histogram_shaded(allobservations(strcmp(clin, 'ahc')),'Alpha',0.3, 'Normalization', 'pdf');hold off


%% %% %% %%%% %%%% %%%% %%%% %%
%% depression vs anxiety vs stress -  structural
%% %% %% %%%% %%%% %%%% %%%% %%
clear; dataimport
ct_ordered=dlmread('D:\Canada_2020\UK_biobank\reports\ordered\ct_ordered.csv');
permutations=1000;
ix=cellfun('isempty', clinical)' | sum(isnan(ct_ordered)')'>0 |isnan(whitebritish_ordered);
ct_ordered(ix,:)=[]; a=age; a(ix)=[];s=sex;s(ix)=[];clin=clinical; clin(ix)=[];
cognitive=cognitive_ordered; cognitive(ix,:)=[]; minimal=minimal_ordered; minimal(ix,:)=[]; MDD_prs=MDD_prs_ordered; MDD_prs(ix)=[]; ANX_prs=ANX_prs_ordered; ANX_prs(ix)=[]; PTSD_prs=PTSD_prs_ordered; PTSD_prs(ix)=[];whitebritish=whitebritish_ordered;whitebritish(ix)=[];

%%%% subset controls
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

ct_ordered(ix,:)=[]; a(ix)=[];s(ix)=[]; clin(ix)=[]; minimal(ix,:)=[];
MDD_prs(ix)=[]; ANX_prs(ix)=[]; PTSD_prs(ix)=[];whitebritish(ix)=[];
1-sum(s(strcmp(clin, 'ahc')))/sum(~strcmp(clin, 'ahc'))
figure; histogram(a(strcmp(clin, 'ahc'))); hold on; histogram(a(~strcmp(clin, 'ahc')))

for i=1:length(ct_ordered(1,:))
        allobservations=ct_ordered(:,i);i
        parfor n = 1:permutations; 
        permutation_index = randperm(length(allobservations));
        randomSample = allobservations(permutation_index,:);
        T=mktblct(a,s,clin', minimal.x54_2_0,   MDD_prs,ANX_prs,PTSD_prs,   randomSample);T=sortrows(T, 3);
        mdl = fitlm(T,'y~x1+x1^2+x1*x2+x2+x3+x4 +x5+x6+x7 ');Td(1,n)=mdl.Coefficients.tStat(5); Tda(1,n)=mdl.Coefficients.tStat(6); Ts(1,n)=mdl.Coefficients.tStat(7); Ta(1,n)=mdl.Coefficients.tStat(4);
        anovamdl=anova(mdl); F(1,n)=anovamdl.F(3);
        end; 
        F_perm(i,:)=F; Tperm_dep(i,:)=Td;Tperm_depan(i,:)=Tda;Tperm_str(i,:)=Ts;Tperm_anx(i,:)=Ta;
  clear T; T=mktblct(a,s,clin', minimal.x54_2_0,     MDD_prs,ANX_prs,PTSD_prs,     allobservations);T=sortrows(T, 3);
  mdl = fitlm(T,'y~x1+x1^2+x1*x2+x2+x3+x4  +x5+x6+x7'); anovamdl=anova(mdl);
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

%%%% distributions plotting
roi=find(Tdep==min(Tdep)); Tdep(roi)
allobservations=ct_ordered(:,roi); Tdepanx(roi)
figure('Color','w'); Tanx(roi)
hold on;Tstr(roi)
subplot(4,1,1);plot_histogram_shaded(allobservations(strcmp(clin, 'ahc')),'Alpha',0.3, 'Normalization', 'pdf');
subplot(4,1,1);plot_histogram_shaded(allobservations(strcmp(clin, 'dep')),'Alpha',0.3, 'Normalization', 'pdf'); xlim([2.2 3.8]);
subplot(4,1,2);plot_histogram_shaded(allobservations(strcmp(clin, 'ahc')),'Alpha',0.3, 'Normalization', 'pdf');
subplot(4,1,2);plot_histogram_shaded(allobservations(strcmp(clin, 'depanx')),'Alpha',0.3, 'Normalization', 'pdf');xlim([2.2 3.8]);
subplot(4,1,3);plot_histogram_shaded(allobservations(strcmp(clin, 'ahc')),'Alpha',0.3, 'Normalization', 'pdf');
subplot(4,1,3);plot_histogram_shaded(allobservations(strcmp(clin, 'anx')),'Alpha',0.3, 'Normalization', 'pdf');xlim([2.2 3.8]);
subplot(4,1,4);plot_histogram_shaded(allobservations(strcmp(clin, 'ahc')),'Alpha',0.3, 'Normalization', 'pdf');
subplot(4,1,4);plot_histogram_shaded(allobservations(strcmp(clin, 'str')),'Alpha',0.3, 'Normalization', 'pdf');xlim([2.2 3.8]);
hold off

%% cognitive_ordered outcomes
cognitive_ordered.x6350_2_0(cognitive_ordered.x6350_2_0> (nanmean(cognitive_ordered.x6350_2_0)+4*nanstd(cognitive_ordered.x6350_2_0)) | cognitive_ordered.x6350_2_0<100)=NaN; %Duration to complete alphanumeric path 
cognitive_ordered.x6348_2_0(cognitive_ordered.x6348_2_0<100 | cognitive_ordered.x6348_2_0> (nanmean(cognitive_ordered.x6350_2_0)+3*nanstd(cognitive_ordered.x6348_2_0)) )=NaN; %Duration to complete numeric/easy path 
cognitive_ordered.tmt_cor=(cognitive_ordered.x6350_2_0 + 5*cognitive_ordered.x6351_2_0);% - (cognitive_ordered.x6348_2_0 +5*cognitive_ordered.x6349_2_0) ;
%ix=(cellfun('isempty', clinical)' | isnan(ica_partial(:,1))); a=age; a(ix)=[];s=sex;s(ix)=[];clin=clinical; clin(ix)=[];
%cognitive=cognitive_ordered; cognitive(ix,:)=[];minimal=minimal_ordered; minimal(ix,:)=[]; MDD_prs=MDD_prs_ordered; MDD_prs(ix)=[];
%% Cogntive data analysis:
ix=cellfun('isempty', clinical)'; a=age; a(ix)=[];s=sex;s(ix)=[];clin=clinical; clin(ix)=[];
cognitive=cognitive_ordered; cognitive(ix,:)=[];minimal=minimal_ordered; minimal(ix,:)=[]; MDD_prs=MDD_prs_ordered; MDD_prs(ix)=[]; ANX_prs=ANX_prs_ordered; ANX_prs(ix)=[]; PTSD_prs=PTSD_prs_ordered; PTSD_prs(ix)=[];

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

a(ix)=[];s(ix)=[]; clin(ix)=[]; minimal(ix,:)=[];  cognitive(ix,:)=[]; 
%%% alphanum trails
T=table(a,s,clin',minimal.x54_2_0, cognitive.tmt_cor); T(isnan(T.Var5),:)=[];T.Var5=zscore(T.Var5);T=sortrows(T, 3);
mdl = fitlm(T,'Var5~a+s+a^2+a*s+Var3+Var4')
tmp=mdl.Coefficients{4:7,:};
%%% fluid intelligence
T=table(a,s,clin', minimal.x54_2_0,cognitive.x20016_2_0); T(isnan(T.Var5),:)=[];T.Var5=zscore(T.Var5);T=sortrows(T, 3);
mdl = fitlm(T,'Var5~a+s+a^2+a*s+Var3+Var4')
tmp=horzcat(tmp, mdl.Coefficients{4:7,:});
%%% word pairs - PAL 
T=table(a,s,clin', minimal.x54_2_0,cognitive.x20197_2_0);T(isnan(T.Var5),:)=[];T.Var5=zscore(T.Var5);T=sortrows(T, 3);
mdl = fitlm(T,'Var5~a+s+a^2+a*s+Var3+Var4')
tmp=horzcat(tmp, mdl.Coefficients{4:7,:});
%%% symbol digit
T=table(a,s,clin',minimal.x54_2_0, cognitive.x23324_2_0);T(isnan(T.Var5),:)=[];T.Var5=zscore(T.Var5);T=sortrows(T, 3);
mdl = fitlm(T,'Var5~a+s+a^2+a*s+Var3+Var4')
tmp=horzcat(tmp, mdl.Coefficients{4:7,:});

%%% pairs matching (Number of incorrect matches in round) - 399.1-2 
T=table(a,s,clin',minimal.x54_2_0, cognitive.x399_2_2+cognitive.x399_2_1); T(isnan(T.Var5),:)=[];T.Var5=zscore(T.Var5);T=sortrows(T, 3);
mdl = fitlm(T,'Var5~a+s+a^2+a*s+Var3+Var4')
tmp=horzcat(tmp, mdl.Coefficients{4:7,:});
%%% reaction time 
T=table(a,s,clin',minimal.x54_2_0, clean(cognitive.x20023_2_0, 3));T(isnan(T.Var5),:)=[];T.Var5=zscore(T.Var5);T=sortrows(T, 3);
mdl = fitlm(T,'Var5~a+s+a^2+a*s+Var3+Var4')
tmp=horzcat(tmp, mdl.Coefficients{4:7,:});
%%% digit span
T=table(a,s,clin',minimal.x54_2_0, clean(cognitive.x4282_2_0,3));T(isnan(T.Var5),:)=[];T.Var5=zscore(T.Var5);T=sortrows(T, 3);
mdl = fitlm(T,'Var5~a+s+a^2+a*s+Var3+Var4')
tmp=horzcat(tmp, mdl.Coefficients{4:7,:});
%%% Raven's matrices
T=table(a,s,clin',minimal.x54_2_0, cognitive.x6373_2_0);T(isnan(T.Var5),:)=[];T.Var5=zscore(T.Var5);T=sortrows(T, 3);
mdl = fitlm(T,'Var5~a+s+a^2+a*s+Var3+Var4')
tmp=horzcat(tmp, mdl.Coefficients{4:7,:});
%%% Tower of london - nr of figures done 21004
T=table(a,s,clin', minimal.x54_2_0, cognitive.x21004_2_0);T(isnan(T.Var5),:)=[];T.Var5=zscore(T.Var5);T=sortrows(T, 3);
mdl = fitlm(T,'Var5~a+s+a^2+a*s+Var3+Var4')
tmp=horzcat(tmp, mdl.Coefficients{4:7,:});

%%% neuroticism
T=table(a,s,clin', minimal.x54_2_0,  cognitive.x20127_0_0);T(isnan(T.Var5),:)=[];T.Var5=zscore(T.Var5);T=sortrows(T, 3);
mdl = fitlm(T,'Var5~a+s+a^2+a*s+Var3+Var4')
tmp=horzcat(tmp, mdl.Coefficients{4:7,:});

FDR=mafdr(reshape(tmp(:,4:4:40),[1 40]), 'BHFDR', true);FDR=reshape(FDR, [4,10])


[RHO,PVAL]=corr([cognitive.tmt_cor, cognitive.x20016_2_0, ... %GF
cognitive.x20197_2_0,...%%% word pairs - PAL 
cognitive.x23324_2_0,...%%% symbol digit
cognitive.x399_2_2+cognitive.x399_2_1,... %%% pairs matching 
cognitive.x20023_2_0,... %%% reaction time 
cognitive.x4282_2_0,...%%% digit span
cognitive.x6373_2_0,...%%% Raven's matrices
cognitive.x21004_2_0,... %%% Tower of london 
cognitive.x20127_0_0], 'rows','complete')%%% neuroticism

%% clinical differences between samples
ix=(cellfun('isempty', clinical)'| isnan(ica_partial(:,1))) ; a=age; a(ix)=[];s=sex;s(ix)=[];clin=clinical; clin(ix)=[];
cognitive=cognitive_ordered; cognitive(ix,:)=[]; minimal=minimal_ordered; minimal(ix,:)=[];  MDD_prs=MDD_prs_ordered; MDD_prs(ix)=[]; ANX_prs=ANX_prs_ordered; ANX_prs(ix)=[]; PTSD_prs=PTSD_prs_ordered; PTSD_prs(ix)=[];
healthydate=healthydate_ordered; healthydate(ix,:)=[]; 
%%%% 2050/2060/2070/2080
PHQ2=(minimal.x2050_2_0 + minimal.x2060_2_0)>=2 ; sum(PHQ2>0)
tmp=table(PHQ2,... %%% PHQ2
    cognitive.x20127_0_0,... %neuroticism
minimal.x2070_2_0>=2,... %%% restless
minimal.x2080_2_0>=2,... %%% tired
tmpssri, tmpsnri, tmptca, tmpmaoi, tmpatypical,...
MDD_prs,ANX_prs, PTSD_prs, minimal.x21003_2_0, minimal.x31_0_0, clin', 'VariableNames',{'PHQ2','neurotic','restless', 'tired', 'ssri', 'snri', 'tca', 'maoi', 'atypical','mdd_prs','ANX_prs', 'PTSD_prs', 'age','sex','Clin'}); 
statsm=grpstats(tmp, 'Clin', {'std'}); statsm=statsm([1 2 5 3 4], :)
[tbl,chi2,p]=crosstab(minimal.x31_0_0, clin')
[tbl,chi2,p]=crosstab(PHQ2, clin')
[tbl,chi2,p]=crosstab(minimal.x2070_2_0>=2, clin')
[tbl,chi2,p]=crosstab(minimal.x2080_2_0>=2, clin')
[tbl,chi2,p]=crosstab(tmpatypical, clin') %tmpsnri, tmptca, tmpmaoi, tmpatypical
anova1(minimal.x21003_2_0,clin')
anova1(MDD_prs,clin')
anova1(ANX_prs,clin')
anova1(PTSD_prs,clin')
clear tmp statsm


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%  %%% time since depression %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
minimal.x20433_0_0(minimal.x20433_0_0<0)=NaN;%Age at first episode of depression exclude <0
minimal.x20434_0_0(minimal.x20434_0_0<0)=NaN;%Age at last episode of depression exclude <0

%firstonset=minimal.x21003_2_0-minimal.x20433_0_0;
lastonset=minimal.x21003_2_0-minimal.x20434_0_0;
%histogram(lastonset(ismember(clin','anx') ))
figure; histogram(-1*lastonset(~ismember(clin','hc') ), 150); xlim([-50 3]); hold on

healthydate.x130906_0_0(healthydate.x130906_0_0<0)=NaN;
healthydate.x130910_0_0(healthydate.x130910_0_0<0)=NaN;
firstonset=minimal.x21003_2_0-(healthydate.x130894_0_0 -healthydate.x34_0_0)
histogram(-1*firstonset(~ismember(clin','hc') ), 20, 'Normalization','probability'); xlim([-60 3]); ylim([0 0.2])
firstonset=minimal.x21003_2_0-(healthydate.x130906_0_0-healthydate.x34_0_0)
hold on; histogram(-1*firstonset(~ismember(clin','hc')), 20, 'Normalization','probability' ); xlim([-60 3]); 
firstonset=minimal.x21003_2_0-(healthydate.x130910_0_0-healthydate.x34_0_0)
hold on; histogram(-1*firstonset(~ismember(clin','hc')), 20, 'Normalization','probability' ); xlim([-60 3]); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%  medication  %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('D:\Canada_2020\UK_biobank\reports\ordered\medication_ordered.mat')
medication_dep_clin=medication_ordered; medication_dep_clin(ix,:)=[]; 
clear medication_class
med_map=readtable('med_map.csv'); medication_class=zeros(length(medication_dep_clin), 10);medication_class(medication_class==0)=NaN;
for i=1:length(medication_dep_clin); clear tmp
    tmp=medication_dep_clin{i};
    if ~istable(tmp)
        medication_class(i)=NaN;
    else 
        tmp=unique(tmp.generic_name);
        for med=1:length(tmp)
            try; meds(med)=unique(med_map.Var4( ~cellfun('isempty', strfind(med_map.Var2, tmp{med})) ));end
        end
        medication_class(i,1:length(unique(meds) ))=sort(unique(meds));clear meds
    end
end
tmpssri=sum((medication_class==1)')'; %%%% SSRI 1; MAOI 3; TCA 2; SNRI 4; Atypical 5  
tmptca=sum((medication_class==2)')'; tmpmaoi=sum((medication_class==3)')'; tmpsnri=sum((medication_class==4)')'; tmpatypical=sum((medication_class==5)')'; 


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



square_P_full = zeros(D); square_P_full(triu(ones(D),1)>0) = p_perm_inter; 
square_T_full = zeros(D); square_T_full(triu(ones(D),1)>0) = Tinteract; square_T_full(square_P_full>0.05)=0;
figure(1);imagesc(square_T_full); set(gca,'XTick',[1:21], 'YTick',[1:21]); %clear square_* 

square_P_full = zeros(D); square_P_full(triu(ones(D),1)>0) = p_perm_edu; 
square_T_full = zeros(D); square_T_full(triu(ones(D),1)>0) = Tedu; square_T_full(square_P_full>0.05)=0;
figure(2);imagesc(square_T_full); set(gca,'XTick',[1:21], 'YTick',[1:21]); %clear square_* 

