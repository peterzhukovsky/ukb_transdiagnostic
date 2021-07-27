%% tabular data
cd D:\Canada_2020\UK_biobank\reports\ordered
selfreport_ordered=readtable('selfreport_ordered.csv');
healthy_ordered=readtable('healthy_ordered.csv');healthydate_ordered=readtable('healthydate_ordered.csv');
minimal_ordered=readtable('minimal_ordered.csv');
ct_ordered=dlmread('ct_ordered.csv');
ms_ordered=dlmread('ms_ordered.csv');
gc_ordered_dem=dlmread('gc_ordered_dem.csv');
gc_ordered=dlmread('gc_ordered.csv');
yeormat_ordered=dlmread('yeormat_ordered.csv');
cognitive_ordered=readtable('cognitive_ordered.csv');
thicknessFS=readtable('thicknessFS_ordered.csv');thicknessFS=thicknessFS{:,:}; thicknessFS(:,1:2)=[];
physical_ordered=readtable('physical_ordered.csv');
MDD_prs_ordered=dlmread('MDD_prs_ordered.csv');y=MDD_prs_ordered;X=[healthydate_ordered.x22009_0_1,healthydate_ordered.x22009_0_2, healthydate_ordered.x22009_0_3, healthydate_ordered.x22009_0_4, healthydate_ordered.x22009_0_5,healthydate_ordered.x22009_0_6,healthydate_ordered.x22009_0_7,healthydate_ordered.x22009_0_8,healthydate_ordered.x22009_0_9,healthydate_ordered.x22009_0_10]; [b,bint,r] = regress(y,X); MDD_prs_ordered=r;
ANX_prs_ordered=dlmread('ANX_prs_ordered.csv');y=ANX_prs_ordered;X=[healthydate_ordered.x22009_0_1,healthydate_ordered.x22009_0_2, healthydate_ordered.x22009_0_3, healthydate_ordered.x22009_0_4, healthydate_ordered.x22009_0_5,healthydate_ordered.x22009_0_6,healthydate_ordered.x22009_0_7,healthydate_ordered.x22009_0_8,healthydate_ordered.x22009_0_9,healthydate_ordered.x22009_0_10]; [b,bint,r] = regress(y,X); ANX_prs_ordered=r;
PTSD_prs_ordered=dlmread('PTSD_prs_ordered.csv');y=PTSD_prs_ordered;X=[healthydate_ordered.x22009_0_1,healthydate_ordered.x22009_0_2, healthydate_ordered.x22009_0_3, healthydate_ordered.x22009_0_4, healthydate_ordered.x22009_0_5,healthydate_ordered.x22009_0_6,healthydate_ordered.x22009_0_7,healthydate_ordered.x22009_0_8,healthydate_ordered.x22009_0_9,healthydate_ordered.x22009_0_10]; [b,bint,r] = regress(y,X); PTSD_prs_ordered=r;

cd D:\Canada_2020\UK_biobank\data
%ica_partial=dlmread('ica_d25_full.csv');
ica_partial=dlmread('ica_d25_par.csv');

cd D:\Canada_2020\UK_biobank\reports
datanames_rs=readtable('subjects_rs.txt');
datanames_fs=readtable('subjects_fs.txt');
datanames_all=unique(vertcat(datanames_fs.Var1, datanames_rs.Var1));
regions=readtable('thicknessFS_names.xlsx');regions=regions.names_short;regions(37)=[];regions(1)=[];%regions(72)=[];regions(69)=[];regions(37)=[];regions(36)=[];regions(33)=[]; regions(1)=[];

minimal_ordered.x2050_2_0=minimal_ordered.x2050_2_0-1; minimal_ordered.x2060_2_0=minimal_ordered.x2060_2_0-1;
minimal_ordered.x2050_2_0(minimal_ordered.x2050_2_0<0)=NaN;minimal_ordered.x2060_2_0(minimal_ordered.x2060_2_0<0)=NaN;
minimal_ordered.x2070_2_0(minimal_ordered.x2070_2_0<0)=NaN; minimal_ordered.x2080_2_0(minimal_ordered.x2080_2_0<0)=NaN; 
llabel_names=readtable('D:\Canada_2020\UK_biobank\data\ROI_names\lh.rsfc_HCP.txt');llabel_names=llabel_names.Var5;
rlabel_names=readtable('D:\Canada_2020\UK_biobank\data\ROI_names\rh.rsfc_HCP.txt');rlabel_names=rlabel_names.Var5;
label_names_all=vertcat(llabel_names,rlabel_names);  label_names_HOA={'Occipital Cortex';'Occipital Cortex';'Occipital Cortex';'Occipital Cortex';'Occipital Cortex';'Occipital Cortex';'Occipital Cortex';'Postcentral Gyrus';'Precentral/postcentral Gyrus';'Precentral Gyrus';'Precentral Gyrus';'Middle Frontal Gyrus';'Occipital Cortex';'Posterior Cingulate Cortex';'Precuneous ';'Occipital Cortex';'Lateral Occipital Cortex ';'Temporal/Occipital Fisiform Cortex';'Lateral Occipital Cortex ';'Lateral Occipital Cortex ';'Lateral Occipital Cortex ';'Lateral Occipital Cortex ';'Lateral Occipital Cortex ';'Heschl''s Gyrus';'Supramarginal Gyrus (posterior)/Angular Gyrus';'Superior Frontal Gyrus';'Precuneous ';'Supramarginal Gyrus (posterior)/Angular Gyrus';'Precuneous ';'Precuneous ';'Precuneous ';'Posterior Cingulate Cortex';'Posterior Cingulate Cortex';'Posterior Cingulate Cortex';'Posterior Cingulate Cortex';'Medial Postcentral Gyrus ';'Medial Postcentral Gyrus ';'Posterior Cingulate Cortex';'Superior Parietal Lobule';'Medial Postcentral Gyrus ';'Medial Postcentral Gyrus ';'Superior Parietal Lobule';'Supplementary Motor Cortex';'Superior Frontal Gyrus';'Precuneous Cortex ';'Precuneous Cortex ';'Superior Parietal Lobule ';'Superior Parietal Lobule ';'Lateral Occipital Cortex ';'Lateral Occipital Cortex ';'Postcentral Gyrus';'Postcentral Gyrus';'Middle Postcentral Gyrus';'Precentral Gyrus';'Precentral Gyrus';'Precentral Gyrus';'Anterior Cingulate Cortex';'Anterior Cingulate Cortex';'Anterior Cingulate Cortex (anterior)';'Anterior Cingulate Cortex (anterior)';'Anterior Cingulate Cortex (ventral, anterior)';'Paracingulate Gyrus';'Medial Superior Frontal Gyrus';'Medial Superior Frontal Gyrus';'Paracingulate Gyrus (anterior)';'Frontal Orbital Cortex';'Middle Frontal Gyrus';'Superior/Middle Frontal Gyrus';'Superior Frontal Gyrus/Frontal Pole';'Frontal Pole';'Frontal Pole (anterior)';'Frontal Pole (anterior)';'Middle Frontal Gyrus/ Inferior Frontal Gyrus';'Inferior Frontal Gyrus, pars opercularis ';'Inferior Frontal Gyrus, pars opercularis ';'Frontal Orbital Cortex';'Frontal Orbital Cortex/Frontal Pole (lateral)';'Inferior Frontal Gyrus, pars opercularis ';'Inferior Frontal Gyrus, pars opercularis ';'Precentral Gyrus/Inferior Frontal Gyrus';'Inferior Frontal Gyrus/Middle Frontal Gyrus';'Inferior Frontal Gyrus, pars triangularis ';'Middle Frontal Gyrus';'Middle Frontal Gyrus';'Middle Frontal Gyrus';'Frontal Pole';'Frontal Pole (dorsal anterior)';'Frontal Medial Cortex ';'Frontal Pole (ventral)';'Frontal Pole (ventral)';'Frontal Orbital Cortex/Frontal Pole (medial)';'Frontal Orbital Cortex (Medial)';'Frontal Orbital Cortex (Medial)';'Frontal Orbital Cortex (Medial)';'Superior Parietal Lobule';'Precentral Gyrus (medial)';'Middle Frontal Gyrus ';'Middle Frontal Gyrus ';'Precentral Gyrus (lateral ventral)';'Central Opercular Cortex';'Parietal Opercular Cortex';'Insular/Parietal Opercular Cortex';'Heschl''s Gyrus';'Parietal Opercular Cortex';'Parietal Opercular Cortex';'Insular Cortex (posterior)';'Insular Cortex (ventral)';'Frontal Operculum Cortex ';'Insular Cortex (dorsal)';'Insular Cortex (ventral)';'Insular Cortex (anterior)/Frontal Orbital Cortex';'Insular Cortex (anterior)';'Central Opercular Cortex ';'Central Opercular Cortex ';'Insular Cortex (dorsal)';'Supramarginal Gyrus (anterior)';'Superior Parietal Lobule';'Parahippocampal Gyrus';'Parahippocampal Gyrus';'Parahippocampal Gyrus';'Posterior Cingulate Cortex (ventral)';'Parahippocampal Gyrus';'Temporal Pole (aTL)';'Planum Temporale (STG)';'Superior Temporal Gyrus (anterior)';'Parahippocampal/Lingual Gyrus';'Temporal Fusiform Cortex';'Middle Temporal Gyrus (anterior)';'Middle Temporal Gyrus (posterior)';'Middle Temporal Gyrus (posterior)';'Temporal Pole';'Middle Temporal Gyrus (anterior)';'Inferior Temporal Gyrus, temporooccipital ';'Inferior Temporal Gyrus, posterior division ';'Temporal Fusiform Cortex, posterior division ';'Inferior Temporal Gyrus, temporooccipital part ';'Middle Temporal Gyrus, temporooccipital part ';'Inferior Temporal Gyrus, temporooccipital part ';'Middle Temporal Gyrus, temporooccipital part ';'Lateral Occipital Cortex, inferior division ';'Lateral Occipital Cortex, superior division ';'Precuneous';'Lateral Occipital Cortex, superior division ';'Supramarginal Gyrus, posterior division ';'Lateral Occipital Cortex, superior division ';'Lateral Occipital Cortex, superior division ';'Postcentral Gyrus (ventral)';'Supramarginal Gyrus, anterior division ';'Angular Gyrus/Supramarginal Gyrus';'Angular Gyrus/Lateral Occipital Cortex';'Angular Gyrus/Lateral Occipital Cortex';'Lateral Occipital Cortex';'Lingual Gyrus';'Temporal Occipital Fusiform Cortex';'Temporal Fusiform Cortex, posterior division ';'Lateral Occipital Cortex, inferior division ';'Lateral Occipital Cortex, inferior division ';'Lateral Occipital Cortex, inferior division ';'Lateral Occipital Cortex, inferior division ';'Lingual Gyrus';'Precuneous Cortex ';'Precuneous Cortex ';'Temporal Occipital Fusiform Cortex ';'Subcallosal Cortex (Subgenual ACC)';'Paracingulate Gyrus (anterior, ventral)';'Subcallosal Cortex (Subgenual ACC)';'Insular Cortex';'Insular Cortex/Parietal Operculum';'Frontal Operculum Cortex ';'Frontal Pole (medial)';'Frontal Pole (lateral)';'Inferior Temporal Gyrus, anterior division ';'Heschl''s Gyrus ';'Planum Temporale (STG)';'Planum Temporale (STG)';'Middle Temporal Gyrus, posterior division';'Middle Temporal Gyrus, posterior division ';'Planum Polare';'Anterior Cingulate';'Anterior Cingulate';'Occipital Cortex';'Occipital Cortex';'Occipital Cortex';'Occipital Cortex';'Occipital Cortex';'Occipital Cortex';'Occipital Cortex';'Postcentral Gyrus';'Precentral/postcentral Gyrus';'Precentral Gyrus';'Precentral Gyrus';'Middle Frontal Gyrus';'Occipital Cortex';'Posterior Cingulate Cortex';'Precuneous ';'Occipital Cortex';'Lateral Occipital Cortex ';'Temporal/Occipital Fisiform Cortex';'Lateral Occipital Cortex ';'Lateral Occipital Cortex ';'Lateral Occipital Cortex ';'Lateral Occipital Cortex ';'Lateral Occipital Cortex ';'Heschl''s Gyrus';'Supramarginal Gyrus (posterior)/Angular Gyrus';'Superior Frontal Gyrus';'Precuneous ';'Supramarginal Gyrus (posterior)/Angular Gyrus';'Precuneous ';'Precuneous ';'Precuneous ';'Posterior Cingulate Cortex';'Posterior Cingulate Cortex';'Posterior Cingulate Cortex';'Posterior Cingulate Cortex';'Medial Postcentral Gyrus ';'Medial Postcentral Gyrus ';'Posterior Cingulate Cortex';'Superior Parietal Lobule';'Medial Postcentral Gyrus ';'Medial Postcentral Gyrus ';'Superior Parietal Lobule';'Supplementary Motor Cortex';'Superior Frontal Gyrus';'Precuneous Cortex ';'Precuneous Cortex ';'Superior Parietal Lobule ';'Superior Parietal Lobule ';'Lateral Occipital Cortex ';'Lateral Occipital Cortex ';'Postcentral Gyrus';'Postcentral Gyrus';'Middle Postcentral Gyrus';'Precentral Gyrus';'Precentral Gyrus';'Precentral Gyrus';'Anterior Cingulate Cortex';'Anterior Cingulate Cortex';'Anterior Cingulate Cortex (anterior)';'Anterior Cingulate Cortex (anterior)';'Anterior Cingulate Cortex (ventral, anterior)';'Paracingulate Gyrus';'Medial Superior Frontal Gyrus';'Medial Superior Frontal Gyrus';'Paracingulate Gyrus (anterior)';'Frontal Orbital Cortex';'Middle Frontal Gyrus';'Superior/Middle Frontal Gyrus';'Superior Frontal Gyrus/Frontal Pole';'Frontal Pole';'Frontal Pole (anterior)';'Frontal Pole (anterior)';'Middle Frontal Gyrus/ Inferior Frontal Gyrus';'Inferior Frontal Gyrus, pars opercularis ';'Inferior Frontal Gyrus, pars opercularis ';'Frontal Orbital Cortex';'Frontal Orbital Cortex/Frontal Pole (lateral)';'Inferior Frontal Gyrus, pars opercularis ';'Inferior Frontal Gyrus, pars opercularis ';'Precentral Gyrus/Inferior Frontal Gyrus';'Inferior Frontal Gyrus/Middle Frontal Gyrus';'Inferior Frontal Gyrus, pars triangularis ';'Middle Frontal Gyrus';'Middle Frontal Gyrus';'Middle Frontal Gyrus';'Frontal Pole';'Frontal Pole (dorsal anterior)';'Frontal Medial Cortex ';'Frontal Pole (ventral)';'Frontal Pole (ventral)';'Frontal Orbital Cortex/Frontal Pole (medial)';'Frontal Orbital Cortex (Medial)';'Frontal Orbital Cortex (Medial)';'Frontal Orbital Cortex (Medial)';'Superior Parietal Lobule';'Precentral Gyrus (medial)';'Middle Frontal Gyrus ';'Middle Frontal Gyrus ';'Precentral Gyrus (lateral ventral)';'Central Opercular Cortex';'Parietal Opercular Cortex';'Insular/Parietal Opercular Cortex';'Heschl''s Gyrus';'Parietal Opercular Cortex';'Parietal Opercular Cortex';'Insular Cortex (posterior)';'Insular Cortex (ventral)';'Frontal Operculum Cortex ';'Insular Cortex (dorsal)';'Insular Cortex (ventral)';'Insular Cortex (anterior)/Frontal Orbital Cortex';'Insular Cortex (anterior)';'Central Opercular Cortex ';'Central Opercular Cortex ';'Insular Cortex (dorsal)';'Supramarginal Gyrus (anterior)';'Superior Parietal Lobule';'Parahippocampal Gyrus';'Parahippocampal Gyrus';'Parahippocampal Gyrus';'Posterior Cingulate Cortex (ventral)';'Parahippocampal Gyrus';'Temporal Pole (aTL)';'Planum Temporale (STG)';'Superior Temporal Gyrus (anterior)';'Parahippocampal/Lingual Gyrus';'Temporal Fusiform Cortex';'Middle Temporal Gyrus (anterior)';'Middle Temporal Gyrus (posterior)';'Middle Temporal Gyrus (posterior)';'Temporal Pole';'Middle Temporal Gyrus (anterior)';'Inferior Temporal Gyrus, temporooccipital ';'Inferior Temporal Gyrus, posterior division ';'Temporal Fusiform Cortex, posterior division ';'Inferior Temporal Gyrus, temporooccipital part ';'Middle Temporal Gyrus, temporooccipital part ';'Inferior Temporal Gyrus, temporooccipital part ';'Middle Temporal Gyrus, temporooccipital part ';'Lateral Occipital Cortex, inferior division ';'Lateral Occipital Cortex, superior division ';'Precuneous';'Lateral Occipital Cortex, superior division ';'Supramarginal Gyrus, posterior division ';'Lateral Occipital Cortex, superior division ';'Lateral Occipital Cortex, superior division ';'Postcentral Gyrus (ventral)';'Supramarginal Gyrus, anterior division ';'Angular Gyrus/Supramarginal Gyrus';'Angular Gyrus/Lateral Occipital Cortex';'Angular Gyrus/Lateral Occipital Cortex';'Lateral Occipital Cortex';'Lingual Gyrus';'Temporal Occipital Fusiform Cortex';'Temporal Fusiform Cortex, posterior division ';'Lateral Occipital Cortex, inferior division ';'Lateral Occipital Cortex, inferior division ';'Lateral Occipital Cortex, inferior division ';'Lateral Occipital Cortex, inferior division ';'Lingual Gyrus';'Precuneous Cortex ';'Precuneous Cortex ';'Temporal Occipital Fusiform Cortex ';'Subcallosal Cortex (Subgenual ACC)';'Paracingulate Gyrus (anterior, ventral)';'Subcallosal Cortex (Subgenual ACC)';'Insular Cortex';'Insular Cortex/Parietal Operculum';'Frontal Operculum Cortex ';'Frontal Pole (medial)';'Frontal Pole (lateral)';'Inferior Temporal Gyrus, anterior division ';'Heschl''s Gyrus ';'Planum Temporale (STG)';'Planum Temporale (STG)';'Middle Temporal Gyrus, posterior division';'Middle Temporal Gyrus, posterior division ';'Planum Polare';'Anterior Cingulate';'Anterior Cingulate'};
ica2yeo7=readtable('ica2yeo7.csv');

minimal_ordered.x21000_0_0; tmp=repmat( {'Missing'} , length(datanames_all),1 );
tmp(minimal_ordered.x21000_0_0==1|minimal_ordered.x21000_0_0==1001|minimal_ordered.x21000_0_0==1002|minimal_ordered.x21000_0_0==1003)={'White'};
tmp(~(minimal_ordered.x21000_0_0==1|minimal_ordered.x21000_0_0==1001|minimal_ordered.x21000_0_0==1002|minimal_ordered.x21000_0_0==1003))={'Non-White'};
minimal_ordered.x21000_0_0=tmp;
tmp=repmat( {'Missing'} , length(datanames_all),1 ); 
tmp(minimal_ordered.x54_2_0==11025)={'Cheadle'}; tmp(minimal_ordered.x54_2_0==11026)={'Reading'}; tmp(minimal_ordered.x54_2_0==11027)={'Newcastle'}; 
minimal_ordered.x54_2_0=tmp; clear tmp;

%%% alternatively use random data: load('./random_data/randomised_data.mat') 

%% %% %% %%%% %%%% %%%% %%%% %%
%% depression vs anxiety vs stress -  rsfmri ICA; setting up the clinical variables
%% %% %% %%%% %%%% %%%% %%%% %%
id_depressed_clin=( (minimal_ordered.x130895_0_0>=30 & minimal_ordered.x130895_0_0<=51) | (healthy_ordered.x130897_0_0 >=30 & healthy_ordered.x130897_0_0<=51) ); sum(id_depressed_clin>0)
id_anxious_clin=(healthy_ordered.x130907_0_0>=30 & healthy_ordered.x130907_0_0<=51); sum(id_anxious_clin>0)
id_stress_clin=(healthy_ordered.x130911_0_0>=30 & healthy_ordered.x130911_0_0<=51); sum(id_stress_clin)
id_stress=(id_stress_clin>0 & id_anxious_clin==0 & id_depressed_clin==0); 
id_anxious=(id_stress_clin==0 & id_anxious_clin>0 & id_depressed_clin==0); 
id_dep=(id_stress_clin==0 & id_anxious_clin==0 & id_depressed_clin>0); 
id_depanx=(id_stress_clin==0 & id_anxious_clin>0 & id_depressed_clin>0); 
id_healthy=sum(~isnan(healthy_ordered{:,3:145})')';
%for i=1:length(id_healthy);
%    if minimal_ordered.x21003_2_0(i)>60 & id_healthy(i)==0
%        id_healthy(i)=round(rand);
%end; end
sum(id_healthy==0)

clinical={'a'}; %clinical(age<60)={[]};
clinical(id_healthy==0)={'hc'};
clinical(id_dep==1)={'dep'};
clinical(id_depanx==1)={'depanx'};
clinical(id_anxious==1)={'anx'};
clinical(id_stress==1)={'str'};clinical(40699)=clinical(40697);
age=minimal_ordered.x21003_2_0;sex=minimal_ordered.x31_0_0; 


%%% preparing data (regressing out confounds), running PLS and using permutation testing to test for significance

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

