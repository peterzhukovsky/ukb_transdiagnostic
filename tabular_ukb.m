%% tabular data
cd D:\Canada_2020\UK_biobank\reports\ordered
selfreport_ordered=readtable('selfreport_ordered.csv');
healthy_ordered=readtable('healthy_ordered.csv');healthydate_ordered=readtable('healthydate_ordered.csv');
minimal_ordered=readtable('minimal_ordered.csv');
ct_ordered=dlmread('ct_ordered.csv');
cognitive_ordered=readtable('cognitive_ordered.csv');
thicknessFS=readtable('thicknessFS_ordered.csv');thicknessFS=thicknessFS{:,:}; thicknessFS(:,1:2)=[];
physical_ordered=readtable('physical_ordered.csv');
MDD_prs_ordered=dlmread('MDD_prs_ordered.csv');y=MDD_prs_ordered;X=[healthydate_ordered.x22009_0_1,healthydate_ordered.x22009_0_2, healthydate_ordered.x22009_0_3, healthydate_ordered.x22009_0_4, healthydate_ordered.x22009_0_5,healthydate_ordered.x22009_0_6,healthydate_ordered.x22009_0_7,healthydate_ordered.x22009_0_8,healthydate_ordered.x22009_0_9,healthydate_ordered.x22009_0_10]; [b,bint,r] = regress(y,X); MDD_prs_ordered=r;
ANX_prs_ordered=dlmread('ANX_prs_ordered.csv');y=ANX_prs_ordered;X=[healthydate_ordered.x22009_0_1,healthydate_ordered.x22009_0_2, healthydate_ordered.x22009_0_3, healthydate_ordered.x22009_0_4, healthydate_ordered.x22009_0_5,healthydate_ordered.x22009_0_6,healthydate_ordered.x22009_0_7,healthydate_ordered.x22009_0_8,healthydate_ordered.x22009_0_9,healthydate_ordered.x22009_0_10]; [b,bint,r] = regress(y,X); ANX_prs_ordered=r;
PTSD_prs_ordered=dlmread('PTSD_prs_ordered.csv');y=PTSD_prs_ordered;X=[healthydate_ordered.x22009_0_1,healthydate_ordered.x22009_0_2, healthydate_ordered.x22009_0_3, healthydate_ordered.x22009_0_4, healthydate_ordered.x22009_0_5,healthydate_ordered.x22009_0_6,healthydate_ordered.x22009_0_7,healthydate_ordered.x22009_0_8,healthydate_ordered.x22009_0_9,healthydate_ordered.x22009_0_10]; [b,bint,r] = regress(y,X); PTSD_prs_ordered=r;

cd D:\Canada_2020\UK_biobank\data
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

                %% permutation testing GLMs
permutations=1000;
ix=(cellfun('isempty', clinical)'| isnan(ica_partial(:,1))) ; a=age; a(ix)=[];s=sex;s(ix)=[];clin=clinical; clin(ix)=[];
cognitive=cognitive_ordered; cognitive(ix,:)=[]; minimal=minimal_ordered; minimal(ix,:)=[]; MDD_prs=MDD_prs_ordered; MDD_prs(ix)=[]; ANX_prs=ANX_prs_ordered; ANX_prs(ix)=[]; PTSD_prs=PTSD_prs_ordered; PTSD_prs(ix)=[];
%%%%% the models can be either ran with PRS as covariates or without - the models without PRS as covariates are shown here, but the the model with PRS as covariates can be ran by unchecking the "%"
for i=1:length(ica_partial(1,:))
        allobservations=ica_partial(:,i);allobservations(ix)=[];i
        parfor n = 1:permutations; 
        permutation_index = randperm(length(allobservations));
        randomSample = allobservations(permutation_index,:);
        T=mktbl(a,s,clin',    minimal.x25741_2_0, minimal.x54_2_0, randomSample); %x25741=motion, x21000=race (1 2 3 4), x54=site, MDD_prs, ANX_prs, PTSD_prs,
        mdl = fitlm(T,'y~x1+x1^2+x1*x2+x2+x3+x4+x5  ');Td(1,n)=mdl.Coefficients.tStat(4); Tda(1,n)=mdl.Coefficients.tStat(5); Ts(1,n)=mdl.Coefficients.tStat(6); Ta(1,n)=mdl.Coefficients.tStat(7);Tmddp(1,n)=mdl.Coefficients.tStat(8);
        anovamdl=anova(mdl); F(1,n)=anovamdl.F(3);%F2(1,n)=anovamdl.F(4);F3(1,n)=anovamdl.F(5);F4(1,n)=anovamdl.F(6);
        %T=mktbl(a,s,clin', MDD_prs, ANX_prs, PTSD_prs,    minimal.x25741_2_0, minimal.x54_2_0, randomSample);
        %mdl = fitlm(T,'y~x1+x1^2+x1*x2+x2+x3+x4+x5  +x6+x7+x8  ');Td(1,n)=mdl.Coefficients.tStat(4); Tda(1,n)=mdl.Coefficients.tStat(5); Ts(1,n)=mdl.Coefficients.tStat(6); Ta(1,n)=mdl.Coefficients.tStat(7);Tmddp(1,n)=mdl.Coefficients.tStat(8);
        end; 
        F_perm(i,:)=F; Tperm_dep(i,:)=Td;Tperm_depan(i,:)=Tda;Tperm_str(i,:)=Ts;Tperm_anx(i,:)=Ta;
        %F_mddprs(i,:)=F2;F_anxprs(i,:)=F3;F_ptsdprs(i,:)=F4; 
        clear Td Tda Ts Ta Tmddp F2 F3 F4
  clear T; T=mktbl(a,s,clin',  minimal.x25741_2_0,minimal.x54_2_0, allobservations); %x25741=motion, x21000=race (1 2 3 4), x54=site
  mdl = fitlm(T,'y~x1+x1^2+x1*x2+x2+x3+x4+x5  '); anovamdl=anova(mdl);
  %T=mktbl(a,s,clin', MDD_prs, ANX_prs, PTSD_prs,    minimal.x25741_2_0, minimal.x54_2_0, allobservations);
  %mdl = fitlm(T,'y~x1+x1^2+x1*x2+x2+x3+x4+x5  +x6+x7+x8  ');
  p_perm_anova(i)=1-sum(anovamdl.F(3)>F_perm(i,:))/permutations; 
  %p_perm_prs(i,1)=1-sum(anovamdl.F(4)> F_mddprs(i,:))/permutations; p_perm_prs(i,2)=1-sum(anovamdl.F(5)> F_anxprs(i,:))/permutations; p_perm_prs(i,3)=1-sum(anovamdl.F(6)> F_ptsdprs(i,:))/permutations; 
  %Tprs(1,i)=mdl.Coefficients.tStat(8);  Tprs(2,i)=mdl.Coefficients.tStat(9);  Tprs(3,i)=mdl.Coefficients.tStat(10);  Pprs(1,i)=mdl.Coefficients.pValue(8);  Pprs(2,i)=mdl.Coefficients.pValue(9);  Pprs(3,i)=mdl.Coefficients.pValue(10);  
  Tage(1,i)=mdl.Coefficients.tStat(2); 
  Tsex(1,i)=mdl.Coefficients.tStat(3);
  Tdep(1,i)=mdl.Coefficients.tStat(4); Pdep(1,i)=mdl.Coefficients.pValue(4);
  Tdepanx(1,i)=mdl.Coefficients.tStat(5); Pdepanx(1,i)=mdl.Coefficients.pValue(5);
  Tstr(1,i)=mdl.Coefficients.tStat(6); Pstr(1,i)=mdl.Coefficients.pValue(6);
  Tanx(1,i)=mdl.Coefficients.tStat(7); Panx(1,i)=mdl.Coefficients.pValue(7);
end
%%% sign testing correlations - potential covariates: 1) simple 2) with three prs scores and 3) Years of education, income, urban vs rural environment, marital status, stressful life events and substance use family history variables for psychiatric disease
r=corr(Tperm_dep, Tperm_anx); r=r(eye(permutations)==1); P(1,1)=1-sum(corr(Tdep', Tanx')>r)/permutations; % column1 for dep vs anx
r=corr(Tperm_dep, Tperm_str); r=r(eye(permutations)==1); P(1,3)=1-sum(corr(Tdep', Tstr')>r)/permutations; % column1 for dep vs anx
r=corr(Tperm_str, Tperm_anx); r=r(eye(permutations)==1); P(2,3)=1-sum(corr(Tstr', Tanx')>r)/permutations; % column1 for dep vs anx
r=corr(Tperm_depan, Tperm_dep); r=r(eye(permutations)==1); P(1,2)=1-sum(corr(Tdepanx', Tdep')>r)/permutations; % column1 for dep vs anx
r=corr(Tperm_depan, Tperm_anx); r=r(eye(permutations)==1); P(2,2)=1-sum(corr(Tdepanx', Tanx')>r)/permutations; % column1 for dep vs anx
r=corr(Tperm_depan, Tperm_str); r=r(eye(permutations)==1); P(3,3)=1-sum(corr(Tdepanx', Tstr')>r)/permutations; % column1 for dep vs anx
r=corr(Tperm_mddprs, Tperm_dep); r=r(eye(permutations)==1); P(1,4)=1-sum(corr(Tmddprs', Tdep')>r)/permutations; % column1 for dep vs anx
r=corr(Tperm_mddprs, Tperm_depan); r=r(eye(permutations)==1); P(3,4)=1-sum(corr(Tmddprs', Tdepanx')>r)/permutations; % column1 for dep vs anx
r=corr(Tperm_mddprs, Tperm_anx); r=r(eye(permutations)==1); P(2,4)=1-sum(corr(Tmddprs', Tanx')>r)/permutations; % column1 for dep vs anx
r=corr(Tperm_mddprs, Tperm_str); r=r(eye(permutations)==1); P(4,4)=1-sum(corr(Tmddprs', Tstr')>r)/permutations; P% column1 for dep vs anx

corr([Tdep', Tanx',Tdepanx',Tstr']) %, -1*Tmci'
                        %% visuals
index=ismember(clin,'hc'); figure(1);histogram(a(index)); hold on;yyaxis right ; %histogram(a(~index)); clear index
Tdep_square = zeros(21); Tdep_square(triu(ones(21),1)>0) = Tanx; %Tanx_square(Tanx_square==0)=NaN;
Panova_square=zeros(21);Panova_square(triu(ones(21),1)>0) = p_perm_anova; %Panx_square(Panx_square==0)=NaN;
Pdep_square=zeros(21);Pdep_square(triu(ones(21),1)>0) = Panx; Pdep_square(Pdep_square==0)=NaN;
tmp=Tdep_square; tmp(Panova_square>0.05 | Pdep_square>0.0125)=0;
figure(2); imagesc(tmp); colorbar; caxis([-5 5]); set(gca, 'XTick', 1:21, 'XTickLabel', ica2yeo7.Yeo7N, 'XTickLabelRotation',90);
set(gca, 'YTick', 1:21, 'YTickLabel', ica2yeo7.Yeo7N);

Tdep_square = zeros(21); Tdep_square(triu(ones(21),1)>0) = Tprs(2,:); %Tanx_square(Tanx_square==0)=NaN;
Panova_square=zeros(21);Panova_square(triu(ones(21),1)>0) = p_perm_prs(:,2); %= p_perm_mddprs; 
tmp=Tdep_square; tmp(Panova_square>0.05)=0;
figure(2); imagesc(tmp); colorbar; caxis([-5 5]); set(gca, 'XTick', 1:21, 'XTickLabel', ica2yeo7.Yeo7N, 'XTickLabelRotation',90);
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



%% %% %% %%%% %%%% %%%% %%%% %%
%% depression vs anxiety vs stress -  structural
%% %% %% %%%% %%%% %%%% %%%% %%
clear T* F*
permutations=1000;
ix=(cellfun('isempty', clinical)' | sum(isnan(ct_ordered)')'>0 ) ; a=age; a(ix)=[];s=sex;s(ix)=[];clin=clinical; clin(ix)=[];
cognitive=cognitive_ordered; cognitive(ix,:)=[]; minimal=minimal_ordered; minimal(ix,:)=[];  MDD_prs=MDD_prs_ordered; MDD_prs(ix)=[]; ANX_prs=ANX_prs_ordered; ANX_prs(ix)=[]; PTSD_prs=PTSD_prs_ordered; PTSD_prs(ix)=[];

for i=1:length(ct_ordered(1,:))
        allobservations=ct_ordered(:,i);allobservations(ix)=[];i
        parfor n = 1:permutations; 
        permutation_index = randperm(length(allobservations));
        randomSample = allobservations(permutation_index,:);
        T=mktblprs(a,s,clin', MDD_prs, ANX_prs, PTSD_prs,    minimal.x54_2_0, randomSample); %x25741=motion, x21000=race (1 2 3 4), x54=site, MDD_prs, ANX_prs, PTSD_prs
        mdl = fitlm(T,'y~x1+x1^2+x1*x2+x2+x3+x4 +x5+x6+x7 ');Td(1,n)=mdl.Coefficients.tStat(4); Tda(1,n)=mdl.Coefficients.tStat(5); Ts(1,n)=mdl.Coefficients.tStat(6); Ta(1,n)=mdl.Coefficients.tStat(7);Tmddp(1,n)=mdl.Coefficients.tStat(8);
        %T=mktblprs(a,s,clin',   minimal.x54_2_0, randomSample); 
        %mdl = fitlm(T,'y~x1+x1^2+x1*x2+x2+x3+x4 ');Td(1,n)=mdl.Coefficients.tStat(4); Tda(1,n)=mdl.Coefficients.tStat(5); Ts(1,n)=mdl.Coefficients.tStat(6); Ta(1,n)=mdl.Coefficients.tStat(7);
        anovamdl=anova(mdl); F(1,n)=anovamdl.F(3);
        F2(1,n)=anovamdl.F(4);F3(1,n)=anovamdl.F(5);F4(1,n)=anovamdl.F(6);
        end; 
        F_perm(i,:)=F; Tperm_dep(i,:)=Td;Tperm_depan(i,:)=Tda;Tperm_str(i,:)=Ts;Tperm_anx(i,:)=Ta;Tperm_mddprs(i,:)=Tmddp;
        F_mddprs(i,:)=F2;F_anxprs(i,:)=F3;F_ptsdprs(i,:)=F4; clear Td Tda Ts Ta Tmddp F2
  clear T; T=mktblprs(a,s,clin',   MDD_prs, ANX_prs, PTSD_prs,     minimal.x54_2_0, allobservations); %x25741=motion, x21000=race (1 2 3 4), x54=site
  mdl = fitlm(T,'y~x1+x1^2+x1*x2+x2+x3+x4 +x5+x6+x7 '); anovamdl=anova(mdl);
  p_perm_anova(i)=1-sum(anovamdl.F(3)>F_perm(i,:))/permutations; 
  p_perm_prs(i,1)=1-sum(anovamdl.F(4)> F_mddprs(i,:))/permutations; p_perm_prs(i,2)=1-sum(anovamdl.F(5)> F_anxprs(i,:))/permutations; p_perm_prs(i,3)=1-sum(anovamdl.F(6)> F_ptsdprs(i,:))/permutations;  % check this line if you're not including PRS as covariates
  Tprs(1,i)=mdl.Coefficients.tStat(8);  Tprs(2,i)=mdl.Coefficients.tStat(9);  Tprs(3,i)=mdl.Coefficients.tStat(10);  
  Tage(1,i)=mdl.Coefficients.tStat(2); 
  Tsex(1,i)=mdl.Coefficients.tStat(3);
  Tdep(1,i)=mdl.Coefficients.tStat(4); Pdep(1,i)=mdl.Coefficients.pValue(4);
  Tdepanx(1,i)=mdl.Coefficients.tStat(5); Pdepanx(1,i)=mdl.Coefficients.pValue(5);
  Tstr(1,i)=mdl.Coefficients.tStat(6); Pstr(1,i)=mdl.Coefficients.pValue(6);
  Tanx(1,i)=mdl.Coefficients.tStat(7); Panx(1,i)=mdl.Coefficients.pValue(7);
end
%figure;imagesc(Tanx);colorbar; colormap default

%%%Create table with summary scores for sign regions for dep (MDD), depanx (MDD+ANX), anx (ANX) or str (STR):
tmp=Pdep; tmp2=Tdep;
table(label_names_all(p_perm_anova<0.05 & tmp<0.0125),...
label_names_HOA(p_perm_anova<0.05 & tmp<0.0125), ...
tmp2(p_perm_anova<0.05 & tmp<0.0125)', tmp(p_perm_anova<0.05 & tmp<0.0125)',...
2*tmp2(p_perm_anova<0.05 & tmp<0.0125)'./sqrt(30419)    ); clear tmp tmp2; sortrows(ans, 2)


%create table with summary scores for 
tmp=3
table(label_names_all(p_perm_prs(:,tmp)<0.05),...
label_names_HOA(p_perm_prs(:,tmp)<0.05),...
Tprs(tmp,p_perm_prs(:,tmp)<0.05)', p_perm_prs(p_perm_prs(:,tmp)<0.05,tmp),...
2*Tprs(tmp,p_perm_prs(:,tmp)<0.05)'./sqrt(28189)); sortrows(ans, 2)

figure(3);subplot(2,3,1); scatter(Tdep', Tanx'); corr(Tdep', Tanx')
subplot(2,3,2);scatter(Tdep', Tstr'); corr(Tdep', Tstr')
subplot(2,3,3);scatter(Tanx', Tstr'); corr(Tanx', Tstr')
subplot(2,3,4); scatter(Tdep', Tdepanx'); corr(Tdep', Tdepanx')
subplot(2,3,5); scatter(Tanx', Tdepanx'); corr(Tanx', Tdepanx')
subplot(2,3,6); scatter(Tstr', Tdepanx'); corr(Tstr', Tdepanx')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Morphometric
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Similarity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% out of iterest
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear T* P*
for i=1:length(ms_ordered(1,:))
        allobservations=ms_ordered(:,i);allobservations(ix)=[];
        %parfor n = 1:permutations; 
        %permutation_index = randperm(length(allobservations));
        %randomSample = allobservations(permutation_index,:);
        %T=mktbl(a,s,clin',MDD_prs,minimal.x54_2_0, randomSample); %x25741=motion, x21000=race (1 2 3 4), x54=site
        %mdl = fitlm(T,'y~x1+x1^2+x1*x2+x2+x3+x4+x5');Td(1,n)=mdl.Coefficients.tStat(4); Tda(1,n)=mdl.Coefficients.tStat(5); Ts(1,n)=mdl.Coefficients.tStat(6); Ta(1,n)=mdl.Coefficients.tStat(7);Tmddp(1,n)=mdl.Coefficients.tStat(8);
        %anovamdl=anova(mdl); F(1,n)=anovamdl.F(3);F2(1,n)=anovamdl.F(4);
        %end; 
        %F_perm(i,:)=F; Tperm_dep(i,:)=Td;Tperm_depan(i,:)=Tda;Tperm_str(i,:)=Ts;Tperm_anx(i,:)=Ta;Tperm_mddprs(i,:)=Tmddp;
        %F_mddprs(i,:)=F2;clear Td Tda Ts Ta Tmddp F2
  clear T; T=mktbl(a,s,clin',MDD_prs,minimal.x54_2_0, allobservations); %x25741=motion, x21000=race (1 2 3 4), x54=site
  mdl = fitlm(T,'y~x1+x1^2+x1*x2+x2+x3+x4+x5'); anovamdl=anova(mdl);
  %p_perm_anova(i)=1-sum(anovamdl.F(3)>F_perm(i,:))/permutations; p_perm_mddprs(i)=1-sum(anovamdl.F(4)> F_mddprs(i,:))/permutations; 
  Tage(1,i)=mdl.Coefficients.tStat(2); Tmddprs(1,i)=mdl.Coefficients.tStat(8);  Pmddprs(1,i)=mdl.Coefficients.pValue(8); 
  Tsex(1,i)=mdl.Coefficients.tStat(3);
  Tdep(1,i)=mdl.Coefficients.tStat(4); Pdep(1,i)=mdl.Coefficients.pValue(4);
  Tdepanx(1,i)=mdl.Coefficients.tStat(5); Pdepanx(1,i)=mdl.Coefficients.pValue(5);
  Tstr(1,i)=mdl.Coefficients.tStat(6); Pstr(1,i)=mdl.Coefficients.pValue(6);
  Tanx(1,i)=mdl.Coefficients.tStat(7); Panx(1,i)=mdl.Coefficients.pValue(7);
end
figure;histogram(Tdep); hold on; histogram(Tdepanx); hold on;histogram(Tanx); hold on;histogram(Tstr);

figure;subplot(2,3,1); scatter(Tdep', Tanx'); corr(Tdep', Tanx')
subplot(2,3,2);scatter(Tdep', Tstr'); corr(Tdep', Tstr')
subplot(2,3,3);scatter(Tanx', Tstr'); corr(Tanx', Tstr')
subplot(2,3,4); scatter(Tdep', Tdepanx'); corr(Tdep', Tdepanx')
subplot(2,3,5); scatter(Tanx', Tdepanx'); corr(Tanx', Tdepanx')
subplot(2,3,6); scatter(Tstr', Tdepanx'); corr(Tstr', Tdepanx')

%% cognitive_ordered outcomes
cognitive_ordered.x6350_2_0(cognitive_ordered.x6350_2_0> (nanmean(cognitive_ordered.x6350_2_0)+4*nanstd(cognitive_ordered.x6350_2_0)) | cognitive_ordered.x6350_2_0<100)=NaN; %Duration to complete alphanumeric path 
cognitive_ordered.x6348_2_0(cognitive_ordered.x6348_2_0<100 | cognitive_ordered.x6348_2_0> (nanmean(cognitive_ordered.x6350_2_0)+3*nanstd(cognitive_ordered.x6348_2_0)) )=NaN; %Duration to complete numeric/easy path 
cognitive_ordered.tmt_cor=(cognitive_ordered.x6350_2_0 + 5*cognitive_ordered.x6351_2_0);% - (cognitive_ordered.x6348_2_0 +5*cognitive_ordered.x6349_2_0) ;
%ix=(cellfun('isempty', clinical)' | isnan(ica_partial(:,1))); a=age; a(ix)=[];s=sex;s(ix)=[];clin=clinical; clin(ix)=[];
%cognitive=cognitive_ordered; cognitive(ix,:)=[];minimal=minimal_ordered; minimal(ix,:)=[]; MDD_prs=MDD_prs_ordered; MDD_prs(ix)=[];
%% Cogntive data analysis (effects of group and PRS):
ix=(cellfun('isempty', clinical)' | isnan(MDD_prs_ordered)); a=age; a(ix)=[];s=sex;s(ix)=[];clin=clinical; clin(ix)=[];
cognitive=cognitive_ordered; cognitive(ix,:)=[];minimal=minimal_ordered; minimal(ix,:)=[]; MDD_prs=MDD_prs_ordered; MDD_prs(ix)=[]; ANX_prs=ANX_prs_ordered; ANX_prs(ix)=[]; PTSD_prs=PTSD_prs_ordered; PTSD_prs(ix)=[];
%%% alphanum trails
T=table(a,s,clin',MDD_prs, PTSD_prs, ANX_prs, minimal.x54_2_0, cognitive.tmt_cor);
mdl = fitlm(T,'Var8~a+s+a^2+a*s+Var3+MDD_prs+ANX_prs+PTSD_prs+Var7')
tmp=mdl.Coefficients{4:10,:};
%%% fluid intelligence
T=table(a,s,clin',MDD_prs, PTSD_prs, ANX_prs, minimal.x54_2_0,cognitive.x20016_2_0);
mdl = fitlm(T,'Var8~a+s+a^2+a*s+Var3+MDD_prs+ANX_prs+PTSD_prs+Var7')
tmp=horzcat(tmp, mdl.Coefficients{4:10,:});
%%% word pairs - PAL 
T=table(a,s,clin', MDD_prs, PTSD_prs, ANX_prs, minimal.x54_2_0,cognitive.x20197_2_0);
mdl = fitlm(T,'Var8~a+s+a^2+a*s+Var3+MDD_prs+ANX_prs+PTSD_prs+Var7')
tmp=horzcat(tmp, mdl.Coefficients{4:10,:});
%%% symbol digit
T=table(a,s,clin',MDD_prs, PTSD_prs, ANX_prs, minimal.x54_2_0, cognitive.x23324_2_0);
mdl = fitlm(T,'Var8~a+s+a^2+a*s+Var3+MDD_prs+ANX_prs+PTSD_prs+Var7')
tmp=horzcat(tmp, mdl.Coefficients{4:10,:});

%%% pairs matching (Number of incorrect matches in round) - 399.1-2 
T=table(a,s,clin',MDD_prs, PTSD_prs, ANX_prs, minimal.x54_2_0, cognitive.x399_2_2+cognitive.x399_2_1);
mdl = fitlm(T,'Var8~a+s+a^2+a*s+Var3+MDD_prs+ANX_prs+PTSD_prs+Var7')
tmp=horzcat(tmp, mdl.Coefficients{4:10,:});
%%% reaction time 
T=table(a,s,clin',MDD_prs, PTSD_prs, ANX_prs, minimal.x54_2_0, clean(cognitive.x20023_2_0, 3));
mdl = fitlm(T,'Var8~a+s+a^2+a*s+Var3+MDD_prs+ANX_prs+PTSD_prs+Var7')
tmp=horzcat(tmp, mdl.Coefficients{4:10,:});
%%% digit span
T=table(a,s,clin',MDD_prs, PTSD_prs, ANX_prs, minimal.x54_2_0, clean(cognitive.x4282_2_0,3));
mdl = fitlm(T,'Var8~a+s+a^2+a*s+Var3+MDD_prs+ANX_prs+PTSD_prs+Var7')
tmp=horzcat(tmp, mdl.Coefficients{4:10,:});
%%% Raven's matrices
T=table(a,s,clin',MDD_prs, PTSD_prs, ANX_prs, minimal.x54_2_0, cognitive.x6373_2_0);
mdl = fitlm(T,'Var8~a+s+a^2+a*s+Var3+MDD_prs+ANX_prs+PTSD_prs+Var7')
tmp=horzcat(tmp, mdl.Coefficients{4:10,:});
%%% Tower of london - nr of figures done 21004
T=table(a,s,clin', MDD_prs, PTSD_prs, ANX_prs, minimal.x54_2_0, cognitive.x21004_2_0);
mdl = fitlm(T,'Var8~a+s+a^2+a*s+Var3+MDD_prs+ANX_prs+PTSD_prs+Var7')
tmp=horzcat(tmp, mdl.Coefficients{4:10,:});

%%% neuroticism
T=table(a,s,clin', MDD_prs, PTSD_prs, ANX_prs, minimal.x54_2_0,  cognitive.x20127_0_0);
mdl = fitlm(T,'Var8~a+s+a^2+a*s+Var3+MDD_prs+ANX_prs+PTSD_prs+Var7')
tmp=horzcat(tmp, mdl.Coefficients{4:10,:});

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

%[~,score,~,~,explained]  = pca([MDD_prs_ordered, ANX_prs_ordered, PTSD_prs_ordered],'NumComponents',3);MDD_prs_ordered=score(:,2);

%% PRS analysis - between PRS correlations
%% PRS analysis:
ix=(cellfun('isempty', clinical)' | isnan(MDD_prs_ordered) | isnan(ANX_prs_ordered) | isnan(PTSD_prs_ordered)); a=age; a(ix)=[];s=sex;s(ix)=[];clin=clinical; clin(ix)=[];
cognitive=cognitive_ordered; cognitive(ix,:)=[];minimal=minimal_ordered; minimal(ix,:)=[]; 
MDD_prs=MDD_prs_ordered; MDD_prs(ix)=[]; ANX_prs=ANX_prs_ordered; ANX_prs(ix)=[]; PTSD_prs=PTSD_prs_ordered; PTSD_prs(ix)=[];
corr([MDD_prs, ANX_prs, PTSD_prs])

index=~(  ismember(clin','hc')   ); %
corr([MDD_prs(index), ANX_prs(index), PTSD_prs(index)])
