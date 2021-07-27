%No. Label Name:				R   G   B   A
%% RIGHT SIDE first
clear full_cell
id=[1:180]';  tstat_roi=Tstr; %tstat_roi=Tprs(3,:);
for roi=1:length(rlabel_names)
    %e.g. t=2; maxt=4 => then color should be half way to the max of the range:
    if tstat_roi(roi+180)>0
        blue_range=0:255;
        green_range=178.5:0.3:255;
    
    tstat_relative=1-tstat_roi(roi+180)/(max(tstat_roi)+0.1);
    color_relative=round(tstat_relative*256);
    
    R(roi,1)=255;
    G(roi,1)=round(green_range(color_relative));
    B(roi,1)=blue_range(color_relative);
    A(roi,1)=0;
    
    full_cell{roi, 1}=id(roi); full_cell{roi, 2}=rlabel_names{roi};full_cell{roi, 3}=R(roi);    full_cell{roi, 4}=G(roi);    full_cell{roi, 5}=B(roi);     full_cell{roi, 6}=A(roi);

    elseif tstat_roi(roi+180)<0;
        red_range=0:255;
        green_range=183.5:0.28:255;
    tstat_relative=1-tstat_roi(roi+180)/(min(tstat_roi)-0.1);
    color_relative=round(tstat_relative*256);
    
    
    R(roi,1)=red_range(color_relative);
    G(roi,1)=round(green_range(color_relative));
    B(roi,1)=255;
    A(roi,1)=0;
    full_cell{roi, 1}=id(roi); full_cell{roi, 2}=rlabel_names{roi};full_cell{roi, 3}=R(roi);    full_cell{roi, 4}=G(roi);    full_cell{roi, 5}=B(roi);     full_cell{roi, 6}=A(roi);
    end
end



%for the blue scales: blue 255 always; green 85:255 red 0:255

%% LEFT SIDE next
clear full_cell
id=[1:180]'; 
for roi=1:length(rlabel_names)
    %e.g. t=2; maxt=4 => then color should be half way to the max of the range:
    if tstat_roi(roi)>0
        blue_range=0:255;
        green_range=178.5:0.3:255;
    
    tstat_relative=1-tstat_roi(roi)/(max(tstat_roi)+0.1);
    color_relative=round(tstat_relative*256);
    
    R(roi,1)=255;
    G(roi,1)=round(green_range(color_relative));
    B(roi,1)=blue_range(color_relative);
    A(roi,1)=0;
    
    full_cell{roi, 1}=id(roi); full_cell{roi, 2}=llabel_names{roi};full_cell{roi, 3}=R(roi);    full_cell{roi, 4}=G(roi);    full_cell{roi, 5}=B(roi);     full_cell{roi, 6}=A(roi);

    elseif tstat_roi(roi)<0;
        red_range=0:255;
        green_range=183.5:0.28:255;
    tstat_relative=1-tstat_roi(roi)/(min(tstat_roi)-0.1);
    color_relative=round(tstat_relative*256);
    
    
    R(roi,1)=red_range(color_relative);
    G(roi,1)=round(green_range(color_relative));
    B(roi,1)=255;
    A(roi,1)=0;
    full_cell{roi, 1}=id(roi); full_cell{roi, 2}=llabel_names{roi};full_cell{roi, 3}=R(roi);    full_cell{roi, 4}=G(roi);    full_cell{roi, 5}=B(roi);     full_cell{roi, 6}=A(roi);
    end
end


%% Age hist

histogram(age(group==4))

%% Boxplot the adjr2

figure
boxplot([adjr1',adjr2', adjr3'],'Notch','on')

boxplot([log(adjr1(adjr1>0 & adjr2>0 & adjr3>0))', log(adjr2(adjr1>0 & adjr2>0 & adjr3>0))', log(adjr3(adjr1>0 & adjr2>0 & adjr3>0))'],'Notch','on')

d2=adjr2'-adjr1'
[h,p,ci,stats]=ttest(d2)
mean(d2)/std(d2)

d3=adjr3'-adjr1'
mean(d3)/std(d3)
[h,p,ci,stats]=ttest(d3)
figure; yyaxis right
scatter(ones(length(d2),1),d2, 'b', 'filled'); hold on
scatter(2*ones(length(d3),1),d3, 'b', 'filled');xlim([0 3]);