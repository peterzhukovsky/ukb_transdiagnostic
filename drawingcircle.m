%https://uk.mathworks.com/matlabcentral/answers/6322-drawing-a-segment-of-a-circle
%drawing a circle
%example P = circle_arc(pi/4,3*pi/3,9,-4,3);
ica2yeo7=readtable('ica2yeo7.csv');


range=[0:0.299:2*pi]
for i=1:21
    try
    P = circle_arc(range(i),range(i+1),9,-4,3);hold on;
    set(P,'edgecolor','w','linewidth',2, 'facecolor', ica2yeo7.Color{i})
    end
end

%radians go between 0 and 2pi
names={'vol0000';'vol0001';'vol0002';'vol0003';'vol0004';'vol0005';'vol0006';'vol0007';'vol0008';'vol0009';'vol0010';'vol0011';'vol0012';'vol0013';'vol0014';'vol0015';'vol0016';'vol0017';'vol0018';'vol0019';'vol0020'};
for i=1:length(names)
tmp=dlmread(strcat('D:\Canada_2020\UK_biobank\data\ICA_d25\Yeo7Mapping\', names{i},'.txt'));
vols(i,:)=tmp(:,1)./(sum(tmp(:,1)));
end

%%%% import the excel manually adjusted

