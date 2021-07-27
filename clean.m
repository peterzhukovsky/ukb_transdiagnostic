function cleandata = clean(input, n)
 upper=nanmean(input)+n*nanstd(input);
 lower=nanmean(input)-n*nanstd(input);
 index= input>upper | input<lower;
 cleandata=input; cleandata(index==1)=NaN;
end
