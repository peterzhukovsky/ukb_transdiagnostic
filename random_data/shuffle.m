 function v=shuffle(v)
    if istable(v)
     v=v(randperm(height(v)),:);
    else
        v=v(randperm(length( v(:,1) )),:);
 end