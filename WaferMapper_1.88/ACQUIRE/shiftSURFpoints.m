function p=shiftSURFpoints(p,shiftx,shifty) 

if ~isempty(p)
    for k=1:length(p)
       p(k).x=p(k).x+shiftx;
       p(k).y=p(k).y+shifty;
    end
end
