function[md5List] = checkSumFolder(TPN);

dirTPN = dir(TPN); dirTPN = dirTPN(3:end);
md5path = pwd;

md5List = cell(length(dirTPN),2);
for i = 1:length(dirTPN)
    for r = 1:3
        nam = dirTPN(i).name;
        s=sprintf('%c%s%s\t%s\n','!',md5path,'\md5.exe -n -otempmd5.txt',[TPN nam]);
        eval(s);
        success = 1;
        
        try
            md5=textread('tempmd5.txt','%s');
        catch err
            delete('tempmd5.txt','%s');
            success = 0;
        end
        
        if success
            md5List{i,1} = nam;
            md5List{i,2} = md5;
            break
        end
    end
    
    
end