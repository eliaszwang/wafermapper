function M=ConvertCompositeAffineTransform(CompTrans)


tdata=CompTrans.tdata;

if length(tdata)>1

    M=tdata(1).tdata.T;
    for i=2:length(tdata)
      M=tdata(i).tdata.T*M;
    end
else
    M=tdata.T;
end


