function [ wx,wy,wz ] = getAvgVort(p)
%GETAVGVELUV reads from output/vel_uv_avg*.bin

for i=1:p.nproc
    
    % Open the file
    fname = ['./output/vort_avg.c',num2str(i-1),'.bin'];
    fid=fopen(fname,'r');
    if (fid < 0) 
        error('getSnap:fname',['Could not open file ',fname]);
    end

    % Determine the interval of the matrix where the data should be stored
    zmin=p.zmin_buf(i);
    zmax=p.zmax_buf(i);
    
    % Scan the data
    N = p.nx*p.ny*p.nz2;
    dummy=fread(fid,N, 'double',p.fmt);
    wx(1:p.nx,1:p.ny,zmin:zmax)=reshape(dummy,p.nx,p.ny,p.nz2);
    dummy=fread(fid,N, 'double',p.fmt); 
    wy(1:p.nx,1:p.ny,zmin:zmax)=reshape(dummy,p.nx,p.ny,p.nz2);
    dummy=fread(fid,N, 'double',p.fmt); 
    wz(1:p.nx,1:p.ny,zmin:zmax)=reshape(dummy,p.nx,p.ny,p.nz2);
    
    fclose(fid);
end

end

