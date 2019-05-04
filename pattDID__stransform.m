function  [A_lpcc]=pattDID__stransform(filename,l)
i=l-1;
% myfile =sprintf('%s%s%d%', 'C:\Users\Sangeet Sagar\Music\Dataset_Train_Test_Same','\',i)
myfile =sprintf('%s%s%d%', 'D:\BTP\EmodB data base\An_ EmoDB_Dataset','\',i)
Files=dir(myfile);
for k=3:length(Files)
    file=Files(k).name;
    fid =  fopen(filename,'a');
    disp(file);
    [signal, fs] = audioread(file);
    [A_lpcc]=lpcc_coeff(signal,fs);
     A_lpcc = A_lpcc(isfinite(A_lpcc(:, 1)), :);
    [m , ~]= size(A_lpcc);
    s=zeros(1,19);
    for p=1:m
        s= A_lpcc(p,:);
        fprintf(fid ,'%f ' , s);
        fprintf(fid, '\n');
    end
    fclose(fid);
end