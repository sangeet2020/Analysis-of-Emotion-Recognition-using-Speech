close all;
clear all;
clc;
tr0=load('test_savee_mfcc_delta_doubledelta.m');
sh_cf=zeros(6,1);%no of emotions
u1=0;u2=0;u3=0;u4=0;u5=0;u6=0;
for i=1:720%6emotions*20frames/utterance*no of samples in one emotion taken for testing=no of samples taken for testing
x=tr0(i,:);
x=x';
k=1;
l=1;
p=load('savee60_gmm_model_mfcc_delta_doubledelta.m');
 for j=1:6% no of emotions
 mus=p(k:k+15,:);
 sigm=p(k+16:k+31,:);
 c=p(k+32,:);
 c=c';
 k=k+33;
 [YM,Y]=lmultigauss_mod(x,mus,sigm,c);
 y(l)=Y;
 l=l+1;
 end

 [log Ind]=max(y);

 if Ind==1
      u1=u1+1;
     sh_cf(1)=u1;
 elseif Ind==2
      u2=u2+1;
     sh_cf(2)=u2;
 elseif Ind==3
      u3=u3+1;
     sh_cf(3)=u3;
 elseif Ind==4
      u4=u4+1;
     sh_cf(4)=u4;
 elseif Ind==5
      u5=u5+1;
     sh_cf(5)=u5;
 elseif Ind==6
      u6=u6+1;
     sh_cf(6)=u6;
end
end
  fprintf(' emotion 1 is identified=') 
 disp(sh_cf(1,1))
 fprintf(' emotion 2 is identified=') 
 disp(sh_cf(2,1))
 fprintf(' emotion 3 is identified=') 
 disp(sh_cf(3,1))
 fprintf(' emotion 4 is identified=') 
 disp(sh_cf(4,1))
 fprintf(' emotion 5 is identified=') 
 disp(sh_cf(5,1))
 fprintf(' emotion 6 is identified=') 
 disp(sh_cf(6,1))

  
if sh_cf(1,1)<120% no of frames taken p=20 (in pattDID_new_2.m)* no of samples in one emotion taken for testing
    a1=120-sh_cf(1,1);
else
    a1=0;
end
if sh_cf(2,1)<120
    a2=120-sh_cf(2,1);
else
    a2=0;
end
 if sh_cf(3,1)<120
    a3=120-sh_cf(3,1);
else
    a3=0;
 end
if sh_cf(4,1)<120
    a4=120-sh_cf(4,1);
else
    a4=0;
end
if sh_cf(5,1)<120
    a5=120-sh_cf(5,1);
else
    a5=0;
if sh_cf(6,1)<120
    a6=120-sh_cf(6,1);
else
    a6=0;
end
end

%identification_rate= [(no of samples taken for tesing-no of mis
%classifcation)/no of samples taken for tesing]*100=
%[(700-(a1+a2+a3+a4+a5+a6+a7))*100]/700=(700-(a1+a2+a3+a4+a5+a6+a7))/7

   
Identification_rate=(720-(a1+a2+a3+a4+a5+a6))/7.2;
disp(Identification_rate);
%      
% %  