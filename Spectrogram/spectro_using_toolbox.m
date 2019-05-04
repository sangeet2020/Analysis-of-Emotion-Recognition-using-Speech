% function [e1, e2,e3] = spectro_using_toolbox(signal,fs,f)

[signal,fs]=audioread('03a02Wb.wav');
f=227;  % f is number of frames
envelope = imdilate(abs(signal), true(1501, 1));
% Find the quiet parts.
quietParts = envelope < 0.06; % Or whatever value you want.
% Cut out quiet parts and plot.
yEdited = signal; % Initialize
yEdited(quietParts) = [];
signal=yEdited(:);

[T,F,B]=spgrambw(signal,fs,'pjcw',200,[0 8000],[-130, -22]); % This function is from Voicebox %Toolkit.
h=imagesc(T,F,10*log10(B'));colormap('jet');
shading interp;
axis xy;
axis off;
saveas(h,'spectro.png');

img = imread('spectro.png'); % Read image
img1 = imcrop(img,[155.5 67.5 931 732]);
img2=imresize(img1,[f f]);
% imshow(img2);
red = img2(:,:,1); % Red channel
green = img2(:,:,2); % Green channel
blue = img2(:,:,3); % Blue channel
a = zeros(size(img2, 1), size(img2, 2));
just_red = cat(3, red, a, a);
just_green = cat(3, a, green, a);
just_blue = cat(3, a, a, blue);

%% For figure plotting 
figure;

subplot(2,3,1); %[55.5 39.5 738 583]
t1=0:1/fs:(length(yEdited)-1)/fs;
plot(t1,yEdited);
set(gca,'FontWeight','bold');
title({'Original angry' ,'speech signal'},'FontSize',18);
xlabel('Time (s)','fontweight','bold','fontsize',16); ylabel('Amplitude','fontweight','bold','fontsize',16);
grid on;

ax1 = subplot(2,3,2);
imagesc(T,F,10*log10(B'));%colormap('jet');
colormap(ax1,jet);
shading interp;
[c1 ,d1,e1] =size(img);
title(['RBG Image Array ',num2str(c1),'x',num2str(d1),' '],'FontSize',18);
set(gca,'FontWeight','bold');
xlabel('Time (s)','fontweight','bold','fontsize',16); 
ylabel('Frequency (Hz)','fontweight','bold','fontsize',16);
axis xy;
% img = imcrop(img1,[55.5 39.5 738 583]);
% imshow(img);


ax2=subplot(2,3,3);
imagesc(T,F,10*log10(B'));%colormap('jet');
colormap(ax2,flipud(gray));
shading interp;
[c1 ,d1,e1] =size(img);
title(['Grayscale Image Array ',num2str(c1),'x',num2str(d1),' '],'FontSize',18);
set(gca,'FontWeight','bold');
xlabel('Time (s)','fontweight','bold','fontsize',16); 
ylabel('Frequency (Hz)','fontweight','bold','fontsize',16);
axis xy;


subplot(2,3,4);
im_yellow=imshow(just_red);
[c ,d] =size(red);
title(['R – Array ',num2str(c),'x',num2str(d)],'FontSize',18);
axis on;
% set(gca,'XTickLabel',{0,0.5,1,1.5});
% % yticks([0 2000 4000 6000 8000])
set(gca,'YTickLabel',{'8000','6000','4000','2000','0'},...
    'YTick',[  1.0000   57.5000  114.0000  170.5000  227.0000],...
     'XTickLabel',{'0.5','1','1.5'},...
    'XTick',[60   120   180],...
    'Layer','top',...
    'TickDir','out',...
    'YDir','reverse',...
    'FontWeight','bold',... 
    'DataAspectRatio',[1 1 1]);

xlabel('Time (s)','fontweight','bold','fontsize',16); 
ylabel('Frequency (Hz)','fontweight','bold','fontsize',16);

subplot(2,3,5);
im_green=imshow(just_green);
[c ,d] =size(green);
title(['G – Array ',num2str(c),'x',num2str(d)],'FontSize',18);
axis on;
set(gca,'YTickLabel',{'8000','6000','4000','2000','0'},...
    'YTick',[  1.0000   57.5000  114.0000  170.5000  227.0000],...
     'XTickLabel',{'0.5','1','1.5'},...
   'XTick',[60   120   180],...
    'Layer','top',...
    'TickDir','out',...
    'YDir','reverse',...
    'FontWeight','bold',...
    'DataAspectRatio',[1 1 1]);
xlabel('Time (s)','fontweight','bold','fontsize',16); 
ylabel('Frequency (Hz)','fontweight','bold','fontsize',16);

subplot(2,3,6);
im_blue=imshow(just_blue);
[c ,d] =size(blue);
title(['B – Array ',num2str(c),'x',num2str(d)],'FontSize',18);
axis on;
set(gca,'YTickLabel',{'8000','6000','4000','2000','0'},...
   'YTick',[  1.0000   57.5000  114.0000  170.5000  227.0000],...
    'XTickLabel',{'0.5','1','1.5'},...
    'XTick',[60   120   180],...
    'Layer','top',...
    'TickDir','out',...
    'FontWeight','bold',...
    'YDir','reverse',...
    'DataAspectRatio',[1 1 1]);
xlabel('Time (s)','fontweight','bold','fontsize',16); 
ylabel('Frequency (Hz)','fontweight','bold','fontsize',16);


%% 
% just_red=double(just_red);
% just_red=just_red(:,:,1);
% normR = just_red - min(just_red(:));
% my_red = normR./ max(normR(:));
% 
% just_green=double(just_green);
% just_green=just_green(:,:,2);
% normG = just_green - min(just_green(:));
% my_green = normG./ max(normG(:));
% 
% just_blue=double(just_blue);
% just_blue=just_blue(:,:,3);
% normB = just_blue - min(just_blue(:));
% my_blue = normB./ max(normB(:));
% 
% %% Eigen Values
% e1=abs(eig(my_red)); 
% e2=abs(eig(my_green));
% e3=abs(eig(my_blue));

%% PCA
% [~, ~, e1]=princomp(my_red'); 
% [~, ~, e2]=princomp(my_green');
% [~, ~, e3]=princomp(my_blue');
