A=imread('cell.tif');                     %读取并显示图像
SE=strel('disk',4,4);                            %定义模板
B=imdilate(A,SE);                             %按模板膨胀
C=imerode(A,SE);                             %按模板腐蚀
figure
subplot(131),imshow(A);
subplot(132),imshow(B);
subplot(133),imshow(C);
