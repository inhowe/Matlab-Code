A=imread('rice.png');                %读入并显示图像
B=fspecial('Sobel');                        %用Sobel算子进行边缘锐化
fspecial('Sobel');
B=B';                                       %Sobel垂直模板
C=filter2(B,A);
figure
subplot(121),imshow(A);                  %显示添加椒盐噪声后的图像
subplot(122),imshow(C);                  %显示平滑处理后图像
