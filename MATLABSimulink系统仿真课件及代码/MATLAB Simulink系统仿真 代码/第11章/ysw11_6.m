clc,clear
close all
n = 0:0.01:2;
for i=1:4
    switch i
    case 1, N=2;
        case 2,N=5;
            case 3,N=10;
                case 4,N=20;
    end
    [z,p,k] = buttap(N);
    [b,a] = zp2tf(z,p,k);
    [H,w] = freqs(b,a,n);
    magH2 = (abs(H).^2);     %传递函数幅值平方
    hold on
    plot(w,magH2,'linewidth',2)
end
xlabel('w/wc');
ylabel('H(jw)|^2');
title('Butterworth低通模拟原型滤波器')
text(1.5,0.18,'N=2')
text(1.3,0.08,'N=5')
text(1.16,0.08,'N=10')
text(0.93,0.98,'N=20')
grid on








