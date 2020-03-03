clc;
clear all;

mode = 0;
varX = 1;
snr = 1;
alpha = 1;
nthm = 2;

% perr = 0.2;
% errv = linspace(0.1,0.9,9);
% na = zeros(nthm,9);
% eps0a = zeros(nthm,9);
% eps1a = zeros(nthm,9);
% kna = zeros(nthm,9);
% for cnt = 1:9
%     err = errv(cnt)
%     [n,eps0,eps1,kn] = get_sample_number(err,perr,snr,varX,alpha,mode);
%     na(:,cnt) = n;
%     eps0a(:,cnt) = eps0;
%     eps1a(:,cnt) = eps1;
%     kna(:,cnt) = kn;
% end
% figure
% plot(errv,log(na)/log(10))
% grid on
% xlabel('err')
% ylabel('ln n')
% legend('Theorem 5','Theorem 6')

perrv = linspace(0.1,0.9,9);
np = size(perrv,2);
err = 0.5;
nb = zeros(nthm,9);
for cnt = 1:9
    perr = perrv(cnt)
    n = get_sample_number(err,perr,snr,varX,alpha,mode);
    nb(:,cnt) = n.';
end
figure
hold on 
plot(perrv,log(nb)/log(10))
hold off
grid on
xlabel('Perr')
ylabel('ln n')
legend('Theorem 5','Theorem 6')

