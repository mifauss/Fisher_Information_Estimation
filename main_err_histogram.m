clear all;
clc;

%% inputs
mode = 1; % 0: Guassian input, 1: Binary input
n = 1e3; % # samples
T = 1e3; % # experiments

% Parameters of estimator: 
a = n^(-1/6); % bandwidth of estimator
kn = log(n);

varX = 1;
snr = 1;


%% FI and MMSE with varing snr
J = FisherInfo(snr,mode); % Fisher info    
Jn = zeros(1,T);
Jr = zeros(1,T);
for cnt = 1:T
    cnt
    
    % samples
    if mode==0
        Xi = randn(1,n); % Gaussian
    else
        Xi = 2*(rand(1,n)>.5) - 1; % Binary
    end
    Zi = randn(1,n); % Gaussian noise
    Yi = sqrt(snr)*Xi + Zi;
    
    % Jn: estimator of J
    Jn(cnt) = EstFisherInfo(Yi,a,a,kn);
    Jr(cnt) = RegEstFisherInfo(Yi,a,a,kn,snr,varX);
end
% MMSE = 1-snr.*J;
% % MMSEg = 1./(1+snr); % MMSE for Gaussian
% MMSEn = 1-snr.*Jn; % estimator of MMSE

errn = abs(J-Jn);
errr = abs(J-Jr);


%% plots
figure
hold on
plot(J,0,'xk','MarkerSize',8,'LineWidth',2)
histogram(Jn,'Normalization','probability')%'BinLimits',[0.1,1],
histogram(Jr,'Normalization','probability')
xlabel('J_n')
ylabel('probability')
legend('J','J1','J2','J3')

figure
hold on
histogram(errn,'Normalization','probability')
histogram(errr,'Normalization','probability')
xlabel('err')
ylabel('probability')
legend('J1','J2','J3')
