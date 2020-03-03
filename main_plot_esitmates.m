clear all;
clc;

%% inputs
mode = 0; % 0: Guassian input, 1: Binary input
varX = 1;

%% Density function, its derivative, and score function with varying y and n but fixed snr
snr0 = 1;

py = 100; % # points of y in the plots of pdf, its derivative, and score function vs y
y = linspace(-10,10,py);   % vector of y
% The functions
[f,df,rho] = ScoreFun(snr0,y,mode);

nv = [100,500,1000,5000,10000];

nn = size(nv,2);
rhon = zeros(nn,py);
fn = zeros(nn,py);
dfn = zeros(nn,py);
for cnt = 1:nn
    n = nv(cnt);
    
    a0 = n^(-1/8); % bandwidth of estimator
    a1 = n^(-1/8); % bandwidth of estimator
    
    % samples
    if mode==0
        Xi = randn(1,n); % Gaussian
    else
        Xi = 2*(rand(1,n)>.5) - 1; % Binary
    end
    Zi = randn(1,n); % Gaussian noise
    Yi0 = sqrt(snr0)*Xi + Zi;
    
    % The estimates
    % fn: estimator of fY, the pdf of Y
    % dfn: estimator of dfY, the derivative of fY
    % rhon: estimator of rho, the score function
    [rhon(cnt,:),fn(cnt,:),dfn(cnt,:)] = EstScoreFun(Yi0,y,a0,a1);

end
figure
hold on
plot(y,f,'-k')
plot(y,fn,'--')
xlabel('y')
ylabel('f_Y')
legend('f_Y','f_n,n=50','f_n,n=100','f_n,n=500','f_n,n=1000','f_n,n=5000')
grid on
figure
hold on
plot(y,df,'-k')
plot(y,dfn,'--')
xlabel('y')
ylabel('f`_Y(y)')
legend('df_Y','df_n,n=50','df_n,n=100','df_n,n=500','df_n,n=1000','df_n,n=5000')
grid on

%% FI and MMSE with varing snr
psnr = 25; % # points of snr in the plots of Fisher information vs snr
snr = linspace(0,10,psnr); % vector of snr
n = 1e4;
% Parameters of estimator: 
a0 = 0.6;%n^(-1/6); % bandwidth of estimator
a1 = 0.6;%n^(-1/8); % bandwidth of estimator
kn = 10;

if mode==0
    Xi = randn(1,n); % Gaussian
else
    Xi = 2*(rand(1,n)>.5) - 1; % Binary
end
Zi = randn(1,n); % Gaussian noise
J = zeros(1,psnr);
Jn_org = zeros(1,psnr);
Jn_reg = zeros(1,psnr);
for cnt = 1:psnr
    cnt
    snri = snr(cnt);        
    
    Yi = sqrt(snri)*Xi + Zi; % psnr*n matrix
    
    J(cnt) = FisherInfo(snri,mode);
    
    % Jn: estimator of J, the Fisher information
    [Jn_org(cnt)] = EstFisherInfo(Yi,a0,a1,kn);
    [Jn_reg(cnt)] = RegularizedEstFI(Yi,a0,a1,kn);

end
% MMSEn: estimator of MMSE
MMSE = (1-J)./snr;
% MMSEg = 1./(1+snr); % MMSE for Gaussian
MMSEn_org = (1-Jn_org)./snr;
MMSEn_reg = (1-Jn_reg)./snr;


figure
hold on
plot(snr,J,'-k')
plot(snr,Jn_org,'--xk')
plot(snr,Jn_reg,'-.ok')
hold off
xlabel('snr')
ylabel('I(snr)')
legend('I','I_n','I^c_n')
grid on

figure
hold on
plot(snr,MMSE,'-k')
plot(snr,MMSEn_org,'--xk')
plot(snr,MMSEn_reg,'-.ok')
hold off
xlabel('snr')
ylabel('MMSE(snr)')
legend('MMSE','MMSE_n','MMSE^c_n')
grid on