function [N, eps0, eps1, kn] = get_sample_number(DI_target, Pout_target, snr, VarX, alpha, mode)
% Inputs
% DI_target: target error
% Pout_target: target outage probability
% snr: signal-noise-ratio
% VarX: variance of X
% alpha: |X| is alpha-sub-Gaussian
% Outputs
% N: number of samples
% eps0: epsilon_0
% eps1: epsilon_1

    N = zeros(2,1);
    eps0 = zeros(2,1);
    eps1 = zeros(2,1);
    kn = zeros(2,1);

    for thm = 5:6
        func = @(log_eps0) get_N_eps0_DI(thm, log_eps0, DI_target, Pout_target, snr, VarX, alpha, mode);
        log_eps0_opt = fminsearch(func, -100);

        N(thm-4) = get_N_eps0_DI(thm, log_eps0_opt, DI_target, Pout_target, snr, VarX, alpha, mode);
        eps0(thm-4) = exp(log_eps0_opt);
        eps1(thm-4) = get_eps1_DI(thm, log_eps0_opt, DI_target, snr, VarX, alpha, mode);
        kn(thm-4) = get_opt_k(thm, eps0(thm-4), eps1(thm-4), snr, VarX, alpha, mode);
    end
    
end


% get sample number to attain targeted DI and Pout given eps0
function N = get_N_eps0_DI(theorem, log_eps0, DI_target, Pout_target, snr, VarX, alpha, mode)
    eps1 = get_eps1_DI(theorem, log_eps0, DI_target, snr, VarX, alpha, mode);
    if isnan(eps1) 
        N = inf;
    else
        eps0 = exp(log_eps0);
        func = @(n) get_Pout(n, eps0, eps1) - Pout_target;
        N = fzero(func, 1e10);
    end
end

% get eps1 such that DI is of desired width for a given eps0
function eps1 = get_eps1_DI(theorem, log_eps0, DI_target, snr, VarX, alpha, mode)   
    fun = @(log_eps1) DI(theorem, log_eps0, log_eps1, snr, VarX, alpha, mode) - DI_target;
    fun_min = fun(-100);
    if isnan(fun_min) || fun_min > 0
        eps1 = NaN;
    else
        log_eps1 = fzero(fun, [-100, 0]);
        eps1 = exp(log_eps1);
    end
end


% get probability of estimate being outside DI for given eps0, eps1
function Pout = get_Pout(n, eps0, eps1)
    n = max(n,1);
    
    c0 = 1/sqrt(2*pi*exp(1));
    c1 = (2/exp(1)+1)/sqrt(2*pi);

    v0 = sqrt(2/pi);
    v1 = sqrt(2/(pi*exp(1)));

    a0 = (2/3)*eps0/c0;
    a1 = (2/3)*eps1/c1;
            
    Pout0 = 2*exp(-2*n*a0^2*(eps0-c0*a0)^2/v0^2);
    Pout1 = 2*exp(-2*n*a1^4*(eps1-c1*a1)^2/v1^2);
    
    Pout = Pout0 + Pout1;
end

% get DI width for given eps0, eps1
function val = DI(theorem, log_eps0, log_eps1, snr, VarX, alpha, mode)
    eps0 = exp(log_eps0);
    eps1 = exp(log_eps1);
    k_opt = get_opt_k(theorem, eps0, eps1, snr, VarX, alpha, mode);
    val = DI_k(theorem, eps0, eps1, snr, VarX, alpha, mode, k_opt);
end

% get optimal k
function k_opt = get_opt_k(theorem, eps0, eps1, snr, VarX, alpha, mode)
    func = @(k) DI_k(theorem, eps0, eps1, snr, VarX, alpha, mode, k);
    if theorem < 6 
        if -log(sqrt(2*pi)*eps0)-snr*VarX <= 0
            k_opt = NaN;
            return
        else
            ub = sqrt(-log(sqrt(2*pi)*eps0)-snr*VarX);
        end
    elseif theorem == 6
        ub = 10*sqrt(VarX);
    else
        error('Argument "theorem" must be 5 or 6');
    end
       
    k_opt = fminbnd(func, 0, ub);
end

% DI as a function of truncation width k
function val = DI_k(theorem, eps0, eps1, snr, VarX, alpha, mode, k)
% f0: maximum of f(t) in [-k,k]
% f1: maximum of f'(t) in [-k,k]
% df: number of modes of f_Y
    f = @(t) fY(snr,t,mode);
    Df = @(t) dfY(snr,t,mode);
    t_df0 = fzero(Df, 0);
    df = size(t_df0,2); % # modes
    f0 = max(f(t_df0));
    
    rhox = rho_max(k, snr, VarX);
    phix = phi(k, snr, VarX);
    ck = c(k, snr, alpha);
    if theorem == 5
        if eps0*phix < 1
            val = (4*eps1*k*rhox + eps0*phix + 2*eps1^2*k*phix)/(1-eps0*phix)+ck;
        else
            val = 1;
        end
    elseif theorem == 6
        psix = psi(0, k, snr, VarX, f0);
        rhox = rho_max(k, snr, VarX);
        val1 = 4*eps1*min((2+df)*psix,rho_max_int(k,snr,VarX));
        val1 = val1 + 2*eps0*min((2+df)*psix*rhox,rho_max_squared_int(k,snr,VarX)) + ck;
%         val1 = 2*(2+df)*psix*(2*eps1+ eps0*rhox) + ck;
        val2 = 2*eps1*rho_max_int(k,snr,VarX)+2*eps0*rho_max_squared_int(k,snr,VarX);
        val = max(val1,val2);
    else
        error('Argument "theorem" must be 5 or 6');
    end
end


% auxiliray functions defined in the paper
function val = rho_max(t, snr, VarX)
    val = sqrt(3*snr*VarX) + 3*abs(t);
end

function val = rho_max_int(t, snr, VarX)
    a = sqrt(3*snr*VarX);
    b = 3;
    val = a*t + (1/2)*b*t^2;
end

function val = rho_max_squared_int(t, snr, VarX)
    a = sqrt(3*snr*VarX);
    b = 3;
    val = a^2*t + a*b*t^2 + (1/3)*b^2*t^3;
end

function val = phi(t, snr, VarX)
    val = sqrt(2*pi)*exp(t^2+snr*VarX);
end

function val = psi(eps0, t, snr, VarX, f0)
    val1 = log(f0+eps0);
    val2 = log(phi(t, snr, VarX)/(1-eps0*phi(t, snr, VarX)));
    val = max(val1, val2);
end

function val = c(t, snr, alpha)
    factor = @(v)2*gamma(v+0.5)^(1/(1+v))/(pi^(1/(2*(1+v))));
    moment = @(v) (2*exp((alpha^2*snr-t^2)/2))^(v/(v+1));
    fun = @(v) factor(v)*moment(v);
    v_opt = fminbnd(fun, 0, 10);
    val = fun(v_opt);
end
