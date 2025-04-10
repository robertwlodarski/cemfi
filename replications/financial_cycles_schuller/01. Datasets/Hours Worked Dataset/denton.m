function [x,lambda] = denton(z,y,a,t)
% PURPOSE: estimate high-frequency series x that are consistent with the 
% reliable low-frequency series y of the same variable and maximally 
% reproducing the movements in the indicator (inconsistent with y) high-
% frequency series z according to the Denton method.
% -------------------------------------------------------------------------
% USAGE: x = denton(z,y,a,t) or [x,lambda] = denton(z,y,a,t), where
% -------------------------------------------------------------------------
% INPUT:
% -> z = indicator/preliminary high-frequency (e.g. sub-annual) series (a 
%        column vector)
% -> y = benchmark or low-frequency (e.g. annual) series (a column vector)
% -> a = a scalar that defines the variants of the Denton procedure:
%       1 - additive level differences 
%       2 - proportional level differences 
%       3 - additive first differences
%       4 - proportional first differences
%       5 - additive second differences
%       6 - proportional second differences
% -> t - a scalar defining the following types of the Denton method:
%       0 - original Denton method with the initial condition of x0=z0
%       1 - modified Denton method without the mentioned initial condition
% -------------------------------------------------------------------------
% OUTPUT:
% -> x = benchmarked or estimated high-frequency series 
% -> lambda = Lagrange multipliers (if needed)
% -------------------------------------------------------------------------
% REFERENCES: 1) Denton, F.T. (1971), Adjustments of monthly or quarterly
% series to annual totals: an approach based on quadratic minimization,
% Journal of the American Statistical Association, 66(333), pp. 99-102.
% 2) Dagum E.B. and Cholette, P.A. (2006), Benchmarking, Temporal
% Distribution, and Reconciliation Methods for Time Series, Lecture Notes
% in Statistics 186, Springer Science+Business Media, LLC, New York.
% 3) Temurshoev, U. (2012), Entropy-based benchmarking methods, GGDC 
% Research Memorandum GD-122, February.
% -------------------------------------------------------------------------
% Written by:   Dr. Umed Temurshoev, 07/10/2010 (updated in 2012)
%               Faculty of Economics and Business
%               University of Groningen, The Netherlands
%               E-mail: utemurshoev@gmail.com

n = length(z);      
m = length(y);      %number of years
k = n/m;            %length of each sub-annual data
B = kron(eye(m),ones(k,1));
r = y-B'*z;         %vector of totals discrepancies
if t == 0
    R = toeplitz([1;zeros(n-1,1)],ones(1,n));  %R'=inv(D)
    RR = R'*R;          %inv(A)=RR
    if a == 1           %additive level differences
        invA = eye(n);
    elseif a == 2       %proportional level differences
        invA = diag(z)*diag(z);
    elseif a == 3       %additive first differences
        invA = RR;
    elseif a == 4       %proportional first differences 
        invA = diag(z)*RR*diag(z);
    elseif a == 5       %additive second differences
        invA = R'*RR*R;
    else                %proportional second differences
        invA = diag(z)*R'*RR*R*diag(z);
    end 
%     invBAB = inv(B'*invA*B);
%     C = invA*B*invBAB;
    BAB = B'*invA*B;
    C = invA*B/BAB;
    x = z+C*r;
%     lambda = -invBAB*r;
    lambda = -BAB\r;
end

if t == 1
    D = toeplitz([-1 zeros(1,n-2)],[-1 1 zeros(1,n-2)]);
    D2 = toeplitz([1 zeros(1,n-3)],[1 -2 1 zeros(1,n-3)]);
    if a == 1           %additive level differences
        A = eye(n);
    elseif a == 2       %proportional level differences
        invz = 1./z;
        invz(z==0) = 1/eps;
        A = diag(invz)*diag(invz);
    elseif a == 3       %additive first differences
        A = D'*D;
    elseif a == 4       %proportional first differences
        DD = D'*D;
        invz = 1./z;
        invz(z==0) = 1/eps;
        A = diag(invz)*DD*diag(invz);
    elseif a == 5       %additive second differences
        A = D2'*D2;
    else                %proportional second differences
        DD = D2'*D2;
        invz = 1./z;
        invz(z==0) = 1/eps;
        A = diag(invz)*DD*diag(invz);
    end 
%     V = inv([A B;B' zeros(m,m)]);
%     s = V*[A zeros(n,m); B' eye(m)]*[z;r];
    V = [A B;B' zeros(m,m)];
    s = V\[A zeros(n,m); B' eye(m)]*[z;r];
    x = s(1:n);
    lambda = s(n+1:end);
end    
