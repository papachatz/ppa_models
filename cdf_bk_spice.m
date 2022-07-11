clc;
clear all;
close all;

myCluster = parcluster('local');
myCluster.NumWorkers = 45;

q = 0.9987;

mu_g_bit = 16.96338323;
mu_p_bit =  9.132502113;
mu_g_1 = 49.62524715;
mu_g_2 = 55.55737201;
mu_g_3 = 57.38386092;
mu_g_4 = 56.83453404;
mu_g_5 = 53.35902084;
mu_g_6 = 49.92165706;
mu_xor = 31.52140507;

sigma_g_bit = 1.69250559906;
sigma_p_bit = 0.729029607527;
sigma_g_1 = 4.60590297822;
sigma_g_2 = 4.90005347143;
sigma_g_3 = 5.09278683566;
sigma_g_4 = 5.07956974211;
sigma_g_5 = 4.53129506003;
sigma_g_6 = 4.49238034182;
sigma_xor = 2.67146025349;

mu_p_bk = [mu_g_1 mu_g_1 mu_g_1 mu_g_1 mu_g_1 mu_g_1 ...
    mu_g_2 mu_g_2 mu_g_2 mu_g_3 mu_g_4+mu_g_5+mu_g_6 ]';

Sigma_p_bk = [sigma_g_1^2 sigma_g_1^2 sigma_g_1^2 sigma_g_1^2 sigma_g_1^2 ...
    sigma_g_1^2 sigma_g_2^2 sigma_g_2^2 sigma_g_2^2 sigma_g_3^2 ...
    sigma_g_4^2+sigma_g_5^2+sigma_g_6^2] .* eye(11)

A_bk = [1 0 0 0 0 0 1 0 0 0 1; ...
        0 1 0 0 0 0 1 0 0 0 1; ...
        0 0 1 0 0 0 0 1 0 1 1; ...
        0 0 0 1 0 0 0 1 0 1 1; ...
        0 0 0 0 1 0 0 0 1 1 1; ...
        0 0 0 0 0 1 0 0 1 1 1];

mu_bk = A_bk*mu_p_bk
Sigma_bk = A_bk * Sigma_p_bk * A_bk'

n = 16;
A2 = kron(eye(6),ones(3,1))
%mu_tot = mu_g_bit + A2 * mu_bk + mu_xor
mu_tot = [ kron(ones(1,n/2-2), [mu_g_bit mu_g_bit mu_p_bit]) ]' + A2*mu_bk + mu_xor

%Sigma_tot = (sigma_g_bit^2)*eye(3*6) ...
%    +  A2 * Sigma_bk * (A2') ...
%    + (sigma_xor^2)
Sigma_tot = [ kron(eye(n/2-2), [sigma_g_bit^2 0 0; 0 sigma_g_bit^2 0; 0 0 sigma_p_bit^2]) ] ...
    + A2*Sigma_bk*(A2') + (sigma_xor^2)

% evaluation of the cdf of the maximum delay 
% from a multivariate gaussian cdf
xrange = 330:0.1:420;
mycdf = zeros(1,length(xrange));
qflag = 1;

parfor i = 1:length(xrange)
   X = xrange(i) * ones(1, 3*6);
   mycdf(i)  = mvncdf(X, mu_tot', Sigma_tot);
end

for i=1:length(xrange)
    if(qflag)
       if(mycdf(i)>=q)
           q1 = xrange(i);
           qflag = 0;
           fprintf("%f quantile: %f\n", q, xrange(i));
       end
    end
end

figure();
plot(xrange, mycdf);

%%
fileID = fopen('BRENT_KUNG_VTHINTRA_hist.txt','r');
A=fscanf(fileID, '%f %f', [2 Inf]);
fclose(fileID);

A = A *1e12;
hold on;
h1 = cdfplot(A(1,:));
legend("CDF of max (proposed model)", "CDF of max (M.C.)");

%estimate 0.95 Quantile of M.C. data
qflag = 1;

for i=1:length(h1.XData)
   if(qflag)
      if(h1.YData(i)>=q)
        q2 = h1.XData(i);
        fprintf("%f quantile of M.C. data: %f\n", q, h1.XData(i));
        qflag = 0;
      end
   end
end

fprintf("Approximation error at %f quantile: %f (perc) \n", q, 100.0*(q2-q1)/q2);


%%
fileID = fopen('brent_kung_16_model.txt', 'w');
mydata = [xrange; mycdf];
fprintf(fileID, '%f %f\n', mydata);
fclose(fileID);

fileID = fopen('brent_kung_16_mc.txt', 'w');
mydata = [h1.XData; h1.YData];
fprintf(fileID, '%f %f\n', mydata);
fclose(fileID);

