clc;
clear all;
close all;

myCluster = parcluster('local');
myCluster.NumWorkers = 10;

n = 16;
q = 0.9987;

mu_g_bit = 15.34840213;
mu_p_bit =  9.132502113;
mu_g_1 = 44.62038555;
mu_g_2 = 47.17094912;
mu_g_3 = 45.2834975;
mu_g_4 = 42.09194804;
mu_g_5 = 50.12835239;
mu_xor = 31.36564507;

sigma_g_bit = 1.51080045527;
sigma_p_bit = 0.729029607527;
sigma_g_1 = 4.09328649379;
sigma_g_2 = 4.3182243412;
sigma_g_3 = 4.20959992913;
sigma_g_4 = 3.93000985737;
sigma_g_5 = 4.65913586135;
sigma_xor = 2.707089662;

Ap1 = [1 0 0 0 0 0 1; 0 1 0 0 0 0 1; 0 0 1 0 0 0 1; 0 0 0 1 0 0 1; ...
    0 0 0 0 1 0 1; 0 0 0 0 0 1 1; 0 0 0 0 0 0 1]
mu_p1 = Ap1 *[mu_g_1 mu_g_1 mu_g_1 mu_g_1 mu_g_1 mu_g_1 mu_g_5]'
Sigma_p1 = Ap1 *([sigma_g_1^2 sigma_g_1^2 sigma_g_1^2 sigma_g_1^2 ...
    sigma_g_1^2 sigma_g_1^2 sigma_g_5^2] .*eye(7)) * Ap1';

mu_ghc = [mu_g_2 mu_g_2 mu_g_2 mu_g_2 mu_g_3 mu_g_3 mu_g_4];
Sigma_ghc = [sigma_g_2^2 sigma_g_2^2 sigma_g_2^2 sigma_g_2^2 ...
    sigma_g_3^2 sigma_g_3^2 sigma_g_4^2] .*eye(7);


A_hc = [ 1 0 0 0; 1 0 0 0; 0 1 0 0; 0 1 0 0; 0 0 1 0; 0 0 1 0; 0 0 0 1]
A_n2 = [1 0 0 0 1 0 1; 0 1 0 0 1 0 1; 0 0 1 0 0 1 1; 0 0 0 1 0 1 1]

mu_hc = mu_p1 + (A_hc * A_n2) * mu_ghc'
Sigma_hc = Sigma_p1 + (A_hc * A_n2) * Sigma_ghc * (A_hc * A_n2)' 


A2 = [ kron(eye(7), ones(3,1))];

mu_tot = [ kron(ones(1,n/2-1), [mu_g_bit mu_g_bit mu_p_bit]) ] + (A2*mu_hc)' + mu_xor
Sigma_tot = [ kron(eye(n/2-1), [sigma_g_bit^2 0 0; 0 sigma_g_bit^2 0; 0 0 sigma_p_bit^2]) ] ...
    + A2*Sigma_hc*(A2') + (sigma_xor^2)


% evaluation of the cdf of the maximum delay 
% from a multivariate gaussian cdf
qflag = 1;
mycdf = [];
xrange = 240:0.1:320;

parfor i = 1:length(xrange)
   X = xrange(i) * ones(1, 3*7);
   mycdf(i)  = mvncdf(X, mu_tot, Sigma_tot);
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
fileID = fopen('HAN_CARLSON_VTHINTRA_hist.txt','r');
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
fileID = fopen('han_carlson_16_model.txt', 'w');
mydata = [xrange; mycdf];
fprintf(fileID, '%f %f\n', mydata);
fclose(fileID);

fileID = fopen('han_carlson_16_mc.txt', 'w');
mydata = [h1.XData; h1.YData];
fprintf(fileID, '%f %f\n', mydata);
fclose(fileID);

