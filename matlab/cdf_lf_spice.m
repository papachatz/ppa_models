clc;
clear all;
close all;

myCluster = parcluster('local');
myCluster.NumWorkers = 10;

n = 16;
q = 0.95;

mu_g_bit = 16.9836353507;
mu_p_bit = 8.408699236;
mu_g_1 = 49.2114329058;
mu_g_2 = 55.2334153707;
mu_g_3 = 59.0434897796;
mu_g_4 = 56.8601707816;
mu_g_5 = 49.8054195792;
mu_xor = 31.5312120842;

sigma_g_bit = 1.55495782602;
sigma_p_bit = 0.708300844205;
sigma_g_1 = 4.48600132244
sigma_g_2 = 4.84948070863;
sigma_g_3 = 5.11619751782;
sigma_g_4 = 5.19056272334;
sigma_g_5 = 4.49846496739;
sigma_xor = 2.7559434048;

mu_p1 = [mu_g_1 mu_g_1 mu_g_1 mu_g_1 mu_g_1 mu_g_1 mu_g_1];
Sigma_p1 = [sigma_g_1^2 sigma_g_1^2 sigma_g_1^2 sigma_g_1^2 sigma_g_1^2 ...
    sigma_g_1^2 sigma_g_1^2] .*eye(n/2 -1)

mu_glf = [0 mu_g_2 mu_g_2 mu_g_2 mu_g_3 mu_g_3 mu_g_4+mu_g_5];
Sigma_glf = [0 sigma_g_2^2 sigma_g_2^2 sigma_g_2^2 sigma_g_3^2 ...
    sigma_g_3^2 sigma_g_4^2+sigma_g_5^2] .*eye(n/2 -1)


A_lf = [ 1 0 0 0; 0 1 0 0; 0 1 0 0; 0 0 1 0; 0 0 1 0; 0 0 0 1; 0 0 0 1]

A = [];
for i=0:log2(n/2)-1
    A = [A kron(eye(2^(log2(n/2)-1-i)), ones(2^i,1)) ];
end

mu_lf = mu_p1' + (A_lf * A) * mu_glf'
Sigma_lf = Sigma_p1 + (A_lf * A) * Sigma_glf * (A_lf * A)' 

A2 = [ kron(eye(n/2 - 1), ones(3,1))];

mu_tot = [ kron(ones(1,n/2-1), [mu_g_bit mu_g_bit mu_p_bit]) ] + (A2*mu_lf)' + mu_xor
Sigma_tot = [ kron(eye(n/2-1), [sigma_g_bit^2 0 0; 0 sigma_g_bit^2 0; 0 0 sigma_p_bit^2]) ] ...
    + A2*Sigma_lf*(A2') + (sigma_xor^2)
%%

% evaluation of the cdf of the maximum delay 
% bya multivariate gaussian cdf
xrange = 275:0.1:380;
mycdf = zeros(1,length(xrange));
qflag = 1;

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
fileID = fopen('LADNER_FISCHER_VTHINTRA_hist.txt','r');
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
fileID = fopen('ladner_fischer_16_model.txt', 'w');
mydata = [xrange; mycdf];
fprintf(fileID, '%f %f\n', mydata);
fclose(fileID);

fileID = fopen('ladner_fischer_16_mc.txt', 'w');
mydata = [h1.XData; h1.YData];
fprintf(fileID, '%f %f\n', mydata);
fclose(fileID);
