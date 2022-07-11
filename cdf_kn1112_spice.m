clc;
clear all;
close all;

myCluster = parcluster('local');
myCluster.NumWorkers = 10;

n = 16;
q = 0.95;

mu_g_bit = 18.3811;
mu_g_1 = 51.8874;
mu_g_2 = 54.5155;
mu_g_3 = 52.8414;
mu_g_4 = 55.5801;
mu_xor = 31.5350;

sigma_g_bit = 1.7273;
sigma_g_1 = 4.8546;
sigma_g_2 = 4.7197;
sigma_g_3 = 4.6617;
sigma_g_4 = 5.5872;
sigma_xor = 2.7220;


fprintf("Number of bits: %d\n",n);

% matrix used in the proposed transformation
A = [];
for i=0:log2(n)-1
    A = [A kron(eye(2^(log2(n)-1-i)), ones(2^i,1)) ];
end

fprintf("Transformation matrix:\n");
A 

fprintf("Transformation matrix 2");
A2 = [ kron(eye(n/2), ones(3,1)) ]

% mean value for prefix network
mu = [mu_g_1 mu_g_1 mu_g_1 mu_g_1 mu_g_1 mu_g_1 mu_g_1 mu_g_1 ...
    mu_g_2 mu_g_2 mu_g_2 mu_g_2 mu_g_3 mu_g_3 mu_g_4 ];

% covariance matrix for prefix network
% no correlation between group generate nodes
Sigma = [sigma_g_1^2 sigma_g_1^2 sigma_g_1^2 sigma_g_1^2 sigma_g_1^2 ...
    sigma_g_1^2 sigma_g_1^2 sigma_g_1^2 sigma_g_2^2 sigma_g_2^2 sigma_g_2^2 ...
    sigma_g_2^2 sigma_g_3^2 sigma_g_3^2 sigma_g_4^2 ] .* eye(n-1);

% path mean delays and covariance matrix 
muprime = (A*mu')';
Sigmaprime = A*Sigma*(A');

% total delay 
mu_tot = mu_g_bit + (A2*muprime')' + mu_xor

%total Sigma 
Sigma_tot = (sigma_g_bit^2)*eye(3*(n/2)) + A2*Sigmaprime*(A2') + (sigma_xor^2)

%%
% evaluation of the cdf of the maximum delay 
% from a multivariate gaussian cdf
xrange = 230:0.1:320;
mycdf = zeros(1,length(xrange));
qflag = 1;

parfor i = 1:length(xrange)
   X = xrange(i) * ones(1, 3*(n/2));
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
fileID = fopen('KNOWLES1112_VTHINTRA_hist.txt','r');
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

fileID = fopen('knowles_1112_16_model.txt', 'w');
mydata = [xrange; mycdf];
fprintf(fileID, '%f %f\n', mydata);
fclose(fileID);

fileID = fopen('knowles_1112_16_mc.txt', 'w');
mydata = [h1.XData; h1.YData];
fprintf(fileID, '%f %f\n', mydata);
fclose(fileID);
