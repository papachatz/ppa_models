clc;
clear all;
close all;

myCluster = parcluster('local');
myCluster.NumWorkers = 40;

n = 32;
q = 0.95;

mu_g_bit = 15.53789664;
mu_p_bit = mu_g_bit;
mu_g_1 = 44.22205226;
mu_g_2 = 41.32753379;
mu_g_3 = 55.55006159;
mu_g_4 = 53.1661757;
mu_g_5 = 52.0548568;
mu_xor = 19.39257639;
 
sigma_g_bit = 1.53899100508;
sigma_p_bit = sigma_g_bit;
sigma_g_1 = 4.16418282958;
sigma_g_2 = 4.03769658148;
sigma_g_3 = 5.00791093559;
sigma_g_4 = 4.72630994533;
sigma_g_5 = 4.66170001163;
sigma_xor = 1.58314851261;


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

%%
% mean value for prefix network

mu = [mu_g_1.*ones(1,16) mu_g_2.*ones(1,8) mu_g_3.*ones(1,4) mu_g_4.*ones(1,2) ...
    mu_g_5.*ones(1,1)]
% covariance matrix for prefix network
% no correlation between group generate nodes
Sigma = [(sigma_g_1^2).*ones(1,16) (sigma_g_2^2).*ones(1,8) (sigma_g_3^2).*ones(1,4) ...
    (sigma_g_4^2).*ones(1,2) (sigma_g_5^2).*ones(1,1)] .* eye(n-1)

% path mean delays and covariance matrix 
muprime = (A*mu')';
Sigmaprime = A*Sigma*(A');

% total delay 
mu_tot = [ kron(ones(1,n/2), [mu_g_bit mu_g_bit mu_p_bit]) ] + (A2*muprime')' + mu_xor

%total Sigma 
SS = [sigma_g_bit^2 0 0; 0 sigma_g_bit^2 0 ; 0 0 sigma_p_bit^2];
Sigma_tot = [ kron(eye(n/2), SS) ] + A2*Sigmaprime*(A2') + (sigma_xor^2)
%%
% evaluation of the cdf of the maximum delay 
% by a multivariate gaussian cdf
xrange = 250:0.1:350;
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

mypdf = [];
for i=6:5:length(xrange)    
   mypdf = [mypdf mycdf(i)-mycdf(i-1)]; 
end
figure();plot(xrange(6:5:end), mypdf);
mydata = [xrange(6:5:end); mypdf];
fileID = fopen('ks_16_pdf.txt','w');
fprintf(fileID,'%f %f\n',mydata);
fclose(fileID);

%%
fileID = fopen('KOGGE_STONE_32_VTHINTRA_hist.txt','r');
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
xrange = 200:0.1:270;
fileID = fopen('kogge_stone_16_model.txt', 'w');
mydata = [xrange; mycdf];
fprintf(fileID, '%f %f\n', mydata);
fclose(fileID);

fileID = fopen('kogge_stone_16_mc.txt', 'w');
mydata = [h1.XData; h1.YData];
fprintf(fileID, '%f %f\n', mydata);
fclose(fileID);
