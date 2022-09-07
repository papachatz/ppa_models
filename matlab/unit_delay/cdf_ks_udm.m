clc;
clear all;
close all;

myCluster = parcluster('local');
myCluster.NumWorkers = 45;

n = 16;     % bit length
q = 0.95;   % quantile point
trees = 1;  % l parameter 
paths = 3;  % prefix paths are extended to #paths adder paths

fanout_coeff = 0.0227;
unit_del = 1;
sigma_coeff = 0.05;

mu_g_bit = unit_del*(1+fanout_coeff*0);
mu_g_1 = unit_del*(1+fanout_coeff*1);
mu_g_2 = unit_del*(1+fanout_coeff*1);
mu_g_3 = unit_del*(1+fanout_coeff*1);
mu_g_4 = unit_del*(1+fanout_coeff*0);
mu_xor = unit_del*(1+fanout_coeff*0);                  

sigma_g_bit = sigma_coeff * mu_g_bit;
sigma_g_1 = sigma_coeff * mu_g_1;
sigma_g_2 = sigma_coeff * mu_g_2;
sigma_g_3 = sigma_coeff * mu_g_3;
sigma_g_4 = sigma_coeff * mu_g_4;
sigma_xor = sigma_coeff * mu_xor;

qflag = 1;
fprintf("Number of bits: %d\n",n);

% matrix used in the proposed transformation
A = [];
for i=0:log2(n)-1
    A = [A kron(eye(2^(log2(n)-1-i)), ones(2^i,1)) ];
end

fprintf("Transformation matrix:\n");
A 

fprintf("Transformation matrix 2");
A2 = [ kron(eye(n/2), ones(paths,1)) ]

% mean value for prefix network
mu = [ones(n-2,1)*mu_g_1; mu_g_4]';

% covariance matrix for prefix network
% no correlation between group generate nodes
Sigma = ((sigma_coeff*mu).^2).* eye(n-1)

% path mean delays and covariance matrix 
muprime = (A*mu')'
Sigmaprime = A*Sigma*(A')

% total delay 
mu_tot = mu_g_bit + (A2*muprime')' + mu_xor
mu_tot = kron(ones(trees,1), mu_tot')'

%total Sigma 
Sigma_tot = kron(eye(trees),(sigma_g_bit^2)*eye(3*(n/2)) ...
    + A2*Sigmaprime*(A2') + (sigma_xor^2))

%%
% evaluation of the cdf of the maximum delay 
% by a multivariate gaussian cdf
xrange = log2(n)+2:0.005:log2(n)+3;
mycdf = zeros(1,length(xrange));
gflag = 1;

parfor i = 1:length(xrange)
   X = xrange(i) * ones(1, paths*trees*n/2);
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

mypdf = []; % construct pdf from consecutive cdf points
for i=2:length(mycdf)    
   mypdf = [mypdf mycdf(i)-mycdf(i-1)]; 
end
figure();
plot(xrange(2:end), mypdf);
legend("PDF of max delay (proposed model)");

mydata = [xrange(2:end); mypdf];
fileID = fopen(sprintf('ks_%d_pdf_%d.txt',n,trees),'w');
fprintf(fileID,'%f %f\n',mydata);
fclose(fileID);

%%
% monte-carlo simulations

mu = 1;
sigma = 0.05;
mcruns = 100000;
dmax = zeros(1, mcruns);
Gmax = zeros(1, mcruns);

parfor mcrun=1:mcruns
    [dmax(mcrun), Gmax(mcrun)] = ks_adder(mu, n, fanout_coeff, sigma_coeff);
end

mudmax = mean(dmax);
sigmadmax = std(dmax);
fprintf("Mean of max delay (M.C.): %f\n", mudmax);
fprintf("Sigma of max delay (M.C.): %f\n", sigmadmax);

mugmax = mean(Gmax);
sigmagmax = std(Gmax);
fprintf("Mean of G%d:0 delay (M.C.): %f\n", n-1, mugmax);
fprintf("Sigma of G%d:0 delay (M.C.): %f\n", n-1, sigmagmax);

figure();
plot(xrange, mycdf,'-o');
hold on;
h1 = cdfplot(dmax);
hold on;
h2 = cdfplot(Gmax);
legend("CDF of max (proposed model)", ...
    "CDF of max (M.C.)", "CDF of G15:0 (M.C.)");

%estimate 0.95 Quantile of M.C. data
xrange = h1.XData;
yrange = h1.YData;
qflag = 1;

for i=1:length(xrange)
   if(qflag)
      if(yrange(i)>=q)
        q2 = xrange(i);
        fprintf("%f quantile of M.C. data: %f\n", q, xrange(i));
       qflag = 0;
      end
   end
end

fprintf("Approximation error at %f quantile: %f (perc) \n", ...
    q, 100.0*(q2-q1)/q2);

%%
% write data to file 
step = 100;

mydata = [h1.XData; h1.YData];
mydata = [mydata(1,[1:step:length(mydata(1,:))]); ...
    mydata(2,[1:step:length(mydata(1,:))])];
fileID = fopen(sprintf('ks_%d_mc_dmax.txt',n),'w');
fprintf(fileID,'%f %f\n', mydata);
fclose(fileID);

mydata = [h2.XData; h2.YData];
mydata = [mydata(1,[1:step:length(mydata(1,:))]); ...
    mydata(2,[1:step:length(mydata(1,:))])];
fileID = fopen(sprintf('ks_%d_mc_gmax.txt',n),'w');
fprintf(fileID,'%f %f\n', mydata);
fclose(fileID);
