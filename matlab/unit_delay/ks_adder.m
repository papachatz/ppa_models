function [dmax, Gmax] = ks_adder(mu, n, fanout_coeff, sigma_coeff)

%dmax: maximum adder delay
%Gmax: maximum delay that end up to G_{n:0} node

stages = log2(n);
dmax = 0;
sigma = mu *sigma_coeff;
% node delays of group propagate nodes
mu_g = mu *(1+fanout_coeff*1);
sigma_g = mu_g *sigma_coeff;
G = normrnd(mu_g, sigma_g, stages-1, n);
mu_g = mu *(1+fanout_coeff*0);
sigma_g = mu_g *sigma_coeff;
G = [G; normrnd(mu_g, sigma_g, 1, n)];
% path delays for group propagate nodes
maxdel = zeros( stages, n);

for st=1:stages
    for i=(2^(st-1))+1:n
    
        if(st==1)
            pg = max( max(normrnd(mu, sigma), normrnd(mu, sigma)), normrnd(mu, sigma));
            maxdel(st, i) = pg + G(st,i); 
        else
            maxdel(st, i) = max(maxdel(st-1,i), maxdel(st-1, i-2^(st-1))) + G(st,i);
        end
        
        if(st==stages)
            maxdel(st, i) = maxdel(st, i) + normrnd(mu, sigma);
        end
        
        dmax = max(dmax, maxdel(st, i));
        
    end    
end

Gmax = maxdel(stages, n);

end

