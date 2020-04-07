function [M,E] = gibbs_sampler_errfun(mu,sigma,n,n1,err_fun)
%gibbs_sampler_params samples the parameter vector x with initial parameter
%vector x_0. Output M is the n x length(mu) matrix of samples. Likelihood
%function is taken to be -ve exponential of err_fun.
%mu is the vector of means for the prior normal distributions.
%sigma is the vector of std devs for the priors.
%n is the total number of iterations.
%n1 is the number of iterations discarded for 'burn in'.
%likelihood_fun is the likelihood function used for the model.

%% Generate initial param vector
l = length(mu);     % number of parameters
M = zeros(n,l);
E = zeros(n,1);
s = 0.01;
x_0 = mu + randn([1,l]).*sigma; % sample initial values

for i = 1:n     % loop for total samples
    if i <= n1
        if i == 1
            x = x_0; % set comparison vec to be initial at start
        end
        x_new = x;
        for j = 1:l     %loop over condtional for each parameter
            x_new(j) = mu(j) + randn(1)*sigma(j); %sample from prior
            marg_like_old = exp(-(x(j)-mu(j))^2/(2*sigma(j)^2))*exp(-err_fun(x)/s)/(sqrt(2*pi*sigma(j)^2));
            marg_like_new = exp(-(x_new(j)-mu(j))^2/(2*sigma(j)^2))*exp(-err_fun(x_new)/s)/(sqrt(2*pi*sigma(j)^2)); %compute marginal likelihood
            r = marg_like_new/marg_like_old; % ratio of marginal likelihoods
            u = rand(1);
            if r >= u % acceptance criteria
                x = x_new;
            end
        end
    end
    for j = 1:l     %loop over condtional for each parameter
            x_new(j) = mu(j) + randn(1)*sigma(j); %sample from prior
            marg_like_old = exp(-(x(j)-mu(j))^2/(2*sigma(j)^2))*exp(-err_fun(x)/s)/(sqrt(2*pi*sigma(j)^2));
            marg_like_new = exp(-(x_new(j)-mu(j))^2/(2*sigma(j)^2))*exp(-err_fun(x_new)/s)/(sqrt(2*pi*sigma(j)^2)); %compute marginal likelihood
            r = marg_like_new/marg_like_old; % ratio of marginal likelihoods
            u = rand(1);
            if r >= u % acceptance criteria
                x = x_new;
            end
            M(i,:) = x;
            E(i) = err_fun(x);
    end
end

