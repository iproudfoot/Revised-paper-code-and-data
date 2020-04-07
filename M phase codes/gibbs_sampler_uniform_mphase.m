function [M,E,L] = gibbs_sampler_uniform_mphase(LB,UB,n,n1,err_fun)
%gibbs_sampler_params samples the parameter vector x with initial parameter
%vector x_0. Output M is the n x length(mu) matrix of samples. Likelihood
%function is taken to be -ve exponential of err_fun.
%LB is the vector of lower bounds for the uniform priors.
%UB is the vector of upper bounds for the uniform priors.
%n is the total number of iterations.
%n1 is the number of iterations discarded for 'burn in'.
%likelihood_fun is the likelihood function used for the model.

%% Generate initial param vector
a = 100;

l = length(LB);     % number of parameters
M = zeros(n-n1,l);
E = zeros(n-n1,1);
L = zeros(n-n1,1);
x_0 = zeros(1,l);
x_0(1:l) = LB + (UB-LB).*rand([1,l]); % sample initial values
% x_0(5) = 95 + randn([1,1])*25;
% x_0(6) = 0.15 + randn([1, 1])*0.05;
for i = 1:n     % loop for total samples
    if i <= n1
        if i == 1
            x = x_0; % set comparison vec to be initial at start
        end
        x_new = x;
        for j = 1:l     %loop over condtional for each parameter
            x_new(j) = LB(j) + rand(1)*(UB(j)-LB(j)); %sample from prior
            marg_like_old = exp(-err_fun(x)/a)/(UB(j)-LB(j));
            marg_like_new = exp(-err_fun(x_new)/a)/(UB(j)-LB(j)); %compute marginal likelihood
            r = (marg_like_new/marg_like_old); % ratio of marginal likelihoods
            u = rand(1);
            if r >= u % acceptance criteria
                x = x_new;
            end
        end
    else
        for j = 1:l     %loop over condtional for each parameter
                x_new(j) = LB(j) + rand(1)*(UB(j)-LB(j)); %sample from prior
                marg_like_old = exp(-err_fun(x)/a)/(UB(j)-LB(j));
                marg_like_new = exp(-err_fun(x_new)/a)/(UB(j)-LB(j)); %compute marginal likelihood
                r = (marg_like_new/marg_like_old); % ratio of marginal likelihoods
                u = rand(1);
                if r >= u % acceptance criteria
                    x = x_new;
                end
                M(i-n1,:) = x;
                E(i-n1) = err_fun(x);
                L(i-n1) = exp(-err_fun(x)/a)/(UB(j)-LB(j));
        end
    end
%     x(5) = 95 + randn([1,1])*25;   % sample next med parameters at this point
%     x(6) = 0.15 + randn([1, 1])*0.05;
end