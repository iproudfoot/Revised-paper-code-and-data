clear all
close all
save_file = 'C:\Users\isaac\Documents\neat_figs3\';
%% data
global t_data conc_data
t_data = [0 4 24];
conc_data = [150 117 86];

Kbar = 781;

%%
k1 = 72;
k2 = 0.11;
n = 4001;
n1 = 1;

mu = [190 0.075];
sigma = [25 0.05];
%%
[M,E] = gibbs_sampler_errfun(mu,sigma,n,n1,@med_err);

%%

meds = [median(M(:,1)),median(M(:,2))];
v1 = meds(1)*ones(1,10);
q1 = linspace(0,3000,10);
figure(1)
hist(M(:,1))
hold on
plot(v1,q1,'k--','linewidth',2)
xlabel('\lambda_{4}')
ylabel('counts')
file_name = [save_file 'eta7hist_normanorm2.eps'];
export_fig(figure(1),file_name)
%%
figure(2)
hist(M(:,2))
v2 = meds(2)*ones(1,10);
hold on
plot(v2,q1,'k--','linewidth',2)
xlabel('\lambda_{5}')
ylabel('counts')
file_name = [save_file 'eta8hist_normanorm2.eps'];
export_fig(figure(2),file_name)
%%
iters = 1:1:n;
figure(3)
plot(iters,E)
xlabel('iteration')
ylabel('error function, E')
file_name = [save_file 'med_err_trace_normanorm2.eps'];
export_fig(figure(3),file_name)
 %%
 figure(4)
 plot(iters,M(:,1))
 xlabel('iteration')
 ylabel('\lambda_{4}')
 file_name = [save_file 'eta7_trace_normanorm2.eps'];
export_fig(figure(4),file_name)
 %%
 figure(5)
 plot(iters,M(:,2))
 xlabel('iteration')
 ylabel('\lambda_{5}')
  file_name = [save_file 'lambda4hist.eps'];
export_fig(figure(5),file_name)
 %%
 w1 = linspace(0, max(M(:,1)),10);
 w2 = linspace(0, max(M(:,2)),10);
 figure(6)
 scatter(M(:,1),M(:,2))
 hold on
 plot(v1,w2,'k--','linewidth',2)
 hold on
 plot(w1,v2,'k--','linewidth',2)
 xlabel('\lambda_{4}')
 ylabel('\lambda_{5}')
  file_name = [save_file 'med_scatter.eps'];
export_fig(figure(6),file_name)
%  %%
%  n2 = 25;
%  V = zeros(2,n2);
%  params = [72 0.11];
%  l = 2;
%  UB = [300 2.5];
%  for k = 1:l
% V(k,:) = linspace(LB(k),UB(k),n2); % matrix of various param values
% end
% 
% [X,Y] = meshgrid(V(1,:),V(2,:));
%         U = zeros(n2,n2);
% for i = 1:l-1
%      for j = 1:n2
%             x = V(i,j);
%                 for m = 1:n2
%                     vec = params;
%                     vec(i) = x;
%                     vec(k) = V(k,m);
%                     U(m,j) =  nocellv2err(vec,t_data,conc_data)/1000;
%                 end
%       end
% end
% %%
%  figure(10)
%         contourf(X,Y,U,[0 0.25 0.5 1.25 2.5 5 10])
%         hold on
%         plot(meds(1),meds(2),'rx','Markersize',10)
%         hold on
%         plot(72,0.11,'wx','Markersize',10)
%         colorbar
% %         xticks([0.5*initial_guess1(i) initial_guess1(i) 1.5*initial_guess1(i)])
% %         yticks([0.5*initial_guess1(k) initial_guess1(k) 1.5*initial_guess1(k)])
%         xlabel('\lambda_{4}','fontsize',18)
%         ylabel('\lambda_{5}','fontsize',18)
%         file_name2 = [save_file 'medparamerrfull.eps'];
%         export_fig(figure(10),file_name2)
 %%
 figure(11)
 q1 = linspace(0,0.35,10);
histogram(M(:,1),'Normalization','probability')
hold on
plot(v1,q1,'k--','linewidth',2)
hold on
x = 0:0.1:300;
m1 = mu(1);
s1 = sigma(1);
plot(x,7.5*gauss_distribution(x,m1,s1),'r','linewidth',2)
xlabel('\lambda_{4}')
ylabel('counts')
file_name = [save_file 'eta7scalehist_normanorm2.eps'];
export_fig(figure(11),file_name)
%%
figure(12)
histogram(M(:,2),'Normalization','probability')
hold on
q1 = linspace(0,0.4,10);
v2 = meds(2)*ones(1,10);
plot(v2,q1,'k--','linewidth',2)
hold on
x2 = 0:0.01:1;
a = 2;
b = 10;
m2 = mu(2);
s2 = sigma(2);
%plot(x2,0.3*beta_dens(x2,a,b)/5,'r')
plot(x2,0.01*gauss_distribution(x2,m2,s2),'r','linewidth',2)
xlabel('\lambda_{5}')
ylabel('counts')
file_name = [save_file 'lambda5hist.eps'];
export_fig(figure(12),file_name)
 %%
 function f = med_err(x)
 global t_data conc_data
 f = nocellv2err(x,t_data,conc_data)/1000;
 end