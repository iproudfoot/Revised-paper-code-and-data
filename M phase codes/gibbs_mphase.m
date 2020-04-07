clear all
close all

%% Fitting to Adrian's Data
mydir = 'C:\Users\isaac\Documents\MATLAB\lab_data\'; % Specify directory
save_file = 'C:\Users\isaac\Documents\neat_figs3\'; % Place to save figures


%%
V_onecell=2.425*10^(-6);
B0=6.93*1000;
T0=100*B0/33;
Bmaxm = 3940;
Kbar = 781;
Tmax = 4*B0;

%%
%% Plotting error function landscapes of parameters pairwise

% str3 = ["$$ \kappa_{1} $$", "$$ \kappa_{2} $$", "$$ \kappa_{3} $$", "$$ \kappa_{4} $$", "$$ a_{2} $$", "$$ a_{3} $$", "$$ \beta_{1} $$", "$$ \beta_{2} $$", "$$ \beta_{3} $$", "$$ \gamma_{2} $$", "$$ \gamma_{3} $$", "$$ \alpha $$"]; 
% str3 = ["$$ \eta_{15} $$", "$$ \eta_{2} $$", "$$ \eta_{16} $$", "$$ \eta_{4} $$", "$$\eta_{10}$$", "$$ \eta_{11} $$", "$$ \eta_{12}$$", "$$ \eta_{17} $$", "$$ \eta_{14} $$", "$$ \eta_{7}$$", "$$\eta_{6}$$", "$$\eta_{5}$$"];
% str4 = str3(~~logic_vec);

%     str3 = str(m); 
%     str5 = sprintf('%g',m);


%%
LB = [0 0];
UB = [0.001 0.01];
scaling_factors = [B0 B0];
n = 1050;
n1 = 50;
tic
[M1,E1,L1] = gibbs_sampler_uniform_mphase(LB,UB,n,n1,@err_fun);
toc
%%
str4 = ["\eta_{18}", "\eta_{19}", "\eta_{20}"];
num_params = 2;
meds = zeros(1,num_params);
for k = 1:num_params
    meds(k) = scaling_factors(k).*median(M1(:,k));
end
for i = 1:num_params
    str1 = sprintf('%g',i);
    for j = i:num_params
        str2 = sprintf('%g',j);
        if j == i
            v = linspace(0,UB(i),10);
            q = meds(i)*ones(1,10);
            figure(10*i+j)
            box on
            hist(scaling_factors(i).*M1(:,i))
            xlabel(str4(i))
            ylabel('counts')
            hold on
            plot(q,v,'k--','linewidth',2)
            file_name = [save_file 'unifhist3' str1 '.eps'];
            export_fig(figure(10*i+j),file_name)
        else
            figure(10*i+j)
            scatter(scaling_factors(i).*M1(:,i),scaling_factors(j).*M1(:,j),'rx')
            hold on
            xticks = get(gca,'Xtick');
            yticks = get(gca,'Ytick');
            v1 = linspace(min(xticks),max(xticks),10);
            q1 = meds(i)*ones(1,10);
            v2 = linspace(min(yticks),max(yticks),10);
            q2 = meds(j)*ones(1,10);
            plot(v1,q2,'k--','linewidth',2)
            hold on
            plot(q1,v2,'k--','linewidth',2)
            box on
            xlabel(str4(i))
            ylabel(str4(j))
            file_name = [save_file 'unifscatter3' str1 str2 '.eps'];
            export_fig(figure(10*i+j),file_name)
        end
    end
end

%%
% iters = 1:1:length(E1);
% 
% figure(50)
% box on
% plot(iters,E1)
% xlabel('iteration')
% ylabel('Error function (E)')
% file_name = [save_file 'uniferr_trace.eps'];
% export_fig(figure(50),file_name)
% 
% %%
% for i = 1:4
%     str3 = sprintf('%g',i);
%     figure(50+i)
%     box on
%     plot(iters, scaling_factors(i).*M1(:,i))
%     xlabel('iteration')
%     ylabel(str4(i))
%     file_name = [save_file 'uniftrace' str3 '.eps'];
%     export_fig(figure(50+i),file_name)
% end
% 
% %%
% mu1 = [202 78.3 10.6 410];
% sigma1 = [60 8 1.1 80];
% % for i = 1:4
% %     str1 = sprintf('%g',i);
% %     V =  linspace(min(scaling_factors(i)*M1(:,i)),max(scaling_factors(i)*M1(:,i)),1001);
% %     figure(60+i)
% %     histogram(scaling_factors(i)*M1(:,i),'Normalization','probability')
% %     hold on
% %     plot(V,scaling_factors(i)*gauss_distribution(V,mu1(i),sigma1(i)),'r','linewidth',2)
% %     xlabel(str4(i))
% %     ylabel('normalized counts')
% %     file_name = [save_file 'unifnormhist' str1 '.eps'];
% %     export_fig(figure(60+i),file_name)
% % end
% %%
% str1 = sprintf('%g',1);
%     V =  linspace(min(scaling_factors(1)*M1(:,1)),max(scaling_factors(1)*M1(:,1)),1001);
%     figure(61)
%     box on
%     histogram(scaling_factors(1)*M1(:,1),'Normalization','probability')
%     hold on
%     plot(V,30*gauss_distribution(V,mu1(1),sigma1(1)),'r','linewidth',2)
%     hold on
%     yticks = get(gca,'Ytick');
%     v = linspace(0,max(yticks),10);
%     q = meds(1)*ones(1,10);
%     plot(q,v,'k--','linewidth',2)
%     xlabel(str4(1))
%     ylabel('probability density')
%     file_name = [save_file 'unifnormhist' str1 '.eps'];
%     export_fig(figure(61),file_name)
%  %%   
%     
%     str1 = sprintf('%g',2);
%     V =  linspace(min(scaling_factors(2)*M1(:,2)),max(scaling_factors(2)*M1(:,2)),1001);
%     figure(62)
%     box on
%     histogram(scaling_factors(2)*M1(:,2),'Normalization','probability')
%     hold on
%     plot(V,2.5*gauss_distribution(V,mu1(2),sigma1(2)),'r','linewidth',2)
%     hold on
%     yticks = get(gca,'Ytick');
%     v = linspace(0,max(yticks),10);
%     q = meds(2)*ones(1,10);
%     plot(q,v,'k--','linewidth',2)
%     xlabel(str4(2))
%     ylabel('probability density')
%     file_name = [save_file 'unifnormhist' str1 '.eps'];
%     export_fig(figure(62),file_name)
%     
%     %%
%     str1 = sprintf('%g',3);
%     V =  linspace(min(scaling_factors(3)*M1(:,3)),max(scaling_factors(3)*M1(:,3)),1001);
%     figure(63)
%     box on
%     histogram(scaling_factors(3)*M1(:,3),'Normalization','probability')
%     hold on
%     plot(V,0.7*gauss_distribution(V,mu1(3),sigma1(3)),'r','linewidth',2)
%     hold on
%     yticks = get(gca,'Ytick');
%     v = linspace(0,0.3,10);
%     q = meds(3)*ones(1,10);
%     plot(q,v,'k--','linewidth',2)
%     xlabel(str4(3))
%     ylabel('probability density')
%     file_name = [save_file 'unifnormhist' str1 '.eps'];
%     export_fig(figure(63),file_name)
%     
%     %%
%     str1 = sprintf('%g',4);
%     V =  linspace(min(scaling_factors(4)*M1(:,4)),max(scaling_factors(4)*M1(:,4)),1001);
%     figure(64)
%     box on
%     histogram(scaling_factors(4)*M1(:,4),'Normalization','probability')
%     hold on
%     plot(V,50*gauss_distribution(V,mu1(4),sigma1(4)),'r','linewidth',2)
%     hold on
%     yticks = get(gca,'Ytick');
%     v = linspace(0,max(yticks),10);
%     q = meds(4)*ones(1,10);
%     plot(q,v,'k--','linewidth',2)
%     xlabel(str4(4))
%     ylabel('probability density')
%     file_name = [save_file 'unifnormhist' str1 '.eps'];
%     export_fig(figure(64),file_name)
% %%
% figure(65)
% box on
% plot(iters,L1)
% xlabel('iteration')
% ylabel('Marginal Likelihood (L)')

%% Credible intervals
for k = 1:4
M2 = zeros(n-n1,4);
M2 = sort(M1(:,k));
LB1 = M2(26)*scaling_factors(k)
UB1 = M2(975)*scaling_factors(k)
end

%%
function f = err_fun(x)
    f = mphase_err(x);
end