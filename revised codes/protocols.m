clear all
close all

save_file = 'C:\Users\isaac\Documents\neat_figs2\';

t7 = -1:0.1:0;
x7 = ones(1,length(t7));


t1 = 0:0.1:24;
t4 = -1:0.1:24;
t2 = -2:0.1:24;
t3 = -3:0.1:24;
t5 = -2:0.1:0;
t6 = 0:0.01:2;
x = 2*ones(1,length(t5));
x1 = ones(1,length(t1));
x2 = 2*ones(1,length(t1));
x3 = 3*ones(1,length(t1));
y1 = ones(1,length(t2));
y2 = 2*ones(1,length(t2));
z1 = ones(1,length(t3));
z2 = 2*ones(1,length(t3));
a1 = ones(1,length(t4));
y = 2*ones(1,length(t6));
%% Binding in medium/uptake
figure(1)
plot(t1,1+x2,'r','linewidth',1.5)
ylim([-0.5 4])
xlim([0 24])
legend({'taxol in the medium','location'},'location','south','fontsize',16)
xlabel('time (hours)','fontsize',16)
set(gca,'fontsize',16)
yticks([1 3])
xticks(0:4:24)
set(gca,'YTickLabel',{'off','on'})
file_name = [save_file 'protocol1.eps'];
export_fig(figure(1),file_name)
%% Nocodazole uptake
figure(2)
plot(t1+1,1+x2,'r','linewidth',1.5)
hold on
plot(t4+1,1+a1,'k','linewidth',1.5)
hold on
plot(t7+1,x7,'r','linewidth',1.5)

xlim([0,25])
ylim([-0.5 4])
legend({'taxol in the medium','nocodazole in the medium'},'location','south','fontsize',16)
xlabel('time (hours)','fontsize',16)
set(gca,'fontsize',16)
yticks([1 2 3])
set(gca,'YTickLabel',{'off','on','on'})
file_name = [save_file 'protocol2.eps'];
export_fig(figure(2),file_name)

%% Nocodazole washout 1
figure(3)
plot(t5+3,1+x,'r','linewidth',1.5)
hold on
plot(t3+3,1+z1,'k','linewidth',1.5)
hold on
plot(t7+1,x7,'r','linewidth',1.5)
hold on
plot(t1+3,x1,'r','linewidth',1.5)

xlim([0,27])
ylim([-0.5 4])
legend({'taxol in the medium','nocodazole in the medium'},'location','south','fontsize',16)
xlabel('time (hours)','fontsize',16)
set(gca,'fontsize',16)
yticks([1 2 3])
set(gca,'YTickLabel',{'off','on','on'})
file_name = [save_file 'protocol3.eps'];
export_fig(figure(3),file_name)

%% Normal washout
figure(4)
plot(t5+2,1+x,'r','linewidth',1.5)
hold on
plot(t1+2,x1,'r','linewidth',1.5)

xlim([0 26])
ylim([-0.5 4])
legend({'taxol in the medium'},'location','south','fontsize',16)
xlabel('time (hours)','fontsize',16)
set(gca,'fontsize',16)
yticks([1 3])
xticks(0:4:28)
set(gca,'YTickLabel',{'off','on'})
file_name = [save_file 'protocol4.eps'];
export_fig(figure(4),file_name)
%% Nocodazole washout 2
figure(5)
plot(t5+2,1+x,'r','linewidth',1.5)
hold on
plot(t1+2,1+x1,'k','linewidth',1.5)
hold on
plot(t1+2,x1,'r','linewidth',1.5)
hold on
plot(t5+2,x-1,'k','linewidth',1.5)

xlim([0 26])
ylim([-0.5 4])
legend({'taxol in the medium','nocodazole in the medium'},'location','south','fontsize',16)
xlabel('time (hours)','fontsize',16)
set(gca,'fontsize',16)
yticks([1 2 3])
set(gca,'YTickLabel',{'off','on','on'})
file_name = [save_file 'protocol5.eps'];
export_fig(figure(5),file_name)

%% 2 hour treatment
figure(6)
plot(t6,1+y,'r','linewidth',1.5)
xlim([0,2])
ylim([-0.5 4])
legend({'taxol in the medium'},'location','south','fontsize',16)
xlabel('time (hours)','fontsize',16)
set(gca,'fontsize',16)
yticks([1 3])
set(gca,'YTickLabel',{'off','on'})
file_name = [save_file 'protocol6.eps'];
export_fig(figure(6),file_name)
