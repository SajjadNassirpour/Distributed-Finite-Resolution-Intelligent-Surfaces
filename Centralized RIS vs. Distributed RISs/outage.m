function Outage_capacity=outage(sum_rate_mat,network_case)


epsilon=0.9:-0.1:0;

sum_each_iter=sum(sum_rate_mat,2);

A=sort(sum_each_iter); 


Outage_capacity=zeros(1,length(epsilon));

for e=1:length(epsilon)
   
   idx_incre=0;
   
   for outage_idx=length(A):-1:1
       
       if round(sum(sum_each_iter>=A(outage_idx))/length(A)-(1-epsilon(e)),4)>=0 && idx_incre==0
            idx_incre=1;
            idx=outage_idx;
            
       end
   end

   Outage_capacity(e)=A(idx);


end


if strcmp(network_case,'center')

    figure(1)
else
    figure(2)
end

hold on
grid on

plot(1-epsilon,Outage_capacity ,'--k','linewidth',2)

xlabel('$1-\gamma$','Interpreter','latex')
ylabel('Outage Capacity','Interpreter','latex')
%title('K=3,  N=4, M=32, $P_{t}=20$dBm, $\sigma^2=-80$dBm, d=50m, d$_1$=5m','Interpreter','latex')
%ylim([2.9,4.3])
xlim([0.1,1])
%xlim([25,55])
box on
grid on
set(gcf,'color','w');
figureHandle = gcf;
set(findall(figureHandle,'type','text'),'Interpreter','latex','fontSize',18)
set(gca,'TickLabelInterpreter','latex','fontSize',20 )
if strcmp(network_case,'center')

    H9 = legend('A centralized RIS','Interpreter','latex')
else
    H9 = legend('K distributed RISs','Interpreter','latex')
end
set(H9,'Interpreter','latex','fontSize',20)
 