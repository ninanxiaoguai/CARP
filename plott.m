clear,clc
close all

load('spea_all.dat');
% nsga_all = spea_all;
load('nsga_all.dat');
% load('mix_all.dat');
% mix_all = nsga_all;
mix_all = spea_all;
figure
set(gcf,'position',[50 50 1200 450])
step = 60;
iter = size(mix_all,1)/step;
idx = 0;
for i = 1:5:iter
    
    ind = ((i-1)*step+1):i*step;
%     subplot(1,3,1)
%     f = spea_all(ind,3) > 0;
%     plot(spea_all(ind,1),spea_all(ind,2),'ro')
%     hold on
%     plot(spea_all(ind(1:60),1),spea_all(ind(1:60),2),'g*')
%     hold on
%     plot(spea_all(ind(f),1),spea_all(ind(f),2),'ks')
%     title('spea')
%     hold off
    
%     subplot(1,2,1)
    
%     f = mix_all(ind,3) > 0;
%     if(1 == 2)
%         plot3(mix_all(ind,1),mix_all(ind,2),mix_all(ind,3),'ro')
%         hold on
%         plot3(mix_all(ind(1:step),1),mix_all(ind(1:step),2),mix_all(ind(1:step),3),'g*')
%         hold on
%         plot3(mix_all(ind(f),1),mix_all(ind(f),2),mix_all(ind(f),3),'ks')
%         name = ['mix: ' num2str(i)];
%     else
%         subplot(1,3,1)
%         plot(mix_all(ind,1),mix_all(ind,2),'ro')
%         hold on
%         plot(mix_all(ind(1:step),1),mix_all(ind(1:step),2),'g*')
%         hold on
%         plot(mix_all(ind(f),1),mix_all(ind(f),2),'ks')
%         hold off
%         
%         subplot(1,3,2)
%         plot(mix_all(ind(1:step),1),mix_all(ind(1:step),2),'g*')
%         hold off
%         name = ['mix: ' num2str(i)];
%         
%         subplot(1,3,3)
%         final = (mix_all(ind(1:step),:));
%         del = unique(final,'rows');
%         plot(final(:,1),final(:,2),'ro')
%         num_clusters = [];
%         for ii = 1:size(del,1)
%             num_clusters =[num_clusters sum(all( final == repmat(del(ii,:),size(final,1),1),2 ))];
%         end
%         for ii = 1:size(del,1)
%             hold on
%             text(del(ii,1)+15,del(ii,2),num2str(num_clusters(ii)));
%         end
%         hold off
%     end
        
        subplot(1,2,1)
        final = (nsga_all(ind(1:step),:));
        del = unique(final,'rows');
        plot(final(:,1),final(:,2),'g*')
        num_clusters = [];
        for ii = 1:size(del,1)
            num_clusters =[num_clusters sum(all( final == repmat(del(ii,:),size(final,1),1),2 ))];
        end
        for ii = 1:size(del,1)
            hold on
            text(del(ii,1)+15,del(ii,2),num2str(num_clusters(ii)));
        end
        hold off
        name = [ 'NSGAII gen: ' num2str(i)];
        title(name)

        
        subplot(1,2,2)
        final = (mix_all(ind(1:step),:));
        del = unique(final,'rows');
        plot(final(:,1),final(:,2),'ro')
        num_clusters = [];
        for ii = 1:size(del,1)
            num_clusters =[num_clusters sum(all( final == repmat(del(ii,:),size(final,1),1),2 ))];
        end
        for ii = 1:size(del,1)
            hold on
            text(del(ii,1)+15,del(ii,2),num2str(num_clusters(ii)));
        end

        name = ['SPEA2 gen: ' num2str(i)];
        title(name)
        hold off
        pic_name = ['pic/' num2str(idx,'%03d') '.jpg'];
        print(gcf,'-djpeg',pic_name)
        idx = idx + 1;
    
%     subplot(1,2,2)
%     f = mix_all(ind,3) > 0;
%     plot3(nsga_all(ind,1),nsga_all(ind,2),nsga_all(ind,3),'ro')
%     hold on
%     plot3(nsga_all(ind(1:60),1),nsga_all(ind(1:60),2),nsga_all(ind(1:60),3),'g*')
%     hold on
%     plot3(nsga_all(ind(f),1),nsga_all(ind(f),2),nsga_all(ind(f),3),'ks')
%     title('nsga')
    hold off
    drawnow
end


% figure
% plot(mix_nondominated(:,1),mix_nondominated(:,2),'ks')
% hold on
% plot(mix_all(ind(1:60),1),mix_all(ind(1:60),2),'g*')


% plot(spea_all(:,1),spea_all(:,2),'r*')
% hold on
% plot(nsga_all(:,1),nsga_all(:,2),'go')
% hold on
% plot(mix_all(:,1),mix_all(:,2),'ks')
% legend('spea2','nsga','mix');

% n = unique(nsga_all,'rows');
% s = unique(spea_all,'rows');
% figure
% plot(s(:,1),s(:,2),'r*')
% hold on
% plot(n(:,1),n(:,2),'go')
% legend('spea2','nsga');


% 
% NSGA2 is done! time_cost: 85.570000
% total cost:  10156
% 
% SPEA2 is done! time_cost: 203.398000
% 
% mix is done! time_cost: 84.403000
% total cost:  10172
% n = load('nsga_finalpop_in_29.mat');
% s = load('spea_finalpop_in_29.mat');
% plot(s.final(:,1),s.final(:,2),'ro')
% hold on
% plot(n.final(:,1),n.final(:,2),'g*')
% legend('SPEA2','NSGA2')




