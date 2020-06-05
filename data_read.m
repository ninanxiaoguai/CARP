clc,clear
problems = 0:80;
runs = 0:9;
show = 0;
%% 收集 PF
res = [];
nsga_best = [];
spea_best = [];
for problem = problems
    res = [];
    res1 = [];
    res2 = [];
    for run = runs
        name = [ num2str(problem), '_mix_nondominated' num2str(run) '.dat'];
        fid = fopen(name);
        C = textscan(fid,'%d%d');
        fclose(fid);
        res = [res; cell2mat(C)];
        res1 = [res1;cell2mat(C)];
        
        name = [ num2str(problem), '_nsga_nondominated' num2str(run) '.dat'];
        fid = fopen(name);
        C = textscan(fid,'%f%f');
        fclose(fid);
        res = [res; cell2mat(C)];
        res2 = [res2;cell2mat(C)];
        
    end
    res = unique(res,'rows');
    
    res1 = [res1(:,2) res1(:,1)];
    res2 = [res2(:,2) res2(:,1)];
    
    res1 = unique(res1,'rows');
    res2 = unique(res2,'rows');
    
    nsga_best = [nsga_best;res1(1,:)];
    spea_best = [spea_best;res2(1,:)];
    [FrontNo,MaxFNo] = NDSort(res,size(res,1));
    PF = double(res(FrontNo==1,:));
    name = ['C:\\Users\\QingQuan\\Desktop\\tt\\res\\pareto\\',num2str(problem), '.mat'];
    save(name, 'PF');
end
%%
% ii = 1;
% for i = problems
%     name = ['C:\Users\QingQuan\Desktop\tt\res\pareto\' num2str(i) '.mat'];
%     temp = load(name);
%     data = temp.PF;
%     subplot(5,5,ii)
%     plot(data(:,1),data(:,2),'ks')
%     title(num2str(i+1))
%     ii = ii + 1;
%     close all
% end
%%
NSGA_i = [];
SPEA_i = [];
MIX_i = [];
NSGA_h = [];
SPEA_h = [];
MIX_h = [];
NSGA_s = [];
SPEA_s = [];
MIX_s = [];
for problem = problems
    res = [];
    name = ['C:\Users\QingQuan\Desktop\tt\res\pareto\' num2str(problem) '.mat'];
    temp = load(name);
    PF = temp.PF;
    ma = max(PF,[],1);
    mi = min(PF,[],1);
    PF = (PF - repmat(mi,size(PF,1),1))./repmat(ma-mi+1,size(PF,1),1);
    mix_i = [];
    spea_i = [];
    nsga_i = [];
    mix_h = [];
    spea_h = [];
    nsga_h = [];
    mix_s = [];
    spea_s = [];
    nsga_s = [];
    for run = runs
        name = [ num2str(problem), '_mix_nondominated' num2str(run) '.dat'];
        fid = fopen(name);
        C = textscan(fid,'%d%d');
        fclose(fid);
        pop1 = double(cell2mat(C));
        pop = (pop1 - repmat(mi,size(pop1,1),1))./repmat(ma-mi + 1,size(pop1,1),1);
        dis = pdist2(PF,pop);
        dis = sort(dis,2);
        mix_i = [mix_i mean(dis(:,1))];
        mix_h = [mix_h HV(pop1,ma)];
        ss = min(pop1,[],1);
        mix_s = [mix_s ss(2)];
        if(show == 1)
            plot(pop(:,1),pop(:,2),'ks')
        end
        
        name = [ num2str(problem), '_nsga_nondominated' num2str(run) '.dat'];
        fid = fopen(name);
        C = textscan(fid,'%d%d');
        fclose(fid);
        pop1 = double(cell2mat(C));
        pop = (pop1 - repmat(mi,size(pop1,1),1))./repmat(ma-mi+1,size(pop1,1),1);
        dis = pdist2(PF,pop);
        dis = sort(dis,2);
        nsga_i = [nsga_i mean(dis(:,1))];
        nsga_h = [nsga_h HV(pop1,ma)];
        ss = min(pop1,[],1);
        nsga_s = [nsga_s ss(2)];
        if(show == 1)
            hold on
            plot(pop(:,1),pop(:,2),'go')
            hold on
            plot(PF(:,1),PF(:,2),'b^')
            legend('mix','nsga','PF')
            title(num2str(problem))
            hold off
            drawnow
            %close all
        end
    end
    MIX_i = [MIX_i;mix_i];
    NSGA_i = [NSGA_i;nsga_i];
    
    MIX_h = [MIX_h;mix_h];
    NSGA_h = [NSGA_h;nsga_h];
    
    MIX_s = [MIX_s;mix_s];
    NSGA_s = [NSGA_s;nsga_s];
end

H_mi = zeros(1,size(MIX_i,1));
H_mh = zeros(1,size(MIX_i,1));
H_ms = zeros(1,size(MIX_i,1));


% 在 IGD 上，统计mix相对于nsga的Wilson rank
for i = 1:size(MIX_i,1)
    [~,h] = ranksum(MIX_i(i,:),NSGA_i(i,:));
    if(h == 1)
        if( mean(MIX_i(i,:)) < mean(NSGA_i(i,:)) )
            H_mi(i) = 1;
        else
            H_mi(i) = -1;
        end
    else
        H_mi(i) = 0;
    end
end

% 在 HV 上，统计mix相对于nsga的Wilson rank
for i = 1:size(MIX_i,1)
    [~,h] = ranksum(-MIX_h(i,:),-NSGA_h(i,:) );
    if(h == 1)
        if( mean(-MIX_h(i,:)) < mean(-NSGA_h(i,:)) )
            H_mh(i) = 1;
        else
            H_mh(i) = -1;
        end
    else
        H_mh(i) = 0;
    end
end

ALL_i = [mean(MIX_i,2) mean(NSGA_i,2)];
ALL_h = [mean(MIX_h,2)  mean(NSGA_h,2)];
ALL_sm = [mean(MIX_s,2) mean(NSGA_s,2)];
ALL_s = [min(MIX_s,[],2) min(NSGA_s,[],2)];
for i = 1:size(MIX_i,1)
    [~,h] = ranksum(MIX_s(i,:),NSGA_s(i,:) );
    if(h == 1)
        if( mean(MIX_h(i,:)) < mean(NSGA_s(i,:)) )
            H_ms(i) = 1;
        else
            H_ms(i) = -1;
        end
    else
        H_ms(i) = 0;
    end
end

fprintf('\n In the IGD, SPEA2 is better in %f , NSGA2 is better in %f\n ',sum([ALL_i(:,1) < ALL_i(:,2)])/length(problems),sum([ALL_i(:,1) > ALL_i(:,2)])/length(problems))
fprintf('In the HV, SPEA2 is better in %f , NSGA2 is better in %f\n ',sum([ALL_h(:,1) > ALL_h(:,2)])/length(problems),sum([ALL_h(:,1) < ALL_h(:,2)])/length(problems))
fprintf('In the mean best, SPEA2 is better in %f , NSGA2 is better in %f\n ',sum([ALL_sm(:,1) < ALL_sm(:,2)])/length(problems),sum([ALL_sm(:,1) > ALL_sm(:,2)])/length(problems))
fprintf('In the best one, SPEA2 is better in %f , NSGA2 is better in %f\n\n ',sum([ALL_s(:,1) < ALL_s(:,2)])/length(problems),sum([ALL_s(:,1) > ALL_s(:,2)])/length(problems))

fprintf('Wilcoxon Rank: 1 means SPEA2 is better; 0 means simliar; -1 means bad\n ')
fprintf('    IGD   HV  mean best \n')

wilcoxon = [H_mi' H_mh' H_ms'];
disp(wilcoxon)












