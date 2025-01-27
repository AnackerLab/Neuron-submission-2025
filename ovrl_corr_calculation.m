function [lag2plot,Ovrl_corr] = ovrl_corr_calculation(signal1,signal2,time,...
    Max_lag,corr_type,show_plot,save_plot,PATH2SAVE,save_name)
% Define lags and index for lags
tmp_t = time - time(1);
dummie = 1:length(time);
dummie = dummie(tmp_t >= Max_lag);
idx_max_lag = dummie(1);
lag_steps = idx_max_lag*-1:idx_max_lag;
lag2plot = tmp_t(1:(idx_max_lag*2)+1)-tmp_t(idx_max_lag+1);

% Compute correlation with lags
Ovrl_corr.r = ones(1,length(lag_steps))*nan;
Ovrl_corr.p = ones(1,length(lag_steps))*nan;
for o = 1:length(lag_steps)
    data1 = signal1;
    data2 = circshift(signal2,lag_steps(o));
    [Ovrl_corr.r(o),Ovrl_corr.p(o)] = corr(data1',data2','Type',corr_type);
end

[~,idx] = findpeaks(abs(Ovrl_corr.r));
if ~isempty(idx)
    tmp = Ovrl_corr.r(idx);
    [~,ix] = max(tmp);
    [~,ix2] = find(abs(Max_lag-lag2plot(idx(ix))) == min(abs(Max_lag-lag2plot(idx(ix)))));
    Ovrl_corr.positive.value = Ovrl_corr.r(idx(ix(ix2)));
    Ovrl_corr.positive.lag.idx = idx(ix(ix2));
    Ovrl_corr.positive.lag.val = lag2plot(idx(ix(ix2)));
    
    [~,ix] = min(tmp);
    [~,ix2] = find(abs(Max_lag-lag2plot(idx(ix))) == min(abs(Max_lag-lag2plot(idx(ix)))));
    Ovrl_corr.negative.value = Ovrl_corr.r(idx(ix(ix2)));
    Ovrl_corr.negative.lag.idx = idx(ix(ix2));
    Ovrl_corr.negative.lag.val = lag2plot(idx(ix(ix2)));
else
    Ovrl_corr.positive.value = [];
    Ovrl_corr.positive.lag.idx = [];
    Ovrl_corr.positive.lag.val = [];
    Ovrl_corr.negative.value = [];
    Ovrl_corr.negative.lag.idx = [];
    Ovrl_corr.negative.lag.val = [];
end

% first_half = abs(Ovrl_corr.r(1:idx_max_lag+1));
% [~,idx] = findpeaks(first_half);
% if ~isempty(idx)
%     idx = idx(end);
%     Ovrl_corr.negative.value = Ovrl_corr.r(idx);
%     Ovrl_corr.negative.lag.idx = idx;
%     Ovrl_corr.negative.lag.val = lag2plot(idx);
% else
%     Ovrl_corr.negative.value = [];
%     Ovrl_corr.negative.lag.idx = [];
%     Ovrl_corr.negative.lag.val = [];
% end
% 
% second_half = Ovrl_corr.r(idx_max_lag+1:end);
% [~,idx] = findpeaks(second_half);
% if ~isempty(idx)
%     idx = idx(1);
%     Ovrl_corr.positive.value = Ovrl_corr.r(idx+idx_max_lag);
%     Ovrl_corr.positive.lag.idx = idx+idx_max_lag;
%     Ovrl_corr.positive.lag.val = lag2plot(idx+idx_max_lag);
% else
%     Ovrl_corr.positive.value = [];
%     Ovrl_corr.positive.lag.idx = [];
%     Ovrl_corr.positive.lag.val = [];
% end
% Plot the result
figure
if show_plot == 0
    set(gcf,'visible','off')
end
plot(lag2plot,Ovrl_corr.r,'k')
% ylim([min(Ovrl_corr.r)-0.1 max(Ovrl_corr.r)+0.1])
ylim([-1 1])
Y_lim = ylim;
xlim([lag2plot(1) lag2plot(end)])
hold on
plot([0 0],Y_lim,'--k')
if ~isempty(Ovrl_corr.negative.lag.val)
    plot([Ovrl_corr.negative.lag.val Ovrl_corr.negative.lag.val],Y_lim,'r')
end
if ~isempty(Ovrl_corr.positive.lag.val)
    plot([Ovrl_corr.positive.lag.val Ovrl_corr.positive.lag.val],Y_lim,'r')
end
xlabel('Lag (s)')
ylabel('r')
title(save_name)
if ~isempty(Ovrl_corr.negative.lag.val)
    if Ovrl_corr.negative.lag.val <= 4 || Ovrl_corr.negative.lag.val >= -4 
        txt = ['- lag = ',num2str(Ovrl_corr.negative.lag.val),'msec \rightarrow'];
        text(Ovrl_corr.negative.lag.val,min(Ovrl_corr.r)-0.05,txt,'HorizontalAlignment','right')
    else
        txt = ['\leftarrow - lag = ',num2str(Ovrl_corr.negative.lag.val),'msec'];
        text(Ovrl_corr.negative.lag.val,min(Ovrl_corr.r)-0.05,txt,'HorizontalAlignment','right')
    end
end
if ~isempty(Ovrl_corr.positive.lag.val)
    if Ovrl_corr.positive.lag.val <= 4 || Ovrl_corr.positive.lag.val >= -4 
        txt = ['+ lag = ',num2str(Ovrl_corr.positive.lag.val),'msec \rightarrow'];
        text(Ovrl_corr.positive.lag.val,max(Ovrl_corr.r)+0.05,txt,'HorizontalAlignment','right')
    else
        txt = ['\leftarrow + lag = ',num2str(Ovrl_corr.positive.lag.val),'msec'];
        text(Ovrl_corr.positive.lag.val,max(Ovrl_corr.r)+0.05,txt,'HorizontalAlignment','right')
    end
end
box off
if save_plot
    saveas(gcf,[PATH2SAVE,'figures\',save_name,'.jpg'])
    saveas(gcf,[PATH2SAVE,'figures\',save_name,'.fig'])
end

end