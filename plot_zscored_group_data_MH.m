function plot_zscored_group_data_Maryam(Stim_data,t_trials,Outcome,Condition,IdChannel,trial_mode,avg_mode,colors2plot,order_crit,limits2plot,Color_scale,show_plot,save_plot,PATH2SAVE)

%% Plot the average responses and the moving correlation
figure
if show_plot == 0
    set(gcf,'visible','off')
end

for o = 1:length(IdChannel)
    if ~isempty(IdChannel{o})
        tmp = Stim_data.(Condition).(IdChannel{o}).(trial_mode);
        tmp_t_trials = t_trials;
        tmp_t_trials(sum(isnan(tmp),2) == size(tmp,2)) = [];
        tmp(sum(isnan(tmp),2) == size(tmp,2),:) = [];
        if ~isempty(tmp)
            if avg_mode == 1
                tmp_avg = mean(tmp,2,'omitnan');
                tmp_error = std(tmp,1,2,'omitnan')./...
                    sqrt(sum(~isnan(tmp(1,:))));
            elseif avg_mode == 2
                tmp_avg = median(tmp,2,'omitnan');
                tmp_error = mad(tmp,1,2,'omitnan');
            end
            error_area(tmp_t_trials,tmp_avg,tmp_error,colors2plot{o},0.25)
        end
        hold on
    end
end
box off
if ~isempty(tmp_t_trials)
    xlim([tmp_t_trials(1) tmp_t_trials(end)])
end
if ~isempty(limits2plot.all)
    ylim(limits2plot.all)
end
xlabel('Time (s)')
ylabel({'Signal size','(A.U.)'})
title(['Avg sensors ',trial_mode,' aligned to ',Condition],'Interpreter', 'none')
ylimits = get(gca,'YLim');
plot([0 0],ylimits,'k')
%     xline(0,'-k');
if save_plot == 1
    saveas(gcf,[PATH2SAVE,Outcome,' ','Zscored Average signals centered on ',Condition,' ',trial_mode,'.jpg'])
    saveas(gcf,[PATH2SAVE,Outcome,' ','Zscored Average signals centered on ',Condition,' ',trial_mode,'.fig'])
end



%% Plot the heatmaps
figure
if show_plot == 0
    set(gcf,'visible','off')
end

for o = 1:length(IdChannel)
    if ~isempty(IdChannel{o})
        if order_crit == 1
            [~,order] = sort(abs(min(Stim_data.(Condition).(IdChannel{o}).(trial_mode)(:,t_trials > 0 & t_trials < 10),[],2)));
            order(isnan(Stim_data.(Condition).(IdChannel{o}).(trial_mode)(order,1))) = [];
        else
            order = 1:size(Stim_data.(Condition).(IdChannel{o}).(trial_mode),2);
        end
        subplot(length(IdChannel),1,o)
        imagesc(t_trials,1,Stim_data.(Condition).(IdChannel{o}).(trial_mode)(:,order)')
        colormap('jet')
        c = colorbar;
        c.Label.String = 'Intensity (A.U.)';
        if ~isempty(Color_scale)
            c_limits = Color_scale;
            caxis(c_limits)
        end
        hold on
        xlim([t_trials(1) t_trials(end)])
        xlabel('Time (s)')
        ylabel('Animal #')
        title([IdChannel{o},' ',trial_mode,' aligned to ',Condition],'Interpreter', 'none')
        box off
        ylimits = get(gca,'YLim');
        plot([0 0],ylimits,'k','LineWidth',2)
        %     xline(0,'-k');
    end
end
if save_plot == 1
    saveas(gcf,[PATH2SAVE,Outcome,' ','Heatmaps of Zscored dFF signals centered on ',Condition,' ',trial_mode,'.jpg'])
    saveas(gcf,[PATH2SAVE,Outcome,' ','Heatmaps of Zscored dFF signals centered on ',Condition,' ',trial_mode,'.fig'])
end

%% Plot the behavior and the individual traces
figure
if show_plot == 0
    set(gcf,'visible','off')
end

for o = 1:length(IdChannel)
    if ~isempty(IdChannel{o})
        subplot(length(IdChannel),1,o)
        plot(t_trials,Stim_data.(Condition).(IdChannel{o}).(trial_mode))
        hold on
        if avg_mode == 1
            plot(t_trials,nanmean(Stim_data.(Condition).(IdChannel{o}).(trial_mode),2),colors2plot{o},...
                'LineWidth',1.5)
        elseif avg_mode == 2
            plot(t_trials,nanmedian(Stim_data.(Condition).(IdChannel{o}).(trial_mode),2),colors2plot{o},...
                'LineWidth',1.5)
        end
        box off
        %                 plot([10 10],limits2plot.dLight,'k')
        xlim([t_trials(1) t_trials(end)])
        if ~isempty(limits2plot.(IdChannel{o}))
            ylim(limits2plot.(IdChannel{o}))
        end
        ylabel({'Signal size','(A.U.)'})
        title([IdChannel{o},' ',trial_mode,' aligned to ',Condition],'Interpreter', 'none')
        ylimits = get(gca,'YLim');
        plot([0 0],ylimits,'k')
        %     xline(0,'-k');
    end
end
if save_plot == 1
    saveas(gcf,[PATH2SAVE,Outcome,' ','Individual Zscored traces centered on ',Condition,' ',trial_mode,'.jpg'])
    saveas(gcf,[PATH2SAVE,Outcome,' ','Individual Zscored traces centered on ',Condition,' ',trial_mode,'.fig'])
end
end
