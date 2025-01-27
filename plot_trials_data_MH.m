function plot_trials_data_Maryam(Stim_data,t_trials,Outcome,Condition,IdChannel,trial_mode,avg_mode,colors2plot,order_crit,limits2plot,LimitXaxis,Color_scale,show_plot,save_plot,PATH2SAVE)

%% Plot the average responses and the moving correlation
figure
if show_plot == 0
    set(gcf,'visible','off')
end
if sum(sum(~isnan((Stim_data.(Condition).corr)))) > 1
    subplot(1,2,1)
    for o = 1:length(IdChannel)
        if avg_mode == 1
            tmp_avg = nanmean(Stim_data.(Condition).dFF.(IdChannel{o}).(trial_mode),1);
            tmp_error = nanstd(Stim_data.(Condition).dFF.(IdChannel{o}).(trial_mode),1,1)./...
                sqrt(sum(~isnan(Stim_data.(Condition).dFF.(IdChannel{o}).(trial_mode)(:,1))));
        elseif avg_mode == 2
            tmp_avg = nanmedian(Stim_data.(Condition).dFF.(IdChannel{o}).(trial_mode),1);
            tmp_error = mad(Stim_data.(Condition).dFF.(IdChannel{o}).(trial_mode),1,1);
        end
        error_area(t_trials,tmp_avg,tmp_error,colors2plot{o},0.25)
        hold on
    end
    % tmp_median = nanmedian(Stim_data.(Condition).dFF.(IdChannel{2}).(trial_mode),1);
    % tmp_mad = mad(Stim_data.(Condition).dFF.(IdChannel{2}).(trial_mode),1,1);
    % error_area(t_trials,tmp_median,tmp_mad,'r',0.25)
    CurrentyLimits = get(gca,'YLim');
    if ~isempty(limits2plot.all)
        plot([0 0],limits2plot.all,'-k')
    else
        plot([0 0],CurrentyLimits,'-k')
    end
    box off
    if ~isempty(LimitXaxis)
        xlim([LimitXaxis(1) LimitXaxis(end)])
    else
        xlim([t_trials(1) t_trials(end)])
    end
    
    xlabel('Time (s)')
    ylabel({'Signal size','(A.U.)'})
    title(['Avg sensors aligned to ',Condition],'Interpreter', 'none')
    subplot(1,2,2)
    if avg_mode == 1
        tmp_avg = nanmean(Stim_data.(Condition).corr,1);
        tmp_error = nanstd(Stim_data.(Condition).corr,1,1)./...
            sqrt(sum(~isnan(Stim_data.(Condition).corr(:,1))));
    elseif avg_mode == 2
        tmp_avg = nanmedian(Stim_data.(Condition).corr,1);
        tmp_error = mad(Stim_data.(Condition).corr,1,1);
    end
    error_area(t_trials,tmp_avg,tmp_error,'k',0.25)
    box off
    if ~isempty(LimitXaxis)
        xlim([LimitXaxis(1) LimitXaxis(end)])
    else
        xlim([t_trials(1) t_trials(end)])
    end
    ylim([-1 1])
    xlabel('Time (s)')
    ylabel('r')
    title('Moving Correlation')
    ylimits = get(gca,'YLim');
    plot([0 0],ylimits,'k')
    %     xline(0,'-k');
    if save_plot == 1
        saveas(gcf,[PATH2SAVE,Outcome,' ','Average signals and moving correlation centered on ',Condition,' ',trial_mode,'.jpg'])
        saveas(gcf,[PATH2SAVE,Outcome,' ','Average signals and moving correlation centered on ',Condition,' ',trial_mode,'.fig'])
    end
else
    for o = 1:length(IdChannel)
        if ~isempty(IdChannel{o})
            if avg_mode == 1
                tmp_avg = nanmean(Stim_data.(Condition).dFF.(IdChannel{o}).(trial_mode),1);
                tmp_error = nanstd(Stim_data.(Condition).dFF.(IdChannel{o}).(trial_mode),1,1)./...
                    sqrt(sum(~isnan(Stim_data.(Condition).dFF.(IdChannel{o}).(trial_mode)(:,1))));
            elseif avg_mode == 2
                tmp_avg = nanmedian(Stim_data.(Condition).dFF.(IdChannel{o}).(trial_mode),1);
                tmp_error = mad(Stim_data.(Condition).dFF.(IdChannel{o}).(trial_mode),1,1);
            end
            error_area(t_trials,tmp_avg,tmp_error,colors2plot{o},0.25)
            hold on
        end
    end
    box off
    if ~isempty(LimitXaxis)
        xlim([LimitXaxis(1) LimitXaxis(end)])
    else
        xlim([t_trials(1) t_trials(end)])
    end
    if ~isempty(limits2plot.all)
        ylim(limits2plot.all)
    end
    xlabel('Time (s)')
    ylabel({'Signal size','(A.U.)'})
    title(['Avg sensors aligned to ',Condition],'Interpreter', 'none')
    ylimits = get(gca,'YLim');
    plot([0 0],ylimits,'k')
%     xline(0,'-k');
    if save_plot == 1
        saveas(gcf,[PATH2SAVE,Outcome,' ','Average signals centered on ',Condition,' ',trial_mode,'.jpg'])
        saveas(gcf,[PATH2SAVE,Outcome,' ','Average signals centered on ',Condition,' ',trial_mode,'.fig'])
    end
end


%% Plot the heatmaps
figure
if show_plot == 0
    set(gcf,'visible','off')
end
% for o = 1:length(IdChannel)
%     if ~isempty(IdChannel{o})
%         [~,order] = sort(abs(min(Stim_data.(Condition).dFF.(IdChannel{o}).(trial_mode)(:,t_trials > 0 & t_trials < 10),[],2)));
%     end
% end
% order(isnan(Stim_data.(Condition).dFF.(IdChannel{o}).(trial_mode)(order,1))) = [];

if sum(sum(~isnan((Stim_data.(Condition).corr)))) > 1
    for o = 1:length(IdChannel)
        if order_crit == 1
            [~,order] = sort(abs(min(Stim_data.(Condition).dFF.(IdChannel{o}).(trial_mode).(Session)(:,t_trials > 0 & t_trials < 10),[],2)));
            order(isnan(Stim_data.(Condition).dFF.(IdChannel{o}).(trial_mode).(Session)(order,1))) = [];
        else
            order = 1:size(Stim_data.(Condition).dFF.(IdChannel{o}).(trial_mode),1);
        end
        subplot(length(IdChannel)+1,1,o)
        imagesc(t_trials,1,Stim_data.(Condition).dFF.(IdChannel{o}).(trial_mode)(order,:))
        colormap('jet')
        c = colorbar;
        c.Label.String = 'Intensity (A.U.)';
        if ~isempty(Color_scale)
            c_limits = Color_scale;
            caxis(c_limits)
        end
        hold on
        if ~isempty(LimitXaxis)
            xlim([LimitXaxis(1) LimitXaxis(end)])
        else
            xlim([t_trials(1) t_trials(end)])
        end
        xlabel('Time (s)')
        ylabel('Trial')
        title([IdChannel{1},' and ',IdChannel{2},' aligned to ',Condition],'Interpreter', 'none')
        box off
        ylimits = get(gca,'YLim');
        plot([0 0],ylimits,'k','LineWidth',2)
        %     xline(0,'-k');
    end
    subplot(length(IdChannel)+1,1,length(IdChannel)+1)
    imagesc(t_trials,1,Stim_data.(Condition).corr(order,:))
    colormap('jet')
    c = colorbar;
    c.Label.String = 'r';
    caxis([-1 1])
    hold on
    if ~isempty(LimitXaxis)
        xlim([LimitXaxis(1) LimitXaxis(end)])
    else
        xlim([t_trials(1) t_trials(end)])
    end
    xlabel('Time (s)')
    ylabel('Trial')
    title(['Correlation aligned to ',Condition],'Interpreter', 'none')
    box off
    ylimits = get(gca,'YLim');
    plot([0 0],ylimits,'k','LineWidth',2)
    %     xline(0,'-k');
    if save_plot == 1
        saveas(gcf,[PATH2SAVE,Outcome,' ','Heatmaps of dFF signals and moving correlation centered on ',Condition,' ',trial_mode,'.jpg'])
        saveas(gcf,[PATH2SAVE,Outcome,' ','Heatmaps of dFF signals and moving correlation centered on ',Condition,' ',trial_mode,'.fig'])
    end
else
    for o = 1:length(IdChannel)
        if ~isempty(IdChannel{o})
            order = 1:size(Stim_data.(Condition).dFF.(IdChannel{o}).(trial_mode),1);
            subplot(length(IdChannel),1,o)
            imagesc(t_trials,1,Stim_data.(Condition).dFF.(IdChannel{o}).(trial_mode)(order,:))
            colormap('jet')
            c = colorbar;
            c.Label.String = 'Intensity (A.U.)';
            if ~isempty(Color_scale)
                c_limits = Color_scale;
                caxis(c_limits)
            end
            hold on
            if ~isempty(LimitXaxis)
                xlim([LimitXaxis(1) LimitXaxis(end)])
            else
                xlim([t_trials(1) t_trials(end)])
            end
            xlabel('Time (s)')
            ylabel('Trial')
            title([IdChannel{o},' aligned to ',Condition],'Interpreter', 'none')
            box off
            ylimits = get(gca,'YLim');
            plot([0 0],ylimits,'k','LineWidth',2)
            %     xline(0,'-k');
        end
    end
    if save_plot == 1
        saveas(gcf,[PATH2SAVE,Outcome,' ','Heatmaps of dFF signals centered on ',Condition,' ',trial_mode,'.jpg'])
        saveas(gcf,[PATH2SAVE,Outcome,' ','Heatmaps of dFF signals centered on ',Condition,' ',trial_mode,'.fig'])
    end
end

%% Plot the behavior and the individual traces
figure
if show_plot == 0
    set(gcf,'visible','off')
end

for o = 1:length(IdChannel)
    if ~isempty(IdChannel{o})
        subplot(length(IdChannel),1,o)
        plot(t_trials,Stim_data.(Condition).dFF.(IdChannel{o}).(trial_mode))
        hold on
        if avg_mode == 1
            plot(t_trials,nanmean(Stim_data.(Condition).dFF.(IdChannel{o}).(trial_mode),1),colors2plot{o},...
                'LineWidth',1.5)
        elseif avg_mode == 2
            plot(t_trials,nanmedian(Stim_data.(Condition).dFF.(IdChannel{o}).(trial_mode),1),colors2plot{o},...
                'LineWidth',1.5)
        end
        box off
        %                 plot([10 10],limits2plot.dLight,'k')
        if ~isempty(LimitXaxis)
            xlim([LimitXaxis(1) LimitXaxis(end)])
        else
            xlim([t_trials(1) t_trials(end)])
        end
        if ~isempty(limits2plot.(IdChannel{o}))
            ylim(limits2plot.(IdChannel{o}))
        end
        ylabel({'Signal size','(A.U.)'})
        title([IdChannel{o},' aligned to ',Condition],'Interpreter', 'none')
        ylimits = get(gca,'YLim');
        plot([0 0],ylimits,'k')
        %     xline(0,'-k');
    end
end
if save_plot == 1
    saveas(gcf,[PATH2SAVE,Outcome,' ','Individual traces centered on ',Condition,' ',trial_mode,'.jpg'])
    saveas(gcf,[PATH2SAVE,Outcome,' ','Individual traces centered on ',Condition,' ',trial_mode,'.fig'])
end
end
