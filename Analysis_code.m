%% BEFORE YOU START THE FOLLOWING SCRIPT REQUIRES THE FLLOWING DATA STRUCTURES
%% X (X POSITION ),Y (Y POSITION),D (DELTA T BETWEEN TWO CONSECUTIVE FRAMES),V (VELOCITY),INFO (INFORMATIONS ABOUT THE ARENA SIZE)
%% each of these variable need to be a cell array with one experimental group per cell 
%% i.e X=1x3 cell if you have 3 experimental group
%%     X{1,1}= N x M cell array where M=number of animals, N= number of experimental sessions\trial performed
%% for savings of the figure set a variable filename= "path where you want to save your figures" and create a folder called figures
%% name_groups (identifyer of each individual experimental group, colorur (color id for each experimental group, used for plotting)
%% i.e. name_groups={'control', 'trained', 'sham'}; colour={'k','r','m'}
%% rule_change = session at which you invert the learning rule, if you do not, set it to 'NO'
time_bin=120000;%2 minutes time bin
base_time=1800000;%take into account only the last 30 minutes of the baseline 
%%CODE IMPLEMENTATION
%%
%%Learning index calculation
switch_rule=rule_change;
for gr=1:max(size(X))
     X_vector=[];Y_vector=[]; D_vector=[]; info_vector=[];
        X_vector=X{1,gr};Y_vector=Y{1,gr}; D_vector=D{1,gr};  info_vector=info{1,gr};
 Learning_index(gr)={LI_with_inversion_of_the_rule_baseline_refered(X_vector,Y_vector,D_vector,info_vector,switch_rule)};
end
%%
%%if you want to select only animals with a positive learning index
session_of_interest=4; % session in which you want to refer to for the learning index
for gr=1:max(size(X))
    LI= Learning_index{gr}(session_of_interest,:);
    animal_to_select=[];animal_to_select=find(LI>0);
     X_vector=[];Y_vector=[]; D_vector=[]; info_vector=[];V_vector=[];
      X_vector=X{1,gr};Y_vector=Y{1,gr}; D_vector=D{1,gr}; info_vector=info{1,gr};V_vector=V{1,gr};
X{1,gr}=X_vector(:,animal_to_select);
Y{1,gr}=Y_vector(:,animal_to_select);
D{1,gr}=D_vector(:,animal_to_select);
info{1,gr}=info_vector(:,animal_to_select);
V(1,gr)={V_vector(:,animal_to_select)};
end

%%
time_window=5; %in minutes
min_size_bin=100;%%min number of frames in a time bin (use usually half the frames in a timebin)
for gr=1:max(size(X))
X_vector=[];Y_vector=[];D_vector=[];info_vector=[];V_vector=[];V_vector=V{1,gr};
X_vector=X{1,gr};Y_vector=Y{1,gr};D_vector=D{1,gr};info_vector=info{1,gr};
[xNew{gr},yNew{gr},dNew{gr},~,vNew{gr}]=get_time_bins(X_vector,Y_vector,D_vector,V_vector,time_window,info_vector,min_size_bin);
end
%%
figure();
titles=["average time of entrance","number of entrance","v of entrance"];
baseline_shift=6;%the first 5 bins of baseline are negletted, check the time_window (those multiplied gives time deleted)
baseline_bins=baseline_shift+[1:6];%i am using a 5 min binn and selection only the last 30 min of the baseline 
for gr=1:max(size(X))
   P=[];H=[];
   IT=[];ITM=[];
   [IT,ITM]=mean_duration_entrance(xNew{1,gr},dNew{1,gr},info{1,gr}(:,:,1),rule_change,vNew{gr});
   %%IT duretion of each single entrance in conditioned zone
   % ITI (:,:,1)=mean time of entrance in seconds
   % ITI (:,:,2)=mean velocity of entrance in cm/s
   % ITI (:,:,3)=number of times the animal enters the conditioned zone
   % ITI (:,:,4)=standard error mean time of entrance
p=0;t=0;
for u=[1,3,4]
    p=p+1;t=t+1;
subplot(3,1,(p))
shadedErrorBar(1:size(ITM,1)-baseline_shift,nanmean(ITM((baseline_shift+1):end,:,u),2),[nanstd(ITM((baseline_shift+1):end,:,u),0,2)./sqrt(size(ITM((baseline_shift+1):end,:,u),2))],'lineprops',colour{gr} , 'patchSaturation', 0.2);
GR_IT(gr,u)={squeeze(ITM((baseline_shift+1):end,:,u))};
for v=13:size(ITM(:,:,u),1)
[P(v,u),H(v,u)]=signrank(mean(ITM(baseline_bins,:,u),1),ITM(v,:,u)); 
end
Pval(gr)={P};
Hval(gr)={H};
hold on, title(titles(t))
end
end
num_of_sessions=block_num;bins_per_session=6;
for u=[1,3,4]
for pv=1:num_of_sessions
    A=[];A=GR_IT{1,u}(((pv-1)*bins_per_session)+1:pv*bins_per_session,:);
    B=[];B=GR_IT{2,u}(((pv-1)*bins_per_session)+1:pv*bins_per_session,:);
[P12(u,pv),H12(u,pv)]=ranksum(A(:),B(:));
end
end
%%

%%
time_window=30; %in minutes we want the hetmaps to have a longer averaged time!
min_size_bin=100;%%min number of frames in a time bin (use usually half the frames in a timebin)
for gr=1:max(size(X))
X_vector=[];Y_vector=[];D_vector=[];info_vector=[];V_vector=[];V_vector=V{1,gr};
X_vector=X{1,gr};Y_vector=Y{1,gr};D_vector=D{1,gr};info_vector=info{1,gr};
[~,~,~,heatmaps{gr},~]=get_time_bins(X_vector,Y_vector,D_vector,V_vector,time_window,info_vector,min_size_bin);
end
%% save averaged and individual heatmaps of exploration patterns, we use the binned time setted above line 33 
sz=get(0,'MonitorPositions');
B=figure();
for gr=1:max(size(X))
pdf_vector=[]; X_vector=[];Y_vector=[]; D_vector=[];  info_vector=[];V_vector=[];
pdf_vector=heatmaps{1,gr}; X_vector=X{1,gr};Y_vector=Y{1,gr}; D_vector=D{1,gr};  info_vector=info{1,gr}; V_vector=V{1,gr};
name=['heat_maps_plot_average',name_groups{gr}]
figure(B)
PDF_mean=plot_heat_maps_averaged([pdf_vector(2:end,:)]);
plot_heat_map(PDF_mean,(gr),block_num,max(size(X)));%
saveas(B,[filename,'figures\',name,'.fig']);saveas(B,[filename,'figures\',name,'.tif']);saveas(B,[filename,'figures\',name,'.svg']); 
counts=0;add=0;A=figure();figure(A)
name_fig=['single_fish_plot_',name_groups{gr}];
for FISH=1:size(pdf_vector,2)
 if (counts==6 )
     saveas(A,[filename,'figures\',name_fig,num2str(FISH-1),'.fig'])
     saveas(A,[filename,'figures\',name_fig,num2str(FISH-1),'.tif'])
     saveas(A,[filename,'figures\',name_fig,num2str(FISH-1),'.svg'])
     close(A) 
     A=figure();
     figure(A)
     counts=0;
     add=add+6;
 end
 counts=counts+1;
   plot_heat_map(pdf_vector(:,add+counts),counts,block_num+1,6)
hold on
      if (FISH==size(pdf_vector,2) )
     saveas(A,[filename,'figures\',name_fig,num2str(FISH),'.fig'])
     saveas(A,[filename,'figures\',name_fig,num2str(FISH),'.tif'])
     saveas(A,[filename,'figures\',name_fig,num2str(FISH),'.svg'])
     close(A) 
     A=figure();
     figure(A)
     counts=0;
     add=add+6;
     end
end
end
%%
%% analyse the time in red for each group individually and compare them together
fig1=figure();
fig2=figure();
for gr=1:max(size(X))
    X_vector=[];Y_vector=[]; D_vector=[];  info_vector=[]; V_vector=[];
    X_vector=X{1,gr};Y_vector=Y{1,gr}; D_vector=D{1,gr};  info_vector=info{1,gr}; V_vector=V{1,gr};
name_fig1=['time binned in red plot individuals'];
name_fig2=['time binned in red plot all'];
binned_vector1=[];Flip_session=rule_change;
[binned_vector1]=new_binned_plot_with_flipping_color(X_vector,Y_vector,D_vector,info_vector,colour{gr},time_bin,base_time,Flip_session);
bin_vector(gr)={binned_vector1};
for p=1:size(binned_vector1,1)
if(p==1)
M=length(binned_vector1{p,1});
else
M= horzcat (M,M(p-1)+length(binned_vector1{p,1}));
end
end
figure(fig1)
subplot(max(size(X)),1,gr)
mean_binned1=[];
mean_binned1=plot_binned(binned_vector1,colour{gr},15); plot_boundary(M , [0,1])
hold on
bin_vector_mean(gr)={mean_binned1};
figure(fig2)
shadedErrorBar(1:size(mean_binned1,2),mean_binned1(1,:),mean_binned1(2,:),'lineprops',colour{gr} , 'patchSaturation', 0.2);
hold on ; plot_boundary(M , [0,1]) ; hold on; ylim([0 1]);ylabel('percentage time in danger')
if(gr==max(size(X)))
name=[];name=name_fig1;
title(name_fig1)
saveas(fig1,[filename,'figures\',name,'.fig']);saveas(fig1,[filename,'figures\',name,'.tif']);saveas(fig1,[filename,'figures\',name,'.svg']); 
name=[];name=name_fig2;
title(name_fig2)
saveas(fig2,[filename,'figures\',name,'.fig']);saveas(fig2,[filename,'figures\',name,'.tif']);saveas(fig2,[filename,'figures\',name,'.svg']); 
end
end
%%
%% Velocity plotting
fig1=figure();name_fig1=[];name_fig1=['binned velocity plot individuals'];
fig2=figure();name_fig2=[];name_fig2=['binned velocity plot all'];
fig3=figure();name_fig3=[];name_fig3=['total distance swimmed'];
fig4=figure();name_fig4=[];name_fig4=['freezing'];
fig5=figure();name_fig5=[];name_fig5=['individual animal dist in red'];
fig6=figure(); 
fig7=figure();
fig10=figure();name_fig10=[];name_fig10=['animal dist in red [cm] or ratio from baseline'];
sec_to_freez=2;dist_to_freez=2;
for gr=1:max(size(X))
    X_vector=[];Y_vector=[]; D_vector=[]; info_vector=[]; V_vector=[];
    X_vector=X{1,gr};Y_vector=Y{1,gr}; D_vector=D{1,gr};  info_vector=info{1,gr}; V_vector=V{1,gr};
    binned_vel1=[];Flip_session=rule_change;   
    normalized_to_baseline=0;
    limity=2;
    Distance1=[];Velocity1=[];binned_vel1=[];Dist_red=[];Dist_raw=[];binned_dist_red=[];binned_dist_all=[];
    [Dist_raw,Distance1,Velocity1,binned_vel1,Dist_red]=binned_plot_velocity(X_vector,Y_vector,D_vector,info_vector,normalized_to_baseline,rule_change,time_bin, base_time,V_vector); %%velocity output is in centimiters per seconds !
    Normalise_to_base=1;%used to normalize distace swim in conditioned zone with the distance swam in conditioned zone during the baseline 
    czcs=0;% boolean 1=ratio distance conditioned over safe distance, 0 only conditioned zone
    [binned_dist_red,binned_dist_all] = dist_in_red_time_bin(X_vector,Y_vector,D_vector,info_vector,normalized_to_baseline,rule_change,time_bin,base_time,Normalise_to_base,czcs);
    D_RED(gr)={binned_dist_red};
    D_ALL(gr)={binned_dist_all};
    Distance_raw(gr)={Dist_raw};
    Distance_ratio(gr)={Dist_red./(Distance1)};
    freez_in_time=[];totfreez=[];
    [freezing_time,freez_raw_time]=freezing_new(X_vector,Y_vector,D_vector,[sec_to_freez,dist_to_freez],info_vector);% dist to freeze in mm , sec to freeze in sec
    freez(1,gr)={freezing_time};
    freez(2,gr)={freez_raw_time};
    mean_binned1=[];all_dist_red=[];
    figure(fig5)
    [mean_binned1,all_dist_red]=plot_binned_dist(binned_dist_red,colour{gr},15);
    title('plot individual animal dist in red')
    hold on
    bin_vector_dist_red(gr)={mean_binned1};
    figure(fig10)
    shadedErrorBar(1:size(mean_binned1,2),mean_binned1(1,:),mean_binned1(2,:),'lineprops',colour{gr}, 'patchSaturation', 0.2);hold on; plot_boundary(M , [0, 2]);hold on;ylim([0 2])
    hold on
    title('animal dist in red [cm] or ratio from baseline')
    D_ratio=[];D_ratio=Distance_ratio{gr};
    figure(fig7)
    shadedErrorBar(1:size(Distance1,1),mean(Distance1,2),std(Distance1,0,2)./sqrt(size(Distance1,2)),'lineprops',colour{gr}, 'patchSaturation', 0.6);
    hold on
    title('total distance swam per session [m]')
    figure(fig6)
    shadedErrorBar(1:size(D_ratio,1),mean(D_ratio,2),std(D_ratio,0,2)./sqrt(size(D_ratio,2)),'lineprops',colour{gr}, 'patchSaturation', 0.2);
    hold on;xlabel('session in protocol ');ylabel('dist in danger/all dist');hold on;
    ylim([0 1])
    figure(fig4)
    totfreez=freezing_time;shadedErrorBar(1:size(totfreez,1),median(totfreez,2),[[median(totfreez,2)-prctile(totfreez,25,2)],[prctile(totfreez,75,2)-median(totfreez,2)]],'lineprops',colour{gr}, 'patchSaturation', 0.2);
    hold on
    xlabel('session ')
    ylabel('freezing time in %')
    hold on   
    bin_vel(gr)={binned_vel1}; velocity(gr)={Velocity1};distance(gr)={Distance1};distance_in_red(gr)={Dist_red};
    for p=1:size(binned_vel1,1)
        if(p==1)
            M=length(binned_vel1{p,1});
        else
            M= horzcat (M,M(p-1)+length(binned_vel1{p,1}));
        end
    end
    figure(fig3)
    pointer=1:1:size(Distance1,1);
    pointer=repmat(pointer,size(Distance1,2),1);
    for t=1:size(Distance1,1)
        scattercloud(pointer(:,t),[Distance1(t,:)]','filled',colour{gr},50)
        hold on
    end
    ylabel('total distance swimmed in meter')
    xlabel('different sessions of protocol')
    for t=1:size(Distance1,2)
        plot(pointer(1,:),[Distance1(:,t)]',colour{gr});
        hold on
    end
    hold on
    figure(fig1)
    subplot(max(size(X)),1,gr)
    mean_vel_binned1=plot_binned_vel(binned_vel1,colour{gr});plot_boundary(M , [0, limity]);hold on ;ylim([0 limity])
    hold on ; title('individual velocity plot')
    bin_vel_mean(gr)={mean_vel_binned1};
    figure(fig2)
    shadedErrorBar(1:size(mean_vel_binned1,2),mean_vel_binned1(1,:),mean_vel_binned1(2,:),'lineprops',colour{gr}, 'patchSaturation', 0.2);hold on; plot_boundary(M , [0, limity]);hold on;ylim([0 limity])
    hold on ; plot_boundary(M , [0,limity]) ; hold on; ylim([0 limity]);ylabel('velocity of animals in cm/s');
    if(gr==max(size(X)))
        name=[];name=name_fig1;title(name_fig1)
        saveas(fig1,[filename,'figures\',name,'.fig']);saveas(fig1,[filename,'figures\',name,'.tif']);saveas(fig1,[filename,'figures\',name,'.svg']);
        name=[];name=name_fig2;title(name_fig2)
        saveas(fig2,[filename,'figures\',name,'.fig']);saveas(fig2,[filename,'figures\',name,'.tif']);saveas(fig2,[filename,'figures\',name,'.svg']);
        name=[];name=name_fig3;title(name_fig3)
        saveas(fig3,[filename,'figures\',name,'.fig']);saveas(fig3,[filename,'figures\',name,'.tif']);saveas(fig3,[filename,'figures\',name,'.svg']);
    end
end
name=[];name=name_fig4;figure(fig4);title(name_fig4)
saveas(fig4,[filename,'figures\',name,'.fig']);saveas(fig4,[filename,'figures\',name,'.tif']);saveas(fig4,[filename,'figures\',name,'.svg']);
name=[];name=name_fig5;figure(fig5);title(name_fig5)
saveas(fig5,[filename,'figures\',name,'.fig']);saveas(fig5,[filename,'figures\',name,'.tif']);saveas(fig5,[filename,'figures\',name,'.svg']);
name=[];name='ratio dist in red';figure(fig6);title(name)
saveas(fig6,[filename,'figures\',name,'.fig']);saveas(fig6,[filename,'figures\',name,'.tif']);saveas(fig6,[filename,'figures\',name,'.svg']);
name=[];name='total dist swimmed';figure(fig7);title(name)
saveas(fig7,[filename,'figures\',name,'.fig']);saveas(fig7,[filename,'figures\',name,'.tif']);saveas(fig7,[filename,'figures\',name,'.svg']);
name=[];name=name_fig10;figure(fig10);title(name)
saveas(fig10,[filename,'figures\',name,'.fig']);saveas(fig10,[filename,'figures\',name,'.tif']);saveas(fig10,[filename,'figures\',name,'.svg']);

%%
fig1=figure();
fig2=figure();
limit=110; %size of distance midline from borders
mean_or_median='median';%'mean'
for gr=1:max(size(X))
    X_vector=[];Y_vector=[]; D_vector=[];  info_vector=[]; V_vector=[];
    X_vector=X{1,gr};Y_vector=Y{1,gr}; D_vector=D{1,gr};  info_vector=info{1,gr}; V_vector=V{1,gr};
name_fig1=['distance from midline individuals'];
name_fig2=['distance from midline plot all'];
binned_pos1=[];Flip_session=rule_change;
[binned_pos1]=binned_mean_position_flipping_color(X_vector,Y_vector,D_vector,info_vector,colour{gr},time_bin,base_time,Flip_session,mean_or_median);
bin_dist_from_midline_vector(gr)={binned_pos1};
for p=1:size(binned_pos1,1)
if(p==1)
M=length(binned_pos1{p,1});
else
M= horzcat (M,M(p-1)+length(binned_pos1{p,1}));
end
end
figure(fig1)
subplot(max(size(X)),1,gr)
mean_binned1=[];
mean_binned1=plot_binned(binned_pos1,colour{gr},15); plot_boundary(M , [-limit,limit])
hold on
title ('individual animal plot')
bin_vector_mean_dist_midline(gr)={mean_binned1};
figure(fig2)
shadedErrorBar(1:size(mean_binned1,2),mean_binned1(1,:),mean_binned1(2,:),'lineprops',colour{gr} , 'patchSaturation', 0.2);
hold on ; plot_boundary(M , [-limit,limit]) ; hold on; ylim([-limit,limit]);ylabel('distance from midline');xlabel('time')
if(gr==max(size(X)))
name=[];name=name_fig1;
title(name_fig1)
saveas(fig1,[filename,'figures\',name,'.fig']);saveas(fig1,[filename,'figures\',name,'.tif']);saveas(fig1,[filename,'figures\',name,'.svg']); 
name=[];name=name_fig2;
title(name_fig2)
saveas(fig2,[filename,'figures\',name,'.fig']);saveas(fig2,[filename,'figures\',name,'.tif']);saveas(fig2,[filename,'figures\',name,'.svg']); 
end
end

%%
%%we now calculate the response of the animal to the stimulation only 
days_num=1;%you can perform it for multiple days experiments if that is the case
division_in_time=1;%set the number of windows in which you want to divide your protocol sessions
for windows_to_analyze=1:division_in_time
    fig1=figure();
    clearvars V_raw_all V_raw;
    switch_rule=rule_change;
    timeWindow=80; %number of frames to take into account for the calculation of the average velocity before and after encountering the midline.
    variability_threshold=0;%%number of pixel the annimal needs to be next to the midline to be considered as an encounter
    time_block=windows_to_analyze;
    for gr=1:max(size(X))
        X_vector=[];Y_vector=[]; D_vector=[];  info_vector=[]; V_vector=[];
        X_vector=X{1,gr};Y_vector=Y{1,gr};D_vector=D{1,gr};info_vector=info{1,gr};V_vector=V{1,gr};
        name_fig1=[];name_fig1=['velocity before and after stimulus'];
        
        V1=[];V_raw=[];Peaks=[];Mean_peaks=[];BASE=[];BASE_raw=[];PRE_Peak_response=[];POST_peak_response=[];
        [Angles,Vel,V_raw,Peaks,Mean_peaks,BASE,BASE_raw,PRE_Peak_response,POST_peak_response]=plot_response_to_only_shock_normalized(V_vector,X_vector,Y_vector,info_vector,timeWindow,switch_rule,timeWindow,variability_threshold,division_in_time,time_block);
        V_response_to_shock(1,gr)={PRE_Peak_response};
        V_response_to_shock(2,gr)={POST_peak_response};
        Angles_avoidance(gr,windows_to_analyze)={Angles};
        V_vector1=[];
        BASE(cellfun('isempty',BASE))={(-0.1).*ones(1,timeWindow)}; %a value of -0.1 is set if the animal never encounts the midline
        Vel(cellfun('isempty',V))={(-0.1).*ones(1,timeWindow)};
        V_raw(cellfun('isempty',V_raw))={(-0.1).*ones(1,timeWindow)};
        BASE_raw(cellfun('isempty',BASE_raw))={(-0.1).*ones(1,timeWindow)};
        V_raw_single_fish(gr)={V_raw};
        Base_raw_single_fish(gr)={BASE_raw};
        for temp1=1:size(V_raw,1)
            for  temp2=1:size(V_raw,2)
                if(temp2==1)
                    V_raw_all(gr,temp1)={[V_raw{temp1,temp2}]};
                    BASE_raw_all(gr,temp1)={[BASE_raw{temp1,temp2}]};
                else
                    V_raw_all(gr,temp1)={ vertcat(V_raw_all{gr,temp1},V_raw{temp1,temp2})};
                    BASE_raw_all(gr,temp1)={ vertcat(BASE_raw_all{gr,temp1},BASE_raw{temp1,temp2})};
                end
            end
        end
        figure(fig1)
        subplot(1,max(size(X)),gr)
        TO_PLOT=[];ones_to_plot=[];counters=1;
        for counting=conditioning_session
            TO_PLOT=[TO_PLOT;PRE_Peak_response(counting,:);POST_peak_response(counting,:)];
            ones_to_plot=[ones_to_plot;ones(1,size(PRE_Peak_response(counting,:),2)).*counters];
            counters=counters+1;
            ones_to_plot=[ones_to_plot;ones(1,size(POST_peak_response(counting,:),2)).*counters];
            counters=counters+1;
            Atemp=[];Btemp=[];
            Atemp=PRE_Peak_response(counting,:);
            Btemp=POST_peak_response(counting,:);
            [ h(gr,counting,days_num),p(gr,counting,days_num)] = signrank(Atemp(Atemp~=-0.1),Btemp(Btemp~=-0.1));
        end
        %if you have more than 2 conditionign sessions you need to adjust this plot accordingly to your necessities
        boxplot(TO_PLOT')
        hold on
        scatter (ones_to_plot(:),TO_PLOT(:));
        hold on
        primo_plot=[];
        second_plot=[];
        first_plot=[];
        sec_plot=[];
        primo_plot=ones_to_plot(1:2,:);
        second_plot=ones_to_plot(3:4,:);
        first_plot=TO_PLOT(1:2,:);
        sec_plot=TO_PLOT(3:4,:);
        plot(primo_plot,first_plot)
        hold on
        plot(second_plot,sec_plot)
        hold on
        ylim([-0.2 5])
        title(['session number : ', num2str(conditioning_session), '; ' , name_groups{gr}])
    end
    sgtitle(name_fig1); hold on
    name=[];name=name_fig1;
    saveas(fig1,[filename,name,'.fig']);saveas(fig1,[filename,name,'.tif']);saveas(fig1,[filename,name,'.svg']);
end
%%
%%
%%Now we calculate the response of the animal to the encounter of the
%%boundary of the arena(safe vs danger)
days_num=1;%you can perform it for multiple days experiments if that is the case
division_in_time=6;%set the number of windows in which you want to divide your protocol sessions
for windows_to_analyze=1:division_in_time
    fig1=figure();
    clearvars V_raw_all V_raw;
    switch_rule=rule_change;
    timeWindow=80; %number of frames to take into account for the calculation of the average velocity before and after encountering the midline.
    variability_threshold=0;%%number of pixel the annimal needs to be next to the midline to be considered as an encounter
    time_block=windows_to_analyze;
    for gr=1:max(size(X))
        X_vector=[];Y_vector=[]; D_vector=[];  info_vector=[]; V_vector=[];
        X_vector=X{1,gr};Y_vector=Y{1,gr};D_vector=D{1,gr};info_vector=info{1,gr};V_vector=V{1,gr};
        name_fig1=[];name_fig1=['velocity before and after stimulus'];
        
        V1=[];V_raw=[];Peaks=[];Mean_peaks=[];BASE=[];BASE_raw=[];PRE_Peak_response=[];POST_peak_response=[];
        [Angles,Vel,V_raw,Peaks,Mean_peaks,BASE,BASE_raw,PRE_Peak_response,POST_peak_response]=plot_avoidance_of_shock(V_vector,X_vector,Y_vector,info_vector,timeWindow,switch_rule,timeWindow,variability_threshold,division_in_time,time_block);
        V_response_to_shock(1,gr)={PRE_Peak_response};
        V_response_to_shock(2,gr)={POST_peak_response};
        Angles_avoidance(gr,windows_to_analyze)={Angles};
        V_vector1=[];
        BASE(cellfun('isempty',BASE))={(-0.1).*ones(1,timeWindow)}; %a value of -0.1 is set if the animal never encounts the midline
        Vel(cellfun('isempty',V))={(-0.1).*ones(1,timeWindow)};
        V_raw(cellfun('isempty',V_raw))={(-0.1).*ones(1,timeWindow)};
        BASE_raw(cellfun('isempty',BASE_raw))={(-0.1).*ones(1,timeWindow)};
        V_raw_single_fish(gr)={V_raw};
        Base_raw_single_fish(gr)={BASE_raw};
        for temp1=1:size(V_raw,1)
            for  temp2=1:size(V_raw,2)
                if(temp2==1)
                    V_raw_all(gr,temp1)={[V_raw{temp1,temp2}]};
                    BASE_raw_all(gr,temp1)={[BASE_raw{temp1,temp2}]};
                else
                    V_raw_all(gr,temp1)={ vertcat(V_raw_all{gr,temp1},V_raw{temp1,temp2})};
                    BASE_raw_all(gr,temp1)={ vertcat(BASE_raw_all{gr,temp1},BASE_raw{temp1,temp2})};
                end
            end
        end
        figure(fig1)
        subplot(1,max(size(X)),gr)
        TO_PLOT=[];ones_to_plot=[];counters=1;
        for counting=conditioning_session
            TO_PLOT=[TO_PLOT;PRE_Peak_response(counting,:);POST_peak_response(counting,:)];
            ones_to_plot=[ones_to_plot;ones(1,size(PRE_Peak_response(counting,:),2)).*counters];
            counters=counters+1;
            ones_to_plot=[ones_to_plot;ones(1,size(POST_peak_response(counting,:),2)).*counters];
            counters=counters+1;
            Atemp=[];Btemp=[];
            Atemp=PRE_Peak_response(counting,:);
            Btemp=POST_peak_response(counting,:);
            [ h(gr,counting,days_num),p(gr,counting,days_num)] = signrank(Atemp(Atemp~=-0.1),Btemp(Btemp~=-0.1));
        end
        %if you have more than 2 conditionign sessions you need to adjust this plot accordingly to your necessities
        boxplot(TO_PLOT')
        hold on
        scatter (ones_to_plot(:),TO_PLOT(:));
        hold on
        primo_plot=[];
        second_plot=[];
        first_plot=[];
        sec_plot=[];
        primo_plot=ones_to_plot(1:2,:);
        second_plot=ones_to_plot(3:4,:);
        first_plot=TO_PLOT(1:2,:);
        sec_plot=TO_PLOT(3:4,:);
        plot(primo_plot,first_plot)
        hold on
        plot(second_plot,sec_plot)
        hold on
        ylim([-0.2 5])
        title(['session number : ', num2str(conditioning_session), '; ' , name_groups{gr}])
    end
    sgtitle(name_fig1); hold on
    name=[];name=name_fig1;
    saveas(fig1,[filename,name,'.fig']);saveas(fig1,[filename,name,'.tif']);saveas(fig1,[filename,name,'.svg']);
end
%%
%%blue dots are aproaching angles (radius = mean velocity of approach)
%%red dots are "escaping" angles (radius = mean velocity of approach)
%%The line represent the circular mean value
divide_window=division_in_time;%set the number of windows in which you want to divide your protocol sessions
for session_windowed=1:divide_window
    tt=1;
    for ses=1:block_num
        figA=figure();
        figure(figA)
        for gr=1:max(size(GROUPS))
            ctemp=[];ctemp=Angles_avoidance{gr,session_windowed};
            A=[];    C=[];    T=[];
            num=0;  DIR_vector=[];
            name_fig=[' mean escape strategy in session ',num2str(ses),' window ' ,num2str(session_windowed)];
            subplot(max(size(GROUPS)),1,gr);
            frequency=zeros(size(ctemp,2),2);
            for fishN=1:size(ctemp,2)
                if(isempty(ctemp{ses,fishN}) || isequal( [size(ctemp{ses,fishN})],[1,1]) )
                    if(ctemp{ses,fishN}==1)
                        polarscatter(0.5,0.01,'g','filled');
                        hold on
                        frequency(fishN,2)=1;
                    else
                        frequency(fishN,2)=0;
                        polarscatter(0.5,0.01,'m','filled');
                        hold on
                    end
                    DIR_vector(fishN,1)=NaN;
                    DIR_vector(fishN,2)=NaN;
                    DIR_vector(fishN,3)=NaN;
                    DIR_vector(fishN,4)=NaN;
                    frequency(fishN,1)=0;
                    continue
                end
                C=ctemp{ses,fishN};
                A=((C(:,1)));
                A= A(A~=361);
                Apost=((C(:,7)));
                Apost= Apost(Apost~=361);
                A=deg2rad(A);
                Apost=deg2rad(Apost);
                T=C(:,6);
                Vpre=(mean((T(T~=-0.1))));
                T=[];T=C(:,12);
                Vpost=(mean((T(T~=-0.1))));
                entrance_after=[];entrance_after=C(:,9);
                entrance_pre=[];entrance_pre=C(:,8);
                if(isempty(entrance_pre))
                    frequency(fishN,1)=0;
                    frequency(fishN,2)=ctemp{ses,fishN};
                else
                    frequency(fishN,1)=max(size(entrance_after));
                    if (strcmp(rule_change,'NO'))
                        frequency(fishN,2)=(max(size((entrance_pre(entrance_after>=entrance_pre)))))./max(size(entrance_after));
                    else
                        if (ses<rule_change)
                            frequency(fishN,2)=(max(size((entrance_pre(entrance_after>=entrance_pre)))))./max(size(entrance_after));
                        else
                            frequency(fishN,2)=(max(size((entrance_pre(entrance_after<=entrance_pre)))))./max(size(entrance_after));
                        end
                    end
                end
                if(isempty(Apost) || isnan(Vpost)|| isempty(A) || isnan(Vpre) )
                    polarscatter(0.3,0.01,'k','filled');
                    hold on
                    DIR_vector(fishN,1)=NaN;
                    DIR_vector(fishN,2)=NaN;
                    DIR_vector(fishN,3)=NaN;
                    DIR_vector(fishN,4)=NaN;
                    continue
                else
                    DIR_vector(fishN,3)=circ_mean((Apost));
                    DIR_vector(fishN,4)=Vpost;
                    DIR_vector(fishN,1)=circ_mean((A));
                    DIR_vector(fishN,2)=Vpre;
                    polarscatter(circ_mean((A)),Vpre,'b');
                    hold on
                    polarscatter(circ_mean((Apost)),Vpost,'r');hold on
                end
                hold on
                Ax = gca; % current axes
                Ax.ThetaTickLabel = [];
            end
            Angles_avoidance_frequency(gr,ses,session_windowed)={frequency};
            DIR_vector( find(isnan(DIR_vector)))=[];
            DIR_vector=reshape(DIR_vector,length(DIR_vector)/4,4);
            figure(figA)
            subplot(max(size(GROUPS)),1,gr);
            polarplot([0,circ_mean( [DIR_vector(:,1)])],[0,(mean(DIR_vector(:,2)))],'b','linewidth',3);
            hold on
            polarplot([0,circ_mean( DIR_vector(:,3))],[0,(mean(DIR_vector(:,4)))],'r','linewidth',3);
            hold on
            hold on
            rlim ([0 max(max(DIR_vector(:,[2,4])))]);
            title( [ name_fig,' ',name_groups{gr}]);
            hold on
        end
        name=[];name=name_fig;
        saveas(figA,[filename,name,'.fig']);saveas(figA,[filename,name,'.tif']);saveas(figA,[filename,name,'.svg']);
    end
end
close(figA)
fig_avoidance=figure();
beseline_reference=1;%session of reference in the baseline to run the statistics
for gr= 1: size(Angles_avoidance_frequency,1)
    FREQ_tot=[];
    for ses =1:size(Angles_avoidance_frequency,2)
        for window = 1:size(Angles_avoidance_frequency,3)
            temp=[];temp=Angles_avoidance_frequency{gr,ses,window};
            freq_avoidance_temp=[];freq_avoidance_temp=temp(:,2);
            if(ses==1 && window==1)
                FREQ_tot=[];
                FREQ_tot=freq_avoidance_temp;
            else
                FREQ_tot=horzcat(FREQ_tot,freq_avoidance_temp);
            end
        end
    end
    shadedErrorBar(1:size(FREQ_tot,2),mean(FREQ_tot,1),std(FREQ_tot,0,1)./sqrt(size(FREQ_tot,1)),'lineprops',colour{gr} , 'patchSaturation', 0.2);
    hold on
    FREQUENCY_GROUPS(gr,days_num)={FREQ_tot};
    base_freq=FREQ_tot(:,beseline_reference);
    for ses_point=7:size(FREQ_tot,2)
        [h_success(gr,ses_point,days_num),p_succes(gr,ses_point,days_num)]=signrank(base_freq,FREQ_tot(:,ses_point));
    end
    GRUPS_AVOIDANCE(gr)={FREQ_tot};
end
name=[];name=['successful avoidance per session ', num2str(variability_threshold)];title(name)
saveas(fig_avoidance,[filename,name,'.fig']);saveas(fig_avoidance,[filename,name,'.tif']);saveas(fig_avoidance,[filename,name,'.svg']);
fig_tot=figure()
for gr= 1: size(Angles_avoidance_frequency,1)
    FREQ_abs=[];
    for ses =1:size(Angles_avoidance_frequency,2)
        for window = 1:size(Angles_avoidance_frequency,3)
            temp=[];temp=Angles_avoidance_frequency{gr,ses,window};
            freq_avoidance_temp=[];freq_avoidance_temp=temp(:,1);
            if(ses==1 && window==1)
                FREQ_abs=freq_avoidance_temp;
            else
                FREQ_abs=horzcat(FREQ_abs,freq_avoidance_temp);
            end
        end  
    end
    shadedErrorBar(1:size(FREQ_abs,2),mean(FREQ_abs,1),std(FREQ_abs,0,1)./sqrt(size(FREQ_abs,1)),'lineprops',colour{gr} , 'patchSaturation', 0.2);
    hold on
    FREQUENCY_TOTAL_APPROACH_GROUPS(gr,days_num)={FREQ_abs};
    base_freq_abs=FREQ_abs(:,beseline_reference);
    for ses_point=7:size(FREQ_abs,2)
        [h_success_abs(gr,ses_point,days_num),p_succes_abs(gr,ses_point,days_num)]=signrank(base_freq_abs,FREQ_abs(:,ses_point));
    end
    GRUPS_TOT_AVOIDANCE(gr)={FREQ_abs};
end
name=[];name='total approach per session ';title(name)
saveas(fig_tot,[filename,name,'.fig']);saveas(fig_tot,[filename,name,'.tif']);saveas(fig_tot,[filename,name,'.svg']);


%%

%%this are all the functions you need to run the above code
function [ ITI,ITImean_std]=mean_duration_entrance(X,D,info,flipping,V)
for i=1:size(X,1)
    for j=1:size(X,2)
        x=[];x=X{i,j};d=[];d=cumsum(D{i,j});
        v=[];v=V{i,j};
        if (strcmp('NO',flipping))
            border=info{3,j}
            x(x<=border)=1;
            x(x>border)=0;
            xdiff=[];xdiff=diff(x);
            ENTRANCE=[];ENTRANCE=find(xdiff==1);
            EXIT=[];EXIT=find(xdiff==-1);
            if (max(size(EXIT))>1)
                if(d(EXIT(1))-d(ENTRANCE(1))<0)
                    EXIT(1)=[];
                end
                if(d(EXIT(end))-d(ENTRANCE(end))<0)
                    ENTRANCE(end)=[];
                end
                vel_tot=[];vel_mean=[];
                for t=1:size(ENTRANCE)
                    vel_tot=[vel_tot,v(ENTRANCE(t):EXIT(t))];
                    vel_mean(t,1)=nanmean(v(ENTRANCE(t):EXIT(t)));
                end
                ITImean_std(i,j,4)=nanmean(vel_tot(:));
                ITI(i,j)={[d(EXIT)'-d(ENTRANCE)']./1000};
                ITImean_std(i,j,1)=nanmean(d(EXIT)'-d(ENTRANCE)')./1000;
                ITImean_std(i,j,2)=nanmean(((d(EXIT)'-d(ENTRANCE)')./1000).*vel_mean);
                ITImean_std(i,j,3)=size([d(EXIT)'-d(ENTRANCE)'],1);
            else
                ITI(i,j)={0};
                ITImean_std(i,j,1)=0;
                ITImean_std(i,j,2)=0;
                ITImean_std(i,j,3)=0;
                ITImean_std(i,j,4)=0;
            end
        else
            if(i<flipping)
                border=info{3,j}
                x(x<=border)=1;
                x(x>border)=0;
            else
                border=info{4,j}
                x(x<border)=0;
                x(x>=border)=1;
            end
            xdiff=[];xdiff=diff(x);
            ENTRANCE=[];ENTRANCE=find(xdiff==1);
            EXIT=[];EXIT=find(xdiff==-1);
            if (max(size(EXIT))>1)
                if(d(EXIT(1))-d(ENTRANCE(1))<0)
                    EXIT(1)=[];
                end
                if(d(EXIT(end))-d(ENTRANCE(end))<0)
                    ENTRANCE(end)=[];
                end
                vel_tot=[];vel_mean=[];
                for t=1:size(ENTRANCE)
                    vel_tot=[A,v(ENTRANCE(t):EXIT(t))];
                    vel_mean(t,1)=nanmean(v(ENTRANCE(t):EXIT(t)));
                end
                ITImean_std(i,j,4)=nanmean(vel_tot(:));
                ITI(i,j)={[d(EXIT)'-d(ENTRANCE)']./1000};
                ITImean_std(i,j,1)=nanmean(d(EXIT)'-d(ENTRANCE)')./1000;
                ITImean_std(i,j,2)=nanmean(((d(EXIT)'-d(ENTRANCE)')./1000).*vel_mean);
                ITImean_std(i,j,3)=size([d(EXIT)'-d(ENTRANCE)'],1);
            else
                ITI(i,j)={0};
                ITImean_std(i,j,1)=0;
                ITImean_std(i,j,2)=0;
                ITImean_std(i,j,3)=0;
                ITImean_std(i,j,4)=0;
            end
        end
    end
end
end
function [index]=LI_with_inversion_of_the_rule_baseline_refered(X_vector,Y_vector,D_vector,info_vector,change_rule)
block=size(X_vector,1);
fish_n=size(X_vector,2);
learning_index=zeros(block,fish_n);
for block_n=1:block
    for i=1: fish_n
        
        if(size(info_vector,1)==4)
            if(block_n <change_rule)
                threshold=info_vector{3,i};
                threshold=cell2mat(info_vector(3,i));
            else
                if(isempty(info_vector{4,i}))
                    threshold=info_vector{3,i};
                    threshold=cell2mat(info_vector(3,i));
                else
                    threshold=info_vector{4,i};
                    threshold=cell2mat(info_vector(4,i));
                end
            end
        else
            threshold=info_vector{3,i};
            threshold=cell2mat(info_vector(3,i));
        end
        X=X_vector{block_n,i};
        Y=Y_vector{block_n,i};
        D=D_vector{block_n,i};
        [timespenttemp,timetotal]=timespent(X,D,threshold);
        if(strcmp(change_rule,'NO'))
            if(block_n==1 )
                timebasetemp(i)=((timespenttemp/timetotal));
                index(block_n,i)=((timespenttemp/timetotal));
            else
                index(block_n,i)=(timebasetemp(i)-((timespenttemp/timetotal)))/timebasetemp(i);
            end
            
        else
            if((block_n==1)) %
                timebasetemp(i)=((timespenttemp/timetotal));
                index(block_n,i)=((timespenttemp/timetotal));
            end
            if(block_n>1)
                if(block_n<change_rule)
                    index(block_n,i)=(timebasetemp(i)-((timespenttemp/timetotal)))/timebasetemp(i);
                else
                    index(block_n,i)=((1-timebasetemp(i))-(((timetotal-(timespenttemp))/timetotal)))/(1-timebasetemp(i));
                    
                end
            end
        end
    end
end
index(index<-1)=-1;
clearvars X D Y threshold
end
function [Tspent, Ttotal]=timespent(X,D,threshold)
timespent2=0;
timetotal2=0;
for j=2:length(X)
    if(X(j)<threshold)
        timespent2 = timespent2 + D(j);
    end
    timetotal2= timetotal2 + D(j);
end
Tspent=timespent2;
Ttotal=timetotal2;
end
function [xnew,ynew,dnew,pdf,vnew]=get_time_bins(X,Y,D,V,time_bin,size_arena,min_size_bin)%%timebin in min
xnew=[];ynew=[];dnew=[];pdf=[];vnew=[];
for i=1:size(X,2)
    size_arena2=max(size_arena{1:2,i});
    xnew2=cell([0]);ynew2=cell([0]);dnew2=cell([0]);pdf2=cell([0]);vnew2=cell([0]);
    index2=1;
    for j=1:(size(X,1))
        x=X{j,i};
        y=Y{j,i};
        v=V{j,i};
        if(size(x,2)==1)
            v=v';
        end
        d=D{j,i};
        if(size(x,2)==1)
            x=x';
        end
        if(size(y,2)==1)
            y=y';
        end
        if(size(d,2)==1)
            d=d';
        end
        dtemp=[];
        dtemp=diff(floor(cumsum(d)./(time_bin*60000)));
        dtemp2=[];dtemp2=[1,find(dtemp)];
        for index=1:size(dtemp2,2)
            if(index<size(dtemp2,2))
                xnew2(index2)={x(dtemp2(index):dtemp2(index+1))};
                ynew2(index2)={y(dtemp2(index):dtemp2(index+1))};
                vnew2(index2)={v(dtemp2(index):dtemp2(index+1))};
                dnew2(index2)={d(dtemp2(index):dtemp2(index+1))};
                pdf2(index2)={get_PDF(x(dtemp2(index):dtemp2(index+1))',y(dtemp2(index):dtemp2(index+1))',size_arena2)};
                index2=index2+1;
            else
                if(max(size(x(dtemp2(index):end)))>(min_size_bin))
                    xnew2(index2)={x(dtemp2(index):end)};
                    ynew2(index2)={y(dtemp2(index):end)};
                    dnew2(index2)={d(dtemp2(index):end)};
                    vnew2(index2)={v(dtemp2(index):end)};
                    pdf2(index2)={get_PDF(x(dtemp2(index):end)',y(dtemp2(index):end)',size_arena2)};
                    index2=index2+1;
                end
            end
        end
        clearvars x d y v
    end
    if(size(xnew,1)>size(xnew2',1) && size(xnew,1)>0 )
        xnew=[xnew,[xnew2';NaN(size(xnew,1)-size(xnew2',1))]];
        ynew=[ynew,[ynew2';NaN(size(xnew,1)-size(xnew2',1))]];
        dnew=[dnew,[dnew2';NaN(size(xnew,1)-size(xnew2',1))]];
        pdf=[pdf,[pdf2';NaN(size(xnew,1)-size(xnew2',1))]];
        vnew=[vnew,[vnew2';NaN(size(xnew,1)-size(xnew2',1))]];
    else
        xnew=[xnew,xnew2'];ynew=[ynew,ynew2'];dnew=[dnew,dnew2'];pdf=[pdf,pdf2'];vnew=[vnew,vnew2'];
    end
end
end
function varargout=shadedErrorBar(x,y,errBar,varargin)
% generatcontinuous error bar area around a line plot
%
% function H=shadedErrorBar(x,y,errBar, ...)
%
% Purpose 
% Makes a 2-d line plot with a pretty shaded error bar made
% using patch. Error bar color is chosen automatically.
%
%
% Inputs (required)
% x - vector of x values [optional, can be left empty]
% y - vector of y values or a matrix of n observations by m cases
%     where m has length(x);
% errBar - if a vector we draw symmetric errorbars. If it has a size
%          of [2,length(x)] then we draw asymmetric error bars with
%          row 1 being the upper bar and row 2 being the lower bar
%          (with respect to y). ** alternatively ** errBar can be a
%          cellArray of two function handles. The first defines which
%          statistic the line should be and the second defines the
%          error bar.
%
% Inputs (optional, param/value pairs)
% 'lineProps' - ['-k' by default] defines the properties of
%             the data line. e.g.:    
%             'or-', or {'-or','markerfacecolor',[1,0.2,0.2]}
% 'transparent' - [true  by default] if true, the shaded error
%               bar is made transparent. However, for a transparent
%               vector image you will need to save as PDF, not EPS,
%               and set the figure renderer to "painters". An EPS 
%               will only be transparent if you set the renderer 
%               to OpenGL, however this makes a raster image.
% 'patchSaturation'- [0.2 by default] The saturation of the patch color.
%
%
%
% Outputs
% H - a structure of handles to the generated plot objects.
%
%
% Examples:
% y=randn(30,80); 
% x=1:size(y,2);
%
% 1)
% shadedErrorBar(x,mean(y,1),std(y),'lineprops','g');
%
% 2)
% shadedErrorBar(x,y,{@median,@std},'lineprops',{'r-o','markerfacecolor','r'});
%
% 3)
% shadedErrorBar([],y,{@median,@(x) std(x)*1.96},'lineprops',{'r-o','markerfacecolor','k'});
%
% 4)
% Overlay two transparent lines:
% clf
% y=randn(30,80)*10; 
% x=(1:size(y,2))-40;
% shadedErrorBar(x,y,{@mean,@std},'lineprops','-r','transparent',1);
% hold on
% y=ones(30,1)*x; y=y+0.06*y.^2+randn(size(y))*10;
% shadedErrorBar(x,y,{@mean,@std},'lineprops','-b','transparent',1);
% hold off
%
%
% Rob Campbell - November 2009



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse input arguments
narginchk(3,inf)

params = inputParser;
params.CaseSensitive = false;
params.addParameter('lineProps', '-k', @(x) ischar(x) | iscell(x));
params.addParameter('transparent', true, @(x) islogical(x) || x==0 || x==1);
params.addParameter('patchSaturation', 0.2, @(x) isnumeric(x) && x>=0 && x<=1);

params.parse(varargin{:});

%Extract values from the inputParser
lineProps =  params.Results.lineProps;
transparent =  params.Results.transparent;
patchSaturation = params.Results.patchSaturation;

if ~iscell(lineProps), lineProps={lineProps}; end


%Process y using function handles if needed to make the error bar dynamically
if iscell(errBar) 
    fun1=errBar{1};
    fun2=errBar{2};
    errBar=fun2(y);
    y=fun1(y);
else
    y=y(:).';
end

if isempty(x)
    x=1:length(y);
else
    x=x(:).';
end


%Make upper and lower error bars if only one was specified
if length(errBar)==length(errBar(:))
    errBar=repmat(errBar(:)',2,1);
else
    s=size(errBar);
    f=find(s==2);
    if isempty(f), error('errBar has the wrong size'), end
    if f==2, errBar=errBar'; end
end

if length(x) ~= length(errBar)
    error('length(x) must equal length(errBar)')
end


%Log the hold status so we don't change
initialHoldStatus=ishold;
if ~initialHoldStatus, hold on,  end

H = makePlot(x,y,errBar,lineProps,transparent,patchSaturation);

if ~initialHoldStatus, hold off, end

if nargout==1
    varargout{1}=H;
end



function H = makePlot(x,y,errBar,lineProps,transparent,patchSaturation)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot to get the parameters of the line
    H.mainLine=plot(x,y,lineProps{:},'linewidth',2);


    % Work out the color of the shaded region and associated lines.
    % Here we have the option of choosing alpha or a de-saturated
    % solid colour for the patch surface.
    mainLineColor=get(H.mainLine,'color');
    edgeColor=mainLineColor+(1-mainLineColor)*0.55;

    if transparent
        faceAlpha=patchSaturation;
        patchColor=mainLineColor;
    else
        faceAlpha=1;
        patchColor=mainLineColor+(1-mainLineColor)*(1-patchSaturation);
    end


    %Calculate the error bars
    uE=y-errBar(1,:);
    lE=y+errBar(2,:);


    %Add the patch error bar



    %Make the patch
    yP=[lE,fliplr(uE)];
    xP=[x,fliplr(x)];

    %remove nans otherwise patch won't work
    xP(isnan(yP))=[];
    yP(isnan(yP))=[];


    H.patch=patch(xP,yP,1,'facecolor',patchColor, ...
                  'edgecolor','none', ...
                  'facealpha',faceAlpha);


    %Make pretty edges around the patch. 
    H.edge(1)=plot(x,lE,'-','color',edgeColor);
    H.edge(2)=plot(x,uE,'-','color',edgeColor);



    uistack(H.mainLine,'top') % Bring the main line to the top
end


end
function [pdf]=get_PDF(X,Y,dim_max)
x= 0 : dim_max/25 : dim_max; % axis x, which you want to see
y = 0 :dim_max/25: dim_max; % axis y, which you want to see
[Xtemp,Ytemp] = meshgrid(x,y); % important for "surf" - makes defined grid
pdf = hist3([X , Y],{x y}); % standard hist3 (calculated for yours axis
pdf_normalize = (((pdf./ length(X)))); % normalization means devide it by length of
pdf_normalize =(pdf_normalize)';% data_x (or data_y)
pdf=pdf_normalize;
end
function [PDF_MEAN] = plot_heat_maps_averaged(PDF)
x=size(PDF{1,1},1);
y=size(PDF{1,1},2);
for i=  1:size(PDF,1)
    for j = 1:size(PDF,2)   
        pdf(1:x,1:y,j)=PDF{i,j}
    end 
    PDF_MEAN(i,1)={mean(pdf(1:end,1:end,:),3)}; 
end
end

function []=plot_heat_map(pdf_mean,y,x,z)

for i=1:size(pdf_mean,1)
subplot(z,x,i+(x*(y-1)))
imagesc(pdf_mean{i,1})
colormap jet 
caxis ([0 0.005])
sum(sum(pdf_mean{i,1}))
axis off
end

end
function [Binned_time]=new_binned_plot_with_flipping_color(X,Y,D,info_vector,color,time_bin,base_time,Flip_time)
for i=1:size(X,2)
    max_size=max(info_vector{3,i},info_vector{2,i});
    for j=1:(size(X,1))
        x=X{j,i};
        y=Y{j,i};
        d=D{j,i};
%         [x_3,y_3]=remove_pos_spike(x, y); you can implement here a
%         function to filter your data if a lot of detection errors are
%         present
%         x=[];  y=[];
%         x=x_3;
%         y=y_3;
        if(size(x,2)==1)
            x=x';
        end
        if(size(y,2)==1)
            y=y';
        end
        if(size(d,2)==1)
            d=d';
        end
        if(size(info_vector,1)==4)
            if( strcmp(Flip_time,'NO'))
                thresh=info_vector(3,:);
                border=thresh{i};
            else
                if(j<Flip_time  )
                    thresh=info_vector(3,:);
                    border=thresh{i};
                else
                    if(isempty(info_vector{4,1}))
                        thresh=info_vector(3,:);
                        border=thresh{i};
                    else
                        thresh=info_vector(4,:);
                        border=thresh{i};
                    end
                end
            end
        else
            thresh=info_vector(3,:);
            border=thresh{i};
        end
        if(j==1)
            time=0;
            for k=max(size(x)):-1:1
                
                if(time<base_time && k>1)
                    time=time+d(k);
                else
                    Binned_time(j,i)={divide_time_in_red(x(k:end),d(k:end),border,time_bin)};
                    break
                end
            end
            clearvars k
        else
            if (strcmp(Flip_time,'NO'))
                Binned_time(j,i)={divide_time_in_red(x,d,border,time_bin)};
            else
                if(j<Flip_time )
                    Binned_time(j,i)={divide_time_in_red(x,d,border,time_bin)};
                else
                    Binned_time(j,i)={divide_time_in_gray(x,d,border,time_bin)};
                end
            end
            clearvars x time  d y
        end
    end
    clearvars time timetot point percentage_time
end
end
function [] = plot_boundary(M, t)
for s = 1:size(M,2)-1
    plot([M(s),M(s)], t, 'k','linewidth', 3);
    hold on
end
hold on
end
function [Mean_std]=plot_binned(Bin,color,binning)
Xtot=NaN(size(Bin,2),size(Bin,1)*binning);
for j=1:size(Bin,2)
    X=[];
    for i=1:size(Bin,1)
        x=NaN(1,binning);
        x=Bin{i,j};
        X=[X,x];
    end
    plot(X)
    hold on
    size(X)
    X_tot(j,1:size(X,2))=X;
end
errorbar(1:size(X_tot,2),nanmean(X_tot,1),nanstd(X_tot,1)./sqrt(size(X_tot,1)),color ,'linewidth', 2 )
xlim([0,size(X_tot,2)+1]);
hold on
Mean_std(1,1:size(X_tot,2))=nanmean(X_tot,1);
Mean_std(2,1:size(X_tot,2))=nanstd(X_tot,1)./sqrt(size(X_tot,1));
end
function [Dist_raw,Distance_sum,Velocity,binned_vel,D_red]=binned_plot_velocity(X,Y,D,info_vector,Normalized,Flip_time,time_bin, base_time,V_vector)
for i=1:size(X,2)
    max_size=max(info_vector{3,i},info_vector{2,i});
    for j=1:(size(X,1))
        x=X{j,i};        y=Y{j,i};        d=D{j,i};
%         [x_3,y_3]=remove_pos_spike(x, y);if you want to filter your data
%         because of tracking errors implement it here
%         x=[];  y=[];    x=x_3;y=y_3;
        if(size(x,2)==1)
            x=x';
        end
        if(size(y,2)==1)
            y=y';
        end
        if(size(d,2)==1)
            d=d';
        end
        if(size(info_vector,1)==4)
           if( strcmp(Flip_time,'NO'))
                thresh=info_vector(3,:);
                border=thresh{i};
            else
                if(j<Flip_time  )
                    thresh=info_vector(3,:);
                    border=thresh{i};
                else
                    if(isempty(info_vector{4,i}))
                        thresh=info_vector(3,:);
                        border=thresh{i};
                    else
                        thresh=info_vector(4,:);
                        border=thresh{i};
                    end
                end
            end
        else
            thresh=info_vector(3,:);
            border=thresh{i};
        end
        if(j==1)
            time=0;
            for k=max(size(x)):-1:1
                if(time<base_time && k>1)
                    time=time+d(k);
                else
                    break
                end
            end
            Distance(j,i)= {(((x(k+1:end)-x(k:end-1)).^2)+((y(k+1:end)-y(k:end-1)).^2)).^(0.5)};
            Dist_raw(j,i)={Distance{j,i}.*(120/max_size)};
            Velocity=V_vector;
            if (Normalized)
                scale_factor=mean(Velocity{1,i});
            else
                scale_factor=1;
            end
            Velocity(j,i)={Velocity{j,i}./scale_factor};
            v_temporary=[];v_temporary=Velocity{j,i};
            binned_vel(j,i)={divide_velocity_bin(v_temporary(k+1:end),d(k+1:end),time_bin)};
            binned_vel(j,i)={binned_vel{j,i}};
            temp_dist_calc=[];
            temp_dist_calc= Distance{j,i};
            D_red(j,i)=sum(temp_dist_calc(x(k+1:end)<border)).*(0.12/max_size);
            Distance_sum(j,i)=sum(Distance{j,i}).*(0.12/max_size);
            clearvars k
        else
            Distance(j,i)= {(((x(2:end)-x(1:end-1)).^2)+((y(2:end)-y(1:end-1)).^2)).^(0.5)};
            Dist_raw(j,i)={Distance{j,i}.*(120/max_size)};
            Velocity(j,i)={Velocity{j,i}./scale_factor};
            v_temporary=[];v_temporary=Velocity{j,i};
            binned_vel(j,i)={divide_velocity_bin(v_temporary,d(2:end),time_bin)};
            binned_vel(j,i)={binned_vel{j,i}};
            if( strcmp(Flip_time,'NO'))
                temp_dist_calc=[];
                temp_dist_calc= Distance{j,i};
                D_red(j,i)=sum(temp_dist_calc(x(1:end-1)<border)).*(0.12/max_size);
            else
                if(j<Flip_time  )
                    temp_dist_calc=[];
                    temp_dist_calc= Distance{j,i};
                    D_red(j,i)=sum(temp_dist_calc(x(1:end-1)<border)).*(0.12/max_size);
                else
                    temp_dist_calc=[];
                    temp_dist_calc= Distance{j,i};
                    D_red(j,i)=sum(temp_dist_calc(x(1:end-1)>border)).*(0.12/max_size);
                end
            end
            Distance_sum(j,i)=sum(Distance{j,i}).*(0.12/max_size);
            clearvars x time  d y
        end
    end
    clearvars time timetot point percentage_time
end

end
function [velocity_bin]=divide_velocity_bin(V,D,time_bin)
time_tot=D(1);
point=1;
counting=1;
while(time_tot>=time_bin)
    velocity_bin(point)= 0;
    time_tot=abs(time_tot-time_bin);
    point=point+1;
end
for n=1:(length (V))
    time_tot=time_tot+D(n);
    if(time_tot<time_bin)
        if(n==length(V) && time_tot > time_bin/20 )
            velocity_bin(point)= mean(V(counting:n));         
            time_tot=0;
            point=point+1;
            counting=n;
        end
    else
        velocity_bin(point)= mean(V(counting:n));
        time_tot/60000
        time_tot=0;
        point=point+1;
        counting=n;
    end  
end
end
function [freezing_percentage,freezing_time]=freezing_new(X,Y,D,parameters,info_vector) % dist to freeze in mm , sec to freeze in sec
dist_to_freeze=parameters(2);
second_to_freeze=parameters(1);
for j=1:size(X,1)
    for i=1:size(X,2)
        tempd=[];tempd= D{j,i};
        tempx= [];
        tempy= [];
        tempx= X{j,i};
        tempy= Y{j,i};
        deltax=[];
        deltay=[];
        dist_tot=[];
        deltax=diff(tempx);
        deltay=diff(tempy);
        dist_tot=sqrt((deltax.^2 )+(deltay.^2 ))./(120/min(info_vector{1:2,i}));
        dist_cumulative=floor(cumsum(dist_tot)./dist_to_freeze);
        dist_cum=diff(dist_cumulative);
        moved=[];moved=[find(dist_cum>0)];
        deltat=[];
        deltat(1)=sum(tempd(1:moved(1))) ;
        for t= 1: max(size(moved))
            if(t< max(size(moved)))
                deltat(t)=sum(tempd((moved(t)+1):moved(t+1))) ;
            else
                deltat(t)=sum(tempd((moved(t)+1):moved(end))) ;
            end
        end
        freezing_time(j,i)=sum(deltat(deltat>=(second_to_freeze*1000)));
        freezing_percentage(j,i)=freezing_time(j,i)./sum(deltat);
    end
end

end
function [percentage_time]=divide_time_in_red(X,D,border,time_bin) %%calculate the time on the left side of the arena
time_tot=0;
time=0;
point=1;
sum(D)/60000
for n=1:(length (X))
    if(time_tot<=time_bin)
        time_tot=time_tot+D(n);
        if(X(n)<border)
            time=time+D(n);
        end
        if(n==length(X) && time_tot > 60000 )
            percentage_time(point)= time/time_tot;
            time=0;
            time_tot=0;
            point=point+1;
            continue
        else
            if(n==length(X) && time_tot <= 60000 && time_tot>30000)
                percentage_time(point)=percentage_time(point-1);
                time=0;
                time_tot=0;
                point=point+1;
                continue
            end
        end
    else
        percentage_time(point)= time/time_tot;
        time_tot/60000
        time=0;
        time_tot=0;
        point=point+1;
    end
end
end
function [binned_dist_red,binned_dist_all] = dist_in_red_time_bin(X,Y,D,info_vector,Normalized,Flip_time,time_bin,base_time,Normalise_to_base,czcs)

for i=1:size(X,2)
    
    max_size=max(info_vector{3,i},info_vector{2,i});
    for j=1:(size(X,1))
        
        x=X{j,i};
        y=Y{j,i};
        d=D{j,i};
%         [x_3,y_3]=remove_pos_spike(x, y); to filter your data
%         x=[];  y=[];
%         
%         x=x_3;
%         y=y_3;
        if(size(x,2)==1)
            x=x';
        end
        if(size(y,2)==1)
            y=y';
        end
        if(size(d,2)==1)
            d=d';
        end
        if(size(info_vector,1)==4)
            if( strcmp(Flip_time,'NO'))
                thresh=info_vector(3,:);
                border=thresh{i};
            else
                if(j<Flip_time  )
                    thresh=info_vector(3,:);
                    border=thresh{i};
                else
                    if(isempty(info_vector{4,1}))
                        thresh=info_vector(3,:);
                        border=thresh{i};
                    else
                        thresh=info_vector(4,:);
                        border=thresh{i};
                    end
                end
            end
        else
            thresh=info_vector(3,:);
            border=thresh{i};
        end
        if(j==1)
            time=0;
            for k=max(size(x)):-1:1
                if(time<base_time && k>1)
                    time=time+d(k);
                else
                    break
                end
            end
            Distance(j,i)= {(((x(k+1:end)-x(k:end-1)).^2)+((y(k+1:end)-y(k:end-1)).^2)).^(0.5)};
            Dist_raw(j,i)={Distance{j,i}.*(12/max_size)};
            [binned_dist_red(j,i),binned_dist_all(j,i)]=divide_distance_bin( d(k+1:end),Dist_raw{j,i},x(k+1:end),time_bin,border,1,czcs);
            if(Normalise_to_base)
                normalise_factor(i)=mean(binned_dist_red{j,i});
            else
                normalise_factor(i)=1;
            end
            binned_dist_red(j,i)={ binned_dist_red{j,i}./normalise_factor(i)};
            clearvars k
        else
            Distance(j,i)= {(((x(2:end)-x(1:end-1)).^2)+((y(2:end)-y(1:end-1)).^2)).^(0.5)};
            Dist_raw(j,i)={Distance{j,i}.*(12/max_size)};
            if( strcmp(Flip_time,'NO'))
                [binned_dist_red(j,i),binned_dist_all(j,i)]=divide_distance_bin(d,Dist_raw{j,i},x(2:end),time_bin,border,1,czcs);
                binned_dist_red(j,i)={ binned_dist_red{j,i}./normalise_factor(i)};
            else
                if(j<Flip_time  )
                    [binned_dist_red(j,i),binned_dist_all(j,i)]=divide_distance_bin(d,Dist_raw{j,i},x(2:end),time_bin,border,1,czcs);
                    TT=[];TT= binned_dist_red(j,i);
                    binned_dist_red(j,i)={ binned_dist_red{j,i}./normalise_factor(i)};
                    if(j==Flip_time-1 && Normalise_to_base)
                        normalise_factor(i)=mean(TT{:}(end-5:end));
                    end
                else
                    [binned_dist_red(j,i),binned_dist_all(j,i)]=divide_distance_bin(d,Dist_raw{j,i},x(2:end),time_bin,border,2,czcs);
                    binned_dist_red(j,i)={ binned_dist_red{j,i}./normalise_factor(i)};
                end
            end
            clearvars x time  d y
        end
    end
    clearvars time timetot point percentage_time
end
end
function [Mean_std,X_tot]=plot_binned_dist(Bin,color,binning)
Xtot=NaN(size(Bin,2),size(Bin,1)*binning);
for j=1:size(Bin,2)
    X=[];
    for i=1:size(Bin,1)
        x=NaN(1,binning);
        x=Bin{i,j};       
     X=[X,x];       
    end
    plot(X)
    hold on 
    size(X)
    X_tot(j,1:size(X,2))=X;
end
errorbar(1:size(X_tot,2),nanmean(X_tot,1),nanstd(X_tot,1)./sqrt(size(X_tot,1)),color ,'linewidth', 2 ) 
xlim([0,size(X_tot,2)+1]);
hold on 
Mean_std(1,1:size(X_tot,2))=nanmean(X_tot,1);
Mean_std(2,1:size(X_tot,2))=nanstd(X_tot,1)./sqrt(size(X_tot,1));
end
function [dist_R,dist_A]=divide_distance_bin(time,D,X,time_bin,threshold,rule,czcs)% if rule is 1 x<threshold else x>threshold
time_tot=0;
point=1;
counting=1;
for n=1:(length (D))
    
    
    
    if(time_tot<time_bin)
        time_tot=time_tot+time(n);
        
        if(n==length(D) && time_tot >(0.5*time_bin ))
            if(rule==1)
                tempD=[];tempD=D(counting:n);
                dist_red(point)= sum(tempD(X(counting:n)<threshold));
                dist_all(point)= sum(D(counting:n));
                
                time_tot=0;
                point=point+1;
                counting=n;
            else
                tempD=[];tempD=D(counting:n);
                dist_red(point)= sum(tempD(X(counting:n)>threshold));
                dist_all(point)= sum(D(counting:n));
                
                time_tot=0;
                point=point+1;
                counting=n;
            end
        end
    else
        if(rule==1)
            tempD=[];tempD=D(counting:n);
            dist_red(point)= sum(tempD(X(counting:n)<threshold));
            dist_all(point)= sum(D(counting:n));
            time_tot/60000
            time_tot=0;
            point=point+1;
            counting=n;
        else
            tempD=[];tempD=D(counting:n);
            dist_red(point)= sum(tempD(X(counting:n)>threshold));
            dist_all(point)= sum(D(counting:n));
            time_tot/60000
            time_tot=0;
            point=point+1;
            counting=n;
        end
    end
end
dist_R= {dist_red};
temp=[];
if (czcs)
    temp=dist_red./(dist_all-dist_red);
    temp(temp>5)=5;
else
    temp=(dist_all-dist_red);
end
dist_A={temp};
end
function [percentage_time]=get_mean_pos_in_timebin(X,D,border,time_bin,mean_or_median)
time_tot=0;
point=1;
sum(D)/60000
x_temporary=[];
point_t=1;
for n=1:(length (X))
    if(time_tot<time_bin)
        time_tot=time_tot+D(n);
       if(n==length(X) && time_tot > 60000 )
            if(strcmp(mean_or_median,'mean'))
                x_temporary= border-X(point_t:n);
                percentage_time(point)= mean(x_temporary);
            else
                                x_temporary= border-X(point_t:n);
                percentage_time(point)= median(x_temporary);
            end
            x_temporary=[];
            point_t=n;
            time_tot=0;
            point=point+1;      
        end
    else
        if(strcmp(mean_or_median,'mean'))
                         x_temporary= border-X(point_t:n);
            percentage_time(point)= mean(x_temporary);
        else
                            x_temporary= border-X(point_t:n);
            percentage_time(point)= median(x_temporary);
        end
        x_temporary=[];
        time_tot=0;
        point=point+1;  
          point_t=n;
    end
end
end
function [Binned_pos]=binned_mean_position_flipping_color(X,Y,D,info_vector,color,time_bin,base_time,Flip_time,mean_or_median)
for i=1:size(X,2)
    max_size=max(info_vector{3,i},info_vector{2,i});
    for j=1:(size(X,1))
        x=X{j,i};
        y=Y{j,i};
        d=D{j,i};
%         [x_3,y_3]=remove_pos_spike(x, y);filter your data in case here
%         x=[];  y=[];
%         x=x_3;
%         y=y_3;
        if(size(x,2)==1)
            x=x';
        end
        if(size(y,2)==1)
            y=y';
        end
        if(size(d,2)==1)
            d=d';
        end
        if(size(info_vector,1)==4)
            if( strcmp(Flip_time,'NO'))
                thresh=info_vector(3,:);
                border=thresh{i};
            else
                if(j<Flip_time  )
                    thresh=info_vector(3,:);
                    border=thresh{i};
                else
                    if(isempty(info_vector{4,1}))
                        thresh=info_vector(3,:);
                        border=thresh{i};
                    else
                        thresh=info_vector(4,:);
                        border=thresh{i};
                    end
                end
            end
        else
            thresh=info_vector(3,:);
            border=thresh{i};
        end
        if(j==1)
            time=0;
            for k=max(size(x)):-1:1
                if(time<base_time && k>1)
                    time=time+d(k);
                else
                    break
                end
            end
            Binned_pos(j,i)={get_mean_pos_in_timebin(x(k:end),d(k:end),border,time_bin,mean_or_median)};  
            clearvars k
        else
            if (strcmp(Flip_time,'NO'))
                Binned_pos(j,i)={get_mean_pos_in_timebin(x,d,border,time_bin,mean_or_median)};
            else
                if(j<Flip_time )      
                    Binned_pos(j,i)={get_mean_pos_in_timebin(x,d,border,time_bin,mean_or_median)};
                else
                    Binned_pos(j,i)={get_mean_pos_in_timebin(x,d,border,time_bin,mean_or_median)};
                end
            end
            clearvars x time  d y
        end
    end
    clearvars time timetot point percentage_time
end
end
function [Mean_std]=plot_binned_vel(Bin,color)
Bin
X_tot=NaN(size(Bin,2),29);
for j=1:size(Bin,2)
    X=[];
    for i=1:size(Bin,1)
        x=Bin{i,j};
     X=[X,x];      
    end
    plot(X)
    hold on 
    size(X)
    X_tot(j,1:size(X,2))=X;
end
errorbar(1:size(X_tot,2),nanmean(X_tot,1),nanstd(X_tot,1)./sqrt(size(X_tot,1)),color ,'linewidth', 2 ) 
hold on 
Mean_std(1,1:size(X_tot,2))=nanmean(X_tot,1);
Mean_std(2,1:size(X_tot,2))=nanstd(X_tot,1)./sqrt(size(X_tot,1));
end
function []=scattercloud(X,Y,filled,color,size)% 1 if you want the scatterplot filled 
data=-0.2:0.02:+0.2;
for i=1:length(X)
    randomnum= datasample(data,1); % returns k observations sampled uniformly at random, with replacement, from the data in data.
   X(i)=X(i)+randomnum; 
end
if(strcmp(filled,'filled'))
scatter(X,Y,size,filled,color);
else
  scatter(X,Y,size,color);
end
end
function [ANGLE_ALL,V,V_raw,Peaks,MeanP,BASE_V,BASE_V_raw,PRE_Peak_response,Peak_response]=plot_response_to_only_shock_normalized(V_input,X,Y,info,window,switch_rule,internal_window,variability_threshold,division_in_time,time_block)
% %%the varianle angles has the following information :
% angles(line,1) angle of aproaching the midline
% angles(line,2) x position of encounter of the midline
% angles(line,3) x position "time_window" frames before meeting the midline
% angles(line,4) y position of encounter of the midline
% angles(line,5) y position "time_window" frames before meeting the midline
% angles(line,6) mean velocity of appraching the midline

% angles(line,7) angle of leaving the midline
% angles(line,8) x position of encounter of the midline
% angles(line,9) x position "time_window" frames after meeting the midline
% angles(line,10) y position of encounter of the midline
% angles(line,11) y position "time_window" frames after meeting the midline
% angles(line,12) mean velocity of leaving the midline

for i=1:size(X,2)
    for j=1:size(X,1)
        line=1;
        x=[];y=[];v=[];
        x=X{j,i};
       if (j==1 && size(x,2)>28000)%% this is a control for the time division of the session which in our case were 30 min long
            x=x(27000:end);
            y=y(27000:end);
            v=v(27000:end);
       end
       x=x(((end/division_in_time)*(time_block-1))+1:((end/division_in_time)*(time_block)));
        v=V_input{j,i};
        y=Y{j,i};
        if( 1~=(strcmp(switch_rule,'NO')))
            if(j>=switch_rule && size(info,1)==4 )
                if(isempty(info{4,i}))
                    tresh=info{3,i}-variability_threshold;
                else
                    tresh=info{4,i}-variability_threshold;
                end
            else
                tresh=info{3,i}+variability_threshold;
                
            end
        else
            tresh=info{3,i}+variability_threshold;
        end
        
        count=window;
        A=[];
        B=[];
        angles=[];
        for t=window:(max(size(x))-(window+1))
            if(t>count && t<=(max(size(v))-(window+1)))
                if( 1~=(strcmp(switch_rule,'NO')))
                    if(j<switch_rule )
                        if (x(t)<tresh && x(t-internal_window)>tresh  && x(t+internal_window)<tresh)
                            count=t+window-1;
                            A(line,1:window)=v(t:t+window-1);
                            if(t>window)
                                HYP=[];CAT1=[];CAT2=[];
                                CAT1=x(t+window)-x(t);
                                CAT2=y(t+window)-y(t);
                                HYP=sqrt(CAT1.^2+CAT2.^2);
                                angles(line,8)=x(t);
                                angles(line,9)=x(t+window);
                                angles(line,10)=y(t);
                                angles(line,11)=y(t+window);
                                angles(line,12)=mean(mean(v(t:(t+window))));
                                if(HYP==0)
                                    angles(line,7)=361;
                                    angles(line,12)=-0.1;
                                else
                                    if(rad2deg(asin(CAT2/HYP))>=0)
                                        angles(line,7)= ( rad2deg(acos(CAT1/HYP)));%post crossing thresh angle
                                    else
                                        angles(line,7)= 360+(-( rad2deg(acos(CAT1/HYP))));%post crossing thresh angle
                                    end
                                end
                                HYP=[];CAT1=[];CAT2=[];
                                CAT1=x(t)-x(t-window);
                                CAT2=y(t)-y(t-window);
                                HYP=sqrt(CAT1.^2+CAT2.^2);
                                angles(line,2)=x(t);
                                angles(line,3)=x(t-window);
                                angles(line,4)=y(t);
                                angles(line,5)=y(t-window);
                                angles(line,6)=mean(mean(v((t-window):t)));
                                if(HYP==0)
                                    angles(line,1)=361;
                                else
                                    if(rad2deg(asin(CAT2/HYP))>=0)
                                        angles(line,1)= ( rad2deg(acos(CAT1/HYP)));%pre crossing thresh angle
                                    else
                                        angles(line,1)= 360+(-( rad2deg(acos(CAT1/HYP))));%pre crossing thresh angle
                                    end
                                end
                                B(line,1:window)=v(t-(window-1):t);
                            end
                            line=line+1;
                        end
                    else
                        if (x(t)>tresh && x(t-internal_window)<tresh )
                            count=t+window-1;
                            A(line,1:window)=v(t:t+window-1);
                            if(t>window)
                                HYP=[];CAT1=[];CAT2=[];
                                CAT1=x(t+window)-x(t);
                                CAT2=y(t+window)-y(t);
                                HYP=sqrt(CAT1.^2+CAT2.^2);
                                angles(line,8)=x(t);
                                angles(line,9)=x(t+window);
                                angles(line,10)=y(t);
                                angles(line,11)=y(t+window);
                                angles(line,12)=mean(mean(v(t:(t+window))));
                                if(HYP==0)
                                    angles(line,7)=361;
                                    angles(line,12)=-0.1;
                                else
                                    if(rad2deg(asin(CAT2/HYP))>=0)
                                        angles(line,7)= ( rad2deg(acos(CAT1/HYP)));%post crossing thresh angle
                                    else
                                        angles(line,7)= 360+(-( rad2deg(acos(CAT1/HYP))));%post crossing thresh angle
                                    end
                                end
                                HYP=[];CAT1=[];CAT2=[];
                                CAT1=x(t)-x(t-window);
                                CAT2=y(t)-y(t-window);
                                HYP=sqrt(CAT1.^2+CAT2.^2);
                                angles(line,2)=x(t);
                                angles(line,3)=x(t-window);
                                angles(line,4)=y(t);
                                angles(line,5)=y(t-window);
                                angles(line,6)=mean(mean(v((t-window):t)));
                                if(HYP==0)
                                    angles(line,1)=361;
                                    angles(line,6)=-0.1;
                                else
                                    if(rad2deg(asin(CAT2/HYP))>=0)
                                        angles(line,1)= ( rad2deg(acos(CAT1/HYP)));%pre crossing thresh angle
                                    else
                                        angles(line,1)= 360+(-( rad2deg(acos(CAT1/HYP))));%pre crossing thresh angle
                                    end
                                end
                                B(line,1:window)=v(t-(window-1):t);
                            end
                            line=line+1;
                        end
                    end
                else
                    if (x(t)<tresh && x(t-internal_window)>tresh && x(t+internal_window)<tresh)
                        count=t+window-1;
                        A(line,1:window)=v(t:t+window-1);
                        if(t>window )
                            HYP=[];CAT1=[];CAT2=[];
                            CAT1=x(t+window)-x(t);
                            CAT2=y(t+window)-y(t);
                            HYP=sqrt(CAT1.^2+CAT2.^2);
                            angles(line,8)=x(t);
                            angles(line,9)=x(t+window);
                            angles(line,10)=y(t);
                            angles(line,11)=y(t+window);
                            angles(line,12)=mean(mean(v(t:(t+window))));
                            t+window
                            if(HYP==0)
                                angles(line,7)=361;
                                angles(line,12)=-0.1;
                            else
                                if(rad2deg(asin(CAT2/HYP))>=0)
                                    angles(line,7)= ( rad2deg(acos(CAT1/HYP)));%post crossing thresh angle
                                else
                                    angles(line,7)= 360+(-( rad2deg(acos(CAT1/HYP))));%post crossing thresh angle
                                end
                            end
                            HYP=[];CAT1=[];CAT2=[];
                            CAT1=x(t)-x(t-window);
                            CAT2=y(t)-y(t-window);
                            HYP=sqrt(CAT1.^2+CAT2.^2);
                            angles(line,2)=x(t);
                            angles(line,3)=x(t-window);
                            angles(line,4)=y(t);
                            angles(line,5)=y(t-window);
                            angles(line,6)=mean(mean(v((t-window):t)));
                            if(HYP==0)
                                angles(line,1)=361;
                                angles(line,6)=-0.1;
                            else
                                if(rad2deg(asin(CAT2/HYP))>=0)
                                    angles(line,1)= ( rad2deg(acos(CAT1/HYP)));%pre crossing thresh angle
                                else
                                    angles(line,1)= 360+(-( rad2deg(acos(CAT1/HYP))));%pre crossing thresh angle
                                end
                            end
                            B(line,1:window)=v(t-(window-1):t);
                        end
                        line=line+1;
                    end
                end
            end
        end
        ANGLE_ALL(j,i)={angles};
        if(size(A,1)>1)
            V(j,i)={mean(A,1)};
            BASE_V(j,i)={mean(B,1)};
            V_raw(j,i)={A};
            BASE_V_raw(j,i)={B};
            A=[];
            B=[];
            pks = mean( V{j,i});
            Peaks(j,i)={pks};
            pks_pre = mean( BASE_V{j,i});
            Peaks_pre(j,i)={pks_pre};
            Peak_response(j,i)=mean(pks);
            PRE_Peak_response(j,i)=mean(pks_pre);
            MeanP(j,i)=(mean(pks)./mean(pks_pre));
            pks=[];locs=[];w=[];p=[];T=[];
            pks_pre=[];locs_pre=[];w_pre=[];p_pre=[];T_pre=[];
            hold on
        else
            if(size(A,1)==1)
                V(j,i)={A};
                BASE_V(j,i)={B};
                V_raw(j,i)={A};
                BASE_V_raw(j,i)={B};
                A=[];
                B=[];
                pks = ( V{j,i});
                Peaks(j,i)={pks};
                pks_pre = ( BASE_V{j,i});
                Peaks_pre(j,i)={pks_pre};
                Peak_response(j,i)=mean(pks);
                PRE_Peak_response(j,i)=mean(pks_pre);
                MeanP(j,i)=(mean(pks)./mean(pks_pre));
                pks=[];locs=[];w=[];p=[];T=[];
                pks_pre=[];locs_pre=[];w_pre=[];p_pre=[];T_pre=[];
                
                hold on
            else
                V(j,i)={-0.1};
                BASE_V(j,i)={-0.1};
                V_raw(j,i)={A};
                BASE_V_raw(j,i)={B};
                A=[];
                B=[];
                pks = ( V{j,i});
                Peaks(j,i)={pks};
                pks_pre = ( BASE_V{j,i});
                Peaks_pre(j,i)={pks_pre};
                Peak_response(j,i)=(pks);
                PRE_Peak_response(j,i)=(pks_pre);
                MeanP(j,i)=((pks)./(pks_pre));
                pks=[];locs=[];w=[];p=[];T=[];
                pks_pre=[];locs_pre=[];w_pre=[];p_pre=[];T_pre=[];
                hold on
            end
        end
    end
end
end
function [ANGLE_ALL,V,V_raw,Peaks,MeanP,BASE_V,BASE_V_raw,PRE_Peak_response,Peak_response]=plot_avoidance_of_shock(V_input,X,Y,info,window,switch_rule,internal_window,variability_threshold,division_in_time,time_block)
% %%the varianle angles has the following information :
% angles(line,1) angle of aproaching the midline
% angles(line,2) x position of encounter of the midline
% angles(line,3) x position "time_window" frames before meeting the midline
% angles(line,4) y position of encounter of the midline
% angles(line,5) y position "time_window" frames before meeting the midline
% angles(line,6) mean velocity of appraching the midline

% angles(line,7) angle of leaving the midline
% angles(line,8) x position of encounter of the midline
% angles(line,9) x position "time_window" frames after meeting the midline
% angles(line,10) y position of encounter of the midline
% angles(line,11) y position "time_window" frames after meeting the midline
% angles(line,12) mean velocity of leaving the midline

for i=1:size(X,2)
    for j=1:size(X,1)
        line=1;
        x=[];y=[];v=[];
        x=X{j,i};
       if (j==1 && size(x,2)>28000)%% this is a control for the time division of the session which in our case were 30 min long
            x=x(27000:end);
            y=y(27000:end);
            v=v(27000:end);
       end
       x=x(((end/division_in_time)*(time_block-1))+1:((end/division_in_time)*(time_block)));
        v=V_input{j,i};
        y=Y{j,i};
        if( 1~=(strcmp(switch_rule,'NO')))
            if(j>=switch_rule && size(info,1)==4 )
                if(isempty(info{4,i}))
                    tresh=info{3,i}-variability_threshold;
                else
                    tresh=info{4,i}-variability_threshold;
                end
            else
                tresh=info{3,i}+variability_threshold;
                
            end
        else
            tresh=info{3,i}+variability_threshold;
        end
        
        count=window;
        A=[];
        B=[];
        angles=[];
        for t=window:(max(size(x))-(window+1))
            if(t>count && t<=(max(size(v))-(window+1)))
                if( 1~=(strcmp(switch_rule,'NO')))
                    if(j<switch_rule )
                        if (x(t)<tresh && x(t-internal_window)>tresh )
                            count=t+window-1;
                            A(line,1:window)=v(t:t+window-1);
                            if(t>window)
                                HYP=[];CAT1=[];CAT2=[];
                                CAT1=x(t+window)-x(t);
                                CAT2=y(t+window)-y(t);
                                HYP=sqrt(CAT1.^2+CAT2.^2);
                                angles(line,8)=x(t);
                                angles(line,9)=x(t+window);
                                angles(line,10)=y(t);
                                angles(line,11)=y(t+window);
                                angles(line,12)=mean(mean(v(t:(t+window))));
                                if(HYP==0)
                                    angles(line,7)=361;
                                    angles(line,12)=-0.1;
                                else
                                    if(rad2deg(asin(CAT2/HYP))>=0)
                                        angles(line,7)= ( rad2deg(acos(CAT1/HYP)));%post crossing thresh angle
                                    else
                                        angles(line,7)= 360+(-( rad2deg(acos(CAT1/HYP))));%post crossing thresh angle
                                    end
                                end
                                HYP=[];CAT1=[];CAT2=[];
                                CAT1=x(t)-x(t-window);
                                CAT2=y(t)-y(t-window);
                                HYP=sqrt(CAT1.^2+CAT2.^2);
                                angles(line,2)=x(t);
                                angles(line,3)=x(t-window);
                                angles(line,4)=y(t);
                                angles(line,5)=y(t-window);
                                angles(line,6)=mean(mean(v((t-window):t)));
                                if(HYP==0)
                                    angles(line,1)=361;
                                else
                                    if(rad2deg(asin(CAT2/HYP))>=0)
                                        angles(line,1)= ( rad2deg(acos(CAT1/HYP)));%pre crossing thresh angle
                                    else
                                        angles(line,1)= 360+(-( rad2deg(acos(CAT1/HYP))));%pre crossing thresh angle
                                    end
                                end
                                B(line,1:window)=v(t-(window-1):t);
                            end
                            line=line+1;
                        end
                    else
                        if (x(t)>tresh && x(t-internal_window)<tresh )
                            count=t+window-1;
                            A(line,1:window)=v(t:t+window-1);
                            if(t>window)
                                HYP=[];CAT1=[];CAT2=[];
                                CAT1=x(t+window)-x(t);
                                CAT2=y(t+window)-y(t);
                                HYP=sqrt(CAT1.^2+CAT2.^2);
                                angles(line,8)=x(t);
                                angles(line,9)=x(t+window);
                                angles(line,10)=y(t);
                                angles(line,11)=y(t+window);
                                angles(line,12)=mean(mean(v(t:(t+window))));
                                if(HYP==0)
                                    angles(line,7)=361;
                                    angles(line,12)=-0.1;
                                else
                                    if(rad2deg(asin(CAT2/HYP))>=0)
                                        angles(line,7)= ( rad2deg(acos(CAT1/HYP)));%post crossing thresh angle
                                    else
                                        angles(line,7)= 360+(-( rad2deg(acos(CAT1/HYP))));%post crossing thresh angle
                                    end
                                end
                                HYP=[];CAT1=[];CAT2=[];
                                CAT1=x(t)-x(t-window);
                                CAT2=y(t)-y(t-window);
                                HYP=sqrt(CAT1.^2+CAT2.^2);
                                angles(line,2)=x(t);
                                angles(line,3)=x(t-window);
                                angles(line,4)=y(t);
                                angles(line,5)=y(t-window);
                                angles(line,6)=mean(mean(v((t-window):t)));
                                if(HYP==0)
                                    angles(line,1)=361;
                                    angles(line,6)=-0.1;
                                else
                                    if(rad2deg(asin(CAT2/HYP))>=0)
                                        angles(line,1)= ( rad2deg(acos(CAT1/HYP)));%pre crossing thresh angle
                                    else
                                        angles(line,1)= 360+(-( rad2deg(acos(CAT1/HYP))));%pre crossing thresh angle
                                    end
                                end
                                B(line,1:window)=v(t-(window-1):t);
                            end
                            line=line+1;
                        end
                    end
                else
                    if (x(t)<tresh && x(t-internal_window)>tresh)
                        count=t+window-1;
                        A(line,1:window)=v(t:t+window-1);
                        if(t>window )
                            HYP=[];CAT1=[];CAT2=[];
                            CAT1=x(t+window)-x(t);
                            CAT2=y(t+window)-y(t);
                            HYP=sqrt(CAT1.^2+CAT2.^2);
                            angles(line,8)=x(t);
                            angles(line,9)=x(t+window);
                            angles(line,10)=y(t);
                            angles(line,11)=y(t+window);
                            angles(line,12)=mean(mean(v(t:(t+window))));
                            t+window
                            if(HYP==0)
                                angles(line,7)=361;
                                angles(line,12)=-0.1;
                            else
                                if(rad2deg(asin(CAT2/HYP))>=0)
                                    angles(line,7)= ( rad2deg(acos(CAT1/HYP)));%post crossing thresh angle
                                else
                                    angles(line,7)= 360+(-( rad2deg(acos(CAT1/HYP))));%post crossing thresh angle
                                end
                            end
                            HYP=[];CAT1=[];CAT2=[];
                            CAT1=x(t)-x(t-window);
                            CAT2=y(t)-y(t-window);
                            HYP=sqrt(CAT1.^2+CAT2.^2);
                            angles(line,2)=x(t);
                            angles(line,3)=x(t-window);
                            angles(line,4)=y(t);
                            angles(line,5)=y(t-window);
                            angles(line,6)=mean(mean(v((t-window):t)));
                            if(HYP==0)
                                angles(line,1)=361;
                                angles(line,6)=-0.1;
                            else
                                if(rad2deg(asin(CAT2/HYP))>=0)
                                    angles(line,1)= ( rad2deg(acos(CAT1/HYP)));%pre crossing thresh angle
                                else
                                    angles(line,1)= 360+(-( rad2deg(acos(CAT1/HYP))));%pre crossing thresh angle
                                end
                            end
                            B(line,1:window)=v(t-(window-1):t);
                        end
                        line=line+1;
                    end
                end
            end
        end
        ANGLE_ALL(j,i)={angles};
        if(size(A,1)>1)
            V(j,i)={mean(A,1)};
            BASE_V(j,i)={mean(B,1)};
            V_raw(j,i)={A};
            BASE_V_raw(j,i)={B};
            A=[];
            B=[];
            pks = mean( V{j,i});
            Peaks(j,i)={pks};
            pks_pre = mean( BASE_V{j,i});
            Peaks_pre(j,i)={pks_pre};
            Peak_response(j,i)=mean(pks);
            PRE_Peak_response(j,i)=mean(pks_pre);
            MeanP(j,i)=(mean(pks)./mean(pks_pre));
            pks=[];locs=[];w=[];p=[];T=[];
            pks_pre=[];locs_pre=[];w_pre=[];p_pre=[];T_pre=[];
            hold on
        else
            if(size(A,1)==1)
                V(j,i)={A};
                BASE_V(j,i)={B};
                V_raw(j,i)={A};
                BASE_V_raw(j,i)={B};
                A=[];
                B=[];
                pks = ( V{j,i});
                Peaks(j,i)={pks};
                pks_pre = ( BASE_V{j,i});
                Peaks_pre(j,i)={pks_pre};
                Peak_response(j,i)=mean(pks);
                PRE_Peak_response(j,i)=mean(pks_pre);
                MeanP(j,i)=(mean(pks)./mean(pks_pre));
                pks=[];locs=[];w=[];p=[];T=[];
                pks_pre=[];locs_pre=[];w_pre=[];p_pre=[];T_pre=[];
                
                hold on
            else
                V(j,i)={-0.1};
                BASE_V(j,i)={-0.1};
                V_raw(j,i)={A};
                BASE_V_raw(j,i)={B};
                A=[];
                B=[];
                pks = ( V{j,i});
                Peaks(j,i)={pks};
                pks_pre = ( BASE_V{j,i});
                Peaks_pre(j,i)={pks_pre};
                Peak_response(j,i)=(pks);
                PRE_Peak_response(j,i)=(pks_pre);
                MeanP(j,i)=((pks)./(pks_pre));
                pks=[];locs=[];w=[];p=[];T=[];
                pks_pre=[];locs_pre=[];w_pre=[];p_pre=[];T_pre=[];
                hold on
            end
        end
    end
end
end
function [mu ul ll] = circ_mean(alpha, w)
%
% mu = circ_mean(alpha, w)
%   Computes the mean direction for circular data.
%
%   Input:
%     alpha	sample of angles in radians
%     [w		weightings in case of binned angle data]
%
%   Output:
%     mu		mean direction
%     ul    upper 95% confidence limit
%     ll    lower 95% confidence limit 
%
% PHB 7/6/2008
%
% References:
%   Statistical analysis of circular data, N. I. Fisher
%   Topics in circular statistics, S. R. Jammalamadaka et al. 
%   Biostatistical Analysis, J. H. Zar
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html

% check vector size
if size(alpha,2) > size(alpha,1)
	alpha = alpha';
end

if nargin<2
  % if no specific weighting has been specified
  % assume no binning has taken place
	w = ones(size(alpha));
else
  if size(w,2) > size(w,1)
    w = w';
  end 
end

% compute weighted sum of cos and sin of angles
r = w'*exp(1i*alpha);

% obtain mean by
mu = angle(r);

% confidence limits if desired
if nargout > 1
  t = circ_confmean(alpha,0.05,w);
  ul = mu + t;
  ll = mu - t;
end
end
function [pdf_vector, pdf_mean, X_vector,Y_vector, D_vector, info_vector, Size_b ,V_vector] = import_data_general(days, month_array, YEAR,user,exp_number, filename, PROTOCOL, blocknum,num_fish_total,roitoloop)

pdf_vector=cell(blocknum,num_fish_total);
pdf_mean=cell(blocknum,1);
Y_vector=cell(blocknum,num_fish_total);
D_vector=cell(blocknum,num_fish_total);
Size_b=zeros(1,num_fish_total);
num_of_days=1;
for protocol_number=1:size(PROTOCOL,2)
if(strcmp(PROTOCOL{protocol_number},'add_day'))
    num_of_days=num_of_days+1;
end
end
info_vector=cell(4,num_fish_total+1,num_of_days);

info_vector(3,end)={'threshold_first_cond'};
info_vector(4,end)={'threshold_secodn_cond'};



line=-1;

for dayloop=1: size(days,2)
    year=YEAR{dayloop};
    exp_num=exp_number{dayloop};
    line=line+1;
    c=days{dayloop}
    c_first_day=c;
    month=month_array{dayloop}
    User=user{dayloop};
    i=0;
    

    for j=roitoloop
         
       clearvars xmaxdim;
        year=YEAR{dayloop};
    exp_num=exp_number{dayloop};
    c=days{dayloop}
    c_first_day=c;
    month=month_array{dayloop}
    User=user{dayloop};

        i=i+1;
        if (strcmp(exp_num,'none'))
            filename_orig=strcat(filename,year,'-',month,'-',c_first_day,'\',year,'-',month,'-',c);
        else
            filename_orig=strcat(filename,year,'-',month,'-',c_first_day,'exp',exp_num,'\',year,'-',month,'-',c)
        end
        ROI=num2str(j)
        conditioning_done=0;
        baseline_done=0;
        after_trial_done=0;
        counter=0;
        count_days=1;
        cnew=[];
        for session=1:size(PROTOCOL,2)
            session_performed=PROTOCOL{session};
            if ((strcmp(session_performed,'add_day')))
                month_new=[];
                cnew=[];
                baseline_done=0;conditioning_done=0;after_trial_done=0;
                count_days=count_days+1
                
                
                if(str2num(c)<9)
                    if (strcmp(exp_num,'none'))
                        filename_orig=strcat(filename,year,'-',month,'-',c_first_day,'\',year,'-',month,'-0',num2str(str2num(c)+1));
                    else
                        filename_orig=strcat(filename,year,'-',month,'-',c_first_day,'exp',exp_num,'\',year,'-',month,'-0',num2str(str2num(c)+1));
                    end
                    if(isempty(cnew))
                        cnew=num2str(str2num(c)+1);
                        c=cnew;
                    end
                    c=(cnew);
                else
                    if((str2num(month)==1 || str2num(month)==3 || str2num(month)==5 ||str2num(month)==7 || str2num(month)==8 || str2num(month)==10) && (str2num(c)==31))
                        if(str2num(month)<9)
                            month_new=['0',num2str(str2num(month)+1)];
                        else
                            month_new=num2str(str2num(month)+1);
                        end
                        cnew='01'; c=cnew;
                    else
                        if(str2num(c)==28 &&  str2num(month)==2)
                            if(str2num(month)<9)
                                month_new=['0',num2str(str2num(month)+1)];
                            else
                                month_new=num2str(str2num(month)+1);
                            end
                            cnew='01'; c=cnew;
                        else
                            if( str2num(c)==30)
                                if(str2num(month)<9)
                                    month_new=['0',num2str(str2num(month)+1)];
                                else
                                    month_new=num2str(str2num(month)+1);
                                end
                                cnew='01'; c=cnew;
                            else
                                if(isempty(month_new))
                                    month_new=month;
                                end
                                if(isempty(cnew))
                                    cnew=num2str(str2num(c)+1); c=cnew;
                                end
                                
                            end
                        end
                    end
                    if (strcmp(exp_num,'none'))
                        filename_orig=strcat(filename,year,'-',month,'-',c_first_day,'\',year,'-',month_new,'-',cnew);
                    else
                        filename_orig=strcat(filename,year,'-',month,'-',c_first_day,'exp',exp_num,'\',year,'-',month_new,'-',cnew);
                    end
                end          
                baseline_done=0;conditioning_done=0;after_trial_done=0;
            end
            if (strcmp(session_performed,'baseline'))
                baseline_done=baseline_done+1;
                counter=counter+1;
                filename1 = strcat(filename_orig,'_',User,'_baseline_exp_',num2str(baseline_done),'_ROI',ROI,'.txt');
                
                [Db,Xb,Yb,size_b] = import_baseline_with_size(filename1);
                if( baseline_done==1 && count_days==1)
                    info_vector(1,line*6+i,count_days)={importheigth_roi_function(filename1,10,10)};
                    info_vector(2,line*6+i,count_days)={importwidth_roi_function(filename1,9,9)};
                    info_vector(1,end,count_days)={'height baseline'};
                    info_vector(2,end,count_days)={'width baseline'};
                    X_vector(1,line*6+i)={[]};
                end
                 info_vector(1,line*6+i,count_days)={importheigth_roi_function(filename1,10,10)};
                    info_vector(2,line*6+i,count_days)={importwidth_roi_function(filename1,9,9)};
                   
                xmaxdim=max(importwidth_roi_function(filename1,9,9),importheigth_roi_function(filename1,10,10));
                [xb,yb,dtempb]=equally_spaced_vector(Xb,Yb,Db);
                sum(dtempb)./60000
                D_vector(counter,line*6+i)={dtempb};
                X_vector(counter,line*6+i)={xb};
                Y_vector(counter,line*6+i)={yb};
                Size_b(1,line*6+i)=mean(size_b(size_b>0));
                temp_base=get_PDF(xb',yb',xmaxdim);
                pdf_vector(counter,line*6+i)={temp_base};
                xb=[];yb=[];dtempb=[];Db=[];Xb=[];Yb=[];
                else
                if(strcmp(session_performed,'conditioning'))
                    conditioning_done=conditioning_done+1;%exp_num%
                    filename2 = strcat(filename_orig,'_',User,'_trial_exp_',num2str(conditioning_done),'_ROI',ROI,'.txt');
                    info_vector(2+conditioning_done,line*6+i,count_days)={importthreshold_function(filename2)};               
%                     info_vector(3,line*6+i,count_days)={importthreshold_function(filename2)};               
                    
[Dt,Xt,Yt,S,size_b] = importtrial_function(filename2);
       S(S~=2)=1;  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
if((Size_b(1,line*6+i))==0)
                                        Size_b(1,line*6+i)=mean(size_b(size_b>0));
                    end
                    [DT,XT,YT,ST]=detection_multiple_period_conditioning(Dt,Xt,Yt,S);
                     check=exist('xmaxdim');
                    if (check)
                    else
                          
                    info_vector(1,line*6+i,count_days)={importheigth_roi_function(filename2,12,12)};
                    info_vector(2,line*6+i,count_days)={importwidth_roi_function(filename2,13,13)};
                    info_vector(1,end,count_days)={'height baseline'};
                    info_vector(2,end,count_days)={'width baseline'};
                                    xmaxdim=max(importwidth_roi_function(filename2,13,13),importheigth_roi_function(filename2,14,14));

                      end
                    for counter_block=1:size(DT,2)
                        
                        xt1=[];
                        dtempt1=[];
                        yt1=[];
                        temp_t1=[];
                        [xt1,yt1,dtempt1]=equally_spaced_vector(XT{counter_block},YT{counter_block},DT{counter_block});
                        if(sum(dtempt1) <= 600000)
                            DT(counter_block)=[];
                            XT(counter_block)=[];
                            YT(counter_block)=[];
                            ST(counter_block)=[];
                            
                            continue
                        end
                        
                        
                        temp_t1=get_PDF(xt1',yt1',xmaxdim);
                        pdf_vector(counter_block+counter,line*6+i)={temp_t1};
                        X_vector(counter_block+counter,line*6+i)={xt1};
                        Y_vector(counter_block+counter,line*6+i)={yt1};
                        D_vector(counter_block+counter,line*6+i)={dtempt1};
                        xt1=[];
                        dtempt1=[];
                        yt1=[];
                        temp_t1=[];
                        
                    end
                    counter=counter+size(DT,2);
                    DT=[];XT=[];YT=[];ST=[];Dt=[];Xt=[];Yt=[];S=[];
                    
                else
                    if(strcmp(session_performed,'after_trial'))
                        
                        after_trial_done=after_trial_done+1;
                        counter=counter+1;
                        
                        filename3 = strcat(filename_orig,'_',User,'_after_trial_exp_',num2str(after_trial_done),'_ROI',ROI,'.txt')
                        
                        [Db,Xb,Yb] = importafter_trial_function(filename3);
                        if(Db(1)<0)
                            Xb(1)=[];
                            Yb(1)=[];
                            Db(1)=[];
                        end
                        
                        if((sum(Db)/60000)<29.5)
                            Xb=[Xb;Xb(end)];
                            Yb=[Yb;Yb(end)];
                            Db=[Db;((29.5)*60000)-(sum(Db))];
                        end
                        
                        [xb,yb,dtempb]=equally_spaced_vector(Xb,Yb,Db);
                        D_vector(counter,line*6+i)={dtempb};
                        X_vector(counter,line*6+i)={xb};
                        Y_vector(counter,line*6+i)={yb};
                        
                        temp_base=get_PDF(xb',yb',xmaxdim);
                        pdf_vector(counter,line*6+i)={temp_base};
                        xb=[];yb=[];dtempb=[];Db=[];Xb=[];Yb=[];
                        
                        
                    end
                    
                end
                
            end
            
        end
        
    end
    
end
% V_vector=get_velocity_new(X_vector,Y_vector,D_vector,info_vector);

V_vector=get_velocity_accurate(X_vector,Y_vector,D_vector,info_vector);
end
function [V]=get_velocity_new(X,Y,D,info_vector)
V=cell (size (D));
for j=1:size(D,1)
    
    for i=1:size(D,2)  
      
        tempd= [];
        tempx= [];
        tempy= [];
        tempd= D{j,i};
        tempx= X{j,i};
        tempy= Y{j,i};    
        V(j)={((((tempx(2:end)-tempx(1:end-1)).^2)+((tempy(2:end)-tempy(1:end-1)).^2)).^(0.5))./tempd(1:end-1) };
        temp_velocity=(((diff(tempx).^2)+(diff(tempy).^2)).^(0.5))./tempd(1:end-1);
        for t=4:(length(temp_velocity)-3)
            if(temp_velocity(t)>=0.4)
                if(temp_velocity(t-3)~=0 | temp_velocity(t+3)~=0 | temp_velocity(t-1)~=0 ||  temp_velocity(t+1)~=0 | temp_velocity(t-2)~=0 |  temp_velocity(t+2)~=0)
                    temp_velocity(t) =mean(temp_velocity([t-3,t-2,t-1,t+1,t+2,t+3]));%% or NaN??????
                else
                    temp_velocity(t) =0;
                end                 
            end
        end
        V(j,i)={temp_velocity.*(12000/min(info_vector{1:2,i}))};
    end
end

end

function [x_test,y_test,d_test] = equally_spaced_vector(x,y,d)
x_temp=x;
y_temp=y;
d_temp=d;
num=1;

        repetition=1;
        if(d_temp(1)<70)
         x_test(num)=x_temp(1);
                y_test(num)=y_temp(1);
                d_test(num)=d_temp(1);
                num=num+1;
        else
            repetition=round(d_temp(1)/67);
                dist_x=(x_temp(2)-x_temp(1));
                dist_y=(y_temp(2)-y_temp(1));
                for rep=1:(repetition)
                    x_test(num)=x_temp(1)+(dist_x/repetition);
                    y_test(num)=y_temp(1)+(dist_y/repetition);
                    d_test(num)=(d_temp(1)./repetition);
                    
                    num=num+1;
                end
        end
            
            
        
        for k=2:max(size(x_temp))
           
            if(d_temp(k)<70)
                x_test(num)=x_temp(k);
                y_test(num)=y_temp(k);
                d_test(num)=d_temp(k);
                num=num+1;
            else
                repetition=round(d_temp(k)/67);
                dist_x=(x_temp(k)-x_temp(k-1));
                dist_y=(y_temp(k)-y_temp(k-1));
                for rep=1:(repetition)
%                     if(rep<=repetition) 
                    x_test(num)=x_temp(k-1)+(dist_x/repetition);
                    y_test(num)=y_temp(k-1)+(dist_y/repetition);
                    d_test(num)=(d_temp(k)./repetition);
                    
                    num=num+1;
% %                     else
%                                                 test_time=((d_temp(k)./66).*);
% 
%                     if(test_time<d_temp(k))       
% 
%                        d_test(num)=(d_temp(k)-test_time);
%                        num=num+1;
%                     end
                end
            end
        end
        [x_temp,y_temp] = remove_pos_spike(x_temp,y_temp);

if(num==1)
%     y_test=(y);
% d_test=(d);
%       x_test=(x);  
else
y_test=round(y_test);
d_test=round(d_test);
      x_test=round(x_test);  
end    
end

function [Da,Xa,Ya] = importafter_trial_function(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as column vectors.
%   [DA,XA,YA] = IMPORTFILE(FILENAME) Reads data from text file FILENAME
%   for the default selection.
%
%   [DA,XA,YA] = IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads data from
%   rows STARTROW through ENDROW of text file FILENAME.
%
% Example:
%   [Da,Xa,Ya] =
%   importfile('2016-03-11_Fabrizio_after_trial_exp_1_ROI6.txt',14, 38650);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2016/03/14 14:28:03

%% Initialize variables.
delimiter = {',',':'};
if nargin<=2
    startRow = 13;
    endRow = inf;
end

%% Format string for each line of text:
%   column2: double (%f)
%	column4: double (%f)
%   column6: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%*s%f%*s%f%*s%f%*s%*s%*s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
textscan(fileID, '%[^\n\r]', startRow(1)-1, 'ReturnOnError', false);
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    textscan(fileID, '%[^\n\r]', startRow(block)-1, 'ReturnOnError', false);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'ReturnOnError', false);
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Allocate imported array to column variable names
D= dataArray{:, 1};
X = dataArray{:, 2};
Y = dataArray{:, 3};
Min_size=min([length(D),length(X),length(Y)]);
Da = D(1:Min_size);
Xa = X(1:Min_size);
Ya = Y(1:Min_size);
end

function [Dt,Xt,Yt,ST]=detection_multiple_period_conditioning (D,X,Y,S)
t=1;
Dt={};Xt={};Yt={};ST={};
previous=1;
for j=1:length(X)-1

%     if(S(j)<2 && S(j+1)==2 || S(j)==2 && S(j+1)==3 )
% if(S(j)<2 && S(j+1)==2 )
 if(S(j) ~= S(j+1) )

    current=j;

if(sum(D(previous:current))/60000 <(29.9) && D(j+1) > ((29.9*60000)-sum(D(previous:current)) ))

    Dt(t)={[D(previous:current);((29.9*60000)-sum(D(previous:current)))]};
  Xt(t)={[X(previous:current);X(current)]};
  Yt(t)={[Y(previous:current);Y(current)]};
  ST(t)={[S(previous:current);S(current)]};
       
D(j+1)=D(j+1)-((29.9*60000)-sum(D(previous:current))); t=t+1;
previous=j+1;
else
  Dt(t)={D(previous:current)};
  Xt(t)={X(previous:current)};
  Yt(t)={Y(previous:current)};
  ST(t)={S(previous:current)};
        t=t+1;
previous=j+1;
end
    end
%    if(S(j)==2 && S(j+1)~=2  )
% current=j;
% if(sum(D(previous:current))/60000 <(29.9) && D(j+1) > ((29.9*60000)-sum(D(previous:current)) ))
% 
%     Dt(t)={[D(previous:current);((29.9*60000)-sum(D(previous:current)))]};
%   Xt(t)={[X(previous:current);X(current)]};
%   Yt(t)={[Y(previous:current);Y(current)]};
%   ST(t)={[S(previous:current);S(current)]};
%        
% D(j+1)=D(j+1)-((29.9*60000)-sum(D(previous:current))); t=t+1;
% previous=j+1;
% else
% 
%   Dt(t)={D(previous:current)};
%   Xt(t)={X(previous:current)};
%   Yt(t)={Y(previous:current)};
%   ST(t)={S(previous:current)};
%         t=t+1;
% previous=j+1;
% 
% end
%    end
if(j==(length(X)-1) && previous < (length(X)-1) )
current=j+1;

  if((current-previous )> 150)
 Dt(t)={D(previous:current)};
  Xt(t)={X(previous:current)};
  Yt(t)={Y(previous:current)};
  ST(t)={S(previous:current)};
  end
end

end

end

function [Dt, Xt, Yt, S, size_b] = importtrial_function(filename, dataLines)
%IMPORTFILE Import data from a text file
%  [DT, XT, YT, S, SIZE_B] = IMPORTFILE(FILENAME) reads data from text
%  file FILENAME for the default selection.  Returns the data as column
%  vectors.
%
%  [DT, XT, YT, S, SIZE_B] = IMPORTFILE(FILE, DATALINES) reads data for
%  the specified row interval(s) of text file FILENAME. Specify
%  DATALINES as a positive scalar integer or a N-by-2 array of positive
%  scalar integers for dis-contiguous row intervals.
%
%  Example:
%  [Dt, Xt, Yt, S, size_b] = importfile("\\forskning.it.ntnu.no\ntnu\mh-kin\yaksi2\bramse\Behavioural_experiments\New_extinction_paradigm\2020-03-06\2020-03-06_Bram_trial_exp_1_ROI5.txt", [18, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 13-Mar-2020 11:09:05

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [18, Inf];
end

%% Setup the Import Options
opts = delimitedTextImportOptions("NumVariables", 10);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = [",", ":"];

% Specify column names and types
opts.VariableNames = ["Var1", "Dt", "Var3", "Xt", "Var5", "Yt", "Var7", "S", "Var9", "size_b"];
opts.SelectedVariableNames = ["Dt", "Xt", "Yt", "S", "size_b"];
opts.VariableTypes = ["string", "double", "string", "double", "string", "double", "string", "double", "string", "double"];
opts = setvaropts(opts, [1, 3, 5, 7, 9], "WhitespaceRule", "preserve");
opts = setvaropts(opts, [1, 3, 5, 7, 9], "EmptyFieldRule", "auto");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
tbl = readtable(filename, opts);

%% Convert to output type
Dt = tbl.Dt;
Xt = tbl.Xt;
Yt = tbl.Yt;
S = tbl.S;
size_b = tbl.size_b;
end

function delta = importthreshold_function(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as column vectors.
%   DELTA = IMPORTFILE(FILENAME) Reads data from text file FILENAME for the
%   default selection.
%
%   DELTA = IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads data from rows
%   STARTROW through ENDROW of text file FILENAME.
%
% Example:
%   delta = importfile('2016-03-11_Fabrizio_trial_exp_1_ROI1.txt',15, 15);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2016/03/14 14:26:40

%% Initialize variables.
delimiter = {',',':'};
if nargin<=2
    startRow = 15;
    endRow = 15;
end

%% Format string for each line of text:
%   column1: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%f%*s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
textscan(fileID, '%[^\n\r]', startRow(1)-1, 'ReturnOnError', false);
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    textscan(fileID, '%[^\n\r]', startRow(block)-1, 'ReturnOnError', false);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'ReturnOnError', false);
    dataArray{1} = [dataArray{1};dataArrayBlock{1}];
end

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Allocate imported array to column variable names
delta = dataArray{:, 1};
end

function [Db,Xb,Yb,Size_b] = import_baseline_with_size(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as column vectors.
%   [DB,XB,YB,SIZE_B] = IMPORTFILE(FILENAME) Reads data from text file
%   FILENAME for the default selection.
%
%   [DB,XB,YB,SIZE_B] = IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads data
%   from rows STARTROW through ENDROW of text file FILENAME.
%
% Example:
%   [Db,Xb,Yb,Size_b] = importfile('2017-11-23_Fabrizio_baseline_exp_1_ROI2.txt',14, 23133);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2017/12/05 14:19:04

%% Initialize variables.
delimiter = {',',' ',':'};
if nargin<=2
    startRow = 14;
    endRow = inf;
end

%% Format string for each line of text:
%   column2: double (%f)
%	column4: double (%f)
%   column6: double (%f)
%	column8: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%*s%f%*s%f%*s%f%*s%f%*s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
textscan(fileID, '%[^\n\r]', startRow(1)-1, 'WhiteSpace', '', 'ReturnOnError', false);
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    textscan(fileID, '%[^\n\r]', startRow(block)-1, 'WhiteSpace', '', 'ReturnOnError', false);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'ReturnOnError', false);
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Allocate imported array to column variable names
Db = dataArray{:, 1};
Xb = dataArray{:, 2};
Yb = dataArray{:, 3};
Size_b = dataArray{:, 4};
end

function Dt = importheigth_roi_function(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as column vectors.
%   DT = IMPORTFILE(FILENAME) Reads data from text file FILENAME for the
%   default selection.
%
%   DT = IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads data from rows
%   STARTROW through ENDROW of text file FILENAME.
%
% Example:
%   Dt = importfile('2016-03-11_Fabrizio_trial_exp_1_ROI1.txt',13, 13);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2016/03/14 14:26:10

%% Initialize variables.
delimiter = {',',':'};
if nargin<=2
    startRow = 13;
    endRow = 13;
end

%% Format string for each line of text:
%   column2: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%*s%f%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
textscan(fileID, '%[^\n\r]', startRow(1)-1, 'ReturnOnError', false);
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    textscan(fileID, '%[^\n\r]', startRow(block)-1, 'ReturnOnError', false);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'ReturnOnError', false);
    dataArray{1} = [dataArray{1};dataArrayBlock{1}];
end

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Allocate imported array to column variable names
Dt = dataArray{:, 1};

end

function Dt = importwidth_roi_function(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as column vectors.
%   DT = IMPORTFILE(FILENAME) Reads data from text file FILENAME for the
%   default selection.
%
%   DT = IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads data from rows
%   STARTROW through ENDROW of text file FILENAME.
%
% Example:
%   Dt = importfile('2016-03-11_Fabrizio_trial_exp_1_ROI1.txt',12, 12);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2016/03/14 14:25:35

%% Initialize variables.
delimiter = {',',':'};
if nargin<=2
    startRow = 12;
    endRow = 12;
end

%% Format string for each line of text:
%   column2: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%*s%f%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
textscan(fileID, '%[^\n\r]', startRow(1)-1, 'ReturnOnError', false);
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    textscan(fileID, '%[^\n\r]', startRow(block)-1, 'ReturnOnError', false);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'ReturnOnError', false);
    dataArray{1} = [dataArray{1};dataArrayBlock{1}];
end

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Allocate imported array to column variable names
Dt = dataArray{:, 1};

end

function [x_test,y_test] = remove_pos_spike(X,Y)
pix_lim=100;
        repetition=1;
        x_test=X;%(1);
        y_test=Y;%(1);
        
        for j=2:(max(size(X))-4)
        
% if(j==1079)
% ra=2
% end
           d_uno=((((X(j)-X(j-1)).^2)+((Y(j)-Y(j-1)).^2)).^0.5); 
           d_due=((((X(j+1)-X(j-1)).^2)+((Y(j+1)-Y(j-1)).^2)).^0.5);
                      d_tre=((((X(j+2)-X(j-1)).^2)+((Y(j+2)-Y(j-1)).^2)).^0.5);
           d_quattro=((((X(j+3)-X(j-1)).^2)+((Y(j+3)-Y(j-1)).^2)).^0.5);
           d_cinque=((((X(j+4)-X(j-1)).^2)+((Y(j+4)-Y(j-1)).^2)).^0.5);

           if( d_uno > pix_lim && d_due < pix_lim  ) %
               x_test(j)= X (j-1)+ (X(j+1)-X(j-1))/2  ;
               y_test(j)= Y (j-1)+ (Y(j+1)-Y(j-1))/2  ;
           else
               if( d_uno > pix_lim && d_due > pix_lim && d_tre < pix_lim ) %
                   x_test(j)= X (j-1)+ (X(j+2)-X(j-1))/3  ;
                   y_test(j)= Y (j-1)+ (Y(j+2)-Y(j-1))/3  ;
                   x_test(j+1)= X (j-1)+ 2*((X(j+2)-X(j-1))/3)  ;
                   y_test(j+1)= Y (j-1)+ 2*((Y(j+2)-Y(j-1))/3 );
               else
                   if( d_uno > pix_lim && d_due > pix_lim && d_tre >pix_lim && d_quattro < pix_lim) %
                       x_test(j)= X (j-1)+ (X(j+3)-X(j-1))/4  ;
                       y_test(j)= Y (j-1)+ (Y(j+3)-Y(j-1))/4  ;
                       x_test(j+1)= X (j-1)+ 2*((X(j+3)-X(j-1))/4)  ;
                       y_test(j+1)= Y (j-1)+ 2*((Y(j+3)-Y(j-1))/4 );
                       x_test(j+2)= X (j-1)+ 3*((X(j+3)-X(j-1))/4)  ;
                       y_test(j+2)= Y (j-1)+ 3*((Y(j+3)-Y(j-1))/4 );
                   else
                       if( d_uno > pix_lim && d_due > pix_lim && d_tre >pix_lim && d_quattro > pix_lim && d_cinque < pix_lim) %
                           x_test(j)= X (j-1)+ (X(j+4)-X(j-1))/5  ;
                           y_test(j)= Y (j-1)+ (Y(j+4)-Y(j-1))/5  ;
                           x_test(j+1)= X (j-1)+ 2*((X(j+4)-X(j-1))/5)  ;
                           y_test(j+1)= Y (j-1)+ 2*((Y(j+4)-Y(j-1))/5 );
                           x_test(j+2)= X (j-1)+ 3*((X(j+4)-X(j-1))/5)  ;
                           y_test(j+2)= Y (j-1)+ 3*((Y(j+4)-Y(j-1))/5 );
                           x_test(j+3)= X (j-1)+ 4*((X(j+4)-X(j-1))/5)  ;
                           y_test(j+3)= Y (j-1)+ 4*((Y(j+4)-Y(j-1))/5 );
                       end
                   end
               end
           end
% 
% else
%                x_test(j)=X(j);
%                y_test(j)=Y(j);
%            end
        end
%         x_test(max(size(X)))=X(max(size(X)));
%         y_test(max(size(X)))=Y(max(size(X)));

end
