%-------------------------------------------------------------------------------
% REPLICATION - ROBERT JACEK WLODARSKI 
% Paper: Financial Cycles: Characterisation and Real-Time Measurement
% Authors: Schuler, Y. S., Hiebert, P. P., and Peltonen, T. A. (2020)
% Journal: Journal of International Money and Finance
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% SECTION 1: Loading data
%-------------------------------------------------------------------------------
%I use the author's naming convention, especially that there are quite many
%functions used. The following variables are used:
%- dY_fin - quarterly growth rates of credit, house, equity, and bond
%- dY_bus - quarterly growth rates of gdp, consumption, investment, and hours
%- fdY_fin - filtered dY_fin
%- fdY_bus - filtered dY_bus
%- for fdY: filter dY with fn_bpass_all(dY_fin_mu,2,200,0,0,0), where dY_fin_mu is
%zero mean
%-real financial variables derived by:   Y_fin = [[credit house equity]./[CPI CPI CPI] (1./(1+BONDYIELD/100))./[CPI]];

%Step 1: Loading the data.
original_data = xlsread("us_excel.xlsx") %column 2: pol_rate, column 3: cpi 
%credit	House	capital	consumption	cpi	gdp	shares
credit_level_nominal = original_data(:,1);
property_real_level_2010 = original_data(:,2);
capital_real_level_2015 = original_data(:,3);
consumption_real_level_2015 = original_data(:,4);
cpi_growth = original_data(:,5);%It's annualised.
gdp_level_2015 = original_data(:,6);
shares_level_real_2015 = original_data(:,7);

%Step 2: Dates.
startDate = datetime('1970-03-01');
endDate = datetime('2019-12-01');
% Generate an array of dates using calmonths
quarters = (startDate:calmonths(3):endDate)';

%Step 3: Expressing real values for credit & property prices. 
%I create the CPI index. 
cpi_index_2015 = zeros(length(cpi_growth),1);
cpi_index_2015(1,1) = 100;
for i=2:1:length(cpi_growth)
    cpi_index_2015(i,1)=cpi_index_2015(i-1,1)*(1+cpi_growth(i-1,1)/400);
end 
%I used 2015:Q1 (obs. 181)
cpi_index_2015 = cpi_index_2015/cpi_index_2015(181,1)*100;%Scale 2015:Q1=100.
%Expressing the variables. 
credit_real_2015 = credit_level_nominal ./ cpi_index_2015;
%For property, they use the average inflation for 2010. For me, it's
%observations 161-4.
property_real_level_2015 = property_real_level_2010*mean(cpi_index_2015(161:164,1),1)/100;

%Step 4: Hours worked. 
hours_level = xlsread("adjusted_hours.xlsx") 

%Step 5: Corporate yields. 
corp =  xlsread("moody_corporate_yields_monthly_excel.xlsx");
corporate_yields_monthly = corp(:,2);
%For the corporate yields, I take values at the end of each quarter.
corporate_yields_nominal = zeros(length(quarters),1);
for i=1:1:length(quarters)
    corporate_yields_nominal(i,1)=corporate_yields_monthly(3*i,1);
end 
%I need to turn them into a real rather than nominal variable. 
corporate_yields_real = 100*(1./(1+corporate_yields_nominal/100))./cpi_index_2015

%Step 6: Setting the variables.
dY_fin = 100*[diff(log(credit_real_2015)) diff(log(property_real_level_2015)) diff(log(shares_level_real_2015)) diff(log(corporate_yields_real))];
dY_bus = 100*[diff(log(gdp_level_2015)) diff(log(consumption_real_level_2015)) diff(log(capital_real_level_2015)) diff(log(hours_level))];
dY_fin_mu = dY_fin-mean(dY_fin);
%mean(dY_fin_mu(:,1)) <- sanity check
dY_bus_mu = dY_bus-mean(dY_bus);
fdY_fin = fn_bpass_all(dY_fin_mu,2,200,0,0,0);
fdY_bus = fn_bpass_all(dY_bus_mu,2,200,0,0,0);

%-------------------------------------------------------------------------------
%SECTION 2: POWER COHESION
%-------------------------------------------------------------------------------
%relevant output:
% T - time dimension
% date_us - cell with dates
% n_co - number of countries
%
% table_years (n_co x 9)    - for FCycle (broad, column 1-3), FCycle (narrow, column 4-6), BCycle (column 7-9): 
%                             [frequency window left, peak, frequency
%                             window right], similar to Table 3 in ESRB paper
%
% FCycle (T x n_co) - band pass filtered financial cycles (broad)
%
%This is the main real-time indicator:
% FCycle_rt (T x n_co) - real time ma filtered financial cycles (broad)
%
% FCycle_n (T x n_co) - band pass filtered financial cycles (narrow)  
% FCycle_n_rt (T x n_co) - real time ma filtered financial cycles (narrow)  
% BCycle (T x n_co) - band pass filtered business cycles   
% BCycle_rt (T x n_co) - real time ma filtered business cycles  
%
% FCycle_unsmooth (T x n_co) - unfiltered financial cycles (broad) (data pre-treated)  
% FCycle_unsmooth_rt (T x n_co) - real time unfiltered financial cycles (broad) (data not pre-treated)  
% FCycle_unsmooth_n (T x n_co) - unfiltered financial cycles (narrow)  (data pre-treated)
% FCycle_unsmooth_n_rt (T x n_co) - real time unfiltered financial cycles (narrow) (data not pre-treated) 
% BCycle_unsmooth (T x n_co) - unfiltered business cycles  (data pre-treated)
% BCycle_unsmooth_rt (T x n_co) - real time unfiltered business cycles  (data not pre-treated)
% 

%Step 1: Prepare Workspace
%1.1:Selecting countries
%1 - CA; 2 - DE, 3 - FR; 4 - IT; 5 - JP; 6 - UK; 7 - US
co = [7];
n_co = length(co);
%1.2: Parameters
q_low = 5;                      % lower bound for spectral density, put 0 for lowest (10=benchmark)
q_high= 200;                    % upper bound for spectral density, put 0 for highest (0 = benchmark)
factor = 8;                     % precision of spectral density estimates
N = 1 ;                         % Number of peaks to select
pct = 0.67;                     % percentage rule for the area of power cohesion to select
corr_rol = 1 ;                  % Composite Cycle: if 1 - rolling corrleations; 0 - linear average
sgn_res = 1;                    % sgn restrictions CISS: emphasises postiviely related movements
idx_fc = [1 1 1 1];
idx_bc = [1 1 1 1];
T=length(quarters)-1;
%1.3: Grid
[grid,step] = fn_grid(q_low,q_high,16);
%1.4: Tables
table_freqband_fcycle = [];
table_freqband_fcycle_n = [];
table_freqband_bcycle = [];
grid_index_total = [];
grid_index_peak = [];
%1.5: Financial Cycle Broad
FCycle = zeros(T,length(co));
FCycle_rt = zeros(T,length(co));
FCycle_unsmooth = zeros(T,length(co));
FCycle_unsmooth_rt = zeros(T,length(co));
series_graph_fc = zeros(T,length(co)*5);
time_varying_weights =  zeros(T,length(co)*4);
%1.6: Financial Cycle Narrow
FCycle_n = zeros(T,length(co));
FCycle_n_rt = zeros(T,length(co));
series_graph_fc_n = zeros(T,length(co)*3);
FCycle_unsmooth_n = zeros(T,length(co));
FCycle_unsmooth_n_rt = zeros(T,length(co));
%1.7: Business Cycle
BCycle = zeros(T,length(co));
BCycle_rt = zeros(T,length(co));
series_graph_bc = zeros(T,length(co)*5);
BCycle_unsmooth = zeros(T,length(co));
BCycle_unsmooth_rt = zeros(T,length(co));
%1.8: Adopting the authors' naming convention
dY_f = fdY_fin;
dY_b = fdY_bus;
dY_f_rt = dY_fin; 
dY_b_rt = dY_bus;
%1.9: Peak space
peak_fc = zeros(1,N);
peak_fc_n = zeros(1,N);
peak_bc = zeros(1,N);

%Step 2: Running the Power Cohesion function.
[pw_coh_fc,f_hats_fc] = pw_cohesion(dY_f,grid,factor,idx_fc);
[pw_coh_fc_n,f_hats_fc_n] = pw_cohesion(dY_f,grid,factor,[1 1 0 0]);
[pw_coh_bc,f_hats_bc] = pw_cohesion(dY_b,grid,factor,idx_bc);

%Step 3: Finding the extrema.
[temp1,~,~] = extrema(pw_coh_fc);
[temp15,~,~] = extrema(pw_coh_fc_n);
[temp2,~,~] = extrema(pw_coh_bc);
    
%Step 4: Identifying the highest N
[~,sortIndex] = sort(pw_coh_fc(temp1(:,1)),'descend');                                                   
maxIndex = sortIndex(1:N); 
peak_fc(1,:) = (pi/2)./grid(temp1(maxIndex,1));
[~,sortIndex] = sort(pw_coh_fc_n(temp15(:,1)),'descend');                                                   
maxIndex = sortIndex(1:N); 
peak_fc_n(1,:) = (pi/2)./grid(temp15(maxIndex,1));
[~,sortIndex] = sort(pw_coh_bc(temp2(:,1)),'descend');                                                   
maxIndex = sortIndex(1:N); 
peak_bc(1,:) = (pi/2)./grid(temp2(maxIndex,1));

%Step 5: Deriving the area and distance interval
[min_distance_fc]=fn_integral_min_dist(pw_coh_fc,peak_fc(1,1),grid,pct);
[min_distance_fc_n]=fn_integral_min_dist(pw_coh_fc_n,peak_fc_n(1,1),grid,pct);
[min_distance_bc]=fn_integral_min_dist(pw_coh_bc,peak_bc(1,1),grid,pct);

%Step 6: Storing the frequency windows and peaks
table_freqband_fcycle = [min_distance_fc(1,3) peak_fc(1,:) min_distance_fc(1,4)];
table_freqband_fcycle_n = [min_distance_fc_n(1,3) peak_fc_n(1,:) min_distance_fc_n(1,4)];
table_freqband_bcycle = [min_distance_bc(1,3) peak_bc(1,:) min_distance_bc(1,4)];

%-------------------------------------------------------------------------------
%SECTION 3: FINANCIAL & BUSINESS CYCLES
%-------------------------------------------------------------------------------
%Step 1: Smoothed values. 
%Note that the functions automatically choose the smoothing parameters
%based on the sample information on the cycles! 
%1.1: Broad Financial Cycle
[series_graph_fc,FCycle,~,~,FCycle_unsmooth,time_varying_weights] = fn_cycle_extract_ecdf_aggr(dY_f,idx_fc,min_distance_fc,corr_rol,sgn_res);
%1.2: Narrow Financial Cycle
[series_graph_fc_n,FCycle_n,~,~,FCycle_unsmooth_n,~] = fn_cycle_extract_ecdf_aggr(dY_f,[1 1 0 0],min_distance_fc_n,0,sgn_res);
%1.3: Business Cycle
[series_graph_bc,BCycle,~,~,BCycle_unsmooth] = fn_cycle_extract_ecdf_aggr(dY_b,idx_bc,min_distance_bc,corr_rol,sgn_res);

%Step 2: Real time financial cycle.
%Like the authors, I use 3 years training sample (12 quarters)
%in fn_fcycle_dy_rt(~,~,~,~,12)
%I later adjust it to higher and lower values to see how this evolves.
%2.1: Broad Financial Cycle
[FCycle_unsmooth_rt,FCycle_rt,dY_ecdf,dY_ecdf_ma,weights] = fn_fcycle_dy_rt(dY_f_rt,[1 1 1 1],1,1,12); %3years training sample
%2.2: Narrow Financial Cycle
[FCycle_unsmooth_n_rt,FCycle_n_rt,~,~,~] = fn_fcycle_dy_rt(dY_f_rt,[1 1 0 0],1,1,12); %3years training sample
%2.3: Business Cycle
[BCycle_unsmooth_rt,BCycle_rt,~,~,~] = fn_fcycle_dy_rt(dY_b_rt,[1 1 1 1],1,1,12); %3years training sample
table_years = [table_freqband_fcycle table_freqband_fcycle_n table_freqband_bcycle];

%-------------------------------------------------------------------------------
%SECTION 4: PLOTTING CROSS SPECTRA
%-------------------------------------------------------------------------------
%Step 1: Preparing the graphs.
freq_low = 2*pi/(2*4);
freq_high= 2*pi/(8*4);
bus_freq = (grid <= freq_low) & (grid >= freq_high);
freq_low_fc =  2*pi/(8*4);
freq_high_fc = 2*pi/(20*4);
fc_freq = (grid <= freq_low_fc) & (grid >= freq_high_fc);
grid_select = grid <= grid(end); %choose window size
      
%Step 1: Naming each variable
Ynam = {'\Deltacr'; '\Deltap_h'; '\Deltap_e'; '\Deltap_b'; '\Deltaq'; '\Deltaco'; '\Deltai'; '\Deltah'}; 

%Step 2: Plotting financial cross spectra   
size_plot = [8/2 8/2];
ii=7;
f_hat_fc_i = f_hats_fc;
figure
hold on
area(grid(grid_select),bus_freq(grid_select)*max(f_hat_fc_i(:)))
colormap([0.69 0.7686 0.87])
area(grid(grid_select),fc_freq(grid_select)*max([f_hat_fc_i(:)]),'FaceColor',[230/255 230/255 250/255])
a1=   plot(grid(grid_select),f_hat_fc_i(grid_select,1) ,'-', 'Color',[0 0 1],'LineWidth',1.5);
a2=   plot(grid(grid_select),f_hat_fc_i(grid_select,2) ,'-','Color',[0 1 0],'LineWidth',1.5);
a3=   plot(grid(grid_select),f_hat_fc_i(grid_select,3) ,'-','Color',[1 0 0],'LineWidth',1.5);
a4=   plot(grid(grid_select),f_hat_fc_i(grid_select,4) ,'-','Color',[1 0 1],'LineWidth',1.5);
a5=   plot(grid(grid_select),f_hat_fc_i(grid_select,5) ,'-','Color',[0 1 1],'LineWidth',1.5);
a6=   plot(grid(grid_select),f_hat_fc_i(grid_select,6) ,'-','Color',[0 0 0],'LineWidth',1.5); 
hold off
lab_x = [grid(18)  grid(56) grid(248)];
set(gca, 'XTick', lab_x , 'XTickLabel', {'20'; '8'; '2'},'FontSize',10)
l=   legend([a1,a2,a3,a4,a5,a6],[Ynam{1} '/' Ynam{2}],[Ynam{1} '/' Ynam{3}],[Ynam{1} '/' Ynam{4}],[Ynam{2} '/' Ynam{3}],[Ynam{2} '/' Ynam{4}],[Ynam{3} '/' Ynam{4}]); %'Location','northoutside','orientation','horizontal'
set(l,'FontSize',10)
legend('boxoff')
ylabel('$|\hat{s}_{x_ix_j}|$','Interpreter','Latex','FontSize',15)
xlabel('Years')
axis tight
set(gcf,'PaperPositionMode','manual')
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[0 0 size_plot]);
figfile = fullfile('figures', ['f_hat_fc_' num2str(ii)]);
%saveas(gcf, figfile , 'epsc'); 

%Step 3: Busiess Cycle Cross Spectra
size_plot = [8/2 8/2];
f_hat_fc_i = f_hats_bc;
figure
hold on
area(grid(grid_select),bus_freq(grid_select)*max(f_hat_fc_i(:)))
colormap([0.69 0.7686 0.87])
area(grid(grid_select),fc_freq(grid_select)*max([f_hat_fc_i(:)]),'FaceColor',[230/255 230/255 250/255])
a1=   plot(grid(grid_select),f_hat_fc_i(grid_select,1) ,'-','Color',[0 0 1],'LineWidth',1.5);
a2=   plot(grid(grid_select),f_hat_fc_i(grid_select,2) ,'-','Color',[0 1 0],'LineWidth',1.5);
a3=   plot(grid(grid_select),f_hat_fc_i(grid_select,3) ,'-','Color',[1 0 0],'LineWidth',1.5);
a4=   plot(grid(grid_select),f_hat_fc_i(grid_select,4) ,'-','Color',[1 0 1],'LineWidth',1.5);
a5=   plot(grid(grid_select),f_hat_fc_i(grid_select,5) ,'-','Color',[0 1 1],'LineWidth',1.5);
a6=   plot(grid(grid_select),f_hat_fc_i(grid_select,6) ,'-','Color',[0 0 0],'LineWidth',1.5);
hold off
lab_x = [grid(18)  grid(56) grid(248)];
set(gca, 'XTick', lab_x , 'XTickLabel', {'20'; '8'; '2'},'FontSize',10)
l=   legend([a1,a2,a3,a4,a5,a6],[Ynam{5} '/' Ynam{6}],[Ynam{5} '/' Ynam{7}],[Ynam{5} '/' Ynam{8}],[Ynam{6} '/' Ynam{7}],[Ynam{6} '/' Ynam{8}],[Ynam{7} '/' Ynam{8}]); %'Location','northoutside','orientation','horizontal'
set(l,'FontSize',10)
legend('boxoff')
ylabel('$|\hat{s}_{x_ix_j}|$','Interpreter','Latex','FontSize',15)
xlabel('Years')
axis tight
set(gcf,'PaperPositionMode','manual')
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[0 0 size_plot]);
figfile = fullfile('figures', ['f_hat_bc_' num2str(ii)]);
%saveas(gcf, figfile , 'epsc'); 

%-------------------------------------------------------------------------------
%SECTION 5: POWER COHESION
%-------------------------------------------------------------------------------


%Step 1: Parameters 
% %print financial and real cycle power cohesion in one graph.
size_plot = [8/2 8/2];
%jip=7;
freq_low = 2*pi/(2*4);
freq_high= 2*pi/(8*4);
bus_freq = (grid <= freq_low) & (grid >= freq_high);
freq_low_fc =  2*pi/(8*4);
freq_high_fc = 2*pi/(20*4);
fc_freq = (grid<= freq_low_fc) & (grid >= freq_high_fc);
grid_s = (grid <= grid(end)) & (grid>= grid(1)); %choose window size

%Step 2: Plotting
figure('Name','Power Cohesion')
hold on
area(grid(grid_s),bus_freq(grid_s)*max([pw_coh_fc(grid_s); pw_coh_bc(grid_s);pw_coh_fc_n(grid_s)]))
colormap([0.69 0.7686 0.87])
area(grid(grid_s),fc_freq(grid_s)*max([[pw_coh_fc(grid_s); pw_coh_bc(grid_s);pw_coh_fc_n(grid_s)]]),'FaceColor',[230/255 230/255 250/255])
a1=plot(grid(grid_s),pw_coh_fc(grid_s,1),'-','Color', [0 0 0],'LineWidth',3); %yearly fc
a2=plot(grid(grid_s),pw_coh_fc_n(grid_s,1),'-','Color', [0 76/255 153/255],'LineWidth',3); %yearly fc
a3=plot(grid(grid_s),pw_coh_bc(grid_s,1),'-','Color',[0.7 0 0],'LineWidth',3); %yearly bc
hold off
lab_x = [grid(18)  grid(56) grid(248)];
set(gca, 'XTick', lab_x , 'XTickLabel', {'20'; '8'; '2'},'FontSize',10)
title("United States",'FontSize', 11);
l=   legend([a1,a2,a3],'Broad financial cycle','Narrow financial cycle','Business cycle');%'Location','northoutside','orientation','horizontal'
set(l,'FontSize',10)
legend('boxoff')
ylabel('PCoh','FontSize',11)
xlabel('Years')
axis tight
set(gcf,'PaperPositionMode','manual')
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[0 0 size_plot]);   
figfile = fullfile('figures\US');
%saveas(gca, figfile , 'epsc'); 
      
%-------------------------------------------------------------------------------
%SECTION 6: OTHER PLOTS AND CALCULATIONS
%-------------------------------------------------------------------------------
%Step 1: Financial Cycles Plot.
%1.1.: Dates
startDate = datetime('1970-06-01');
endDate = datetime('2019-12-01');
% Generate an array of dates using calmonths
quarters = (startDate:calmonths(3):endDate)';
%1.2.: Elements of the plot
median = 0.5*ones(length(quarters),1);
%1.3: Broad Financial Cycles
figure(1001)
a1 = plot(quarters, FCycle,'-','Color', [0 0 0],'LineWidth',3)
hold on
a2 = plot(quarters, FCycle_unsmooth,"--",'Color', [0 0 0],'LineWidth',1)
plot(quarters, median)
recessionplot
hold off
ylabel('Deviation from the Historical Median Growth')
xlabel('Time')
l=   legend([a1,a2],'Filtered Broad Financial Cycle','Unfiltered Broad Financial Cycle')
title("Broad Financial Cycle")
%1.3: Narrow Financial Cycle
figure(1002)
a1 = plot(quarters, FCycle_n,'-','Color', [0 0 0],'LineWidth',3)
hold on
a2 = plot(quarters, FCycle_unsmooth_n,"--",'Color', [0 0 0],'LineWidth',1)
plot(quarters, median)
recessionplot
hold off
ylabel('Deviation from the Historical Median Growth')
xlabel('Time')
l=   legend([a1,a2],'Filtered Narrow Financial Cycle','Unfiltered Narrow Financial Cycle')
title("Narrow Financial Cycle")
%1.4: Comparing the types of cycles.
figure(1003)
a1 = plot(quarters, FCycle_n,'-','Color', [0 76/255 153/255],'LineWidth',3)
hold on
a2 = plot(quarters, FCycle,"-",'Color', [0 0 0],'LineWidth',3)
a3 = plot(quarters, BCycle,"-",'Color', [0.7 0 0],'LineWidth',3)
plot(quarters, median)
recessionplot
hold off
ylabel('Deviation from the Historical Median Growth')
xlabel('Time')
l=   legend([a1,a2, a3],'Narrow Financial Cycle','Broad Financial Cycle','Business Cycle')
title("Filtered Cycles")
%1.5: Filtered financial variables. %APPENDIX C ANALYSIS
fdY_fin_f2 = fn_bpass_all(dY_fin_mu,4.9*4,33.6*4,0,0,0);
%- dY_fin_f2 - detrended quarterly growth rates of credit, house, equity, and bond
figure(1004)
a1 = plot(quarters,fdY_fin_f2(:,1),"-",'Color',[0 0 1],'LineWidth',2)
hold on 
a2 = plot(quarters,fdY_fin_f2(:,2),"-",'Color',[1 0 0],'LineWidth',2)
a3 = plot(quarters,fdY_fin_f2(:,3),"-",'Color',[0 1 0],'LineWidth',2)
a4 = plot(quarters,fdY_fin_f2(:,4),"-",'Color',[0.5 0 0.5],'LineWidth',2)
recessionplot
hold off
ylabel('Per Cent')
xlabel('Time')
l=   legend([a1,a2, a3, a4],'Credit','House Prices','Equity Prices','Bond Prices')
title("Cyclical Growth Rates")
%1.6: Real Time Cycles.
figure(1006)
a1 = plot(quarters, FCycle_n_rt,'-','Color', [0 76/255 153/255],'LineWidth',3)
hold on
a2 = plot(quarters, FCycle_rt,"-",'Color', [0 0 0],'LineWidth',3)
a3 = plot(quarters, BCycle_rt,"-",'Color', [0.7 0 0],'LineWidth',3)
plot(quarters, median)
recessionplot
hold off
ylabel('Deviation from the Historical Median Growth')
xlabel('Time')
l=   legend([a1,a2, a3],'Narrow Financial Cycle','Broad Financial Cycle','Business Cycle')
title("Real Time Filtered Cycles")

%-------------------------------------------------------------------------------
%SECTION 7: COINCIDENT & EARLY WARNING INDICATORS - PSEUDO OUT-OF-SAMPLE
%(PART 1: CALCULATIONS & GRAPHS)
%-------------------------------------------------------------------------------
%Step 1: Adjusting the time (I use only 1970:Q1 through 1999:Q4)
%1-119
experiment_dY_f_rt=dY_f_rt;
experiment_dY_b_rt=dY_b_rt;

%Step 2: Real time financial cycle.
%Like the authors, I use 10 years training sample (39 quarters; one quarter is lost)
%I later adjust it to higher and lower values to see how this evolves.
%2.1: Broad Financial Cycle
[FCycle_unsmooth_rt_ex,FCycle_rt_ex,dY_ecdf_ex,dY_ecdf_ma_ex,weights_ex] = fn_fcycle_dy_rt(experiment_dY_f_rt,[1 1 1 1],1,1,39); %10years training sample
%2.2: Narrow Financial Cycle
[FCycle_unsmooth_n_rt_ex,FCycle_n_rt_ex,~,~,~] = fn_fcycle_dy_rt(experiment_dY_f_rt,[1 1 0 0],1,1,39); %10years training sample
%2.3: Business Cycle
[BCycle_unsmooth_rt_ex,BCycle_rt_ex,~,~,~] = fn_fcycle_dy_rt(experiment_dY_b_rt,[1 1 1 1],1,1,39); %10years training sample

%Step 3: Variables for regression. 
% Collected variables. 
%- dY_fin - quarterly growth rates of credit, house, equity, and bond
cr_ex = experiment_dY_f_rt(40:119,1);
cr_pp_ex = experiment_dY_f_rt(40:119,1:2);
all_fin_ex = experiment_dY_f_rt(40:119,:);
broad_fin_cycle=FCycle_unsmooth_rt_ex(40:119,1);
narrow_fin_cycle=FCycle_unsmooth_n_rt_ex(40:119,1);
business_cycle=BCycle_unsmooth_rt_ex(40:119,1);
%Lagged version (I use only one lag due to a small number of observations)
cr_ex_l1= lagmatrix(cr_ex,1);
cr_pp_ex_l1=[lagmatrix(cr_pp_ex(:,1),1) lagmatrix(cr_pp_ex(:,1),1)];
all_fin_ex_l1=[lagmatrix(all_fin_ex(:,1),1) lagmatrix(all_fin_ex(:,1),1) lagmatrix(all_fin_ex(:,1),1) lagmatrix(all_fin_ex(:,1),1)];
broad_fin_cycle_l1=lagmatrix(broad_fin_cycle,1);
narrow_fin_cycle_l1=lagmatrix(narrow_fin_cycle,1);
business_cycle_l1=lagmatrix(business_cycle,1);

%Step 4: Preparing the experiment.
startDate = datetime('1980-03-01');
endDate = datetime('1999-12-01');
quarters_experiment = (startDate:calmonths(3):endDate)';
onset_dummy=ones(length(quarters_experiment),1);
onset_dummy(31,1) = 2;
early_dummy=zeros(length(quarters_experiment),1)+1;
early_dummy(27:30,1) = 2;
super_early_dummy = zeros(length(quarters_experiment),1)+1;
super_early_dummy(22:29,1)=2;

%Step 4: Regressions for narrow financial cycle (1).
%4.1. Onset, one lag. 
X=[narrow_fin_cycle(2:end,1) narrow_fin_cycle_l1(2:end,1)];
X_warning = X([1:33,40:79],:)
Y=onset_dummy(2:end,1);
onset_regression_1 = mnrfit(X,Y);
%4.2. Early warning, one lag. 
%I follow Anundsen et al. (2016) in removing the onset of the crisis as well as 
%6 quarters thereafter from the samples used in the vulnerability-centred samples
Y=early_dummy(2:end,1);
Y=Y([1:33,40:79]);
early_regression_1 = mnrfit(X_warning,Y);
%4.3. Super early warning, one lag. 
Y=super_early_dummy(2:end,1);
Y=Y([1:33,40:79]);%Anundsen et al. (2016) 
super_early_regression_1 = mnrfit(X_warning,Y);

%Step 5: Regressions for broad financial cycle (2).
%5.1. Onset, one lag. 
X=[broad_fin_cycle(2:end,1) broad_fin_cycle_l1(2:end,1)];
X_warning = X([1:33,40:79],:);
Y=onset_dummy(2:end,1);
onset_regression_2 = mnrfit(X,Y);
%5.2. Early warning, one lag. 
Y=early_dummy(2:end,1);
Y=Y([1:33,40:79]);%Anundsen et al. (2016)
early_regression_2 = mnrfit(X_warning,Y);
%5.3. Super early warning, one lag. 
Y=super_early_dummy(2:end,1);
Y=Y([1:33,40:79]);%Anundsen et al. (2016) 
super_early_regression_2 = mnrfit(X_warning,Y);

%Step 6: Regressions for business cycle (3).
%6.1. Onset, one lag. 
X=[business_cycle(2:end,1) business_cycle_l1(2:end,1)];
X_warning = X([1:33,40:79],:);
Y=onset_dummy(2:end,1);
onset_regression_3 = mnrfit(X,Y);
%6.2. Early warning, one lag. 
Y=early_dummy(2:end,1);
Y=Y([1:33,40:79]);%Anundsen et al. (2016)
early_regression_3 = mnrfit(X_warning,Y);
%6.3. Super early warning, one lag. 
Y=super_early_dummy(2:end,1);
Y=Y([1:33,40:79]);%Anundsen et al. (2016) 
super_early_regression_3 = mnrfit(X_warning,Y);

%Step 7: Regressions for credit (4).
X=[cr_ex(2:end,1) cr_ex_l1(2:end,1)];
X_warning = X([1:33,40:79],:);
Y=onset_dummy(2:end,1);
onset_regression_4 = mnrfit(X,Y);
%6.2. Early warning, one lag. 
Y=early_dummy(2:end,1);
Y=Y([1:33,40:79]);%Anundsen et al. (2016)
early_regression_4 = mnrfit(X_warning,Y);
%6.3. Super early warning, one lag. 
Y=super_early_dummy(2:end,1);
Y=Y([1:33,40:79]);%Anundsen et al. (2016) 
super_early_regression_4 = mnrfit(X_warning,Y);

% %Step 8: Regressions for credit & housing prices (5).
% X=[cr_pp_ex(2:end,:) cr_pp_ex_l1(2:end,:)];
% X_warning = X([1:33,40:79],:);
% Y=onset_dummy(2:end,1);
% onset_regression_5 = mnrfit(X,Y);
% %6.2. Early warning, one lag. 
% Y=early_dummy(2:end,1);
% Y=Y([1:33,40:79]);%Anundsen et al. (2016)
% early_regression_5 = mnrfit(X_warning,Y);
% %6.3. Super early warning, one lag. 
% Y=super_early_dummy(2:end,1);
% Y=Y([1:33,40:79]);%Anundsen et al. (2016) 
% super_early_regression_5 = mnrfit(X_warning,Y);
% 

%Step 9: Fitting. 
%Variables
f_cr_ex = experiment_dY_f_rt(120:end,1);
f_cr_pp_ex = experiment_dY_f_rt(120:end,1);
f_all_fin_ex = experiment_dY_f_rt(120:end,1);
f_broad_fin_cycle=FCycle_unsmooth_rt_ex(120:end,1);
f_narrow_fin_cycle=FCycle_unsmooth_n_rt_ex(120:end,1);
f_business_cycle=BCycle_unsmooth_rt_ex(120:end,1);
%Lagged version (I use only one lag due to a small number of observations)
f_cr_ex_l1= lagmatrix(f_cr_ex,1);
f_cr_pp_ex_l1=lagmatrix(f_cr_pp_ex,1);
f_all_fin_ex_l1=lagmatrix(f_all_fin_ex,1);
f_broad_fin_cycle_l1=lagmatrix(f_broad_fin_cycle,1);
f_narrow_fin_cycle_l1=lagmatrix(f_narrow_fin_cycle,1);
f_business_cycle_l1=lagmatrix(f_business_cycle,1);
%9.1. Onset regressions
v_ones = ones(79,1);
X_obs = [v_ones f_narrow_fin_cycle(2:end,1) f_narrow_fin_cycle_l1(2:end,1)];
fitted_onset_regression_1= X_obs*onset_regression_1;
fitted_onset_regression_1=1-1./(1+exp(-fitted_onset_regression_1));
X_obs = [v_ones f_broad_fin_cycle(2:end,1) f_broad_fin_cycle_l1(2:end,1)];
fitted_onset_regression_2= X_obs*onset_regression_2;
fitted_onset_regression_2=1-1./(1+exp(-fitted_onset_regression_2));
X_obs = [v_ones f_business_cycle(2:end,1) f_business_cycle_l1(2:end,1)];
fitted_onset_regression_3= X_obs*onset_regression_3;
fitted_onset_regression_3=1-1./(1+exp(-fitted_onset_regression_3));
X_obs = [v_ones f_cr_ex(2:end,1) f_cr_ex_l1(2:end,1)];
fitted_onset_regression_4= X_obs*onset_regression_4;
fitted_onset_regression_4=1-1./(1+exp(-fitted_onset_regression_4));
%9.2. Early regressions
v_ones = ones(79,1);
X_obs = [v_ones f_narrow_fin_cycle(2:end,1) f_narrow_fin_cycle_l1(2:end,1)];
fitted_early_regression_1= X_obs*early_regression_1;
fitted_early_regression_1=1-1./(1+exp(-fitted_early_regression_1));
X_obs = [v_ones f_broad_fin_cycle(2:end,1) f_broad_fin_cycle_l1(2:end,1)];
fitted_early_regression_2= X_obs*early_regression_2;
fitted_early_regression_2=1-1./(1+exp(-fitted_early_regression_2));
X_obs = [v_ones f_business_cycle(2:end,1) f_business_cycle_l1(2:end,1)];
fitted_early_regression_3= X_obs*early_regression_3;
fitted_early_regression_3=1-1./(1+exp(-fitted_early_regression_3));
X_obs = [v_ones f_cr_ex(2:end,1) f_cr_ex_l1(2:end,1)];
fitted_early_regression_4= X_obs*early_regression_4;
fitted_early_regression_4=1-1./(1+exp(-fitted_early_regression_4));
%9.3. Early regressions
v_ones = ones(79,1);
X_obs = [v_ones f_narrow_fin_cycle(2:end,1) f_narrow_fin_cycle_l1(2:end,1)];
fitted_super_early_regression_1= X_obs*super_early_regression_1;
fitted_super_early_regression_1=1-1./(1+exp(-fitted_super_early_regression_1));
X_obs = [v_ones f_broad_fin_cycle(2:end,1) f_broad_fin_cycle_l1(2:end,1)];
fitted_super_early_regression_2= X_obs*super_early_regression_2;
fitted_super_early_regression_2=1-1./(1+exp(-fitted_super_early_regression_2));
X_obs = [v_ones f_business_cycle(2:end,1) f_business_cycle_l1(2:end,1)];
fitted_super_early_regression_3= X_obs*super_early_regression_3;
fitted_super_early_regression_3=1-1./(1+exp(-fitted_super_early_regression_3));
X_obs = [v_ones f_cr_ex(2:end,1) f_cr_ex_l1(2:end,1)];
fitted_super_early_regression_4= X_obs*super_early_regression_4;
fitted_super_early_regression_4=1-1./(1+exp(-fitted_super_early_regression_4));

%Step 10: Plotting the fitted values. 
startDate = datetime('2000-06-01');
endDate = datetime('2019-12-01');
quarters_fitted = (startDate:calmonths(3):endDate)';
%10.1: Onset of a banking crisis dummy 
figure(1007)
a1 = plot(quarters_fitted, fitted_onset_regression_1,'-','Color', [0 76/255 153/255],'LineWidth',2)
hold on
a2 = plot(quarters_fitted, fitted_onset_regression_2,"-",'Color', [0 0 0],'LineWidth',2)
a3 = plot(quarters_fitted, fitted_onset_regression_3,"-",'Color', [0.7 0 0],'LineWidth',2)
a4 =plot(quarters_fitted, fitted_onset_regression_4,"-",'Color', [0.7 0 1],'LineWidth',2)
recessionplot
hold off
ylabel('Probability')
xlabel('Time')
l=   legend([a1,a2, a3, a4],'Narrow Financial Cycle','Broad Financial Cycle','Business Cycle','Total Credit')
title("Fitted Probability of the Onset of a Systemic Banking Crisis")
%10.2: Early warning indicator of a banking crisis dummy 
figure(1008)
a1 = plot(quarters_fitted, fitted_early_regression_1,'-','Color', [0 76/255 153/255],'LineWidth',2)
hold on
a2 = plot(quarters_fitted, fitted_early_regression_2,"-",'Color', [0 0 0],'LineWidth',2)
a3 = plot(quarters_fitted, fitted_early_regression_3,"-",'Color', [0.7 0 0],'LineWidth',2)
a4 =plot(quarters_fitted, fitted_early_regression_4,"-",'Color', [0.7 0 1],'LineWidth',2)
recessionplot
hold off
ylabel('Probability')
xlabel('Time')
l=   legend([a1,a2, a3, a4],'Narrow Financial Cycle','Broad Financial Cycle','Business Cycle','Total Credit')
title("Fitted Early Indicator of a Systemic Banking Crisis (1-5 Quarters)")
%10.3: Early warning indicator of a banking crisis dummy 
figure(1009)
a1 = plot(quarters_fitted, fitted_super_early_regression_1,'-','Color', [0 76/255 153/255],'LineWidth',2)
hold on
a2 = plot(quarters_fitted, fitted_super_early_regression_2,"-",'Color', [0 0 0],'LineWidth',2)
a3 = plot(quarters_fitted, fitted_super_early_regression_3,"-",'Color', [0.7 0 0],'LineWidth',2)
a4 =plot(quarters_fitted, fitted_super_early_regression_4,"-",'Color', [0.7 0 1],'LineWidth',2)
recessionplot
hold off
ylabel('Probability')
xlabel('Time')
l=   legend([a1,a2, a3, a4],'Narrow Financial Cycle','Broad Financial Cycle','Business Cycle','Total Credit')
title("Fitted Very Early Indicator of a Systemic Banking Crisis (5-12 Quarters)")

%-------------------------------------------------------------------------------
%SECTION 8: COINCIDENT & EARLY WARNING INDICATORS - PSEUDO OUT-OF-SAMPLE
%(PART 2: SUMMARY STATISTICS)
%-------------------------------------------------------------------------------

%Step 1: Optimal threshold. 
% onset_vectors = [fitted_onset_regression_1 fitted_onset_regression_2 fitted_onset_regression_3 fitted_onset_regression_4];
% onset_thresholds = zeros(4,1);
% for i=1:1:4
%     z_function = @(x) sum_of_errors(x,onset_vectors(:,i));
%     xyz=min(z_function,[0,1]);
%     onset_thresholds(i,1)=xyz(1,1);
% end 
% % z_function = @(x) sum_of_errors(x,fitted_onset_regression_1)
% % x = fminsearch(z,[0,1])
% early_vectors = [fitted_early_regression_1 fitted_early_regression_2 fitted_early_regression_3 fitted_early_regression_4];
% early_thresholds = zeros(4,1);
% for i=1:1:4
%     z_function = @(x) sum_of_errors(x,early_vectors(:,i));
%     xyz=fminsearch(z_function,[0,1]);
%     early_thresholds(i,1)=xyz(1,1);
% end 
% super_early_vectors = [fitted_super_early_regression_1 fitted_super_early_regression_2 fitted_super_early_regression_3 fitted_super_early_regression_4];
% super_early_thresholds = zeros(4,1);
% for i=1:1:4
%     z_function = @(x) sum_of_errors(x,super_early_vectors(:,i));
%     xyz=fminsearch(z_function,[0,1]);
%     super_early_thresholds(i,1)=xyz(1,1);
% end 
% z_function = @(x) sum_of_errors(x,super_early_vectors(:,1));
% fplot(z_function,[0,1])
%Values are very weird. I will just use 0.3 for everything. 

%Step 2: Calculating the statistics. 
%Predicted values for threshold 0.5
%Vectors.
onset_vectors = [fitted_onset_regression_1 fitted_onset_regression_2 fitted_onset_regression_3 fitted_onset_regression_4];
early_vectors = [fitted_early_regression_1 fitted_early_regression_2 fitted_early_regression_3 fitted_early_regression_4];
super_early_vectors = [fitted_super_early_regression_1 fitted_super_early_regression_2 fitted_super_early_regression_3 fitted_super_early_regression_4];
%Real events.
true_onset = zeros(79,1);
true_onset(30,1)= 1;
true_early = zeros(79,1);
true_early(26:29,1) = 1;
true_very_early = zeros(79,1);
true_very_early(18:25) = 1
%Fitted events.
fitted_onset_values = (onset_vectors>0.3);
fitted_early_values = (early_vectors>0.3);
fitted_very_early_values = (super_early_vectors>0.3);
% Calculating key values
TP_onset = zeros(4,1);
FP_onset= zeros(4,1);
TN_onset= zeros(4,1);
FN_onset= zeros(4,1);
T1_onset = zeros(4,1);
T2_onset = zeros(4,1);
usefulness_onset = zeros(4,1);
noise_to_signal_onset = zeros(4,1);
auc_onset= zeros(4,1);
for i=1:1:4
    [TP_onset(i,1), FP_onset(i,1), TN_onset(i,1), FN_onset(i,1)] = calError(true_onset, fitted_onset_values(:,i));
    T1_onset(i,1)= FN_onset(i,1)/(FN_onset(i,1)+TP_onset(i,1));
    T2_onset(i,1)= FP_onset(i,1)/(FP_onset(i,1)+TN_onset(i,1));
    usefulness_onset(i,1) = (0.5-T1_onset(i,1)-T2_onset(i,1))/0.5;
    noise_to_signal_onset(i,1) = T2_onset(i,1)/(T1_onset(i,1));
    auc_onset(i,1)=TP_onset(i,1)/(FP_onset(i,1)+TP_onset(i,1))
end
TP_early = zeros(4,1);
FP_early= zeros(4,1);
TN_early= zeros(4,1);
FN_early= zeros(4,1);
T1_early = zeros(4,1);
T2_early = zeros(4,1);
usefulness_early = zeros(4,1);
noise_to_signal_early = zeros(4,1);
auc_onset= zeros(4,1);
for i=1:1:4
    [TP_early(i,1), FP_early(i,1), TN_early(i,1), FN_early(i,1)] = calError(true_early, fitted_early_values(:,i));
    T1_early(i,1)= FN_early(i,1)/(FN_early(i,1)+TP_early(i,1));
    T2_early(i,1)= FP_early(i,1)/(FP_early(i,1)+TN_early(i,1));
    usefulness_early(i,1) = (0.5-T1_early(i,1)-T2_early(i,1))/0.5;
    noise_to_signal_early(i,1) = T2_early(i,1)/(T1_early(i,1));
    auc_early(i,1)=TP_early(i,1)/(FP_early(i,1)+TP_early(i,1))
end
TP_very_early = zeros(4,1);
FP_very_early= zeros(4,1);
TN_very_early= zeros(4,1);
FN_very_early= zeros(4,1);
T1_very_early = zeros(4,1);
T2_very_early = zeros(4,1);
usefulness_very_early = zeros(4,1);
noise_to_signal_very_early = zeros(4,1);
auc_onset= zeros(4,1);
for i=1:1:4
    [TP_very_early(i,1), FP_very_early(i,1), TN_very_early(i,1), FN_very_early(i,1)] = calError(true_very_early, fitted_very_early_values(:,i));
    T1_very_early(i,1)= FN_very_early(i,1)/(FN_very_early(i,1)+TP_very_early(i,1));
    T2_very_early(i,1)= FP_very_early(i,1)/(FP_very_early(i,1)+TN_very_early(i,1));
    usefulness_very_early(i,1) = (0.5-T1_very_early(i,1)-T2_very_early(i,1))/0.5;
    noise_to_signal_very_early(i,1) = T2_very_early(i,1)/(T1_very_early(i,1));
    auc_very_early(i,1)=TP_very_early(i,1)/(FP_very_early(i,1)+TP_very_early(i,1))
end

%Step 3: Recreating Schuller et al (2020)'s Table.
row_onset = [TP_onset FP_onset TN_onset FN_onset T1_onset T2_onset usefulness_onset noise_to_signal_onset auc_onset];
row_early = [TP_early FP_early TN_early FN_early T1_early T2_early usefulness_early noise_to_signal_early auc_early];
row_very_early = [TP_very_early FP_very_early TN_very_early FN_very_early T1_very_early T2_very_early usefulness_very_early noise_to_signal_very_early auc_very_early];
big_table=[row_onset;row_early;row_very_early]
