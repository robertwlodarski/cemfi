% -----------------------------------------------------------------
% Hours Worked - for the US
% -----------------------------------------------------------------
% This is based on Ohanian and Raffo (2012)
% CONTENT
% 1 - Loading the data series. 
% 2 - Removing outliers from the ILO series (Iglewicz & Hoaglin, 1993)
% 3 - Econometric model of the level of the official series
% 4 - Adjusting the quarterly indicators so that they match the official
% series. 
% 5 - Comparing the series with the OR dataset covering up to 2014. 

% -----------------------------------------------------------------
% SECTION 1: LOADING DATA. 
% -----------------------------------------------------------------
% Datasets
% (I) ILS Average Monthly Hours Worked (1994Q1 to 2019Q4).*
% (II) BLS Average Hours Per Week Worked (1970:1 to 2019:12)
% (III) TED Hours Worked Per Year (1950-2023)
% (IV) Ohanian and Raffo (2012) Dataset Obtained From Schuller
% (1960:Q1-2014) - Only US Annual Hours Worked 
% *: The most recent obervations are at the top. 

%Step 1: ILS Quarterly Data. 
ils_data = xlsread("ils_hours_worked_excel.xlsx");
ils_data = flip(ils_data);
%Data Vector.
startDate = datetime('1994-03-01');
endDate = datetime('2019-12-01');
quarters_ils = (startDate:calmonths(3):endDate)';

%Step 2: BLS Quarterly Data.
bls_data = xlsread("bls_hours_per_week_excel.xlsx"); 
bls_data = bls_data(:,2:end)*4.3;%It's given in hours per week. I multiply for hours per month.
bls_vector = reshape(bls_data', size(bls_data,1)*size(bls_data,2),1);
%Panel data (rows: year, columns: month)
startDate = datetime('1970-03-01');
endDate = datetime('2019-12-01');
quarters_bls = (startDate:calmonths(3):endDate)';
%I sum the monthly observations. 
bls_quarterly = zeros(length(quarters_bls),1);
for i=1:1:length(bls_quarterly)
    bls_quarterly(i,1)=sum(bls_vector((3*i-2):(3*i),1));
end 

%Step 3: TED Data.
ted_data = xlsread("ted_1950_2023_excel.xlsx"); 
startDate = datetime('1950-01-01');
endDate = datetime('2023-01-01');
years_ted = (startDate:calmonths(12):endDate)';
ted_data = ted_data(21:70,1);

%Step 4: The author's dataset. 
authors = xlsread("from_schuller_just_us_excel.xlsx"); 
startDate = datetime('1960-03-01');
endDate = datetime('2013-12-01');
years_author = (startDate:calmonths(3):endDate)';
authors = authors(41:end,1)

% -----------------------------------------------------------------
% SECTION 2: REMOVING OUTLIES FROM THE ILO SERIES (IGLEWICZ & HOAGLIN, 1993) 
% -----------------------------------------------------------------
%Step 1: Modified Z-score statistics.
numerator_M = ils_data - median(ils_data);
denominator_M = median(abs(numerator_M));
M_t = 0.6745*(numerator_M/denominator_M);
M_t = abs(M_t);

%Step 2: Identifying outliers. 
outlier = zeros(length(M_t),1);
for i=1:1:length(outlier)
    if M_t(i,1) > 3.48
    outlier(i,1)=1;
    else
    outlier(i,1)=0;
    end
end
%No data are removed from this series. 
%The authors use longer series that has only 4 outliers.
%This makes my results understandable. 

% -----------------------------------------------------------------
% SECTION 3: ECONOMETRIC MODEL OF THE OFFICIAL SERIES 
% -----------------------------------------------------------------
%Ohanian & Raffo (2012) do not apply this estimation to the official series
%of the US and UK. I (happily) follow their lead.

% -----------------------------------------------------------------
% SECTION 4: ADJUSTING THE QUARTERLY INDICATORS (DENTON, 1971)
% -----------------------------------------------------------------
%I use a matlab function found online.
%Source:  mathworks.com/matlabcentral/fileexchange/47568-denton-benchmarking-method?focused=3833101&tab=function
a=1; %Proportional first difference
int=0; %Initial condition
[adjusted_bls,lambda] = denton(bls_quarterly,ted_data,a,int);
writematrix(adjusted_bls,'adjusted_hours.xlsx')

% -----------------------------------------------------------------
% SECTION 5: COMPARING WITH DATASET FROM SCHULLER
% -----------------------------------------------------------------
%Step 1: Cut the Adjusted to Compare
adjusted_cut =adjusted_bls(1:176,1);

%Step 2: Difference b/n Datasets
%Author uses annualised version.
delta_per_cent = 100*(4*adjusted_cut-authors)./authors;
mean(delta_per_cent)
%What about growth rates?
delta_per_per_cent = diff(log(4*adjusted_cut-authors)) - diff(log(authors));
mean(delta_per_per_cent)


%My results are on average 6.42 per cent larger than the authors. 
%I am satisfied with the result. 
%My BLS dataset is much shorter than the one of the authors.
%They include the 1960s that have much lower unemployment rate. 

%Step 3: Plot.
startDate = datetime('1970-03-01');
endDate = datetime('2013-12-01');
years_author = (startDate:calmonths(3):endDate)';
figure(1)
a1 = plot(years_author, 4*adjusted_cut,'-','Color', [0 76/255 153/255],'LineWidth',2)
hold on
a2 = plot(years_author, authors,"-",'Color', [0 0 0],'LineWidth',2)
recessionplot
hold off
ylabel('Hours')
xlabel('Time')
l=   legend([a1,a2],'My Results','Ohanian and Raffo (2012) Sample')
title("Differences between My Results and That of Ohanian and Raffo (2012)")

