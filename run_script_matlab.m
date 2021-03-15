clear
close all

%This script can be used to calculculate combined reference metabolite
%ratios (R_CRef) or concentrations (C_CRef) from single subject data. See:
%Quantification of Magnetic Resonance Spectroscopy data using
% a combined reference: Application in typically developing infants
%Calculated values of R_CRef and C_CRef are compared with normative data to
%estimate a percentile for each metabolite. The percentile calculation is
%performed by comparing R_CRef or C_CRef with the expected value for that
%age, and by assuming a normal distribution. The standard deviation of the
%assumed normal distribution is calculated from the normative data for each
%metabolite by first removing post-menstrual age dependnece using a second
%order polynomial

%ssd is a structure that holds all variables related to the single-subject
%analysis.
%ald is a structure that holds data from the normative data described in
%the paper

%ssd.met_names lists the metabolites to be analyzed.
%metabolies should be ordered to be consistent with ordering in the normative data, which is: 
%'GABA','tNAA','NAA','tCr','Cho','Ins','Glx','Glu','GSH','Tau'

ssd.met_names = {'GABA','tNAA','NAA','tCr','Cho','Ins','Glx','Glu','GSH','Tau'}

%ssd.sgn contains signal values for each of the listed metabolites. Signal
%values should already be adjusted for number of protons, such that signal
%values are is proportional to the number of molecules. Data from LCModel
%have already received this adjustment. If you don't have your own data
%there is an option below,after the normative data are loaded into ssd, to
%obtain the test subject from the normative data.
%Default data: from infant #1 of normative data set.
ssd.sgn = [2.1468   10.6573    8.9165    7.3482    2.5725    8.2330    6.7045    6.0122    2.8335    2.8844];

%the reference metabolites.   
ssd.ref_metabolites = {'GABA','tNAA','tCr','Cho','Ins','Glx','GSH'};
%ssd.ref_metabolites = {'tCr'};

%PMA is chronological age + gestational age, in days.
%Default data: from infant #1 of normative data set.
ssd.PMA = 399;

%calculation options.
ssd.is_weighting = 'y'; %using weighting factors (normalization by weighted SD values). if no, all weighting factors are one.  y if yes, n if no, default: y.
ssd.is_self = 'y';  %allow self-referencing? y if yes, n if no, default: n. 
if length(ssd.ref_metabolites) == 1
    ssd.is_self = 'y'; 
end

%results directory
results_dir = fullfile('results','test1'); %directory for results
mkdir(results_dir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load the normative data:
load(fullfile('InfantData','App_data.mat'))

%identify the indices of the metabolites within the ald
[ssd.aldi, ordno] = ismember(ald.met_names,ssd.met_names);

%data checks:

%metabolies should be ordered to be consistent with ordering in the normative data, which is: 
%{'GABA'}    {'tNAA'}    {'NAA'}    {'tCr'}    {'Cho'}    {'Ins'}    {'Glx'}    {'Glu'}    {'GSH'}    {'Tau'}
if nnz(find(sign(diff(find(ordno)))==1)) ~= length(ordno)-1
    error('metabolite order error')
end

%all of the reference metabolites must be included in the metabolites to be
%analyzed, ssd.met_names
if nnz(ismember(ssd.ref_metabolites,ssd.met_names )) ~= length(ssd.ref_metabolites)
   error('include all ref metabolites in list of mets to be analyzed')  
end

%bring the needed normative data from the data structure ald to ssd.
%only the metabolites to be analyzed are brought over.
ssd.metwSD = ald.metwSD(ssd.aldi); %standard deviation values used to calculate weighting factors, calculated from PMA-corrected water-scaled data from the normative data
ssd.metw = ald.metw(:,ssd.aldi); %water-scaled data from normative data
ssd.metfw = ald.metfw(:,ssd.aldi); % concentration values from normative data: water-scaled data from normative data, corrected for from partial-volume corrected water-scaled data
ssd.metfwE = ald.metfwE(:,ssd.aldi); %expected values of concentration
ssd.metfw_ps = ald.metfw_ps(ssd.aldi,:); %fitting parameters to describe expected values as functions of PMA
ssd.mu = ald.mu; %normalization parameters used tofit expected values as functions of PMA
ssd.all_PMA = ald.PMA; %post-menstrual age of normative data.

% %%%%%%%%%%%%%%%%%%%%%%%%
% %option:  use data from one of the infants in the normatived data:
% infno =1;
% %get signals from the normative data
% ssd.sgn = ssd.metw(infno,:);
% 
% %get PMA from normative data:
% ssd.PMA = ssd.all_PMA(infno);
% %%%%%%%%%%%%%%%%%%%%%%%

%numbers of metabolites, subjects
ssd.Nmets = length(ssd.met_names);
ssd.Nsubs = size(ssd.metw,1);

%Identify which of the available metabolites are reference metabolites:
ssd.my_refs = ismember(ssd.met_names,ssd.ref_metabolites);
ssd.refi = find(ssd.my_refs);

%If there is one reference metabolite, no results for that metaoblite will
%be reported
%variables with "rpt" in the name refer to metabolites that will be in the
%report.
if length(ssd.ref_metabolites) == 1
    ssd.Nmets_rpt = ssd.Nmets-1;       
    ssd.met_names_rpt = ssd.met_names(~strcmp(ssd.met_names,ssd.ref_metabolites));
else
    ssd.Nmets_rpt = ssd.Nmets;
    ssd.met_names_rpt = ssd.met_names;    
end
ssd.rpt_inds = find(ismember(ssd.met_names,ssd.met_names_rpt));

%ssd.mat_refs is a matrix of reference metabolites; each row give the reference
%metabolites for that metabolite.
ssd.mat_refs = NaN*ones(ssd.Nmets_rpt,ssd.Nmets);
for jj = 1:ssd.Nmets_rpt
    %calculat a metabolite-specific mask
    ssd.mat_refs(jj,:) = ssd.my_refs;
    if strcmp(ssd.is_self,'n') %no self referencing
        ssd.mat_refs(jj,jj) = 0; %no self-referencing
        if strcmp(ssd.met_names{jj},'Glu') %we don't want to use Glx to reference Glu
            ssd.mat_refs(jj,strcmp(ssd.met_names,'Glx')) = 0;
        end
        if strcmp(ssd.met_names{jj},'NAA') %we don't want to use tNAA to reference NAA
            ssd.mat_refs(jj,strcmp(ssd.met_names,'tNAA')) = 0;
        end
        if strcmp(ssd.met_names{jj},'Glx') %we don't want to use Glx to reference Glu
            ssd.mat_refs(jj,strcmp(ssd.met_names,'Glu')) = 0;
        end
        if strcmp(ssd.met_names{jj},'tNAA') %we don't want to use tNAA to reference NAA
            ssd.mat_refs(jj,strcmp(ssd.met_names,'NAA')) = 0;
        end
    elseif strcmp(ssd.is_self,'y') %allow self referencing, nothing to be done
    else
        error('allow self referencing? y or n')
    end
end

% %calculate weighting factors:
ssd.metw_wn = NaN*ones(ssd.Nmets_rpt,ssd.Nmets);
allmets = ones(1,ssd.Nmets);
for jj = 1:ssd.Nmets_rpt
    msk = logical(ssd.mat_refs(jj,:));
    har_mean = 1/mean(1./ssd.metwSD(msk));
    if strcmp(ssd.is_weighting,'y') %if weighiting factors are to be used:
        ssd.metw_wn(jj,msk) = har_mean./ssd.metwSD(msk)/nnz(msk);
    elseif strcmp(ssd.is_weighting,'n')
        ssd.metw_wn(jj,msk) = allmets(msk);
    else
        error('use weighting factors? y or n')
    end
end

%For the single subject, calculate expected concentration values using 
%polynomial constants obtained from the normative data:
for jj = 1:ssd.Nmets
    ssd.C_E(jj) = polyval(ssd.metfw_ps(jj,:),ssd.PMA,[],ssd.mu)';
end

%for the single subject, calculate R_CRef and C_CRef:
for jj = 1:ssd.Nmets_rpt      
    msk = logical(ssd.mat_refs(jj,:));  %identify reference metabolties
    kk = ssd.rpt_inds(jj);  %identify the index number of each metabolite in the signal vector
    %R_CReF:
    ssd.R_CRef(jj) = ssd.sgn(kk)/sum(ssd.metw_wn(jj,msk).*ssd.sgn(msk));
    %C_CRef:
    ssd.C_CRef(jj) = ssd.R_CRef(jj).*sum(ssd.metw_wn(jj,msk).*ssd.C_E(msk));
end

%For the normative data, calculate R_CRef and C_CRef
for kk = 1:ssd.Nsubs
    for jj = 1:ssd.Nmets_rpt
        msk = logical(ssd.mat_refs(jj,:));        
        mm = ssd.rpt_inds(jj);  %identify the index number of each metabolite in the metw.
        %R_CReF:
        ssd.all_R_CRef(kk,jj) = ssd.metw(kk,mm)/sum(ssd.metw_wn(jj,msk).*ssd.metw(kk,msk));
        %C_CRef:
        ssd.all_C_CRef(kk,jj) = ssd.all_R_CRef(kk,jj).*sum(ssd.metw_wn(jj,msk).*ssd.metfwE(kk,msk));
    end
end

%for normative data, calcuate average PMA dependence:
for jj = 1:ssd.Nmets_rpt
    
    kk = ssd.rpt_inds(jj);  %identify the index number of each metabolite in the expectation value parameters.
    [r,~,~] = find(~isnan(ssd.all_R_CRef(:,jj))); %identify NaNs
    %R_CRef
    [ps] = polyfit((ssd.all_PMA(r)-ssd.mu(1))/ssd.mu(2),ssd.all_R_CRef(r,jj),2);
    ssd.all_R_CRef_ps(jj,:) = ps;
    ssd.all_R_CRef_E(:,jj) = polyval(ps,ssd.all_PMA,[],ssd.mu); %expected values for normatived data.
    ssd.R_CRef_E(jj) = polyval(ps,ssd.PMA,[],ssd.mu); %expected value for the single subject.
    
    %C_CRef
    [ps] = polyfit((ssd.all_PMA(r)-ssd.mu(1))/ssd.mu(2),ssd.all_C_CRef(r,jj),2);
    ssd.all_C_CRef_ps(jj,:) = ps;
    ssd.all_C_CRef_E(:,jj) = polyval(ps,ald.PMA,[],ssd.mu); %expected values for normatived data.
    ssd.C_CRef_E(jj) = polyval(ps,ssd.PMA,[],ssd.mu); %expected value for the single subject.
end

%for normative data, calculate SD values of C_CRef and R_CRef
for jj = 1:ssd.Nmets_rpt
    [r,~,~] = find(~isnan(ssd.all_R_CRef(:,jj))); %identify NaNs
    ssd.R_CRef_SD(jj) = std(ssd.all_R_CRef(r,jj) - ssd.all_R_CRef_E(r,jj));
    ssd.C_CRef_SD(jj) = std(ssd.all_C_CRef(r,jj) - ssd.all_C_CRef_E(r,jj));
end

%for normative data, calculate error estimates via comparison of estimated concentrations of normative data with "ground truth" :corrected-water scaled values from a subset of infants
for jj = 1:ssd.Nmets_rpt
  ssd.C_Cref_error(jj) =nanstd(ssd.all_C_CRef(:,jj) - ssd.metfw(:,ssd.rpt_inds(jj) ));
end

%calculate percentiles of the single subject, based on group statistics.
for jj = 1:ssd.Nmets_rpt
    pct_R_CRef(jj) = round(100*normcdf(ssd.R_CRef(jj),ssd.R_CRef_E(jj),ssd.R_CRef_SD(jj)),0);
    pct_C_CRef(jj) = round(100*normcdf(ssd.C_CRef(jj),ssd.C_CRef_E(jj),ssd.C_CRef_SD(jj)),0);
end

%create a table and spreadsheet showing the results
restab = table;
restab.metabolites = ssd.met_names_rpt';
restab.R_CRef = ssd.R_CRef';
restab.percentile_R_CRef = pct_R_CRef';
restab.C_CRef = ssd.C_CRef';
restab.percentile_C_CRef = pct_C_CRef';

% save results spreedsheet and the ssd structure to the results directory.
results_spreadsheet = 'results.csv';
delete(fullfile(results_dir,results_spreadsheet));
writetable(restab, fullfile(results_dir,results_spreadsheet));

save(fullfile(results_dir,'ssd'),'ssd')

%plot data
xs = 260:1:440;% PMA range for continuous plots
for jj = 1:ssd.Nmets_rpt
    
    %construct a title line describing which references to use.
    myrefs = ssd.met_names(logical(ssd.mat_refs(jj,:)))
    line2 = ['refs: '];
    for pp = 1:length(myrefs)
        if pp == 1
            line2 = [line2 myrefs{pp}]
        else
            line2 = [line2 ', ' myrefs{pp}]
        end
    end
    
    kk = ssd.rpt_inds(jj);  %identify the index number of each metabolite in the expectation value parameters.
    
    figure(jj) %R_CRef
    %plot the single subject:
    plot(ssd.PMA*12/365.25,ssd.R_CRef(jj),'ro','markersize',9,'linewidth',2,'markerfacecolor','r')
    set(gca,'YMinorTick','on')
    set(gca,'XMinorTick','on')
    hold on;
    
    %plot the curve fits
    ys= polyval(ssd.all_R_CRef_ps(jj,:),xs,[],ssd.mu);
    plot(xs*12/365.25,ys,'k-','linewidth',2);
    plot(xs*12/365.25,ys+ssd.R_CRef_SD(jj),'k--','linewidth',2);
    plot(xs*12/365.25,ys-ssd.R_CRef_SD(jj),'k--','linewidth',2);
    
    %option: plot individuals from the normative data
    plot(ald.PMA*12/365.25,ssd.all_R_CRef(:,jj),'ks','markersize',3)
    
    set(gca,'TickLength',[0.03 0.03])
    set(gca,'Fontsize',11,'linewidth',1.2)
    
    title({'R_{CRef}';strcat(ssd.met_names_rpt{jj},' percentile: ', num2str(pct_R_CRef(jj)),'%') ;line2},'fontsize',16)
    ylabel(ssd.met_names_rpt{jj},'fontsize',18)
    xlabel('PMA (months)','fontsize',18)
    xlim([260 440]*12/365.25)
    
    set(gca,'color','w');
    set(gca,'Fontsize',16,'linewidth',2)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    figure(100+jj) %C_CRef
    
    %plot the single subject:
    errorbar(ssd.PMA*12/365.25,ssd.C_CRef(jj),ssd.C_Cref_error(jj),'ro','markerfacecolor','r','markersize',9,'linewidth',2)
    set(gca,'YMinorTick','on')
    set(gca,'XMinorTick','on')
    hold on;
    
    %plot the curve fits
    ys= polyval(ssd.all_C_CRef_ps(jj,:),xs,[],ssd.mu);
    plot(xs*12/365.25,ys,'k-','linewidth',2);
    plot(xs*12/365.25,ys+ssd.C_CRef_SD(jj),'k--','linewidth',2);
    plot(xs*12/365.25,ys-ssd.C_CRef_SD(jj),'k--','linewidth',2);
    
    %option: plot individuals from the normative data
    plot(ald.PMA*12/365.25,ssd.all_C_CRef(:,jj),'ks','markersize',3)
    
%     %option: plot the expected concentration for this metabolite.
%     ys = polyval(ssd.metfw_ps(jj,:),xs,[],ssd.mu)';
%     plot(xs*12/365.25,ys,'m--','linewidth',1.5);
    
    set(gca,'TickLength',[0.03 0.03])
    set(gca,'Fontsize',11,'linewidth',1.2)
    
    title({'C_{CRef}';strcat(ssd.met_names_rpt{jj},' percentile: ', num2str(pct_C_CRef(jj)),'%') ;line2},'fontsize',16)
    ylabel(strcat(ssd.met_names_rpt{jj},' (mmol/kg)'),'fontsize',18)
    xlabel('PMA (months)','fontsize',18)
    xlim([260 440]*12/365.25)
    
    set(gca,'color','w');
    set(gca,'Fontsize',16,'linewidth',2)
    
    print(strcat('-f',num2str(jj)), '-dtiff',fullfile(results_dir,strcat('met',num2str(jj),'_',ssd.met_names_rpt{jj},'_R_CRef','.tif')));
    print(strcat('-f',num2str(jj+100)), '-dtiff',fullfile(results_dir,strcat('met',num2str(jj),'_',ssd.met_names_rpt{jj},'_C_CRef','.tif')));
end