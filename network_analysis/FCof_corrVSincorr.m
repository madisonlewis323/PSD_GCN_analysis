% Load outcomes
ICD_outcome=readtable('/booboo_workspace/mlewis/UKB_HCP_Neocortex/Code/GCN_SCZ_Classification/ICD_outcomes_test5.csv');
true_ICDvcont=table2array(ICD_outcome(:,'true'));
pred=table2array(ICD_outcome(:,'prediction'));

%Load phenotypes
load('/booboo_workspace/mlewis/UKB_HCP_Neocortex/UKB_HCP_919.mat')

numICD=sum(true_ICDvcont);
numcont=size(pred,1)-numICD;

I_c=1;
I_i=1;
C_c=1;
C_i=1;

%% Retrieve networks and diagnosis
for i=1:length(pred)
    subject=char(table2array(ICD_outcome(i,'subjects')));
    subID=str2num(strtok(subject, 'sub-'));
    index919(i,1)=find(Subject_EIDs(:,1)== subID);
    phenos(i,1)=pred(i,1);
    phenos(i,2)=phenotypes_psychosis(index919(i,1),1);
    phenos(i,3)=phenotypes_control(index919(i,1),1);
    phenos(i,4)=phenotypes_age(index919(i,1),1);
    phenos(i,5)=phenotypes_sex(index919(i,1),1);
    phenos(i,6)=phenotypes_bipolar(index919(i,1),1);
    phenos(i,7)=phenotypes_manic(index919(i,1),1);
    phenos(i,8)=phenotypes_recurrantDepression(index919(i,1),1);
    phenos(i,9)=phenotypes_schizoaffective(index919(i,1),1);
    phenos(i,10)=phenotypes_schizophrenia(index919(i,1),1);
    if true_ICDvcont(i,1)==1 %ICD
        if pred(i,1)==1 %correct
            ICD_correct{I_c,1}=subject;
            MSTD_ICDcorr(:,:,I_c)=readmatrix(['/booboo_workspace/mlewis/UKB_HCP_Neocortex/FCrest/MaxTreeDens/' subject '_FCrest_MaxTreeDenser.csv']);
            I_c=I_c+1;
        else %incorrect
            ICD_incorrect{I_i,1}=subject;
            MSTD_ICDincorr(:,:,I_i)=readmatrix(['/booboo_workspace/mlewis/UKB_HCP_Neocortex/FCrest/MaxTreeDens/' subject '_FCrest_MaxTreeDenser.csv']);
            I_i=I_i+1;
        end
    else %Control
        if pred(i,1)==0 %correct
            CONT_correct{C_c,1}=subject;
            MSTD_CONTcorr(:,:,C_c)=readmatrix(['/booboo_workspace/mlewis/UKB_HCP_Neocortex/FCrest/MaxTreeDens/' subject '_FCrest_MaxTreeDenser.csv']);
            C_c=C_c+1;
        else %incorrect
            CONT_incorrect{C_i,1}=subject;
            MSTD_CONTincorr(:,:,C_i)=readmatrix(['/booboo_workspace/mlewis/UKB_HCP_Neocortex/FCrest/MaxTreeDens/' subject '_FCrest_MaxTreeDenser.csv']);
            C_i=C_i+1;
        end
    end   
end

%% Look at nets
%Controls Correct
for x=1:size(CONT_correct,1)
    BC_CONTcorr(:,x)=betweenness_bin(MSTD_CONTcorr(:,:,x));
    mean_BC_CONTcorr(x)=mean(BC_CONTcorr(:,x));
    std_BC_CONTcorr(x)=std(abs(BC_CONTcorr(:,x)));      
    
    % Degree
    [deg_CONTcorr(:,x)] = degrees_und(MSTD_CONTcorr(:,:,x));
    mean_deg_CONTcorr(x)=mean(deg_CONTcorr(:,x));
    std_deg_CONTcorr(x)=std(deg_CONTcorr(:,x));
     
     % Eigenvector Centrality
     EC_CONTcorr(:,x)= eigenvector_centrality_und(MSTD_CONTcorr(:,:,x));
     mean_EC_CONTcorr(x)=mean(EC_CONTcorr(:,x));
     std_EC_CONTcorr(x)=std(EC_CONTcorr(:,x));
     
     CC_CONTcorr(:,x)=clustering_coef_bu(MSTD_CONTcorr(:,:,x));
     mean_CC_CONTcorr(x)=mean(CC_CONTcorr(:,x));
     std_CC_CONTcorr(x)=std(CC_CONTcorr(:,x));
     
     %Hubs
     Hublist_EC_CONTcorr(:,x)=zeros(length(MSTD_CONTcorr(:,:,x)),1);
     Hublist_BC_CONTcorr(:,x)=zeros(length(MSTD_CONTcorr(:,:,x)),1);
     Hublist_deg_CONTcorr(:,x)=zeros(length(MSTD_CONTcorr(:,:,x)),1);
     for n=1:length(MSTD_CONTcorr(:,:,x))
        if EC_CONTcorr(n,x) > mean_EC_CONTcorr(x) + (2*std_EC_CONTcorr(x))
            Hublist_EC_CONTcorr(n,x)=n;
        end
        if BC_CONTcorr(n,x) > mean_BC_CONTcorr(x) + (2*std_BC_CONTcorr(x))
            Hublist_BC_CONTcorr(n,x)=n;
        end
        if deg_CONTcorr(n,x) > mean_deg_CONTcorr(x) + (2*std_deg_CONTcorr(x))
            Hublist_deg_CONTcorr(n,x)=n;
        end
     end
     
     %Char path and ecc
     [~,dist_CONTcorr]=reachdist(MSTD_CONTcorr(:,:,x));
     [path_CONTcorr(x),~,ECC_CONTcorr(:,x),~,~]=charpath(dist_CONTcorr,0,0);
     
     %Assortativity
     assort_CONTcorr(x) = assortativity_bin(MSTD_CONTcorr(:,:,x),0);
     
     %Modularity
     [Ci_CONTcorr(:,x),Q_CONTcorr(x)] = modularity_und(MSTD_CONTcorr(:,:,x));
end

%Controls Incorrect
for x=1:size(CONT_incorrect,1)
    BC_CONTincorr(:,x)=betweenness_bin(MSTD_CONTincorr(:,:,x));
    mean_BC_CONTincorr(x)=mean(BC_CONTincorr(:,x));
    std_BC_CONTincorr(x)=std(abs(BC_CONTincorr(:,x)));      
    
    % Degree
    [deg_CONTincorr(:,x)] = degrees_und(MSTD_CONTincorr(:,:,x));
    mean_deg_CONTincorr(x)=mean(deg_CONTincorr(:,x));
    std_deg_CONTincorr(x)=std(deg_CONTincorr(:,x));
     
     % Eigenvector Centrality
     EC_CONTincorr(:,x)= eigenvector_centrality_und(MSTD_CONTincorr(:,:,x));
     mean_EC_CONTincorr(x)=mean(EC_CONTincorr(:,x));
     std_EC_CONTincorr(x)=std(EC_CONTincorr(:,x));
     
     CC_CONTincorr(:,x)=clustering_coef_bu(MSTD_CONTincorr(:,:,x));
     mean_CC_CONTincorr(x)=mean(CC_CONTincorr(:,x));
     std_CC_CONTincorr(x)=std(CC_CONTincorr(:,x));
     
     %Hubs
     Hublist_EC_CONTincorr(:,x)=zeros(length(MSTD_CONTincorr(:,:,x)),1);
     Hublist_BC_CONTincorr(:,x)=zeros(length(MSTD_CONTincorr(:,:,x)),1);
     Hublist_deg_CONTincorr(:,x)=zeros(length(MSTD_CONTincorr(:,:,x)),1);
     for n=1:length(MSTD_CONTincorr(:,:,x))
        if EC_CONTincorr(n,x) > mean_EC_CONTincorr(x) + (2*std_EC_CONTincorr(x))
            Hublist_EC_CONTincorr(n,x)=n;
        end
        if BC_CONTincorr(n,x) > mean_BC_CONTincorr(x) + (2*std_BC_CONTincorr(x))
            Hublist_BC_CONTincorr(n,x)=n;
        end
        if deg_CONTincorr(n,x) > mean_deg_CONTincorr(x) + (2*std_deg_CONTincorr(x))
            Hublist_deg_CONTincorr(n,x)=n;
        end
     end
     
     %Char path and ecc
     [~,dist_CONTincorr]=reachdist(MSTD_CONTincorr(:,:,x));
     [path_CONTincorr(x),~,ECC_CONTincorr(:,x),~,~]=charpath(dist_CONTincorr,0,0);
     
     %Assortativity
     assort_CONTincorr(x) = assortativity_bin(MSTD_CONTincorr(:,:,x),0);
     
     %Modularity
     [Ci_CONTincorr(:,x),Q_CONTincorr(x)] = modularity_und(MSTD_CONTincorr(:,:,x));
end

%ICD Correct
for x=1:size(ICD_correct,1)
    BC_ICDcorr(:,x)=betweenness_bin(MSTD_ICDcorr(:,:,x));
    mean_BC_ICDcorr(x)=mean(BC_ICDcorr(:,x));
    std_BC_ICDcorr(x)=std(abs(BC_ICDcorr(:,x)));      
    
    % Degree
    [deg_ICDcorr(:,x)] = degrees_und(MSTD_ICDcorr(:,:,x));
    mean_deg_ICDcorr(x)=mean(deg_ICDcorr(:,x));
    std_deg_ICDcorr(x)=std(deg_ICDcorr(:,x));
     
     % Eigenvector Centrality
     EC_ICDcorr(:,x)= eigenvector_centrality_und(MSTD_ICDcorr(:,:,x));
     mean_EC_ICDcorr(x)=mean(EC_ICDcorr(:,x));
     std_EC_ICDcorr(x)=std(EC_ICDcorr(:,x));
     
     CC_ICDcorr(:,x)=clustering_coef_bu(MSTD_ICDcorr(:,:,x));
     mean_CC_ICDcorr(x)=mean(CC_ICDcorr(:,x));
     std_CC_ICDcorr(x)=std(CC_ICDcorr(:,x));
     
     %Hubs
     Hublist_EC_ICDcorr(:,x)=zeros(length(MSTD_ICDcorr(:,:,x)),1);
     Hublist_BC_ICDcorr(:,x)=zeros(length(MSTD_ICDcorr(:,:,x)),1);
     Hublist_deg_ICDcorr(:,x)=zeros(length(MSTD_ICDcorr(:,:,x)),1);
     for n=1:length(MSTD_ICDcorr(:,:,x))
        if EC_ICDcorr(n,x) > mean_EC_ICDcorr(x) + (2*std_EC_ICDcorr(x))
            Hublist_EC_ICDcorr(n,x)=n;
        end
        if BC_ICDcorr(n,x) > mean_BC_ICDcorr(x) + (2*std_BC_ICDcorr(x))
            Hublist_BC_ICDcorr(n,x)=n;
        end
        if deg_ICDcorr(n,x) > mean_deg_ICDcorr(x) + (2*std_deg_ICDcorr(x))
            Hublist_deg_ICDcorr(n,x)=n;
        end
     end
     
     %Char path and ecc
     [~,dist_ICDcorr]=reachdist(MSTD_ICDcorr(:,:,x));
     [path_ICDcorr(x),~,ECC_ICDcorr(:,x),~,~]=charpath(dist_ICDcorr,0,0);
     
     %Assortativity
     assort_ICDcorr(x) = assortativity_bin(MSTD_ICDcorr(:,:,x),0);
     
     %Modularity
     [Ci_ICDcorr(:,x),Q_ICDcorr(x)] = modularity_und(MSTD_ICDcorr(:,:,x));
end

%ICD Incorrect
for x=1:size(ICD_incorrect,1)
    BC_ICDincorr(:,x)=betweenness_bin(MSTD_ICDincorr(:,:,x));
    mean_BC_ICDincorr(x)=mean(BC_ICDincorr(:,x));
    std_BC_ICDincorr(x)=std(abs(BC_ICDincorr(:,x)));      
    
    % Degree
    [deg_ICDincorr(:,x)] = degrees_und(MSTD_ICDincorr(:,:,x));
    mean_deg_ICDincorr(x)=mean(deg_ICDincorr(:,x));
    std_deg_ICDincorr(x)=std(deg_ICDincorr(:,x));
     
     % Eigenvector Centrality
     EC_ICDincorr(:,x)= eigenvector_centrality_und(MSTD_ICDincorr(:,:,x));
     mean_EC_ICDincorr(x)=mean(EC_ICDincorr(:,x));
     std_EC_ICDincorr(x)=std(EC_ICDincorr(:,x));
     
     CC_ICDincorr(:,x)=clustering_coef_bu(MSTD_ICDincorr(:,:,x));
     mean_CC_ICDincorr(x)=mean(CC_ICDincorr(:,x));
     std_CC_ICDincorr(x)=std(CC_ICDincorr(:,x));
     
     %Hubs
     Hublist_EC_ICDincorr(:,x)=zeros(length(MSTD_ICDincorr(:,:,x)),1);
     Hublist_BC_ICDincorr(:,x)=zeros(length(MSTD_ICDincorr(:,:,x)),1);
     Hublist_deg_ICDincorr(:,x)=zeros(length(MSTD_ICDincorr(:,:,x)),1);
     for n=1:length(MSTD_ICDincorr(:,:,x))
        if EC_ICDincorr(n,x) > mean_EC_ICDincorr(x) + (2*std_EC_ICDincorr(x))
            Hublist_EC_ICDincorr(n,x)=n;
        end
        if BC_ICDincorr(n,x) > mean_BC_ICDincorr(x) + (2*std_BC_ICDincorr(x))
            Hublist_BC_ICDincorr(n,x)=n;
        end
        if deg_ICDincorr(n,x) > mean_deg_ICDincorr(x) + (2*std_deg_ICDincorr(x))
            Hublist_deg_ICDincorr(n,x)=n;
        end
     end
     
     %Char path and ecc
     [~,dist_ICDincorr]=reachdist(MSTD_ICDincorr(:,:,x));
     [path_ICDincorr(x),~,ECC_ICDincorr(:,x),~,~]=charpath(dist_ICDincorr,0,0);
     
     %Assortativity
     assort_ICDincorr(x) = assortativity_bin(MSTD_ICDincorr(:,:,x),0);
     
     %Modularity
     [Ci_ICDincorr(:,x),Q_ICDincorr(x)] = modularity_und(MSTD_ICDincorr(:,:,x));
end

%% organize hubs
% Percolation
%Active
[hubs_CONTcorr, hubs_CONTincorr] = hubs_organize(Hublist_BC_CONTcorr, Hublist_BC_CONTincorr, Hublist_EC_CONTcorr, Hublist_EC_CONTincorr, ... 
    Hublist_deg_CONTcorr, Hublist_deg_CONTincorr);
[hubs_ICDcorr, hubs_ICDincorr, Hublist_ICDcorr, Hublist_ICDincorr] = hubs_organize(Hublist_BC_ICDcorr, Hublist_BC_ICDincorr, Hublist_EC_ICDcorr, Hublist_EC_ICDincorr, ... 
    Hublist_deg_ICDcorr, Hublist_deg_ICDincorr);

%% Hub connections

% Hub L_A4 is node 174/358 so 193/377
for i=1:size(MSTD_ICDcorr,3)
    hub_L_A4_cons_ICDcorr(i,:)=MSTD_ICDcorr(193,:,i);
end
for i=1:size(MSTD_ICDincorr,3)
    hub_L_A4_cons_ICDincorr(i,:)=MSTD_ICDincorr(193,:,i);
end

% Hub R_A4 is node 353/358 so 372/377
for i=1:size(MSTD_ICDcorr,3)
    hub_R_A4_cons_ICDcorr(i,:)=MSTD_ICDcorr(372,:,i);
end
for i=1:size(MSTD_ICDincorr,3)
    hub_R_A4_cons_ICDincorr(i,:)=MSTD_ICDincorr(372,:,i);
end
%% Compare degree between sections

%Load Section labels
sectionlabel=readtable('/booboo_workspace/mlewis/HCPcortical_sectionlabels.xlsx', 'Sheet', 'Cortical Only');
sectionnums=table2array(sectionlabel(:,'Section'));

section_nodenames=cell(32,max(sectionnums));
nodedeg_ICDcorr=zeros(size(MSTD_ICDcorr,3),32,max(sectionnums));
nodedeg_ICDincorr=zeros(size(MSTD_ICDincorr,3),32,max(sectionnums));
nodedeg_CONTcorr=zeros(size(MSTD_CONTcorr,3),32,max(sectionnums));
nodedeg_CONTincorr=zeros(size(MSTD_CONTincorr,3),32,max(sectionnums));

secdeg_ICDcorr=zeros(size(MSTD_ICDcorr,3),22);
secdeg_ICDincorr=zeros(size(MSTD_ICDincorr,3),22);
secdeg_CONTcorr=zeros(size(MSTD_CONTcorr,3),22);
secdeg_CONTincorr=zeros(size(MSTD_CONTincorr,3),22);

for s=1:max(sectionnums)
    indx=find(sectionnums(:,1)~=s);
    %ICD correct
    temp_net=MSTD_ICDcorr;
    temp_net(1:19,:,:)=[];
    temp_net(:,1:19,:)=[];
    temp_net(indx,:,:)=[];
    temp_net(:,indx,:)=[];
    for sub=1:size(temp_net,3)
        nodedeg_ICDcorr(sub,1:size(temp_net,1),s)=degrees_und(temp_net(:,:,sub));
        secdeg_ICDcorr(sub,s)=sum(sum(tril(temp_net(:,:,sub),-1)));
    end
    %ICD incorrect
    temp_net=MSTD_ICDincorr;
    temp_net(1:19,:,:)=[];
    temp_net(:,1:19,:)=[];
    temp_net(indx,:,:)=[];
    temp_net(:,indx,:)=[];
    for sub=1:size(temp_net,3)
        nodedeg_ICDincorr(sub,1:size(temp_net,1),s)=degrees_und(temp_net(:,:,sub));
        secdeg_ICDincorr(sub,s)=sum(sum(tril(temp_net(:,:,sub),-1)));
    end
    %CONT correct
    temp_net=MSTD_CONTcorr;
    temp_net(1:19,:,:)=[];
    temp_net(:,1:19,:)=[];
    temp_net(indx,:,:)=[];
    temp_net(:,indx,:)=[];
    for sub=1:size(temp_net,3)
        nodedeg_CONTcorr(sub,1:size(temp_net,1),s)=degrees_und(temp_net(:,:,sub));
        secdeg_CONTcorr(sub,s)=sum(sum(tril(temp_net(:,:,sub),-1)));
    end
    %CONT incorrect
    temp_net=MSTD_CONTincorr;
    temp_net(1:19,:,:)=[];
    temp_net(:,1:19,:)=[];
    temp_net(indx,:,:)=[];
    temp_net(:,indx,:)=[];
    for sub=1:size(temp_net,3)
        nodedeg_CONTincorr(sub,1:size(temp_net,1),s)=degrees_und(temp_net(:,:,sub));
        secdeg_CONTincorr(sub,s)=sum(sum(tril(temp_net(:,:,sub),-1)));
    end
    %Get node names
    x=1;
    for n=1:358
        if sectionnums(n) == s
            section_nodenames{x,s}=NodeNames{n+19,1};
            x=x+1;
        end
    end
end
%% Compare sections

for s=1:22
    [~, pvalue_ContCorrVContIncorr(s), ~, stats_ContCorrVContIncorr(s)] = ttest2(secdeg_CONTcorr(1:end,s), secdeg_CONTincorr(1:end,s));
    [~, pvalue_ContCorrVICDIncorr(s), ~, stats_ContCorrVICDIncorr(s)] = ttest2(secdeg_CONTcorr(1:end,s), secdeg_ICDincorr(1:end,s));
    [~, pvalue_ContCorrVICDCorr(s), ~, stats_ContCorrVICDCorr(s)] = ttest2(secdeg_CONTcorr(1:end,s), secdeg_ICDcorr(1:end,s));
    [~, pvalue_ContIncorrVICDcorr(s), ~, stats_ContIncorrVICDcorr(s)] = ttest2(secdeg_CONTincorr(1:end,s), secdeg_ICDcorr(1:end,s));
    [~, pvalue_ICDcorrVICDincorr(s), ~, stats_ICDcorrVICDIncorr(s)] = ttest2(secdeg_ICDcorr(1:end,s), secdeg_ICDincorr(1:end,s));
    %Classified as contol vs classified as psychosis
    [~, pvalue_sec_CONTVICD(s), ~, stats_sec_CONTVICD(s)] = ttest2([secdeg_CONTcorr(1:end,s); secdeg_ICDincorr(1:end,s)], [secdeg_ICDcorr(1:end,s); secdeg_CONTincorr(1:end,s)]);    
end
%% Graph measures 
% Average BC
[~, pvalue_BC_ICDcvICDi, ~, stats_BC_ICDcvICDi] = ttest2(mean_BC_ICDcorr, mean_BC_ICDincorr);
[~, pvalue_BC_ICDcvCONTi, ~, stats_BC_ICDcvCONTi] = ttest2(mean_BC_ICDcorr, mean_BC_CONTincorr);
[~, pvalue_BC_ICDcvCONTc, ~, stats_BC_ICDcvCONTc] = ttest2(mean_BC_ICDcorr, mean_BC_CONTcorr);
[~, pvalue_BC_CONTcvCONTi, ~, stats_BC_CONTcvCONTi] = ttest2(mean_BC_CONTcorr, mean_BC_CONTincorr);
[~, pvalue_BC_CONTcvICDi, ~, stats_BC_CONTcvICDi] = ttest2(mean_BC_CONTcorr, mean_BC_ICDincorr);
% Average EC
[~, pvalue_EC_ICDcvICDi, ~, stats_EC_ICDcvICDi] = ttest2(mean_EC_ICDcorr, mean_EC_ICDincorr);
[~, pvalue_EC_ICDcvCONTi, ~, stats_EC_ICDcvCONTi] = ttest2(mean_EC_ICDcorr, mean_EC_CONTincorr);
[~, pvalue_EC_ICDcvCONTc, ~, stats_EC_ICDcvCONTc] = ttest2(mean_EC_ICDcorr, mean_EC_CONTcorr);
[~, pvalue_EC_CONTcvCONTi, ~, stats_EC_CONTcvCONTi] = ttest2(mean_EC_CONTcorr, mean_EC_CONTincorr);
[~, pvalue_EC_CONTcvICDi, ~, stats_EC_CONTcvICDi] = ttest2(mean_EC_CONTcorr, mean_EC_ICDincorr);
% Modularity
[~, pvalue_Q_ICDcvICDi, ~, stats_Q_ICDcvICDi] = ttest2(Q_ICDcorr, Q_ICDincorr);
[~, pvalue_Q_ICDcvCONTi, ~, stats_Q_ICDcvCONTi] = ttest2(Q_ICDcorr, Q_CONTincorr);
[~, pvalue_Q_ICDcvCONTc, ~, stats_Q_ICDcvCONTc] = ttest2(Q_ICDcorr, Q_CONTcorr);
[~, pvalue_Q_CONTcvCONTi, ~, stats_Q_CONTcvCONTi] = ttest2(Q_CONTcorr, Q_CONTincorr);
[~, pvalue_Q_CONTcvICDi, ~, stats_Q_CONTcvICDi] = ttest2(Q_CONTcorr, Q_ICDincorr);

