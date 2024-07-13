%% This compares DEGs from 3 groups. It outputs a triple Venn diagram and lists of DEGs. The input should be the filename and sheet to read 7 columns: 
% 1) Gene names
% 2) group 1 log2FC
% 3) group 1 padj
% 4) group 2 log2FC
% 5) group 2 padj
% 6) group 3 log2FC
% 7) group 3 padj
%
%Then just specify the excel sheet file name, sheet name, and the range,
%which should be A2:E26000 ish
%
%Be sure to make all the dependent scripts available before calling the function, which should be something like:
%addpath('C:\Data\Wainger_Lab\Matlab\Scripts\Universal_functions','C:\Data\Wainger_Lab\Matlab\Scripts\Graphs'; 'C:\Data\Wainger_Lab\Matlab\Scripts\RNAseq')


function [data] = Triple_DEG_comparison(file,sheet)

[~,names,~] = xlsread(file,sheet,'A2:A17605');
data = xlsread(file,sheet,'B2:G17605'); %data should be organized as group1_log2fc, group1_padj, group2_log2fc, group2_padj, group2_log2fc, group2_padj

sig1=data(:,2)<.05; %sig genes for group1
sig2=data(:,4)<.05; %sig genes for group2
sig3=data(:,6)<.05; %sig genes for group2
gp1_up=data(:,1)>0;
gp1_down=data(:,1)<0;
gp2_up=data(:,3)>0;
gp2_down=data(:,3)<0;
gp3_up=data(:,5)>0;
gp3_down=data(:,5)<0;

%significant genes
gp1_sigup=sig1.*gp1_up;
gp1_sigdown=sig1.*gp1_down;
gp2_sigup=sig2.*gp2_up;
gp2_sigdown=sig2.*gp2_down;
gp3_sigup=sig3.*gp3_up;
gp3_sigdown=sig3.*gp3_down;

%triple intersection
all_up = gp1_sigup.*gp2_sigup.*gp3_sigup;
all_down = gp1_sigdown.*gp2_sigdown.*gp3_sigdown;

%double intersection
gp1_gp2_sigup=[gp1_sigup.*gp2_sigup]-all_up;
gp1_gp2_sigdown=[gp1_sigdown.*gp2_sigdown]-all_down;

gp1_gp3_sigup=[gp1_sigup.*gp3_sigup]-all_up;
gp1_gp3_sigdown=[gp1_sigdown.*gp3_sigdown]-all_down;

gp2_gp3_sigup=[gp2_sigup.*gp3_sigup]-all_up;
gp2_gp3_sigdown=[gp2_sigdown.*gp3_sigdown]-all_down;

%unique to each
unique_gp1_up=gp1_sigup-all_up-gp1_gp3_sigup-gp1_gp2_sigup;
unique_gp1_down=gp1_sigdown-all_down-gp1_gp3_sigdown-gp1_gp2_sigdown;
unique_gp2_up=gp2_sigup-all_up-gp1_gp2_sigup-gp2_gp3_sigup;
unique_gp2_down=gp2_sigdown-all_down-gp1_gp2_sigdown-gp2_gp3_sigdown;
unique_gp3_up=gp3_sigup-all_up-gp1_gp3_sigup-gp2_gp3_sigup;
unique_gp3_down=gp3_sigdown-all_down-gp1_gp3_sigdown-gp2_gp3_sigdown;

%lists and sig values
allup_names=names(logical(all_up));
allup_values=data(logical(all_up),:);
alldown_names=names(logical(all_down));
alldown_values=data(logical(all_down),:);

gp1_gp2_up_names=names(logical(gp1_gp2_sigup));
gp1_gp2_up_values=data(logical(gp1_gp2_sigup),:);
gp1_gp2_down_names=names(logical(gp1_gp2_sigdown));
gp1_gp2_down_values=data(logical(gp1_gp2_sigdown),:);

gp1_gp3_up_names=names(logical(gp1_gp3_sigup));
gp1_gp3_up_values=data(logical(gp1_gp3_sigup),:);
gp1_gp3_down_names=names(logical(gp1_gp3_sigdown));
gp1_gp3_down_values=data(logical(gp1_gp3_sigdown),:);

gp2_gp3_up_names=names(logical(gp2_gp3_sigup));
gp2_gp3_up_values=data(logical(gp2_gp3_sigup),:);
gp2_gp3_down_names=names(logical(gp2_gp3_sigdown));
gp2_gp3_down_values=data(logical(gp2_gp3_sigdown),:);

gp1_up_names=names(logical(unique_gp1_up));
gp1_up_values=data(logical(unique_gp1_up),:);
gp1_down_names=names(logical(unique_gp1_down));
gp1_down_values=data(logical(unique_gp1_down),:);

gp2_up_names=names(logical(unique_gp2_up));
gp2_up_values=data(logical(unique_gp2_up),:);
gp2_down_names=names(logical(unique_gp2_down));
gp2_down_values=data(logical(unique_gp2_down),:);

gp3_up_names=names(logical(unique_gp3_up));
gp3_up_values=data(logical(unique_gp3_up),:);
gp3_down_names=names(logical(unique_gp3_down));
gp3_down_values=data(logical(unique_gp3_down),:);

%write out data
xlswrite(file,allup_names,sheet,'i2')
xlswrite(file,allup_values,sheet,'j2')
xlswrite(file,alldown_names,sheet,'q2')
xlswrite(file,alldown_values,sheet,'r2')
xlswrite(file,{'All_Up', 'Group1_log2FC','Group1_padj','Group2_log2FC','Group2_padj', 'Group3_log2FC', 'Group3_padj'},sheet,'i1');
xlswrite(file,{'All_Down', 'Group1_log2FC','Group1_padj','Group2_log2FC','Group2_padj', 'Group3_log2FC', 'Group3_padj'},sheet,'q1');

xlswrite(file,gp1_gp2_up_names,sheet,'y2')
xlswrite(file,gp1_gp2_up_values,sheet,'z2')
xlswrite(file,gp1_gp2_down_names,sheet,'ag2')
xlswrite(file,gp1_gp2_down_values,sheet,'ah2')
xlswrite(file,{'Group1_Group2_Up', 'Group1_log2FC','Group1_padj','Group2_log2FC','Group2_padj', 'Group3_log2FC', 'Group3_padj'},sheet,'y1');
xlswrite(file,{'Group1_Group2_Down', 'Group1_log2FC','Group1_padj','Group2_log2FC','Group2_padj', 'Group3_log2FC', 'Group3_padj'},sheet,'ag1');

xlswrite(file,gp1_gp3_up_names,sheet,'ao2')
xlswrite(file,gp1_gp3_up_values,sheet,'ap2')
xlswrite(file,gp1_gp3_down_names,sheet,'aw2')
xlswrite(file,gp1_gp3_down_values,sheet,'ax2')
xlswrite(file,{'Group1_Group3_Up', 'Group1_log2FC','Group1_padj','Group2_log2FC','Group2_padj', 'Group3_log2FC', 'Group3_padj'},sheet,'ao1');
xlswrite(file,{'Group1_Group3_Down', 'Group1_log2FC','Group1_padj','Group2_log2FC','Group2_padj', 'Group3_log2FC', 'Group3_padj'},sheet,'aw1');

xlswrite(file,gp2_gp3_up_names,sheet,'be2')
xlswrite(file,gp2_gp3_up_values,sheet,'bf2')
xlswrite(file,gp2_gp3_down_names,sheet,'bm2')
xlswrite(file,gp2_gp3_down_values,sheet,'bn2')
xlswrite(file,{'Group2_Group3_Up', 'Group1_log2FC','Group1_padj','Group2_log2FC','Group2_padj', 'Group3_log2FC', 'Group3_padj'},sheet,'be1');
xlswrite(file,{'Group2_Group3_Down', 'Group1_log2FC','Group1_padj','Group2_log2FC','Group2_padj', 'Group3_log2FC', 'Group3_padj'},sheet,'bm1');

xlswrite(file,gp1_up_names,sheet,'bu2')
xlswrite(file,gp1_up_values,sheet,'bv2')
xlswrite(file,gp1_down_names,sheet,'cc2')
xlswrite(file,gp1_down_values,sheet,'cd2')
xlswrite(file,{'Group1_Up_Unique', 'Group1_log2FC','Group1_padj','Group2_log2FC','Group2_padj', 'Group3_log2FC', 'Group3_padj'},sheet,'bu1');
xlswrite(file,{'Group1_Down_Unique', 'Group1_log2FC','Group1_padj','Group2_log2FC','Group2_padj', 'Group3_log2FC', 'Group3_padj'},sheet,'cc1');

xlswrite(file,gp2_up_names,sheet,'ck2')
xlswrite(file,gp2_up_values,sheet,'cl2')
xlswrite(file,gp2_down_names,sheet,'cs2')
xlswrite(file,gp2_down_values,sheet,'ct2')
xlswrite(file,{'Group2_Up_Unique', 'Group1_log2FC','Group1_padj','Group2_log2FC','Group2_padj', 'Group3_log2FC', 'Group3_padj'},sheet,'ck1');
xlswrite(file,{'Group2_Down_Unique', 'Group1_log2FC','Group1_padj','Group2_log2FC','Group2_padj', 'Group3_log2FC', 'Group3_padj'},sheet,'cs1');

xlswrite(file,gp3_up_names,sheet,'da2')
xlswrite(file,gp3_up_values,sheet,'db2')
xlswrite(file,gp3_down_names,sheet,'di2')
xlswrite(file,gp3_down_values,sheet,'dj2')
xlswrite(file,{'Group3_Up_Unique', 'Group1_log2FC','Group1_padj','Group2_log2FC','Group2_padj', 'Group3_log2FC', 'Group3_padj'},sheet,'da1');
xlswrite(file,{'Group3_Down_Unique', 'Group1_log2FC','Group1_padj','Group2_log2FC','Group2_padj', 'Group3_log2FC', 'Group3_padj'},sheet,'di1');


%venn diagrams made using file exchange respository found here: https://www.mathworks.com/matlabcentral/fileexchange/22282-venn
%venn diagram up
venn([numel(gp1_up_names), numel(gp2_up_names), numel(gp3_up_names), numel(gp1_gp2_up_names),numel(gp1_gp3_up_names),numel(gp2_gp3_up_names),numel(allup_names)])
hold on
ylim([-40 60])
txt= strcat('z1=',int2str(numel(gp1_up_names)),'; z2=',int2str(numel(gp2_up_names)),'; z3=',int2str(numel(gp3_up_names)),'; z12=',int2str(numel(gp1_gp2_up_names)),'; z13=',int2str(numel(gp1_gp3_up_names)),'; z23=',int2str(numel(gp2_gp3_up_names)),'; z123=',int2str(numel(allup_names)));
text(0,60,txt)
hold off
saveas(gcf,'venn_up.pdf')
close

%venn diagram down
venn([numel(gp1_down_names), numel(gp2_down_names), numel(gp3_down_names), numel(gp1_gp2_down_names),numel(gp1_gp3_down_names),numel(gp2_gp3_down_names),numel(alldown_names)])
hold on
ylim([-40 60])
txt= strcat('z1=',int2str(numel(gp1_down_names)),'; z2=',int2str(numel(gp2_down_names)),'; z3=',int2str(numel(gp3_down_names)),'; z12=',int2str(numel(gp1_gp2_down_names)),'; z13=',int2str(numel(gp1_gp3_down_names)),'; z23=',int2str(numel(gp2_gp3_down_names)),'; z123=',int2str(numel(alldown_names)));
text(0,60,txt)
hold off
saveas(gcf,'venn_down.pdf')
close
