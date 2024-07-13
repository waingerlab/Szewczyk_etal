%% This function reads in two sets of DEGs from an excel sheet and compares them. The format in the excel sheet should be:
%column1: Gene names
%column2: DEG set 1 log2FC
%column3: DEG set 1 padj
%column4: DEG set 2 log2FC
%column5: DEG set 2 padj
%
%Then just specify the excel sheet file name, sheet name, and the range,
%which should be A2:E26000 ish
%
%Be sure to make all the dependent scripts available before calling the function, which should be something like:
%addpath('C:\Data\Wainger_Lab\Matlab\Scripts\Universal_functions','C:\Data\Wainger_Lab\Matlab\Scripts\Graphs'; 'C:\Data\Wainger_Lab\Matlab\Scripts\RNAseq')


function [data] = DEG_comparison(file,sheet)


[~,names,~] = xlsread(file,sheet,'A2:A17605');
data = xlsread(file,sheet,'B2:E17605'); %data should be organized as group1_log2fc; group1_padj; group2_log2fc; group2_padj

sig1=data(:,2)<.05; %sig genes for group1
sig2=data(:,4)<.05; %sig genes for group2
gp1_up=data(:,1)>0;
gp1_down=data(:,1)<0;
gp2_up=data(:,3)>0;
gp2_down=data(:,3)<0;

%significant genes
gp1_sigup=sig1.*gp1_up;
gp1_sigdown=sig1.*gp1_down;
gp2_sigup=sig2.*gp2_up;
gp2_sigdown=sig2.*gp2_down;

%intersection
both_sigup=gp1_sigup.*gp2_sigup;
both_sigdown=gp1_sigdown.*gp2_sigdown;

%unique to each
unique_gp1_up=gp1_sigup-both_sigup;
unique_gp1_down=gp1_sigdown-both_sigdown;
unique_gp2_up=gp2_sigup-both_sigup;
unique_gp2_down=gp2_sigdown-both_sigdown;

%lists and sig values
bothup_names=names(logical(both_sigup));
bothup_values=data(logical(both_sigup),:);
bothdown_names=names(logical(both_sigdown));
bothdown_values=data(logical(both_sigdown),:);

gp1_up_names=names(logical(unique_gp1_up));
gp1_up_values=data(logical(unique_gp1_up),:);
gp1_down_names=names(logical(unique_gp1_down));
gp1_down_values=data(logical(unique_gp1_down),:);

gp2_up_names=names(logical(unique_gp2_up));
gp2_up_values=data(logical(unique_gp2_up),:);
gp2_down_names=names(logical(unique_gp2_down));
gp2_down_values=data(logical(unique_gp2_down),:);

%write out data
xlswrite(file,bothup_names,sheet,'g2')
xlswrite(file,bothup_values,sheet,'h2')
xlswrite(file,bothdown_names,sheet,'m2')
xlswrite(file,bothdown_values,sheet,'n2')
xlswrite(file,{'Both_Up', 'Group1_log2FC','Group1_padj','Group2_log2FC','Group2_padj'},sheet,'g1');
xlswrite(file,{'Both_Down', 'Group1_log2FC','Group1_padj','Group2_log2FC','Group2_padj'},sheet,'m1');

xlswrite(file,gp1_up_names,sheet,'s2')
xlswrite(file,gp1_up_values,sheet,'t2')
xlswrite(file,gp1_down_names,sheet,'y2')
xlswrite(file,gp1_down_values,sheet,'z2')
xlswrite(file,{'Group1_Up_Unique', 'Group1_log2FC','Group1_padj','Group2_log2FC','Group2_padj'},sheet,'s1');
xlswrite(file,{'Group1_Down_Unique', 'Group1_log2FC','Group1_padj','Group2_log2FC','Group2_padj'},sheet,'y1');

xlswrite(file,gp2_up_names,sheet,'ae2')
xlswrite(file,gp2_up_values,sheet,'af2')
xlswrite(file,gp2_down_names,sheet,'ak2')
xlswrite(file,gp2_down_values,sheet,'al2')
xlswrite(file,{'Group2_Up_Unique', 'Group1_log2FC','Group1_padj','Group2_log2FC','Group2_padj'},sheet,'ae1');
xlswrite(file,{'Group2_Down_Unique', 'Group1_log2FC','Group1_padj','Group2_log2FC','Group2_padj'},sheet,'ak1');

%venn diagram up
venn([numel(gp1_sigup) numel(gp2_sigup)],numel(bothup_names),'FaceColor',{'r','k'})
hold on
ylim([-40 40])
txt= strcat('both=',int2str(numel(bothup_names)),'; group1=',int2str(numel(gp1_up_names)),'; group2=',int2str(numel(gp2_up_names)))
text(0,40,txt)
hold off
saveas(gcf,'venn_up.pdf')
close

%venn diagram down
venn([numel(gp1_sigdown) numel(gp2_sigdown)],numel(bothdown_names),'FaceColor',{'r','k'})
hold on
ylim([-40 40])
txt= strcat('both=',int2str(numel(bothdown_names)),'; group1=',int2str(numel(gp1_down_names)),'; group2=',int2str(numel(gp2_down_names)))
text(0,40,txt)
hold off
saveas(gcf,'venn_down.pdf')
close

%correlation plot- significant in either
corrgenes=logical(sig1.*sig2);
gp1_corrgenes=data(corrgenes,1);
gp2_corrgenes=data(corrgenes,3);
[r,p]= corrcoef(gp1_corrgenes,gp2_corrgenes)

txt=strcat('R=',num2str(r(1,2)),'; P=',num2str(p(1,2)),'; gene number=',int2str(sum(corrgenes)))
scatter(gp1_corrgenes,gp2_corrgenes,3,'k','filled')
l=lsline
l.Color='r';
text(-2,max(gp2_corrgenes)+.05,txt)
saveas(gcf,'correlation.pdf')
close