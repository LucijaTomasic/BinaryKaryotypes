%% TCGA ASCAT
% loads summary.ascatv3TCGA.penalty70.hg38.tsv and segments downloaded from
% ASCAT github, filters the data, calculates chromosome copy num as mode
% of segmental copy num and plots the data
clear
%%
path = strcat('/Documents/ .. /ASCAT'); % path to folder where summary.ascatv3TCGA.penalty70.hg38.tsv file is located
cd (path)
Files = dir('*.tsv');
clinical_data=readtable(Files.name,"FileType","text");
names= clinical_data.name( (clinical_data.purity<1.000) & (clinical_data.purity>0.4)); %ALL DATA
%names= clinical_data.name( (clinical_data.purity<1.000) & (clinical_data.purity>0.4) & (clinical_data.WGD<1) ); %ALL DATA WGD-

%%
path = strcat('/Documents/ ... /ASCAT/segments 3'); % path to folder where segment txt files are located
cd (path)
Files = dir('*.txt');
numfiles = length(names);

patient_table = cell(1, numfiles);
karyotype=zeros(numfiles, 23);
%% calculate chromosome copy num as statistical  mode of segmental copy numbers
for s = 1:numfiles  
    
    txtname = strcat(names(s),'*.txt');
    File = dir(string(txtname));
    patient_table{s} = readtable(File.name);
    patient_table{s}.("size") = patient_table{s}.endpos - patient_table{s}.startpos; 
    patient_table{s}.("copy_num") = patient_table{s}.nMajor + patient_table{s}.nMinor; 
    
    for chr_i= 1:22
        if ~isa(patient_table{s}.chr,'double')
            patient_table{s}.chr=str2double(patient_table{s}.chr);
        end
        idx_chr_i = (patient_table{s}.chr==chr_i);
        opcije_copy_num=unique(patient_table{s}.copy_num(idx_chr_i));
        ukupni_interval=zeros(length(opcije_copy_num), 1);
        i=1;
        for opcije_i = opcije_copy_num'
           ukupni_interval(i)=sum( patient_table{s}.size( ((patient_table{s}.copy_num == opcije_i) & idx_chr_i) ) ) ;
            i=i+1;
        end
        [~,i_max]=max(ukupni_interval);
        karyotype(s,chr_i) = opcije_copy_num(i_max);
    end
    
end

%% Total chr num 1N 2N 3N 4N 
chr_x= 23-1;

total_chr_copy_num = sum(karyotype,2);
total_chr_copy_num = sort(total_chr_copy_num);


figure('color','white');
hold on;
histogram(total_chr_copy_num,'BinMethod', 'integers');%,'BinLimits',[55,77])
histogram(total_chr_copy_num,'BinMethod', 'integers','BinLimits',[55,77],'FaceColor','r')
xlabel('Ploidy')
ylabel('Number of samples')
xticks((1:6)*chr_x)
xticklabels({'1N','2N','3N','4N','5N','6N'})
ylim([0,1120])
set(gca,'FontSize',19)
%breakyaxis([220 1100]);
hold off;
%% Gain/Loss plot compared to fixed base ploidy
chr_x= 23-1;
gain = zeros(1, chr_x);
loss = zeros(1, chr_x);
base_ploidy=2;
%base_ploidy=mode(karyotype(:,1:chr_x),2);

%%%%%%%%

for k= 1:chr_x
    gain(k) = sum( karyotype(:,k) > base_ploidy );
    loss(k) = sum( karyotype(:,k) < base_ploidy );
end
%%%%%%%%

num_dipl= sum(sum(karyotype(:,1:22)==2*ones(1,chr_x),2)==chr_x); %number of diploid samples
num_total =  size(karyotype ,1) - num_dipl ; % remove diploid


figure();
hold on;
bar(1:22,gain/num_total,'FaceColor','r')
bar(1:22,-loss/num_total,'FaceColor','b')
hold off
xlabel('Chromosome')
ylabel('Chr. Gains/losses')
legend('Trisomy and up','Monosomy')
xticks(1:chr_x)
ylim([-1,1])
txt_title = strcat('gain/loss from ploidy=',num2str(base_ploidy));
title(txt_title)
set(gca,'FontSize',19)

%% Macro-Karyotype plot

chr_x= 23-1;
p = zeros(chr_x+1,chr_x+1);

xi = zeros(size(karyotype,1),2);
xi(1:size(karyotype,1),1) = sum(karyotype==2,2);
xi(1:size(karyotype,1),2) = sum(karyotype==3,2);

for i2 = 2:(chr_x+2)
    for i3 = 2:(chr_x-i2+4)
        x = [i2-2,i3-2];
        p(x(1)+1,x(2)+1) = sum(sum(x==xi,2)==2);
    end
end
for x2= 2:chr_x+1
    for x3= (chr_x+1-x2+2):(chr_x+1)
        p(x2,x3)=NaN;
    end
end

baseColor = [0 0.4470 0.7410]; %blue
%baseColor =[1 0 1]; %magenta

m = 256;
v = linspace(0.1, 1, m)';
map = baseColor.*v + (1-v);


p(p==0)=NaN;
figure();
set(gcf, 'Position',  [0, 0, 869, 1801]) %square axis
%h = heatmap(flip(p.'), 'MissingDataColor', [1 1 1],'ColorScaling','log','ColorLimits',[0 4]);
%h = heatmap(flip(p.'), 'MissingDataColor', [1 1 1],'ColorLimits',[0 0.001*max(p,[],'all')]);
h = heatmap(flip(p.'), 'MissingDataColor', [1 1 1],'ColorScaling','log');
%caxis([1 100]);
caxis(log([1 100]));

h.GridVisible = 'off';
%h.CellLabelColor = h.Colormap(1,:);
h.CellLabelColor = [0 0 0];
set(gca,'FontSize',7)
h.YDisplayLabels(:) = num2cell(chr_x:-1:0);
h.XDisplayLabels(:) = num2cell(0:chr_x);
%HS = struct(h);
%HS.Axes.Box = 'off';
hHeatmap = struct(h).Heatmap;
hGrid = struct(hHeatmap).Grid;
hGrid.ColorData = uint8([238;238;238;125]);
colormap(map);
colorbar
% xlabel('x_2')
% ylabel('x_3')
%% distance from binary 
xi = zeros(size(karyotype,1),7);
xi(1:size(karyotype,1),1) = sum(karyotype==1,2);
xi(1:size(karyotype,1),2) = sum(karyotype==2,2);
xi(1:size(karyotype,1),3) = sum(karyotype==3,2);
xi(1:size(karyotype,1),4) = sum(karyotype==4,2);
xi(1:size(karyotype,1),5) = sum(karyotype==5,2);
xi(1:size(karyotype,1),6) = sum(karyotype==6,2);
xi(1:size(karyotype,1),7) = sum(karyotype>6,2);
num_dipl=sum(xi(:,2)==22);

[xii,sortIdx]= sort(xi,2,'descend');
udaljenost= 22-sum(xii(:,1:2),2);

figure('color','white');
n=zeros(12,1);
for idx=0:11
n(idx+1)=(sum(udaljenost==idx) );
end
n(1)=n(1)-num_dipl;
b=bar(0:11, n/(size(karyotype,1)-num_dipl));
%b=bar(0:11, n);

b.FaceColor = [1 0.3 0.5];
b.EdgeColor = [1 0.3 0.5];
box off;
xlabel('Distance from Binary karyotypes')
ylabel('Percentage')
set(gca,'FontSize',19)


%% 3N peak in MK
p = zeros(chr_x+1,chr_x+1);
total_chr_copy_num = sum(karyotype,2);
kar2=karyotype(((total_chr_copy_num>55) & (total_chr_copy_num<77)),:);
xi = zeros(size(kar2,1),2);
xi(1:size(kar2,1),1) = sum(kar2==2,2);
xi(1:size(kar2,1),2) = sum(kar2==3,2);

for i2 = 2:(chr_x+2)
    for i3 = 2:(chr_x-i2+4)
        x = [i2-2,i3-2];
        p(x(1)+1,x(2)+1) = sum(sum(x==xi,2)==2);
    end
end
for x2= 2:chr_x+1
    for x3= (chr_x+1-x2+2):(chr_x+1)
        p(x2,x3)=NaN;
    end
end
%baseColor = [0 0.4470 0.7410]; %blue
baseColor = [255 192 203]/255; %pink


p(p==0)=NaN;

m = 256;
v = linspace(0.1, 1, m)';
map = baseColor.*v + (1-v);


figure();
set(gcf, 'Position',  [0, 0, 869, 1801]) %square axis
%h = heatmap(flip(p.'), 'MissingDataColor', [1 1 1],'ColorScaling','log','ColorLimits',[0 4]);
%h = heatmap(flip(p.'), 'MissingDataColor', [1 1 1],'ColorLimits',[0 0.001*max(p,[],'all')]);
h = heatmap(flip(p.'), 'MissingDataColor', [1 1 1],'ColorScaling','log');
%caxis([1 100]);
caxis(log([1 150]));

h.GridVisible = 'off';
%h.CellLabelColor = h.Colormap(1,:);
h.CellLabelColor = [0 0 0];
set(gca,'FontSize',7)
h.YDisplayLabels(:) = num2cell(chr_x:-1:0);
h.XDisplayLabels(:) = num2cell(0:chr_x);
%HS = struct(h);
%HS.Axes.Box = 'off';
hHeatmap = struct(h).Heatmap;
hGrid = struct(hHeatmap).Grid;
hGrid.ColorData = uint8([238;238;238;125]);
colormap(map);
colorbar

% xlabel('x_2')
% ylabel('x_3')
