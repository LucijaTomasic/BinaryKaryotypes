%% TCGA TP53
%%% load oncoprint file downloaded from cBioPortal
clear
%% copy num data
path = strcat('/Documents/ ... /ASCAT'); % path to folder where summary.ascatv3TCGA.penalty70.hg38.tsv file is located
cd (path)
Files = dir('*.tsv');
clinical_data=readtable(Files.name,"FileType","text");
names= clinical_data.name( (clinical_data.purity<1.000) & (clinical_data.purity>0.4)); %ALL DATA

%%
path = strcat('/Documents/ .. /ASCAT/segments 3'); % path to folder where segment txt files are located
cd (path)
Files = dir('*.txt');
numfiles = length(names);

patient_table = cell(1, numfiles);
karyotype=zeros(numfiles, 23);
%% calculate chromosome copy num as mode of segmental copy numbers
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

%% P53

path = strcat('/Documents/ .. /TP53.tsv'); % path to oncoprint file downloaded from cBioPortal

opts = detectImportOptions(path,'FileType','delimitedtext'); 
opts.VariableNamingRule= "preserve";
TP53data=readtable(path,opts);

TP53_alternation =zeros(numfiles, 1);

for s = 1:numfiles
    if ismember(names{s}, TP53data.Properties.VariableNames)
        if isa(TP53data.(names{s})(2:4),'cell') 
            
            if any(~cellfun('isempty',  TP53data.(names{s})(5:9))) 
                TP53_alternation(s)=1;    
            end

        elseif any(~isnan(TP53data.(names{s})(2:4))) 
            error('Error. \n check var type.')
        else
            TP53_alternation(s)=NaN;
        end
    else
        TP53_alternation(s)=NaN;
    end
end

%% MK plot

chr_x= 23-1;
p = zeros(chr_x+1,chr_x+1);

kar2=karyotype(TP53_alternation==1,:);

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



baseColor = [0 0.4470 0.7410]; %blue
%baseColor = [1 0.3412 0.2]; % orange

m = 256;
v = linspace(0.1, 1, m)';
map = baseColor.*v + (1-v);


p(p==0)=NaN;
figure();
set(gcf, 'Position',  [0, 0, 869, 1801]) %square axis
h = heatmap(flip(p.'), 'MissingDataColor', [1 1 1],'ColorScaling','log');
caxis(log([1 200]));


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
%% MK plot of percent of samples with TP53 mutations
chr_x= 23-1;
p1 = zeros(chr_x+1,chr_x+1);
p2 = zeros(chr_x+1,chr_x+1);

kar1=karyotype((TP53_alternation==1),:);
xi = zeros(size(kar1,1),2);
xi(1:size(kar1,1),1) = sum(kar1==2,2);
xi(1:size(kar1,1),2) = sum(kar1==3,2);

for i2 = 2:(chr_x+2)
    for i3 = 2:(chr_x-i2+4)
        x = [i2-2,i3-2];
        p1(x(1)+1,x(2)+1) = sum(sum(x==xi,2)==2);
    end
end

kar2=karyotype((TP53_alternation==0),:);
xi = zeros(size(kar2,1),2);
xi(1:size(kar2,1),1) = sum(kar2==2,2);
xi(1:size(kar2,1),2) = sum(kar2==3,2);

for i2 = 2:(chr_x+2)
    for i3 = 2:(chr_x-i2+4)
        x = [i2-2,i3-2];
        p2(x(1)+1,x(2)+1) = sum(sum(x==xi,2)==2);
    end
end

p_sum=p1+p2;
%p=p2-p1;
p=p1./p_sum;
%p=1-p;

for x2= 2:chr_x+1
    for x3= (chr_x+1-x2+2):(chr_x+1)
        p(x2,x3)=NaN;
    end
end


baseColor1 = [0 0.4470 0.7410];
baseColor2 = [1 0.3412 0.2];
m = 256;
v = linspace(0.1, 1, m)';
map1 = baseColor1.*v + (1-v);
map2 = baseColor2.*v + (1-v);
map=[map1(end:-1:1,:);map2];

figure();
set(gcf, 'Position',  [0, 0, 869, 1801]) %square axis
h = heatmap(round(flip(p.'),2)*100, 'MissingDataColor', [1 1 1]);
caxis([0 100]);


h.GridVisible = 'off';
%h.CellLabelColor = h.Colormap(1,:);
h.CellLabelColor = [0 0 0];
set(gca,'FontSize',10)
h.YDisplayLabels(:) = num2cell(chr_x:-1:0);
h.XDisplayLabels(:) = num2cell(0:chr_x);

hHeatmap = struct(h).Heatmap;
hGrid = struct(hHeatmap).Grid;
hGrid.ColorData = uint8([238;238;238;125]);

colormap(map);

% xlabel('x_2')
% ylabel('x_3')
%% MK Monosomi-disomy
chr_x= 23-1;
p1 = zeros(chr_x+1,chr_x+1);
p2 = zeros(chr_x+1,chr_x+1);

kar1=karyotype((TP53_alternation==1) ,:);
xi = zeros(size(kar1,1),2);
xi(1:size(kar1,1),1) = sum(kar1==1,2);
xi(1:size(kar1,1),2) = sum(kar1==2,2);

for i2 = 2:(chr_x+2)
    for i3 = 2:(chr_x-i2+4)
        x = [i2-2,i3-2];
        p1(x(1)+1,x(2)+1) = sum(sum(x==xi,2)==2);
    end
end

kar2=karyotype((TP53_alternation==0),:);
xi = zeros(size(kar2,1),2);
xi(1:size(kar2,1),1) = sum(kar2==1,2);
xi(1:size(kar2,1),2) = sum(kar2==2,2);

for i2 = 2:(chr_x+2)
    for i3 = 2:(chr_x-i2+4)
        x = [i2-2,i3-2];
        p2(x(1)+1,x(2)+1) = sum(sum(x==xi,2)==2);
    end
end

p_sum=p1+p2;
%p=p2-p1;
p=p1./p_sum;
%p=1-p;

for x2= 2:chr_x+1
    for x3= (chr_x+1-x2+2):(chr_x+1)
        p(x2,x3)=NaN;
    end
end

baseColor1 = [0 0.4470 0.7410]; %plava
baseColor2 = [1 0.3412 0.2];
% 
m = 256;
v = linspace(0.1, 1, m)';
map1 = baseColor1.*v + (1-v);
map2 = baseColor2.*v + (1-v);
map=[map1(end:-1:1,:);map2];


figure();
set(gcf, 'Position',  [0, 0, 869, 1801]) %square axis
h = heatmap(round(flip(p.'),2)*100, 'MissingDataColor', [1 1 1]);

caxis([0 100]);


h.GridVisible = 'off';
%h.CellLabelColor = h.Colormap(1,:);
h.CellLabelColor = [0 0 0];
set(gca,'FontSize',10)
h.YDisplayLabels(:) = num2cell(chr_x:-1:0);
h.XDisplayLabels(:) = num2cell(0:chr_x);
%HS = struct(h);
%HS.Axes.Box = 'off';
hHeatmap = struct(h).Heatmap;
hGrid = struct(hHeatmap).Grid;
hGrid.ColorData = uint8([238;238;238;125]);

% colorbar
colormap(map);

% xlabel('x_2')
% ylabel('x_3')

%% dist from binary 
xi = zeros(size(karyotype,1),7);
xi(1:size(karyotype,1),1) = sum(karyotype==1,2);
xi(1:size(karyotype,1),2) = sum(karyotype==2,2);
xi(1:size(karyotype,1),3) = sum(karyotype==3,2);
xi(1:size(karyotype,1),4) = sum(karyotype==4,2);
xi(1:size(karyotype,1),5) = sum(karyotype==5,2);
xi(1:size(karyotype,1),6) = sum(karyotype==6,2);
xi(1:size(karyotype,1),7) = sum(karyotype>6,2);
% xi(1:size(karyotype,1),7) = sum(karyotype==7,2);
% xi(1:size(karyotype,1),8) = sum(karyotype==8,2);
% xi(1:size(karyotype,1),9) = sum(karyotype>8,2);

[xii,sortIdx]= sort(xi,2,'descend');
udaljenost= sum(xii(:,3:end),2);


%exclude diploid and tetraploid from stat
p53_alternation1=TP53_alternation;
udaljenost1=udaljenost;
p53_alternation1((xi(:,4)==22)|(xi(:,2)==22))=[];
udaljenost1((xi(:,4)==22)|(xi(:,2)==22))=[];

figure('color','white');
y=zeros(11,1);
n=zeros(11,1);
for idx=0:10
y(idx+1)=sum(p53_alternation1(udaljenost1==idx)==1)/(sum(udaljenost1==idx)-sum( isnan(p53_alternation1(udaljenost1==idx)) ) );
n(idx+1)=(sum(udaljenost1==idx)-sum( isnan(p53_alternation1(udaljenost1==idx)) ) );
end
b=bar(0:10, y);
% b.FaceColor = [1 0.3412 0.2];
% b.EdgeColor = [1 0.3412 0.2];
box off;
xlabel('Distance from Binary karyotypes')
ylabel('Percentage of TP53 altered samples')
set(gca,'FontSize',19)


%% postptak binarnih mutiranih
%dijag=sum(p53_alternation( (sum(xi(:,2:3),2)==22) )==0 ,"all" )/(sum(p53_alternation( (sum(xi(:,2:3),2)==22)  )==0 ,"all" ) + sum(p53_alternation( (sum(xi(:,2:3),2)==22) )==1 ,"all" ));
dijag=sum(TP53_alternation( (sum(xi(:,2:3),2)==22) &(xi(:,3)>0) )==0 ,"all" )/(sum(TP53_alternation( (sum(xi(:,2:3),2)==22) &(xi(:,3)>0) )==0 ,"all" ) + sum(TP53_alternation( (sum(xi(:,2:3),2)==22) &(xi(:,3)>0))==1 ,"all" ));
horiz12=sum(TP53_alternation( (sum(xi(:,1:2),2)==22) &(xi(:,1)>0) )==0 ,"all" )/(sum(TP53_alternation( (sum(xi(:,1:2),2)==22) &(xi(:,1)>0) )==0 ,"all" ) + sum(TP53_alternation( (sum(xi(:,1:2),2)==22) &(xi(:,1)>0) )==1 ,"all" ));
ver34=sum(TP53_alternation( (sum(xi(:,3:4),2)==22) &(xi(:,3)>0) )==0 ,"all" )/(sum(TP53_alternation( (sum(xi(:,3:4),2)==22) &(xi(:,3)>0) )==0 ,"all" ) + sum(TP53_alternation( (sum(xi(:,3:4),2)==22)&(xi(:,3)>0) )==1 ,"all" ));
wgd24=sum(TP53_alternation( (sum(xi(:,[2,4]),2)==22)&(xi(:,2)>0) &(xi(:,4)>0)  )==0 ,"all" )/(sum(TP53_alternation((sum(xi(:,[2,4]),2)==22)&(xi(:,2)>0)&(xi(:,4)>0)   )==0 ,"all" ) + sum(TP53_alternation( (sum(xi(:,[2,4]),2)==22)&(xi(:,2)>0)&(xi(:,4)>0)   )==1 ,"all" ));
tetraploid=sum(TP53_alternation( (sum(xi(:,4),2)==22) )==0 ,"all" )/(sum(TP53_alternation( (sum(xi(:,4),2)==22)  )==0 ,"all" ) + sum(TP53_alternation( (sum(xi(:,4),2)==22) )==1 ,"all" ));
diploid=sum(TP53_alternation( (sum(xi(:,2),2)==22) )==0 ,"all" )/(sum(TP53_alternation( (sum(xi(:,2),2)==22)  )==0 ,"all" ) + sum(TP53_alternation( (sum(xi(:,2),2)==22) )==1 ,"all" ));
wgd45=sum(TP53_alternation( (sum(xi(:,[4,5]),2)==22)&(xi(:,5)>0) )==0 ,"all" )/(sum(TP53_alternation((sum(xi(:,[4,5]),2)==22)&(xi(:,5)>0)   )==0 ,"all" ) + sum(TP53_alternation( (sum(xi(:,[4,5]),2)==22)&(xi(:,5)>0)  )==1 ,"all" ));
% wgd56=sum(p53_alternation( (sum(xi(:,[5,6]),2)==22) )==0 ,"all" )/(sum(p53_alternation((sum(xi(:,[5,6]),2)==22)  )==0,"all" ) + sum(p53_alternation( (sum(xi(:,[5,6]),2)==22)  )==1 ,"all" ));
% wgd48=sum(p53_alternation( (sum(xi(:,[4,8]),2)==22)&(xi(:,4)>0) &(xi(:,8)>0)  )==0 ,"all" )/(sum(p53_alternation((sum(xi(:,[4,8]),2)==22)&(xi(:,4)>0)&(xi(:,8)>0)   )==0 ,"all" ) + sum(p53_alternation( (sum(xi(:,[4,8]),2)==22)&(xi(:,4)>0)&(xi(:,8)>0)   )==1 ,"all" ));


%Distance_from_Binary_karyotype= sum(xii(:,3:end),2);
Distance_from_Binary_karyotype=udaljenost;
mix=sum(TP53_alternation(Distance_from_Binary_karyotype>0)==0 ,"all" )/(sum(TP53_alternation(Distance_from_Binary_karyotype>0)==0 ,"all" ) + sum(TP53_alternation(Distance_from_Binary_karyotype>0)==1 ,"all" ));
mix_with_x1=sum(TP53_alternation((Distance_from_Binary_karyotype>0)&(xi(:,1)>0))==0 ,"all" )/(sum(TP53_alternation((Distance_from_Binary_karyotype>0)&(xi(:,1)>0))==0 ,"all" ) + sum(TP53_alternation((Distance_from_Binary_karyotype>0)&(xi(:,1)>0))==1 ,"all" ));


figure();
y=[1-diploid,diploid;1-tetraploid,tetraploid;1-mix, mix;1-dijag, dijag;1-wgd24, wgd24;1-horiz12,horiz12; 1-ver34, ver34];
bar( y,'stacked' ) 

xticklabels({'diploid','tetraploid','Mix','dijag','2&4','horiz12','ver34'})
