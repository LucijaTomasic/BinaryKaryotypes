%% Mitelman Cytoconverted data
clear
%% load MItelman data

path = strcat('/Documents/ .. /data'); %%% path to folder where downloaded data is located
cd (path)
Files = dir('*.TXT.DATA');
opts = detectImportOptions('CYTOCONVERTED.TXT.DATA','FileType','delimitedtext'); 
opts = setvartype(opts,{'CaseNo'},'char');
opts_Clinical = detectImportOptions('CYTOGEN.TXT.DATA','FileType','delimitedtext'); 
opts_Clinical = setvartype(opts_Clinical,{'CaseNo'},'char');
opts_original = detectImportOptions('KCLONE.TXT.DATA','FileType','delimitedtext'); 
opts_original = setvartype(opts_original,{'CaseNo'},'char');

converted_data=readtable('CYTOCONVERTED.TXT.DATA',opts);
error_log=tdfread('CYTOCONVERTED_LOG.TXT.DATA','tab');
Clinical_data= readtable('CYTOGEN.TXT.DATA',opts_Clinical);
original_kar= readtable('KCLONE.TXT.DATA',opts_original );
Koder_data= readtable('KODER.TXT.DATA',"FileType",'delimitedtext');

%% definira Topo Morph
%RefCaseNo=Clinical_data((Clinical_data.Topo==305)&(Clinical_data.Morph==3131),1:2); %Bladder Transitional cell carcinoma 
%RefCaseNo=Clinical_data((Clinical_data.Topo==305),1:2); %Bladder all
%RefCaseNo=Clinical_data((Clinical_data.Topo==801)&(Clinical_data.Morph==5101),1:2); %Brain Glioma
%RefCaseNo=Clinical_data((Clinical_data.Topo==801)|(Clinical_data.Topo==806),1:2); %Brain and brain stem
%RefCaseNo=Clinical_data((Clinical_data.Topo==401)&((Clinical_data.Morph==3111)|(Clinical_data.Morph==3102)),1:2); %breast carcinoma and adenocarcinoma
%RefCaseNo=Clinical_data((Clinical_data.Topo==505)|(Clinical_data.Topo==507)|(Clinical_data.Topo==509),1:2); %Uterus and Vagina all
%RefCaseNo=Clinical_data((Clinical_data.Topo==222),1:2); %Gallbladder/Biliary system
%RefCaseNo=Clinical_data(((Clinical_data.Topo==225)|(Clinical_data.Topo==227))&((Clinical_data.Morph==3111)|(Clinical_data.Morph==3102)),1:2); %LS intestine carcinoma
%RefCaseNo=Clinical_data((Clinical_data.Topo>200)&(Clinical_data.Topo<213),1:2); %oral
%RefCaseNo=Clinical_data((Clinical_data.Morph==1820),1:2); % Diffuse large B-cell lymphoma
%RefCaseNo=Clinical_data((Clinical_data.Morph==8810),1:2); % Mesothelioma
%RefCaseNo=Clinical_data((Clinical_data.Topo==501),1:2); %Ovary
%RefCaseNo=Clinical_data((Clinical_data.Topo==218),1:2); %Pancreas
%RefCaseNo=Clinical_data((Clinical_data.Topo==703),1:2); %Adrenal
%RefCaseNo=Clinical_data((Clinical_data.Topo==602),1:2); %
%RefCaseNo=Clinical_data((Clinical_data.Morph==6011),1:2); % Melanoma
%RefCaseNo=Clinical_data((Clinical_data.Topo==216),1:2); %Stomach
%RefCaseNo=Clinical_data((Clinical_data.Topo==601),1:2); %Testis
%RefCaseNo=Clinical_data((Clinical_data.Topo==1301),1:2); %Thymus
%RefCaseNo=Clinical_data((Clinical_data.Topo==704),1:2); %Thyroid
%RefCaseNo=Clinical_data((Clinical_data.Morph>=8000)&(Clinical_data.Morph<=8899),1:2);%Sarcoma

%%%%%%%%%%%%%%%%%%%%%%          Leuk.           %%%%%%%%%%%%%%%%%%%%%%% 
%RefCaseNo=Clinical_data((Clinical_data.Morph>=1000)&(Clinical_data.Morph<=1099),1:2);%AUBL
%RefCaseNo=Clinical_data((Clinical_data.Morph>=1100)&(Clinical_data.Morph<=1199),1:2);%AML
%RefCaseNo=Clinical_data((Clinical_data.Morph>=1200)&(Clinical_data.Morph<=1299),1:2);%CML
%RefCaseNo=Clinical_data((Clinical_data.Morph>=1600)&(Clinical_data.Morph<=1699),1:2);%ALL 
%RefCaseNo=Clinical_data((Clinical_data.Morph==1802),1:2);%CLL
%RefCaseNo=Clinical_data((Clinical_data.Morph==1708)|(Clinical_data.Morph==1804)|(Clinical_data.Morph==1818)|(Clinical_data.Morph==1822)|(Clinical_data.Morph==1902)|(Clinical_data.Morph==1904)|(Clinical_data.Morph==1914)|(Clinical_data.Morph==1926),1:2);

%RefCaseNo=Clinical_data((Clinical_data.Morph==1708),1:2);%PCL
%RefCaseNo=Clinical_data((Clinical_data.Morph==1804),1:2);%B-PLL
%RefCaseNo=Clinical_data((Clinical_data.Morph==1818),1:2);%HCL
%RefCaseNo=Clinical_data((Clinical_data.Morph==1822),1:2);%BL
%RefCaseNo=Clinical_data((Clinical_data.Morph==1902),1:2);%T-PLL
%RefCaseNo=Clinical_data((Clinical_data.Morph==1904),1:2);%T-GLL
%RefCaseNo=Clinical_data((Clinical_data.Morph==1914),1:2);%ANL
%RefCaseNo=Clinical_data((Clinical_data.Morph==1926),1:2);%ATL

%%%%%%%%%%%%%%%%%%%%%%          lymphoma            %%%%%%%%%%%%%%%%%%%%%%% 
%RefCaseNo=Clinical_data((Clinical_data.Morph>=2000)&(Clinical_data.Morph<=2005),1:2);%Hodgkin disease (all subtypes)
%RefCaseNo=Clinical_data((Clinical_data.Morph>=1800)&(Clinical_data.Morph<=1899),1:2);%B-cell(all subtypes)
%RefCaseNo=Clinical_data((Clinical_data.Morph>=1900)&(Clinical_data.Morph<=1999),1:2);%T-cell(all subtypes)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RefCaseNo=Clinical_data((Clinical_data.Morph>=5100)&(Clinical_data.Morph<=5299),1:2);%
%RefCaseNo=Clinical_data((Clinical_data.Morph>=9000)&(Clinical_data.Morph<=9099),1:2);% Nonneoplastic disorders


%RefCaseNo=Clinical_data((Clinical_data.Morph>=3000)&(Clinical_data.Morph<=3099),1:2);%pretumor
%RefCaseNo=Clinical_data(Clinical_data.Topo==501,1:2); % Ovarian
%RefCaseNo=Clinical_data(Clinical_data.Topo==214,1:2); % Oesophagus/ oesophageal tumors 
%RefCaseNo=Clinical_data((Clinical_data.Topo<600)&(Clinical_data.Topo>500),1:2); % All female genital organs 
%RefCaseNo=Clinical_data(Clinical_data.Topo==401,1:2); % Breas
%RefCaseNo=Clinical_data((Clinical_data.Topo==225)|(Clinical_data.Topo==227)|(Clinical_data.Topo==216)|(Clinical_data.Topo==230),1:2); % intestine stomac and colorectal
%RefCaseNo=Clinical_data((Clinical_data.Topo==227)&((Clinical_data.Morph>3100)&(Clinical_data.Morph<3900)),1:2); % carcinoma
%RefCaseNo=Clinical_data((Clinical_data.Topo==227)&((Clinical_data.Morph>3000)&(Clinical_data.Morph<3100)),1:2); % adenoma
%RefCaseNo=Clinical_data((Clinical_data.Topo<300)&(Clinical_data.Topo>200)&((Clinical_data.Morph>3000)&(Clinical_data.Morph<3100)),1:2); % adenoma digestive
%RefCaseNo=Clinical_data((Clinical_data.Topo<300)&(Clinical_data.Topo>200)&((Clinical_data.Morph>3100)&(Clinical_data.Morph<3900)),1:2); % carcinoma digestive
%RefCaseNo=Clinical_data(Clinical_data.Morph==1602,1:2); % ALL leukemia
%RefCaseNo=Clinical_data((Clinical_data.Morph<1204)|(Clinical_data.Morph==1602),1:2); %all data Leukemia
%RefCaseNo=Clinical_data(~((Clinical_data.Morph<1204)|(Clinical_data.Morph==1602)),1:2); %all EXCEPT Leukemia
%RefCaseNo=Clinical_data((Clinical_data.Morph>2900)&(Clinical_data.Morph<3200),1:2); %Malignant and tumor epithelial neoplasms
%RefCaseNo=Clinical_data(strcmp('Y' ,Clinical_data.PrevTum),1:2); 
%RefCaseNo=Clinical_data((Clinical_data.Topo==227)&( ~((Clinical_data.Morph==3003)|(Clinical_data.Morph==3099))),1:2);
%RefCaseNo=Clinical_data((((Clinical_data.Morph==3001)|(Clinical_data.Morph==3003)|(Clinical_data.Morph==3011)|(Clinical_data.Morph==3011)|(Clinical_data.Morph==3041)|(Clinical_data.Morph==3055)|(Clinical_data.Morph==3099))),1:2); %Benign

RefCaseNo=Clinical_data(:,1:2); % ALL

%% FILTER

%%%%%%%%%%%%%%%%%%%
FilterOut=original_kar((endsWith(original_kar.CloneShort,'inc')),1:4); %svi koji zavrsavaju na inc, ukupno 3340
FilterOut3 = original_kar((original_kar.ChromoMin~=original_kar.ChromoMax),1:4); % svi koji nemaju tocno odreden ukupni br kromosoma, ukupno  12355
FilterOut2=FilterOut; % svi koji za koje cytoconvereter javlja warning, ukupno 2718
i=0;
for jj=1:size(error_log.RefNo,1)
    if strcmp(error_log.Message(1,:),error_log.Message(jj,:))
        i=i+1;
        FilterOut2.RefNo(i)=error_log.RefNo(jj);
        FilterOut2.CaseNo(i)={strtrim(error_log.CaseNo(jj,:))};
        FilterOut2.InvNo(i)=error_log.InvNo(jj);
        FilterOut2.CloneNo(i)=error_log.Clone(jj);
    end
end
FilterOut2(i+1:end,:)=[];
%size(FilterOut2)
FilterOut=[FilterOut;FilterOut2;FilterOut3];
%FilterOut=[FilterOut;FilterOut2];
%size(FilterOut)
% ostavlja samo odredene clinical data
i=zeros(height(FilterOut),1);
for jj=1:height(RefCaseNo)
   i=i + (( FilterOut.RefNo==RefCaseNo.RefNo(jj))&( strcmp(RefCaseNo.CaseNo(jj),FilterOut.CaseNo)  )); 
end
FilterOut(~i,:)=[];
% size(FilterOut)

%% Karyotype

total_chr_length=[248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468]; 	
resolution=1000;
numfiles=100000; 
karyotype=zeros(numfiles, 22); 
sample_i=1;
preskoci=0;
idxx_for_table=FilterOut;
idx=1;
for j=1:height(RefCaseNo)
    %RefCaseNo(j,:)
    convdata_i=converted_data(converted_data.RefNo==RefCaseNo.RefNo(j),:);
    convdata_i=convdata_i(strcmp(RefCaseNo.CaseNo(j),convdata_i.CaseNo),:);
    
    %     [~,idx]=sort(convdata_i.ChrOrd);
    %     convdata_i=convdata_i(idx,:);
    %     [~,idx]=sort(convdata_i.Clone);
    %     convdata_i=convdata_i(idx,:);
    %     [~,idx]=sort(convdata_i. CaseNo);
    %     convdata_i=convdata_i(idx,:);
    %     convdata_i
    InvNoClone=unique( convdata_i(:,3:4), 'rows');
    
    for InvNoClone_i= 1:height(InvNoClone)
        sample=convdata_i(convdata_i.InvNo==InvNoClone.InvNo(InvNoClone_i),:);
        sample=sample(sample.Clone==InvNoClone.Clone(InvNoClone_i),:);

        if any((FilterOut.RefNo==sample.RefNo(1))&(strcmp(sample.CaseNo(1),FilterOut.CaseNo))&(FilterOut.InvNo==sample.InvNo(1))&(FilterOut.CloneNo==sample.Clone(1)))
            preskoci=preskoci+1;
        else
            for chr_i=1:22
                bins=2*ones(resolution,1);
                chr_data=sample(sample.ChrOrd==chr_i,:);
                for i=1:height(chr_data) 
                    if strcmp('Loss',chr_data.Type(i))
                        s =  round( chr_data.Start(i)*resolution/total_chr_length(chr_i) ) +1;
                        e =  round( chr_data.End(i)*resolution/total_chr_length(chr_i) ) ;
                        if e>resolution
                            e=resolution;
                        end
                        bins(s:e) = bins(s:e) - 1;
                    end
                    if strcmp('Gain',chr_data.Type(i))
                        s =  round( chr_data.Start(i)*resolution/total_chr_length(chr_i) ) +1;
                        e =  round( chr_data.End(i)*resolution/total_chr_length(chr_i) ) ;
                        if e>resolution
                            e=resolution;
                        end
                        bins(s:e) = bins(s:e) + 1;
                    end
                end
                 karyotype( sample_i ,chr_i) = mode(bins);

            end
            idxx_for_table.RefNo(idx)=sample.RefNo(1);
            idxx_for_table.CaseNo(idx)=sample.CaseNo(1);
            idxx_for_table.InvNo(idx)=sample.InvNo(1);
            idxx_for_table.CloneNo(idx)=sample.Clone(1);
            
            idx=idx+1;
            if( (sum(karyotype(sample_i ,:)==3*ones(1,22),2)==19))
                msg = strcat('refNo_i=',num2str(RefCaseNo.RefNo(j)),'\n caseNo_i=',RefCaseNo.CaseNo{j,:},'\n');
                fprintf(msg)
                
            end
            sample_i = sample_i +1;
        end
       
    end
end

karyotype(sample_i:end,:)=[];


%% Total chr num 1N 2N 3N 4N 
chr_x= 23-1;

total_chr_copy_num = sum(karyotype,2);
total_chr_copy_num = sort(total_chr_copy_num);

figure();
hold on;
histogram(total_chr_copy_num,1000)
hold off
xlabel('Ploidy')
ylabel('Number of samples')
xticks((1:6)*chr_x)
xticklabels({'1N','2N','3N','4N','5N','6N'})
xlim([20,96])
ylim([0,300])

set(gca,'FontSize',19)
%% Gain/Loss plot 
chr_x= 23-1;
gain = zeros(1, chr_x);
loss = zeros(1, chr_x);
%base_ploidy=mode(karyotype,2);
base_ploidy=2;
%%%%%%%% 
for k= 1:chr_x
    gain(k) = sum( karyotype(:,k) > base_ploidy );
    loss(k) = sum( karyotype(:,k) < base_ploidy );
end
%%%%%%%%

num_dipl= sum(sum(karyotype(:,1:22)==2*ones(1,chr_x),2)==chr_x); %number of diploid samples
num_total =  size(karyotype ,1) - num_dipl ; % remove diploid


figure();
set(gcf, 'Position',  [0, 0, 1000, 350])
hold on;
bar(1:22,gain/num_total,'FaceColor','r')
bar(1:22,-loss/num_total,'FaceColor','b')
hold off
xlabel('Chromosome')
ylabel('Chr. Gains/losses')
legend('Gain','Loss')
xticks(1:chr_x)
ylim([-0.15,0.25])
set(gca,'FontSize',19)
%% Distance from Binary
xi = zeros(size(karyotype,1),7);
xi(1:size(karyotype,1),1) = sum(karyotype==1,2);
xi(1:size(karyotype,1),2) = sum(karyotype==2,2);
xi(1:size(karyotype,1),3) = sum(karyotype==3,2);
xi(1:size(karyotype,1),4) = sum(karyotype==4,2);
xi(1:size(karyotype,1),5) = sum(karyotype==5,2);
xi(1:size(karyotype,1),6) = sum(karyotype==6,2);
xi(1:size(karyotype,1),7) = sum(karyotype>6,2);
%xi(xi(:,2)==22,:)=[]; % remove diploid

xi= sort(xi,2,'descend');

udaljenost_ukupno= sum(xi(:,3:end),2);
figure('color','white');
h=histogram(udaljenost_ukupno,'BinMethod', 'integers');
box off;
xlim([-0.99 5])
xlabel('Distance from Binary karyotypes')
ylabel('Number of samples')
set(gca,'FontSize',19)


%% Macro-Karyotype plot 
chr_x= 23-1;
p = zeros(chr_x+1,chr_x+1);

xi = zeros(size(karyotype,1),1);
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

baseColor = [0 0.4470 0.7410]; 


m = 256;
v = linspace(0.1, 1, m)';
map = baseColor.*v + (1-v);


p(p==0)=NaN;
figure();
set(gcf, 'Position',  [0, 0, 869, 1801]) 
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

%% seq random 400 stanica i klasifikacija na 3 ruba i bulk

color= generatecolormapthreshold([0 1 2 3 4 5 6 7 8],[0 0.3 0.85;0 0.4470 0.7410; 0.9020 0.9059 0.9020;1 0.302 0.302;0.6350 0.0780 0.1840;0.2863 0 0;0.2863 0 0;0.2863 0 0]); 
color2= generatecolormapthreshold([-1 0 1 2 3 4 5 6 7 8],[1 1 1;0 0.3 0.85;0 0.4470 0.7410; 0.9020 0.9059 0.9020;1 0.302 0.302;0.6350 0.0780 0.1840;0.2863 0 0;0.2863 0 0;0.2863 0 0]); 
karyotype(:,23)=0;
chr=23;
chr_x= 23-1;
[~,N_idx]=sort(sum(karyotype,2));
karyotype =karyotype(N_idx,:);

br_plot= 400; 
Ctm= karyotype(~(sum(karyotype(:,1:chr_x)==2*ones(1,chr_x),2)==chr_x),:); % remove diploid
Ctm = Ctm(randperm(size(Ctm,1),br_plot),:);
[~,N_idx]=sort(sum(Ctm,2));
Ctm =Ctm(N_idx,:);
%%%
%%%%%%%%%%%%%%%%%
distances = pdist(Ctm, 'euclidean');
linkage_tree = linkage(distances, 'ward');  
leaf_order = optimalleaforder(linkage_tree, distances);
Ctm = Ctm(leaf_order, :);
%%%%%%%%%%%%%%%%%%
%%%%%% definira debljinu stupaca
rowWidths = 10*ones(1,br_plot);  
chrSize= [113 110 90 87 82 78 72 67 62 62 60 61 52 49 46 42 37 37 27 29 21 23 72];
colWidths = (chrSize+1)/10;  %mm, one for each col of data; applied left to right
%%%%%%%
fig = figure();

set(gcf, 'Position',  [20, 20, 1500, 400])
ax = axes(fig);
[xg,yg] = meshgrid([0,cumsum(colWidths)], [0,cumsum(rowWidths)]);
sh = surf(xg, yg, zeros(size(xg)),Ctm, 'MeshStyle', 'column');

colormap(ax,color)
cb = colorbar(ax);
caxis([0 8]) % color limit colormap
view(2)
set(gca,'FontSize',14)
xlim([0 sum(colWidths)])
xticksPosition = cumsum(colWidths)-colWidths/2;
xticks(xticksPosition)

xlabel('')
ylabel('')
yticks([])
kk=["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X"];
xticklabels(kk)
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
chr_x= 23-1;
xi = zeros(size(Ctm,1),2);
xi(1:size(Ctm,1),1) = sum(Ctm==2,2);
xi(1:size(Ctm,1),2) = sum(Ctm==3,2);

Ctm23 = Ctm( (sum(xi,2)==chr_x),:);
Ctm34 = Ctm( (xi(:,1)==0),:); % x2=0 vertikala
Ctm12 = Ctm( ((xi(:,2)==0)&(xi(:,1)~=0)),:); 
Ctm_bulk = Ctm((~(sum(xi,2)==chr_x))&(xi(:,1)~=0)&(xi(:,2)~=0),:);
%%%%%%%%% Sorting and clustering
distances = pdist(Ctm23, 'euclidean');
linkage_tree = linkage(distances, 'ward');  
leaf_order = optimalleaforder(linkage_tree, distances);
Ctm23 = Ctm23(leaf_order, :);

distances = pdist(Ctm34, 'euclidean');
linkage_tree = linkage(distances, 'ward');  
leaf_order = optimalleaforder(linkage_tree, distances);
Ctm34 = Ctm34(leaf_order, :);

distances = pdist(Ctm12, 'euclidean');
linkage_tree = linkage(distances, 'ward');  
leaf_order = optimalleaforder(linkage_tree, distances);
Ctm12 = Ctm12(leaf_order, :);

distances = pdist(Ctm_bulk, 'euclidean');
linkage_tree = linkage(distances, 'ward');  
leaf_order = optimalleaforder(linkage_tree, distances);
Ctm_bulk = Ctm_bulk(leaf_order, :);

b=-1*ones(5,23);
Ctmm=[Ctm23;b;Ctm12;b;Ctm34;b;Ctm_bulk];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% definira debljinu stupaca
rowWidths = 10*ones(1,size(Ctmm,1));  
chrSize= [113 110 90 87 82 78 72 67 62 62 60 61 52 49 46 42 37 37 27 29 21 23 72];
colWidths = (chrSize+1)/10;  
%%%%%%%
fig = figure();

set(gcf, 'Position',  [20, 20, 1500, 400])
ax = axes(fig);
[xg,yg] = meshgrid([0,cumsum(colWidths)], [0,cumsum(rowWidths)]);
surf(xg, yg, zeros(size(xg)),Ctmm(end:-1:1,:), 'MeshStyle', 'column');;

colormap(ax,color2)
cb = colorbar(ax);
caxis([-1 8]) % color limit colormap
view(2)
set(gca,'FontSize',14)
xlim([0 sum(colWidths)])
xticksPosition = cumsum(colWidths)-colWidths/2;
xticks(xticksPosition)

xlabel('')
ylabel('')
yticks([])
kk=["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X"];
xticklabels(kk)


%% Percentage of binary kar
% Ctm= karyotype(~(sum(karyotype(:,1:chr_x)==2*ones(1,chr_x),2)==chr_x),:); % remove diploid
% xi1 = zeros(size(Ctm,1),2);
% xi1(1:size(Ctm,1),1) = sum(Ctm==2,2);
% xi1(1:size(Ctm,1),2) = sum(Ctm==3,2);
% sum((sum(xi1,2)==chr_x),1)%/size(Ctm,1) %23
% 
% xi2 = zeros(size(Ctm,1),2);
% xi2(1:size(Ctm,1),1) = sum(Ctm==3,2);
% xi2(1:size(Ctm,1),2) = sum(Ctm==4,2);
% sum((sum(xi2,2)==chr_x),1)%/size(Ctm,1) %34
% 
% xi3 = zeros(size(Ctm,1),2);
% xi3(1:size(Ctm,1),1) = sum(Ctm==1,2);
% xi3(1:size(Ctm,1),2) = sum(Ctm==2,2);
% sum((sum(xi3,2)==chr_x),1)%/size(Ctm,1) %12
% 
% xi3 = zeros(size(Ctm,1),2);
% xi3(1:size(Ctm,1),1) = sum(Ctm==4,2);
% xi3(1:size(Ctm,1),2) = sum(Ctm==5,2);
% sum((sum(xi3,2)==chr_x),1)%/size(Ctm,1) %45


%% colormap fun

function color_value = generatecolormapthreshold(threshold,colour)
    arguments
        threshold (:,:) {mustBeNumeric}
        colour (:,:) {mustBeNumeric, mustBeGreaterThanOrEqual(colour,0), mustBeLessThanOrEqual(colour,1)}
    end
    validateattributes(threshold,{'numeric'},{'size', [1 NaN]},mfilename,"Threshold",2)
    validateattributes(colour,{'numeric'},{'size', [NaN 3]},mfilename,"Colour",2)
    if length(threshold)-1 ~= size(colour,1)
        error("threshold and colour are not match. Either more colour or more threshold is assigned")
    else
        color_value=[];
        for i = 1:1:length(threshold)-1
            createrange = threshold(i):threshold(i+1);
            for w = 1:1:length(createrange)-1
               color_value = [color_value;colour(i,:)];
            end
        end
    end
end