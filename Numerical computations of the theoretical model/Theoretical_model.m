%%% Theoretical model numerical solutions 
%%% Euler method

clear
chr=22;
kor = 1/(2*log(2));


p_miss= 0.002*kor;
pgd = @(x) 0.005;
p_MP= 0.002;

apop = @(x) 0.15-0.15*( ((x(3)+x(4))==chr) || ((x(3)+x(2))==chr) || ((x(2)+x(1))==chr)  )*( 1-0.10*( x(1)>0 ) ) ;
%apop = @(x) 0;

beta=@(x) log(2);

T =500;
dt = 0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = 0;
time = [0:10:100,T];   %trenutci t u kojima se spremju vrijednosti n i ntot
%time=0:10:T;

%%% Initialisation
nnew = zeros(chr+3,chr+3,chr+3,chr+3);
nold = zeros(chr+3,chr+3,chr+3,chr+3);

%%%% INITIAL CONDITIONS
nold(2,2+chr,2,2) = 1; %%%Diploid
%nold(2,2,2,2+chr) = 1; %%%Tetraploid


ii = 1;

p1 = zeros(chr+1,chr+1);

while t<=T
    for i1 = 2:(chr+2)
        for i2 = 2:(chr-i1+4)
            for i3 = 2:(chr-i1-i2+6)
                for i4 = 2:(chr-i1-i2-i3+8)
                    x = [i1-2,i2-2,i3-2,i4-2,chr-i4-i3-i2-i1+8];
                    
                    if ( x(1) == 0 && x(3) ==0 && x(5) ==0 ) 
                        ngd = nold(i2, i4, 2, 2); 
                    else
                        ngd = 0;
                    end
                    
                    if ( x(5)==0 ) %MP triploid
                        
                        p_multi =3* p_MP *(factorial(chr)/factorial(x(1))/factorial(x(2))/factorial(x(3))/factorial(x(4)) * (8*(1/3)^(4))^(x(1)) * (24*(1/3)^(4))^(x(2)) * (32*(1/3)^(4))^(x(3)) * ((2/3)^(4))^(x(4)) );
                        
                    else
                        p_multi = 0;
                    end
                    
                    if  ( x(2)==chr )  
                        MP=p_MP;
                    else
                        MP=0;
                    end
                    
                    F = (1 - 2*(sum((1:5).*x)*p_miss  + apop(x) +pgd(x) -MP  ) )*beta(x)*nold(i1,i2,i3,i4) ...
                        +(  (x(1)+1)*beta([x(1)+1,x(2)-1,x(3),x(4),x(5)])*nold(i1+1,i2-1,i3,i4) ...
                        + 2*(x(2)+1)* ( beta([x(1)-1,x(2)+1,x(3),x(4),x(5)])*nold(i1-1,i2+1,i3,i4) + beta([x(1),x(2)+1,x(3)-1,x(4),x(5)])*nold(i1,i2+1,i3-1,i4) )...
                        + 3*(x(3)+1)* ( beta([x(1),x(2)-1,x(3)+1,x(4),x(5)])*nold(i1,i2-1,i3+1,i4) + beta([x(1),x(2),x(3)+1,x(4)-1,x(5)])*nold(i1,i2,i3+1,i4-1) )...
                        + 4*(x(4)+1)* ( beta([x(1),x(2),x(3)-1,x(4)+1,x(5)])*nold(i1,i2,i3-1,i4+1) + beta([x(1),x(2),x(3),x(4)+1,x(5)-1])*nold(i1,i2,i3,i4+1) )...
                        + 5*(x(5)+1)*   beta([x(1),x(2),x(3),x(4)-1,x(5)+1])*nold(i1,i2,i3,i4-1)  )*p_miss ...
                        + pgd([x(2),x(4),0,0,0]) * beta([x(2),x(4),0,0,0]) * ngd...
                        + p_multi *beta([0,chr,0,0,0])*nold(2,2+chr,2,2);
                    
                    
                    nnew(i1,i2,i3,i4) = nold(i1,i2,i3,i4) + dt*F;
                    

                end
            end
        end
    end
    
    t = t + dt;
    
    if t >= time(ii)    %save for plot
        p = zeros(chr+1,chr+1);
        for ii2 = 2:(chr+2)
            for ii3 = 2:(chr-ii2+4)
                p(ii2-1,ii3-1) = sum(  nnew(:,ii2,ii3,:), 'all');
            end
        end
        p=p/sum(p, 'all');
        p1=p1+p;
        
        ii = ii + 1;
    end
    nold = nnew;
        
end


%% MK x2 x3

p2=p1;

p2=p2/sum(p2, 'all')*64417; 

for x2= 2:chr+1
    for x3= (chr+1-x2+2):(chr+1)
        p2(x2,x3)=NaN;
    end
end
%%
figure();
h = heatmap(round(flip(p2.')), 'MissingDataColor', [1 1 1],'ColorScaling','log','ColorLimits',[-0.2 8]);
set(gca,'FontSize',4)
set(gcf, 'Position',  [0, 0, 869, 1801])
HS = struct(h);
HS.Axes.Box = 'off'
h.CellLabelColor = [0 0 0];
hHeatmap = struct(h).Heatmap;
hGrid = struct(hHeatmap).Grid;
hGrid.ColorData = uint8([238;238;238;125]);
h.YDisplayLabels(:) = num2cell(chr:-1:0);
h.XDisplayLabels(:) = num2cell(0:chr);
%xlabel('x_2')
%ylabel('x_3')


%% udaljenost od bi-modalnih
udaljenost = zeros(12,1);
for x2 = 0:(chr)
    for x3 = 0:(chr-x2)
        xi=[x2,x3,chr-x2-x3];
        xi= sort(xi,'descend');
        udaljenost(xi(:,3)+1)=udaljenost(xi(:,3)+1)+ p2(x2+1,x3+1);
    end
end
udaljenost(1)= udaljenost(1) - p2(chr+1, 1);

%%% exp
broj_stanica=[ 36423,8107,1832,571,194,98,56,20,11,5,4,2]; %  Mitelman all data 
%broj_stanica=[19361 3847 763 172 59 33 14 6 1 0 0 0]; % Mitelman Leukemias
gain = broj_stanica/sum(broj_stanica);

%%%%%%%%%%%%%%%%%%%%
figure('color','white');
hold on
plot(0:11,udaljenost/sum(udaljenost), 'LineWidth',2,'Color','red'); %theory
plot(0:11, gain,'.','MarkerSize',25); %exp
hold off
box off;
xlim([-0.5 5])
xlabel('Distance from Binary karyotypes')
ylabel('Fraction')
set(gca,'FontSize',19)


