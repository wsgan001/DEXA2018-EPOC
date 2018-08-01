% color1 = [255,215,0]/255;%³È
% color4 = [192,80,77]/255;%×Ï
% color7 = [150,200,50]/255;%ÂÌ
% color8 = [230,0,0]/255;%»Æ
% color5 = [75,172,198]/255;%Ç³À¶
% color6 = [150,80,30]/255;%À¶
% color3 = [57,173,72]/255;%ÂÌ2 
% color2 = [255,140,0]/255;%ºì

color7 = [228,26,28]/255;%³È
color4 = [55,126,184]/255;%×Ï
color8 = [77,175,74]/255;%ÂÌ
color1 = [152,78,163]/255;%»Æ
color5 = [255,127,0]/255;%Ç³À¶
color6 = [255,255,51]/255;%À¶
color3 = [166,86,40]/255;%ÂÌ2 
color2 = [247,129,191]/255;%ºì


figure(1)
%% k-coverage
subplot(224)
SVR = [0.51875	0.4875	0.55	0.596875	0.690625	0.721875];
Linear = [0.73125	0.778125	0.7625	0.778125	0.809375	0.840625];
PreWhether = [0.631372549	0.650980392	0.729411765	0.749019608	0.768627451	0.807843137];
Seismic = [0.471875	0.45625	0.471875	0.51875	0.4875	0.55		];
Sansnet = [0.729411765	0.768627451	0.807843137	0.82745098	0.827058824	0.846666667];
EPOC = [0.690196078	0.768888889	0.82744186	0.847058824	0.866666667	0.864736842];
plot(Linear,'-p','color',color1,'linewidth',2,'Markersize',3,'MarkerEdgeColor',color1,'MarkerFaceColor',color1);
hold on ;
plot(SVR,'-o','color',color6,'linewidth',2,'Markersize',3,'MarkerEdgeColor',color6,'MarkerFaceColor',color6);
plot(PreWhether,'-^','color',color3,'linewidth',2,'Markersize',3,'MarkerEdgeColor',color3,'MarkerFaceColor',color3);
plot(Seismic,'-h','color',color4,'linewidth',2,'Markersize',3,'MarkerEdgeColor',color4,'MarkerFaceColor',color4);
plot(Sansnet,'->','color',color8,'linewidth',2,'Markersize',3,'MarkerEdgeColor',color8,'MarkerFaceColor',color8);
plot(EPOC,'-s','color',color2,'linewidth',2,'Markersize',3,'MarkerEdgeColor',color2,'MarkerFaceColor',color2);
% plot(B_ESP,'-p','color',color7,'linewidth',2,'Markersize',3,'MarkerEdgeColor',color7,'MarkerFaceColor',color7);
% plot(FB_ESP,'-h','Color',color5,'linewidth',2,'Markersize',3,'MarkerEdgeColor',color5,'MarkerFaceColor',color5);
set(gca,'xtick',1:1:5)
set(gca, 'GridLineStyle' ,'-')
set(gca,'linewidth',1);
xlabel( "Time (hour)");
ylabel("K-Coverage")
grid on


%% cost
subplot(223)
SVR = [0.124573379	0.129692833	0.107508532	0.099829352	0.094709898	0.086177474];
Linear = [0.090443686	0.099829352	0.103242321	0.103242321	0.087030717	0.074232082];
PreWhether = [0.133959044	0.111774744	0.104095563	0.09556314	0.091296928	0.078498294];
Seismic = [0.163447099	0.172832765	0.178805461	0.167713311	0.169419795	0.151428571];
Sansnet = [0.103242321	0.087030717	0.074232082	0.077645051	0.075085324	0.069112628];
EPOC = [0.087030717	0.08112628	0.072525597	0.06996587	0.064705882	0.065862069];
plot(Linear,'-p','color',color1,'linewidth',2,'Markersize',3,'MarkerEdgeColor',color1,'MarkerFaceColor',color1);
hold on ;
plot(SVR,'-o','color',color6,'linewidth',2,'Markersize',3,'MarkerEdgeColor',color6,'MarkerFaceColor',color6);
plot(PreWhether,'-^','color',color3,'linewidth',2,'Markersize',3,'MarkerEdgeColor',color3,'MarkerFaceColor',color3);
plot(Seismic,'-h','color',color4,'linewidth',2,'Markersize',3,'MarkerEdgeColor',color4,'MarkerFaceColor',color4);
plot(Sansnet,'->','color',color8,'linewidth',2,'Markersize',3,'MarkerEdgeColor',color8,'MarkerFaceColor',color8);
plot(EPOC,'-s','color',color2,'linewidth',2,'Markersize',3,'MarkerEdgeColor',color2,'MarkerFaceColor',color2);
% plot(B_ESP,'-p','color',color7,'linewidth',2,'Markersize',3,'MarkerEdgeColor',color7,'MarkerFaceColor',color7);
% plot(FB_ESP,'-h','Color',color5,'linewidth',2,'Markersize',3,'MarkerEdgeColor',color5,'MarkerFaceColor',color5);
set(gca,'xtick',1:1:5)
set(gca, 'GridLineStyle' ,'-')
set(gca,'linewidth',1);
xlabel( "Time (hour)");
ylabel("Cost")
grid on

legend('Linear','SVR','PreWhether','SEISMIC','Sansnet','EPOC','Orientation','horizontal');
