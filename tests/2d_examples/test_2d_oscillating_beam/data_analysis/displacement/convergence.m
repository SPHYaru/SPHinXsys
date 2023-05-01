clear;
clc;
dpi='-r600'; r=400; paper_size_m=[2000 1860]/r; page_size_m=[0 0 paper_size_m];
paper_size_p=[3000 1500]/r; page_size_p=[0 0 paper_size_p];
LegendFontSize = 10;
 flag_p_pdf=0;
string_path_root='D:\SPHinXsys\CompositeMaterial5\Source\tests\2d_examples\test_2d_oscillating_beam\data_analysis\displacement\'; 
stringpre=[string_path_root,'one-displacement-convergence'];
data1=load('D:\SPHinXsys\CompositeMaterial5\Source\tests\2d_examples\test_2d_oscillating_beam\data_analysis\displacement\numerical-20.dat');
data2=load('D:\SPHinXsys\CompositeMaterial5\Source\tests\2d_examples\test_2d_oscillating_beam\data_analysis\displacement\numerical-30.dat');
data3=load('D:\SPHinXsys\CompositeMaterial5\Source\tests\2d_examples\test_2d_oscillating_beam\data_analysis\displacement\numerical-40.dat');
data4=load('D:\SPHinXsys\CompositeMaterial5\Source\tests\2d_examples\test_2d_oscillating_beam\data_analysis\displacement\Ana_num_material_1.txt');


x1=data1(:,1);
y1=data1(:,3);

x2=data2(:,1);
y2=data2(:,3);

x3=data3(:,1);
y3=data3(:,3);

x4=data4(:,1);
y4=data4(:,2);



plot(x1,y1,'b-','LineWidth',1.0);
hold on
plot(x2,y2,'k','LineWidth',1.0);
hold on
plot(x3,y3,'r--','LineWidth',1.0);
hold on
plot(x4,y4,'-.','color',[0.13 0.55 0.13],'LineWidth',1.0);


set(gca,'XLim',[0 1]);%X轴的数据显示范围
set(gca,'YLim',[-0.045,0.03]);%X轴的数据显示范围
% set(gca, 'YTick', -0.025:0.1:0.005);
xlabel('{t (s)}','FontSize',11,'Fontname','Times new roman');
ylabel('{\eta (m)}','FontSize',11,'Fontname','Times new roman');
axis on;  %设置坐标轴开启
set(gca,'FontSize',11,'Fontname','Times new roman');
box on;
h1=legend("Present-ds = t/20 m","Present-ds = t/30 m","Present-ds = t/40 m","Analytical solution");
set(h1,'FontSize',11,'Fontname','Times new roman','position', [0.4,0.19,0.3,0.1])
%set(gcf,'unit','centimeters','position',[3 5 10 6])
%============================输出图片设置=================================
 set(gcf, 'PaperUnits', 'inches', 'PaperPosition', page_size_p);
  if flag_p_pdf==0
     string=[stringpre, '.eps']; print('-depsc',dpi,string)
%    string=[stringpre, string_num, '.tif']; print('-dtiff',dpi,string)
  elseif flag_p_pdf==1
     set(gcf, 'PaperSize', paper_size_p);
     string=[stringpre, '.pdf']; print('-dpdf',dpi,string)
  end
 
