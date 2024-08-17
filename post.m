% 使用readmatrix函数读取data.txt文件  
data = readmatrix('shakecpu_walls.txt');  
  
% data现在是一个包含txt文件中所有数据的矩阵

time = data(:,1)';
contact = data(:,2)';


plot(time,contact,'-o')

labs=28;
ls=1.2
fonts=36
set(gca,'FontSize',labs);
set(gca,'linewidth',ls);
xlabel('time','Interpreter','latex','fontsize',fonts,'fontweight','bold');
ylabel('contact number','Interpreter','latex','fontsize',fonts,'fontweight','bold');