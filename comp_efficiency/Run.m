figure(1)
PlotTheory([1, 0.9, 0.71, 0.5, 0.2, 0], 1000)

data = importdata('m1000_insertion.dat', '\t', 2);

gammabar = data.data(:,1);
m = max(gammabar);
gammabar = gammabar/m;

vg = data.data(:,7);
comptime = data.data(:,8);

efficiency = 1./(vg.*comptime);
efficiency = efficiency / efficiency(end);

p = 0.711; % Power fit found with R
[theory,gb] = ComputationalEfficiency(p,m,1000*m);

figure(2)
plot(gammabar, efficiency, 'o', gb, theory,'m-', 'LineWidth',2,'MarkerSize',8)
hold on
plot([0,1],[1,1],'k')
xlabel('\gamma/m', 'FontSize',16)
ylabel('Relative Computational Efficiency', 'FontSize',16)
legend('NVD with Insertion Sort, m=1000', 'Theoretical, p=0.71', 'Location','SouthEast')

figure(3)
CompareNonoverlapping([1, 0.9, 0.71, 0.5, 0.2, 0], 1000, 1000)

