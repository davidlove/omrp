function CompareNonoverlapping( p, m, k )
%COMPAREOVERLAPPING compares computational time in overlapping and
%nonoverlapping to get equivalent variance reduction

if nargin < 3
    k = 30;
    if nargin < 2
        m = 1000;
    end
end

plot_options = {'MarkerSize',8, 'LineWidth',2};
label_options = {'FontSize', 16};
LineDesign = {'b-', 'g--', 'm-', 'r--', 'c-', ...
    'y--', 'kp:', 'bh:', 'g+:', 'r*:', 'cx:', 'mo:'};
legend_text = cell(1,length(p));

gammabar = (1:m)/m;

for ii = 1:length(p)
T = gammabar.^p(ii);
overlap = T.*NumBatches(gammabar,m,k*m);

plot(gammabar,(overlap/k-1)*100,LineDesign{ii}, plot_options{:})
legend_text{ii} = ['Overlapping, p = ' num2str(p(ii))];
hold on    
end

nb_effective = k./RelativeVariance(gammabar);
plot(gammabar,(nb_effective/k-1)*100,'k-', plot_options{:})

xlabel('\gamma/m', label_options{:})
ylabel('Work increase over base MRP (%)', label_options{:})
ylim([0,200])

legend(legend_text{:}, 'Nonoverlapping', 'Location','NorthEast')

hold off

end

