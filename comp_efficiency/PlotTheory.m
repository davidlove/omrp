function PlotTheory(p,m,n)

if nargin < 2
    m = 1000;
end

if nargin < 3
    n = 1000*m;
end

plot_options = {'MarkerSize',8, 'LineWidth',2};
label_options = {'FontSize', 16};
LineDesign = {'b-', 'g-', 'm-', 'r-', 'c-', ...
    'y<:', 'kp:', 'bh:', 'g+:', 'r*:', 'cx:', 'mo:'};
legend_text = cell(1,length(m));

for ii = 1:length(p)
    legend_text{ii} = ['p = ' num2str(p(ii))];
    [cd,gb] = ComputationalEfficiency(p(ii),m,n);
    
    plot(gb, cd, LineDesign{ii}, plot_options{:})
    hold on
end
plot([0,1],[1,1],'k')
xlabel('\bar{\gamma}', label_options{:})
ylabel('Relative Computational Efficienty', label_options{:})
legend(legend_text{:}, 'Location','Best')

hold off
