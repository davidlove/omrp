function plot_rdk( k, file_start, file_end, varargin )

plot_options = {'MarkerSize',8, 'LineWidth',2};
label_options = {'FontSize', 16};
LineDesign = {'kx:', 'bo:', 'gd:', 'm>:', 'rs:', 'cx:', ...
    'y<:', 'kp:', 'bh:', 'g+:', 'r*:', 'cx:', 'mo:'};

% NOTE: This assumes that m = 200.  Can be changed in the future.
m = 200;

for ii = 1:length(k)
    clear full;
    full = load([file_start num2str(k(ii)) file_end]);
    full = full';
    gamma = full(1,:);
    cover = full(2,:);
    ci = full(3,:);
    avgvar = full(4,:);
    varvar = full(5,:);
    g = gamma / m;
    
    legend_text{ii} = ['n/m = ' num2str(k(ii))];
    
    for jj=1:size(varargin,2)
        sw_var = varargin{jj};
        switch sw_var
            case 'cover'
                figure(2)
                plot(g,cover, LineDesign{ii}, ...
                    plot_options{:} );
                ylim([0.8,1]);
                xlabel('\gamma/m', label_options{:})
                ylabel('Coverage Probability of CI', label_options{:})
                hold on
            case 'ci'
                figure(3)
                plot(g,ci, LineDesign{ii}, ...
                    plot_options{:} );
                xlabel('\gamma/m', label_options{:})
                ylabel('Confidence Interval Width', label_options{:})
                hold on
            case 'avgvar'
                figure(5)
                plot(g,avgvar, LineDesign{ii}, ...
                    plot_options{:} );
                xlabel('\gamma/m', label_options{:})
                ylabel('E[VG]', label_options{:})
                hold on
            case 'varvar'
                figure(6)
                plot(g,varvar/varvar(end), LineDesign{ii}, ...
                    plot_options{:} );
                xlabel('\gamma/m', label_options{:})
                ylabel('Normalized Var(VG)', label_options{:})
                hold on
        end
    end
end

for jj=1:size(varargin,2)
    sw_var = varargin{jj};
    switch sw_var
        case 'cover'
            figure(2)
            x = [0,1];
            y = [.9,.9];
            plot(x, y, 'k-', ...
                plot_options{:} );
            legend(legend_text{:},'Planned coverage', 'Location','SouthWest');
            hold off
        case 'ci'
            figure(3)
            legend(legend_text{:}, 'Location','SouthWest');
            hold off
        case 'avgvar'
            figure(5)
            legend(legend_text{:}, 'Location','Best');
            hold off
        case 'varvar'
            figure(6)
            N = 1:1000;
            theo = (2*N.^2 + 1) ./ (3*N.^2);
            plot(1./N, theo, 'k-', ...
                plot_options{:} );
            legend(legend_text{:},'Theoretical', 'Location','NorthWest');
            hold off
    end
end
