function plot_files( m, file_base, file_end )

plot_options = {'MarkerSize',8, 'LineWidth',2};
label_options = {'FontSize', 16};
LineDesign = {'kx:', 'bo:', 'gd:', 'm>:', 'rs:', 'cx:', ...
    'y<:', 'kp:', 'bh:', 'g+:', 'r*:', 'cx:', 'mo:'};
legend_text = cell(1,length(m));

numfiles = length(m);

cell_data = cell(1,numfiles);
for ii = 1:length(m)
    cell_data{ii} = importdata([file_base num2str(m(ii)) file_end]);
    assert( isequal( cell_data{1}.colheaders, cell_data{ii}.colheaders ) )
end
assert( isequal( cell_data{1}.colheaders{1}, 'gamma' ) )

numcols = length(cell_data{1}.colheaders);

for ii = 1:numfiles
    gamma  = cell_data{ii}.data(:,1);
    g = gamma / m(ii);
    
    legend_text{ii} = ['m = ' num2str(m(ii))];
    
    for jj=2:numcols
        figure(jj)
        y_data = cell_data{ii}.data(:,jj);
        
        colhead = cell_data{1}.colheaders{jj};
        switch colhead
            case 'var_of_var'
                plot(g,y_data/y_data(end), LineDesign{ii}, ...
                    plot_options{:} );
                hold on
            otherwise
                plot(g,y_data, LineDesign{ii}, ...
                    plot_options{:} );
                hold on
        end
    end
end

for jj=2:numcols
    figure(jj)
    xlabel('\gamma/m', label_options{:})
    
    colhead = cell_data{1}.colheaders{jj};
    switch colhead
        case 'cover_prob'
            x = [0,1];
            y = [.9,.9];
            plot(x, y, 'k-', ...
                plot_options{:} );
            ylim([0.8,1]);
            ylabel('Coverage Probability of CI', label_options{:})
            legend(legend_text{:},'Planned coverage', 'Location','SouthWest')
        case 'ci_width'
            ylabel('Confidence Interval Width', label_options{:})
            legend(legend_text{:}, 'Location','SouthWest')
        case 'avg_of_var'
            ylabel('E[VG]', label_options{:})
            legend(legend_text{:}, 'Location','Best')
        case 'var_of_var'
            gg = 0.001:0.001:1;
            cg = ceil(1./gg);
            vt1 = -gg + 2*gg.^2.*cg - gg.^3.*cg.^2;
            vt2 = 2*(gg.*cg - gg.^2.*cg.^2);
            vt3 = 2/3*gg.^3.*cg.^3 + 1/3*gg.^3.*cg;
            theo = vt1 + vt2 + vt3;
            plot(gg, theo, 'k-', ...
                plot_options{:} );
            ylabel('Normalized Var(VG)', label_options{:})
            legend(legend_text{:},'Theoretical', 'Location','NorthWest')
        case 'comp_time'
            ylabel('Computation Time (s)', label_options{:})
            legend(legend_text{:}, 'Location','NorthEast')
        otherwise
            warning(['Column type ' colhead ' not recognized'])
            ylabel(['Unknown label "' colhead '"'], 'Interpreter','none',label_options{:})
            legend(legend_text{:}, 'Location','Best')
    end
    hold off
end
