function plot_rd( m, problem, varargin )

plot_options = {'MarkerSize',8, 'LineWidth',2};
label_options = {'FontSize', 16};
LineDesign = {'kx:', 'bo:', 'gd:', 'm>:', 'rs:', 'cx:', ...
    'y<:', 'kp:', 'bh:', 'g+:', 'r*:', 'cx:', 'mo:'};
legend_text = cell(1,length(m));

% mysql parameters
outfile = '/tmp/matlab_mysql.txt';
system_call  = 'mysql omrp -Be ''';
mysql_select =  'SELECT gamma,cover_prob,ci_width,avg_of_var,var_of_var';
mysql_from   =  ' FROM run_data INNER JOIN ';
mysql_latest = [' ( SELECT a.name,a.problem,a.batch_size ' ...
                '  FROM run_desc AS a WHERE NOT EXISTS ' ...
                '  ( SELECT * FROM run_desc AS b ' ...
                '   WHERE a.problem=b.problem ' ...
                '    AND a.batch_size=b.batch_size ' ...
                '    AND b.date>a.date ' ...
                ' ) ) AS latest ' ];
mysql_on     =  ' ON run_data.name = latest.name';
mysql_where  = [' WHERE latest.problem="' problem '"' ...
                ' AND latest.batch_size = '];
mysql_order  =  ' ORDER BY gamma'''; % Must be ascending order
system_redir = [' > ' outfile];

G = 1;
C = 2;
W = 3;
A = 4;
V = 5;

for ii = 1:length(m)
    system([ system_call ...
        mysql_select mysql_from mysql_latest mysql_on ...
        mysql_where num2str(m(ii)) mysql_order ...
        system_redir]);
    
    clear temp;
    temp   = importdata(outfile);
    
    assert( strcmp( temp.colheaders{G}, 'gamma' ) )
    assert( strcmp( temp.colheaders{C}, 'cover_prob' ) )
    assert( strcmp( temp.colheaders{W}, 'ci_width' ) )
    assert( strcmp( temp.colheaders{A}, 'avg_of_var' ) )
    assert( strcmp( temp.colheaders{V}, 'var_of_var' ) )
    
    gamma  = temp.data(:,G);
    cover  = temp.data(:,C);
    ci     = temp.data(:,W);
    avgvar = temp.data(:,A);
    varvar = temp.data(:,V);
    g = gamma / m(ii);
    
    legend_text{ii} = ['m = ' num2str(m(ii))];
    
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
