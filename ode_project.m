
classdef ode_project
    % repos of functions used for the project
   methods (Static)
       function [t, Y] = adam_bashforth2(f, tspan, X0, n)
           % Implementation of the Adams-Bashforth 2nd Order method 
            % INPUTS:
            %   f - function handle or RHS
            %   tspan - time interval 
            %   X0 - the initial state
            %   n - number of points or nodes >= 2
            % OUTPUTS:
            %   t - time 
            %   Y - the solution
            
           
           dt = (tspan(2)-tspan(1))/n;
           t = tspan(1):dt:tspan(2);
           
           Y(:,1)=X0(:);
           
           % obtain X1 via RK4
           k1 = f(t(1), Y(:,1));
           k2 = f(t(1)+dt/2, Y(:,1)+dt/2*k1);
           k3 = f(t(1)+dt/2, Y(:,1)+k2*dt/2);
           k4 = f(t(1)+dt, Y(:,1)+dt*k3);
           Y(:,2)=Y(:,1)+dt*(k1/6+(k2+k3)/3+k4/6);
           
           % the remaining points
           for i=2:n
              Y(:, i+1) = Y(:,i)+(1/2)*dt*(3*f(t(i), Y(:,i))-f(t(i-1),Y(:,i-1))); 
           end
       end
       
       function [t, Y] = adam_bashforth4(f, tspan, X0, n)
           % Implementation of Adams-Bashhforth 4th Order method => requires n>=4
           % INPUTS:
            %   f - function handle or RHS
            %   tspan - time interval 
            %   X0 - the initial state
            %   n - number of points or nodes >= 4
           % OUTPUTS:
            %   t - time steps
            %   Y - the solution
           dt = (tspan(2)-tspan(1))/n;
           t = tspan(1):dt:tspan(2);
           
           Y(:, 1) = X0(:);
          
           % obtain X2to4 via RK4
           for m=1:3
               k1 = f(t(m), Y(:, m));
               k2 = f(t(m)+dt/2, Y(:, m)+(dt/2)*k1);
               k3 = f(t(m)+dt/2, Y(:, m)+k2*dt/2);
               k4 = f(t(m)+dt, Y(:, m)+dt*k3);
               Y(:, m+1)=Y(:,m)+dt*(k1/6+(k2+k3)/3+k4/6);
           end
           
           % the remaining points
           for i=4:n
              Y(:, i+1) = Y(:,i)+(1/24)*dt*(55*f(t(i), Y(:,i))-59*f(t(i-1), Y(:,i-1)) + 37*f(t(i-2), Y(:,i-2))-9*f(t(i-3), Y(:,i-3))); 
           end
           
       end
       
        function [t, Y] = RK2(f, tspan, X0, n)
            % Implementation of the Runge-Kutta 2nd Order Method 
            % INPUTS:
            %   f - function handle or RHS
            %   tspan - time interval 
            %   X0 - the initial state
            %   n - number of points or nodes
            % OUTPUTS:
            %   t - time steps 
            %   Y - the solution
            dt = (tspan(2)-tspan(1))/n;
            t = tspan(1):dt:tspan(2);
            Y(:, 1)=X0(:);
            for i=1:n
                Y(:, i+1) = Y(:, i)+dt*f(t(i)+dt/2, Y(:,i)+(dt/2)*f(t(i), Y(:,i)));
            end
            
        end

        function [t, Y]=RK3(f, tspan, X0, n)
            % Implementation of Runge-Kunta 3rd Order Methhod
            % INPUTS:
            %   f - function handle or RHS
            %   tspan - time interval 
            %   X0 - the initial state
            %   n - number of points or nodes
            % OUTPUTS:
            %   t - time steps 
            %   Y - the solution
            dt = (tspan(2)-tspan(1))/n;
            t = tspan(1):dt:tspan(2);
            Y(:, 1)=X0(:);
            for i=1:n
                k1=f(t(i), Y(:,i));
                k2=f(t(i)+dt/3, Y(:,i)+dt/3*k1);
                k3=f(t(i)-dt/3, Y(:, i)-dt*k1+2*dt/3*k2);
                Y(:,i+1)=Y(:,i)+dt*((-2)*k1+9/4*k2+3/4*k3);
            end
            
        end

        function [t, Y]=RK4(f, tspan, X0, n)
            % IMplementation of Runge-Kutta 4th Order Method
            % INPUTS:
            %   f - function handle or RHS
            %   tspan - time interval 
            %   X0 - the initial state
            %   n - number of points or nodes
            % OUTPUTS:
            %   t - time steps 
            %   Y - the solution
            dt = (tspan(2)-tspan(1))/n;
            t = tspan(1):dt:tspan(2);
            Y(:, 1)=X0(:);
            for i=1:n
                k1 = f(t(i), Y(:,i));
                k2 = f(t(i)+dt/2, Y(:,i)+dt/2*k1);
                k3 = f(t(i)+dt/2, Y(:, i)+k2*dt/2);
                k4 = f(t(i)+dt, Y(:,i)+dt*k3);
                Y(:,i+1)=Y(:,i)+dt*(k1/6+(k2+k3)/3+k4/6);
            end
 
        end

        % Euler Method
        function [t, Y] = Euler_ode(f, X0, tspan, n)
            % Implementation of Euler method 
            % INPUTS:
            %   f - function handle or RHS
            %   tspan - time interval 
            %   X0 - the initial state
            %   n - number of points or nodes
            % OUTPUTS:
            %   t - time steps 
            %   Y - the solution 
            Y(:, 1) = X0(:);
            dt = (tspan(2)-tspan(1))/n;
            t = tspan(1):dt:tspan(2);           
            for i=1:n
                Y(:, i+1) = Y(:, i) + f(t(i), Y(:, i))*dt;
            end
           
        end
        
        
        function []=portrait_plot(f, lim)
            % plot phase portrait and solution
            x1 = linspace(-lim, lim, 20);
            x2 = linspace(-lim, lim, 20);
            [x, y] = meshgrid(x1, x2);
            size(x); 
            size(y);

            u = zeros(size(x)); % value of the derivative of X(1) at each point x, y
            v = zeros(size(x)); % value of the derivative of X(2) at each point x, y

            % at t=0
            t=0;
            for i=1:numel(x)
                X_prime = f(t, [x(i); y(i)]);
                u(i) = X_prime(1);
                v(i) = X_prime(2);
            end

            % plot
            figure;
            quiver(x, y, u, v, 'black')
            xlabel('x(t)')
            ylabel('y(t)')
            axis tight equal
            title('Phase Portrait')
        end

        function []=plot_solution(S, color, annotation)
            % plots the solution
            %INPUTS:
            % S - the solution obtained from the nunmerical methods
            % color - the color of the plot
            % annnotation - name of the method used
            %OUTPUTS:
            % plot
            hold on
            plot(S(1,1), S(2, 1), 'gx', 'MarkerSize', 15, 'LineWidth', 5)
            plot(S(1, 2:end), S(2, 2:end), color, 'LineWidth', 2);
            legend('phase portrait','Starting point',annotation)
            axis tight equal
            hold off
        end

        function [t, lin_alg_subs]=linear_exact_solution(A, X0, tspan, n)
            % implementation of the analytical method for the linear system 
            % using linear algebra
            %INPUTS:
            % A - the matrix corresponding to the system
            % X0 - the initial state
            % tspan - time interval
            % n - thhe number of nodes
            %OUTPUTS:
            % t - time steps
            % lin_alg_subs - the solution
            
            dt = (tspan(2)-tspan(1))/n;
            t = tspan(1):dt:tspan(2);
            [V, D] = eig(A);
            V1 = V(:, 1);
            V2 = V(:, 2);
            lamdas = diag(D); 
            lam1 = lamdas(1);
            lam2 = lamdas(2);
            C = [V1 V2]\X0; c1 = C(1); c2 = C(2);
            lin_alg = @(t) real(c1*V1*exp(lam1*t)+c2*V2*exp(lam2*t))+imag(c1*V1*exp(lam1*t)+c2*V2*exp(lam2*t));
            lin_alg_subs=lin_alg(t);
        end    

        function []=plot_xy(results)
            % plots of the indepenndent results with respect to the
            % independable variable, here t.
            %INPUTS:
            % results - metrix of the results from the numerical methods
            %OUTPUTS:
            % plots
            
            figure;
            ax1 = subplot(2, 1, 1);
            ax2 = subplot(2, 1, 2);
            % x vs time
            x = results(1, :);
            plot(ax1, x)
            title(ax1, 'x vs time')
            ylabel(ax1, 'x(t)')
            xlabel(ax1, 'time')
            % y vs time
            y = results(2, :);
            plot(ax2, y)
            title(ax2, 'y vs t')
            ylabel(ax2, 'y(t)')
            xlabel(ax2, 'time')
        end

        function []=plot_comparison(exact, sol, tspan, n, comparison)
            % comparivative plot of the solution and approximation with
            % respect to the inndedpendable variable, here t.
            %INPUTS:
            % exact - analytical solution connsidered or "true" solution
            % sol - results from the numerical methods
            % n - number of nodes
            % comparison - name of the method, e.g 'Euler method'
            %OUTPUTS:
            % plot
            
            figure;
            ax1 = subplot(2, 1, 1);
            ax2 = subplot(2, 1, 2);
            dt = (tspan(2)-tspan(1))/n;
            x_axis = tspan(1):dt:tspan(2);
            plot(ax1, x_axis, exact(1, :), 'g', x_axis, sol(1, :), 'r')
            legend(ax1, 'Analytical solution', comparison)
            ylabel(ax1, 'x(t)')
            plot(ax2, x_axis, exact(2,:), 'g', x_axis, sol(2,:), 'r')
            legend(ax2, 'Analytical solution', comparison)
            ylabel(ax2, 'y(t)')
            xlabel('time')
        end
        
        function [] = plot_linear_vs_nlinear(lin_sol, nlin_sol, tspan, n)
            % comparivative plot of the solutions from both linear and nonlinear 
            % system with respect to the inndedpendable variable, here t.
            %INPUTS:
            % lin_sol - analytical solution of the linear system
            % nlin_sol - analytical solution of the nonlinear system
            % tspan - time interval
            % n - number of nodes
            %OUTPUTS:
            % plot
            
            figure;
            ax1 = subplot(2, 1, 1);
            ax2 = subplot(2, 1, 2);
            dt = (tspan(2)-tspan(1))/n;
            x_axis = tspan(1):dt:tspan(2);
            plot(ax1, x_axis, lin_sol(1,:), 'g', x_axis, nlin_sol(1,:), 'r')
            legend(ax1, 'Solution to linear system', 'Solution to non linear system')
            ylabel(ax1, 'x(t)')
            plot(ax2, x_axis, lin_sol(2,:), 'g', x_axis, nlin_sol(2,:), 'r')
            legend(ax2, 'Solution to linear system', 'Solution to non linear system')
            ylabel(ax2, 'y(t)')
            xlabel('time')
        end

        function [err]=error(analytical_res, sol)
            % computes the norm of the error between the analytical solution
            % and the approximation from the numerical methods
            %INPUTS:
            % analytical_res - analytical solution or "true" solution
            % sol - results from the numerical methods
            %OUTPUTS:
            % err - the norm error  
            
            format long
            % norm of (analytical_res - sol)./analytical_res =>relative
            % error. Uncomment the line bellow to consider it.
            % err = norm((analytical_res-sol)./analytical_res);
            
            err = norm(analytical_res-sol);
        end

        function [errs]=err_on_linearSystem_vs_n(A, X0, f, tspan, ns, method)
            % computes the errors of with respect to the number of nodes -
            % case: solutions of the linear system
            %INPUTS:
            % A - the matrix corresponding to the system
            % X0 - the initial state
            % tspan - time interval
            % ns - array of number of nodes
            % methods - numerical method to use. e.g ("RK2" for Runge-Kutta 2nd)
            %OUTPUTS:
            % errs - array of errors with respect to theh number of nodes.
            
            errs = [];
            for n=ns
                [~, analytical_sol] = ode_project.linear_exact_solution(A, X0, tspan, n);
                if strcmp(method, 'Euler')
                    [~, approx] = ode_project.Euler_ode(f, X0, tspan, n);
                    anno = 'Loglog plot of the error of Euler method';
                    
                
                elseif strcmp(method, 'RK2')
                    [~, approx] = ode_project.RK2(f, tspan, X0, n);
                    anno = 'LogLog plot of the Error of RK 2nd order Method on linear system';
                
                
                elseif strcmp(method, 'RK3')
                    [~, approx] = ode_project.RK3(f, tspan, X0, n);
                    anno = 'LogLog plot of the Error of RK 3rd order Method on linear system';
               
                elseif strcmp(method, 'ab2')
                    [~, approx] = ode_project.adam_bashforth2(f, tspan, X0, n);
                    anno = 'LogLog plot of the Error of Adam Bashforth 2nd Order Method on linear system';
                    
                elseif strcmp(method, 'ab4')
                    [~, approx] = ode_project.adam_bashforth4(f, tspan, X0, n);
                    anno = 'LogLog plot of the Error of Adam Bashforth 4th Order Method on linear system';
                
                else
                    method = 'RK4';
                    [~, approx] = ode_project.RK4(f, tspan, X0, n);
                    anno = 'LogLog plot of the Error of RK 4th order Method on linear system';
                end
                err = ode_project.error(analytical_sol, approx);
                errs = [errs; err];
            end
            figure;
            loglog(ns, errs)
            ylabel('error')
            xlabel('n - number of nodes')
            title(anno)
        end
        
        function [t, sol] = fixed_step_ode45(f, tspan, X0, n)
            %Implementation of the analytical method for the nonlinear
            %system using ode45 for fixed step size
            % INPUTS:
            %   f - function handle or RHS
            %   tspan - time interval 
            %   X0 - the initial state
            %   n - number of points or nodes
            % OUTPUTS:
            %   t - time steps 
            %   sol - the solution
           dt = (tspan(2)-tspan(1))/n;
           t = tspan(1):dt:tspan(2);
           sol(:, 1) = X0(:);
           for i = 2:length(t)
              [~, x_ode] = ode45(f, [0, dt], sol(:, end));
              sol(:, i) = x_ode(end, :)';
           end
        end
        
        function [errs]=err_nonlinearSystem_vs_n(f, X0, tspan, ns, method)
            % computes the errors of with respect to the number of nodes -
            % case: solutions of the nonlinear system
            %INPUTS:
            % f - the nonlinear function handle or RHS
            % X0 - the initial state
            % tspan - time interval
            % ns - array of number of nodes
            % methods - numerical method to use. e.g ("RK2" for Runge-Kutta 2nd)
            %OUTPUTS:
            % errs - array of errors with respect to theh number of nodes.
            
            errs = [];
            for n=ns
                [~, analytical_nl_res] = ode_project.fixed_step_ode45(f, tspan, X0, n);
                        
                if strcmp(method,'Euler')
                    [~, approx] = ode_project.Euler_ode(f, X0, tspan, n);
                    
                    anno = 'LogLog plot of the Error of Euler Method on Non-linear system';
               
                
                elseif strcmp(method, 'RK2')
                    [~, approx] = ode_project.RK2(f, tspan, X0, n);
                    anno = 'LogLog plot of the Error of RK 2nd order Method on Non-linear system';
                
                
                elseif strcmp(method, 'RK3')
                    [~, approx] = ode_project.RK3(f, tspan, X0, n);
                    anno = 'LogLog plot of the Error of RK 3rd order Method on Non-linear system';
                
                elseif strcmp(method, 'ab2')
                    [~, approx] = ode_project.adam_bashforth2(f, tspan, X0, n);
                    anno = 'LogLog plot of the Error of Adam Bashforth 2nd Order Method on Non-linear system';
                
                elseif strcmp(method, 'ab4')
                    [~, approx] = ode_project.adam_bashforth4(f, tspan, X0, n);
                    anno = 'LogLog plot of the Error of Adam Bashforth 4th Order Method on Non-linear system';
                
                else
                    method = 'RK4';
                    [~, approx] = ode_project.RK4(f, tspan, X0, n);
                    anno = 'LogLog plot of the Error of RK 4th order Method on Non-linear system';
                end
                err = ode_project.error(analytical_nl_res, approx);
                errs = [errs; err];
            end
            figure;
            loglog(ns, errs)
            ylabel('error')
            xlabel('n - number of nodes')
            title(anno)
        end
        
        %%% TODO => future work: implementation of adaptive step size
        %%% methods

   end
end
