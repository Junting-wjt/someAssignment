%-------------------------------------------------
% MAFTDSP Matlab Assignment 2 - Part Four
% 
% General Lagrange interpolation
% 
% Junting Wang 10/12/23
%-------------------------------------------------


function output_matrix = MA2_s2327978_Wang_Linterp(N, Q, fmode)


% restrict N to be even ---------------------------------------------------
    if mod(N, 2) ~= 0
        error('N must be even.'); 
    end


% Initialize the alpha range for estimate ---------------------------------
    alpha_q = linspace(-0.5, 0.5, Q); 


% Generate alpha values for N interpolation points ------------------------
    alpha = linspace(- (N - 1) / 2, (N - 1) / 2, N);


% Choose different mode by the value of fmode -----------------------------
    % lookup table generation mode
    if fmode == 1
        % preallocate a matrix for the lookup table
        lookup_table = zeors(Q, N);

        % for each evaluation point
        for i = 1:Q 
            % current evaluation point alpha value
            a_q = alpha_q(i);

            % for each interpolation point
            for j = 1:N 
                % current interpolation point alpha value
                a_j = alpha(j);           
                P = 1; 

                % compute Lagrange Basis Polynomials
                for q = alpha 
                    if q ~= a_j 
                        P = P * (a_q - q) / (a_j - q);
                    end
                end

                % store the result in the appropriate location of the lookup table
                lookup_table(i, j) = P;
            end
        end
        output_matrix = lookup_table;


    % polynomial plotting mode
    elseif fmode == 2
        % generate alpha values for Q evaluation points
        alpha_q = linspace(- (N - 1) / 2, (N - 1) / 2, Q);
        figure;

        % for each interpolation point
        for j = 1:N 
            a_j = alpha(j);
            P_n_alpha = zeros(1, Q); 

            % for each evaluation point
            for i = 1:Q 
                  a_q=alpha_q(i);
                  P = 1;

                % compute Lagrange Basis Polynomials
                for q = alpha
                    if q ~= a_j 
                    P = P * (a_q - q) / (a_j - q);
                    end
                end
                P_n_alpha(i) = P;
            end

% plot the polynomial -----------------------------------------------------
            plot(alpha_q, P_n_alpha); 
            hold on;
        end
        hold off; 
        output_matrix = [];
        title('Lagrange Basis Polynomials'); 
        xlabel('\alpha'); 
        ylabel('P_n(\alpha)'); 
        legend_labels = arrayfun(@(m) sprintf('P_{%.2f}^{%.2f}(\\alpha)', alpha(m), N), 1:N, 'UniformOutput', false);
        legend(legend_labels, 'Location', 'best');
        grid on;
        hold off;

    else
        % Throw an error if fmode is not 1 or 2
        error('Invalid fmode. Use 1 for lookup table or 2 for plotting polynomials.'); 
    end
end