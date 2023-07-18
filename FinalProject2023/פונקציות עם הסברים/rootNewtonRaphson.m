function [x,n] = rootNewtonRaphson(f,df, initialGuess, tolerance, maxIterations)
    % Newton-Raphson method for finding the root of a function
    %INPUT: f- the function
    %       df- the derivative
    %       initial guess
    %       tolerance- criterion of convergence
    %       max iterations
    %OUTPUT: x- X value of the function
    %        n- number of iterations

    % Set the initial guess
    x = initialGuess;
    
    % Iterate until the maximum number of iterations is reached
    for n = 1:maxIterations
        % Calculate the function value and its derivative at the current point
        y =         f(x)  ; %num
        yprime =    df(x) ; %den
        
        % Check if the derivative is close to zero (division by zero)
        if abs(yprime) < 100*eps
            error('Derivative is close to zero. Cannot continue iteration.');
        end
        
        % Update the guess using Newton-Raphson formula
        x_new = x - y / yprime;
        clear y yprime;
        
        % Check if the tolerance condition is satisfied
        % desire stop only over here.
        if abs(x_new - x) < tolerance
            return; % Stop when the result is within the desired tolerance
        end
        
        % Update the guess for the next iteration
        x = x_new;
    end
    
    % If the maximum number of iterations is reached without convergence
    error('Maximum number of iterations reached without convergence.');
end
