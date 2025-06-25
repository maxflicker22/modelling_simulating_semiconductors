% Utilization of method of nested intervals to find the root of a function
% defined on a onedimensional vector
%
% usage : 
%   chemical_potential = 
%   FindRootNestedIntervals(@(E) chargeNeutrality(E,E_C,E_V,m_n_eff,m_p_eff,T),...
%                           energies, (E_C-E_V)/2, 1d-3, 20);
%
%   or
%
%   fh = @(E) chargeNeutralityIntrinsic(E,E_C,E_V,m_n_eff,m_p_eff,T);
%   chemical_potential = FindRootNestedIntervals(fh,energies, ...
%                        (E_C-E_V)/2, 1d-3, 20); 
%
% requires: definition of function FEval
%
%           e.g., FuncEval(x) = chargeNeutralityIntrinsic(energies) 
%           checking for charge neutrality
% input:   
%     Interval        ..  vector containing search interval 
%                         {x} = [xmin,xmax] of target quantity x_target 
%     FuncEval        ..  string quantity f(x) dependent on x
%     initial_guess   ..  initial guess for x_target
%     iter_threshold  ..  threshold how close to f(x)==0 will be iterated 
%                         (until |f(x)| < iter_threshold)
%     max_iter_steps  ..  upper limit of number of interval splitting steps 
%
% output:
%     x_target        ..  one value x_target within interval {x} satisfying 
%                         the condition f(x)==0
%                         x_target is NOT necessarily member of the vector x!
%                       


function [ x_target, num_iter, error ] = ...
       FindRootNestedIntervals( FuncEval,  Interval, initial_Guess, ...
       iter_threshold, max_iter_steps )


x_max = max(Interval);
x_min = min(Interval);


x_target = 0.;

iter =  1;
error = 1;

    while ((abs(error)> iter_threshold) && (iter < max_iter_steps ))
        
        if (iter == 1)
            x_center = initial_Guess;
        else
           x_center = (x_max - x_min)/2 + x_min;
        end

        % error here is an absolute error :(
        error = FuncEval(x_center) ;


        if (error > 0)
             x_max = x_center;
        else
             x_min = x_center;
        end; % if


        iter = iter + 1;
    end; % while

    x_target = x_center;
    num_iter = iter;
end
