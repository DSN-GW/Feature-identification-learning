% GR adaptations
% test this; it takes smaller and smaller intervals at the end
% does that make sense?

function out = slidingavg_nan(in, N)
nx   = length(in);  

if (isempty(in)) | (N<=0)                                              % If the input array is empty or N is non-positive,
 disp(sprintf('SlidingAvg: (Error) empty input data or N null.'));     % an error is reported to the standard output and the
 return;                                                               % execution of the routine is stopped.
end % if

if (N==1)                                                              % If the number of neighbouring points over which the sliding 
 out = in;                                                             % average will be performed is '1', then no average actually occur and
 return;                                                               % OUTPUT_ARRAY will be the copy of INPUT_ARRAY and the execution of the routine
end % if                                                               % is stopped.

out = zeros(size(in));       

   for i=1:nx,                                              
    if ((i <= 1) & ((i+N) <= nx))                          
       out(i) = nanmean(in(1:N));
    elseif (i >= 1) & ((i+N) <= nx)
       out(i) = nanmean(in(i:N+(i-1)));                                                         % then we proceed to evaluate the mean on 2*m elements centered on the i-th position.
    elseif ((i >= 1) & ((i+N) >= nx))                          % If not enough points are available on the rigth of the i-th element..
       out(i) = nanmean(in(i:nx));                                % then we proceed to evaluate the mean from the element (i - m)-th to the last one.
   % elseif ((i-N) < 1) & ((i+N) >= nx)                          % If not enough points are available on the left and on the rigth of the i-th element..
    %   out(i) = nanmean(in(1:nx));                                  % then we proceed to evaluate the mean from the first element to the last.
     end % if
    end % for i
end

 

