function y = make_random_1D_line(pts,numPoints,smoothness,upper_bound)
start_pt = pts(1);
end_pt = pts(2);
% Define the original signal.
% x = linspace(0, end_pt-start_pt, numPoints); % XX
x = linspace(0, 1, numPoints);
y = x;%zeros(1, numPoints); % Can be whatever you want.  I picked all zeros.

% Define number of iterations, spacing, and how much noise amplitude falls off with iteration number.
numIterations = 50;

% Can't make it more than numPoints because the most we can alter is every one of the numPoints elements.
if numIterations > numPoints
	numIterations = numPoints;
end
% Set the smoothness factor - how fast the noise amplitude dies away with iteration in the exponential decay.
% Smaller numbers give noisier signal.  Bigger numbers give smoother signal.
% smoothness = 2;
% Make linearly decreasing amplitude.
% amplitude = linspace(10, 1, numIterations)
max_amplitude = 10;%10*abs(end_pt-start_pt);
amplitude = linspace(max_amplitude, 1, numIterations);
% Now make it exponential decay to decrease faster.
amplitude = amplitude .* exp(-smoothness .* (1:length(amplitude)));
% Determine how many locations to adjust the y at for each iteration.
loc = round(linspace(5, numPoints, numIterations));

% Now here are the iterations where we actually add noise to the original y signal.
for k = 1 : numIterations
	% We don't add noise to the entire y vector at all locations, just some of them.  
	% An increasing number of locations, meaning a decreasing spacing.
	% Find the exact few locations we will add noise to.
	all_locations = ceil(linspace(1, numPoints, loc(k)));
	locations = all_locations; 
    y_prev = y;
    
	% Add noise.  Noise decreases in amplitude and spacing as iterations progress.
    while ~isempty(locations)
        noise = amplitude(k) * (rand(1, length(locations)) - 0.5);
        y(locations) = y(locations) + noise;
        
        % convert y vals to corresponding volume position 
        y_check = (y(locations)-y(1))*(end_pt-start_pt)+start_pt; 
        % check if any points fall outside volume bounds after noise
        oob = y_check < 1 | y_check > upper_bound; 
        % update locations to any points still out of bounds 
        locations = locations(oob);
        % remove noise from out of bounds points
        y(locations) = y_prev(locations);

    end
    
    
% 	plot(x, y, '-', 'LineWidth', 2);
% 	grid on;
	% Interpolate y between those new y values.
	y = interp1(x(all_locations), y(all_locations), x);
	

end
% Now the ends may not be the original ends.  So let's rotate it down so that y(1) and y(end) are the original values.
% coefficients = polyfit([x(1), x(end)], [y(1), y(end)], 1);
% rampSignal = polyval(coefficients, x);
% y = y - rampSignal;
% y = y+start_pt-y(1);  % XX
% y = [linspace(x(1),y(1),20), y, linspace(y(end),x(end),20)];
y = (y-y(1))*(end_pt-start_pt)+start_pt; 

y = y'; 

