signal = randn(10,1);

gaussFilter = gausswin(5)
gaussFilter = gaussFilter / sum(gaussFilter); % Normalize.

% Do the blur.
smoothedVector = conv(vector, gaussFilter,'same');

%find laplacian


figure;
hold on; 
plot(y1, 'r', 'linewidth', 3); 
plot(y2, 'b'); 
%plot(y3, 'g', 'linewidth', 3);