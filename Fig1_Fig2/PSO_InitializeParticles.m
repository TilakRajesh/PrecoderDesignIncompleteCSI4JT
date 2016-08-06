function [position, velocity] = PSO_InitializeParticles(numberOfParticles, numberOfVariables, xmax, xmin, alpha, deltaT)
% position = zeros(numberOfParticles,numberOfVariables);
% velocity = zeros(numberOfParticles,numberOfVariables);
% for i = 1:numberOfParticles
%     for j = 1:numberOfVariables
%         position(i,j) = xmin + rand*(xmax - xmin);
%         velocity(i,j) = (alpha/deltaT)*((-(xmax-xmin)/2) + rand*(xmax-xmin));
%     end
% end

position = xmin + rand(numberOfParticles,numberOfVariables).*(xmax - xmin);
velocity = (alpha/deltaT)*((-(xmax-xmin)/2) + rand(numberOfParticles,numberOfVariables).*(xmax-xmin));

end
