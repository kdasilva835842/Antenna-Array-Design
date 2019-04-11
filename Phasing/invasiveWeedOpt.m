coeffs = [0.3595, 0.2528, 0.2958, 0.3072, 0.3827, 0.5186, 0.5869, 0.6160, 0.7092, 0.7786, 0.8304, 0.8904, 0.9424, 0.9575, 1.0000, 0.9988, 0.9969, 1.0000, 0.9757, 0.9049, 0.8868, 0.8317, 0.7572, 0.7039, 0.6628, 0.5290, 0.4582, 0.4497, 0.3931, 0.2831, 0.1987, 0.3559];
coeffsInterp = interp(coeffs,2);

nPerRow = 32;
distribution = ones(nPerRow,nPerRow);


for r=1:nPerRow
    for c=1:nPerRow
        distribution(r,c) =  coeffs(1,c)*coeffs(1,r);
    end
end

nPerRow = 30;
x = linspace(-1.6,1.6,nPerRow);
sincFunc = sinc(x);
sincDist = 0.3*sincFunc +0.7026; %sincDist = 0.6*sincFunc +0.4026;
                                 %sincDist = 0.7*sincFunc +0.3026;

sincDistMat = ones(nPerRow,nPerRow);
for r=1:nPerRow
    for c=1:nPerRow
        sincDistMat(r,c) =  sincDist(1,c)*sincDist(1,r);
    end
end

singleColumn = reshape(sincDistMat,[1,nPerRow^2]);
phase = zeros(1,nPerRow^2);

importVector = [singleColumn;phase];

fileID = fopen('sincDistribution.txt','w');
fprintf(fileID,'%f %f\r\n',importVector); % '%6.2f %12.8f\r\n'
fclose(fileID);

figure(1)
surf(distribution)

figure(2)
plot(coeffs)

figure
plot(sincDist)

