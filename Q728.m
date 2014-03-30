function average = SB_average(set)
    if ~isvector(set)
        error('Stupid human, the input parameter must be a vector')
    end
    average=sum(set)/length(set)
end

grid=zeros(27,0)
grid=[352.7 344.1 347.9 329.6 328.7 333.6 345.2 343.8 328.3 331.8 335.8 329.2 332.3 337.6 324.3 328.3 330 347.4 335.4 339.4 335 339 333.2 347 334.1 338.9 341.6]
mean=SB_average(grid)
% calculate a poitn estimate of the pull off force of all connections in the population
% and more
