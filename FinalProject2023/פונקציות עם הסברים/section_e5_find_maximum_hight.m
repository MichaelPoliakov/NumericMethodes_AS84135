function section_e5_find_maximum_hight(y)
% find the max height of the rocket
% INPUT: y- the state vector we derived
max = y(1,3);
for z = y(:,3)'
    if z>max
        max = z;
    end
end
fprintf('Maximum hight rocket reaches:\n % .4f[km] = % .4f[feet]. \n', ...
    max*1e-3,convlength(max,'m','ft'))
end