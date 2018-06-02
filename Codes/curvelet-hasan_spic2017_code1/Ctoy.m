function [y,C_array]=Ctoy(C)
y=0;
% Array to store sizes in C to use them if needed with inv.curvelet.:
C_array=zeros(1,4);
for i=1:length(C)
    for j=1:length(C{i})
        C_array=[C_array;i length(C{i}) size(C{i}{j})];
        y=[y reshape(C{i}{j},[1 numel(C{i}{j})])];
    end
end


y=y(2:end)';


end