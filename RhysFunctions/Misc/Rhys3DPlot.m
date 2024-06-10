function Rhys3DPlot(FigNum,Z,r)
% % % % Use this for plotting at the end
Z = r.data.N.mean;
param1 = r.data.param1;
param2 = r.data.param2;
% % % Sort the data
[Param1Sorted Index1] = sort(param1);
[Param2Sorted Index2] = sort(param2);
sorted_Z_rows = Z(Index1, :);
sorted_Z = sorted_Z_rows(:, Index2);


%%
[X Y] = meshgrid(param2,param1);

figure(FigNum);clf;hold on
for ii = 1:1:size(Z,2)
% for ii = [1 2 3 size(N,2)]
    scatter3(param1, param2(ii) * ones(size(param1)), Z(:,ii), 'LineWidth', 2);
    hold on
    plot3(Param1Sorted, Param2Sorted(ii) * ones(size(param1)), sorted_Z(:,ii),'k', 'LineWidth', 1);
end
view(3); % Set the default 3D view

% clear Z param1 param2 Param1Sorted Index1 Param2Sorted Index2 sorted_Z_rows Z sorted_Z X Y ii
end