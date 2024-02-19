files = get_correct_image_files('E:\SpatialFringes_data\tides-take2',2,6943 - 500,1);

% MZ_data = zeros(6943,2);
MZ2_data = zeros(size(files,2),2);
for nn = 1:size(files,2)
    nn
    img = Abs_Analysis_V2('files',files(:,nn));
%     pause(10)
    MZ2_data(nn,:) = img.get('Nsum');
end
%% Old method
% MZ_data = zeros(6943,2);
% for i = 500:6942
%     img = Abs_Analysis_V2('last',6943-i);
%     MZ_data(i+1,:)=[img.clouds(1).Nsum, img.clouds(2).Nsum];
%     i
% end