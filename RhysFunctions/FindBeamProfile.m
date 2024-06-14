% % %inputs
directory = 'E:\labview-images';
FigNum = 10;


% % % Grab/plot data 
raw = BinaryImageData.loadImageSets('directory',directory,'rotation',90,'datatype','mono8');
%
figure(FigNum);clf
surf(raw.images(:,:,5))
title('DDS Channel 1')

figure(FigNum + 1);clf
surf(raw.images(:,:,6))
title('DDS Channel 2')

%%
[Beam1_cy,Beam1_cx] = find_beam_position(raw.images(:,:,5))
[Beam2_cy,Beam2_cx] = find_beam_position(raw.images(:,:,6))

% amp, pos, size, offset


%%
x = 1:numel(raw.images(:,1,5));
y = 1:numel(raw.images(1,:,5));
z = raw.images(:,:,3);

figure(FigNum + 2)
subplot(3,1,1)
plot(x)
subplot(3,1,2)
plot(y)
subplot(3,1,3)
surf(z)
