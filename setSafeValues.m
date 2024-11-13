function setSafeValues(sq)
    
for nn = 1:sq.numChannels
    sq.channels(nn).set(0);
end

sq.find('2D MOT coils').set(1);

% sq.find('DDS TTL').set(1);
