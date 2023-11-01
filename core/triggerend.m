%This function is used in the gravimeter control interface to trigger the
%end of r.make.urun(@triggerend,@Abs_Analysis_GUI); the idea is that this
%changes the run statue to cycleneded and the times wihtin the app that
%checks the status will detect the end of the run! 

function triggerend()
assignin('base', 'CycleEnded', true); %this is for the app to turn off the running!
end