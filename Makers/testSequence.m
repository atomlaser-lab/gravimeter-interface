function varargout = testSequence(varargin)

sq = initSequence;
trig_channel = sq.find('87 cam trig');
sq.analog(4).setName('Test');
sq.find('3Dmot').set(0);
sq.find('2Dmot').set(0);
sq.find('87 repump').set(0);
sq.find('CD bit 1').set(0);
sq.find('CD2').set(0);
sq.find('DKC TTL').set(1);
sq.find('87 imag amp').set(8);
sq.find('87 imag freq').set(FtoV('image',0));
sq.find('87 imag amp').set(7);
sq.camDelay = 2;
sq.delay(2);
sq.find('DKC power').set(0.3);
sq.delay(1);

time_at_drop = sq.time;
sq.find('dkc power').set(2);
sq.delay(20e-3);
sq.find('dkc power').set(0);

sq.anchor(time_at_drop + 35e-3);
trig_channel.set(1);
sq.delay(20e-6);
sq.find('87 imag').set(1).after(20e-6,0);
trig_channel.set(0);

% sq.delay(50e-3);
% trig_channel.set(1);
% sq.delay(20e-6);
% sq.find('87 imag').set(1).after(20e-6,0);
% trig_channel.set(0);

sq.delay(0.5);
trig_channel.set(1);
sq.delay(20e-6);
sq.find('87 imag').set(1).after(20e-6,0);
trig_channel.set(0);

sq.delay(0.5);
trig_channel.set(1);
sq.delay(20e-6);
trig_channel.set(0);

varargout{1} = sq;
