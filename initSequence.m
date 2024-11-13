function sq = initSequence
    sq = TimingSequence(32,24);
    
    %% Name digital channels
    sq.channels(1).setName('2DMOT','A0').setDefault(1);
    sq.channels(2).setName('3DMOT','A1').setDefault(1);
    sq.channels(3).setName('87 repump','A2').setDefault(1);
    sq.channels(4).setName('87 push','A3').setDefault(1);
    sq.channels(5).setName('87 imag','A4').setDefault(0);
    sq.channels(6).setName('85 repump','A5').setDefault(0);
    sq.channels(7).setName('85 push','A6').setDefault(0);
    sq.channels(8).setName('85 imag','A7').setDefault(0);
    sq.channels(9).setName('CD bit 0','B0').setDefault(0);
    sq.channels(10).setName('CD bit 1','B1').setDefault(0);
    sq.channels(11).setName('RF atten','B2').setDefault(0);
    sq.channels(12).setName('Repump shutter','B3').setDefault(1);
    sq.channels(13).setName('87 cam trig','B4').setDefault(0);
    sq.channels(14).setName('ND cam trig','B5').setDefault(0);
    sq.channels(15).setName('RedPower TTL','B6').setDefault(0);
    sq.channels(16).setName('Probe','B7').setDefault(0);
    sq.channels(17).setName('MW Switch','C0').setDefault(0);
    sq.channels(18).setName('H-Bridge Quad','C1').setDefault(1);
    sq.channels(19).setName('H-Bridge Helm','C2').setDefault(0);
    sq.channels(20).setName('MOT Bias','C3').setDefault(0);
    sq.channels(21).setName('Repump Switch','C4','Inverted').setDefault(0);
    sq.channels(22).setName('LG Shutter','C5').setDefault(0);
    sq.channels(23).setName('C6 - N/C','C6').setDefault(0);
    sq.channels(24).setName('RF Switch','C7').setDefault(0);
    sq.channels(25).setName('2D MOT Coils','D0','INVERTED').setConversionFunction(@(x) ~x).setDefault(1);
    sq.channels(26).setName('DDS TTL','D1').setDefault(0);
    sq.channels(27).setName('Feedback Laser TTL','D2').setDefault(0);
    sq.channels(28).setName('Raycus TTL','D3').setDefault(1);
    sq.channels(29).setName('DKC TTL','D4').setDefault(1);
    sq.channels(30).setName('vertical cam trig','D5').setDefault(0);
    sq.channels(31).setName('n/c','D6','DO NOT USE').setDefault(0);
    sq.channels(32).setName('n/c','D7','DO NOT USE').setDefault(0);
    
    %% Name analog channels
    sq.analog(1).setName('RF Frequency','AO/0')...
        .setConversionFunction(@RFtoV,'MHz')...
        .setBounds([0,20]).setDefault(20);
    sq.analog(2).setName('3DMOT Freq','AO/1')...
        .setConversionFunction(@(x) FtoV('trap',x),'MHz').setDefault(26);
    sq.analog(3).setName('87 repump freq','AO/2')...
        .setConversionFunction(@(x) FtoV('repump',x),'MHz').setDefault(0);
    sq.analog(4).setName('2DMOT amp','AO/3')...
        .setConversionFunction(@(x) x,'V').setDefault(0);
    sq.analog(5).setName('87 imag freq','AO/4')...
        .setConversionFunction(@(x) FtoV('image',x),'MHz').setDefault(0);
    sq.analog(6).setName('2DMOT Freq','AO/5','')...
        .setConversionFunction(@(x) FtoV('2d',x),'MHz').setDefault(14);
    sq.analog(7).setName('Raycus CW','AO/6','3.0V MAXIMUM')...
        .setDefault(0).setBounds([-0.1,3.5]);
    sq.analog(8).setName('Push Freq','AO/7')...
        .setConversionFunction(@(x) FtoV('push',x),'MHz').setDefault(15);
    sq.analog(9).setName('RedPower CW','B0/0')...
        .setDefault(-0.1).setBounds([-0.1,7]); %was setBounds([-0.1,2.6])
    sq.analog(10).setName('3DMOT amp','BO/1')...
        .setConversionFunction(@(x) TrapPtoV('trap',x),'').setDefault(1);
    sq.analog(11).setName('87 repump amp','BO/2')...
        .setConversionFunction(@(x) TrapPtoV('repump',x),'').setDefault(1);
    sq.analog(12).setName('MOT bias coil','BO/3')...
        .setDefault(0);
    sq.analog(13).setName('87 imag amp','BO/4')...
        .setConversionFunction(@(x) x,'V').setDefault(8);
    sq.analog(14).setName('Bias E/W','BO/5')...
        .setConversionFunction(@(x) x,'V').setDefault(0.4);
    sq.analog(15).setName('Variable Wave Plate','BO/6')...
        .setDefault(-3.4);
    sq.analog(16).setName('85 imag amp','BO/7')...
        .setConversionFunction(@(x) x,'V').setDefault(0);
    sq.analog(17).setName('CD3','CO/0')...
        .setConversionFunction(@(x) dBtoV('normal',x),'G/cm').setDefault(0);
    sq.analog(18).setName('CD2','CO/1')...
        .setConversionFunction(@(x) dBtoV('normal',x),'G/cm').setDefault(0);
    sq.analog(19).setName('CD1','CO/2')...
        .setConversionFunction(@(x) dBtoV('normal',x),'G/cm').setDefault(0);
    sq.analog(20).setName('CD0 Fast','CO/3')...
        .setConversionFunction(@(x) dBtoV('normal',x),'G/cm').setDefault(25);
    sq.analog(21).setName('CD Fine/Fast','CO/4')...
        .setConversionFunction(@(x) dBtoV('fine',x),'G/cm').setDefault(0);
    sq.analog(22).setName('Bias N/S','CO/5')...
        .setConversionFunction(@(x) x,'V').setDefault(1.4);
    sq.analog(23).setName('DPAOM Waveplate','CO/6')...
        .setDefault(0);
    sq.analog(24).setName('Bias U/D','CO/7')...
        .setConversionFunction(@(x) x,'V').setDefault(6.06);

    %% DDS channels
%     sq.dds(1).rfscale = 3;
%     sq.dds(2).rfscale = 3;
    calib_data = load('C:\Users\admin\Desktop\matlab-control\raman-aom-data.mat');
    sq.dds(1).calibrationData = calib_data.data(1);
    sq.dds(2).calibrationData = calib_data.data(2);
    sq.dds(1).powunits = DDSChannel.POW_UNITS_HEX;
    sq.dds(2).powunits = DDSChannel.POW_UNITS_HEX;
    sq.dds(1).setName('DDS 1').setDefault([110,0,0]);
    sq.dds(2).setName('DDS 2').setDefault([110,0,0]);
        
end