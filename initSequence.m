function sq = initSequence
    sq = TimingSequence(32,24);
    
    %% Name digital channels
    sq.channels(1).setName('2DMOT','A0').setDefault(1);
    sq.channels(2).setName('3DMOT','A1').setDefault(1);
    sq.channels(3).setName('87 repump','A2').setDefault(1);
    sq.channels(4).setName('87 push','A3').setDefault(1);
    sq.channels(5).setName('87 imag','A4').setDefault(0);
    sq.channels(6).setName('85 repump','A5').setDefault(1);
    sq.channels(7).setName('85 push','A6').setDefault(0);
    sq.channels(8).setName('85 imag','A7').setDefault(0);
    sq.channels(9).setName('CD bit 0','B0').setDefault(0);
    sq.channels(10).setName('CD bit 1','B1').setDefault(1);
    sq.channels(11).setName('RF atten','B2').setDefault(0);
    sq.channels(12).setName('B3 - N/C','B3').setDefault(0);
    sq.channels(13).setName('87 cam trig','B4').setDefault(0);
    sq.channels(14).setName('B5 - N/C','B5').setDefault(0);
    sq.channels(15).setName('RedPower TTL','B6').setDefault(0);
    sq.channels(16).setName('DDS TTL','B7').setDefault(1);
    sq.channels(17).setName('C0 - N/C','C0').setDefault(0);
    sq.channels(18).setName('H-Bridge Quad','C1').setDefault(1);
    sq.channels(19).setName('H-Bridge Helm','C2').setDefault(0);
    sq.channels(20).setName('MOT Bias','C3').setDefault(1);
    sq.channels(21).setName('Repump Switch','C4','Inverted').setDefault(0);
    sq.channels(22).setName('C5 - N/C','C5').setDefault(0);
    sq.channels(23).setName('C6 - N/C','C6').setDefault(0);
    sq.channels(24).setName('RF Switch','C7').setDefault(0);
    sq.channels(25).setName('2D MOT Coils','D0').setDefault(1);
    sq.channels(26).setName('Probe','D1').setDefault(0);
    sq.channels(27).setName('Scope','D2').setDefault(0);
    sq.channels(28).setName('Control','D3').setDefault(0);
    sq.channels(29).setName('Stark','D4').setDefault(0);
%     sq.channels(30).setDefault(0);
    sq.channels(31).setDefault(1);
    sq.channels(32).setDefault(0);
    
    %% Name analog channels
    sq.analog(1).setName('RF Frequency','AO/0').setDefault(0).setConvert(@(x) RFtoV(x));
    sq.analog(2).setName('3DMOT Freq','AO/1').setDefault(7.1).setConvert(@(x) FtoV('trap',x));
    sq.analog(3).setName('87 repump freq','AO/2').setDefault(4.565).setConvert(@(x) FtoV('repump',x));
    sq.analog(4).setName('Keopsys MO','AO/3','3.9V MAXIMUM').setDefault(0).setBounds([0,3.9]);
    sq.analog(5).setName('87 imag freq','AO/4').setDefault(7.61).setConvert(@(x) FtoV('image',x));
    sq.analog(6).setName('85 Repump freq','AO/5').setDefault(4.64);
    sq.analog(7).setName('Keopsys FA','AO/6','3.5V MAXIMUM').setDefault(0).setBounds([0,3.5]).setConvert(@(x) DipolePtoV('keopsys',x));
    sq.analog(8).setName('85 imag freq','AO/7').setDefault(8.354);
    sq.analog(9).setName('RedPower CW','B0/0').setDefault(0).setBounds([0,7.5]).setConvert(@(x) DipolePtoV('redpower',x));
    sq.analog(10).setName('3DMOT amp','BO/1').setDefault(8).setConvert(@(x) TrapPtoV('trap',x));
    sq.analog(11).setName('87 repump amp','BO/2').setDefault(8).setConvert(@(x) TrapPtoV('repump',x));
    sq.analog(12).setName('MOT Bias Coil','BO/3').setDefault(2.9);
    sq.analog(13).setName('87 imag amp','BO/4').setDefault(8);
    sq.analog(14).setName('85 repump amp','BO/5').setDefault(8);
    sq.analog(15).setName('BO6 - N/C','BO/6').setDefault(0);
    sq.analog(16).setName('85 imag amp','BO/7').setDefault(3.6);
    sq.analog(17).setName('CD3','CO/0').setDefault(0).setConvert(@(x) dBtoV('normal',x));
    sq.analog(18).setName('CD2','CO/1').setDefault(0.8).setConvert(@(x) dBtoV('normal',x));
    sq.analog(19).setName('CD1','CO/2').setDefault(0.314).setConvert(@(x) dBtoV('normal',x));
    sq.analog(20).setName('CD0 Fast','CO/3').setDefault(0).setConvert(@(x) dBtoV('normal',x));
    sq.analog(21).setName('CD Fine/Fast','CO/4').setDefault(0).setConvert(@(x) dBtoV('fine',x));
    sq.analog(22).setName('CO/5 - N/C','CO/5').setDefault(0);
    sq.analog(23).setName('CO/6 - N/C','CO/6').setDefault(0);
    sq.analog(24).setName('CO/7 - N/C','CO/7').setDefault(0);

    %% Fake DDS channels
    sq.dds(1).rfscale = 3;
    sq.dds(2).rfscale = 3;
    sq.dds(1).setName('DDS 1').setDefault([110,0,0]);
    sq.dds(2).setName('DDS 2').setDefault([110,0,0]);
        
end