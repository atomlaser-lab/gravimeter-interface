# Table of Contents
  1. [Introduction](#introduction)
  2. [Initializing a sequence](#initializing-a-sequence)
  3. [Defining channel updates](#defining-channel-updates)
  4. [Building a multi-channel sequence](#building-a-multi-channel-sequence)
  5. [Compiling data](#compiling-data)
  6. [Uploading and running a single sequence](#uploading-and-running-a-single-sequence)
  7. [Visualizing sequences](#visualizing-sequences)
  8. [Executing multiple runs and parameter scans](#executing-multiple-runs-and-parameter-scans)
  9. [Additional features and tips and tricks](#additional-features-and-tips-and-tricks)

# Introduction

This set of MATLAB classes and functions is for creating and running sequences on the gravimeter and experiments with related control apparatuses.  It is meant to replace the run builder LabVIEW program, but it still needs the LabVIEW control program to handle writing data to the National Instruments hardware.  The idea behind this project is to replace defining sequences (or runs) using a graphical interface that needs lots of pointing and clicking with a set of MATLAB commands in text form.  This allows for easier inspection, modification, and automation of sequences.

The basic organization of this interface is as follows.  The classes

  - TimingSequence
  - TimingControllerChannel
  - DigitalChannel
  - AnalogChannel
  - DDSChannel

together are used to define a `TimingSequence` object which consists of an array of `TimingControllerChannels` (the channels), each of which is an instance of sub-classes `DigitalChannel`, `AnalogChannel`, or `DDSChannel`.  Each channel contains an array of values and the times at which to output those values.  These channels can be programmed independently.  The TimingSequence object has a compilation function which collates all the channel times and values and creates a set of arrays that are understood by the LabVIEW control program.  

The function `makeSequence` is the default program used for creating a specific `TimingSequence` object that runs a given experiment.  The idea is that `makeSequence` has a variable argument list, and the user of the program can edit the function so that different arguments do different things.  So `makeSequence` should be called as
```
sq = makeSequence(varargin);
```
where `sq` is the `TimingSequence` object and `varargin` is a list of input arguments.  `makeSequence` internally calls `initSequence` which defines the names of channels and their default values.

Finally, the `RemoteControl` class is used to manage communication with the LabVIEW control program.  It uses the TCP/IP protocol to send and receive data to the LabVIEW program.  The class can also be used to automate the running of multiple sequences in series, with input parameters to `makeSequence` varying in an arbitrary way, and then the subsequent analysis of that data.

# Initializing a sequence

If one wants to create a sequence from scratch, the first thing to do is to create a `TimingSequence` object `sq` by calling
```
sq = TimingSequence(numDigitalChannels,numAnalogChannels)
```
where the input arguments are the number of digital and analog channels, respectively.  For the gravimeter these are 32 and 24, respectively (as of November 2020).  The number of digital channels should not exceed 32.  This command creates the `TimingSequence` object `sq` with the appropriate channels, and each channel can be accessed using
```
sq.channels(idx); %For accessing any channel - the list starts with digital channels and then goes to analog channels
sq.digital(idx);  %For accessing the idx'th digital channel.  Call without any arguments to get an array of all digital channels
sq.analog(idx);   %For accessing the idx'th analog channel.  Call without any arguments to get an array of all analog channels
sq.dds(idx);      %For accessing the idx'th DDS channel.  Call as an array without the index to return all DDS channels
```
Remembering the index of each channel is inconvenient, so channels can be given names, ports, and descriptions.  For instance, one might want to call the 14th channel "Cam Trig".  This can be done using
```
sq.channels(14).setName('Cam Trig','B5','The camera trigger');
```
where the second string is a port number that is not used internally but may be useful for matching up the named channels to labels on breakout boards.  The last string is a description that can be useful for understanding what the channel does.  

Although the actual outputs of the channels are in volts, it is often more useful to have the channel outputs specified in terms of the physical things they change in the apparatus.  For instance, one analog voltage might change the detuning of the 3D MOT beams, while another sets the gradient of the MOT coils.  We define these values (detuning, gradients, relative powers, etc) as the "semantic" values, while the output voltages are the "physical" values.  To allow users to specify the sequence in terms of semantic values instead of physical values, channels have a conversion function that converts semantic values to physical values.  These can be specified using
```
sq.channels(14).setConversionFunction(@(x) x*1.6/8,'G/cm');
```
where the first argument is a function handle that does the conversion, and the second argument gives the unit.  Upper and lower bounds can be set using the `setBounds()` function, invoked as
```
sq.channels(14).setBounds([minBound,maxBound]);
```
where the bounds are in semantic units.  The `setBounds()` function can be appended after the `setName()` function to form a single line
```
sq.channels(37).setName('Some amplitude','AO/5','An amplifier amplitude').setBounds([minBound,maxBound]);
```
The function `sq = initSequence` should be used for defining all the channel names and default values (using the method setDefault()).  Once a channel is named it can be accessed using
```
ch = sq.find('Cam trig');
```
where the input argument for `find` is case-insensitive.  In the above command, `ch` is the `TimingControllerChannel` object corresponding to the channel with name 'Cam trig'.

# Defining channel updates

Sequences are created by telling channels what to do and when to do it -- these are called *updates*.  The fundamental method for this process is the `at(time,value)` method with input arguments `time` and `value`.  You can also use `on(time,value)` -- it's an alias.  The values for digital channels can only be 0 or 1, while for analog values they can be anything between the `bounds` property of the channel.  Times are always in seconds.  For instance, the command
```
sq.find('cam trig').at(6.05,1);
```
tells the camera trigger channel to ouput a value of 1 at 6.05 s from the start of the sequence.  If a value at that time already exists, it is overwritten.  The `at` command can also be used with an array of times and values
```
sq.find('cam trig').at([1,2,3,4],[1,0,1,0]);
```
or with an array of times and a function handle
```
sq.find('cam trig').at([1,2,3,4],@(x) mod(x,2)+1);
```
A series of `at` commands for a channel do not need to be in order
```
sq.find('cam trig').at(5,0);
sq.find('cam trig').at(1,1);
sq.find('cam trig').at(2.5,0);
```
When the sequence is compiled the updates are placed in chronological order.  If you want to force sorting of the updates, use the `TimingControllerChannel.sort()` method.

An important property associated with the `TimingControllerChannel` class is the `lastTime` property.  This value is set to the most recently set update time -- note that this does not mean that it is the latest time in a chronologically ordered set of times.  For instance, in the following set of commands the `lastTime` property is set as indicated.
```
sq.find('cam trig').at(5,0);    %lastTime = 5
sq.find('cam trig').at(1,1);    %lastTime = 1
sq.find('cam trig').at(2.5,0);  %lastTime = 2.5
```
If you call the sort method on a channel `ch` using `ch.sort()` then the `lastTime` property will be set to the latest time in the chronologically ordered set of times.

The `lastTime` property is used for certain semantically useful methods for defining channel updates.  There are four such methods: `anchor(time)`, `set(values)`, `before(times,values)`, and `after(times,values)`.  The `anchor(time)` method sets the value of `lastTime` to the input argument `time` without adding an update.  The method `set(value)` sets the value to `value` at the time corresponding to `lastTime`.  The methods `before(times,values)` and `after(times,values)` can be used to create updates that occur either before or after `lastTime`.  The methods `set(value)`, `before(times,values)`, and `after(times,values)` internally invoke `at(lastTime,value)`, `at(lastTime-times,values)` or `at(lastTime+times,values)` for the before and after methods, respectively, so they work with the same inputs as `at`, namely array inputs.  As a result, the following commands have the results indicated in the comments
```
sq.find('cam trig').at(0,0);                    %Sets the value to 0 a 0 s. lastTime is 0 s
sq.find('cam trig').at(3,1);                    %Sets the value to 1 at 3 s.  lastTime is 3 s.
sq.find('cam trig').after(50e-3,0);             %Sets the value to 0 at 3.05 s (50 ms after the last update).  lastTime is 3.05 s
sq.find('cam trig').anchor(10);                 %Sets lastTime to 10 s. No update is added
sq.find('cam trig').before(10e-3,1);            %Sets the value to 1 at 9.99 s (10 ms before 10 s). lastTime is set to 9.99 s.
sq.find('cam trig').after(50e-6,0);             %Sets the value to 0 at 9.99005 s (50 us after 9.99 s). lastTime is set to 9.99005 s.
sq.find('cam trig').at(15:20,@(x) mod(x,2)+1);  %Sets values to the output of the function at times 15 s through 20 s in 1 s steps.  lastTime is set to 20 s.
sq.find('cam trig').before(1e-3,1);             %Sets the value to 1 at 19.999 s (1 ms before 20 s). lastTime is set to 19.999 s.
```
Note that all of these commands return the channel object, which means that they can be chained together. The command
```
sq.find('cam trig').at(0,0).at(3,1).after(50e-3,0).anchor(10).before(10e-3,1).after(50e-6,0);
```
is equivalent to the first 6 lines of the last set of commands.  

## DDS Channels
An update to a DDS channel involves changing 3 parameters at the given time: the frequency, amplitude/power, and the phase.  The parent class `TimingControllerChannel` defines all methods for adding updates, and it can handle passing value arrays that are Nx3 where N is the number of updates.  Additionally, updates can be added using a multiple-argument syntax.  For instance, the following commands are equivalent:
```
t = linspace(-100e-3,100e-3,100)';  %Define a time vector. Note the transpose operator to make it a column vector!
f = 110 + t/200e-3;                 %Define a frequency ramp in MHz
p = exp(-t.^2/50e-6^2);             %Define a Gaussian pulse. The unit of power is a normalized optical power
ph = zeros(size(t));                %Define a phase array in radians
%The following commands are equivalent
sq.dds(1).after(t,[f,p,ph]);        %Pass the values as a single array
sq.dds(1).after(t,f,p,ph);          %Pass the values as multiple arguments
sq.dds(1).after(t,f,p,0);           %Any argument with only one element will be expanded to size(t)
```

Note that the units of the power are a normalized optical power that is converted internally into an RF power.  When initializing a `DDSChannel`, you should set the associated `rfscale` parameter as
```
sq.dds(1).setName('DDS 1').setDefault([110,0,0]);
sq.dds(1).rfscale = 2.38;
```
This sets the `rfscale` parameter to 2.38 W.  The conversion between normalized optical power and RF power is
```
rf = (asin((P).^0.25)*2/pi).^2*rfscale;
```

# Building a multi-channel sequence

Updates for each channel are added independently of the others, which makes it quite easy to create parallel sets of updates.  For instance, suppose one wants to create a basic imaging sequence with digital channels named 'imaging shutter ttl', 'imaging aom ttl', 'repump aom ttl', and 'cam trig'.  The imaging AOM and the camera trigger should start at the same time of 6 s and last for 30 us.  The repump AOM should turn on for 30 us before the camera trigger and last for the same duration.  The imaging shutter should be raised 2.5 ms before the imaging pulse starts to give it time to open and should close at the same time that the camera trigger goes low.  A set of commands that would work for this purpose is below:
```
sq.find('imaging aom ttl').at(6,1).after(30e-6,0);
sq.find('cam trig').at(6,1).after(30e-6,0);
sq.find('repump aom ttl').at(6,0).before(30e-6,1);
sq.find('imaging shutter ttl').anchor(6).before(2.5e-3,1).at(sq.find('cam trig').lastTime,1);
```
At the same time that these updates are occuring, one can set updates for other channels completely independently of the above set of commands.

For many experiments, however, it may be easier to consider the sequence as a set of sequential commands where the values for all relevant channels are set at each update time before a delay is added and new updates added at the delayed time.  A sequentially built process for the above behaviour might look like
```
sq.anchor(0);
sq.delay(6-2.5e-3);
sq.find('imaging shutter ttl').set(1);
sq.delay(2.5e-3-30e-6);
sq.find('repump aom ttl').set(1);
sq.delay(30e-6);
sq.find('repump aom ttl').set(0);
sq.find('imaging aom ttl').set(1);
sq.find('cam trig').set(1);
sq.delay(30e-6);
sq.find('imaging aom ttl').set(0);
sq.find('cam trig').set(0);
```
In this set of commands, I have started with the `TimingSequence.anchor(time)` command which sets the `lastTime` property for *every* channel to `time`.  Similarly, the method `TimingSequence.delay(delayTime)` (or `TimingSequence.wait(delayTime)`) advances the internal `TimingSequence.time` property by `delayTime` and sets the `lastTime` property of *every* channel to `TimingSequence.time`.  

One can also mix-and-match the two paradigms (parallel and sequential) to harness the power of both.  One could write the above sequence as
```
sq.find('cam trig').at(6,1);
sq.find('imaging aom ttl').at(sq.find('cam trig').last,1);
sq.find('imaging shutter ttl').anchor(sq.find('cam trig').last).before(2.5e-3,1);
sq.find('repump aom ttl').at(sq.find('cam trig').last,0).before(30e-6,1);
sq.delay(30e-6);
sq.find('imaging aom ttl').set(0);
sq.find('imaging shutter ttl').set(0);
sq.find('cam trig').set(0);
```

A slightly different function called `TimingSequence.waitFromLatest(delayTime)` first finds the latest update (chronogically speaking) and then sets that `TimingSequence.time` to that value plus `delayTime`; all channels' `lastTime` properties are set to the same value.  This is equivalent to `sq.anchor(sq.latest).wait(delayTime)`.  The difference between `wait` and `waitFromLatest` is important and best illustrated by an example with ramps. Consider the analog channel 'amp' and the sequence
```
sq.anchor(0);
sq.find('amp').at(0:10,@(x) x);
sq.wait(10);
sq.find('amp').set(0);
```
This sequence ramps 'amp' from 0 to 10 in 10 s with an update every second.  At the end of the ramp (10 s from its start), the value of 'amp' is set to 0.  Now consider the same sequence but using `waitFromLatest` (and assuming only 'amp' exists).
```
sq.anchor(0);
sq.find('amp').at(0:10,@(x) x);
sq.waitFromLatest(10);
sq.find('amp').set(0);
```
Since the `lastTime` value for 'amp' when the ramp is set is 10 s, the third command waits *another* 10 s after the end of the ramp before setting 'amp' to 0.  So this sequence takes a total of 20 s where the value of 'amp' for the last 10 s is 10.

## Using functions to simplify sequences

You can also create your own stand-alone functions that can simplify the creation of certain sequences that are fixed except for certain parameters.  A good example of this would be an imaging sequence, where you don't need to see all the details of the imaging sequence but you do want an easy way of changing the, for instance, imaging pulse duration.  You can create a function that modifies the whole `TimingSequence` object, or even just a couple of channels.  As a simple example, suppose we want to create a very simple absorption imaging sequence with two images.  We might define a function as
```
function makeImagingSequence(sq,expTime,camLoopTime)
%Creates two absorption images starting at the current sequence time

sq.anchor(sq.latest);           %Anchor the sequence time at now
%
% Creates the first absorption image
%
sq.find('imaging ttl').set(1);
sq.find('cam trig').set(1);
sq.delay(expTime);
sq.find('imaging ttl').set(0);
sq.find('cam trig').set(0);

%Wait for the camera to read out the image
sq.delay(camLoopTime);

%
% Creates the second absorption image
%
sq.find('imaging ttl').set(1);
sq.find('cam trig').set(1);
sq.delay(expTime);
sq.find('imaging ttl').set(0);
sq.find('cam trig').set(0);

end
```
Note that due to the way MATLAB deals with so-called 'handle' classes, you do not need to return the sequence `sq` as an output argument for the function.  MATLAB passes 'handle' objects by *reference*, which means that it passes to the function the memory location at which the object is stored: no copy of the object is made.  This is distinct from the way it passes all other variables, which is by *value* where the value of the variable is copied to a new location in memory which is discared when the function terminates.  What this means is that if you call `makeImagingSequence` in your `makeSequence` function as so
```
function makeSequence(varargin)

sq = initSequence;

% Bunch of sequence stuff with sq

% Call the imaging sequence function
% sq is modified by this function!
makeImagingSequence(sq,30e-6,50e-3);

end
```
Then the sequence object `sq` is modified by the function `makeImagingSequence`.

# Compiling data

Once a sequence is defined it has to be transformed into a form that can be sent to the LabVIEW control program.  Use the method
```
sq.compile;
```
to compile the data.  While the updates for each channel are programmed separately, the output from the hardware is such that any time there is an update from a channel, *all* channels must provide a value.  So compilation is the process by which the independent updates from each channel are combined into a set times and values such that values for inactive channels are held constant while active channels change.  For instance, suppose we have the analog channels '3D MOT Freq' and '3D MOT Amp' that are programmed as
```
sq.find('3D MOT Freq').at(0,6.8);
sq.find('3D MOT Amp').at([0,1,2,3],[8,7,6,5]);
```
The 'Freq' channel has only one update while the 'Amp' channel has four.  During compilation, the 'Freq' channel updates are populated such that it has updates at times 0, 1, 2, and 3 seconds with all values being 6.8.

Compiled data is stored in the `TimingSequence` property `TimingSequence.data` which is a structure with fields `t`, `d`, and `a`.  The `t` field is an Nx1 array of double precision values where N is the total number of compiled updates and represents all the update times.  The `d` field is an Nx1 array of unsigned 32-bit integers that represents the output of up to 32 digital channels at each update time.  The `a` field is an NxM array of double-precision values with M the number of analog channels and represents the analog values.

DDS channels do not run through the same hardware as the digital and analog channels, so their compilation process is a little different.  First off, the DDS will start only when it receives a falling edge trigger from the National Instruments box, so a property of the `TimingSequence` called `TimingSequence.ddsTrigDelay` has to be set to the time at which this edge occurs.  This property is used to shift the times of the `DDSChannel` objects from being referenced to when the whole sequence starts to being referenced to when the falling edge occurs.  DDS channel data is stored in the field `TimingSequence.data.dds`.

There are two additional properties in the `TimingSequence` class that interface with our LabVIEW control systems.  The first is the `camDelay` property, which is used to tell the LabVIEW camera VI when to start waiting for an image.  This is necessary because otherwise it will always time out, or one needs ridiculously long timeout values.  Typically, set this value to 1-3 s before you expect the first image to occur.  The second property is the `waitForImage` property, which tells the LabVIEW control interface to report the recorded image number before it tells MATLAB that the run is done.  This is used for saving copies of the run builder file and the sequence options, so that an image (or other data) file can be associated with the exact sequence and parameters used for its generation.  These flags are copied in the `TimingSeqence.data` structure.

The compiled data can be easily stored in a MATLAB data file and opened on a computer that does not have the interface classes and functions installed.  Additionally, the method `TimingSequence.loadCompiledData(data)` can convert a compiled data structure into a `TimingSequence`.

# Uploading and running a single sequence

The `RemoteControl` class is used for uploading sequences created with the `TimingController` class to the gravimeter control program and running it.  It can also be used to automate scanning through parameters for optimization in addition to on-line data analysis.  The `RemoteControl` class communicates with the LabVIEW control program using TCP/IP: the LabVIEW program listens for TCP/IP connections on port 6666 as the host, and the `RemoteControl` class connects to that host using the MATLAB `tcpip` class.  The `open()` and `stop()` can be used to connect to and disconnect from the LabVIEW host.  

The `RemoteControl` class has the property `sq` which is used for storing a sequence to be written to the device.  Additionally, there is a method `make()` that can be called as
```
r = RemoteControl;  %Create RemoteControl object
r.make(varargin);
```
which internally calls the function stored in `RemoteControl.makerCallback` as `r.sq = r.makerCallback(varargin)`.  If there is no callback specified in `makerCallback` then it throws an error.  If you wanted to specify a different function, say `myfunc`, then set it using
```
r.makerSequence = @myfunc;
```
The reason for specifying a separate `make()` method for the `RemoteControl` class is so that it can be chained together with the `upload` and `run` methods to enact a single-line make, compile, upload, and run command:
```
r.make(varargin).upload.run;
```
where the `upload` method called without an input argument uses the compiled data from the internal `sq` property.  The `run` method tells the LabVIEW program to execute the currently stored program.  Alternatively, given compiled data `data` you can upload that using `r.upload(data)`.  For those who are even lazier, there is a method called `urun` that calls the upload and run methods as above so the following two commands are the same
```
r.make(varargin).upload.run;
r.make(varargin).urun;
```
Both `run` and `urun` can accept an input argument that is a callback function which will be executed when the run is finished.  This is useful, for instance, when analyzing data for each run.  Suppose the callback function is an absorption image analysis program called `Abs_Analysis`.  Then running the command
```
r.make(varargin).upload.run(@Abs_Analysis);
```
(or `urun`) will cause `Abs_Analysis` to be called when the run completes.

If the program needs to be looped continuously while analyzing data, which might occur when aligning optics while monitoring atom number of temperature, then the user can use the `loop()` method:
```
r.make(varargin).upload.loop(@Abs_Analysis);
```
which will create, upload, and then run that sequence forever while calling the analysis function after every run.  Stop the infinite loop using `r.stop`.

Data destined for the National Instruments box is sent the LabVIEW control interface VI over TCP/IP.  Data for the DDS is converted into a series of commands for the MOGLabs ARF box and sent asynchronously.  This is necessary because MOGLabs designed a very stupid controller in the box itself which cannot handle more than one command at a time.  As a result, a set of commands cannot be sent together as a single block of text, which would cut done enormously on I/O time; instead, each command has to be sent separately.  To upload 1000 instructions (total) for two channels takes about 7 s.  If this uploading is done synchronously, in that it blocks the command line and prevents the sequence from running, then the cycle time of the experiment takes an extra 7 s.  Instead, the data is sent asynchronously and a message is printed on the command line when the upload is finished.  The user needs to ensure that the upload is complete before the DDS is triggered.

## Storing run builder files and options

If your sequence has set the `TimingSequence.waitForImage` property to `true`, then our current LabVIEW and MATLAB functions will save a copy of the run builder stored in `makerCallback` as well as the options used to create that sequence.  A copy of the run builder, with options embedded, is stored in the `RemoteControl.MAKER_STORAGE_DIRECTORY` as a `.m` file, and the sequence options are stored in the same directory as a MATLAB data file `.mat`.  This needs to wait for images to be saved so that it stores the files with the appropriate image number.

# Visualizing sequences

The update times and values of each channel can be accessed through the `TimingControllerChannel.times` and `TimingControllerChannel.values` properties, and these can be plotted in whatever way the user wishes.  To simplify matters, the `TimingControllerChannel.plot()` method has been included.  `plot()` takes a variable argument list in name/value pairs with valid names being:
  * 'offset': plot the channel values plus the given vertical offset.
  * 'finaltime': plot the channel values up to the given final time.
  * 'plotidx': plot the channel values corresponding to the given indices.

If the user wants to plot all the channels, they can use the `TimingSequence.plot()` function.  This has a single input argument which is the incremental vertical offset at which to plot each channel's values.

Finally, there is a GUI that can be used to display sequence data.  Given a `RemoteControl` object `r` in the base workspace that has a sequence in the field `sq`, so that `r.sq` is a TimingSequence object, the GUI can be started using
```
DisplayGUI(r);
```
Multiple channels can be plotted by selecting multiple channels on the right hand pane using either Shift-Click or Ctrl-Click.  The data in the GUI is automatically updated whenever a new sequence is made using the `r.make()` method.

# Executing multiple runs and parameter scans

The true power of this interface is that it allows for easy scanning of parameters in a sequence for optimization and data taking purposes.  This is accomplished by implementing a finite-state machine in the RemoteControl object with states 'initialize', 'set', and 'analyze'.  The actions undertaken in each of these states is governed by a callback function stored in the `RemoteControl.callback` property.  The callback function should have the basic structure:
```
function mycallback(r)  %r here is the RemoteControl object
  if r.isInit()
    %Initialize the set of runs. This stage is only called once.
    %Use it to specify the values of the parameter you want to change
    %and the number of runs to take
  elseif r.isSet()
    %This is called at the start of every run of the experiment.
    %Use it to create a new sequence to upload based on the current
    %parameter, and upload that sequence.  Do not use the r.run()
    %method here!
  elseif r.isAnalyze()
    %This is called after the LabVIEW control program is finished
    %Use it to analyze data generated by the experiment
    %You can store processed data in the experiment as fields in 
    %the property r.data
  end
end
```
To facilitate the automated collection of data there is the property `RemoteControl.data`, and different kinds of data can be stored as fields in this property.  For instance, one could store the number of atoms as `r.data.N` and the temperature as `r.data.T` with both `N` and `T` being vectors.

Keeping track of the total number of runs and the current run is handled by the property `RemoteControl.c` which is an instance of the `RolloverCounter` class.  `RolloverCounter` is a counter with multiple indices that increments using modular arithmetic -- it functions very much like a car odometer.  The first index has a particular maximum value, and when that index exceeds its maximum value it rolls back to its start value and the next index increments by 1.  Supposing that one creates a `RolloverCounter` object `c`, one can define a set of three counters with maximum values 3, 4, and 5, as
```
c = RolloverCounter([3,4,5]);   %Creates and initializes a RolloverCounter object with 3 indicies with maximum values 3, 4, and 5
fprintf(1,'Total number of runs is %d\n',c.total());  %Displays the total number of runs. c.total() returns the total number
while ~c.done
fprintf(1,'Index 1 %d/%d, Index 2 %d/%d, Index 3 %d/%d, Total counts %d/%d\n',c.i(1),c.final(1),c.i(2),c.final(2),c.i(3),c.final(3),c.now(),c.total());
c.increment();
end
```
The while loop shows how the counters increment.  The function `c.now()` (or `c.current()`) returns the current total index; i.e. the total number of increments that have occurred (minus 1).  The `done(idx)` method is true if the counter has reached the end of the current counter, specified by `idx`, or all counters if no `idx` is given.  You can reset the counter using the method `RolloverCounter.reset()`, and you can reuse a `RolloverCounter` object with different index ranges using the `RolloverCounter.setup()` method:
```
c.setup([10,2,5]);
```

Back to running multiple instances.  When the `RemoteControl` object is first created using `r = RemoteControl`, the value of `r.c` is set to a `RolloverCounter` object where the total number of runs is infinite.  After the successful completion of each run, which is signaled by the LabVIEW program sending a 'ready' word to MATLAB using TCP/IP, the value of `r.c.now()` is checked against `r.c.total()` and, if it is smaller, `r.c` is incremented using `r.c.increment()`.  If `r.c.now() == r.c.total()` then the set of runs is considered finished and the `r.stop()` method is called.  

A set of runs is started by using `r.start()`.  As long as `r.c.now() == 1` it will set the internal state of `r` to 'initialize' so when the callback function is executed it will execute the case corresponding to `r.isInit() == true`.  Use this case to define the parameters of interest and also the number of runs.  Note that `r.start()` **does not** reset `r.c`; this behaviour is so that if you can resume a sequence of runs in case of errors.  Use `r.reset` to reset the run counter and clear the `r.data` property.  From here, the state switches to 'set' and executes the callback case `r.isSet() == true`.  Use this case to create a sequence to upload based on the current parameter.  **Do not** use the `r.run()` method in the callback, as it is automatically called once the callback returns and is in the 'set' state.  When the LabVIEW control program indicates that it is done and ready for a new sequence, `r` moves to the 'analyze' state and executes the case `r.isAnalyze() == true`.  Use this to analyze the data generated by the sequence that just finished.  Pretty much anything can be placed into this section to do nearly any kind of analysis.  The data resulting from this analysis can then be stored as fields in the `r.data` property.

Let's consider an example of a very simple multiple run where we want to change the time-of-flight for the atoms to measure their temperature.  Let us suppose that we have set up our sequence creating function to be `makeSequence(tof)` where `tof` is the time of flight of the atoms.  We want to run through several times-of-flight and analyze the resulting absorption images.  Suppose that we have a function called `Abs_Analysis` that returns a structure with the *x* and *y* widths from the last absorption image.  A potential callback might look like
```
function MeasureTemperature(r)
  if r.isInit
    r.data.tof = 10e-3:2e-3:30e-3;    %Set the times of flights to scan through
    r.c.setup('var',r.data.tof);      %This is a method of setting the counter up just by passing the keyword 'var' and the parameters to loop through
  elseif r.isSet
    r.make(r.data.tof(r.c.now));       %Create the sequence
    r.upload;                                       %Upload the sequence
    %Print something to the command line so that we know how far along we are
    fprintf(1,'Run %d/%d, TOF: %.1f ms\n',r.c.now,r.c.total,r.data.tof(r.c.now)*1e3);
  elseif r.isAnalyze
    nn = r.c.now;                           %Make a shorter variable name
    img = Abs_Analysis;                     %Analyze the absorption image, return structure img
    r.data.xw(nn,1) = img.clouds.gaussWidth(1);             %Store the x width
    r.data.yw(nn,1) = img.clouds.gaussWidth(1);             %Store the y width

    %Plot the results
    figure(1);clf;
    plot(r.data.tof(1:nn),r.data.xw,'o');
    hold on;
    plot(r.data.tof(1:nn),r.data.yw,'o');
  end
end
```
This measurement can then be started by using the commands
```
r = RemoteControl;  %Assuming that there is no variable r yet
r.callback = @MeasureTemperature;
r.reset;r.start;
```
and the program will run through all the times of flight and plot the widths on the fly.  The stored data can be saved to disk, if desired, by using
```
data = r.data;
save('path/to/file.mat','data');
```
It is recommended that you save only simple data structures and not classes to disk because then loading those data files does not require the correct class definitions to be in your path.  Additionally, other programs such as Julia, Python, and Mathematica can easily read MATLAB .mat files when the variables that are saved are simple scalars, strings, arrays, and structures, but not when they are classes.

A more complicated example might be optimizing the temperature of a cloud after PGC by changing the cooling frequency.  To do this properly several times-of-flight need to be acquired for each cooling frequency, and the widths of the clouds need to be fitted to a model of ballistic expansion to get the temperature.  Suppose `makeSequence` is set up with input arguments `makeSequence(freq,tof)`.  A callback function that would accomplish this task might be
```
function TemperatureMeasurement(r)

  if r.isInit()
    %Initialize run
    r.data.tof = 14e-3:3e-3:29e-3; 
    r.data.freq = linspace(6.5:0.1:9);
    
    r.c.setup('var',r.data.tof,r.data.freq);
  elseif r.isSet()
    
    r.make(r.data.freq(r.c(1)),r.data.tof(r.c(2)));   %You can also use r.c.i(1) and r.c.i(2). The subscript indexing function has been redefined to allow the behaviour shown here
    r.upload;
    r.data.sq(r.c.now,1) = r.sq.data;  %This stores the sequence in the data property in case you need to go back and figure out what changed.
    fprintf(1,'Run %d/%d, Freq: %.3f V, TOF: %.1f ms\n',r.c.now,r.total,...
        r.data.freq(r.c(1)),r.data.tof(r.c(2))*1e3);
  elseif r.isAnalyze()
    pause(0.1);
    img = Abs_Analysis;
    if ~img.raw.status.ok()
        %
        % Checks for an error in loading the files (caused by a missed
        % image) and reruns the last sequence
        %
        r.c.decrement;
        return;
    end
    r.data.N(r.c(1),r.c(2)) = img.clouds.N;
    r.data.xw(r.c(1),r.c(2)) = img.clouds.gaussWidth(1);
    r.data.yw(r.c(1),r.c(2)) = img.clouds.gaussWidth(2);

    Ntof = numel(r.data.tof);
    if r.c.done(1)
      %After recording the desired times-of-flight, analyze data according to ballistic expansion model
      xfit = r.data.xw(:,r.c(2));
      yfit = r.data.yw(:,r.c(2));
      %Insert fitting routines for widths vs times of flight here to get xtemp and ytemp
      r.data.Tx(r.c(2),1) = xtemp;
      r.data.Ty(r.c(2),1) = ttemp;

      figure(1);clf;
      plot(r.data.freq(1:r.c(2),r.data.Tx,'o');
      hold on;
      plot(r.data.freq(1:r.c(2)),r.data.Ty,'sq');
    end
  end
end
```
This example has two counters, one for the time of flight and one for the frequency, and increments them in such a way to set the frequency, run through the times of flight, analyze the data when all the times of flight have been collected, and then go to the next frequency.  You will need to insert your own method for extracting the temperature from the widths vs times of flight.  As with the previous example, you can start it by using
```
r = RemoteControl;
r.callback = @TemperatureMeasurement;
r.reset;r.start;
```
and go get yourself a coffee as it automatically runs through 156 different sequences.  

This sequence also has an error checking section (the 'if' statement) which uses the imaging-analysis error RawImageData error checking to detect when an absorption image has not been taken properly.  If an error has occurred, the method uses the ''RolloverCounter.decrement()' function to go back one step, and then immediately returns from the callback function.  This re-runs the current iteration when an error occurs.  For long runs where lots of time might be wasted if an error occurs and is unchecked, it is suggested that you add error checking functionality.


# Additional features and tips and tricks

## Persistent data for multiple runs

Sometimes you want to keep certain data around between scans of parameters, such as MATLAB objects representing instruments.  The `RemoteControl.devices` property is a persistent property of the `RemoteControl` object that is *not* reset when `r.reset()` is called.

## Restarting multiple runs

Sometimes a parameter scan will fail because of a programming error and you will want to restart it from that point.  Suppose the error occurs somewhere in the analysis phase.  You can re-run the analysis using `r.analyze()`, and then manually increment the counter and start the process using
```
r.c.increment;
r.start;
```
If you need to go back and re-do a run, use the `r.c.decrement()` method instead.

Additionally, if you want to add more parameter values to your run while it is scanning, you can simply update the parameters and the counter.  Suppose that you are scanning over the parameters in the field `r.data.params`, and you want to add more points.  Use a set of commands such as
```
r.data.params = [r.data.params;(5:5:20)'];
r.c.final(1) = numel(r.data.params);
```
while the sequence is running to update the values.  If the sequence has already finished, use the above commands and then
```
r.c.increment;
r.start;
```

## Sequence options

It can be convenient to store sequence options that are commonly changed, such as the time of flight or imaging detuning, in a structure that is accepted as an input by the `makeSequence` function.  A class called `SequenceOptions` has been created for this purpose.  Currently, it has properties for `tof`, `detuning`, and more.  There is also a generic `params` field for parameters that you may wish to vary.  `SequenceOptions` is a handle class, so it gets passed by reference.  Use the `SequenceOptions` as follows
```
opt = SequenceOptions('tof',20e-3,'detuning',0,'load_time',10); %TOF = 20 ms, Detuning = 0 MHz, Loading time = 10 s
r.make(opt).urun; %Creates and uploads the sequence using opt for options
% You can also change options
opt.tof = 10e-3;        %Change tof to 10 ms
opt.set('detuning',8):  %Change detuning to 8 MHz
r.make(opt).urun;
% You can also change sequence parameters in the make() call
r.make(opt,'tof',30e-3).urun; %TOF is updated to be 30 ms, all other parameters are kept the same
```

Use the function `parse_maker_variable_argument_list(varargin{:})` to parse the variable argument list.

`SequenceOptions` is a subclass of the abstract superclass `SequenceOptionsAbstract`, which has common methods for setting, replacing, and printing values.  You can instantiate other subclasses of `SequenceOptionsAbstract` as properties in the main `SequenceOptions` class in order to better organize more complex option handling.  For instance, options related to non-destructive imaging are combined under the `FeedbackOptions` class and instantiated as `SequenceOptions.nd`.  

`SequenceOptions` are stored in the `RemoteControl.MAKER_STORAGE_LOCATION` directory as `.mat` files, labelled by the image number they generated.

### StageSequenceOptions

One particularly useful subclass is the `SequenceOptions.stage` property, which can be used to control which stages of the sequence are activated.  This is also a subclass of `SequenceOptionsAbstract` called `StageSequenceOptions`.  This class has the properties `mot`, `cmot`, `pgc`, `pump`, `mag`, `evap_map`, `dipoles`, and `evap_dipoles`, along with an `override` property.  Set these properties, and then in your particular make sequence property use the functions `use_mot()`, `use_cmot`, etc to determine if those sequences are activated.  Since the later stages invariably require the earlier stages, the output of `use_dipoles()` is only true if `mot`, `cmot`, `pgc`, `pump`, `mag`, `evap_mag`, and `dipoles` are true.  See the default function `makeSequence` for how to use the `StageSequenceOptions`.



