%% Analyzing Neural Time Series Data :: a tutorial in Matlab
% This "tutorial" goes through the 'Analyzing Neural Time Series' (ANTSD)
% book, chapter by chapter, expanding on some of the examples provided in
% the code by M.X. Cohen and also doing the exercises at the ends of some
% of the chapters. This tutorial assumes you are already familiar with
% Matlab (e.g., creating and manipulating arrays, structures, and so on).
%
% I recommend that you use this 'tutorial' alongside the book and M.X.
% Cohen's own chapter-by-chapter scripts recreating figures in the book.
% Cohen's scripts will also often jump ahead to later topics to demonstrate
% various things early on (e.g., Morlet wavelets in Chapter 2's scripts,
% when wavelets aren't discussed until about 10 chapter later). This is to
% demonstrate the multidimensionality of discrete time-varying signals
% (like EEG) right away. But it can be confusing if you're going through
% scripts chapter-by-chapter (you kinda just have to take his word for
% what's happening). I'll try to avoid doing that, but you can still turn
% to his scripts to see how all plots are made in each chapter. So, again:
% this code is to complement his work, not to substitute his work.
%
% @hbt7


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                        %
%%                               CHAPTER 1.                               %
%%                                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Cognitive electrophysiology is the study of how cognitive functions are
% underpinned by electrical activity of populations of neurons. Our goal is
% to understand the conceptual, mathematical, and implementational (e.g.,
% in Matlab) bases of time-, time-frequency, and synchronization-based
% analyses of time series data such as electroencephalography (EEG), 
% magnetoencephalography (MEG), or local field potential (LFP).
%           As long as the data contain a discretely sampled time-varying 
% signal, the following tutorial could prove relevant in analyzing those 
% that. Importantly, this requires an understanding of the math and logic 
% behind the various steps involved in analyzing such data and acquiring 
% that understanding is the goal of this book.

% Throughout the tutorial, we'll be working with an example EEG dataset.
% Let's take a quick look at it, so that we know what typical EEG data look
% like.

% sampleEEGdata :: 
% This dataset is one condition, from one subject, from the study Cohen & 
% Ridderinkhof (2013). The data have already been cleaned (manual trial
% rejection and independent components analysis).

%% 1.1 Change directory and load example data
% We can load the sampleEEGdata as follows in Matlab. First, change the
% directory (using the command 'cd') to go to the folder that contains the
% sampleEEGdata (which will be the folder downloaded from Cohen's github.
ANTSD_dir = '/Users/hamid/Desktop/Research/+Github/MikeXCohen/AnalyzingNeuralTimeSeries-main';
cd(ANTSD_dir)

load sampleEEGdata

% Double click on the 'EEG' structure in the Workspace in Matlab. You'll
% see that it's a 1x1 structure with 41 'fields'. For instance, there is an
% 'nbchan' field with the value 64, a 'trials' field with the value 99.
% This means that we have an EEG dataset of 99 trials, recorded with 64
% channels. The 'srate', i.e. sample rate, is 256. This is in Hertz (Hz),
% meaning that they recorded a 'snapshot' of the electrical activity 256
% times per second. This is confirmed when we look at the field 'times'. If
% you click on 'times' in the 'EEG' structure, you'll be shown a bunch of
% number starting from -1000, then -996.0938, and so on. This vector of
% numbers represents time points around the trial onsets, where 0 is the
% moment at which, e.g., a stimulus appeared on screen. So, -1000 means
% 1000 milliseconds before the stimulus onset. Remember that the srate is
% 256 Hz, so 1/256 = .0039, meaning that they recorded data every .0039
% seconds, which is once every (.0039 * 1000) 3.9 milliseconds. Indeed,
% when you run (1/EEG.srate)*1000, Matlab returns 3.9062, which is the
% difference between each of the time points in EEG.times (i.e., -1000 to
% -996.0938 is a difference of 3.9062). And there are 640 'pnts', meaning
% 640 samples recorded per trial and, therefore, that there are (640 *
% 3.9062) 2.5 seconds represented in each trial, which in turn means that
% the trial goes from -1 to +1.5 seconds (well, technically, 1500-3.9062,
% which is 1.4961e+03, which is indeed the last entry of EEG.times).

%% 1.2 Plot a single trial's data
% We can plot the data from a single channel, for a single trial. Looking
% at EEG.data, it's a 64 x 640 x 99 field. Thus, the first dimension is the
% channel, the second refers to the time points in the trial, and the
% last to trials. So, let's plot the entirety of trial 50 on channel 1 for 
% this one participant's data. We'll also plot an x=0 line, that tells us 
% where a trial started/stimulus appeared.
plot(EEG.times,EEG.data(1,:,50)); xline(0,"LineWidth",2);

% Basically, we're seeing 640 data points, all connected by a line,
% representing the recorded electrical data of some population of neurons.
% Note, we don't have continuous data, but Matlab shows us these points
% as though it were a continuous line. They are, instead, discrete data,
% sampled over time, and - as we can see - varying over time. Thus, it is a 
% 'discretely sampled time-varying signal' - the kind of signal we will
% learn to analyze.

% The question now is: is there something in data like these (i.e.,
% discretely sampled time-varying) that reflects the structure of our
% experiment? In other words, if we have two conditions, A and B, do the
% data look different between those conditions? How can we tell?

% End of Chapter 1
clear all; close all; clc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                        %
%%                               CHAPTER 2.                               %
%%                                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Neuroimaging techniques have to make trade-offs; some techniques are
% higher in spatial resolution (telling us /where/ a process is taking
% place), other in temporal resolution (telling us /when/ a process is
% taking place). EEG, which is going to be the primary technique through
% which we'll be learning about ANTSD, falls in the latter group. EEG has
% high temporal resolution, measures the neural activity of neurons at the
% population level directly (as in, we're measuring voltage fluctuations of
% neurons, rather than oxygenation of blood such as with fMRI), and EEG is
% multidimensional (e.g., time, space, frequency, power, phase, each
% providing different insights). There are also bad sides to using EEG:
% it's not suited for localizing the activity or for processes that slow
% and unfold with a very uncertain and variable time course.

% The unit of EEG measurement is volts (typically, microvolts, or µV).
% These µV values are relative to the reference, thus should not be
% interpreted in an absolute sense. Indeed, the values will change
% depending on your processing pipeline (e.g., based on your choice of 
% reference) and across participants (e.g. based on skull shape and
% thickness, whether they washed their hair that morning, the exact 
% placement of the electrodes/EEG cap relative to the underlying cortical
% morphology). That's why many analyses involve scale transformations to be
% able to compare results across electrodes, participants, recording
% equiptment, and publications.

% One particularly common analysis approach is to compute event-related
% potentials (ERPs). They are easy to compute, require few assumptions or
% parameters, offer high temporal precision and accuracy, can be
% contexualized in a decades-long literature, and - because of their
% simplicity to compute - allow for straightforward data quality checks.
% However, there are also limitations: null results are difficult to
% interpret and they don't easily allow for linking to physiological
% mechanisms.

%% 2.1 ERPs

load sampleEEGdata

% Recall that all trials in this sample are from the same condition. An ERP
% is simply the average of all the time series data for that condition on a
% given channel. To do so, we need to first extract the same window of data 
% around each stimulus onset (i.e., the moments you've defined as time=0).
% Fortunately, that's already been done for us. So, we can just go ahead 
% and average the data.

channel = 1; % Which of the 64 channels we'll plot from
n.trials = 6; % Number of trials to plot

% Plot a few of the individual trials
figure
title('The first n trials')
for t = 1:n.trials
    subplot(n.trials,1,t)
    plot(EEG.times,EEG.data(channel,:,t))
    xline(0,"LineWidth",2);
    ylabel('µV') 
end
xlabel('Time (ms)') 

% These are the first 6 trials in our data on channel 1 (unless you changed
% the above parameters). Now we can average them together and see what
% happens when we do so.

% To make the code more accessible, we can first extract the data for a
% single channel, as follows:
data = squeeze(EEG.data(channel,:,:))';
% Squeeze() drops the unused 3rd dimension, but the trials then end up in
% columns with data points in rows. We flip this back by using Matlab's
% transpose function which is the ' at the end. Chan_data now contains 99 
% rows (trials) and 640 columns (time points).

% Matlab performs averaging column-wise (unless there's only one row, or 
% unless you specify a different dimension explicitly). So, since all our 
% time points are along the columns, we are now averaging all 99 values 
% (one for each trial) that were recorded at time point -1000 ms. Then, we 
% average all 99 values representing -996.0938 ms, and so on. And.. that's 
% it! That's our ERP, the electrical potential related to our event. 
% Easy-peasy, data squeezy.
ERP = mean(data);
figure
plot(EEG.times,ERP); xline(0,"LineWidth",2);
title('The ERP'); xlabel('Time (ms)'); ylabel('Signal (µV)') 

% The ERP shows.. fluctuations, again. But, comparing it to our previously
% plotted single-trial data, the ERP seems to have fewer fast squiggles and
% more slower moving fluctuations. This should be intuitive: the more
% trials we average, the more we are smoothing out the fast fluctuations
% and the more we are emphasizing slower changes that are structurally
% present in the data. Averaging, essentially, works like a low-pass filter
% (i.e., rolling off high frequencies and letting the slower frequencies
% "pass through" the filter). Naturally, if there are no fluctuations (fast
% or slow) that are present in the data structurally (meaning, that they
% result from the stimulus and unfold more or less the same way on each
% trial), then averaging trials together will just result in a more and
% more flat line. This is shown in Figure 2.1 in the book, where averaging
% simulated data results in a flat line. However, recall that EEG data is
% multidimensional. So, just because something structural doesn't show up
% on the ERP, doesn't imply there's no structure in /all/ dimensions of the
% EEG data. Figure 2.1 shows that there can still be consistent alpha power
% effects, despite no clear ERP. Figure 2.2 shows another example of
% averaged trials with no structure, but a time-frequency plot that shows
% structure in frequencies over time.

% Indeed, another thing one can do is a time-frequency analysis. As the
% name suggests, this looks at frequency as a function of time. These
% analyses (relative to ERPs) are easier to link to physiology (of neural 
% oscillations), allow for better cross-species comparisons, and they are 
% perhaps more likely to uncover statistical structure in the data. However,
% they also decrease temporal precision and this kind of analysis is more
% complicated than a simple ERP. We'll learn more about time-frequency
% analyses in subsequent chapters.

% End of Chapter 2
clear all; close all; clc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                        %
%%                               CHAPTER 3.                               %
%%                                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Discrete time-varying signals like EEG contain rhythms reflecting neural
% oscillations, which - in turn - reflect fluctuations in the excitability
% of neural populations. These rhythms are described by three features:
% frequency, power, and phase. These features can be extracted with
% signal-processing techniques that we'll be learning about in this book.

% Let's quickly remind ourselves of what each of those are.
% Add the path that contains ANTSD_animatedwaves.m then run the animation.
addpath '/Users/hamid/Desktop/Research/+Tutorials/ANTSD'
ANTSD_animatedwaves(1)

% You can change the number passed to the function to change the number of
% rotations around the circle. Note that the angle, theta, is in radians.

% In the top left corner, we see a dot moving around the unit circle. This
% dot thereby creates an angle between the positive part of the x-axis, the
% origin of the system (0,0), and itself. This angle is referred to as 
% theta. So, when the dot is all the way at the top, theta = 90 degrees.
% All the way on the left, it's theta = 180 degrees, at the bottom theta =
% 270 degrees. When it's all the way at the rightest part of the circle,
% it's theta = 360 degrees (or simply 0 degrees, because it's all the way
% back to where it started. I hope your basic trigonometry knowledge is now
% coming back!

% In the top right, we see the sine of each angle. Look at the horizontal
% red lines in the top left figure and the top right one. The sine of a
% given angle theta is simply how far up the dot is on the y-axis.

% Similarly, the cosine of theta (pictured in the bottom left), shows us
% where the dot is along the x-axis.

% The bottom right shows both the sine and cosine waves plotted alongside
% each other as the dot moves along the circle. Essentially, they are the
% same wave, but shifted along the x-axis.

% If you're having trouble with the animation, try this video:
% https://www.youtube.com/watch?v=z82I6u4DFTo


% Note that, because we're dealing with the unit circle (a circle with
% radius one), sin and cos fluctuate between -1 and +1. That's the
% amplitude of the wave, 1 (the maximal sway of the wave away from its
% central point, in this case 0). But it doesn't need to be that way.

figure; subplot(711)  % Initialize a figure
x=linspace(0,1,1000); % A thousand linearly spaced points between 0 and 1
plot(x,   sin(2*pi*x))% Plot the sine of all those 1000 points
subplot(712)
plot(x, 2*sin(2*pi*x)) % Plot the sine, but time 2.
% Now we've doubled the amplitude. Look at the y-axis! But we can do more!

subplot(713)
plot(x,   sin(2*pi*x) + 1) % This shifts the entire wave up on the y-axis

subplot(714)
plot(x,   sin(2*pi*x - 1)) % This shifts the entire wave to the RIGHT
subplot(715)
plot(x,   sin(2*pi*x + 1)) % This shifts the entire wave to the LEFT
% These shifts may feel unintuitive at first. But adding something to the
% value that gets given to the sin() function, means we are reaching a
% given output of sin() faster, so the entire wave moves to the LEFT,
% whereas subtracting something means we're reaching a given output later,
% so everything shifts towards the right.
 
% For now, think of where we are on the curve as being called the phase. 
% And so, these shifts to the left and right of the curve are changes in 
% its phase relative to x = 0 (it's a bit more complicated, but we'll get 
% to that later).

% We can also multiply, in addition to addition and subtraction. In the
% other plots, the wave goes up and down once. But if we multiply the input
% by - say - 3, here's what happens
subplot(716)
plot(x,   sin(2*pi*x * 3))
% Now the wave does its 'once-up-once-down' 3 times before reaching the end
% of the plot. Each 'once-up-once-down' is called a period, which is why
% these functions are also called periodic functions. The fact that the
% wave has 3 periods means that its frequency = 3. So, the previous plots
% all had a frequency of 1.

% If we now take the x-axis to refer to time, in seconds, that means that
% that the first 5 plots have a frequency of 1 Hertz (Hz) and the last one
% 3 Hz. Hertz is the unit that refers to 'how often something happens in a
% given second'.

subplot(717)
plot(x,   sin(2*pi*x) + sin(2*pi*x * 3) + .7 * sin(2*pi*x * 36 - .5) + 1)
% Finally, waves can also be added together to form a new curve. When we do
% so, you see a bunch of squiggles that start to look a bit like the data
% that we get from an EEG device. The big question is, if we start with the
% actual data itself, can we somehow figure out what the original waves
% were that were all added together to produce the final data? In other
% words, in the wave above, can we figure out that one of the waves had a
% frequency of 3 Hz and the other 36 Hz? What about that the 36 Hz wave had
% an amplitude of 0.7 and was shifted in time? Turns out, yes we can!

% The thing to remember is that our discrete time-varying data is a single
% curve (e.g., the data from a single electrode) that represents a 
% summation of multiple waves, with different amplitudes, frequencies, and 
% phases. And, using signal processing techniques, we can actually
% reconstruct what the original waves were that summed to the final single
% curve we see in our data.

% This series of waves - of varying frequencies, amplitudes, and phases -
% which constitute a given signal is known as a 'Fourier series'.

% End of Chapter 3
clear all; close all; clc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                        %
%%                               CHAPTER 4.                               %
%%                                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If you are indeed already familiar with programming in Matlab, you can
% probably skip this chapter. If not, try the exercises at the end of the
% chapter and see how far you get.

%% Exercises
% Note that there are many ways to solve any of the exercises. My solutions
% are aimed at making the solution easy to understand and not necessarily
% 'elegant and/or efficient'.

addpath '/Users/hamid/Desktop/Research/+Tutorials/ANTSD'

% 4.7.1
n_rows = 4; n_cols = 8;
ANTSD_randmatrix(n_rows,n_cols)

% 4.7.2
% Show the image of Amsterdam
img = imread('amsterdam.bmp');
imshow(img)

% Draw the red line
line([360 365], [300 80], 'Color', 'red', 'LineWidth', 8)

% Plot a 'star' (his exercise solutions show something simple like this,
% perhaps something drawn with multiple lines using the line() function -
% don't overthink it)
text(340,520,"*",'Color','magenta','FontSize',80)

% Find the max value for R, G, and B. If there are multiple pixels with
% that max value, select a random pixel.
R=img(:,:,1); maxR=max(R); rand_maxR=maxR(randi(length(maxR)));
G=img(:,:,2); maxG=max(G); rand_maxG=maxG(randi(length(maxG)));
B=img(:,:,3); maxB=max(B); rand_maxB=maxB(randi(length(maxB)));
% Then, draw 3 circles with each of those R, G, B values.
theta = 0 : 0.01 : 2*pi;
radius = 20;
x1 = radius * cos(theta) + 100; y1 = radius * sin(theta) + 100;
x2 = radius * cos(theta) + 300; y2 = radius * sin(theta) + 500;
x3 = radius * cos(theta) + 600; y3 = radius * sin(theta) + 600;
hold on; plot(x1, y1, 'Color', [rand_maxR 0 0], 'LineWidth', 3);
hold on; plot(x2, y2, 'Color', [0 rand_maxG 0], 'LineWidth', 3);
hold on; plot(x3, y3, 'Color', [0 0 rand_maxB], 'LineWidth', 3);

% 4.7.3
matrix = ANTSD_randmatrix(32,3);
writematrix(matrix, "output_exercise473.txt", 'delimiter', '\t')

% End of Chapter 4
clear all; close all; clc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                               CHAPTER 5.                                %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% EEG reflects (mainly) the summation of excitatory and inhibitory
% postsynaptic potentials at the dendrites of ensembles of neurons with
% parallel geometric orientation. In particular, the signal is dominated by
% 10,000 - 50,000 neurons in the superficial cortical layers. Therefore,
% EEG does not and cannot measure most neural events. Then again, this is
% the case for all neuroimaging techniques, so no need to feel
% disheartened. But for scalp EEG specifically, it's difficult to measure
% subcortical areas (like thalamus, hippocampus etc.). Very slow
% fluctuations (<.01 Hz) are also difficult to pick up, since most modern
% EEG systems have built-in high-pass filters that attenuate very slow
% fluctuations (because those can cause amplifier saturation). If you're
% looking to measure below 1 Hz, some DC-coupled amplifiers may work. Alas,
% very fast fluctuations (>100 Hz) are also hard to measure, since they
% tend to have less power (meaning they are tougher to distinguish from
% noise).

% Figure 5.2 shows us that ERPs show up when structural fluctuations are
% aligned in time (e.g., elicitated by a stimulus and then aligned to that
% stimulus onset during data analysis) and aligned in their phase. This
% should be intuitive, even if things are aligned in time, if the data are
% out of phase across trials, averaging them all together will result in a
% flat(ish) line. Therefore, because ERPs can indeed be seen in EEG data,
% this must mean that there is a neural mechanism that allows for them to
% appear. Some proposals are that stimuli produce a phase reset (thereby
% aligning the phases in the moments following a stimulus, meaning that an
% ERP can appear). You can go back to Chapter 3's code, make various in and
% out of phase simulated data and average them to see for yourself how
% misaligned phases will result in a signal disappearing across averaging.
% Various other proposals about ERPs suggest that an elicited signal is
% additive, meaning that there is an ongoing background noise and a
% stimulus produces a signal that 'adds on top' of the noise; across many
% averaged trials, the noise gets averaged out and the signal remains.
% Another possiblity is amplitude assymmetry/baseline shift, where changes
% in power - elicited by a stimulus - produce assymetry in the ongoing EEG,
% which results in ERPs when averaged across trials. This debate is still
% ongoing, so it remains unclear where the ERP is exactly coming from.
%          What we can hope for, at least, is that the electrical signals
% are indeed somehow involved in cognition. There are several lines of
% evidence that reassure us that this is the case. But then, even if they
% weren't involved causally - and that things such as ERPs were mere
% epiphenomena - they may still be informative.

% End of Chapter 5


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                               CHAPTER 6.                                %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% It should be clear by now what you need, in terms of an EEG lab and an
% experimental design, to run a good EEG experiment. More resources are
% available, I can recommend Luck's (2005) ERP book myself. But, in short,
% you need high quality data and a well-designed experiment
% (unsurprisingly). You need event markers in the data so that you can line
% up your trials. Intertrial intervals can be jittered or not, depending on
% your needs. More trials per experimental condition are better in theory,
% but in practice your participants will get tired, meaning that
% eventually, an added block of trials with sleepy performance may be
% adding noise to your data instead of signal. More electrodes is better in
% theory, but in practice it adds time to set up, time to clean up, and
% storage space & analysis time to produce results. But above all, you need
% high quality data: there's NOTHING you can do in your data preprocessing
% that can save low quality data. High quality data with a poor design
% might, potentially, maybe, perhaps produce an interesting and robust
% result. But low quality data (riddled with noise from bad electrodes
% recording data from a sleepy participant) won't be saved by a perfect
% experimental design.

% End of Chapter 6


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                               CHAPTER 7.                                %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Okay, let's say you managed to get your high quality data - the bare
% minimum you needed. Then what? Analyses? Well, hold on! First, you gotta
% preprocess your data. The main thing you need to do is.. keep a record
% of everything you do. Which trials did you reject? Which electrodes were
% interpolated? Which independent components were removed from the data?
% All of it. Everything. That way, you can always go back to the raw data
% and recreate all steps of the preprocessing so that you can reproduce
% your results, but also because at each stage of the prepro, you might end
% up making choices that introduce problems/biases/errors into the data
% that you will want to reverse and subsequently avoid.

load sampleEEGdata

% By now it should be obvious that one thing we're gonna need to do is
% create so-called epochs (when we have task-based experiments, we don't
% necessarily need epochs for resting-state data). Epoching data
% essentially means that we're taking out the those parts of the data
% around trials. Take another look at the sample EEG data. Before prepro,
% each of the 64 channels (i.e. one for each electrode) had a continuous
% stream of data recorded from the moment the experimenter hit record until
% the moment that they stopped the recorded. An entire hour, perhaps two
% could have gone by. Using event markers, the entire recording is cut into
% trials (in this case, 99 trials). If you recall, these windows around the
% event markers were 1 second before up to 1.5 seconds after each marker.
% With your own data, you'll need to do that too. This also means that
% there needs to be sufficient time available before and after each marker
% for you to do so. If your trials are too fast, you won't be able to
% window out sufficiently lengthy trials. The time before the stimulus
% onset will be useful to do e.g. baselining. But all this is the minimum
% for ERPs. If you wanna do time-frequency analyses (i.e., where you
% analyze the various frequency bands over time), you'll need even longer
% windows of time to have sufficient data available to resolve information
% about the various frequencies. If you don't have long enough epochs, your
% results will be contaminated with edge artifacts like the ones seen in
% Figure 7.2. As a rule of thumb, you should have at least 3 full periods
% of data of the lowest frequency you're interested. So, if you're
% interested in 2 Hz activity, that requires at least 1500 ms per epoch (2 
% Hz means that a full period happens twice per second, so once per 500 ms,
% meaning 3 times per [3 * 500 ms =] 1500 ms). This extra buffer is needed
% with e.g. Morlet wavelet convolution or a Hilbert filter, all of which
% we'll be learning about in later chapters. For some approaches, this
% buffer isn't as important (we'll also learn about those), but allowing
% for (and making use of) longer epochs is always helpful.

% In addition to the above, you will also need to make decisions about
% filters to apply, trials to reject, a reference to pick, and so on. None
% of this is easy, because there isn't a right answer to any of it.

% But one thing is always, always, always true: prepro can turn good data
% into a great data, but it can't turn low quality data into great data
% (and barely into okay-ish data). You can always go back and make
% different prepro choices about filters and references, but you can't go
% back and make your low quality data better. So, always, always, always
% ensure that your raw data is of high quality. That should be your main
% priority during data collection.

% End of Chapter 7
clear all; close all; clc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                               CHAPTER 8.                                %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Even high quality data will still need preprocessing. All data -
% including those of high quality - will contain artifacts that you will
% likely want to remove before running analyses. These artifacts include
% eye blinks, muscle movements, transient amplifier saturations, line
% noise, and more (not to mention various forms of 'cognitive' artifacts).
% However, keep in mind that no measure is noise-free - there is always
% going to be noise in your data that you cannot get rid of and that's ok.
% If you have high quality data and many trials per condition, any
% remaining noise will mostly cancel out.

% Artifacts can be removed in several ways. In the context of EEG, this
% typically involves subtracting out noise fluctuations (recall, the final
% data is a summation of many different signals and noise, so we can try to
% identify the noise and simply subtract that out) or removing data points
% altogether (such as dropping trials that are overly contaminated with
% noise artifacts).
%          One can try to identify noise, for example, with independent
% components analysis (ICA), which will split your data up into so-called
% 'components' which capture different patterns in the variance of your
% data. In total, you can get as many components as you have electrodes
% (although you don't necessarily need that many components to perform the
% data cleaning). Each component (i.e. separated-out time series) will
% contain signal and noise, but some components may look primarily like
% noise. Once you've identified those components, you can subtract them out
% of your data. However, what constitutes noise is a judgment call.
%          Eye blinks are some of the most salient forms of noise in EEG
% data, as they produce large shifts in the EEG. There are ways to
% attenuate those shifts (e.g. with ICA or regression). However, even when
% blinks can easily be removed from the data, you may want to reject the
% data for that trial anyway, because blinks could indicate fatigue or lack
% of focus in some other way. To help minimize blinks in critical points of
% the data (i.e., during trials / times you plan to epoch), you can
% instruct the participant to minimize blinking during those times and give
% them sufficient breaks along the way to 'get their blinks out'.
% Naturally, depending on your participant, this may not be easily done
% (i.e., if you work with infants, kids, or older adults). When it comes to
% cognitive artifacts like fatigue or lack of focus, you can again offer
% sufficient breaks, interact with the participant during those breaks,
% ensure that the testing room is comfortable and cleaned up. Luck (2005)
% proposes that you play some music during breaks (or even during the
% experiment) that the participant picked out and accept the trade-off
% between noise from music with gain in focus and task performance.
%          Similar things can be said about oculomotor activity, like
% (micro)saccades. Participants will make more eye movements if they're
% distracted and/or tired. But even if they're not, saccades still occur.
% You can minize them by having experimental designs with visual stimuli
% centered on the monitor. Additionally, if participants are to give
% responses on a keyboard, make sure that they can find the buttons without
% having to look down. Ideally, you'd have a dedicated button box. But if
% you don't, simplify things by having people respond with the space-bar
% rather than a specific key that may be harder to find (and easier to
% lose). If you have an eye tracker available, consider using it to track
% participant's eye movements and blinks, to help you clean the data (or as
% second source of data for analysis, of course!).
%          Another source of noise is.. the rest of the body, which can
% also move! The odds of the IRB allowing you to strap down your
% participant are practically zero, so we also have to deal with
% electromyogram (EMG) noise, which results from physical movement of
% basically any kind and it's deleterious for EEG data above 15 Hz. Some
% EMG noise can be detected with ICA, some is restricted to certain
% electrodes, some is correlated with aspects of the experimental design
% (movement on some conditions versus others, if some conditions involve
% salient, scary, and/emotional stimuli for example). If possible, consider
% recording EMG seperately to help you identify what parts of the data
% contain movement noise. It can also help identify partial errors on
% trials.
%          If you thought that was it, think again! Remember those
% cognitive artifacts? Yup, there's more! Cognitive noise isn't just lack
% of focus, it's also errors on trials. You'll probably want to split out
% error trials and analyze those seperately (or simply remove them). Also
% take a look at response times, if responses needed to be given, are they
% trials with response times that are far away from typical? What does that
% mean? Was the participant not paying attention on that trial? It's a
% judgement call.

% Regardless, the point is again: high quality data is king.

% End of Chapter 8


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                               CHAPTER 9.                                %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Okay, let's get back to time-domain analyses and ERPs. Recall that data
% contains signal and noise and that, for ERPs, we epoch the data around an
% event (such as a stimulus onset) that will presumably introduce
% structural signal into our data that we can subsequently see by averaging
% a lot of trials together. If all you're doing is using ERPs as quality
% assessments of your data, that's all you need to know. But if you wanna
% make inferences about behavior, cognition, and/or the brain, have a look
% at Luck (2005) for a book on analyzing ERPs specifically.

% Also recall that time-domain signal averaging - the thing we do to create
% the ERPs - is like a low-pass filtering. When signals are out of phase,
% averaging across them tends to cancel things out and signals above 15 Hz
% tend to be non-phase-locked. Therefore, when we create ERPs, we end up
% with fluctuations up to 15 Hz. Some researchers will apply additional
% low-pass filters when computing ERPs, but this practice is still a topic
% of debate; filters can also produce artifacts (like ripple patterns) into
% your data, when the filter is poorly constructed. This means that,
% instead of interpreting neural signals, you end up trying to make
% inferences about an artifact that you yourself introducted. Typically,
% these additional filters have cutoffs around 20 or 30 Hz (but sometimes
% as low as 5 to 10 Hz). Another concern is that filters can shift the ERPs
% across time, making it seems like occurred earlier/later than they
% actually did. Thus, filters can reduce temporal precision of your
% results.

% Here's an example of how time-domain signal averaging is like a low-pass
% filter. We'll plot all the trials from a given electrode as well as its
% ERP.

load sampleEEGdata

% Instead of picking the channel number, let's pick by its location label.
% We'll use the function strcmpi() - 'string compare' - to find 'fcz' in
% the list of all labels under EEG.chanlocs.labels.
chanlabel = 'fcz';
channel = strcmpi(chanlabel,{EEG.chanlocs.labels});
plot(EEG.times,squeeze(EEG.data(channel,:,:)))
set(gca,'xlim',[-200 1000],'ydir','reverse') % Note that we reversed the y-axis
xline(0,"LineWidth",2)
title(strcat("All trials from electrode ",chanlabel))
ylabel('Voltage (µV)'); xlabel('Time (ms)')
hold on % Hold on, we're gonna plot the ERP too!

% Now we'll overlay the ERP
ERP = mean(squeeze(EEG.data(channel,:,:))');
plot(EEG.times,ERP,'black','LineWidth',5)

% There we are! As you can see, the wide variability across all the trials,
% once they're averaged into the ERP, is reduced -- as though it was
% low-pass filtered.

% Also note the reversed y-axis. So, why is that? Historical convention, as
% far as I can tell. When EEG first began being used, different labs had
% their own way of presenting data and.. one particular lab with the
% negative-upward convention won out. That's it. But that doesn't mean that
% you have to do so, too. Since then, a lot of research has been published
% with positive-up figures. The main take-away is that you always need to
% inspect /how/ things are plotted when reading others' work and also make
% it absolutely clear how you yourself are plotting things.

% So what about if we plot all the ERPs from each electrode - what would
% that look like? That's called a butterfly plot and here's what that looks
% like.
figure
subplot(311)
plot(EEG.times,squeeze(mean(EEG.data,3)))
set(gca,'xlim',[-200 1000],'ydir','reverse')
xline(0,"LineWidth",2)
xlabel('Time (ms)'), ylabel('Voltage (µV)')
title('Butterfly plot: ERP from all sensors')
% Take another look at how the ERPs are created, we're averaging along the
% 3rd dimension of the data. EEG.data is 64 x 640 x 99, so the 3rd
% dimension is all the individual trials of the 64 electrodes with 640 time
% points. So, when we collapse across that 3rd dimension, we're left with
% an array of size 64 x 640 (after squeezing) and we can immediately plot
% all those 64 ERPs in one go!

% As you can see, there's quite a bit of variance among the ERPs, although
% there is some structure there that is a bit more visible now. To see some
% more of that variance, try this:
subplot(312)
plot(EEG.times,squeeze(var(mean(EEG.data,3))))
set(gca,'xlim',[-200 1000])
xline(0,"LineWidth",2)
xlabel('Time (ms)'), ylabel('Voltage (µV)')
title('Topographical variance')
% Don't flip the y-axis here, that wouldn't make sense. We're looking at
% variation, so that's a non-negative number by definition. Same goes for
% the standard deviation in the next plot.

% Instead of var(), you can get the std() for global field power
subplot(313)
plot(EEG.times,squeeze(std(mean(EEG.data,3))))
set(gca,'xlim',[-200 1000])
xline(0,"LineWidth",2)
xlabel('Time (ms)'), ylabel('Voltage (µV)')
title('Global field power')

% Perhaps, at this point, you're curious about seeing what things like ERPs
% look like /spatially/, across the scalp. Some ERPs have strong
% inflections upward or downward (depending on how you plotted the y-axis)
% in the first second or so after stimulus onset. To investigate the
% topography of your data, you can make topoplots. We'll do that using free
% software: EEGLab.
%          As you make topoplots for various time points, you may notice
% that the overall pattern can stay (somewhat) stable for some period of
% time (e.g. for 200 ms, but also shorter or longer). These are called
% 'microstates' and such microstates have been linked to various
% (cognitive) processes.

%% Download and install EEGLAB
% EEGLAB is the software we'll be using to create topographical plots
% of our analyses. Fortunately, downloading and installing EEGLAB
% is easy!
% Step 1. Download EEGLAB: https://sccn.ucsd.edu/eeglab/download.php
%         Enter your info and you'll be taken to a page where you can
%         immediately download EEGLAB as a compressed file (e.g. zip).
% Step 2. Once downloaded, unzip the file and place it in a convenient
%         folder. For instance, near your Matlab installation.
% Step 3. Add the path to EEGLAB as follows (change the path here to
%         be the folder on your own computer):
addpath(genpath('/Users/hamid/matlab/eeglab2023.1'))

% And that's it! Now we can use all the functions in EEGLAB!

%% Exercises

load sampleEEGdata % If you haven't already

% 9.8.1
% First, let's compute the ERP at each electrode. Take a look again at
% the EEG structure in Matlab. Recall we have 64 channels, 99 trials,
% and 640 time points in our 99 trials ("epochs"). Under EEG.chanlocs,
% we can find the labels of the electrodes (which correspond to locations
% on the scalp).

% Since we have 64 electrodes, we'll end up with 64 ERPs after computing
% the average of 99 trials per electrode.
ERPs = zeros(EEG.nbchan,EEG.pnts); % This initializes our final ERP array.
for elec = 1:EEG.nbchan % For each electrode..

    % ..extract all the epochs, turn it into a 99x640 array
    % by squeezing and transposing..
    all_epochs = squeeze(EEG.data(elec,:,:))';

    % ..and then average across the 99 trials for each of the 640 time
    % points, which is the ERP at that electrode. Since time points are
    % in the columns, this is easily done, because Matlab averages
    % column-wise. Once computed, put this ERP aside (in the array that 
    % will hold all of our 64 ERPs.
    ERPs(elec,:) = mean(all_epochs);

end

% Now that we have all the ERPs, select five time points at which to show
% topographical plots from 0 to 400 ms in 100-ms steps. Let's find those 5
% time points first. Looking at EEG.times, we see that the sampled time
% points don't neatly map onto 100-ms steps from 0 to 400. So, we need to
% find the time points that are closest to our desired time-points-to-plot.

% From EEG.times, we can get the column info as follows. This isn't a
% 'standard'/'typical' procedure exactly. We're subtracting the value we're
% looking for from all the values in EEG.times and then simply finding the
% index for where our sought-for value had the lowest difference with
% values in EEG.times.
[~, idx] = min(abs(EEG.times-0)) % 257
[~, idx] = min(abs(EEG.times-100)) % 283
[~, idx] = min(abs(EEG.times-200)) % 308
[~, idx] = min(abs(EEG.times-300)) % 334
[~, idx] = min(abs(EEG.times-400)) % 359
idx = [257 283 308 334 359]; % Our column indices
% This means that the 0 ms data can be found in column 257 of our ERP data.
% Time point 100 ms is in column 283, and so on.

% Now we use EEGLAB's topoplot() function to make our topographical plots
% and we'll be done!
figure
cbar_range = [-10 10]; % Range of our colorbar
subplot(321); topoplot(ERPs(:,idx(1)),EEG.chanlocs); colorbar; caxis(cbar_range); title("ERP around 0ms")
subplot(322); topoplot(ERPs(:,idx(2)),EEG.chanlocs); colorbar; caxis(cbar_range); title("ERP around 100ms")
subplot(323); topoplot(ERPs(:,idx(3)),EEG.chanlocs); colorbar; caxis(cbar_range); title("ERP around 200ms")
subplot(324); topoplot(ERPs(:,idx(4)),EEG.chanlocs); colorbar; caxis(cbar_range); title("ERP around 300ms")
subplot(325); topoplot(ERPs(:,idx(5)),EEG.chanlocs); colorbar; caxis(cbar_range); title("ERP around 400ms")
% That's it! Compare your output to the exercise solutions by M.X. Cohen.
% The plotted patterns should look (essentially) the same.

% To increase the signal-to-noise ratio, make each plot show the average of
% activity from 20 ms before until 20 ms after each time point. Recall that
% each time step is 3.9062 ms. So, 20/3.9062 = 5.1201 time points before
% and after the point of interest. Let's say 5 points.
figure
cbar_range = [-10 10]; % Range of our colorbar
TOIs = [0 100 200 300 400]; % Time points of interest
% Did you know you can loop through subplots, by the way? Here's how:
for p = 1:5 % For each of the 5 plots..
    subplot(3,2,p) % Initialize the p-th subplot
    % Calculate the mean from -5 to +5 time points
    % Mean(x,2) computes means row-wise, which is what we need, since each
    % row is an electrode and we want to do this for each electrode
    % independently. The 2 signifies the '2nd dimension' of the array,
    % where the 1st is the columns.
    topoplot(mean(ERPs(:,idx(p)-5:idx(p)+5),2),EEG.chanlocs)
    colorbar; caxis(cbar_range);
    title(strcat("ERP around ",num2str(TOIs(p)),"ms")) % Indicate center time point
end
% The plots now look essentially like his solutions. The 'hotter' areas on
% the topoplots (i.e., those with higher values) have mellowed out, because
% those values were averaged with neighboring time points where values were
% lower.

% 9.8.2
% Loop through each electrode and find the peak time of the ERP between 100
% and 400 ms. Make a topo plot of those peak times. Recall from 9.8.1 that
% the indices of the 100ms and 400ms time points are 283 and 359,
% respectively.
peaks = zeros(EEG.nbchan,1); % This initializes our final peak array.
for elec = 1:EEG.nbchan % For each electrode..

    % ..find the index of max (i.e., peak) value of the ERP in the 100-400
    % ms time window..
    [~, idx] = max(ERPs(elec,283:359));
    % Note now that idx is relative to the window from 283:359. So, if the
    % peak is in column 283, this will return idx = 1. So, to correct for
    % that, we add 283-1 = 282.
    idx = idx+282;

    % ..and store the final peak time
    peaks(elec) = EEG.times(idx);

end

% And now we can plot those peak times
figure
cbar_range = [100 400];
topoplot(peaks,EEG.chanlocs);
colorbar; caxis(cbar_range)
ylabel(colorbar,'Latency (ms)','FontSize',15,'Rotation',270,'VerticalAlignment',"bottom")
title("ERP peak latency",'FontSize',18)

% Areas showing the earliest peak responses are more frontal, but also
% posterior-medial above the precuneus. The latest peaks are occurring more
% occipitally and temporally.

% End of Chapter 9
clear all; close all; clc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                        %
%%                               CHAPTER 10.                              %
%%                                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Back to some math - sorry! But don't worry, this is gonna turn out to be
% quite straightforward. It's also gonna turn out to be very important. I
% cannot understate how important this topic is. Not only does the
% following math support our Fourier transform calculations (which are
% needed to do frequency-based, power-based, and phase-based analyses),
% the math in this chapter actually stands at the heart of.. everything we
% do as scientists interested in behavior and the brain. No, seriously.
% It's the same mathematical operation that we perform when we do things
% like linear regression. This is math that we have to know: the dot
% product.

%% The dot product
% The book mentions there are several interpretations of the dot product,
% but let's look at the actual mathematical operation first.

% Let's say we have four numbers, paired in twos.
set1 = [1 2];
set2 = [3 4];

% To compute the dot product, first, we multiply each corresponding entry 
% in the sets with each other. So, multiply the first number of set1 with
% the first number in set2. Then, the second ones. If there were more
% numbers in each set, we would simply continue multiplying the entries --
% the 3rd entries, the 4th entries, and so on. This also means that, in
% order to compute a dot product, we need sets of numbers that are equal in
% length.

% In this case, that results in the numbers 3 and 8.
% 1 * 3 = 3
% 2 * 4 = 8

% Next, we.. just add all those numbers together and we're done!
% 3 + 8 = 11
% And so, the dot product between set1 and set2 is 11.

% You might think that, for the multiplication step, we can just do:
set1 * set2
% Try it out, you'll get an error saying 'to perform elementwise
% multiplication, use '.*' '
% In other words, don't just write * (i.e. 'the product of set1 and set2').
% Instead, write 'dot' and then 'product': .*
set1 .* set2 % 'dot product'
% Okay, I lied. It doesn't give you the actual dot product, it only does
% the multiplication, but not the summation. However, you'll see .* in many
% scripts and it's obviously related to the actual dot product operation.
% So, you should be aware of what .* does in Matlab.
% Regardless, to get the actual dot product in this case, we need:
sum(set1 .* set2)

% You can, technically, just use the * symbol, if the first set is a row
% and the second set is a column of numbers. Just use transpose:
set1 * set2'

% Or, simply use dot()
dot(set1,set2)

% This should remind you of the linear functions we make when we do
% regressions and so on:
% y ~ b_1*x_1 + b_2*x_2 + b_3*x_3 + constant + error
% What we have there is a dot product:
% a summation of multiplications (b_1*x_1 + b_2*x_2 + b_3*x_3)!
% That's why it's so important to also understand this concept outside of 
% the context of this book. We use it all the time.

clear all; close all; clc

%% But the question now is, what does that dot product ~*"mean"*~?!
%% 1. The dot product is 'the sum of one vector weighted by another'
% The first interpretation is as simple as it gets: it's the sum.. of the
% elements in one vector weighted by the other vector. It's really about
% the weighting of one vector by another. If one of the sets of numbers
% was all ones, those 'weights' wouldn't change the values of the other set
% of numbers. But for any number that isn't one, we'd get a change of the
% value in the corresponding value in the other set. This is why the beta
% coefficients in the linear model above are also called beta 'weights'.
% Coefficients are weights, weights are coefficients, beta's are weights
% are coefficients. It's all just different terms for the same things.
% Often the connection between all of these concepts is missed because it's
% made obscure by the difference in jargon. But the underlying math shows
% that it's all about the same thing. This perspective on the dot product
% is 'algebraic', it's focused simply on the algebra of the operation.

%% 2. The dot product is 'the covariance or similarity between vectors'
% Covariance is the measure of 'joint variability' of two random variables.
% It's similar to correlation, which is the covariance divided by the
% product of the standard deviations of both variables, putting the
% variables on the same scale (and producing a value from -1 to +1).
% But to compute the covariance, we compute the deviations of each variable
% from their own means, multiply the deviations elementwise and sum them
% all together (!!!!!), and then divide by n-1. Let's look:
var1 = rand(1,50);
var2 = rand(1,50);
% Two variables (var1, var2) with 50 random values between 0 and 1.
cov(var1,var2) % produces the covariance matrix
%    0.0751    0.0151
%    0.0151    0.0961
% In my case, I got the above. You will probably get something different,
% given the random draw of numbers for var1 and var2.

% Subtract out the mean, leaving the deviations from the mean. And then,
% compute the dot product.
var1_dev = var1-mean(var1);
var2_dev = var2-mean(var2);
dot_devs = dot(var1_dev,var2_dev);
% And finally, divide by n-1, which is (50-1=) 49.
dot_devs/49 % 0.0151 (again, your output is likely different)
% You'll see this value appear twice in the covar. matrix above, the bottom
% left and upper right. This is because cov() computes self-covariance too.
% This is the matrix with labels:
%        var1      var2
% var1   0.0751    0.0151
% var2   0.0151    0.0961
% Indeed, if I ran cov(var1), it returns 0.0751.

% What about cov(var1, -1 * var2)?
cov(var1,-var2)
%    0.0751   -0.0151
%   -0.0151    0.0961
% We've flipped the sign of one of the variables and that flips the sign of
% the covariance. The 'similarity' (or lack there of) interpretation is
% that covariance (and correlation) say something about whether two vectors
% behave 'similarly' or not.

% Let's say that we have two variables with means of zero: both the
% covariance and the dot() will be zero.
var3=[3 -3 0]; var4=[4 4 -8];
cov(var3,var4)
%     9     0
%     0    48
dot(var3,var4)
%     0

% When the means of both variables are positive, dot() is positive
dot(var1,var2)
%    14.5051
% When both means are negative, dot() is also positive
dot(-var1,-var2)
%   14.5051

% When one mean is positive, one is negative, dot() is negative
dot(var1,-var2)
dot(-var1,var2)
%   -14.5051

% Thus, if the sign of means of the variables matches, the sign of dot() is
% positive. If they don't match, it's negative. If they're both zero, dot()
% is zero. So again, the covariance (and correlations) - and hence, the dot 
% product - say something about 'similarity' in behavior of two vectors.

clear all; close all; clc

%% 3. The dot product is 'a (geometric) mapping between two vectors'
% This is, potentially, the most relevant way to think about the dot
% product when it comes to the content of this book. Here, the dot product
% is the multiplication of two vectors scaled by the cosine angle between
% them. Let's see how that works.

% You may have heard that a set of (however many) numbers can be seen as a
% vector. Let's start simple, with a set of two numbers.
vec1 = [1 2];

% Plot the vector. Vectors can be thought of as a number (i.e. a
% 'magnitude') with a direction. We'll plot it with an arrow to indicate
% that. Matlab's quiver() can help us out, it plots the starting point of
% the arrow with the first two numbers (starting coordinates) and its end
% point with the next two numbers (i.e., end-point coordinates). We add 
% 'off' to turn off scaling of the vector, that's why that's there -- but 
% it's not relevant for our purposes. See below:
figure
quiver(0,0,vec1(1),vec1(2),'off','linewidth',5) % Starting from the origin to coordinates (1,2)
xlim([-1 3]); ylim([-1 3]) % Adjust the range of the axes
hold on; text(vec1(1)+.2,vec1(2),'Vector 1') % Add a label

% Great, now let's add another vector. Recall, for the dot(), we need
% vectors of the same length (in this case, 2 elements).
vec2 = [2 0];
hold on; quiver(0,0,vec2(1),vec2(2),'off','linewidth',5)
hold on; text(vec2(1)+.2,vec2(2),'Vector 2')

% Okay, the dot() between these two vectors is:
dot(vec1,vec2) % 2

% So what does the 2 mean here? And where does the part about 'scaled
% cosine angle' come in? Well, right now! We're going to 'project' the
% first vector 'onto' the second one. This effectively comes down to
% drawing a line from the end-point of the first one 'down' to the second
% one to create a right angle with that second vector.
hold on; line([1 1],[2 0]) % Projection line
hold on; line([.8 .8],[0 0.2]); hold on; line([.8 1],[.2 0.2]); % right angle
hold on; text(vec1(1)+.2,vec1(2)/2,'Perpendicular line')
% A more fun way to think about it is that we place a light above vector 1
% and point it directly at vector 2. The perpendicular line we just drew is
% like a ray of light. So then what's the projection? Well, it's the shadow
% of the first vector on the second one. We've "projected" that first
% vector onto the second one!
hold on; quiver(0,0,1,0,'off','linewidth',5)
hold on; text(.15,-.2,'Projection of Vector 1')

% As you can see, the projection of vec1 lands on coordinates (1,0) and so,
% it has a magnitude of 1. The magnitude of vec2 was 2. Their
% multiplication? 1*2 = 2.. which is the dot product!

% In other words, the dot product between two vectors is the multiplication
% of their magnitude, once one of the vectors has been projected onto the
% other. Think about this for a moment. What's the range of values that the
% projection can take on? Its longest possible length is the length of the
% original vector itself (if it was lying right on top of the other vector,
% the projection would be exactly the length of the original vector
% itself). If vec1 was pointing exactly perpendicular to vec2, the
% projection of vec1 would be zero. If the angle between them was greater
% than 90 degrees, it would project onto the line that vec2 is on, but the
% projection would be pointing the other way, like so:
figure
quiver(0,0,-vec1(1),vec1(2),'off','linewidth',5)
xlim([-1.5 3]); ylim([-1.5 3])
hold on; text(-vec1(1)+.2,vec1(2),'Vector 1')
hold on; quiver(0,0,vec2(1),vec2(2),'off','linewidth',5)
hold on; text(vec2(1)+.2,vec2(2),'Vector 2')
hold on; quiver(0,0,-1,0,'off','linewidth',5)
hold on; text(-1,-.2,'Projection of Vector 1')

% Thus, the more these two vectors overlap, the larger their eventual 
% multiplication. They are like 'forces' with directions being multiplied 
% together and that final number that comes out of the multiplication 
% represents how much these vectors are 'working together' or, better yet, 
% how much they 'REPRESENT A SIMILAR SIGNAL'! And, of course, the number
% that comes out of this multiplication-after-projection just so happens to
% match the dot product (because, of course, it IS the dot product).
 
% So, if the vectors are a little off, like in our
% first example, the value that represents their 'working-togetherness' (or
% their 'similarity') gets attenuated (they're considered less similar). If
% the vectors are at a right angle (i.e. they're orthogonal), they're not
% working together at all - one vector is pushing things in one direction, 
% the other pushing things in an orthogonal direction. They're not similar.
% So, to be clear, they are orthogonal, NOT opposing - they have no effect 
% on each other. But once the angle is greater than 90 degrees (like in our 
% second example), the vectors are actually working AGAINST each other. So,
% they are "similar", in that they are closer to lying on the same infinite
% line/plane, just in opposing directions. And the more the angle is closer 
% to 180 degrees, the more the vectors are 'similar-but-opposing').
%          Taking a step back, what's our eventual goal again? We have to
% figure out how to decompose an EEG signal into its constituent
% frequencies. So, we ask 'How similar is our EEG signal to a wave of.. 1
% Hz? And of 2 Hz? And of 3 Hz?' And so on. Well, the dot product between
% our EEG data and each of those frequencies is gonna give us a measure of
% similarity.

% By the way, in case you didn't notice, the point where vec1's projection
% falls onto vec2 is simply the x-coordinate of vec1. So, since:
% vec1 = [1 2] and vec2 = [2 0]
% it seems like we could've just taken the x-coordinate of vec1 (i.e. 1)
% and multiplied it with the x-coordinate of vec2 (i.e. 2). Well.. yeah!
% Like I said, it's no coincidence that projecting a vector and then
% multiplying the projection with another vector.. results in the same
% thing as elementwise multiplication followed by summation. They are
% different ways to perform the SAME mathematical operation.  It is exactly 
% what we're doing in the first (algebraic) perspective on
% the dot product; if we see each of the elements of the vectors as a
% 'dimension' or 'component', we are computing the product of the
% magnitudes in each dimension and then summing them all together. When we
% do sum(vec1 .* vec2) we are asking 'What's the similarity in
% dimension/component 1? And in 2? And in 3? ... And in n? Okay, now let's
% sum all those similarities together'. And that final dot product is like
% a 'sum of similarities' (that's my term for it, it's copyrighted). Heads
% up: you may remember that when we get things like a 'sum of squares', we
% divide it with some value, to scale it down -- we're gonna be doing that
% with the dot product later too.

% Because of the projection, another way to look at it is to see vec2 as a 
% number line and we're linearly transforming vec1 to place it on the
% number line that is vec2. In the case of the example, vec2 is already on
% our typical real number line (the current x-axis). But you can also think
% vec1 being a new number line and us transforming vec2 to be placed onto
% this angled number line instead. Indeed, this is akin to what we do with
% principal components analysis; we're taking our original dimensions
% ("components") and rotating our current number lines to better match the
% (direction of the) variance in the data. But let's not get too
% philosophical.

% Hold on, what about cosine?! What happened to that? Let's return to our
% vector figures. Recall from our basic trigonometry review in Chapter 3,
% once we have the projection of one vector onto another, we can use cosine 
% to calculate the length of the projection. After all, the cosine is the 
% adjacent side (of a triangle) divided by the hypotenuse:
%      cos(angle) = adjacent / hypotenuse
% This also means that (if we multiply both sides by the hypotenuse, we can
% define the adjacent side:
%      hypotenuse * cos(angle) = adjacent
% Therefore, if we have the length of the hypotenuse and cos(angle), we
% know the length of the adjacent side (which is the projection).

% In the unit circle, the hypotenuse is 1, but we're not on the unit cirlce, 
% so let's use the Pythagorean theorem on vec1=[1 2] 
% sqrt(1*1 + 2*2) = 2.2361. So, the length of the vec1 is 2.2361
%           The angle between the vec1 and vec2, in our example, is 63.4349 
% degrees and cosd() allows us to enter the angle in degrees (as opposed to 
% radians with cos()) and compute the cosine. So:
cosd(63.4349)
%     0.4472

% So, to get the length of the projection:
2.2361 * cosd(63.4349) %    1.0000

% The dot product is then defined as:
% dot(vec1,vec2) = ||vec1|| * cos(angle between vec1,vec2) * ||vec2||
% where || || means the 'magnitude' of the vector, meaning the length of
% the arrow from the origin to the end. So, we take the length of vec1
% (which is ||vec1||), scale it down (i.e., get the projection,
% which is to say, multiply ||vec1|| by cos(angle)), and then multiply the 
% scaled-down-version-of-vec1 by vec2.

% This magnitude is also known as the 'norm' of a vector and is easily
% computed with norm().
norm(vec1) % = 2.2361
norm(vec2) % = 2
norm(vec1) * norm(vec2) * cosd(63.4349) % = 2
% And there you have it, it equals 2. The order of the terms above is
% typically stated with the cos(angle) term at the end:
% dot(vec1,vec2) = ||vec1|| * ||vec2|| * cos(angle)
% because it doesn't matter whether we project vec1 onto vec2 or vec2 onto
% vec1, both result in the same single number that represents their
% similarity.

% Although I gave you the angle between the vectors, it can be computed as
% follows, in case you're wondering:
cos_of_angle = dot(vec1,vec2) / ( norm(vec1) * norm(vec2) ); % 0.4472
angle_in_degrees = real(acosd(cos_of_angle)); % 63.4379
% In other words, if you don't know the angle between vectors.. you can
% always compute it by starting with the dot product and then dividing that
% by the product of the vector magnitudes. Speaking of which..

%% Bonus: Cosine similarity
% At this point, you might be wondering: "Hey, isn't there a thing called
% cosine similarity? Does it have anything to do with this?". Well guess
% what, buddy, it has EVERYTHING to do with it.

% Cosine similarity is.. the cosine of the angle between two vectors. The
% cosine similarity between vec1 and vec2? It's cosd(63.4349) = 0.4472.
% Like I said, if you don't know the angle between vectors, you can easily 
% compute it now:
% cos(angle) = dot(vec1,vec2) / ( norm(vec1) * norm(vec2) )

% What's nice about this measure is that it puts everything into a range
% from -1 to +1. Anyway, boom - you now know what 'cosine similarity' is
% and how to compute it.

clear all; close all; clc

%% Convolution

% Okay, enough of that for now. Let's move on to convolution! What's that?
% Well.. it's essentially a /sliding/ dot product. Remember how we took a
% step back earlier to look at the bigger picture? We have an EEG signal
% consisting of a summation of various frequencies and we can compare that
% to pure waves of frequencies to compute how similar our EEG data is to a
% frequency of - say - 47 Hz. But also 68 Hz. Or 3 Hz. Any frequencies
% we're interested in, really. To do so, we can compute that similarity
% (i.e. the dot product) at a given time point between our data (which
% we'll call the 'signal' [cause that's where our signal is thought to be])
% and a wave of a frequency we're comparing it to (the 'kernel', which
% comes from the word for 'seed', with the signal potentially containing
% features that "stem/sprouted from" the kernel - I don't know if that's 
% the 'seed' of the term, but I'm rolling with it. 'Signal' comes from 
% the Latin 'signum' meaning mark or token, so it bears the mark of the 
% kernel. Makes sense to me! Don't @ me.).

% Anyway, we compare the signal to the kernel at a given time point and 
% then.. slide the kernel to next time point and compute dot() again.. move 
% it over and compute it again.. and so on.. checking the extent to which 
% the signal bears the mark of the kernel (i.e. their 'similarity') each
% time. Let's start by looking at the basic mechanics of the computation:
x = 1:100;
signal = sind(2*pi*x);
figure; plot(signal)
kernel = [1 2 3 4 5];
% Our signal is a sine wave, our kernel is a line. Notice that the kernel
% and the signal don't have the same length. I said earlier that, in order
% to compute the dot(), we need vectors of the same length. That still
% holds: we're just gonna compute the dot product between the signal and
% kernel /where they overlap/ as we're sliding the kernel along. The books
% mentions zero-padding, which simply results in the non-overlapping parts
% being weighted by zero (and thus being canceled out). Either way, the
% result is the same.

% Imagine we start by multiplying the first entry of our signal with that
% of our kernel:
kernel(1) * signal(1) %     0.1094
% Next, the kernel moves up so that we get:
kernel(1) * signal(2) + kernel(2) * signal(1) %     0.4365
% Make sure that you understand what's happening there.

% First time point, only the first elements overlap:
% Kernel: [5   4    3   2   1]
%                           ^
% Signal:                 [0.1094  0.2175  0.3230  0.4247  0.5212 ...]
% So, 1 * 0.1094 = 0.1094

% Second time point, the first two elements overlap:
% Kernel: [5     4      3      2      1]
%                              ^      ^
% Signal:                 [0.1094  0.2175  0.3230  0.4247  0.5212 ...]
% So, 1*0.2175 + 2*0.1094  = 0.4365

% And so on and so forth, until we have reached the last time point where
% only the last entry of the kernel (the 5) overlaps with the last entry of
% the signal. Importantly, NOTICE THAT WE FLIPPED THE KERNEL AROUND WHEN
% VISUALIZING THE COMPUTATION! This should make sense, I hope.

% Matlab has an easy function that performs the entire convolution, at all
% time points, for us:
conv_signal = conv(signal,kernel);
% Confirm:
conv_signal(1) %    0.1094
conv_signal(2) %    0.4365
kernel(end)*signal(end) %  -4.9978 (last kernel entry * last signal entry)
conv_signal(end) %  -4.9978 (last convolution entry, matches previous line)

% But note that length(conv_signal) = 104. Our signal has length 100, our
% kernel has 5. So what happened? The first convolution is between the
% first entries of each vector. By the time kernel(1) reaches signal(100),
% a hundred dot()'s have been computed.. but then there are still several
% trailing entries in the kernel that still need to slide forward. As a
% matter of fact, the number of trailing entries in the kernel is the
% length(kernel) - 1 (there's 5 entries, meaning 4 trailing after the first
% entry). And so, we get 4 more dot products:
% length(convolution) = length(signal) + length(kernel) - 1
%           We need to fix that! Simply, we trim the ends of the
% convolution with the floor(length(kernel)/2)). In this case, 2.
trim = floor(length(kernel)/2);
conv_signal = conv_signal(trim+1:end-trim);
% This way, we end up with a convolution of the same length as the signal.

% You can also request that Matlab does this trimming for you
% automatically, as follows:
conv_signal_trimmed = conv(signal,kernel,'same');
% This trims to convolution, with the same approach as above, so that it's
% the same length as the first argument passed to conv(), in this case
% signal. So, always make sure that you give conv() the data first and the
% kernel second.

% But what does it all "mean", man!?!? Now that you know that dot products
% are about similarity, you should be able to guess that one interpretation
% of a convolution is that it's.. merely a sliding similarity measure.
figure; plot(signal); hold on; plot(conv_signal)
% The signal fluctuations between -1 and +1, the convolution between -15
% and +15. When we move a positive vector (which is what our kernel is)
% across the signal, the convolution is telling us that the kernel is more
% similar to the signal.. when the signal is also positive. Our convolution
% is strongly negative, when.. the positive kernel is convolved ("dotted")
% with a signal that is negative. Think back to the earlier points on the
% sign of the dot product when looking at whether the sign of the averages
% of vectors match or not.

% Okay, so what about this scaling the convolution thing? With all this
% multiplying and summing going on, the dot() can quickly get very big. So,
% we scale it down, to more easily compare it to the original signal that
% we convolved. To do so, we divide it by the sum of the kernel, as such:
conv_scaled = conv_signal_trimmed / sum(kernel);

% The scaling works also if the kernel is normalized to have a sum of 1.
kern_scaled = kernel / norm(kernel,1);
conv_kernscaled = conv(signal,kern_scaled,'same');

% Let's see the results
figure
subplot(411); plot(signal); hold on; plot(conv_signal_trimmed)
legend('Signal','Trimmed Convolution'); title('Unscaled')
subplot(412); plot(signal); hold on; plot(conv_scaled)
legend('Signal','Trimmed Convolution'); title('Scaled with /sum(kernel)')
subplot(413); plot(signal); hold on; plot(conv_kernscaled)
legend('Signal','Trimmed Convolution'); title('Scaled with kern-scaled')

% As you can see, the signal and its convolution with the (scaled) kernel,
% are now back on the same range on the y-axis, so we can more compare
% them. Now flip the kernel's sign and see what happens.
kern_flipped = -kern_scaled;
conv_kernflipped = conv(signal,kern_flipped,'same');
subplot(414); plot(signal); hold on; plot(conv_kernflipped)
legend('Signal','Trimmed Convolution'); title('Scaled with kern-scaled')

% The signal and the original kernel (which was a vector with positive,
% increasing numbers) is most similar to the signal when that signal is
% also steadily increasing. When we flip the sign of the kernel, the kernel
% becomes a steadily decreasing vector, which is most similar to parts of
% the signal that are also steadily decreasing (so the convolution is
% highest, when the signal is going down).

clear all; close all; clc

%% Exercises

% 10.6.1
% Create two kernels for convolution: inverted U and a decay function. No
% need to be too sophisticated, numerical approximations are fine.
kern_decay = [1 .9 .8 .7 .5 .3 .2 .1 0];
kern_invU = [1 .8 .3 .1 0 .1 .3 .8 1];
figure; plot(kern_invU); hold on; plot(kern_decay)

% 10.6.2
% Convolve these kernels with 50 time points of EEG data from one
% electrode. The 'snippet of EEG data' he plots is channel 47, time points 
% 100:149, trial 10.
load sampleEEGdata
signal = squeeze(EEG.data(47,100:149,10));
conv_invU = conv(signal,kern_invU,'same')/sum(kern_invU);
conv_decay = conv(signal,kern_decay,'same')/sum(kern_invU);

% Plot the kernels, EEG data, and the result of the convolution.
figure
% Kernels
subplot(311); plot(kern_invU,'-o'); hold on; plot(kern_decay,'-o')
legend('invertedU','decay')
title('Convolution kernels')
xlabel('Points')
hold on
% EEG data
subplot(312); plot(signal);
title('Snippet of EEG data (EEG.data(47.100:149,10))')
xlabel('Time points'); ylabel('Voltage µV')
hold on
% Convolution result
subplot(313); plot(signal, 'Color', 'black'); hold on % EEG again
plot(conv_invU,'Color','blue','LineWidth',3); hold on % Convolution 1
plot(conv_decay,'Color','red','LineWidth',3); hold on % Convolution 2
legend('EEG data','invertedU-convolved','decay-convolved')
title('EEG data before and after convolution')
xlabel('Time points'); ylabel('Voltage µV')

% Based on visual inspection, what is the effect of convolving the EEG data
% with these two kernels?
% The convolved signals appear low-pass filtered ('smoothed out'). The
% inverted-U kernel allows for more troughs and peaks in the final result, 
% compared to the decay-convolved signal. Thus, the features of the given 
% kernel are allowed to 'pass through', if you will, if they match those of
% the signal (locally). Note that at the starting and ending time points,
% the convolved results appear inverted (troughs became peaks and vice
% versa), especially in the case of the inverted-U convolution (but also
% with the decay-convolution at the end).

% End of Chapter 10
clear all; close all; clc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                        %
%%                               CHAPTER 11.                              %
%%                                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We're almost ready to learn about Fourier transform, the backbone of
% (most) EEG analyses. It works by computing the dot product between the
% EEG data (the signal) and sine waves of a range of frequencies (the 
% kernels). Recall from Chapter 3 that EEG data have frequency, power,
% and phase. The result of a Fourier transform is a three-dimensional
% representation of the frequency, power, and phase similarity between the
% signal and the kernel.

% The remaining thing to understand is the concepts of imaginary and
% complex numbers. Recall from the Chapter 10 how we briefly talked about
% number lines. Traditionally, we have one number line, visualized as.. a
% horizontal line centered on 0, with positive numbers off to the right and
% negative numbers to the left:

%     <<< ... -3  -2  -1  0  +1  +2  +3 ... >>>

% But we can, in principle, place a number line wherever we want. You might
% think 'Well, that's like Cartesian coordinates, isn't it? An x-axis and a
% y-axis? Two number lines?" Well, not quite. In that case, we're
% visualizing how a function takes us from one point on the traditional
% number line to another point on the same traditional number line. But
% here, we're talking about a completely separate number line altogether.
% But since we're familiar with Cartesian coordinates, let's /imagine/ that
% we have a completely different number line at a right angle to our first
% number line:

%                       +3i
%
%                       +2i
%
%                       +1i      o
% 
%            -3  -2  -1  0  +1  +2  +3
%
%                       -1i
%
%                       -2i
%
%                       -3i


% This second number line is called the 'imaginary number line', with the
% numbers along it called 'imaginary numbers'. The traditional number line 
% is then called the 'real number line' with 'real numbers'. In my opinion, 
% these are a TERRIBLE names, just the absolute worst, but it is what it is. 
% Regardless, the imaginary number line has a unique property where the 'i' 
% represents the 'imaginary unit: square root of -1'. The point is, now we 
% can express numbers in relation to these two number lines (as opposed to 
% in relation to just the real number line).

% The number 2 now gets written as '2 + 0i', meaning that it's 2 on the
% real number line and 0 on the imaginary line (hence the i, to indicate
% the imaginary line). The expression of a number is now a bit more complex
% as a you can see. Well, guess what -- we call these numbers 'complex 
% numbers'. The complex number 2+1i is shown above with a dot. And the 
% 2-dimensional plane that these complex numbers lie on is called the 
% 'complex plane'. Again, just 0 out of 10, two-thumbs-down on the naming 
% convention. Point is, complex numbers come in the form 'a + bi', where
% the 'a' term is our "real" component and the 'b*i' term is our
% "imaginary" component. Apparently, it was René Descartes who came up with
% the name 'imaginary numbers', because he - and many other mathematicians
% at the time - thought the idea was ludicrous and that those numbers were
% pointless.

% Now note that we can also use trigonometry to draw vectors, circles, and
% - hence - periodic functions, but.. on the complex plane! So, whereas
% previously, the cos(angle) was on the x-axis and the sin(angle) on the
% y-axis, we now have cos on the real line and sin on the imaginary line.
% Writing complex numbers in aforementioned way [ a + bi ] and plotting
% them in a 2-dimensional complex plane is of course similar to the x and y
% Cartesian coordinate system. And so, 'a + bi' is known as the 'Cartesian
% form' of a complex number.

% In Matlab, we can define complex numbers several ways. Let's recreate the
% commented image above:
z = complex(2,1); % Real part in the first spot, imaginary in the second
plot(real(z),imag(z),'o')
grid on
xlim([-3 3]); ylim([-3 3]) % The x-axis now is Real, y-axis Imaginary
xlabel("Real(z)"); ylabel("Imag(z)")

% A complex number can also be written as follows:
z = -2 - 1i;
hold on; plot(real(z),imag(z),'+')

% Alright, now let's plot a unit circle on the complex plane
radius = 1;
theta = 0: pi/100 : 2*pi;
z = complex(radius * cos(theta), radius * sin(theta));
close all % close all other figures
figure; plot(real(z),imag(z))
grid on
xlim([-2 2]); ylim([-2 2])
xlabel("Real(z)"); ylabel("Imag(z)")

% From this, you should be able to see that any /point/ along this unit
% circle can be expressed as follows:
% radius( cos(angle) + i * sin(angle) )
% Thus, notice that the radius term disappears (because it equals 1 for the
% unit circle, and the i term is introduced to indicate that sin is on the
% imaginary axis.
%           Importantly, notice how the cosine and sine components map onto
% our earlier point about complex numbers (a + bi) -- the cosine maps onto
% the real axis and the sine onto the imaginary axis.
%           By the way, remember how - when creating a sine or cosine wave
% Chapter 3 - we wrote something like the following:
% .7 * sin(2*pi*x * 36 - .5) + 1
% That amplitude of 0.7 is effectively the radius. But then also notice
% that the imaginary 1i is the radius - and therefore amplitude - of the
% sine wave in the earlier expression [ i * sin(angle) ]. That amplitude
% can itself of course also be scaled [e.g., 5i * sin(angle) ].

% Okay, we're now going to be rewriting cos(angle) + i * sin(angle)
% to produce an expression that is more commonly used, using some rules from
% calculus. However, this other expression is - ultimately - the SAME as
% the one above. So, conceptually, nothing changes and all you need to
% understand is that expression above. If you understand that, you're all
% set. What happens next is NOT changing the above expression in any way in
% terms of what it's doing, it's just 'writing it differently'.
%           We're gonna show that  e^(i*angle) = cos(angle) + i*sin(angle).
% First, set the right hand side equal to a value 'z' and take the
% derivative, where theta is now the angle: 
%     z = cos(theta) + i*sin(theta)
%     dz/dtheta = -sin(theta) + i*cos(theta)
% Next, rearrange the derivative and use the definition of i, which is the
% imaginary 'square root of -1' to rewrite the '-' in front of sin(theta)
%     dz/dtheta = i*cos(theta) + i^2*sin(theta)
% Now we factor out the i:
%     dz/dtheta = i*(cos(theta) + i*sin(theta))
% The expression in the parentheses is now equal to the z we started with:
%     dz/dtheta = i*z
% This can be rearranged with algebra (multiply both sides by dtheta and
% divide by z):
%     dz/z = i*dtheta
% Now, with what we have on the left-hand side, we can integrate, which
% gives us:
%     ln(z) = i*theta + C (some constant)
% Which we can turn into:
%     e^(ln(z)) = i*theta + C
%     z = e^(i*theta) * C
%     cos(theta) + i*sin(theta) = e^(i*theta) * C
% When we put input 0 as our angle theta, this results in the constant C
% being shown to equal 1 (and so, it disappears).
% Therefore:
%     e^(i * theta) = cos(theta) + i*sin(theta)
% Thus, although we'll be continuing with the left-hand side as a way to
% express where we are on the unit circle in the complex plane, you can
% always interpret that expression as what we're seeing on the right-hand
% side, because they are the same thing (and the right-hand side is perhaps
% easier to grasp if you're not familiar with the alternative expression).
% Either way, this alternative form [ e^(i*theta) ] is known as the 'Polar
% form' of a complex number, because we're expressing the number in terms
% of a radius and magnitude from some central point (like we're on the
% North pole and trying to figure out in which direction and what distance
% one's home is). The 'pole' is a reference point (akin to the origin (0,0)
% in the Cartesian form) and the 'polar' is the line to/from a point on the
% plane that, in relation to an axis of the plane, creates an angle.

theta = 0: pi/100 : 2*pi; % In radians
z = exp(1i*theta); % To use the number e, we use exp().
figure; plot(real(z),imag(z))
grid on
xlim([-2 2]); ylim([-2 2])
xlabel("Real(z)"); ylabel("Imag(z)")
% As you can see, we get the same result as before

% Okay, but now look at this. Instead of plotting it in 2 dimensions, we'll
% plot our 'circle' in 3 dimensions:
figure
radius=1;
plot3(radius*real(z), radius*imag(z), theta,'-','linew',3)
grid on
xlabel('Real(z)')
ylabel('Imag(z)')
zlabel('Angle (radians)')
title('3-D representation of our complex unit circle')

% So it's not really a circle, so much as it is a corkscrew. You can use
% drag the 3d image around (left-mouse click on the figure, hold down, and
% drag around). Look at the axis labels and drag it around till you see the
% angle as a function of the real numbers (Real at the bottom, x-axis;
% Angle plotted on the y-axis). Doesn't this look like the line drawn by
% the cosine in our animation from Chapter 3? It should, because that IS the
% same curve. If you're having trouble rotating it, run the following:
view([0 90]) % x: Real, y: Imaginary

% Here is the animation again:
addpath '/Users/hamid/Desktop/Research/+Tutorials/ANTSD'
figure; ANTSD_animatedwaves(1)

% Similarly, if you rotate the 3d plot to put Imag(z) on the x-axis and
% Angle on the y-axis, you'll see the sine wave from our animation (albeit
% rotated, because the animation has the axes flipped). Indeed, the cosine
% maps onto the real values, the sine onto the imaginary values.
view([90 0]) % x: Imaginary, y: Angle

% Similar to our 2d animation, we can also create a 3d animation
freq = 6;
z = exp(1i*theta*freq);
figure; comet3(radius*real(z), radius*imag(z), theta)
xlabel('Real(z)')
ylabel('Imag(z)')
zlabel('Angle (radians)')
% You can enter any number you like as the frequency and the animated
% corkscrew will consist of that many 'cycles'. Notice that our animation
% is rotating upward (in a positive direction). Rotate the 3d plot again so
% that real is on the x-axis, the imaginary line on the y-axis (and keep in
% mind how the line rotated during the animation). You'll see that it
% rotated in a counter-clockwise manner.

% You may have guessed: we can also make the corkscrew move in the
% opposite, clockwise direction. All we need to do, is flip the sign in the
% expenonent as so:
z = exp(-1i*theta*freq); % Note the minus sign at the start in the exponent
figure; comet3(radius*real(z), radius*imag(z), theta)
xlabel('Real(z)')
ylabel('Imag(z)')
zlabel('Angle (radians)')
% Yes, the animation is still moving upward; the point here is about the
% /direction/ of the dot being clockwise or counter-clockwise and that
% moving clockwise /is like/ moving in the downward direction of the first
% animation (because going downward would equal going in the clockwise
% direction, when viewed from the perspective produced by view([90 0)). And
% we're, of course, not changing anything about the frequency. And that's
% critical, because we can move the negative sign toward the back, like
% this:
z = exp(1i*theta* -freq); % "Negative" frequency
comet3(radius*real(z), radius*imag(z), theta)
% Same result, because it doesn't matter where we flip the sign in the
% expression. This all might seem trivial, but it's actually rather
% important. In fact, we're going to stick with this label of 'negative
% frequency' moving forward, because the counter-clockwise (positive) vs
% clockwise (negative) motion is gonna prove critical.

% To be complete, we can - of course - also move the negative sign to be in
% front of the angles. Since angles are on the z-dimension, we're kinda 
% heading into 'negative' angles.
z = exp(1i* -theta *freq);
comet3(radius*real(z), radius*imag(z), theta)
% Same clockwise motion as before. If you then insist, you can also flip
% the orientation of the z-axis itself, so that things do move downward
% now.
comet3(radius*real(z), radius*imag(z), -theta)

% In any case, the critical matter is the direction the dot is moving:
% counter-clockwise [ exp( 1i*theta*freq) ] or.. 
% clockwise         [ exp(-1i*theta*freq) ]
close all


% Let's take another look at the 2d unit circle animation.
ANTSD_animatedwaves(1)
% Look at the bottom right panel, it shows the two waves (sine and cosine)
% plotted in the same figure. The cosine starts at 1, whereas the sine
% starts at 0. After one period, they of course are back where they
% started. Now as yourself, 'How much out of phase are they?'. If you had
% to start with the sine function [ sin() ] and had to make a cosine, what
% would you do? Well, when the cos() is at 1, sin() is at 0. When sin() is
% at 1, cos() is at 0. Interesting! If you look at the points on the circle
% where those peaks (at 1) happen, they are 90 degrees apart. They are
% /orthogonal/. This should make sense, cos() gets its value on the real
% line and sin() on the imaginary line.. and we defined those lines earlier
% as.. lines that were orthogonal. We'll come back to the importance of
% this in a moment.





clear all; close all; clc

% And, finally, as with the other periodic functions, we can now 'speed
% up', 'slow down', and 'shift' that function around by introducing a
% multiplier that represents the frequency and introducing added/subtracted
% terms that represent phase, and so on (see sine waves in Chapter 2).
%          Figure 11.2 is another example of multiple individual sine waves
% being summed to create a 'summed sine wave', comparable to the example in
% Chapter 2. The summed sine wave looks much more like an actual EEG signal
% than any of the individual waves. Fourier transform can take that summed
% sine wave and 'de-sum/decompose/deconstruct/demix' that wave from Figure
% 11.2B and give us back (more or less) the individual waves we see in
% Figure 11.2A. Even when there's random noise added to the summed wave, as
% in Figure 11.3. Let's FINALLY see how that works.

%% Fourier Transform
% Yes, it's finally time -- we are now ready for Fourier transforms.
% Specifically, what we'll be doing is /discrete time/ Fourier transform.
% It is discrete time, as opposed to continuous, because we have our time
% domain in discrete steps/time points. If we'd recorded the EEG signal
% continously, rather than with samples every e.g. 100 ms, it wouldn't be
% discrete -- then it would be called continuous-time Fourier transform,
% which is simply 'the' Fourier transform. However, in practice, we're
% unlikely to have access to continuous data and are only able to take
% samples of that continuous signal with a certain sample rate.

% By the way, what is a 'transform'? It's akin to a function, which takes
% an input value and produces an output value. A transform takes a /set/ of
% input values and produces a /set/ of output values. The way a function
% takes you from a number to a number, a transform takes you from one
% function to another function. I mean, okay - this still feels very vague 
% and imprecise, but that definition is sufficient, in case you were
% wondering why it's Fourier /transform/ instead of Fourier /function/.

% Let's look at the code on page 125, where a Fourier transform is shown.
N = 10; % length of sequence
data = rand(1,N); % random numbers
% initialize Fourier coefficients
fourier = zeros(size(data));
time = (0:N-1)/N; % time starts at 0; dividing by N normalizes to 1
% Fourier transform
for fi=1:N % We're testing all frequencies from 0 Hz to N-1 Hz against the data
    % create sine wave (well.. 'sinusoid' but okay)
    sine_wave = exp(-1i*2*pi*(fi-1).*time); % Note the -1i (negative!)
    % compute dot product between sine wave and data
    fourier(fi) = sum(sine_wave .* data);
end
% In the next example, we'll also divide the sums by the number of data
% points as follows (the example in the book doesn't do that yet):
% fourier=fourier/N
% Why is that done? In the last line of the loop, the dot product, we can
% see that each Fourier coefficient is the result of a sum of N data
% points. After all.. that's what the dot product is. This means that, with
% every additional data point, the dot product will keep growing.. and
% growing.. and growing. What we're ultimately interested in, of course, is
% more akin to an 'average' measure of similarity for each data point. And
% so, we divide it by N.

% So, we have 10 data points (defined by N and data) and will be resolving 
% 10 different frequencies (in variable 'fourier') with a sine wave that is
% represented by frequencies from 0 to 9 (i.e., fi-1) unfolding over 1
% second ((0:N-1)/N).

% Why are we subtracting 1 from each fi? fi is defined as numbers from 1:N,
% where N=10. So: 1, 2, 3 ... 10. But, by subtracing 1 each time, we get:
% 0, 1, 2 ... 9. In doing so, the first iteration is with frequency=0,
% which is the "direct current" component (i.e., DC component). The DC
% component is a flat line and thus captures the overal mean offset of the
% signal (i.e., the bias or 'constant' shift up or down of the entire
% signal).

% You might also be wondering, why is it -1i instead of 1i in the
% exponent? Simple answer is that 'terms in the math of the transform will
% end up canceling out' if we don't do it this way. So, for now, keep in
% mind that this is what needs to be done. In fact, when we do a so-called
% /inverse/ Fourier transform, you'll see that we rotate things 
% counter-clockwise again. You might see people state that this negative
% sign is 'mere convention', but there's a mathematical reason for it.

% Now that you know how to do the transform, let's see how you can actually
% view the results. We'll recreate Figure 11.4.
srate=500;
time=-1:1/srate:1;

% Create three sine waves
s1 = sin(2*pi*3*time); % Freq = 3, Amp = 1
s2 = 0.5*sin(2*pi*8*time); % Freq = 8, Amp = 0.5
s3 = s1+s2; % The summed wave

hz = linspace(0,srate/2,floor(length(time)/2)+1);
% With srate=500, we have a Nyquist of (srate/2=) 250. And because we
% want to include the DC component, we go from 0 to the Nyquist frequency 
% in (n_datapoints/2)+1 steps.

% plot the sine waves
figure
for i=1:3
    subplot(2,3,i) % Top row of figure
    
    % plot sine waves, using the eval command (evaluate the string)
    eval([ 'plot(time,s' num2str(i) ')' ]);
    set(gca,'ylim',[-1.6 1.6],'ytick',-1.5:.5:1.5)
    
    % plot power
    subplot(2,3,i+3) % Bottom row of figure
    f  = eval([ 'fft(s' num2str(i) ')/length(time)' ]);
    % Note that eval() takes text and then evaluates the expression.
    % For instance, with i=1, this becomes:
    % f = fft(s1)/length(time);
    % So, it returns the fast Fourier transform for wave s1

    % NOTE:
    % Only plotting the so-called 'positive' frequencies!
    % We are also doubling (some of) the coefficients we got!
    % More on both these choices in a moment.
    f(2:length(hz)) = 2 * f(2:length(hz));
    bar(hz,abs(f(1:length(hz))))
    set(gca,'xlim',[0 11],'xtick',0:10,'ylim',[0 1.2])
end
% As you can see in the bar graphs, the first result tells us that s1 has
% an amplitude of 1 at 3 Hz (which is indeed correct!) and that s2 has an
% amplitude of .5 at 8 Hz (also correct!). The final panel shows us that
% even when the waves are added, the different features can still be
% demixed accurately. Cool!

% If you're a particularly smart cookie, you might now ask yourself "But
% how do I know if this is the 'true' decomposition (if I didn't have that
% information to start with, like with actual data)? After all, each
% decomposed sinusoid can itself also be represented by a summation of
% other sinusoids! And so on and so forth, sinusoids all the way down!
% And actually, why sinusoids?! Why not square waves or triangular waves?!"
% And, well, you'd be right in asking any of those questions. There are an
% infinite number of ways to decompose any given signal; our goal here
% is to find something that is easier to manage, analyze, and interpret
% than the original signal. (If you're interested in why we use sinusoids
% specifically, look up 'sinusoidal fidelity').


% Okay, so when we were plotting amplitudes, we made two important choices.
% First, we only plotted /some/ of the Fourier coefficients (i.e., the
% output of the transform as produced by fft(). In fact, we're plotting
% only the first half of the entire output. Second, the half that we're
% plotting gets multiplied by 2. Why?
%           The frequency spectrum of a the signal is symmetrical around
% the DC (0 Hz): thus, there are positive frequencies and negative
% frequencies. Let's see the FULL power spectrum of the summed signals in
% the figure above. Also, let's NOT double our coefficients yet.
figure
subplot(311)
fullhz = linspace(0,srate,length(time));
bar(fullhz,abs( f )) % Note, no doubling for now!
ylim([0 1.5])
% When we plot the entirety of our fft() output, we see that, towards the
% end of our output, the pattern at the start is mirrored. Importantly,
% it's mirrored around a particular frequency. If you look at the x-axis,
% you can probably guess what frequency that is. Here, let's plot it:
xline(srate/2,"LineWidth",5);
% That's right -- it's the Nyquist rate! Once we get above the Nyquist
% frequency, the power spectrum is mirrored. This is also referred to as
% Hermitian symmetry.

% We can use ffshift() to shift our Nyquist to the center of the spectrum.
subplot(312)
bar(fullhz,fftshift(abs(f)))
% If you hover your mouse over the tops of the 2 shorter graphed bars,
% you'll see that one is plotted at 258 and the other at 242 -- both
% equally far away from 250 (i.e. +/- 8).

% In our first plot, the shorter bar was plotted at +8 Hz, so let's shift
% our center point to be 0 Hz - i.e., the DC - instead.
subplot(313)
mirroredhz = -(srate/2):.5:(srate/2);
bar(mirroredhz,fftshift(abs(f)))
% This should, hopefully, help you understand where labels like positive
% and negative frequencies come from. The negative frequencies are the
% mirrored versions of the positive frequencies. But.. that doesn't
% actually explain what those negative frequencies /are/ (just why they may
% be called that).

% Before we get to that, though, let's go back to the step where we double
% the output of the transform to get the amplitude back. The reason we do
% that is because half of the amplitude of the decomposed signal is
% represented in the positive frequencies and the other half in the
% negative frequencies -- we can see that now in the previous figures. And
% so, when we DON'T plot a double-sided spectrum (i.e., the full spectrum
% with the mirror) and only plot the single-sided spectrum.. we need to
% multiply the coefficients by 2 to get the 'full picture'.
xlim([-10 10])
% You can zoom in on the x-axis in the bottom plot to get a closer look.
% Hover your mouse over the top of the bars again and look at the y-values.

% Back to negative frequencies, what they are and where they come from. 
% Negative frequencies result from us working with discretized signals.
% Because our data isn't continuous, the individually sampled data points
% in our signal can be fit equally well by multiple frequencies.
%           For instance, say you have a 2 Hz signal sampled at 16 Hz. The 
% Nyquist frequency, therefore, is 8 Hz. The Fourier transform might then
% test frequencies from 0 to 16 Hz in 1 Hz steps, testing the similarity to
% the signal using a cosine and inverted sine wave.
%           Recall: cos(angle) + i * sin(angle), which then gets rotated in
% a CLOCKWISE manner, meaning the sine wave starts out heading in a negative
% direction rather than positive. There was a reason for that, if you 
% remember, which is that 'some math terms' cancel out that way and that 
% that's a good thing. Well, get ready for that to be demystified shortly.
%           As the transform is testing the various frequencies, each
% cosine and sine component of the kernel get a 'similarity' score. The
% problem with discretized signals, as mentioned, is that multiple
% frequencies can fit the discretized points equally well. In our example
% of a 2 Hz signal sampled at 16 Hz, a kernel of 2 Hz gets the same
% similarity score as one that is 14 Hz, because we don't know what is
% going on in between the sampled points. The same is true is for 7 and 9
% Hz, and for 6 and 10 Hz, and so on. This is known as 'aliasing'.
%           Okay, still not quite explaining what is 'negative' about any
% of this. Take another look at this animation:
addpath '/Users/hamid/Desktop/Research/+Tutorials/ANTSD'
ANTSD_animatedwaves(1)
% As we move around the cirle in a COUNTER-CLOCKWISE fashion, the sine and
% cosine wave both start out positive. But now imagine moving in the
% CLOCKWISE direction. This doesn't change anything for the cosine, but the
% sine wave is now flipped - starting out negative. The sine wave is now..
% a mirrored image of itself when moving in the counter-clockwise
% direction.
%           When we perform a Fourier transform, you'll have noticed, that
% we end up with complex coefficient. That's the numbers in the form a + bi
% where the 'a' represents the cosine component and the 'bi' represents the
% sine component. Here's the thing: we DID move along the unit circle in
% the clockwise manner, creating a mirrored sine wave. Moving in this
% direction, at a given frequency, is called the 'negative' frequency,
% whereas moving around the circle in a typical, counter-clockwise manner
% is called a 'positive' frequency. It's like we're taking the sine wave,
% at a certain frequency, and then flipping its sign - hence, negative
% frequency. But here's the brilliance in that: to create a sinusoid, we
% need to add the cosine and (now mirrored) sine components together,
% (see, they're summed:  cos(angle) + i * sin(angle) )
% meaning that complex value of the cosine component takes on the form
% a + bi, but the complex value of the sine component takes on the form
% a - bi. And so, when we add them together, a + a + bi - bi, and the bi
% terms cancel out and we're just left with 2a. Badda-bing, badda-boom.
% These so-called complex conjugates cancel out and we're left with a real
% number (the 'a') that gets doubled. We can go through the whole
% procedure, add the components together, let them cancel out, and so on..
% or we can just plot the positive frequencies and double the coefficient
% and get on with our lives. I suggest the latter.

% There's one point left to make!
% !!! WARNING:  Since we are mirroring things around the 0 frequency and 
%               Nyquist frequency, you DON'T double their amplitudes (only
%               double positive frequencies /between/ 0 Hz and Nyquist).

% If some of this is still unclear, I can recommend these videos:
% https://www.youtube.com/watch?v=IIofPiVVC64
% https://www.youtube.com/watch?v=_3-qntJ12q4
% In fact, the whole channel is great.

% If you simply want to move on, all you need to know is that there are
% things called negative frequencies that we're not particularly interested
% in (for our purposes here). So, we (1) show the power spectrum of the
% positive frequencies, and we need to (2) double the absolute value of the
% Fourier coefficient we get to get back the amplitude of the signal.


% Once you have the Fourier coefficients, you can then reconstruct the
% signal (an /inverse/ Fourier transform). We'll take a look at that in a
% moment. But first, you might be wondering what fft(), i.e. the so-called
% Fast Fourier Transform is doing and how it relates to the code we saw on
% page 125.

% For instance, you might be wondering, how does the Fourier transform know
% which frequencies to compare against the signal? There's an infinite
% number of them. That process would go on forever. So, it has to limit
% itself, in order to get the job done /fast/. To do so, the Fourier
% transform limits itself to a total number of frequencies that is equal to
% the length of the number of data points in your data: if you have 1000
% data points, fft() will return 1000 coefficients, where the first entry
% is the coefficient for 0 Hz (the DC) and every subsequent frequency is a
% constant interval above the previous one (generally speaking, 1): e.g., 0
% Hz, 1 Hz, 2 Hz, ... N Hz (where N is the number of data points). To
% comfortably measure all the frequencies present in the signal, you should
% therefore have at least 2 times the number of data points in your data as
% your Nyquist frequency. And indeed, the code on page 125 also tested N
% frequencies, where N equalled the number of data points available. So,
% what happens when we have an explicit set of frequencies we want to test?
srate=500;
time=-1:1/srate:1;

% Create three sine waves
s1 = sin(2*pi*3*time); % Freq = 3, Amp = 1
s2 = 0.5*sin(2*pi*8*time); % Freq = 8, Amp = 0.5
s3 = s1+s2; % The summed wave

%hz = linspace(0,srate/2,floor(length(time)/2)+1);
hz = 0:1:10; % Same signal as before, but now only testing these 11 freqs

for fi=1:length(hz) % Note that it's not 1:n_datapoints, but 1:n_freqs
    sine_wave   = exp(-1i*2*pi*(fi-1).*time);
    %plot(sine_wave)
    fourier(fi) = sum(sine_wave.*s3);
end
fourier=fourier/length(time);

figure
bar(hz,abs(fourier))
ylim([0 2])
% Compare this output to the output of fft(), the bars at 3 Hz and 8 Hz are
% exactly half as big as when we run fft() and double the coefficients. So
% far, so good -- we simply didn't double.
bar(hz,abs(fourier*2))
% When we do double, everything looks the same as running fft() and doubling.
% So, it doesn't matter whether we compare N (n_datapoints) frequencies
% /or/ only a small number of frequencies, right? We can just double the
% output and - as long as our smaller set contains the right frequencies in
% it - everything is fine, right?!

% Well.. not quite. Let's reconstruct the original signal using the same
% shortened set of test frequencies.
reconstructed_data = zeros(1,length(time)); % Reconstructing the signal, so same length as the original
for fi=1:length(hz)
    % scale sine wave by fourier coefficient
    sine_wave = fourier(fi)*exp(1i*2*pi*(fi-1).*time); % Note the counter-clockwise rotation (to reconstruct)!
    % sum sine waves together (take only real part)
    reconstructed_data = reconstructed_data + real(sine_wave);
end
% And now let's plot the original signal and overlay its reconstruction:

figure
plot(s3,'-o')
hold on
plot(reconstructed_data,'r-*')
legend({'original data';'inverse Fourier transform data'})
% Okay.. that doesn't look right. I mean, it's /almost/ there. It's got the
% general shape of the signal, but it's somehow scaled down. Look at the
% y-axis, it's as if the reconstructed signal is exactly half the size of
% the original. In fact...
plot(reconstructed_data * 2,'r-*') % Double the reconstructed
% ..and there it is! There somehow was a factor of 2 missing. You might
% already know where this is going. That factor was hiding, somewhere, in
% the negative frequencies, which we wouldn't capture unless we are running
% comparisons between kernel and signal where the kernels also include
% negative frequencies (i.e., ones above Nyquist). So, if our sample rate
% is 500, our Nyquist is 250, we should probably run at least 501 kernels
% (0 Hz to 500 Hz).

% So let's run it with 501 frequencies and plot the reconstructed signal
% WITHOUT doubling. It should come out matching the original signal right
% away, right?
n_freqs=501;
for fi=1:n_freqs
    sine_wave   = exp(-1i*2*pi*(fi-1).*time);
    %plot(sine_wave)
    fourier(fi) = sum(sine_wave.*s3);
end
fourier=fourier/length(time);

reconstructed_data = zeros(1,length(time)); % Reconstructing the signal, so same length as the original
for fi=1:n_freqs
    % scale sine wave by fourier coefficient
    sine_wave = fourier(fi)*exp(1i*2*pi*(fi-1).*time); % Note the counter-clockwise rotation (to reconstruct)!
    % sum sine waves together (take only real part)
    reconstructed_data = reconstructed_data + real(sine_wave);
end

figure
plot(s3,'-o')
hold on
plot(reconstructed_data,'r-*')
legend({'original data';'inverse Fourier transform data'})
% Now when we plot the reconstructed signal (notice, WITHOUT having to
% double it), we indeed immediately get our original signal back.


% What if we miss some of the negative frequencies?
n_freqs=350;
for fi=1:n_freqs
    sine_wave   = exp(-1i*2*pi*(fi-1).*time);
    %plot(sine_wave)
    fourier(fi) = sum(sine_wave.*s3);
end
fourier=fourier/length(time);
reconstructed_data = zeros(1,length(time)); % Reconstructing the signal, so same length as the original
for fi=1:n_freqs
    % scale sine wave by fourier coefficient
    sine_wave = fourier(fi)*exp(1i*2*pi*(fi-1).*time); % Note the counter-clockwise rotation (to reconstruct)!
    % sum sine waves together (take only real part)
    reconstructed_data = reconstructed_data + real(sine_wave);
end
figure
plot(s3,'-o')
hold on
plot(reconstructed_data,'r-*')
legend({'original data';'inverse Fourier transform data'})
% Again, it's a scaled down version of the original.

% Given what you know about mirroring of the spectrum, and the fact that we
% need to sum together the amplitude of the positive frequency with its 
% corresponding, mirrored negative frequency, you should be able
% to figure which negative frequencies you can and cannot miss. Keep lowering 
% the n_freqs number and see if you can spot where the drop in amplitude for
% the reconstructed signal happens. Can you figure out why it happens where
% it does?
% HINT: our original signal has power at 3 Hz (which is 247 Hz away from
% the Nyquist of 250) and 8 Hz (which is 242 Hz away).
n_freqs=492;
for fi=1:n_freqs
    sine_wave   = exp(-1i*2*pi*(fi-1).*time);
    %plot(sine_wave)
    fourier(fi) = sum(sine_wave.*s3);
end
fourier=fourier/length(time);
reconstructed_data = zeros(1,length(time)); % Reconstructing the signal, so same length as the original
for fi=1:n_freqs
    % scale sine wave by fourier coefficient
    sine_wave = fourier(fi)*exp(1i*2*pi*(fi-1).*time); % Note the counter-clockwise rotation (to reconstruct)!
    % sum sine waves together (take only real part)
    reconstructed_data = reconstructed_data + real(sine_wave);
end
figure
plot(s3,'-o')
hold on
plot(reconstructed_data,'r-*')
legend({'original data';'inverse Fourier transform data'})
% The answer is in the section below:

% The mirrored, negative frequencies of 3 Hz and 8 Hz are 497 Hz and 492 Hz, 
% respectively. So, if you enter those numbers as n_freqs, the for-loop will 
% run from, e.g., 0 to 497-1 (and so, it will miss 497 Hz). At that point,
% you will see that the amplitude of teh figure has come down by 0.5 (half
% of 1). If you then set it to 492, it will miss the 492 Hz frequency and
% the amplitude of the reconstructed signal will come down again, but this
% time by .25 (which is half of the amplitude of for the 8 Hz wave.

% What's the takeaway? It's best to just run the Fourier transform with N
% sinusoids, starting from 0 Hz to n_datapoints-1 Hz.


% We'll do a few more examples and see what else we can do. It's also going
% to get a bit more complicated than what I've outlined above. Sorry!
% This next example, Figure 11.5, is a more complex:
N       = 10;         % length of sequence
data    = randn(1,N); % random numbers
srate   = 200;        % sampling rate in Hz
nyquist = srate/2;    % Nyquist frequency -- the highest frequency you can measure in the data

% initialize Fourier output matrix
fourier = zeros(size(data)); % Equal to the number of data points

% These are the actual frequencies in Hz that will be returned by the
% Fourier transform. The number of unique frequencies we can measure is
% exactly 1/2 of the number of data points in the time series plus 1 (the DC).
% That's because the we need at least half of the data points to account
% for the negative frequencies.
frequencies = linspace(0,nyquist,N/2+1);
time = ((1:N)-1)/N;

% Fourier transform is dot-product between sinusoid and the data
for fi=1:N
    sinusoid   = exp(-1i*2*pi*(fi-1).*time);
    %plot(sinusoid)
    fourier(fi) = sum(sinusoid.*data);
end
fourier=fourier/N; % Note that we're scaling down the dot product this time

figure
subplot(221)
plot(data,'-o')
set(gca,'xlim',[0 N+1])
title('Time domain representation of the data')

subplot(222)
plot3(frequencies,angle(fourier(1:length(frequencies))),abs(fourier(1:length(frequencies))).^2,'-o','linew',3)
grid on
xlabel('Frequency (Hz)')
ylabel('Phase')
zlabel('power')
title('3-D representation of the Fourier transform')
view([20 20])

subplot(223)
bar(frequencies,abs(fourier(1:length(frequencies))).^2)
set(gca,'xlim',[-5 105])
xlabel('Frequency (Hz)')
ylabel('Power')
title('Power spectrum derived from discrete Fourier transform')

subplot(224)
bar(frequencies,angle(fourier(1:length(frequencies))))
set(gca,'xlim',[-5 105])
xlabel('Frequency (Hz)')
ylabel('Phase angle')
set(gca,'ytick',-pi:pi/2:pi)
title('Phase spectrum derived from discrete Fourier transform')
% This shows us that we can recover more than amplitude. Indeed, we can get
% the power (which is, for all intents and purposes here, amplitude
% squared) at each frequency (of interest) and phase.



% We'll bring back this summed sine wave example from Chapter 2 and see if 
% the Fourier transform can return the components from the signal (okay, I
% changed it a little -- but it's like the one we saw in Chapter 2).
tot_dur = 10;
srate = 200;
% I want 10 seconds of data with a sample rate of 200 Hz
n_data = tot_dur*srate; % ..meaning 2000 data points
time = linspace(0,tot_dur,n_data); % Our discretized time vector
% Create signal data and add some noise
data = sin(2*pi*time) + sin(2*pi*time * 3 + (pi/4)) + .7 * sin(2*pi*time * 36 - (pi/2)) + 1.3;
data = data+(randn(size(data))*1.5);
plot(time,data)
% In other words, we're looking for signals coming through with a frequency
% of 1 Hz (the first sine), 3 Hz (the second sine), and 36 Hz (the third
% sine). Note that I changed the bias (i.e., the zero frequency) to 1.3.
% The data look kinda noisy now, so let's see what happens!

frequencies = 0:1:40; % All frequencies we're gonna compare against data
n_freqs = length(frequencies);
fourier = zeros(1,n_freqs);

for fi=1:n_freqs
    sine_wave = exp(-1i*2*pi*(fi-1).*time);
    %plot(sine_wave)
    fourier(fi) = sum(sine_wave.*data);
end
% Recall, each Fourier coefficient is the result of the sum of N data points, 
% so we take the average /per/ data point by dividing by n_data:
fourier=fourier/n_data;

figure
subplot(221)
plot(data,'-')
set(gca,'xlim',[0 n_data+1])
title('Time domain representation of the data')

subplot(222)
plot3(frequencies,angle(fourier),abs(fourier).^2,'-o','linew',3)
grid on
xlabel('Frequency (Hz)')
ylabel('Phase')
zlabel('Power')
title('3-D representation of the Fourier transform')
view([20 20])

subplot(223)
bar(frequencies,abs(fourier).^2)
set(gca,'xlim',[-1 40], 'ylim',[0 2])
xlabel('Frequency (Hz)')
ylabel('Power')
title('Power spectrum derived from discrete Fourier transform')

subplot(224)
bar(frequencies,angle(fourier))
set(gca,'xlim',[0 40])
xlabel('Frequency (Hz)')
ylabel('Phase angle')
set(gca,'ytick',-pi:pi/2:pi)
title('Phase spectrum derived from discrete Fourier transform')

% Things are noisy, but we're seeing some hints of the original signal! We
% see power at the 0 Hz, 1 Hz, 3 Hz, and 36 Hz! The power of the 0 Hz
% component is roughly the square of 1.3 (=1.69). Indeed, we set the 0 Hz 
% (DC component) to be 1.3 in the signal. In my case, I got values
% a little above 0.2 power for 1 Hz and 3 Hz. At 36 Hz, power was ~0.1.
% Your results will differ, because we're introducing random noise.
% Nevertheless, it recovered the relative strengths of the amplitudes of
% the signals with different frequencies. Impressive! The signal is also
% shifted in phase for two of the signals out of three. Did it recover that
% information? It's showing positive shifts in lower frequencies, negative
% shifts in frequencies between 30-40 Hz. But it's not very precise.










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                               CHAPTER 12.                               %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




























%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                               CHAPTER 13.                               %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Exercises

load sampleEEGdata

%% 13.11.1
% Create a family of complex morlet wavelets, ranging in frequencies
% from 2 Hz to 30 Hz in five steps.

n_wav = 5;                      % number of frequencies in family
freqs = linspace(2,30,n_wav);   % family of frequencies:  2     9    16    23    30
srate = EEG.srate;              % data sampling rate in Hz
time  = -2:1/srate:2;           % time vector, from -2 to 2 second in steps of 1/srate
n_cycs = 6;                     % number of cycles

% Make Morlet wavelets
for w = 1:n_wav
    s = n_cycs/(2*pi*freqs(w));  % spread of the gaussian
    sinewave = exp(2*pi*1i*freqs(w).*time);
    gauswindow = exp(-time.^2./(2*s^2));
    wavelet{w} = sinewave .* gauswindow;
end

% Plot a few of the complex wavelets
figure
subplot(321)
plot3(time,real(wavelet{1}),imag(wavelet{1}),'m')
xlabel('Time (ms)'), ylabel('real axis')
title('Wavelet 2 Hz')
view(0,90)

subplot(322)
plot3(time,real(wavelet{1}),imag(wavelet{1}),'m')
xlabel('Time (ms)'), ylabel('real axis')
title('Wavelet 2 Hz')

subplot(312)
plot3(time,real(wavelet{3}),imag(wavelet{3}),'m')
xlabel('Time (ms)'), ylabel('real axis')
title('Wavelet 16 Hz')
view(0,90)

subplot(313)
plot3(time,real(wavelet{5}),imag(wavelet{5}),'m')
xlabel('Time (ms)'), ylabel('real axis')
title('Wavelet 30 Hz')
view(0,90)

%% 13.11.2
% Convolve each wavelet with EEG data from all electrodes and from one trial.

% Extract a given trial, all channels
which_trial = 34;
data = squeeze(EEG.data(:,:,which_trial));

% Convolve with data
l_wav = length(wavelet{1});
halfwaveletsize = ceil(l_wav/2);
n_conv = l_wav + EEG.pnts - 1;

for c = 1:EEG.nbchan
    chan_d = data(c,:);
    for w = 1:n_wav
        fft_w = fft(wavelet{w},n_conv);
        fft_d = fft(chan_d,n_conv);
        s = 6/(2*pi*freqs(w));
        ift   = ifft(fft_d.*fft_w,n_conv)*sqrt(s)/10;
        %wavelet_conv_data{w}(c,:) = real(ift(halfwaveletsize:end-halfwaveletsize+1));
        wavelet_conv_data{w}(c,:) = ift(halfwaveletsize:end-halfwaveletsize+1);
    end
end

%% 13.11.3
% Extract power and phase from the convolutions.
for c = 1:EEG.nbchan
    for w = 1:n_wav
        power{w}(c,:) = abs(wavelet_conv_data{w}(c,:)).^2;
        phase{w}(c,:) = angle(wavelet_conv_data{w}(c,:));
    end
end

%% 13.11.4
% Make topo plots of power and phase at 180 ms at all frequencies.

% From EEG.times, we can see that column 303 represents ~180 ms
[~, idx] = min(abs(EEG.times-180));
EEG.times(idx) % 179.6875

figure
% Power
subplot(251); topoplot(double(squeeze(power{1}(:,idx))),EEG.chanlocs)
subplot(252); topoplot(double(squeeze(power{2}(:,idx))),EEG.chanlocs)
subplot(253); topoplot(double(squeeze(power{3}(:,idx))),EEG.chanlocs)
subplot(254); topoplot(double(squeeze(power{4}(:,idx))),EEG.chanlocs)
subplot(255); topoplot(double(squeeze(power{5}(:,idx))),EEG.chanlocs)
% Phase
subplot(256); topoplot(double(squeeze(phase{1}(:,idx))),EEG.chanlocs)
subplot(257); topoplot(double(squeeze(phase{2}(:,idx))),EEG.chanlocs)
subplot(258); topoplot(double(squeeze(phase{3}(:,idx))),EEG.chanlocs)
subplot(259); topoplot(double(squeeze(phase{4}(:,idx))),EEG.chanlocs)
subplot(2,5,10); topoplot(double(squeeze(phase{5}(:,idx))),EEG.chanlocs)

%% 13.11.5
% Make topo plots of power and phase at 360 ms at all frequencies.

% From EEG.times, we can see that column 303 represents ~180 ms
[~, idx] = min(abs(EEG.times-360));
EEG.times(idx) % 179.6875

figure
% Power
subplot(251); topoplot(double(squeeze(power{1}(:,idx))),EEG.chanlocs)
subplot(252); topoplot(double(squeeze(power{2}(:,idx))),EEG.chanlocs)
subplot(253); topoplot(double(squeeze(power{3}(:,idx))),EEG.chanlocs)
subplot(254); topoplot(double(squeeze(power{4}(:,idx))),EEG.chanlocs)
subplot(255); topoplot(double(squeeze(power{5}(:,idx))),EEG.chanlocs)
% Phase
subplot(256); topoplot(double(squeeze(phase{1}(:,idx))),EEG.chanlocs)
subplot(257); topoplot(double(squeeze(phase{2}(:,idx))),EEG.chanlocs)
subplot(258); topoplot(double(squeeze(phase{3}(:,idx))),EEG.chanlocs)
subplot(259); topoplot(double(squeeze(phase{4}(:,idx))),EEG.chanlocs)
subplot(2,5,10); topoplot(double(squeeze(phase{5}(:,idx))),EEG.chanlocs)

%% 13.11.6

%% 13.11.7


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                               CHAPTER 14.                               %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We're gonna take a look at another method for time-frequency
% decomposition: bandpass filtering and the Hilbert transform (the
% filter-Hilbert method). The result of this method, just like the result
% of the complex wavelet convolution, is called the 'analytic signal'.
% Analytic, in math, refers to the ability to start with e.g. one set of 
% values (like a signal) and arrive at another set of values purely through
% algebraic manipulation (following rules of arithmatic, trig, and so on -
% as opposed to estimating those values by starting with guesses and 
% closing in on a solution or as opposed to a graphical method where you 
% plot things and then look at a function's characteristics).
%           The advantage of filter-Hilbert is that it allows control over
% frequency characteristics of the filter, whereas the shape of a Morlet
% wavelet is always Gaussian. The disadvantage is that filter construction
% , at least in Matlab, requires the Signal Processing Toolbox and that
% bandpass filtering is slower than wavelet convolution.

% Unprocessed EEG signal consists of real values (as in, numbers on the
% real number line). The Hilbert transform extracts a complex signal from a 
% signal with only a real part, similar to wavelet convolution (because a 
% Morlet wavelet contains an imaginary part).
%          We'll perform a Hilbert transform step by step below, following
% the guide on p.176 and Cohen's scripts. Let's bring in our trusted
% example signal from previous chapters:

% Signal
n_sec = 5;     % Number of seconds of data to simulate
srate=1000;    % Sample rate (i.e. data points per seconds)
n_datapoints=srate*n_sec; % Total number of data points
time=linspace(0,n_sec,srate); % Create discretized time line
nyquist=srate/2; % Nyquist frequency
% Simulate signal and visualize
signal = sin(2*pi*time) + sin(2*pi*time * 3 + (pi/4)) + .7 * sin(2*pi*time * 36 - (pi/2)) + 1;
plot(time,signal)

%% Hilbert transform step-by-step (p.176, Chapter 14)
% "First, compute the Fourier transform of a signal.."
fourier = fft(signal);

% "..and create a copy of the Fourier coefficients that have been
% multiplied by the complex operator (i). This turns the M*cos(2*pi*f*t)
% into i*M*cos(2*pi*f*t*)."
f_complex = 1i*fourier;

% Next, identify the positive and negative frequencies. The positive
% frequencies are those between but not including the zero and the
% Nyquist, .."
posF = 2:1:(nyquist-1); % Note we're starting from 2
% "..and the negative frequencies are those above the Nyquist frequency
% (throughout the Hilbert transform, the zero and Nyquist frequencies are
% left untouched)."
negF = (nyquist+1):1:srate;

% "The next step is to convert the iMcos(2πft) to iMsin(2πft). Remember
% that cosine and sine are related to each other by one-quarter cycle;
% thus, to convert a cosine to sine, you rotate the positive-frequency
% coefficients one-quarter cycle counter-clockwise in complex space (-90°
% or -π/2) (think about the complex plane: rotating from the positive real
% axis to the positive imaginary axis turns a cosine into a sine). To
% convert a cosine to sine in negative frequencies, you rotate the
% negative-frequency coefficients one-quarter cycle clockwise (90° or π/2).

% (note 1: this works by computing the iMsin(2πft) component, i.e., the phase quadrature)
% (note 2: positive frequencies are rotated counter-clockwise; negative frequencies are rotated clockwise)
fourier(posF) = fourier(posF) + -1i*f_complex(posF);
fourier(negF) = fourier(negF) +  1i*f_complex(negF);
% The next two lines are an alternative and slightly faster method. 
% The book explains why this is equivalent to the previous two lines.
% fourier(posF) = fourier(posF)*2;
% fourier(negF) = fourier(negF)*0;

% "The final step is to take the inverse Fourier transform of the modulated
% Fourier coefficients. The results is the analytic signal, which can be
% used in the same way that you use the result of complex Morlet wavelet
% convolution."

% Take inverse FFT
MyFirstHilbert = ifft(fourier);

% Using hilbert() in Matlab
MatlabHilbert = hilbert(signal);

% Plot results
figure
subplot(211)
plot(abs(MatlabHilbert),'b-')
hold on
plot(abs(MyFirstHilbert),'ro')
legend({'Matlab Hilbert';'My Hilbert'})
title('Magnitude of Hilbert transform')

subplot(212)
plot(angle(MatlabHilbert),'b-')
hold on
plot(angle(MyFirstHilbert),'ro')
legend({'Matlab Hilbert';'My Hilbert'})
title('Phase of Hilbert transform')

% As you can see, this results in the same analytic signal.

% The hilbert() function can be executed on a matrix, to perform the
% transform on many trials or many electrodes simultaneously. You just have
% to be sure that the matrix is formatted properly: the time dimension
% needs to be in the rows (i.e., the first row is time point 1, the
% second row is time point 2, third row is time point 3, and so on),
% meaning that the columns represent trials or electrodes. For example, in
% our sampleEEGdata, we have 64 electrodes and 640 time points per trial.
% So, your final matrix should be a 640x64 matrix. The output of hilbert()
% is going to be transposed, meaning that it returns a 64x640 matrix, in
% our example case. And so, time is in the first dimension when going into
% hilbert() and the second dimension coming out.


%% Bandpass filtering before applying a Hilbert transform

% Although bandpass filtering isn't necessary, it's recommended to do so
% prior to applying the Hilbert transform. All frequencies in the data
% contribute to the final analytic signal following a Hilbert transform,
% meaning that the resultant signal can be difficult to interpret. So, by
% bandpass filtering first, the analytic signal can be interpreted in a
% frequency-band-specific manner. And that leads to the main advantage over
% complex wavelet convolution: we can control the characteristics of our
% filter with Hilbert transforms (whereas it's always a Gaussian with
% wavelets). The remainder of the chapter focuses on filters.

% Filters can be classified as finite impulse response (FIR) or infinite
% impulse response (IIR). For filter-Hilbert methods to decompose
% time-frequency features of EEG data, FIRs are preferred, because they are
% more stable/less likely to introduce nonlinear phase distortion.

% Filters come in four forms: bandpass, band-stop, high-pass, and low-pass.
% Respectively, they are designed to 'let certain frequencies through',
% 'stop certain frequencies', 'let frequencies higher than some threshold
% through', or 'let frequencies below a certain threshold through'. Thus, a
% Morlet wavelet is, for instance, essentially a bandpass filter.

% Broadly speaking, creating a filter kernel requires specifying the
% desired frequencies to keep or cut, which then defines the overall shape 
% of the filter. Im Matlab (and the Signal Processing Toolbox
% specifically), there are a few functions to construct a filter:
% firls, fir1, fir2, firrcos, guassfir, and firpm. You can even go about
% constructing your own from scratch, if you want.

% You should always check your filters before applying them. The goodness
% of a filter can be quantified as the similarity between the actual
% frequency characteristics of the filter and the characteristics you
% specified (obviously - the difference between what you wanted and what
% you ended up getting). This can be measured as the sum of squared errors
% between the ideal/specified filter and the actual/resultant filter. The
% SSE should be very close to zero. If it's above one, don't use it.

% An important note about filters is that they are 'causal', in the sense
% that they are impacted by preceding values. So, filters are applied
% 'directionally': forwards through the data, or backwards. But this
% results in a problem where the underlying signal gets pushed out of
% phase. Fortunately, the solution is simply: always apply the filter
% twice. Once forward (which pushes things out of phase) and once backwards
% (which pushes things back), or vice versa. The function filtfilt() from
% the Signal Processing Toolbox takes care of this for us (as implied by
% the function name). If the toolbox is unavailable, apply Matlab's
% built-in function filter() twice, once in each direction. Note that
% wavelets are noncausal, as they are impacted by preceding as well as
% upcoming values. So they don't have this problem.

% Despite all the knowledge we have nowadays about how filters affects
% signals, the issue of filtering remains debated - especially when it
% comes to EEG data. For further info, take a look at Luck (2005, chapter
% 5), the references in section 9.2 of the book, and also Acunzo,
% Mackenzie, & Van Rossum 2012; Rousselet, 2012).


%% Exercises

load sampleEEGdata
nyquist = EEG.srate/2;
lower_freq_bound = 4;
bound_multiplier = 3;
filter_order = round(bound_multiplier*(EEG.srate/lower_freq_bound));

%% 14.13.1
% Pick two frequencies (e.g., 5 Hz and 25 Hz) and one electrode
freqs = [5 25]; n_freqs = length(freqs);
elec = 5;

% EEG data
data = squeeze(EEG.data(elec,:,:))'; % 99 trials x 640 time points from chosen electrode
dims_data = size(data);
n_trials = dims_data(1);

% Perform complex Morlet wavelet convolution and filter-Hilbert using those
% two frequencies as the peak/center frequencies for all trials.
for fi = 1:n_freqs
        
    % Parameters
    center_freq = freqs(fi); % Chosen center frequency
    filter_spread = 3;       % +/- Hz around center
    wavelet_spread = 5;      % Number of cycles

    % Make filter
    transition_width = 0.2;
    ffrequencies   = [ 0 (1-transition_width)*(center_freq-filter_spread) (center_freq-filter_spread) (center_freq+filter_spread) (1+transition_width)*(center_freq+filter_spread) nyquist ]/nyquist;
    idealresponse  = [ 0 0 1 1 0 0 ];
    filtweights_LS = zscore(firls(filter_order,ffrequencies,idealresponse));
    filtweights_1  = zscore(fir1(filter_order,[center_freq-filter_spread center_freq+filter_spread]./nyquist));

    % Make wavelet
    center_wavelettime = filter_order/2;
    time = (-center_wavelettime/EEG.srate):(1/EEG.srate):(center_wavelettime/EEG.srate);
    s = wavelet_spread/(2*pi*center_freq);
    sinewave = exp(2*pi*1i*center_freq.*time);
    gauswindow = exp(-time.^2./(2*s^2));
    wavelet{fi} = zscore(sinewave .* gauswindow);

    % Power spectrum (normalized)
    fft_LS = abs(fft(filtweights_LS)); filtkern_LS{fi} = fft_LS./max(fft_LS);
    fft_1  = abs(fft(filtweights_1));  filtkern_1{fi}  = fft_1./max(fft_1);

end

% Plot wavelet and filter weights
figure; hold on
for fi = 1:length(freqs)
    subplot(2,1,fi)
    plot(real(wavelet{fi})./max(real(wavelet{fi})),'b'); hold on
    plot(filtweights_LS./max(filtweights_1),'r'); hold on
    plot(filtweights_1./max(filtweights_1),'g'); hold on
    legend({'wavelet';'firls';'fir1'})
    set(gca,'xlim',[0 filter_order+1],'ylim',[-2 2])
    xlabel('Time')
end

% Filter kernels
figure; hold on
subplot(411); plot(filtkern_LS{1},'r'); hold on
subplot(412); plot(filtkern_LS{2},'r'); hold on
subplot(413); plot(filtkern_1{1},'g'); hold on
subplot(414); plot(filtkern_1{2},'g'); hold off

% Filter
filtered_firLS = zeros(n_trials,EEG.pnts,n_freqs);
filtered_fir1 = zeros(n_trials,EEG.pnts,n_freqs);
for k = 1:n_freqs
    for t = 1:n_trials
        filtered_firLS(t,:,k) = filtfilt(filtkern_LS{k},1,double(data(t,:)));
        filtered_fir1(t,:,k) = filtfilt(filtkern_1{k},1,double(data(t,:)));
    end
end

% Plot single and averaged trials for each filter
trial2plot = 10; ylims = [-4000 4000];

% 5 Hz
figure
subplot(311); plot(EEG.times,data(trial2plot,:,1)); hold on;
subplot(312); plot(EEG.times,filtered_firLS(trial2plot,:,1)); hold on; ylim(ylims)
subplot(313); plot(EEG.times,filtered_fir1(trial2plot,:,1)); ylim(ylims)

% 25 Hz
figure
subplot(311); plot(EEG.times,data(trial2plot,:)); hold on;
subplot(312); plot(EEG.times,filtered_firLS(trial2plot,:,2)); hold on; ylim(ylims)
subplot(313); plot(EEG.times,filtered_fir1(trial2plot,:,2)); ylim(ylims)

% Average across trials
for k = 1:n_freqs
    avg_filtered_firLS{k} = mean(filtered_firLS(:,:,k));
    avg_filtered_fir1{k} = mean(filtered_fir1(:,:,k));
end

% 5 Hz
figure
subplot(211); plot(EEG.times,avg_filtered_firLS{1}); hold on; ylim([-350 350])
subplot(212); plot(EEG.times,avg_filtered_fir1{1}); hold on; ylim([-350 350])

% 25 Hz
figure
subplot(211); plot(EEG.times,avg_filtered_firLS{2}); hold on; ylim([-350 350])
subplot(212); plot(EEG.times,avg_filtered_fir1{2}); hold on; ylim([-350 350])

% Convolution
conv_wavelet = zeros(n_trials,EEG.pnts,n_freqs);
conv_firLS = zeros(n_trials,EEG.pnts,n_freqs);
conv_fir1 = zeros(n_trials,EEG.pnts,n_freqs);
for k = 1:n_freqs
    for t = 1:n_trials
        conv_wavelet(t,:,k) = conv(data(t,:),wavelet{k},'same');
    end
    conv_firLS(:,:,k) = hilbert(filtered_firLS(:,:,k)')';
    conv_fir1(:,:,k) = hilbert(filtered_fir1(:,:,k)')';
end

% Plot the resulting power and bandpass-filtered signal (that is, the real
% component of the analytic signal) from each method.
figure
subplot(311); plot(EEG.times,real(conv_wavelet(trial2plot,:,1))); hold on
subplot(312); plot(EEG.times,real(conv_firLS(trial2plot,:,1))); hold on
subplot(313); plot(EEG.times,real(conv_fir1(trial2plot,:,1))); hold on






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                        %
%%                               CHAPTER 15.                              %
%%                                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% There is an alternative to the FFT, which uses a large chunk of data, that
% we've discussed: the short-time FFT (stFFT). The stFFT helps address 2
% limitations of the Fourier transform identified in Chapter 11: it (1)
% obscures time-varying changes in frequency structure and (2) assumes data
% are stationary. The idea behind stFFT is simply: rather than using an
% entire time series, we use brief segments of the data to extract the
% frequency structure.

% Step by step, the stFFT works as follows: (1) take a small segment of the
% data (a few hundred milliseconds), (2) apply a windowing taper to
% minimize edge artifacts and detrend the data, (3) Fourier transform the
% detrended and tapered segment. The stFFT, just like any other Fourier
% transform, returns as many frequencies between DC and Nyquist as there
% are time points in the segment. The resulting power spectrum of the 
% segment is then considered to be the power spectrum at the point in time
% which corresponds to the center time point of the segment.

% The length of the chosen segment depends on the frequencies you're
% interested in: longer time segments are necessary to accurately decompose
% lower frequencies. Additionally, to boost signal-to-noise ratio (SNR),
% one can average power of neighboring frequencies together or apply a
% distance-weighted measure (like a Gaussian).
%           Longer time segments, of course, provide better frequency
% precision and resolution. But, shorter time segments allow for more
% temporal precision. One way to manage this tradeoff nicely is to not use
% the same segment length for all frequencies, but instead to shorten the
% segment with each subsequent higher frequency sinusoid being tested in
% the FFT.

% Each successive time segment can have some overlap with the neighboring
% ones. A 50-90% overlap per segment is typical/considered acceptable.

%% Exercises

load sampleEEGdata

% 15.6.1
% Compute the short-term FFT at each electrode and make topographical maps
% of theta-band (around 6 Hz) power and alpha-band (around 10 Hz) power at
% 150 ms and 700 ms.

do.average = 1;
t = 89; % Trial to use (solutions appear to be trial 89 or 69)
stfft_output = zeros(EEG.nbchan,4); % Initialize final array (4 columns: 150ms-6Hz, 150ms-10Hz, 700ms-6Hz, 700ms-10Hz)
counter = 0; 
for ttime = [150 700] % target time
    for tfreq = [6 10] % target frequency
        warning('off','all')

        % Counter for output array column index
        counter = counter + 1;

        for elec = 1:EEG.nbchan

            % Recall, lower frequencies require longer segments, where higher
            % frequencies allow for shorter segments.
            % For example, if we want to estimate 6 Hz power, the window should at least be (1000/6
            % =) 166.7 ms long and ideally 2-3 times that (x3 = 500 ms)    
            timewin = round((1e3/tfreq)*3,-2); % Round to nearest hundred

            timewinidx = round(timewin/(1e3/EEG.srate)); % So I need 128 data points to get a 500 ms window
            [~,targettimeidx] = min(abs(EEG.times-ttime)); % The time point 150 ms is found at index 295
            starttime = targettimeidx-(timewinidx/2); % If we put the target time in the center of the window, this is where the window starts
            endtime = targettimeidx+(timewinidx/2); % This is where the window ends
            % This should now constitute a ~500 ms window, with our target time of 150
            % ms post stimulus onset at the center of that window

            % Taper
            taper = .5*(1-cos(2*pi*(0:timewinidx-1)/(timewinidx-1))); % Hann
            %taper = .54 - .46*cos(2*pi*(0:timewinidx-1)/(timewinidx-1)); % Hamming
            %taper = exp(-.5*(2.5*(-timewinidx/2:timewinidx/2-1)/(timewinidx/2)).^2); % Gauss

            % Detrend the chosen trial at current electrode
            d = detrend(EEG.data(elec,:,t));

            % Tapered, detrended EEG data
            dtap = d(starttime:(endtime-1)).*taper;

            % Plot
            %timevec = floor(EEG.times(starttime)):1e3/EEG.srate:ceil(EEG.times(endtime-1));
            %plot(timevec,d(starttime:(endtime-1))); xlim([floor(EEG.times(starttime)) ceil(EEG.times(endtime-1))])
            %hold on; plot(timevec,dtap,'r'); xline(ttime,"LineWidth",3); title('One short-time window of data, windowed')

            % Run FFT on the short-time dtap
            dfft = fft(dtap)/timewinidx;

            % Frequencies returned by fft()
            f = linspace(0,EEG.srate/2,floor(length(taper)/2)+1);
            [~,targetfreqidx] = min(abs(f-tfreq)); % The target frequency is at f(4)

            % The entire power spectrum
            %plot(f,abs(dfft(1:floor(length(taper)/2)+1)).^2,'.-');
            %title('power spectrum from that time window'); xline(tfreq)

            % Save power result for the current electrode
            if (do.average) % Average results with neighboring frequencies (can boost SNR)?
                stfft_output(elec,counter) = mean(abs(dfft(targetfreqidx-1:targetfreqidx+1)).^2);
            else
                stfft_output(elec,counter) = abs(dfft(targetfreqidx)).^2;
            end
        end
    end
end
warning('on','all')

% Topoplots
figure
subplot(221)
cbar_range = [-5 3]; % Range of our colorbar
topoplot(stfft_output(:,1),EEG.chanlocs); colorbar; caxis(cbar_range); title("150 ms, 6 Hz")

subplot(222)
cbar_range = [-3 3]; % Range of our colorbar
topoplot(stfft_output(:,2),EEG.chanlocs); colorbar; caxis(cbar_range); title("150 ms, 10 Hz")

subplot(223)
cbar_range = [-8 5]; % Range of our colorbar
topoplot(stfft_output(:,3),EEG.chanlocs); colorbar; caxis(cbar_range); title("700 ms, 6 Hz")

subplot(224)
cbar_range = [-5 3]; % Range of our colorbar
topoplot(stfft_output(:,4),EEG.chanlocs); colorbar; caxis(cbar_range); title("700 ms, 10 Hz")


channel = 1; % Which of the 64 channels we'll plot from
n.trials = 6; % Number of trials to plot

% Plot a few of the individual trials
figure
title('The first n trials')
for t = 1:n.trials
    subplot(n.trials,1,t)
    plot(EEG.times,EEG.data(channel,:,t))
    xline(0,"LineWidth",2);
    ylabel('µV') 
end
xlabel('Time (ms)') 

% These are the first 6 trials in our data on channel 1 (unless you changed
% the above parameters). Now we can average them together and see what
% happens when we do so.

% To make the code more accessible, we can first extract the data for a
% single channel, as follows:
data = squeeze(EEG.data(channel,:,:))';
% Squeeze() drops the unused 3rd dimension, but the trials then end up in
% columns with data points in rows. We flip this back by using Matlab's
% transpose function which is the ' at the end. Chan_data now contains 99 
% rows (trials) and 640 columns (time points).

% Matlab performs averaging column-wise (unless there's only one row, or 
% unless you specify a different dimension explicitly). So, since all our 
% time points are along the columns, we are now averaging all 99 values 
% (one for each trial) that were recorded at time point -1000 ms. Then, we 
% average all 99 values representing -996.0938 ms, and so on. And.. that's 
% it! That's our ERP, the electrical potential related to our event. 
% Easy-peasy, data squeezy.
ERP = mean(data);
figure
plot(EEG.times,ERP); xline(0,"LineWidth",2);
title('The ERP'); xlabel('Time (ms)'); ylabel('Signal (µV)') 


% Same stFFT, but instead of a single trial, compute it on the ERP
do.average = 0;
stfft_output = zeros(EEG.nbchan,4); % Initialize final array (4 columns: 150ms-6Hz, 150ms-10Hz, 700ms-6Hz, 700ms-10Hz)
counter = 0; 
for ttime = [150 700] % target time
    for tfreq = [6 10] % target frequency
        warning('off','all')

        % Counter for output array column index
        counter = counter + 1;

        for elec = 1:EEG.nbchan

            % Recall, lower frequencies require longer segments, where higher
            % frequencies allow for shorter segments.
            % For example, if we want to estimate 6 Hz power, the window should at least be (1000/6
            % =) 166.7 ms long and ideally 2-3 times that (x3 = 500 ms)    
            timewin = round((1e3/tfreq)*3,-2); % Round to nearest hundred

            timewinidx = round(timewin/(1e3/EEG.srate)); % So I need 128 data points to get a 500 ms window
            [~,targettimeidx] = min(abs(EEG.times-ttime)); % The time point 150 ms is found at index 295
            starttime = targettimeidx-(timewinidx/2); % If we put the target time in the center of the window, this is where the window starts
            endtime = targettimeidx+(timewinidx/2); % This is where the window ends
            % This should now constitute a ~500 ms window, with our target time of 150
            % ms post stimulus onset at the center of that window

            % Taper
            taper = .5*(1-cos(2*pi*(0:timewinidx-1)/(timewinidx-1))); % Hann
            %taper = .54 - .46*cos(2*pi*(0:timewinidx-1)/(timewinidx-1)); % Hamming
            %taper = exp(-.5*(2.5*(-timewinidx/2:timewinidx/2-1)/(timewinidx/2)).^2); % Gauss

            % Detrend the ERP at current electrode
            ERP = mean(squeeze(EEG.data(elec,:,:))');
            d = detrend(ERP);

            % Tapered, detrended EEG data
            dtap = d(starttime:(endtime-1)).*taper;

            % Plot
            %timevec = floor(EEG.times(starttime)):1e3/EEG.srate:ceil(EEG.times(endtime-1));
            %plot(timevec,d(starttime:(endtime-1))); xlim([floor(EEG.times(starttime)) ceil(EEG.times(endtime-1))])
            %hold on; plot(timevec,dtap,'r'); xline(ttime,"LineWidth",3); title('One short-time window of data, windowed')

            % Run FFT on the short-time dtap
            dfft = fft(dtap)/timewinidx;

            % Frequencies returned by fft()
            f = linspace(0,EEG.srate/2,floor(length(taper)/2)+1);
            [~,targetfreqidx] = min(abs(f-tfreq)); % The target frequency is at f(4)

            % The entire power spectrum
            %plot(f,abs(dfft(1:floor(length(taper)/2)+1)).^2,'.-');
            %title('power spectrum from that time window'); xline(tfreq)

            % Save power result for the current electrode
            if (do.average) % Average results with neighboring frequencies (can boost SNR)?
                stfft_output(elec,counter) = mean(abs(dfft(targetfreqidx-1:targetfreqidx+1)).^2);
            else
                stfft_output(elec,counter) = abs(dfft(targetfreqidx)).^2;
            end
        end
    end
end

% Topoplots
figure
subplot(221)
cbar_range = [-1 .5]; % Range of our colorbar
topoplot(stfft_output(:,1),EEG.chanlocs); colorbar; caxis(cbar_range); title("150 ms, 6 Hz")

subplot(222)
cbar_range = [-1 1]; % Range of our colorbar
topoplot(stfft_output(:,2),EEG.chanlocs); colorbar; caxis(cbar_range); title("150 ms, 10 Hz")

subplot(223)
cbar_range = [-.08 .05]; % Range of our colorbar
topoplot(stfft_output(:,3),EEG.chanlocs); colorbar; caxis(cbar_range); title("700 ms, 6 Hz")

subplot(224)
cbar_range = [-.08 .05]; % Range of our colorbar
topoplot(stfft_output(:,4),EEG.chanlocs); colorbar; caxis(cbar_range); title("700 ms, 10 Hz")

clear all; close all; clc

% 15.6.2
% Select one electrode and one frequency and compute power over time at
% that electrode using complex wavelet convolution, filter-Hilbert, and the
% short-time FFT.
load sampleEEGdata
elec = 1; tfreq = 12; % chosen electrode and target frequency
ERP = mean(squeeze(EEG.data(elec,:,:))');
d = detrend(ERP);
nyquist = EEG.srate/2;
do.average=0;

% Create wavelet
srate = EEG.srate;              % data sampling rate in Hz
time  = -1:1/srate:1;           % time vector, from -2 to 2 second in steps of 1/srate
n_cycs = 6;                     % number of cycles
s = n_cycs/(2*pi*tfreq);         % spread of the gaussian
sinusoid = exp(2*pi*1i*tfreq.*time);
gauswindow = exp(-time.^2./(2*s^2));
wavelet = sinusoid .* gauswindow;
%plot(time,real(wavelet))

% Create filter
spread = 4;
transition_width = 0.2;
ffrequencies   = [ 0 (1-transition_width)*(tfreq-spread) (tfreq-spread) (tfreq+spread) (1+transition_width)*(tfreq+spread) nyquist ]/nyquist;
idealresponse  = [ 0 0 1 1 0 0 ];
lower_freq_bound = 10;
bound_multiplier = 3;
filter_order = round(bound_multiplier*(EEG.srate/lower_freq_bound));

% We can't use /all/ time points available to us, because we assign the
% power to the center time point in our window. So, some portion of the
% available data on either end will be available as part of the 'data
% segment' but not as 'assignable'. Trial data run from -1000ms to +1500ms.
% We'll run it from -200ms to +800, leaving plenty of room for the first
% and last time segment.
[~,firstbin] = min(abs(EEG.times--200)); [~,lastbin] = min(abs(EEG.times-800));
power.wavelet = zeros(length(firstbin:lastbin),1);
power.hilbert = zeros(length(firstbin:lastbin),1);
power.stfft = zeros(length(firstbin:lastbin),1);
count = 0;
for tpoint = firstbin:lastbin  % At each time point still available
            
            % Update counter
            count = count + 1;

            % Segment 
            timewin = round((1e3/tfreq)*3,-2); % Segment size in ms
            timewinidx = round(timewin/(1e3/EEG.srate)); % Number of data points for the segment
            starttime = tpoint-floor(timewinidx/2); % Start of segment
            endtime = tpoint+ceil(timewinidx/2); % End of segment

            % Taper
            taper = .5*(1-cos(2*pi*(0:timewinidx-1)/(timewinidx-1))); % Hann
            %taper = .54 - .46*cos(2*pi*(0:timewinidx-1)/(timewinidx-1)); % Hamming
            %taper = exp(-.5*(2.5*(-timewinidx/2:timewinidx/2-1)/(timewinidx/2)).^2); % Gauss

            % Frequency to be extracted
            f = linspace(0,EEG.srate/2,length(taper)); % All freqs
            [~,fidx] = min(abs(f-tfreq)); % Index of freq of interest

            % Detrend the ERP at chosen electrode
            d = detrend(ERP);

            % Tapered, detrended EEG data
            dtap = d(starttime:(endtime-1)).*taper;
            

            %%%%% %%%%% Wav %%%%% %%%%%
            halfwaveletsize = ceil(length(wavelet)/2);
            n_conv = length(wavelet) + timewinidx - 1;
            fft_w = fft(wavelet,n_conv);
            fft_d = fft(ERP(starttime:(endtime-1)),n_conv);
            ift   = ifft(fft_d.*fft_w,n_conv)*sqrt(s)/10;
            wavelet_conv_data = ift(halfwaveletsize:end-halfwaveletsize+1);
            powerspectrum = abs(wavelet_conv_data).^2;
            power.wavelet(count) = powerspectrum(fidx);
            %%%%% %%%%% Wav %%%%% %%%%%


            %%%%% %%%%% fHilb %%%%% %%%%%
            filterweights  = firls(filter_order,ffrequencies,idealresponse);
            filter_result = filtfilt(filterweights,1,double(ERP(starttime-80:(endtime-1+80))));
            convol_result = conv(ERP(starttime:(endtime-1)),filterweights,'same'); % could also use ifft(fft(dtap...
            powerspectrum = abs(hilbert(convol_result)).^2;
            power.hilbert(count) = powerspectrum(fidx);
            %%%%% %%%%% fHilb %%%%% %%%%%


            %%%%% %%%%% stFFT %%%%% %%%%%
            dfft = fft(dtap);
            power.stfft(count) = abs(dfft(fidx)).^2;
            %%%%% %%%%% stFFT %%%%% %%%%%

end

% Plot the results of these three time-frequency decomposition methods.
tvec = -200:1000/EEG.srate:800; % Time vector
plot(tvec,normalize(power.stfft), 'Color', 'blue'); hold on
plot(tvec,normalize(power.wavelet), 'Color', 'green'); hold on
plot(tvec,normalize(power.hilbert), 'Color', 'red');
xline(0,'LineWidth',2)
legend('stFFT','wavelet','filter-Hilbert')

% End of Chapter 15
clear all; close all; clc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                        %
%%                               CHAPTER 16.                              %
%%                                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% As we saw in the previous chapter, one way to manage risks of edge
% artifacts is to apply tapers that roll down the signal (fully or mostly)
% to zero. Another thing one can do is use /multitapers/. They are
% especially useful in scenarios where boosting signal-to-noise is
% critical, like when looking at higher-frequency activity (> 30 Hz, say).

% The procedure is simple: instead of using one taper, we multiply copies
% of the time segment by multiple tapers and take their Fourier transform.
% The resulting power spectra are then averaged together. That's it!

% The multitapers are called 'discrete prolate spheroidal sequences' and,
% at least in Matlab, can be obtained with the function dpss() in the
% Signal Processing Toolbox. They are also often called Slepian tapers
% (Slepian, 1978). The tapers are orthogonal to each other and have
% different frequency characteristics.



%% Exercises

load sampleEEGdata

% 16.6.1
% Pick one electrode and compute a time-frequency map of power using both
% multitaper method 

elec = strcmpi('P7',{EEG.chanlocs.labels});
data = squeeze(EEG.data(elec,:,:));
smoothing = 3;
timewin = 500; % in ms
timewinidx = round(timewin/(1e3/EEG.srate));
[~,firstbin] = min(abs(EEG.times--200)); [~,lastbin] = min(abs(EEG.times-1000));
out.multitaper = zeros(EEG.trials,length(firstbin:lastbin),65);

% Time-frequency plot
for t = 1:EEG.trials
    counter = 0;
    for ttime = firstbin:lastbin
        counter = counter + 1;

        % Segment
        [~,targettimeidx] = min(abs(EEG.times-ttime));
        starttime = targettimeidx-(timewinidx/2);
        endtime = targettimeidx+(timewinidx/2);
        segment = data(starttime:(endtime-1),t);

        % Detrend
        d = detrend(segment);

        %%%%% %%%%% Multitaper %%%%% %%%%%
        tapers = dpss(timewinidx,smoothing); % note that in practice, you'll want to set the temporal resolution to be a function of frequency
        f = linspace(0,EEG.srate/2,floor(timewinidx/2)+1); % returned frequency info
        taperpow = zeros(floor(timewinidx/2)+1,1); % initialize power vector (over tapers)

        % loop through tapers
        for tapi = 1:size(tapers,2)-1 % Recall, we don't really need the last taper

            % window and taper data, and get power spectrum
            dtap      = d.*tapers(:,tapi);
            pow       = fft(dtap,timewinidx)/timewinidx;
            pow       = pow(1:floor(timewinidx/2)+1,:);
            taperpow  = taperpow + mean(pow.*conj(pow),2);
        end

        % Save result
        out.multitaper(t,counter,:) = taperpow/tapi;
        %%%%% %%%%% Multitaper %%%%% %%%%%

    end
end
tf.multitaper = squeeze(mean(out.multitaper,1))';

% plot full TF map
figure
tvec=linspace(-200,1000,length(firstbin:lastbin));
contourf(tvec,f,tf.multitaper,'linecolor','none')
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title([ 'Power via multitaper from channel P7' ])
xline(0,'LineWidth',3)
















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                               NOTES                                     %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





% Plot the star using points on a small circle and large circle,
% surrounding the defined star_center.
% sind() and cosd() take the angles in degrees, rather than radians
% Inner and outer points, in degrees
inn_p = [50 130 200 270 340]; out_p = [10 90 170 260 280];
inn_ycoord = sind(inn_p); inn_xcoord = cosd(inn_p);
out_ycoord = sind(out_p); out_xcoord = cosd(out_p);
% Point adjustment and star magnitude
adj = .42; mag = 300;
% Center coordinate of star on the image
star_center = 100;
h=drawpolygon('Position',[mag*star_center+out_xcoord(1) mag*star_center+out_ycoord(1); ...
                        mag*star_center+inn_xcoord(1)-adj mag*star_center+inn_ycoord(1)-adj; ...
                        mag*star_center+out_xcoord(2) mag*star_center+out_ycoord(2); ...
                        mag*star_center+inn_xcoord(2)+adj mag*star_center+inn_ycoord(2)-adj; ...
                        mag*star_center+out_xcoord(3) mag*star_center+out_ycoord(3); ...
                        mag*star_center+inn_xcoord(3)+adj mag*star_center+inn_ycoord(3);...
                        mag*star_center+out_xcoord(4)-adj mag*star_center+out_ycoord(4);...
                        mag*star_center+inn_xcoord(4) mag*star_center+inn_ycoord(4)+adj;...
                        mag*star_center+out_xcoord(5)+adj mag*star_center+out_ycoord(5);...
                        mag*star_center+inn_xcoord(5)-adj mag*star_center+inn_ycoord(5)]);






hold on
% Plot the star using points on a small circle and large circle,
% surrounding the defined star_center.
% sind() and cosd() take the angles in degrees, rather than radians
% Inner and outer points, in degrees
inn_p = [50 130 200 270 340]; out_p = [10 90 170 260 280];
inn_ycoord = sind(inn_p); inn_xcoord = cosd(inn_p);
out_ycoord = sind(out_p); out_xcoord = cosd(out_p);
% Point adjustment and star magnitude
adj = .42; mag = 10;
% Center coordinate of star on the image
star_center = 100;
plot(polyshape([mag*star_center+out_xcoord(1) mag*star_center+inn_xcoord(1)-adj mag*star_center+out_xcoord(2) mag*star_center+inn_xcoord(2)+adj mag*star_center+out_xcoord(3) ...
                mag*star_center+inn_xcoord(3)+adj mag*star_center+out_xcoord(4)-adj mag*star_center+inn_xcoord(4) mag*star_center+out_xcoord(5)+adj mag*star_center+inn_xcoord(5)-adj], ...
               [mag*star_center+out_ycoord(1) mag*star_center+inn_ycoord(1)-adj mag*star_center+out_ycoord(2) mag*star_center+inn_ycoord(2)-adj mag*star_center+out_ycoord(3) ...
                mag*star_center+inn_ycoord(3) mag*star_center+out_ycoord(4) mag*star_center+inn_ycoord(4)+adj mag*star_center+out_ycoord(5) mag*star_center+inn_ycoord(5)]))
hold off








