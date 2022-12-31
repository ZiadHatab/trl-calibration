# TRL Calibration

## About this repo

The purpose of this repository is to support the tutorial I wrote on TRL (Thru-Reflect-Line) calibration: [https://ziadhatab.github.io/posts/trl-calibration](https://ziadhatab.github.io/posts/trl-calibration) (at the moment still working on it!).

The code in this repo is not contained in a class, they are collection of functions that you apply at every frequency point. Therefore, it can be in convenience to use if you just want to apply TRL calibration to some measurements you have. For such cases, you are better off with the package [`scikit-rf`](https://github.com/scikit-rf/scikit-rf), as they have TRL implemented and you can do even more.

### Content of the scripts

* `TRL.py` contain a collection of functions that are necessary for the implementation of TRL algorithm. You should think of `TRL.py` as a library of functions that you can load in your main script. You only need here the package [`numpy`](https://github.com/numpy/numpy) as dependency.

* `main.py` is an example to demonstrate how to use the functions in `TRL.py` to apply TRL calibration on actual measurements. Becaus I loaded `.s2p` files and plotted some figures, you additionally need both [`scikit-rf`](https://github.com/scikit-rf/scikit-rf) and [`matplotlib`](https://github.com/matplotlib/matplotlib) installed in your python environment.

* `find_line_length.py` is a script that allow you to design the length of the line standard in a TRL-kit. It has two main functions; one to compute the length of the line standard given frequency limits, while the second is to do the reverse. You only need here [`numpy`](https://github.com/numpy/numpy).


## Homework for you ;)

I personally think the usage of TRL calibration is straightforward easy. I mean, you just measure 3 standards and run your script and you get the answer. Even better, most VNAs will do the computation for you. However, the design of the calibration standards is often something you need to do. One aspect of TRL design that is often overlooked is the intuition on how to choose the length of the line standard. I think it is reasonable to ask questions like: given a length, what are the frequency limits? Can my kit work on multiple frequency bands? Maybe you want to know which max/min length you can use, while still satisfying your target frequency?

So, in the script `main.py` I gave measurements of a TRL-kit based on 50 ohm microstrip line. The line standard has a length of 15 mm. The effective relative permittivity of the microstrip is approximately 2.6. If you are interested in broadening you intuition on TRL calibration, try tackle below tasks:

1. The length of the line is 15 mm and the measurements were done from 0.1 up to 14 GHz. Do you think this TRL-kit is appropriate for this frequency band? If not, then which length would have been more appropriate for this frequency band? With your new length, can you still get a phase margin above 20 degrees? If not, which frequency limit and length would you recommend? Hint: use the script `find_line_length.py`. 

2. If you run the script `main.py` and look at the S21 of calibrated line standard you will notice that it looks pretty good, but if you look at S11 of the calibrated reflect standard (open), you will notice some weirds spikes at around 6.2 GHz and 12.4 GHz. Why do you think S21 of the calibrated line has no spikes, while S11 of the calibrated reflect standard does? Can these spikes be explained? Also, why exactly at those specific frequencies? If the measurements were to continue in frequency, do you think you would see more spikes? If yes, can you predict the frequencies at which the spikes would occur? Hint: use the script `find_line_length.py`, but now use a phase margin of zero.

If you have any questions, feel free to ask (you can post an issue).