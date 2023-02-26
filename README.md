# TRL Calibration

## About this repo

The purpose of this repository is to support the tutorial I wrote on TRL (Thru-Reflect-Line) calibration: <https://ziadhatab.github.io/posts/trl-calibration>.

The code in this repo is not contained in a class, it is a collection of functions that you apply at each frequency point. The code is written this way for the sake of the tutorial. If you just want to apply TRL, it may be more convenient to use the [`scikit-rf`](https://github.com/scikit-rf/scikit-rf) package instead, as it has TRL implemented and can do even more.

### Content of the scripts

* [`TRL.py`](TRLpy) contains a collection of functions necessary for implementing the TRL algorithm. You should think of [`TRL.py`](TRLpy) as a library of functions that you can load into your main script. You only need the [`numpy`](https://github.com/numpy/numpy) package as a dependency.

* [`main.py`](mainpy) is an example of how to use the functions in [`TRL.py`](TRLpy) to apply TRL calibration to real measurements. Since I have loaded `.s2p` files and plotted some figures, you will also need both [`scikit-rf`](https://github.com/scikit-rf/scikit-rf) and [`matplotlib`](https://github.com/matplotlib/matplotlib) installed in your python environment.

* [`find_line_length.py`](find_line_lengthpy) is a script that allows you to determine the length of the line standard in a TRL kit. It has two main functions; one is to compute the length of the line standard given frequency limits, while the second is to do the reverse. All you need is [`numpy`](https://github.com/numpy/numpy).

## Homework for you ;)

Personally, I think the use of TRL calibration is straightforward and easy. I mean, you just measure 3 standards and run your script and you get the answer. Even better, most VNAs will do the calculation for you. However, the design of the calibration standards is often something you have to do. One aspect of TRL design that is often overlooked is the intuition about how to choose the length of the line standard. I think it is reasonable to ask questions like: Given a length, what are the frequency limits? or Can my TRL kit work at multiple frequency bands? Maybe you want to know what max/min length you can use and still meet your target frequency?

So in the script [`main.py`](mainpy) I have given measurements of a TRL kit based on a 50 Ohm microstrip line. The standard line has a length of 15 mm. The effective relative permittivity of the microstrip is about 2.6. If you are interested in broadening your intuition on TRL calibration, try tackling the following tasks:

1. The length of the line is 15 mm and the measurements were made from 0.1 to 14 GHz. Do you think this TRL kit is appropriate for this frequency band? If not, what length would have been more appropriate for this frequency band? Can you still get a phase margin above 20 degrees with your new length? If not, what frequency limit and length would you recommend? Hint: Use the script [`find_line_length.py`](find_line_lengthpy).

2. If you run the script [`main.py`](mainpy) and look at S21 of the calibrated line standard, you will notice that it looks pretty good, but if you look at S11 of the calibrated reflect standard (open), you will notice some strange spikes at about 6.2 GHz and 12.4 GHz. Why do you think S21 of the calibrated line has no spikes while S11 of the calibrated reflect standard does? Can you explain these spikes? And why at these specific frequencies? If the measurements were continued in frequency, do you think you would see more spikes? If so, can you predict the frequencies at which the spikes would occur? Hint: Use the [`find_line_length.py`](find_line_lengthpy) script, but now use a phase margin of zero.

If you have any questions, feel free to ask (you can post an issue).

[TRLpy]: https://github.com/ZiadHatab/trl-calibration/blob/main/TRL.py

[mainpy]: https://github.com/ZiadHatab/trl-calibration/blob/main/main.py

[find_line_lengthpy]: https://github.com/ZiadHatab/trl-calibration/blob/main/find_line_length.py
