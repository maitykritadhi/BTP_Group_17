## BTP_Group_17 By Ayush Gupta (1901ME13) & Kritadhi Maity (1901ME36)
## Topic : Advance Prediction of Boiling Crisis Through Acoustic Signal Analysis

Here,<br />
1. **User Interface** - We developed an user-interface in order to notice the effect of changing Window-size and Window-overlap on predictive ability of a particular feature. This takes window-size and window overlp input from user and plots a particular feature graph along with original data graph. Window-szie and Window-overlap are also mentioned on top of output graph.

2. **data_feature_std_mean** - This code gives us 4 plots in output. First one is data vs time, second one feature vs time, third one is standard deviation of feature v time and fourth one is mean of feature vs time. Different lines demarcate different boiling regimes. Window size and overlap mentioned in the code.

3. **feature_vs_feature** - This code gives us video of two features plotted at a time.

4. **MeanF1_vs_meanF2** - This code gives us video of mean of two features plotted at a time. There are different color combinations for different boiling regimes (blue for nucleate boiling, red for transition boiling and green again for nucleate boiling - CHF happens at interface of blue and red).All the demarcation in time domain for different boiling regimes are done using unitary method. Alteration can be done to code for plotting standard deviation of two features instead of mean.

5. **color_3D** - This code gives us 3-D video plotting of three features taken at a time. In our code we have F6, F53 and F59 plotted together. There are different color combinations for different boiling regimes (blue for nucleate boiling, red for transition boiling and green again for nucleate boiling - CHF happens at interface of blue and red).All the demarcation in time domain for different boiling regimes are done using unitary method.

6. **distance_3D** - This code gives us the plot of distance vs time of a 3D plot. Distance is calculated using distance formula that is distance between two points is square root of the sum of squares of differences between corresponding cordinates and then is plotted against time.There are different color combinations for different boiling regimes.

7. **Moment** - This is the most important code in our analysis. This code takes input from paper "Deep learning the sound of boiling" - input is a dataset. This code gives you nine different output which is written directly to an excel sheet specified in the code. The nine different output are value of percentage heat flux, which is calculated as heat flux at a particular time divided by heat flux at CHF(maximum heat flux) and multiplied by 100, at nine different thresholds for a particular moment. The code is first used to calculate percentage heat flux for all the moments (2 to 10) - it is simply done by varying the formula of moment in the code- and then same is performed for all other datasets. 
