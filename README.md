# Multi-Spectral-Imaging--Diffusion-Experiments

This is a project that examines the rate of diffusion of 4 different fluorescent dyes namely, IRDye800, IRDye680, Alexa Four 750(AF750), and 
Alexa Four 700(AF700), conjugated to unique antibodies in a 5ml thick of Agar. The rate of diffusion is gauged by the diffusion co-efficient that is 
calculated using the Eistein's equation born of the Eistein's theory of Brownian motion.

      x^2 = 4Dt -> where is the x the radius, D is the diffusion co-efficient of the dye and t is time in seconds.
 
To find D, diffusion co-efficient, I did the following;

I had 18/17 acquistions for each run(experiment), therefore, I had 18/17 folders of unmixed dyes. The first folder contains an unmixed files of 
an empty well (imaged without dye) to account for autofluorescence(background noise) that is accounted for in all other the acquisitions that are imaged
with a dye/that have fluorescence due to the dye.

For each run, I used a suitable acquisition to create 4 regions-of-interest (ROIs) over 4 different cross-sections of the well (horizontal, left diagonal,
vertical line and right diagonal). I used the same ROI generated on all other acquistion to get a consistent result. For all acquistions, I plotted a graph of
fluorescence Gray Scale/Pixel vs distance (cm): (distance converted from pixel to cm: 1cm - 142.003pixels), and found the distance between the peaks that's the x-intercept
on either sides of peak as the diameter. Using this diameter I found the radius, squared it and stored it in an array.

For the time, I imaged every 5 minutes from acqiuistion 2 to acquistion 13, and 30 minutes apart for other acquisitions after the 13th one. Since 
there are four cross-sections, hence 4 x^2, I had 4 slots for each time point. The two arrays, x^2 and t(time), had 68 elements because there are 
17 acquisitions and 4 cross-sections recorded for each time point -> 17*4.

To fit the data to a linear plot, I used ployfit that returns the 2*1 matrix which contains the slope and the y-intercept respectively. I used the 
linspace function to provide the x-values to plot against what ployfit returns. Using the linear plot, I found the diffusion co-efficient.
