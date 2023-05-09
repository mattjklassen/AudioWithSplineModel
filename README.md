# AudioWithSplineModel
audio with spline modeling

To do a test on this version, as of May 8, 2023:

1. open A450.wav then push "graph signal" button to see graph
2. play back sound file with space bar (or push "play sound" button)
3. set frequency guess to 450 (should be already default)
4. push "shade cycles" button (computes cycles with zero crossings based on frequency guess)
5. there should be 604 cycles verified by scrolling all the way to the right (note: bottom of cycle is number of samples)
6. go back near the beginning and zoom in on one cycle with two finger spread (on Mac trackpad) or use left/right arrows, or mouse wheel
7. when one cycle takes up about half of window you can see samples are graphed as dots
8. use amplitude slider to scale the height of the graph for more detail (useful near the end of the audio graph)
9. left click anywhere on shaded cycle and select "graph cycle spline" (graphs cubic spline of that cycle with default k=30 subintervals)
10. change to k=10 subintervals with slider and select "graph cycle spline" (twice: once off, once to regraph)
11. click "plot targets" to see n=k+3 dots which are the spline interpolation points
12. left click and select "graph cycle spline" again to turn off, then select drop down "Regular Cycle Interpolation" (default)
13. put k=12 (or around 10-15) and m=5 (default) to do cycle interpolation model, then click "compute model"
14. click "Graph Model" to see blue graph which is not very accurate with k around 12, but is more accurate in middle of audio graph
15. every fifth cycle is pale green with spline model computed, the rest in between are computed with B-spline coefficient interpolation
16. click "Play Model" to hear distorted audio reproduction but still recognizable
17. note cycle shape discontinuity near cycle 360 (explained in paper: https://azrael.digipen.edu/research/SplineModelingAudioSignals.pdf)
18. click "Graph Model" to turn off graph, then select "No Cycle Interp" 
19. check box for "use Delta Model" (this allows cycle lengths to be recomputed without zero crossings)
20. select k=30 on subintervals slider (which uses n=k+3=33 B-splines per cycle)
21. click "Compute Model" (this takes a few minutes as delta model is optimizing for cycle shape continuity)
22. click "Graph Model" and slide over to cycle 360 to see that model preserves cycle shape and is not forced to use zero crossings
23. click "Play Model" (should sound close to original, using about one third of the samples so no surprise)
24. select "Regular Cycle Interp" (and leave "use Delta Model" selected, which is in memory)
25. select m=20 on key cycle mult slider (will use every 20th cycle as key, from Delta Model)
26. click "Compute Model" (this is basically instant, forming interpolated cycles with B-splines in between key cycles from Delta model)
27. click "Graph Model" twice (once to erase old model, once more to graph new)
28. Note: key cycles are still accurate and are shaded bright green 
29. Note: intermediate cycles are only approximating the shape of original cycles, move around and zoom in and out to see more
30. click "Play Model" (this is pretty close to the original despite only using about 2.5% of the original amount of data)
    The data reduction is approximate, based on 1/2 of samples per cycle, and 1/20 of cycles, so 1/40 = 0.025 of original data.

Here is a test to load instrument B-spline coefficient files:

1. Run the program and when it opens push the button one back from the right side called "load Bcoeffs". 
2. This reads in a bunch of binary files of floats [instrument].keyBcoeffs and produces and loads a graph.
3. Play the output file by hitting the spacebar or clicking Play.
4. You should hear a chromatic scale starting on middle C, which starts out sounding like a guitar and then morphs into a flute.
5. Note: these sounds are interpolated by simply using the key cycles.
