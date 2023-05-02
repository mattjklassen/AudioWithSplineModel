# AudioWithSplineModel
audio with spline modeling

To do a simple test on this version, as of April 1, 2022:

1. open A450.wav
2. play back sound file with space bar
3. set frequency guess to 450 (should be already default)
4. push "graph signal" and then "shade cycles" buttons
5. there should be 604 cycles, and scrolling all the way to the right should not crash the app now (ha ha)
6. select "No Cycle Interp" (should be default)
7. click "use Delta Model" (this allows cycle lengths to be recomputed without endpoint necessarily at zero crossings)
8. select k=50 on subintervals slider (should be default)
9. click "Compute Model" (this takes a while, under 4 minutes, but now has a complete model for each cycle.) 
10. (use k=20 for much faster but less accurate model)
11. click "Graph Model" (graphs entire model in blue, covering the grey signal graph)
12. click "Play Model" (should sound the same as original, using about half the samples so no surprise)
13. select "Regular Cycle Interp" (and leave "use Delta Model" selected, which is in memory)
14. select m=20 on key cycle mult slider (will use every 20th cycle as key, from Delta Model)
15. click "Compute Model" (this is basically instant, forming interpolated cycles with B-splines in between key cycles from Delta model)
17. click "Graph Model" twice (once to erase old model, once more to graph new)
18. Note: key cycles are still accurate and are shaded bright green 
19. Note: intermediate cycles are only approximating the shape of original cycles, move around and zoom in and out to see more
20. click "Play Model" (this is pretty close to the original despite only using about 2.5% of the original amount of data)
21. the data reduction is approximate, based on 1/2 of samples per cycle, and 1/20 of cycles, so 1/40 = 0.025 of original data

As of May 17, 2022, I added new instrument Bcoeff files, so here is a test to run:

1. Run the program and when it opens push the button one back from the right side called "load Bcoeffs". 
2. This reads in a bunch of binary files of floats [instrument].keyBcoeffs and produces and loads a graph.
3. Play the output file by hitting the spacebar or clicking Play.
4. You should hear a chromatic scale starting on middle C, which starts out sounding like a french horn and then morphs into a cello.
