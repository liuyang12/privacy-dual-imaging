# Privacy Dual Imaging
Imaging Privacy Threats from an Ambient Light Sensor

by [Yang Liu, MIT CSAIL](https://liuyang12.github.io/)


## Usage
0. Save the recorded raw data `.csv` file in `./rawdata/csv/` and rename it following the format `touch16x16_22cm_20fps_NNN.csv` by simply change the last three digits accordingly.
1. Go to `./tests/` and run `screenCamReconCS.m`.
2. Change the file number according to the filename, i.e., `scfile=[NNN];`
3. Run a single pass of the first block or the whole block of `screenCamReconCS.m` to see whether the red line (threshold or `params.measthresh`) hit right the middle of the measurements. Adjust that number accordingly, so the recognized number of measurements equal the number of valid frames on the screen. In that case, the following step will get you the reconstruction of the image, as the demo shows.
4. Two cases it might fail:
   - Cannot get good segmentation because the sensor fails to gather enough data points per screen frame. *Solution:* Use a lower framerate of the video clip.
   - The sensor gets enough data points per screen frame, but the threshold is not correct. Adjust that threshold in `params.measthresh` for better segmentation and re-run the script. 
5. If everything (screen frame rate and measurement threshold) seems correct and it still fails to get any sore of image reconstruction, just acquire another set of raw data and try it again. 
