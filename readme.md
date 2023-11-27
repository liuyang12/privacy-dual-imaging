# Privacy Dual Imaging (PDI)
This repository contains the MATLAB code for the paper **Imaging Privacy Threats from an Ambient Light Sensor** (***preprint [submitted]*** 2023) by [Yang Liu](https://yangliu.mit.edu/), [Gregory W. Wornell](http://allegro.mit.edu/~gww/), [William T. Freeman](https://billf.mit.edu/), [Fredo Durand](https://people.csail.mit.edu/fredo/).
[[project]](https://yangliu.mit.edu/publication/privacy_dual_imaging/ "project page")   [[github]](https://github.com/liuyang12/privacy-dual-imaging "github repository")


## Usage
### Prerequisites
- [MatConvNet](https://www.vlfeat.org/matconvnet/) for [FFDNet](https://github.com/cszn/FFDNet) image denoiser. We highly recommend compiling with GPU support. *Note that you might find [this issue](https://github.com/vlfeat/matconvnet/issues/1143 "Problem Compiling with GPU Support on MATLAB 2018a") useful if you compile it on MATLAB 2018a and above.*  

Note: If you do not have MatConvNet setup already, set `USE_MATCONVNET=false` in `test_privacy_dual_imaging.m` and it will use TV denoiser instead.

### Run PnP-QCS on experimental data
1. Test proposed inversion algorithm PnP-QCS on measured pointing hand data (Figure 3) in `test_privacy_dual_imaging.m`.

2. All touch gesture figures in the paper can be reproduced by running corresponding scripts in `figures/` directory. All data are measured from the Samsung View2 tablet (17.3-inch) using .

### Acquire your own data
- [Android Devices] We use the [AndroSensor](https://play.google.com/store/apps/details?id=com.fivasim.androsensor&hl=en_US&gl=US) app made by [Fiv Asim](https://play.google.com/store/apps/developer?id=Fiv+Asim) to acquire the raw data from the ambient light sensor. We set the "Update Interval" to "Very Fast" to maximize the sensor sampling rate. The raw data is stored in `.csv` files in the `rawdata/` directory.
- [Mac Devices] We write a new app `AmbientLight` to acquire the raw data from the ambient light sensor on MacOS (Macbooks and iMacs). Private API of accessing ambient light sensors on macOS hiden by Apple T1/T2 security chips is borrowed and appreciated from [DarkModeBuddy](https://github.com/insidegui/DarkModeBuddy) by [Guilherme Rambo](https://github.com/insidegui). More details about the app can be found in the [AmbientLight](https://github.com/liuyang12/AmbientLight/) repo.

## Structure of directories

| directory  | description  |
| :--------: | :----------- | 
| `AmbientLight` | Light sensor acquisition app on MacOS | 
| `rawdata`    | Raw `.csv` files acquired from ambient light sensors on Android tablet or Mac laptops / desktops. |
| `dataset`    | Processed data from raw datapoints for each measurement |
| `figures`    | MATLAB scripts to reproduce the results and figures in the paper |
| `packages`   | Proposed PnP-QCS algorithm and other adapted state-of-art algorithms |
| `results`    | results of reconstruction (after reconstruction) |
| `utils`      | utility functions |

## Citation
```
@article{Liu23PDI,
   author    = {Liu, Yang and Wornell, Gregory W. and Freeman, William T. and Durand, Fr{\'e}do},
   title     = {Imaging Privacy Threats from an Ambient Light Sensor},
   journal   = {preprint [submitted]},
   year      = {2023},
   type      = {Journal Article}
}
```

## Contact
[Yang Liu, MIT CSAIL](mailto:yliu@csail.mit.edu "Yang Liu, MIT CSAIL") 