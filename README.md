# SHARP-track pipline for HirschLab

## Setup

- download [allenCCF tools](https://github.com/cortex-lab/allenCCF/tree/master) into MATLAB [`userpath`](https://www.mathworks.com/help/matlab/ref/userpath.html)
- download [npy-matlab](https://github.com/kwikteam/npy-matlab) into MATLAB `userpath`
- download [allenCCF dataset](http://data.cortexlab.net/allenCCF/), and update the `data_path` accordingly.
- download this repository into MATLAB `userpath`

## Guide

- create a dedicated working directory (just to keep the original data intact)
- copy histology images into the working directory
- run `pipeline_*` scripts for each task

### `pipeline_probe_CCF` for coordinates of each shank

- register histology
- label probe tracks
    - the very first point should be the tip of each track
- export probe track coordinates

## Notes

- Coordinates of probe tracks
    - each mat file for probe tracks has
        - `probe_tip`: the CCF coordinate of the probe tips
- For visualization
    - check out Yinan's code on the lab servers
    - check out the [brainrender](https://github.com/brainglobe/brainrender) python package