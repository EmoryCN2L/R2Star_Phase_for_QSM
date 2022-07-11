# A Probabilistic Bayesian Approach to Recover R2* map and Phase Images for QSM
A nonlinear approximate message passing (AMP) framework is proposed to recover R2* map and phase images from undersampled measurements. The proposed AMP-PE approach combines information from both the sparse prior and mono-exponential decay model, it treats the parameters as unknown varialbes and automatically estimates them.

* If you use this code and find it helpful, please cite the above paper. Thanks :smile:
```
@ARTICLE{Bayesian_R2Star_2022,
    author={Shuai Huang and James J. Lah and Jason W. Allen and Deqiang Qiu},
    journal={Magnetic Resonance in Medicine},
    title={A Probabilistic Bayesian Approach to Recover R2* map and Phase Images for Quantitative Susceptibility Mapping},
    pages={1-19},
    year={2022},
}
```



## Summary
```
    ./AMP_MRI_T2Star_Phase_3D	-- This folder contains MATLAB files to perform 3D reconstruction using the AMP-PE approach.
    ./AMP_MRI_T2Star_Phase_2D	-- This folder contains MATLAB files to perform 2D reconstruction using the AMP-PE approach.
    ./L1_MRI_T2Star_Phase_3D	-- This folder contains MATLAB files to perform 3D reconstruction using the l1-norm regularization (L1) approach.
    ./L1_MRI_T2Star_Phase_2D	-- This folder contains MATLAB files to perform 2D reconstruction using the l1-norm regularization (L1) approach.
    ./data			-- The data folder contains simulated data "Sim1" derived from the 2019 QSM challenge 2.0
    ./result			-- The result folder
```

## Dataset Preparation
The datasets need to be placed in the data folder. We first define the following notations:
```
	sx	-- the number of samples along the x direction (readout direction)
	sy	-- the number of samples along the y direction
	sz	-- the number of samples along the z direction
	Ne	-- the number of echoes
	Nc	-- the nunber of channels (sensitivity coils)
```

* 3D datasets
```
	1) Multi-echo datasets are ordered according to the echo indices. Suppose there are 4 echoes, there will be 4 datasets named: "echo_1_3d.mat", "echo_2_3d.mat", "echo_3_3d.mat", "echo_4_3d.mat". Take the dataset from the 1st echo "echo_1_3d.mat" for example, it is a MAT-file containing a 4-dimensional array named "data". The size of "data" is "sx by sy by sz by Nc".
	2) Sensitivity map dataset is named as "sensitivity_map_3d.mat". It is a MAT-file containing a 4-dimensional array named "maps_3d". The size of "maps_3d" is "sx by sy by sz by Nc".
	3) The full-sampling location mask is named as "full_sampling_loc.mat". It is a MAT-file containing a 2-dimensional array named "full_sampling_loc". The size of "full_sampling_loc" is "sy by sz". The 0-1 array consists of 0s and 1s, it is a mask showing the sampling locations in the full-sampling case.
```

* 2D datasets
```
	1) Multi-echo dataset is saved in the file "slice_2d.mat" , it is a MAT-file containing a 4-dimensional array named "data". The size of "data" is "sy by sz by Ne by Nc".
	2) Sensitivity map dataset is named as "sensitivity_map_2d.mat". It is a MAT-file containing a 3-dimensional array named "maps". The size of "maps" is "sy by sz by Nc".
	3) The sampling location mask is named as "full_sampling_loc.mat". It is a MAT-file containing a 2-dimensional array named "full_sampling_loc". The size of "full_sampling_loc" is "sy by sz". The 0-1 array consists of 0s and 1s, it is a mask showing the sampling locations in the full-sampling case.
```

## Usage
You can follow the following steps to run the program. 

`Detailed comments` are within each demo file.

Open `MATLAB` and type the following commands into the console:

* Step 0) Generate simulated measurements based on the Sim1 dataset from the 2019 QSM challenge 2.0
```
    >> cd ./data/Sim1
    >> Sim1_generate_measurements
```

* Step 1) Reconstruction using the AMP-PE approach. The parameters are adaptively and automatically estimated from the data.
```
    >> cd ./AMP_MRI_T2Star_Phase_3D	% perform 3D reconstruction
    >> addpath('.')
    >> mri_wave_amp_3d_x0_r2star_nw_sep_dict	% this is the demo file, look inside for detail comments
    >>
    >> cd ../AMP_MRI_T2Star_Phase_2D	% perform 2D reconstruction
    >> addpath('.')
    >> mri_wave_amp_2d_x0_r2star_nw_sep_dict	% this is the demo file, look inside for detail comments
```
* Step 2) Reconstruction using the L1 approach. The regularization parameters need to be `tuned manually` for different types of datasets.
```
    >> cd ./L1_MRI_T2Star_Phase_3D	% perform 3D reconstruction
    >> addpath('.')
    >> mri_rec_wave_l1_reg_3d	% this is the demo file, look inside for detail comments
    >>
    >> cd ../L1_MRI_T2Star_Phase_2D	% perform 2D reconstruction
    >> addpath('.')
    >> mri_rec_wave_l1_reg_2d	% this is the demo file, look inside for detail comments
```
* Step 3) Reconstruction using the least square approach
```
    >> cd ./L1_MRI_T2Star_Phase_3D	% perform 3D reconstruction
    >> addpath('.')
    >> mri_rec_wave_lsq_3d	% this is the demo file, look inside for detail comments
    >>
    >> cd ../L1_MRI_T2Star_Phase_2D	% perform 2D reconstruction
    >> addpath('.')
    >> mri_rec_wave_lsq_2d	% this is the demo file, look inside for detail comments
```
