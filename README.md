![image](https://user-images.githubusercontent.com/5902974/59154328-fe0c6300-8a67-11e9-9261-4d79fcf8ee94.png)

# AERONET_multimode

Algorithm to fit multi-lognormal distributions (1-8 modes) and perform statistical hypothesis testing to determine the best fit to a AERONET aerosol volume size distribution (AVSD).

The theory is published in:

Taylor, M., Kazadzis, S., and Gerasopoulos, E.: Multi-modal analysis of aerosol robotic network size distributions for remote sensing applications: dominant aerosol type cases, Atmos. Meas. Tech., 7, 839-858, https://doi.org/10.5194/amt-7-839-2014, 2014. 

[doi.org/10.5194/amt-7-839-2014](https://doi.org/10.5194/amt-7-839-2014)

## Contents

* `multimode.m` - main script to be run with Matlab
* `160128_160128_Exeter_MO.siz` - AERONET size distribution test data
* `160128_10128_Exeter_MO.gif` - AERONET plot of size distribution test data
* `GMM1.jpg` - fit to the test data with a single lognormal mode
* `GMM2.jpg` - fit to the test data with a two lognormal modes
* `GMM3.jpg` - fit to the test data with a three lognormal modes
* `GMM4.jpg` - fit to the test data with a four lognormal modes

The first step is to clone the latest AERONET_multimode code and step into the check out directory: 

    $ git clone https://github.com/patternizer/AERONET_multimode.git
    $ cd AERONET_multimode
    
### Matlab

The code should run with [Matlab](https://www.mathworks.com/products/matlab.html) although a version will be soon released for Python 3.7+.

The code is designed to take a 22-element vector Y containing the AVSD as input and is run in Matlab with:

    >> [GMM]=multimode(Y)

To run the default case (AERONET Aerosol Inversions L1.5 Version 2 for EXETER_MO, 28/01/2016, 10:29:08 GMT) issue:

   >>  [GMM]=multimode([])

Alternatively, issue the commands:
	 
   >>  Y = [0.000262,0.001774,0.005377,0.007652,0.005960,0.003439,0.002175,0.002049,0.003043,0.005746,0.009857,0.013233,0.016114,0.019276,0.022546,0.024126,0.022326,0.016942,0.009724,0.003899,0.001043,0.000184];
   >> [GMM]=multimode(Y)
       
## License

The code is distributed under terms and conditions of the [MIT license](https://opensource.org/licenses/MIT).

## Contact information

* [Michael Taylor](https://patternizer.github.io) 
