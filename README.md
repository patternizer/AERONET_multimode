![image](https://user-images.githubusercontent.com/5902974/62958465-a505d600-bdee-11e9-8704-81183f95daad.gif)

# AERONET_multimode

Algorithm to fit multi-lognormal distributions (1-8 modes) and perform statistical hypothesis testing to determine the best fit to a AERONET aerosol volume size distribution (AVSD).

The theory is published in:

Taylor, M., Kazadzis, S., and Gerasopoulos, E.: Multi-modal analysis of aerosol robotic network size distributions for remote sensing applications: dominant aerosol type cases, Atmos. Meas. Tech., 7, 839-858, https://doi.org/10.5194/amt-7-839-2014, 2014. 

## Contents

* `multimode.m` - main script to be run with Matlab
* `20160202_20160202_Exeter_MO.siz` - AERONET size distribution test data
* `20160202_20160202_Exeter_MO.gif` - AERONET plot of size distribution test data
* `GMM1.jpg` - fit to the test data with a single lognormal mode
* `GMM2.jpg` - fit to the test data with a two lognormal modes
* `GMM3.jpg` - fit to the test data with a three lognormal modes
* `BEST.jpg` - best statistical fit to the test data with a four lognormal modes
* `RUN.mat` - statistics arrays

The first step is to clone the latest AERONET_multimode code and step into the check out directory: 

    $ git clone https://github.com/patternizer/AERONET_multimode.git
    $ cd AERONET_multimode
    
### Matlab

The code should run with [Matlab](https://www.mathworks.com/products/matlab.html) although a version will be soon released for Python 3.7+.

The code is designed to take a 22-element vector Y containing the AVSD as input and is run in Matlab with:

    >> [GMM]=multimode(Y)

To run the default case (AERONET Almucantar Level 1.5 Version 3: 2 Feb 2016) issue the command:

   >>  [GMM]=multimode([])
       
## License

The code is distributed under terms and conditions of the [MIT license](https://opensource.org/licenses/MIT).

## Contact information

* [Michael Taylor](https://patternizer.github.io) 
