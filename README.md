PUP
===

PET Unified Pipeline

PUP is a package for processing PET images. It is developed by Yi Su, PhD (ysu001@gmail.com)

The package includes the following subfolders of different components:

4dfp
----
A package developed by Dr. Abraham Z. Snyder (azsnyder@wustl.edu) at Washington University School of Medicine in St. Louis for image analysis and quantification.

src
----
source code for PUP that needs compilation.

scripts
---- 
PUP scripts

AVtemplate
---- 
Florbetapir templates for PET templates-based processing (not required for standard PUP processing).

IDAIF
----
Experimental code for image-derived arterial input function extraction (not required for standard PUP processing). 

matlabcode
----
Matlab code that are used by some of the scripts (not required for standard PUP processing).


INSTALLATION
=============

1. Following the included instructions for the 4dfp package for its installation. This should be done first, and the package is required by PUP. If the user already has 4dfp installed compare the included tools and make sure the existing installation is not missing any components.

2. Compile source code in the src folder using the included makefile using the release target. It will put the binaries at the same location as the 4dfp $RELEASE folder.

3. Copy the content of the scripts folder to desired location or add its location to your path.

______________________
This README is prepared by Yi Su, PhD (ysu001@gmail.com) on 2018/10/16 


