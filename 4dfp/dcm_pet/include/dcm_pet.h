/*
          Copyright (C) 2010 Washington University

*/
/* Copyright marker.  Copyright will be inserted above.  Do not remove */
/*
**		     Neuro Imaging Laboratory
**		   Mallinckrodt Institute of Radiology
**		Washington University School of Medicine
** 		
** Last Update:		$Author:  $
** Source File:		$RCSfile:   $
** Revision:		$Revision: 1.1 $
*/
/*
* PT Series Attributes Not Structured:
* 0018 1100 reconstruct_dia:	diameter in mm of the region used for reconstruction
* 0054 1104 detector_lines:	detector lines of response used for tomographic reconstruction
* 0018 0073 acq_start:		how the data collection was started 
* 0018 0071 acq_start_data:	change or trigger to cause data collection to start
* 0018 0000 acq_term:		how the data collection was stopped
* 0018 0075 acq_term_data:	change or trigger to cause data collection to stop
* 0018 1147 fov_shape:		the shape of the PET camera FOV
* 0018 1149 fov_dimensions:	in mm transverse detector followed by axial width
* 0018 1121 gantry_slew:	directional movement of gantry/detector
* 0018 1120 gantry_tilt:	tilt of gantry/detector in degrees
* 0054 0202 detector_motion:	type of detector motion during acq
* 0018 1181 collimator_type:	NONE or RING
* 0018 1180 collimator_name:	label describing the collimator used
* 0054 1200 axial_acceptance:	maximum axial angle accepted
* 0054 1201 axial_mash:		number of adjacent axial lines of response mashed together
* 0054 1202 trans_mash:		number of adjacent transverse lines of response mashed
* 0054 1203 detector_size:	the size of an individual detector element in mm
* 0054 1210 coincidence_width:	the width of the coincidence timing window
* 0054 0013 energy_window_seq:	the sequence of repeating items that describe the energy windows
* 0054 0014 energy_window_low:	the lower limit of the energy window in KeV
* 0054 0015 energy_window_up:	the upper limit of the energy window in KeV
* 0054 1220 secondary_counts:	array defining the type of additional counts
*/ 

typedef struct {
  void *reserved[2];
  int	part10Flag;
  U16	rows;					/* 0028 0010 */
  U16	columns;				/* 0028 0011 */
  U16	bitsAllocated;				/* 0028 0100 */
  U16	bitsStored;				/* 0028 0101 */
  U16	highBit;				/* 0028 0102 */
  U16	pixelRepresentation;			/* 0028 0103 */
  char	seriesUID[DICOM_UI_LENGTH + 1];		/* 0020 000e */
  U16	relativeImgNumber;			/* 0020 0013 */
  int	seriesNumber;				/* 0020 0011 */
  char	photometricInterpretation[DICOM_CS_LENGTH+1]; /* 0028 0004 */
  
  char	pixelSpacing[DICOM_DS_LENGTH * 2 + 2];	/* 0028 0030 */
  char	sliceThickness[DICOM_DS_LENGTH + 1];	/* 0018 0050 */
  
  float positionX;				/* 0020 0032 */
  float positionY;
  float positionZ;
  
  float orientrowX;				/* 0020 0037 */
  float orientrowY;
  float orientrowZ;
  float orientcolX;
  float orientcolY;
  float orientcolZ;
  
  U16	planarConfiguration;			/* 0028 0006 */
  char	studyDate[DICOM_DA_LENGTH + 1];		/* 0008 0020 */
  char	seriesDescription[DICOM_PN_LENGTH + 1];	/* 0008 103e */
  char	protocolName[DICOM_LO_LENGTH + 1];	/* 0018 1030 */
  char	patientName[DICOM_LO_LENGTH + 1];	/* 0010 0010 */
  char	patientID[DICOM_PN_LENGTH + 1];		/* 0010 0020 */
  
  U16	imgLargestPix;
  U16	imgSmallestPix;
  
  char	acqSequenceName[DICOM_LO_LENGTH + 1];	/* 0018 0024 */
  char	idManufacturer[DICOM_LO_LENGTH + 1];	/* 0008 0070 */
  char	idModel[DICOM_LO_LENGTH + 1];		/* 0008 1090 */
  char	idModality[32 + 1];			/* 0008 0060 */
  char	metaTransferUID[DICOM_UI_LENGTH + 1];	/* 0002 0010 */
  char	idInstitution[DICOM_LO_LENGTH + 1];	/* 0008 0080 */
  char	sliceSpacing[DICOM_DS_LENGTH + 1];	/* 0018 0088 */
  int	instanceNumber;
  
  /* PET */
  char	nmi_series_type[DICOM_CS_LENGTH + 1];	/* 0054 1000 NMI Series Type*/
  char	image_correction[DICOM_LO_LENGTH + 1];	/* 0028 0051 IMG Corrected Image*/
  U16	nmi_image_index;			/* 0054 1330 NMI Image Index*/
  U16	nmi_number_of_time_slices;		/* 0054 0101 NMI Number of Time Slices*/
  char	counts_source[DICOM_CS_LENGTH + 1];	/* 0054 1002 NMI Counts Source*/
  U16	nmi_slice_max;				/* 0054 0081 NMI Number of Slices*/
  
  int	img_number_frames;			/* 0028 0008 IMG Number of Frames*/
  
  char  atten_correction[DICOM_LO_LENGTH + 1];	/* 0054 1101 NMI Attenuation Correction Method*/
  char	scatt_correction[DICOM_LO_LENGTH + 1];	/* 0054 1105 NMI Scatter Correction Method*/
  char	decay_correction[DICOM_CS_LENGTH + 1];	/* 0054 1102 NMI Decay Correction*/
  char	rad_volume[DICOM_DS_LENGTH + 1];	/* 0018 1071 ACQ Radiopharmaceutical Volume*/
  char	inject_time[DICOM_TM_LENGTH + 1];	/* 0018 1072 ACQ Radiopharmaceutical Start Time*/
  char	rad_stop[DICOM_TM_LENGTH + 1];		/* 0018 1073 ACQ Radiopharmaceutical Stop Time*/
  char	rad_dose[DICOM_DS_LENGTH + 1];		/* 0018 1074 ACQ Radionuclide Total Dose*/
  char	half_life[DICOM_DS_LENGTH + 1];		/* 0018 1075 ACQ Radionuclide Half Life*/
  char	rad_fraction[DICOM_DS_LENGTH];		/* 0018 1076 ACQ Radionuclide Positron Fraction*/
  float	frame_duration;				/* 0018 1242 ACQ Actual Frame Duration*/
  char	rad[DICOM_LO_LENGTH + 1];	  	/* 0018 0031 ACQ Radiopharmaceutical*/
  
  char	series_time[DICOM_TM_LENGTH + 1];	/* 0008 0031 ID Series Time*/
  char	acq_start_time[DICOM_TM_LENGTH + 1];	/* 0008 0032 ID Acquisition Time*/
  
  float	frame_ref;				/* 0054 1300 NMI Frame Reference Time*/
  char	rescale_intercept[DICOM_AS_LENGTH];	/* 0028 1052 IMG Rescale Intercept=0 for PT*/
  char	rescale_slope[DICOM_DS_LENGTH + 1];	/* 0028 1053 IMG Rescale Slope*/
  char	lossy_comp[DICOM_AS_LENGTH];		/* 0028 2110 IMG Lossy Image Compression=1 or 0*/
  char	nmi_units[DICOM_CS_LENGTH];		/* 0054 1001 NMI Units */
  char	atten_method[DICOM_LO_LENGTH + 1];	/* 0054 1101 NMI Attenuation Correction Method */
  char	decay_factor[DICOM_DS_LENGTH + 1];	/* 0054 1321 NMI Decay Factor */
  char	reconstruct_method[DICOM_LO_LENGTH + 1];/* 0054 1103 Method used for reconstruction */
  char	randoms_method[DICOM_CS_LENGTH];	/* 0054 1100 Randomization correction method */
  
  char	id_code_value[DICOM_SH_LENGTH + 1];	/* 0008 0100 ID Code Value */
  char	id_code_scheme[DICOM_SH_LENGTH + 1];	/* 0008 0102 ID Code Scheme */
  char	id_code_meaning[DICOM_LO_LENGTH + 1];	/* 0008 0104 ID Code Meaning */

  char	fileName[1024];
} PARAMS;

void* extractPixels(DCM_OBJECT** obj, PARAMS *p);
