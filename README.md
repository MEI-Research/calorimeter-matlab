# Calorimeter Matlab processing scripts

Matlab scripts designed to process calorimeter data from CalRQ

- Push_Human

Process human study data from push-style calorimeters

- Push_Infusion_Pull_Equivalent

Process infusion study data from push-style calorimeters, using a simulated transform to a pull calorimeter

For simplicity, this is the recommended processing method.

- Push_Infusion_Corrected_Derivative

Process infusion study data from push-style calorimeters, correcting the derivative terms to account for differences between nitrogen displacement and oxygen consumption

This is not the recommended processing method due to complexity.

- SaveReport

Save result from processing with above scripts as an Excel file

- ShrinkFile

Retime a CalRQ file to minute-by-minute data