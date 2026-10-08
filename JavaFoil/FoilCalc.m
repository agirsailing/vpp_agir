function [cl cd] = FoilCalc(Foil,Re)

%==========================================================================
% FoilCalc
% 
% FoilCalc is a tool for calculating lift and drag coefficients for 2-d
% general 2-d foils. The script uses JavaFoil and detailed documentation 
% of JavaFoil may be found here: http://www.mh-aerotools.de/airfoils/javafoil.htm
%
%==========================================================================

% Define JavaFoil path
javapath = [pwd,'/JavaFoil/'];

% Generate JavaFoil inputfile
filepath = [javapath 'script.jfscript'];
genJava(Re, Foil ,filepath);

% Run analysis
system(['java -cp "',javapath,'mhclasses.jar" -jar "',javapath,'javafoil.jar" Script="',filepath,'"']);

% Read results
[cl, cd] = readJava(filepath);


