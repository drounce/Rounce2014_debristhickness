%% OVERVIEW: DEBRIS THICKNESS ESTIMATES FROM LANDSAT 7 THERMAL IMAGERY
% Objective: derive debris thickness based on surface temperature from
% Landsat 7 thermal band satellite imagery and an energy balance model.
% This code and details on its methods and results have been published in
% The Cryosphere.  Any use of the code, or portions of the code, should 
% cite:
%
% Rounce, D.R. and McKinney, D.C.: Debris thickness of glaciers in the
% Everest area (Nepal Himalaya) derived from satellite imagery using a
% nonlinear energy balance model, The Cryosphere, 8, 1317-1329,
% doi:10.5194/tcd-8-1317-2014, 2014.

% Inputs required to run model:
%  - Debris properties
%  - Digital elevation model (raster)
%  - Slope (raster)
%  - Aspect (raster)
%  - Hillshade (raster)
%  - Surface temperature (raster)
%  - Meteorological data

% Output:
%  - geotiff of debris thickness [m]

% Assumptions:
%  - The latent heat flux and heat flux supplied by rain are zero based on
%    the assumption that the debris is dry and it is not raining.
%  - The thermal conductivity is constant through the debris layer.
%  - The nonlinear temperature gradient within the debris can be
%    approximated using the nonlinear approximation factor, Gratio.
%  - Sensible heat flux assumes neutral atmosphere (no corrections made for
%    an unstable atmosphere)

% Notes on geotiffs:
% -The .tif files were exported from ArcGIS using the spatial analyst
%  tools. Make sure that all of the .tif files are properly co-registered 
%  and have the same number of rows and columns. Hillshade is computed
%  using the Hillshade tool in ArcGIS with 'model shadows' selected.  The 
%  position of the sun can be calculated using the NOAA solar calcs script 
%  or NOAA Solar Calculator online 
%  (https://www.esrl.noaa.gov/gmd/grad/solcalc/azel.html).
% -The output does not differentiate between on-glacier and off-glacier
%  pixels.  This should be done by the user in post-processing. 

% Note on code development:
%  In the spirit of providing open-source code, users are encouraged to
%  share any modifications to the code via github. If you have any 
%  questions regarding the use of the code or improvements, please contact
%  David Rounce (david.rounce@utexas.edu)

%% INPUT DATA:
% Debris properties:
Albedo = 0.3;           % debris-cover albedo [-] (Nicholson and Benn, 2012)
z0_Benn = 0.016;        % surface roughness [m] based on very rough ice and satrugi (Brock, 1997) and previous
                        %  calculations over debris-covered surfaces (Takeuchi and others, 2000)
k_eff = 0.96;           % effective thermal conductivity [W m-1 K-1] based on temeprature profiles
                        %  from Imja-Lhotse Shar Glacier
G_ratio_UTC_4_5 = 2.7;  % nonlinear approximation factor [-] based on temperature profiles 
                        %  from Imja-Lhotse Shar Glacier
                        
% Filenames of .tifs:
filename_elevation = 'MaskTR_DEM.tif';
filename_slope = 'MaskTR_Slope.tif';
filename_aspect = 'MaskTR_Aspect.tif';
filename_hillshade = 'L7_20021004_hillshade.tif';
filename_surfacetemp = 'Landsat7_20021004_degK_MaskTR.tif';
% Output filename
date_start = datestr(now, 'YYYYmmdd');
output_filename_debristhickness = ['DebrisThickness_Example_L7_20021004_generated',num2str(date_start),'.tif'];
output_refcode = 32645;
%  output reference coordinate system
%  (http://geotiff.maptools.org/spec/geotiff6.html#6.3.3.1)
                       
% AWS site data:
% Elevation [m.a.s.l.]   
Elev_PyrStat = 5035;     
% Slope [deg] (assume AWS flat on top of ridge, i.e., Sky View Factor = 1)
Slope_PyrStat_deg = 0;  
Slope_PyrStat_rad = Slope_PyrStat_deg.*(pi/180);
% Aspect [deg clockwise from N]
Aspect_PyrStat_deg = 0;
Aspect_PyrStat_rad = Aspect_PyrStat_deg.*(pi/180);
% Latidue [deg N]
lat_deg = 27.9;
lat_rad = lat_deg.*(pi/180);  
% Longitude [deg E]
long_deg = 86.9;
long_rad = long_deg.*(pi/180);
% Time zone [hr] (Nepal is +5:45)
timezone = 5.75;
        
% Meteorological data (example data from Pyramid Station)
za = 2;                         % Height of air temperature instrument [m]
Sin_PyrStat = 935;              % Incoming Shortwave Radiation [W/m2] 
Lin_PyrStat = 196.465;          % Incoming Longwave Radiation [W/m2]
Tair_PyrStat = 5.57+273.15;     % Air Temperature [deg K]
RH_PyrStat = 42.2;              % Relative Humidity [%]
RH_PyrStat = RH_PyrStat/100;    % Divide by 100 to convert to decimal
u_PyrStat = 3.195;              % Wind speed [m/s]        
q_PyrStat = 4.28;               % Specific Humidity [g/kg]
q_PyrStat = q_PyrStat * 10^-3;  % convert to kg/kg

% Date and time of satellite image (example from Landsat 7 October 4, 2002)
Julian_Day_of_year = 277;       % Julian day of year
Month = 10;                     % Month
Day   = 4;                      % Day
Year  = 2002;                   % Year
UTC_Time = 4.5;                 % Time according to UTC
Time = UTC_Time + timezone;

% Net energy threshold [W m-2] used to avoid generating unrealistically 
% high debris thickness from low values of positive energy 
netenergy_threshold = 10;

% Constants:
I0 = 1368;              % Solar constant (W/m2)
transmissivity=0.75;    % Vertical atmospheric clear-sky transmissivity (assumed)
emissivity = 0.95;      % Emissivity of debris surface (Nicholson and Benn, 2006)
P0 = 101325;            % Standard Atmospheric Pressure (Pa)       
density_air_0 = 1.29;   % Density of air (kg/m3)     
LapseRate = 0.0065;     % Temperature Lapse Rate (K/m)
Kvk = 0.41;             % Von Karman's Constant
Lv = 2.49*10^6;         % Latent Heat of Vaporation of water
cA = 1010;              % Specific Heat Capacity of Air (J kg-1 K-1)
R = 461;                % Gas constant for water vapor (J K-1 kg-1)

%% Import .tif files
% Elevation [m.a.s.l]
Elevation = imread(filename_elevation);
Elevation = double(Elevation);
    nrows = size(Elevation,1);
    ncols = size(Elevation,2);
% Slope [deg]
Slope_deg = imread(filename_slope);
Slope_deg = double(Slope_deg);
Slope_rad = Slope_deg.*(pi/180);
% Aspect [deg clockwise from North]
Aspect_deg = imread(filename_aspect); % deg CW from N
Aspect_deg = double(Aspect_deg);
Aspect_rad = Aspect_deg.*(pi/180);
% Hillshade [value of 0 is a shadow]
Hillshade = imread(filename_hillshade);
Hillshade = double(Hillshade);
% Surface temperature [deg K]
Ts = imread(filename_surfacetemp);
Ts = double(Ts);

% Remove bad values from .tifs
Slope_deg(Slope_deg < 0) = 0;
Slope_rad(Slope_deg < 0) = 0;
Aspect_deg(Aspect_deg < 0) = 0;
Aspect_rad(Aspect_deg < 0) = 0;
Elevation(Elevation < 0) = 0;
                
%% Additional parameters needed for script
% Pressure at Pyramid Station [Pa]
P_PyrStat = P0*exp(-0.0289644*9.81*Elev_PyrStat/(8.31447*288.15));
% Pressure at each pixel [Pa]
P = P0*exp(-0.0289644*9.81*Elevation/(8.31447*288.15));
P(Elevation < 0) = 0;
% Adjust wind speed [m/s] from 5 m to 2 m accounting for surface roughness
u = u_PyrStat*(log(2/z0_Benn)/(log(5/z0_Benn)));
% Dimensionless transfer coefficient [-]
A_Benn = Kvk^2/(log(za/z0_Benn))^2;
% Adjust air temperature [deg K] based on elevation using the lapse rate
T_air = Tair_PyrStat - LapseRate*(Elevation - Elev_PyrStat);
T_air(Elevation <= 0) = 0;
% Vapor pressure
eZ_Saturated = 611*exp(-Lv/R*(1./T_air-1/273.15));
eZ = eZ_Saturated.*RH_PyrStat;
        
%% Compute the position of the sun (NOAA Solar Calculator)
% Source: Astronomical Algorithms used by NOAA
% Julian Day 
Year_NOAA = Year;
%  change name to Year_NOAA, so value of 'Year' is not altered for other scripts
if Month <= 2
    Month = Month + 12;
    Year_NOAA = Year_NOAA - 1;
end
A = floor(Year_NOAA/100);
B = -13;
JulianDay_no_time_zone = floor(365.25*(Year_NOAA+4716)) + floor(30.6001*(Month+1)) + (Day + Time/24) + B - 1524.5;
%  Time zone not accounted for, i.e., in Greenwich Time Zone
%  floor() rounds numbers down
% Julian Day based on local time
JulianDay = JulianDay_no_time_zone - timezone/24;
% Julian Century
JulianCentury = (JulianDay-2451545)/36525;
% Geom Mean Long Sun (deg)
%  mod returns the remainder after division
GeomMeanLongSun_deg = mod(280.46646+JulianCentury.*(36000.76983+JulianCentury*0.0003032),360);
GeomMeanLongSun_rad = GeomMeanLongSun_deg.*(pi/180);
% Geom Mean Anom Sun (deg)
GeomMeanAnomSun_deg = 357.52911+JulianCentury.*(35999.05029-0.0001537*JulianCentury);
GeomMeanAnomSun_rad = GeomMeanAnomSun_deg.*(pi/180);
% Eccent Earth Orbit
EccentEarthOrbit = 0.016708634-JulianCentury.*(0.000042037+0.0000001267*JulianCentury);
% Sun Eq of Ctr
SunEqofCtr = sin(GeomMeanAnomSun_rad).*(1.914602-JulianCentury.*(0.004817+0.000014*JulianCentury))+sin(2*GeomMeanAnomSun_rad).*(0.019993-0.000101*JulianCentury)+sin(3*GeomMeanAnomSun_rad)*0.000289;
% Sun True Long (deg)
SunTrueLong_deg = GeomMeanLongSun_deg + SunEqofCtr;
% Sun True Anom (deg)
SunTrueAnom_deg = GeomMeanAnomSun_deg + SunEqofCtr;
SunTrueAnom_rad = SunTrueAnom_deg.*(pi/180);
% Sun Rad Vector (AUs)
SunRadVector = (1.000001018*(1-EccentEarthOrbit.*EccentEarthOrbit))./(1+EccentEarthOrbit.*cos(SunTrueAnom_rad));
% Sun App Long (deg)
SunAppLong_deg = SunTrueLong_deg-0.00569-0.00478*sin((125.04-1934.136*JulianCentury).*(pi/180));
SunAppLong_rad = SunAppLong_deg.*(pi/180);
% Mean Obliq Ecliptic (deg)
MeanObliqEcliptic_deg = 23+(26+((21.448-JulianCentury.*(46.815+JulianCentury.*(0.00059-JulianCentury*0.001813))))/60)/60;
MeanObliqEcliptic_rad = MeanObliqEcliptic_deg.*(pi/180);
% Obliq Corr (deg)
ObliqCorr_deg = MeanObliqEcliptic_deg+0.00256*cos((125.04-1934.136*JulianCentury).*(pi/180));
ObliqCorr_rad = ObliqCorr_deg.*(pi/180);
% Sun Rt Ascen (deg)
SunRtAscen_deg = 180/pi*atan((cos(ObliqCorr_rad).*sin(SunAppLong_rad))./cos(SunAppLong_rad));
% Sun Declin (deg)
SunDeclin_deg = 180/pi*asin(sin(ObliqCorr_rad).*sin(SunAppLong_rad));
SunDeclin_rad = SunDeclin_deg.*(pi/180);
% VarY
VarY = tan(ObliqCorr_deg/2.*(pi/180)).*tan(ObliqCorr_deg/2.*(pi/180));
% Eq of Time (min)
EqofTime = 4*180/pi*(VarY.*sin(2*GeomMeanLongSun_rad)-2*EccentEarthOrbit.*sin(GeomMeanAnomSun_rad)+4*EccentEarthOrbit.*VarY.*sin(GeomMeanAnomSun_rad).*cos(2*GeomMeanLongSun_rad)-0.5*VarY.*VarY.*sin(4*GeomMeanLongSun_rad)-1.25*EccentEarthOrbit.*EccentEarthOrbit.*sin(2*GeomMeanAnomSun_rad));
% True Solar Time (min)
TrueSolarTime = mod((Time*60*1440+Time*60+EqofTime+4*long_deg-60*timezone),1440);
% Hour Angle (deg)
if TrueSolarTime/4 < 0
    HourAngle_deg = TrueSolarTime/4+180;
else
    HourAngle_deg = TrueSolarTime/4-180;
end
HourAngle_rad = HourAngle_deg.*(pi/180);
% Solar Zenith Angle (deg)
SolarZenithAngle_deg = 180/pi*acos(sin(lat_rad).*sin(SunDeclin_rad)+cos(lat_rad).*cos(SunDeclin_rad).*cos(HourAngle_rad));
SolarZenithAngle_rad = SolarZenithAngle_deg.*(pi/180);
% Solar Elevation Angle (deg)
SolarElevationAngle_deg = 90-SolarZenithAngle_deg;
SolarElevationAngle_rad = SolarElevationAngle_deg.*(pi/180);
% Approx Atmospheric Refraction (deg)  
if SolarElevationAngle_deg > 85
    ApproxAtmosRefrac_deg = 0;
elseif SolarElevationAngle_deg > 5
    ApproxAtmosRefrac_deg = 58.1./tan(SolarElevationAngle_rad)-0.07./((tan(SolarElevationAngle_rad)).^3)+0.000086./((tan(SolarElevationAngle_rad)).^5);
elseif SolarElevationAngle_deg > -0.575
    ApproxAtmosRefrac_deg = 1735+SolarElevationAngle_deg.*(-518.2+SolarElevationAngle_deg.*(103.4+SolarElevationAngle_deg.*(-12.79+SolarElevationAngle_deg*0.711)));
else
    ApproxAtmosRefrac_deg = -20.772./tan(SolarElevationAngle_rad);
end
ApproxAtmosRefrac_deg = ApproxAtmosRefrac_deg/3600;
% Solar Elevation Correct for Atm Refraction (deg)
SolarElevationAngleCorr_deg = SolarElevationAngle_deg + ApproxAtmosRefrac_deg;
% Solar Zenith Angle Corrected for Atm Refraction (deg)
SolarZenithAngleCorr_deg = 90 - SolarElevationAngleCorr_deg;
SolarZenithAngleCorr_rad = SolarZenithAngleCorr_deg.*(pi/180);
% Solar Azimuth Angle (deg CW from N)    
if HourAngle_deg > 0
    SolarAzimuthAngle_deg = ((180/pi*(acos(((sin(lat_rad).*cos(SolarZenithAngle_rad))-sin(SunDeclin_rad))./(cos(lat_rad).*sin(SolarZenithAngle_rad))))+180)/360-floor((180/pi*(acos(((sin(lat_rad).*cos(SolarZenithAngle_rad))-sin(SunDeclin_rad))./(cos(lat_rad).*sin(SolarZenithAngle_rad))))+180)/360))*360;
else
    SolarAzimuthAngle_deg = ((540-180/pi*(acos(((sin(lat_rad).*cos(SolarZenithAngle_rad))-sin(SunDeclin_rad))./(cos(lat_rad).*sin(SolarZenithAngle_rad)))))/360-floor((540-180/pi*(acos(((sin(lat_rad).*cos(SolarZenithAngle_rad))-sin(SunDeclin_rad))./(cos(lat_rad).*sin(SolarZenithAngle_rad)))))/360))*360;
end  
SolarAzimuthAngle_rad = SolarAzimuthAngle_deg.*(pi/180);
% Distance from sun based on eccentricity of orbit (dbar/d)^2
J = 0:364; 
phi_d = (J.*2*pi)./365; %Julian Date radians (eq. A.5 Hartmann 1994)
d0 = 1.000110.*cos(0*phi_d);
d1 = 0.034221.*cos(1*phi_d) + 0.001280.*sin(1.*phi_d);
d2 = 0.000719.*cos(2*phi_d) + 0.000077.*sin(2.*phi_d);
d = d0 + d1 + d2; %(dbar/d)^2-one value for each day of year
% (rm/r)^2 for specific day of satellite image
rm_r2 = d(Julian_Day_of_year);
  
%% Calculate thermal resistance
% Use Benn model with latent heat flux assumed to be zero.  If using a high
% resolution DEM, then you need to resample the grids such that the 
% incoming solar radiation agrees with the resolution of the surface
% temperature raster.

theta = zeros(nrows,ncols);
I = zeros(nrows,ncols);
Sin = zeros(nrows,ncols);
Rn_Benn = zeros(nrows,ncols);
H_Benn = zeros(nrows,ncols);
LE_Benn = zeros(nrows,ncols);
EnergyFlux_Benn = zeros(nrows,ncols);
Qc_Benn = zeros(nrows,ncols);
TR_Benn = zeros(nrows,ncols);
NetEnergy_Benn = zeros(nrows,ncols);

% Angle of Incidence b/w normal to grid slope and solar beam at PyrStat
theta_PyrStat = acos(cos(Slope_PyrStat_rad).*cos(SolarZenithAngleCorr_rad)+sin(Slope_PyrStat_rad).*sin(SolarZenithAngleCorr_rad).*cos(SolarAzimuthAngle_rad-Aspect_PyrStat_rad));

% Potential Clear-Sky Solar Radiation at Pyramid Station
I_PyrStat = I0*rm_r2*transmissivity.^(P_PyrStat./(P0*cos(SolarZenithAngleCorr_rad))).*cos(theta_PyrStat);
I_PyrStat(I_PyrStat < 0) = 0;
%  incoming shortwave radiation cannot be negative
            
% Compute incoming solar radiation (Sin) for each pixel
for r=1:nrows
    for c=1:ncols
        if Elevation(r,c) < 0
            TR_Benn(r,c) = -1;
            Sin(r,c) = 0;
        else
            % Angle of Incidence b/w normal to grid slope and solar beam
            theta(r,c) = acos(cos(Slope_rad(r,c)).*cos(SolarZenithAngleCorr_rad)+sin(Slope_rad(r,c)).*sin(SolarZenithAngleCorr_rad).*cos(SolarAzimuthAngle_rad-Aspect_rad(r,c)));
            % Potential Clear-Sky Solar Radiation
            I(r,c) = I0*rm_r2*transmissivity^(P(r,c)/(P0*cos(SolarZenithAngleCorr_rad))).*cos(theta(r,c));
            if SolarZenithAngleCorr_rad > (pi/2)
                I(r,c) = 0;
            end
            % Effect of Shading (from Hock and Noetzli, 1997)
            if I_PyrStat == 0
                Sin = 0;
            else
                if Hillshade(r,c) == 0
                    Sin(r,c) = 0.15*I_PyrStat;
                else
                    Sin(r,c) = I(r,c).*Sin_PyrStat./I_PyrStat;
                end
            end
        end
    end
end
Sin(Sin < 0) = 0;
%  incoming shortwave radiation cannot be negative

% Solve Energy Balance for Thermal Resistance
for r = 1:nrows
    for c = 1:ncols
        % Avoid unrealistic values due to no surface temperature data
        if Ts(r,c) > 200      
            Rn_Benn(r,c) = Sin(r,c)*(1-Albedo) + emissivity*(Lin_PyrStat-(5.67*10^-8)*Ts(r,c)^4);
            H_Benn(r,c) = density_air_0*(P(r,c)/P0)*cA*A_Benn*u*(T_air(r,c)-Ts(r,c));
            LE_Benn(r,c) = 0;
        else
            Rn_Benn(r,c) = 0;
            H_Benn(r,c) = 0;
            LE_Benn(r,c) = 0;
        end
        EnergyFlux_Benn(r,c) = Rn_Benn(r,c)+H_Benn(r,c)+LE_Benn(r,c);
        % If net energy is positive, then compute the thermal resistiance
        if EnergyFlux_Benn(r,c) >= netenergy_threshold & Ts(r,c) > 273.15
            TR_Benn(r,c) = G_ratio_UTC_4_5*(Ts(r,c)-273.15)/EnergyFlux_Benn(r,c);
        else
            TR_Benn(r,c) = 0;
        end
        % Compute Qc and check that energy balance equals zero (if positive energy flux)
        if TR_Benn(r,c) > 0
            Qc_Benn(r,c) = -G_ratio_UTC_4_5*(Ts(r,c)-273.15)/TR_Benn(r,c);
        else
            Qc_Benn(r,c) = 0;
        end
        NetEnergy_Benn(r,c) = Rn_Benn(r,c)+H_Benn(r,c)+LE_Benn(r,c)+Qc_Benn(r,c);
    end
end

% Debris thickness [m]
DebrisThickness = TR_Benn * k_eff;

% Export geotiff
reference_tif = worldfileread(getworldfilename(filename_elevation),'planar',size(Elevation));
geotiffwrite(output_filename_debristhickness, DebrisThickness, reference_tif, 'CoordRefSysCode',output_refcode)