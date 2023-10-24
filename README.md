
# cosmosoc project

Using python and its repositories extracting data from fits file

### SOFTWARE REQUIREMENT
Python  
Matplotlab  
Astropy package       
Numpy     
Ipykernel   
vs code   

note:  
 must have python,pylance and jupyter installed


```python
#for installing ipykernel

for mac 
pip3 install ipykernel

for windows 
pip install ipykernel



#for installing astropy

for mac 
pip3 install astropy

for windows 
pip install astropy

# for installing "numpy", "panda" and "matplotlab" similar as above
```
## to display the contents present in fits file
```python
from astropy.io import fits
test1=fits.open('lightcurve_1.fits')
test1.info()
print(test1)
```
output:  
Filename: lightcurve_1.fits   
No.   Name   Ver    Type      Cards   Dimensions   Format     
0  PRIMARY    1 PrimaryHDU      55   ()      
1  LIGHTCURVE 1 BinTableHDU    156   3862R x 20C   [D, E, J, E, E, E, E, E, E, J, D, E, D, E, D, E, D, E, E, E]   
  2  APERTURE      1 ImageHDU        49   (10, 9)   int32   
[<astropy.io.fits.hdu.image.PrimaryHDU object at 0x128d4a7e0>, <astropy.io.fits.hdu.table.BinTableHDU object at 0x10e7bc200>, <astropy.io.fits.hdu.image.ImageHDU object at 0x128d57380>]

```python
k=test1[1].data.columns
print(k)
```
output:  
    name = 'TIME'; format = 'D'; unit = 'BJD - 2454833';  disp = 'D14.7'.  
    name = 'TIMECORR'; format = 'E'; unit = 'd'; disp = 'E13.6'  
    name = 'CADENCENO'; format = 'J'; disp = 'I10'    
    name = 'SAP_FLUX'; format = 'E'; unit = 'e-/s'; disp = 'E14.7'   
    name = 'SAP_FLUX_ERR'; format = 'E'; unit = 'e-/s'; disp = 'E14.7'.  
    name = 'SAP_BKG'; format = 'E'; unit = 'e-/s'; disp = 'E14.7'. 
    name = 'SAP_BKG_ERR'; format = 'E'; unit = 'e-/s'; disp = 'E14.7'       
    name = 'PDCSAP_FLUX'; format = 'E'; unit = 'e-/s'; disp = 'E14.7'     
    name = 'PDCSAP_FLUX_ERR'; format = 'E'; unit = 'e-/s'; disp = 'E14.7'     
    name = 'SAP_QUALITY'; format = 'J'; disp = 'B16.16'    
    name = 'PSF_CENTR1'; format = 'D'; unit = 'pixel'; disp = 'F10.5'    
    name = 'PSF_CENTR1_ERR'; format = 'E'; unit = 'pixel'; disp = 'E14.7'     
    name = 'PSF_CENTR2'; format = 'D'; unit = 'pixel'; disp = 'F10.5'      
    name = 'PSF_CENTR2_ERR'; format = 'E'; unit = 'pixel'; disp = 'E14.7'      
    name = 'MOM_CENTR1'; format = 'D'; unit = 'pixel'; disp = 'F10.5'      
    name = 'MOM_CENTR1_ERR'; format = 'E'; unit = 'pixel'; disp = 'E14.7'       
    name = 'MOM_CENTR2'; format = 'D'; unit = 'pixel'; disp = 'F10.5'       
    name = 'MOM_CENTR2_ERR'; format = 'E'; unit = 'pixel'; disp = 'E14.7'        
    name = 'POS_CORR1'; format = 'E'; unit = 'pixels'; disp = 'E14.7'         
    name = 'POS_CORR2'; format = 'E'; unit = 'pixels'; disp = 'E14.7'         
)





# scatterplot of sap_flux
```python
from astropy.io import fits
import matplotlib.pyplot as plt
hdulist=fits.open('lightcurve_1.fits')
sap_flux=hdulist[1].data['SAP_FLUX']
time=hdulist[1].data['Time']
plt.figure(figsize=(10,4))
plt.plot(time,sap_flux, 'b.')
plt.xlabel("Time")
plt.ylabel("Flux")
plt.show()
hdulist.close()
```
output:        
you will get flux vs time graph but will have negative value too 
can change it to pdcsap instead and sap quality>0 will get only positive data    
 
```python
from astropy.io import fits
import matplotlib.pyplot as plt

hdulist = fits.open('lightcurve_1.fits')

pdc_sap_flux = hdulist[1].data['PDCSAP_FLUX']

sap_quality=hdulist[1].data['PDCSAP_FLUX']

time=hdulist[1].data['Time']

#filter data 
pdc_sap_flux_filtered = pdc_sap_flux[sap_quality>0]
time_filtered = time[sap_quality>0]

#plot
plt.figure(figsize=(10,4))
plt.plot(time,pdc_sap_flux, 'r.',label='All Data')
plt.plot(time_filtered,pdc_sap_flux_filtered, 'c',label='sap_quality>0')
plt.title('PDCSAP Flux')
plt.xlabel("Time")
plt.ylabel("Flux")
plt.show()

hdulist.close()
```

# Display the array of pixels as an image
```python
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np

hdulist = fits.open('lightcurve_1.fits')

image_data=hdulist[2].data
hdulist.close()
#image_data = image_data.astype(int)

plt.figure(figsize=(8,8))
plt.imshow(image_data, cmap='gray', origin='lower', interpolation='none')
plt.colorbar()
plt.title('fits image data')
plt.show()

plt.imshow(image_data, cmap='plasma', origin='lower', interpolation='none')
plt.colorbar()
plt.title('fits image data')
plt.show()

px=np.ndarray(shape=(7,7), dtype=float)
plt.imshow(px,cmap='plasma')
plt.show()
```
output:     
three images


# BONUS QUESTION

## to display the contents present in fits file

```python
from astropy.io import fits
bonus=fits.open('bonus_image_1.fits')
bonus.info()
print(bonus)
```
output:
Filename: bonus_image_1.fits
No.    Name      Ver    Type      Cards   Dimensions   Format
  0  PRIMARY       1 PrimaryHDU     161   (891, 893)   int16   
  1  er.mask       1 TableHDU        25   1600R x 4C   [F6.2, F6.2, F6.2, F6.2] 
  [<astropy.io.fits.hdu.image.PrimaryHDU object at 0x16c826390>, <astropy.io.fits.hdu.table.TableHDU object at 0x16c5c2de0>]

### image
```python

k=bonus[0].data
import matplotlib.pyplot as plt
plt.figure()
plt.imshow(k, origin='lower')
plt.show()

plt.imshow(k, origin='lower',cmap='gray')
plt.show()
```

output:    
images one color and other gray


```python
from astropy.io import fits
bonus=fits.open('bonus_image_1.fits')

k3=bonus[1].data

import matplotlib.pyplot as plt
plt.show()
print(k3)

h=bonus[1].data.columns
print(h)
```
output:      
[(-3.12, -3.12, 0.09, 0.04) (-2.96, -3.12, 0.02, 0.07)
 (-2.8, -3.12, -0.07, 0.15) ... (2.8, 3.12, 0.0, 0.0)
 (2.96, 3.12, 0.0, 0.0) (3.12, 3.12, 0.0, 0.0)]
ColDefs(
    name = 'XI'; format = 'F6.2'; unit = 'DEGREES'; start = 1      
    name = 'ETA'; format = 'F6.2'; unit = 'DEGREES'; start = 7       
    name = 'XI_CORR'; format = 'F6.2'; unit = 'ARCSEC'; start = 13      
    name = 'ETA_CORR'; format = 'F6.2'; unit = 'ARCSEC'; start = 19
)







