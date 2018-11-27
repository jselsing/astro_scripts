#!/usr/local/anaconda3/envs/py36 python
# -*- coding: utf-8 -*-

# Imports
import numpy as np
from astroquery.irsa_dust import IrsaDust
import astropy.coordinates as coord
import astropy.units as u
from dust_extinction.parameter_averages import F04
import glob


def get_ebv(ra, dec):
    """Query IRSA dust map for E(B-V) value and returns reddening array
    ----------

    ra : float
        Right Ascencion in degrees
    dec : float
        Declination in degrees

    Returns
    -------
    ebv : float

    Notes
    -----
    For info on the dust maps, see http://irsa.ipac.caltech.edu/applications/DUST/
    """


    C = coord.SkyCoord(ra*u.deg, dec*u.deg, frame='fk5')

    dust_table = IrsaDust.get_query_table(C, section='ebv', timeout=60)
    ebv = dust_table["ext SandF ref"][0]

    return ebv


def main():

    ra = 197.449524
    dec = -23.38097

    ebv = get_ebv(ra, dec)

    bands = np.sort(glob.glob("../data/engrave/*"))
    band_out, central_wave, A_lambda = [0]*len(bands), [0]*len(bands), [0]*len(bands)
    for ii, band in enumerate(bands):
        band_name = band.split("/")[-1].split(".")[0]

        filt = np.genfromtxt(band)
        lamb_T, T = filt[:, 0], filt[:, 1]
        cwav = np.mean(lamb_T[T > max(T)/100])

        # initialize the model
        ext = F04(Rv=3.1)

        try:
            reddening = 1/ext.extinguish(cwav*u.angstrom, Ebv=0.32)
        except ValueError:
            print("%s outside wavelength coverage of extinction law. Replacing with nan"%band_name)
            reddening = np.nan

        band_out[ii], central_wave[ii], A_lambda[ii] = band_name, cwav, -2.5*np.log10(reddening)

    np.savetxt("../output/GW170817.dat", list(zip(band_out, np.around(central_wave, 2), np.around(A_lambda, 3))), fmt='%s %s %s', header="Extinction magnitudes for E(B-V) = %s"%ebv)


if __name__ == '__main__':
    main()
