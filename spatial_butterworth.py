# Code for denoise butterworth filter
# Sec. 5 in Menezes et al. 2014 presentes a description of the main ideas behind this filter and its application in images
  
import numpy as np
from scipy import interpolate, fft
from astropy.io import fits
import math

class butterworth:
    def __init__(self, data,image_type='cube', order=1.0, cutoff_freq=[0.5,0.5], npad=3, nexp=15,filter_shape='elipse2'):
        '''
            data: input data, cube or single image
            order: filtering order. For SPLUS, the best value is 1. Rings start to appear at order = 1.0
            cutoff_freq: cutout frequency of filtering. For S-PLUS, a good value is 0.5
            image_type: 'cube' or 'single', choose if input data is a cube or a single image
            npad: the original image will be enlarger by 2*npad. Default: npad=3
            nexp: second increment (2*nexp) to original image. Default: nexp=15
            filter_shape: final filter shape is built using a combination of functions. Available options are:
                'elipse': uses only one elipse
                'elipse2': combines two elipses
                'elipsevssquare': combines one elipse to a function that is a product of two other elipses
                Default:'elipse2'
        '''
        self.data = data
        self.order = order
        self.cutoff_frequency_x = cutoff_freq[0]
        self.cutoff_frequency_y = cutoff_freq[1]
        self.image_type = image_type
        self.npad = npad
        self.nexp = nexp
        self.filter_shape = filter_shape

    def apply(self):
        if self.image_type=='cube':
            BWflux = np.zeros_like(self.data)
            for i_l in range(self.data.shape[0]):
                print('Applying Butterworth to filter #', i_l, '...')
                zero_image = self.creates_zeroimage(self.data[i_l],self.npad)
                image_expand = self.creates_secondimage(zero_image,self.nexp)
                image_filtered = self.creates_filterimage(image_expand, self.order, self.filter_shape, self.npad + self.nexp)
                BWflux[i_l] = image_filtered
            return BWflux
        if self.image_type=='single':
            print('Applying Butterworth...')
            zero_image = self.creates_zeroimage(self.data, self.npad)
            image_expand = self.creates_secondimage(zero_image,self.nexp)
            image_filtered = self.creates_filterimage(image_expand, self.order, self.filter_shape, self.npad + self.nexp)
            return image_filtered

    def creates_zeroimage(self,data,npad):
        '''
        Creates zero_image
        '''

        # zero image with larger size than the original, it receives the information from the original image
        ny, nx = data.shape
        zero_image = np.zeros((ny + 2 * npad, nx + 2 * npad))

        for i in range(0, nx):
            for j in range(0, ny):
                zero_image[npad + j, npad + i] = data[j, i]

        # mimic image information in the borders (left side)
        for i in range(0, npad):
            for j in range(0, ny):
                zero_image[npad + j, i] = data[j, 0]

        #  mimic image information in the borders (right side)
        for i in range(0, npad):
            for j in range(0, ny):
                zero_image[npad + j, npad + nx + i] = data[j, nx - 1]

        # mimic image information in the borders (bottom)
        for i in range(0, nx):
            for j in range(0, npad):
                zero_image[j, npad + i] = data[0, i]

        # mimic image information in the borders (top)
        for i in range(0, nx):
            for j in range(0, npad):
                zero_image[npad + ny + j, npad + i] = data[ny - 1, i]

        # it remains 4 corners of the image with zeros that need to filled
        # bottom left corner
        for i in range(0, npad):
            for j in range(0, npad):
                zero_image[j, i] = zero_image[j, 3]

        # bottom right corner
        for i in range(0, npad):
            for j in range(0, npad):
                zero_image[j, npad + i + nx] = zero_image[j, npad + nx - 1]

        # top left corner
        for i in range(0, npad):
            for j in range(0, npad):
                zero_image[npad + ny + j, i] = zero_image[npad + ny + j, npad]

        # top right corner
        for i in range(0, npad):
            for j in range(0, npad):
                zero_image[npad + ny + j, npad + nx + i] = zero_image[npad + ny + j, npad + nx - 1]

        return zero_image

    def creates_secondimage(self,imagezero,nexp):
        '''
        Creates an expanded image
        '''
        y1, x1 = imagezero.shape

        # create larger image with borders containing zero and with zero_image at the center
        image_expanded = np.zeros((y1 + 2 * nexp, x1 + 2 * nexp))
        image_expanded[nexp:y1 + nexp, nexp:x1 + nexp] = imagezero

        interpol_values = np.zeros((2, 2))
        interpol_values1 = np.zeros((2, 2))
        interpol_values2 = np.zeros((2, 2))
        interpol_values3 = np.zeros((2, 2))

        temporary_column = np.zeros((y1 + 2, 2))
        temporary_column[y1 + 1, 0] = y1 + 2 * nexp
        temporary_column[1:y1 + 1, 0] = range(nexp, y1 + nexp)

        # image borders must have values dropped to zero.
        interpol_values[1, 0] = nexp
        for i in range(nexp, x1 + 2 * nexp):
            interpol_values[1, 1] = image_expanded[nexp, i]
            interpolate_function = interpolate.interp1d(interpol_values[0:2, 0], interpol_values[0:2, 1], kind='linear')
            image_expanded[0:nexp + 1, i] = interpolate_function(np.arange(nexp + 1))

        interpol_values1[0, 0] = y1 + nexp - 1
        interpol_values1[1, 0] = y1 + 2 * nexp
        interpol_values1[1, 1] = 0.0
        for i in range(nexp, x1 + 2 * nexp):
            interpol_values1[0, 1] = image_expanded[y1 + nexp - 1, i]
            interpolate_function1 = interpolate.interp1d(interpol_values1[0:2, 0], interpol_values1[0:2, 1],
                                                         kind='linear')
            image_expanded[y1 + nexp - 1:y1 + 2 * nexp, i] = interpolate_function1(
                y1 + nexp - 1 + np.arange(nexp + 1))

        interpol_values2[0, 0] = 0.0
        interpol_values2[1, 0] = nexp
        interpol_values2[0, 1] = 0.0
        for j in range(0, y1 + nexp):
            interpol_values2[1, 1] = image_expanded[j, nexp]
            interpolate_function2 = interpolate.interp1d(interpol_values2[0:2, 0], interpol_values2[0:2, 1],
                                                         kind='linear')
            image_expanded[j, 0:nexp + 1] = interpolate_function2(np.arange(nexp + 1))

        interpol_values3[0, 0] = x1 + nexp - 1
        interpol_values3[1, 0] = x1 + 2 * nexp
        interpol_values3[1, 1] = 0.0
        for j in range(0, y1 + 2 * nexp):
            interpol_values3[0, 1] = image_expanded[j, x1 + nexp - 1]
            interpolate_function3 = interpolate.interp1d(interpol_values3[0:2, 0], interpol_values3[0:2, 1],
                                                         kind='linear')
            image_expanded[j, x1 + nexp - 1:x1 + 2 * nexp] = interpolate_function3(
                x1 + nexp - 1 + np.arange(nexp + 1))

        return image_expanded

    def creates_filterimage(self, expandfinalimage, order, elipse_option, ntot):
        '''
        Creates final image after applying the filter
        '''
        yf, xf = expandfinalimage.shape

        # create filter with double size of expandfinalimage, which will have the size of the Fourier transform
        filtro = np.zeros((2 * yf, 2 * xf))

        aa = self.cutoff_frequency_x * xf
        bb = self.cutoff_frequency_y * yf
        cc = self.cutoff_frequency_x * xf
        dd = self.cutoff_frequency_y * yf

        imfilter1 = np.zeros_like(filtro)
        imfilter2 = np.zeros_like(filtro)
        xm = xf
        ym = yf

        for i in range(0, 2 * xf):
            for j in range(0, 2 * yf):
                conta = math.sqrt(((i - xm) / aa) ** 2 + ((j - ym) / bb) ** 2)
                imfilter1[j, i] = 1.0 / (1.0 + conta ** (2 * order))

        for i in range(0, 2 * xf):
            for j in range(0, 2 * yf):
                conta1 = abs(i - xm) / cc
                conta2 = abs(j - ym) / dd
                imfilter2[j, i] = (1.0 / (1.0 + conta1 ** (2 * order))) * (1.0 / (1.0 + conta2 ** (2 * order)))

        if elipse_option == 'elipse':
            filtro = imfilter1
        elif elipse_option == 'elipse2':
            filtro = imfilter1 * imfilter1
        elif elipse_option == 'elipsevssquare':
            filtro = imfilter1 * imfilter2

        # image that will receive the filtering
        outimage = expandfinalimage
        npadd = np.zeros((2 * yf, 2 * xf))

        # filtering
        imagein = expandfinalimage[0:yf, 0:xf]
        imageflag = npadd
        # padding - procedure to multiply certain image pixels by -1, so that Fourier transform ends up with
        # frequency = 0 at the center.
        # padding for each image
        for i in range(0, 2 * xf):
            for j in range(0, 2 * yf):
                if (i <= xf - 1) and (j > yf - 1):
                    if ((i + j) % 2) != 0:
                        flag = 1
                    else:
                        flag = -1
                    imageflag[j, i] = flag * imagein[j - yf, i]
                else:
                    imageflag[j, i] = 0.0

        # apply fourier transform
        twodfft = fft.fft2(imageflag)
        # multiply by the filter
        filtragem = twodfft * filtro
        # apply inverse fourier transform
        twodifft = fft.ifft2(filtragem)
        invfft = np.real(twodifft)
        paddout = npadd
        imageout = imagein

        # remove o padding
        for i in range(0, 2 * xf):
            for j in range(0, 2 * yf):
                if ((i + j) % 2) != 0:
                    flag = 1
                else:
                    flag = -1
                paddout[j, i] = flag * invfft[j, i]
                if (i <= xf - 1) and (j > yf - 1):
                    imageout[j - yf, i] = paddout[j, i]
        outimage[0:yf, 0:xf] = imageout

        # creates the final image by removing the borders
        final_image = outimage[ntot:yf - ntot, ntot:xf - ntot]

        return final_image
