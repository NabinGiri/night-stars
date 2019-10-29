import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.convolution import convolve, Gaussian2DKernel, Tophat2DKernel, Box2DKernel
from astropy.utils.data import download_file


def image(url,position,size):
    image_file = download_file(url, cache=True, show_progress='True')
    hdu_lists = fits.open(image_file, memmap='True')
    hdu_lists.info()
    image_data = hdu_lists[0].data
    print(
        "Above is the list of Header Data Unit of FITS file. \nWe will be working on PRIMARY Block. \nThe shape of the array is as below:")
    print((type(image_data)))
    print(image_data.shape)
    hdu_lists.close()
    print("Viewing Image..\nStatistics of the Image: ")
    plt.imshow(image_data, cmap='gray')
    plt.colorbar()
    plt.savefig('Original_Image.png')
    plt.show()
    print('Min:', np.min(image_data))
    print('Max:', np.max(image_data))
    print('Mean:', np.mean(image_data))
    print('Stdev:', np.std(image_data))
    crop = Cutout2D(image_data, position, size)
    plt.imshow(crop.data, origin='upper', cmap='gray')
    plt.savefig('Cropped_Image.png')
    plt.show()

    gauss_kernel = Gaussian2DKernel(0.5)
    smoothed_data_gauss = convolve(crop.data, gauss_kernel)
    plt.imshow(smoothed_data_gauss, cmap='gray')
    plt.savefig('Gaussian Smoothing Filter')
    plt.show()

    tophat_kernel = Tophat2DKernel(2)
    smoothed_data_tophat = convolve(crop.data, tophat_kernel)
    plt.imshow(smoothed_data_tophat, cmap='gray')
    plt.savefig('Tophat Smoothing Filter')
    plt.show()

    box_kernel = Box2DKernel(5)
    smoothed_data_box = convolve(crop.data, box_kernel)
    plt.imshow(smoothed_data_box, cmap='gray')
    plt.savefig('Box_kernel Smoothing Filter')
    plt.show()


if __name__ == '__main__':
    url = 'http://data.astropy.org/tutorials/FITS-images/HorseHead.fits'
    position = (500, 450)
    size = (350, 400)
    image(url,position,size)
    print("Files Saved")

