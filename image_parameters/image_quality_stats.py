import re
from os import path
from math import ceil
import logging

import pylab
import numpy as np
from optparse import OptionParser

logging.basicConfig(level=logging.INFO)

BACKGROUD_SPHERE_CENTER = [(71.349, 50.031), (-0.0, 76.), (-71.349, 50.031), (-95.146, 22.169),
                           (-108.574, -11.922), (-110.173, -48.527), (61., -84.), (93.787, -71.414),
                           (109.732, -40.122), (0.0, -84.), (-50., -84.), (-84.409, -73.597)]
# 'name': (x, y, z, r)
ROI = {'sphere22in': (-28.6, -49.54, 37., 11.0),
       'sphere17in': (-57.2, 0.0, 37., 8.5),
       'sphere13in': (-28.6, 49.54, 37., 6.5),
       'sphere10in': (28.6, 49.54, 37., 5.0)}

parser = OptionParser()

parser.add_option("-i", "--input",
                  action="store", dest="input", metavar="PATH", type="string",
                  help="Path out file png")
parser.add_option("-o", "--out",
                  action="store", dest="output", metavar="PATH", type="string",
                  help="Path input header interfile")


def interfile_parser(file_name):
    """
    Parse interfile header to dict, get image size, scaling factor (mm/pixel)...
    :param file_name: srt interfile path
    :return:
    """
    f = open(file_name, 'r')
    param = {}
    for line in f.readlines():
        matchObj = re.match(r'(.*) := (.*)', line, re.M | re.I)
        if matchObj:
            param[matchObj.group(1)] = matchObj.group(2)
    try:
        param['size'] = (int(param['!matrix size [3]']), int(param['!matrix size [2]']),
                         int(param['!matrix size [1]']))
    except KeyError:
        raise Exception("Bad parsing Matrix size")

    if param['!number format'] == 'float':
        if param['!number of bytes per pixel'] == '4':
            param["type"] = np.float32
        else:
            raise Exception("Bad number format")
    else:
        raise Exception("Bad number format")
    try:
        param['scaling_factor_xy'] = float(param['scaling factor (mm/pixel) [1]'])
        param['scaling_factor_z'] = float(param['scaling factor (mm/pixel) [3]'])
    except:
        raise Exception("Bad parsing scaling_factor")
    try:
        param['offset_z'] = float(param['first pixel offset (mm) [3]'])
    except:
        raise Exception("Bad parsing z offset")
    param['path_to_data_file'] = path.join(path.dirname(file_name), param['name of data file'])
    return param


def interfile2array(param):
    """
    conveft interfile to numpy array
    :param param: dict contens interfile poarametrs
    :return:
    """

    f = open(param['path_to_data_file'], 'r')
    v_list = np.fromfile(f, dtype=param['type'])
    f.close()
    resh_arr = np.asarray(v_list).reshape(param['size'])
    return resh_arr[:, ::-1, ::-1]  # rotate


def txt2array(file_names, pram_dict, filip=False):
    """
    conveft interfile to numpy array
    :param param: dict contens interfile poarametrs
    :return:
    """
    arr = np.empty((pram_dict["size_z"], pram_dict["size_xy"], pram_dict["size_xy"]))
    i = 0
    for file in file_names:
        f = open(file, 'r')
        arr[i, :, :] = np.loadtxt(f)
        f.close()
        i += 1
    if flip:
        arr = arr[:, ::-1, ::-1]
    return arr


def txt2interfile(file_names, pram_dict, ouput_name='test.v'):
    """
    Convert interfile to numpy array
    :param param: dict contens interfile poarametrs
    :return:
    """
    image = txt2array(file_names, pram_dict)
    image.astype(np.float32).tofile(ouput_name)


def pos2arry_index(pos, scaling_factor, size):
    """
    Convert position in cm to array index
    """
    return int((size / 2.0 + (pos * (1 / scaling_factor))))


def get_sphere_measure(image, x0, y0, z0, r, scaling_factor_xy, scaling_factor_z, size_xy,
                       size_z, debug=False):
    """
    Return mean value in  sphere, center (x0,y0,z0) with radius r
    """

    idx_x0 = pos2arry_index(x0, scaling_factor_xy, size_xy)
    idx_y0 = pos2arry_index(y0, scaling_factor_xy, size_xy)
    idx_z0 = pos2arry_index(z0, scaling_factor_z, size_z)


    radius_xy = (r * (1 / scaling_factor_xy))
    idx_radius = int(round(radius_xy))

    # build metric matrix, value is distance from center of cell to cener of matrix
    size = 2*int(round(radius_xy))+1 
    metric_array = np.empty((size,size))
    delta = int(size/2)
    for i in range(size):
        for j in range(size):
            metric_array[i,j] = ((i-delta)**2+(j-delta)**2)
    # mask
    index = metric_array<radius_xy**2 

    roi_values_list = image[idx_z0, idx_y0 - idx_radius:idx_y0 + idx_radius + 1,
               idx_x0 - idx_radius:idx_x0 + idx_radius + 1][
        index].flatten()

    if debug:
        # DEBUG draw ROI mask
        logging.debug(" x:%f y:%f z:%f", x0, y0, z0)
        logging.debug("   %d   %d   %d", idx_x0, idx_y0, idx_z0)
        image[idx_z0, idx_y0 - idx_radius:idx_y0 + idx_radius + 1,
        idx_x0 - idx_radius:idx_x0 + idx_radius + 1][index] = -0.01
        logging.debug("numer of pixels:%d", len(roi_values_list))

    return np.mean(roi_values_list)


def show_image(_array, xy_lim, norm=1, zoom=1, title=None):
    """Show image. If title exist, save figure to *title*.png"""
    pylab.figure()
    # linear normalize
    _array = (_array - min(_array.flatten()))
    _array = _array / norm

    pylab.xlabel('X [mm]')
    pylab.ylabel('Y [mm]')
    x_size, y_size = np.shape(_array)
    center_x = int(x_size / 2)
    center_y = int(y_size / 2)
    delta_x = int(center_x / (zoom))
    delta_y = int(center_y / (zoom))
    temp_array = _array[center_x - delta_x:center_x + delta_x,
                 center_y - delta_y:center_y + delta_y]
    flip_matrix = temp_array[::-1, :]  # rotate
    pylab.imshow(flip_matrix, cmap="gray",
                 extent=[-xy_lim / zoom, xy_lim / zoom, -xy_lim / zoom, xy_lim / zoom])
    pylab.colorbar()
    if title != None:
        logging.debug("save image as %s", title)
        pylab.savefig(title)
    else:
        pylab.show()


def CRC(roi_mean, background_mean, activity_ratio=4):
    """ Contrast Recovery Coefficient"""
    return ((roi_mean / background_mean - 1) / (activity_ratio - 1))


def BV(background_std, background_mean):
    """Background Variability"""
    return background_std / background_mean


def CNR(roi_mean, background_mean, background_std):
    """Contrast-to-Noise Ratio"""
    return abs(roi_mean - background_mean) / background_std


def measure(image, params, volume_name, title=None):
    "Salculates NEMA statistics, save image as *title*"

    parameters = params
    x, y, z, r = ROI[volume_name]
    z = z - params['offset_z']

    scaling_factor_xy = parameters['scaling_factor_xy']
    scaling_factor_z = parameters['scaling_factor_z']
    size_xy = parameters['size'][1]
    size_z = parameters['size'][0]

    background_values = []
    for delta_z in [-20, -10, 0, 10, 20]:  # in mm
        for _x, _y in BACKGROUD_SPHERE_CENTER:
            background_values.append(
                get_sphere_measure(image, _x, _y, z + delta_z, r, scaling_factor_xy,
                                   scaling_factor_z, size_xy, size_z))
    background_mean = np.mean(background_values)
    background_std = np.std(background_values)

    roi_mean = get_sphere_measure(image, x, y, z, r, scaling_factor_xy, scaling_factor_z, size_xy,
                                  size_z)

    show_image(image[pos2arry_index(z, scaling_factor_z, size_z), :, :],
               xy_lim=scaling_factor_xy * (size_xy / 2),
               norm=background_mean, zoom=1, title=title)
    logging.debug("bg_std : %.4f    bg_mean :%.4f    mean_roi %.4f", background_std,
                  background_mean, roi_mean)
    logging.debug("bv_std : %.4f    crc_:%.4f    cnr %.4f", BV(background_std, background_mean),
                  CRC(roi_mean, background_mean), CNR(roi_mean, background_mean, background_std))
    return BV(background_std, background_mean), CRC(roi_mean, background_mean)


if __name__ == "__main__":
    (options, args) = parser.parse_args()
    interfile_header = options.input
    if interfile_header:
        params = interfile_parser(interfile_header)

        arr = interfile2array(params)
        print("ROI   BV    CRC")
        for roi in ROI:
            print(roi + " %.4f %.4f" % measure(arr, params, roi, title=options.output))
    else:
        print("Error select  input header interfile")
