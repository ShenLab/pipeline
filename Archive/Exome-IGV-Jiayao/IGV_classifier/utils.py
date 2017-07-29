from PIL import Image
import numpy as np
import os
import re

from keras.models import Sequential
from keras.layers.core import Dense, Dropout, Activation, Flatten
from keras.layers.advanced_activations import PReLU
from keras.layers.convolutional import Convolution2D, MaxPooling2D,ZeroPadding2D
from keras.optimizers import SGD, Adadelta, Adagrad, Adam
from keras.utils import np_utils, generic_utils
from keras.regularizers import l2
from keras import backend as K
from six.moves import range
import sys
import random
import matplotlib.pyplot as plt



def process_one_image(image_adr, n_block, new_width, new_height):

    img = Image.open(image_adr)
    width, height = img.size   # Get dimensions
    img_grey = img.convert('L')
    image_data = np.array(img_grey)

    # find left, right, up down black boundary
    row1 = image_data[0,:]
    for index, pixel in enumerate(row1):
        black_left = index + 4
        if all(row1[index: index+4] != [0, 0, 0, 0]):
            break
    black_right = next(index  for index, pixel in enumerate(row1[::-1]) if pixel != 0) - 1
    black_right = width - black_right
    col_mid = image_data[:,black_left + black_right / 2 ]
    for index, pixel in enumerate(col_mid):
        if all(col_mid[index: index+4] == [0, 0, 0, 64]):
            top = index
            break
    rev_col = col_mid[::-1]
    for index, pixel in enumerate(rev_col):
        if all(rev_col[index: index+4] == [64 , 0, 0, 0]):
            bottom = height - index
            break
    # crop black border
    img = img.crop((black_left, top, black_right, bottom))

    # crop to center
    width, height = img.size   # Get dimensions
    img_grey = img.convert('L')
    image_data = np.array(img_grey)
    center = width / 2
    block_size = int(round((width)/ 41.0)) # all images has 41bp, this can be modified(using the coverage line)
    left = center  - int(block_size * (n_block + 0.5))
    right = center + int(block_size * (n_block + 0.5)) - 4
    new_left, new_right = 0, 0 
    for index in range(int(center - block_size * 0.6) , int(center - block_size * 0.4)):
        pixles = image_data[:,index]
        for pos, pixel in enumerate(pixles[:-20]):
            if all(pixles[pos: pos+20] == [100, 0, 0, 0] * 5) or all(pixles[pos: pos+20] == [100, 0] * 10):
                new_left = index
                break
    for index in range(int(center + block_size * 0.4) , int(center + block_size * 0.6)):
        pixles = image_data[:,index]
        for pos, pixel in enumerate(pixles[:-20]):
            if all(pixles[pos: pos+20] == [100, 0, 0, 0] * 5) or all(pixles[pos: pos+20] == [100, 0] * 10):
                new_right = index
                break          
    if new_left != 0 and  new_right != 0:
        center = (new_left + new_right) / 2
        block_size = new_right - new_left
        left = new_left  - int(block_size * n_block)
        right = new_right + int(block_size * n_block + 0.5)
        # change indels into snp like
        img_rgb = img.convert('RGB')
        img_rgb_data = np.array(img).transpose((2, 0, 1))
        for pos, pixel in enumerate(image_data[:,new_right - 1]):
            if pixel == 74:
                img_rgb_data[0][pos, new_left+1:new_right] = [118] * (new_right - new_left - 1)
                img_rgb_data[1][pos, new_left+1:new_right] = [24] * (new_right - new_left - 1)
                img_rgb_data[2][pos, new_left+1:new_right] = [220] * (new_right - new_left - 1)
        img = Image.fromarray(img_rgb_data.transpose((1, 2, 0)))
    img = img.crop((left, 0, right, height))

    # remove white in the middle
    while True:
        width, height = img.size
        img_grey = img.convert('L')
        image_data = np.array(img_grey)
        middle_str  = ','.join( map(str,image_data[:,width / 2 ]))
        white = r'((250,){20,})'
        m = re.search(white,middle_str)
        if m == None:
            break
        else:
            white_before, _ =  m.span()
            white_before = len(middle_str[:white_before].split(','))
            white_after = white_before + len(m.group().split('250,')) 
            upper = img.crop((0, 0, width, white_before))
            lower = img.crop((0, white_after, width, height))
            imgs_comb = np.vstack( (upper, lower ) )
            img = Image.fromarray( imgs_comb)
    # resize the images 
    img = img.resize((new_width, new_height), Image.BILINEAR)
    return img

def load_data(folder, n_block, image_width, image_height, color_channel = 3):
    '''
    Load the data from folder.
    
    Parameters
    ----------
    :type folder: str 
    :param folder: image folder address

    :type n_block: int
    :param n_block: the blocks to keep around base pair of interest
    
    :type image_width: int
    :param image_width: new image width
 
    :type image_height: int
    :param image_height: new image height

    returns
    ----------
    :type dataset: np.uint8
        image after process, ndarray of dimension: (# of images, channel , width, height)
    '''
    image_files = os.listdir(folder)
    dataset = np.ndarray(shape=(len(image_files), color_channel, image_width, image_height),
                         dtype=np.uint8)  
    image_index = 0
    for image in image_files:
        image_file = os.path.join(folder, image)
        #print ' x',image_file
        img = process_one_image(image_file, n_block, image_width, image_height)
        if image_index % 500 == 0:
            print ' x', image_index, image
            img.show()
        image_data = np.array(img).transpose((2, 0, 1))
        dataset[image_index, :, :, :] = image_data
        image_index += 1
    return dataset


def maybe_npy(data_folders, n_block, image_width, image_height, force=False):
    '''
    save the data to npy.
    
    Parameters
    ----------
    :type data_folders: list 
    :param data_folders: list of folders of images

    :type n_block: int
    :param n_block: the blocks to keep around base pair of interest
    
    :type image_width: int
    :param image_width: new image width
 
    :type image_height: int
    :param image_height: new image height
    
    :type force: boolean
    :param force: whether overwrite
    

    returns
    ----------
    :type dataset_names: list
        saved file name
    '''
        
    dataset_names = []
    for folder in data_folders:
        set_filename = folder + '.npy'
        dataset_names.append(set_filename)
        if os.path.exists(set_filename) and not force:
            # You may override by setting force=True.
            print('%s already present - Skipping npying.' % set_filename)
        else:
            print('npying %s.' % set_filename)
            dataset = load_data(folder, n_block, image_width, image_height)
            try:
                np.save(set_filename, dataset)
            except Exception as e:
                print('Unable to save data to', set_filename, ':', e)
    return dataset_names

def get_minst_model(weights_path=None):

    nb_classes = 6
    # input image dimensions
    img_rows, img_cols = 128, 128
    # number of convolutional filters to use
    nb_filters = 32
    # size of pooling area for max pooling
    nb_pool = 2
    # convolution kernel size
    nb_conv = 3

    model = Sequential()

    model.add(Convolution2D(nb_filters, nb_conv, nb_conv,
                            border_mode='valid',
                            input_shape=(3, img_rows, img_cols)))
    model.add(Activation('relu'))
    model.add(Convolution2D(nb_filters, nb_conv, nb_conv))
    model.add(Activation('relu'))
    model.add(MaxPooling2D(pool_size=(nb_pool, nb_pool)))
    model.add(Dropout(0.25))
    
    model.add(Convolution2D(nb_filters, nb_conv, nb_conv))
    model.add(Activation('relu'))
    model.add(Convolution2D(nb_filters, nb_conv, nb_conv))
    model.add(Activation('relu'))
    model.add(MaxPooling2D(pool_size=(nb_pool, nb_pool)))
    model.add(Dropout(0.25))
    
    model.add(Convolution2D(nb_filters, nb_conv, nb_conv))
    model.add(Activation('relu'))
    model.add(Convolution2D(nb_filters, nb_conv, nb_conv))
    model.add(Activation('relu'))
    model.add(MaxPooling2D(pool_size=(nb_pool, nb_pool)))
    model.add(Dropout(0.25))
    

    model.add(Flatten())
    model.add(Dense(128))
    model.add(Activation('relu'))
    model.add(Dropout(0.5))
    model.add(Dense(nb_classes))
    model.add(Activation('softmax'))

    model.compile(loss='categorical_crossentropy', optimizer='Adam') # Adam adadelta
    
    if weights_path:
        model.load_weights(weights_path)
        
    return model