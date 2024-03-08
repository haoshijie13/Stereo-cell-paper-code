import stardist
import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import skimage.io as io
import argparse
import PIL
PIL.Image.MAX_IMAGE_PIXELS = 933120000

parser = argparse.ArgumentParser(
                    prog='runStarDist',
                    description='detect cells in a dapi staining image.',
                    epilog='')
parser.add_argument('-i', '--input',dest="input",help="The input dapi staining image.")      # option that takes a value
parser.add_argument('-o', '--output',dest="output",help="The output labeling table.")      # option that takes a value
parser.add_argument('-m', '--mask',dest="mask",help="The output labeling image.")      # option that takes a value

args = parser.parse_args()

im=io.imread(args.input)

from stardist.models import StarDist2D

# prints a list of available models
StarDist2D.from_pretrained()

# creates a pretrained model
model = StarDist2D.from_pretrained('2D_versatile_fluo')

from stardist.plot import render_label
from csbdeep.utils import normalize
import matplotlib.pyplot as plt

img = im #[5000:7000,5000:7000]
#img = im[1:7000,1:7000]

labels, _ = model.predict_instances_big(normalize(img),axes='YX', block_size=2048,
                                                   sparse=True, min_overlap=192, context=128, n_tiles=(8,8))

tab=pd.DataFrame(labels).reset_index().melt('index')
tab.columns=("y","x","value")
tab=tab[tab.value!=0]
tab.to_csv(args.output,index=False)
#io.imsave(arr=labels,fname=args.mask,compression="zlib")
