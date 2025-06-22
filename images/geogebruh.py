# cropping whitespace border out of images

if __name__ == '__main__':
    import argparse
    import numpy as np
    import os
    from PIL import Image

    PARSER = argparse.ArgumentParser()
    PARSER.add_argument('--color', type=int, nargs=3, required=False, default=[255, 255, 255],
                        help='definition of \'white\' pixel')
    PARSER.add_argument('--border-size', type=int, default=1, required=False,
                        help='whitespace to add to/keep on border')
    PARSER.add_argument('--img-file', nargs='*', default=[], required=False,
                        help='file(s) to remove whitespace border of')
    args = PARSER.parse_args()

    files = args.img_file
    border_size = args.border_size
    white = np.array(args.color, dtype=np.uint8)
    for file in files:
        print('cropping',file)
        img = np.asarray(Image.open(file))

        for idx in (0, -1):
            while np.all(img[idx, :, :3] == white.reshape(1, 1, -1)):
                img = np.concatenate((img[:idx], img[(idx%len(img)) + 1:]), axis=0)

            while np.all(img[:, idx, :3] == white.reshape(1, 1, -1)):
                img = np.concatenate((img[:, :idx], img[:, (idx%len(img[0])) + 1:]), axis=1)
        if border_size > 0:
            out_img = np.ones((img.shape[0] + border_size*2, img.shape[1] + border_size*2, img.shape[2]),
                              dtype=img.dtype)*255
            out_img[:, :, :3] = white.reshape(1, 1, -1)
            out_img[border_size:-border_size, border_size:-border_size, :] = img
        else:
            out_img = img
        os.remove(file)
        Image.fromarray(out_img).save(file)
