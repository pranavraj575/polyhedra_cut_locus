# cropping whitespace border out of images

if __name__ == '__main__':
    import argparse
    import numpy as np
    import os
    from PIL import Image

    PARSER = argparse.ArgumentParser()
    PARSER.add_argument('--color', type=int, nargs="+", required=False, default=[255, 255, 255],
                        help='definition of \'white\' pixel, "-1" means ignore this channel')
    PARSER.add_argument('--replacement', type=int, nargs="+", required=False, default=None,
                        help='pixel to set border as (defaults to --color)')
    PARSER.add_argument('--border-size', type=int, default=1, required=False,
                        help='whitespace to add to/keep on border')
    PARSER.add_argument('--img-file', nargs='*', default=[], required=False,
                        help='file(s) to remove whitespace border of')
    args = PARSER.parse_args()

    files = list(args.img_file)
    border_size = args.border_size
    channels = [i for i in range(len(args.color)) if args.color[i] >= 0]
    white = np.array([c for c in args.color if c >= 0],
                     dtype=np.uint8)

    img_files = []
    i = 0
    suffixes = ['.png', '.jpg', '.jpeg']
    while i < len(files):
        thing = files[i]
        if os.path.isdir(thing):
            for c in os.listdir(thing):
                files.append(os.path.join(thing, c))
        elif any(thing.endswith(suffix) for suffix in suffixes):
            img_files.append(thing)
        i += 1

    for thing in img_files:
        print('cropping', thing)
        img = np.asarray(Image.open(thing))
        for idx in (0, -1):
            while np.all(img[idx, :][:, channels] == white.reshape(1, -1)):
                img = np.concatenate((img[:idx], img[(idx%len(img)) + 1:]), axis=0)

            while np.all(img[:, idx][:, channels] == white.reshape(1, 1, -1)):
                img = np.concatenate((img[:, :idx], img[:, (idx%len(img[0])) + 1:]), axis=1)
        if border_size > 0:
            out_img = np.ones((img.shape[0] + border_size*2, img.shape[1] + border_size*2, img.shape[2]),
                              dtype=img.dtype)*255

            if args.replacement is not None:
                replacement = np.array(args.replacement, dtype=img.dtype)
            else:
                replacement = np.zeros(img.shape[2], dtype=img.dtype)
                replacement[channels] = white
            out_img[:, :, :] = replacement.reshape(1, 1, -1)
            out_img[border_size:-border_size, border_size:-border_size, :] = img
        else:
            out_img = img
        os.remove(thing)
        Image.fromarray(out_img).save(thing)
