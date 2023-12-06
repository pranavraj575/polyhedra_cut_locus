def hex2(d):
    out = hex(d)[2:]
    while len(out) < 2:
        out = '0' + out
    return out


def pix_hex(r, g, b):
    return '#' + hex2(r) + hex2(g) + hex2(b)


def rb_gradient(n):
    if n == 0:
        return []
    if n == 1:
        return [pix_hex(127, 0, 127)]
    step = 255/(n - 1)
    stuff = [round(step*i) for i in range(n)]
    return [pix_hex(255 - b, 0, b) for b in stuff]