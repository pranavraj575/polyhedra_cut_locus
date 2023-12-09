# cut locus plot for a point on S^2

from matplotlib import pyplot as plt
import numpy as np
import os

TOL = .003


def line_into_segs(a, b):
    pass


def plot_bounds(nine=False):
    low, high = (-1, 2) if nine else (0, 1)
    linewidth = 4
    alpha = .4
    color = 'black'
    for x in (0, 1):
        plt.plot((x, x), (low, high), '--', color=color, linewidth=linewidth, alpha=alpha)
        plt.plot((low, high), (x, x), '--', color=color, linewidth=linewidth, alpha=alpha)
    length = .05
    ht = length/2

    plt.plot((.5 - length/2, .5 + length/2), (1 - ht, 1), '-', color=color, linewidth=linewidth, alpha=alpha)
    plt.plot((.5 - length/2, .5 + length/2), (1 + ht, 1), '-', color=color, linewidth=linewidth, alpha=alpha)

    plt.plot((.5 - length/2, .5 + length/2), (ht, 0), '-', color=color, linewidth=linewidth, alpha=alpha)
    plt.plot((.5 - length/2, .5 + length/2), (-ht, 0), '-', color=color, linewidth=linewidth, alpha=alpha)

    offset = length/2

    plt.plot((1 - ht, 1), (.5 - length/2 + offset, .5 + length/2 + offset), '-', color=color, linewidth=linewidth, alpha=alpha)
    plt.plot((1 + ht, 1), (.5 - length/2 + offset, .5 + length/2 + offset), '-', color=color, linewidth=linewidth, alpha=alpha)

    plt.plot((1 - ht, 1), (.5 - length/2 - offset, .5 + length/2 - offset), '-', color=color, linewidth=linewidth, alpha=alpha)
    plt.plot((1 + ht, 1), (.5 - length/2 - offset, .5 + length/2 - offset), '-', color=color, linewidth=linewidth, alpha=alpha)

    plt.plot((-ht, 0), (.5 - length/2 + offset, .5 + length/2 + offset), '-', color=color, linewidth=linewidth, alpha=alpha)
    plt.plot((+ht, 0), (.5 - length/2 + offset, .5 + length/2 + offset), '-', color=color, linewidth=linewidth, alpha=alpha)

    plt.plot((-ht, 0), (.5 - length/2 - offset, .5 + length/2 - offset), '-', color=color, linewidth=linewidth, alpha=alpha)
    plt.plot((+ht, 0), (.5 - length/2 - offset, .5 + length/2 - offset), '-', color=color, linewidth=linewidth, alpha=alpha)


def geodesic_counts_plain(points, consider, tol=TOL):
    # points is nx2
    # points is points to calc the number of geodesics
    # consider is 9?x2
    # consider is points to consider (i.e., if it is (S1)^2 torus, copy the unit square over
    # returns n x 1 number of geodesics each point has from the consideration points
    bol = which_shortest_paths(points, consider, tol)
    return np.expand_dims(np.sum(bol, axis=1), 1)


def geodesic_counts_complicated(points, consider, consider_obs, obs_r, tol=TOL):
    # points is n x 2
    # consider is m x 2
    # consider_obs is 9 x 2
    # assumes every point intersects a obstacle
    # geodesic counts, redoes the complicated ones

    easy_dists = paths_dist(points, consider)
    # n x m

    shortest_paths = which_shortest_paths(points, consider, tol)
    # n x m

    points = np.expand_dims(points, 1)
    # n x 1 x 2
    consider = np.expand_dims(consider, 0)
    # 1 x m x 2

    n, m = easy_dists.shape
    # n x m

    obs_that_each_intersects = which_paths_intersect_which_obs(points, consider, consider_obs, obs_r)
    # n x m x 9
    redo_indices = np.where(obs_that_each_intersects)
    # each of these is a k array
    # k is the number of times (point,consider_pt) intersects with an object
    # TODO: we also assume each (pt,consider_pt) intersects w at most 1 obj

    k = len(redo_indices[0])

    redo_pts = points[redo_indices[0], :, :].reshape((k, 2))
    # k x 2
    # redo_pts_abs=np.linalg.norm(redo_pts,axis=-1)

    redo_consider = consider[:, redo_indices[1], :].reshape((k, 2))
    # k x 2
    # redo_consider_abs=np.linalg.norm(redo_consider,axis=-1)

    redo_obs = consider_obs[redo_indices[2], :]
    # k x 2

    to_pt = redo_pts - redo_obs
    to_pt_dist = np.linalg.norm(to_pt, axis=-1)
    to_cons = redo_consider - redo_obs
    to_cons_dist = np.linalg.norm(to_cons, axis=-1)

    alphas = np.arccos(np.matmul(np.expand_dims(to_pt, -2),
                                 np.expand_dims(to_cons, -1)).reshape(k)/(
                               to_pt_dist*to_cons_dist
                       ))
    worse_alphas = 2*np.pi - alphas
    # k dimensional smallest angles between points and obstacle center

    thetas = np.arcsin(obs_r/to_pt_dist)
    theta_primes = np.arcsin(obs_r/to_cons_dist)
    # k dimensional angle to tangent points of circle
    arc_angles = alphas - np.pi + thetas + theta_primes

    worse_arc_angles = worse_alphas - np.pi + thetas + theta_primes

    dists = np.cos(thetas)*to_pt_dist
    dists_prime = np.cos(theta_primes)*to_cons_dist
    # k dimensional distance from pt to tangeant of obstacle

    new_dists = dists + dists_prime + arc_angles*obs_r
    # k dimensional shortest distances

    worse_dists = dists + dists_prime + worse_arc_angles*obs_r
    temp = np.full((n, m), np.inf)
    temp[redo_indices[0], redo_indices[1]] = worse_dists
    # n x m, should have at least one non inf value per row
    temp = np.min(temp, axis=-1)
    # n array

    easy_dists[redo_indices[0], redo_indices[1]] = new_dists

    complicated_dists = np.concatenate((easy_dists, np.expand_dims(temp, 1)), axis=1)
    # n x (m+1), appends the extra 'longer' way to go around circle
    diffs = np.expand_dims(np.min(complicated_dists, axis=1), 1)

    bol = (complicated_dists - diffs) <= tol
    return np.expand_dims(np.sum(bol, axis=1), 1)


def which_to_reconsider(points, consider, consider_obs, obs_r):
    # points n x 2
    # returns which points to considre (n boolean array)
    stuff = np.logical_and(which_intersect_obs(points, consider, consider_obs, obs_r),
                           which_shortest_paths(points, consider, tol=TOL))
    return np.sum(stuff, axis=1) > 0


def which_paths_intersect_which_obs(a, b, consider_obs, obs_r):
    # b is "consder", should be 1 x m x 2
    # a is grid points, should be n x 1 x 2
    # a,b broadcast to n x m x 2 (m probably 9)
    # consider obs is 9 x 2
    # nxmx9 array of which a->b overlaps with which obs
    a = np.expand_dims(a, 2)
    b = np.expand_dims(b, 2)
    # broadcast to n x m x 1 x 2

    consider_obs = np.expand_dims(consider_obs, 0)
    consider_obs = np.expand_dims(consider_obs, 0)
    # 1 x 1 x 9 x 2

    pt_to_obs = consider_obs - b
    # 1 x b x 9 x 2

    d_obs_to_pt = np.linalg.norm(pt_to_obs, axis=-1)
    c = np.sqrt(d_obs_to_pt**2 - obs_r**2)*np.cos(np.arcsin(obs_r/d_obs_to_pt))
    # 1 x m x 9

    blast = np.matmul(np.expand_dims(pt_to_obs/np.expand_dims(np.linalg.norm(pt_to_obs, axis=-1), -1), -2),
                      np.expand_dims(a - b, -1))
    blast = blast.reshape(blast.shape[:3])
    # n x m x 9

    # middle=np.matmul(np.expand_dims(consider_obs-a,-2),np.expand_dims(b-consider_obs,-1))
    # middle=middle.reshape(middle.shape[:3])
    # n x m x 9

    res = np.abs(np.cross(b - a, a - consider_obs, axis=-1))/np.linalg.norm(b - a, axis=-1)
    # n x m x 9
    # intersets= np.logical_and(res<obs_r, middle>0)
    intersets = np.logical_and(res < obs_r, blast > c)
    return intersets


def which_paths_intersect_obs(a, b, consider_obs, obs_r):
    # b is "consder", should be 1 x m x 2
    # a is grid points, should be n x 1 x 2
    # a,b broadcast to n x m x 2 (m probably 9)
    # consider obs is 9 x 2
    # nxm array of which a->b overlaps with obs
    intersects = which_paths_intersect_which_obs(a, b, consider_obs, obs_r)
    return np.sum(intersects, axis=-1)


def which_intersect_obs(points, consider, consider_obs, obs_r):
    # points is nx2
    # consider is 9x2
    # consider obs is 9x2
    # nx9 array of which lines overlap with obs
    points = np.expand_dims(points, 1)
    # n x 1 x 2
    consider = np.expand_dims(consider, 0)
    # 1 x 9 x 2
    intersects = which_paths_intersect_obs(points, consider, consider_obs, obs_r)
    return intersects


def paths_dist(points, consider):
    # points is nx2
    # points is points to calc the number of geodesics
    # consider is mx2
    # consider is points to consider (i.e., if it is (S1)^2 torus, copy the unit square over
    # nxm array of closest straight dist to each point
    points = np.expand_dims(points, 1)
    # n x 1 x 2
    consider = np.expand_dims(consider, 0)
    # 1 x m x 2
    actual = np.linalg.norm(points - consider, axis=2)
    # n x m
    return actual


def which_shortest_paths(points, consider, tol=TOL):
    # points is nx2
    # points is points to calc the number of geodesics
    # consider is mx2
    # consider is points to consider (i.e., if it is (S1)^2 torus, copy the unit square over
    # nxm array of which consideration points have the closest distance to each point
    actual = paths_dist(points, consider)
    diffs = np.expand_dims(np.min(actual, axis=1), 1)
    return actual - diffs <= tol


def break_into_segs(a, b):
    def break_i(a, b, i):
        if a[i] > b[i]:
            aa = b
            bb = a
        else:
            aa = a
            bb = b
        # this ensures aa[i]<=bb[i]
        ai = aa[i]
        bi = bb[i]
        if bi <= np.floor(ai) + 1:  # if ai and bi are on the same integer interval
            yield (aa, bb)
            return
        v = bb - aa
        v = v/v[i]  # defined since bi>ai, so v[i] non-zero
        partition = [j for j in range(int(np.floor(ai) + 1), int(np.ceil(bi)))]  # how to partition this range into unit intervals
        partition.append(bi)
        last_val = ai
        last_point = aa
        for val in partition:
            diff = val - last_val
            point = last_point + v*diff
            yield (last_point, point)
            last_val = val
            last_point = point

    for (aa, bb) in break_i(a, b, 0):
        for (aaa, bbb) in break_i(aa, bb, 1):
            yield (aaa, bbb)


def plot_line(a, b, color='red'):
    for (aa, bb) in break_into_segs(a, b):
        shift = np.array((0, 0))
        for v in (aa, bb):
            vv = v + shift
            for i in range(2):
                while vv[i] < 0:
                    vv[i] += 1
                    shift[i] += 1
                while vv[i] > 1:
                    vv[i] -= 1
                    shift[i] -= 1
        aa = aa + shift
        bb = bb + shift
        plt.plot([aa[0], bb[0]], [aa[1], bb[1]], color=color)


p = np.array([0.97021987, 0.082961])
obstruction_pt = np.array([0.12157324, 0.40209816])
# obstruction_pt = np.array([0.221157324, 0.40209816])
# obstruction_pt = np.array([0.25, 0.40209816])
# obstruction_pt = np.array([0.27, 0.40209816])
# obstruction_pt = np.array([0.28, 0.40209816])
# obstruction_pt = np.array([0.29, 0.40209816])
# obstruction_pt = np.array([0.295, 0.40209816])
p = np.array([.3, .3])
obstruction_pt = np.array([0.5, 0.5])

# obstruction_pt = np.array([0.52157324, 0.40209816])
# p = np.array([0.2, 0.1])
# obstruction_pt = np.array([0.2, 0.2])

# p = np.random.random(2)
# obstruction_pt = np.random.random(2)
obstruction_r = .15
while np.linalg.norm(obstruction_pt - p) < obstruction_r:
    obstruction_pt = np.random.random(2)
print(p)
print(obstruction_pt)

fig, ax = plt.subplots()
ax.add_patch(plt.Circle(obstruction_pt, obstruction_r))

consider_obs = []
consider = []
for x in range(-1, 2):
    for y in range(-1, 2):
        consider.append(p + [x, y])
        consider_obs.append(obstruction_pt + [x, y])

consider = np.array(consider)
consider_obs = np.array(consider_obs)

nx, ny = (500, 500)

# 9 x 2
x = np.linspace(0, 1, nx)
y = np.linspace(0, 1, ny)
xv, yv = np.meshgrid(x, y)
xv = xv.flatten()
yv = yv.flatten()
vv = np.stack((xv, yv), axis=1)
vv = vv + (1/(nx + ny))*(np.random.random(vv.shape) - .5)
bb_dist = np.min(np.linalg.norm(np.expand_dims(vv, 1) - np.expand_dims(consider_obs, 0), axis=-1), axis=-1)
# n , distance to closest obstacle copy
vv = vv[np.where(bb_dist >= obstruction_r)[0], :]

annoying_indices_bool = which_to_reconsider(vv, consider, consider_obs, obstruction_r)
nice_grid = vv[np.where(np.logical_not(annoying_indices_bool))[0]]
gg = geodesic_counts_plain(nice_grid, consider, tol=TOL)

largest = np.max(gg)

cols = ['purple', 'blue', 'green', 'red']
seen = []
for k in range(2, largest + 1):
    vals = nice_grid[np.where(gg == k)[0], :]
    if len(vals) > 0:
        label = None
        if not k in seen:
            seen.append(k)
            label = str(k) + ' geodesics'
        plt.scatter(vals[:, 0], vals[:, 1], color=cols[(k - 2)%len(cols)], label=label)
plt.scatter(p[:1], p[1:], s=40, color='black')

annoying_vals = vv[np.where(annoying_indices_bool)[0], :]
max_geodesic = None
if len(annoying_vals > 0):
    gg2 = geodesic_counts_complicated(annoying_vals, consider, consider_obs, obstruction_r)
    largest2 = np.max(gg2)
    for k in range(2, largest2 + 1):
        vals = annoying_vals[np.where(gg2 == k)[0], :]
        if len(vals) > 0:
            label = None
            max_geodesic = vals[0], k
            if not k in seen:
                seen.append(k)
                label = str(k) + ' geodesics'
            plt.scatter(vals[:, 0], vals[:, 1], color=cols[(k - 2)%len(cols)], label=label)

# plt.scatter(annoying_vals[:,0],annoying_vals[:,1],color='red',alpha=.5)

plot_bounds()
plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=4, )
plot_line(p, max_geodesic[0] - [1, 0])

plt.savefig(os.path.join("images", "hole_torus", "pq.png"))
plt.show()
