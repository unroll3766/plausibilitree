"""
treesim.py: Functions for tree model simulation and sampling
Author: Mitch Bryson

Requires NumPy and Trimesh
"""

import os, sys
from math import cos, pi, sin, sqrt, pow
import trimesh
import numpy as np


def _seed_trimesh_sampling():
    """Make trimesh's surface sampling follow numpy's global seed.

    trimesh >= 5 samples mesh surfaces with its own unseeded generator, which
    ignores ``np.random.seed`` -- so the same seed would build the same tree
    skeleton but scatter different surface points on every run. Wrap
    ``Trimesh.sample`` to draw its ``seed`` from the global RNG, which
    the pipeline seeds per tree. treesim's own geometry code stays on
    ``np.random`` throughout, so with this wrapper a (model, seed) pair
    always produces the same cloud. Idempotent.
    """
    orig = trimesh.Trimesh.sample
    if getattr(orig, "_ssf_seeded", False):
        return

    def seeded(self, count, *args, **kw):
        # Derive the sampler seed from the global MT19937 state WITHOUT
        # consuming a draw: treesim samples the stem mesh before it draws
        # branch geometry off the same global stream, so a consuming draw
        # here would shift the branch RNG and change skeletons. The
        # state-derived seed (buffer word at the current position) still
        # differs call to call, keeping surface sampling deterministic and
        # varied while skeleton bits stay exactly as they were.
        if "seed" not in kw:
            s = np.random.get_state()
            kw["seed"] = (int(s[1][int(s[2]) % len(s[1])]) & 0x7FFFFFFF) or 1
        return orig(self, count, *args, **kw)

    seeded._ssf_seeded = True
    trimesh.Trimesh.sample = seeded


_seed_trimesh_sampling()

# get_spline_params: get spline parameters from points
def get_spline_params(p0, p1, p2):
    A = np.array([ [2/(p1[0]-p0[0]), 1/(p1[0]-p0[0]), 0], \
        [1/(p1[0]-p0[0]), 2*((1/(p1[0]-p0[0]))+(1/(p2[0]-p1[0]))), 1/(p2[0]-p1[0])], \
        [0, 1/(p2[0]-p1[0]), 2/(p2[0]-p1[0])] ])
    b = np.array([3*(p1[1]-p0[1])/pow(p1[0]-p0[0],2), \
        3*( (p1[1]-p0[1])/pow(p1[0]-p0[0],2) + (p2[1]-p1[1])/pow(p2[0]-p1[0],2) ), \
        3*(p2[1]-p1[1])/pow(p2[0]-p1[0],2)])
    k = np.linalg.solve(A, b)
    return k

# get_spline_points: get points along a spline defined by points and parameters kx, ky
def get_spline_points(p_spline,kx,ky,z):
	
	x0 = p_spline[0][0]
	y0 = p_spline[0][1]
	z0 = p_spline[0][2]
	x1 = p_spline[1][0]
	y1 = p_spline[1][1]
	z1 = p_spline[1][2]
	x2 = p_spline[2][0]
	y2 = p_spline[2][1]
	z2 = p_spline[2][2]
	
	a01x = kx[0]*(z1-z0) - (x1-x0);
	a12x = kx[1]*(z2-z1) - (x2-x1);
	a01y = ky[0]*(z1-z0) - (y1-y0);
	a12y = ky[1]*(z2-z1) - (y2-y1);
	
	b01x = -kx[1]*(z1-z0) - (x1-x0);
	b12x = -kx[2]*(z2-z1) - (x2-x1);
	b01y = -ky[1]*(z1-z0) - (y1-y0);
	b12y = -ky[2]*(z2-z1) - (y2-y1);
	
	t01 = (z-z0)/(z1-z0)
	t12 = (z-z1)/(z2-z1)
	
	x01out = x0*(1-t01) + x1*t01 + np.multiply(np.multiply(t01,(1-t01)),(a01x*(1-t01) + b01x*t01))
	x12out = x1*(1-t12) + x2*t12 + np.multiply(np.multiply(t12,(1-t12)),(a12x*(1-t12) + b12x*t12))
	y01out = y0*(1-t01) + y1*t01 + np.multiply(np.multiply(t01,(1-t01)),(a01y*(1-t01) + b01y*t01))
	y12out = y1*(1-t12) + y2*t12 + np.multiply(np.multiply(t12,(1-t12)),(a12y*(1-t12) + b12y*t12))
	
	xout = np.multiply(x01out,(z < z1)) + np.multiply(x12out,(z >= z1))
	yout = np.multiply(y01out,(z < z1)) + np.multiply(y12out,(z >= z1))
	
	return (xout,yout)

# _split_stem: truncate a stem at a drawn height and grow two daughter stems
# from the cut -- the generalisation of the original single-fork construction.
# The split height is drawn (split_height_range fraction) over the stem's own
# span; the daughters separate around the parent stem's tip by 1.5-3.5 m
# scaled by `sep_scale` (1.0 = the original single-fork distance); one
# daughter keeps the parent's tip height and the other is shortened; both
# bend outward through mid control points.
# Returns (truncated_parent, daughter1, daughter2, split_height).
def _split_stem(stem, splhr, sep_scale=1.0):
	(xc, yc, zc, p_spline, kx, ky, zrange) = stem
	zbase, ztop = zrange
	frac = splhr[0] + (splhr[1]-splhr[0])*np.random.random()
	frac = min(max(frac, 0.05), 0.95)
	z_split = zbase + frac*(ztop - zbase)

	# truncated parent: same spline, cut at z_split
	zc_t = np.linspace(zbase, z_split, num=20)
	(xc_t, yc_t) = get_spline_points(p_spline, kx, ky, zc_t)
	truncated = [xc_t, yc_t, zc_t, p_spline, kx, ky, [zbase, z_split]]

	# the parent stem's own tip (from its spline) anchors the daughters
	(tipx, tipy) = get_spline_points(p_spline, kx, ky, np.array([ztop]))
	tip = [float(tipx[0]), float(tipy[0]), ztop]

	dr = (2.0*np.random.random()+1.5)*sep_scale
	theta = 2*pi*np.random.random()
	dx = dr*sin(theta)
	dy = dr*cos(theta)
	topxy1 = [tip[0]+dx, tip[1]+dy, ztop]
	l2 = (0.3*np.random.random()+0.7)*(ztop - z_split)+z_split
	topxy2 = [tip[0]-dx, tip[1]-dy, l2]
	bottomxy = [float(xc_t[-1]), float(yc_t[-1]), z_split]

	bend = 0.3+0.1*np.random.random()
	mid1 = [bottomxy[0]+0.3*(topxy1[0]-bottomxy[0]),
	        bottomxy[1]+0.3*(topxy1[1]-bottomxy[1]),
	        bend*(topxy1[2]-z_split)+z_split]
	mid2 = [bottomxy[0]+0.3*(topxy2[0]-bottomxy[0]),
	        bottomxy[1]+0.3*(topxy2[1]-bottomxy[1]),
	        bend*(topxy2[2]-z_split)+z_split]

	out = []
	for (mid_d, top) in ((mid1, topxy1), (mid2, topxy2)):
		ps = np.array([bottomxy, mid_d, top])
		nkx = get_spline_params( (ps[0][2],ps[0][0]), (ps[1][2],ps[1][0]), (ps[2][2],ps[2][0]))
		nky = get_spline_params( (ps[0][2],ps[0][1]), (ps[1][2],ps[1][1]), (ps[2][2],ps[2][1]))
		zc_d = np.linspace(z_split, top[2], num=20)
		(xdc, ydc) = get_spline_points(ps, nkx, nky, zc_d)
		out.append([xdc, ydc, zc_d, ps, nkx, nky, [z_split, top[2]]])
	return truncated, out[0], out[1], z_split

# gen_single_simtree: generate synthetic trees
def gen_single_simtree(Nstem=300, Nfol=1000, model_params=None, tx=None, ty=None,
                       return_skeleton=False):
    
    # class IDs used for different point types
    class_id = {}
    class_id['folliage'] = 0
    class_id['stem'] = 1
    
    # Generate main stem parameters
    if tx is None:
        tx = 0
    if ty is None:
        ty = 0
    tz = 0
    hr = model_params['height_range']
    h = hr[0] + (hr[1]-hr[0])*np.random.random()
    dr = model_params['diam_range']
    r = 0.5*dr[0] + 0.5*(dr[1]-dr[0])*np.random.random()
    
    # randomise fork status. split_prob gates whether the tree forks at all;
    # num_splits (range [lo, hi], default [1, 1]) is how many cascading
    # splits a forked tree gets: the trunk splits, then one of the live stems
    # splits again, and so on. 1 split = the original single-fork tree.
    rn = np.random.random()
    n_splits = 0
    if rn < model_params['split_prob']:
        snr = model_params.get('num_splits', [1, 1])
        lo, hi = int(snr[0]), int(snr[1])
        if hi > lo:
            lo = lo + int(np.random.random()*(hi - lo + 1))
        n_splits = max(1, lo)
    fork = 2 if n_splits >= 1 else 1   # legacy field: 1 = single stem, 2 = forked

    stem_data = []

    # Initial spline for main stem
    p0 = [tx,ty,tz]
    sweeptopx = 2*model_params['tree_top_dist']*np.random.random()-model_params['tree_top_dist']
    sweeptopy = 2*model_params['tree_top_dist']*np.random.random()-model_params['tree_top_dist']
    tmdist = model_params['tree_mid_dist']
    p1 = [tx+2*tmdist*np.random.random()-tmdist+0.2*sweeptopx,ty+2*tmdist*np.random.random()-tmdist+0.2*sweeptopy,0.5*h+(0.167*h*np.random.random())]
    p2 = [tx+sweeptopx,ty+sweeptopy,h]
    p_spline = np.array([p0, p1, p2])
    kx = get_spline_params( (p_spline[0][2],p_spline[0][0]), (p_spline[1][2],p_spline[1][0]), (p_spline[2][2],p_spline[2][0]))
    ky = get_spline_params( (p_spline[0][2],p_spline[0][1]), (p_spline[1][2],p_spline[1][1]), (p_spline[2][2],p_spline[2][1]))

    split_log = []
    zc = np.linspace(tz,h,num=20)
    (xc,yc) = get_spline_points(p_spline,kx,ky,zc)
    stem_data.append([xc,yc,zc,p_spline,kx,ky,[0,h]])

    if n_splits >= 1:
        # Split/forked stem: n_splits cascading splits. Each split truncates
        # one live stem at a height drawn over that stem's own span (the same
        # split_height_range fraction) and grows two daughters from the cut.
        splhr = model_params['split_height_range']
        # The main fork separates its daughters by the full 1.5-3.5 m. Later
        # cascading splits are subordinate structure: their separation scales
        # with the stem's remaining span relative to the first split AND is
        # damped (0.6), so the divergence slope stays well under the main
        # fork's rate. The first split's scale is exactly 1.0, so single-fork
        # trees are unchanged draw for draw.
        ref_span = None
        for k in range(n_splits):
            # Leader cascade: always split the tallest live stem, so the
            # forks climb the dominant leader and every crown forms in the
            # upper canopy, interlocked with its neighbours. (Splitting a
            # random or low stem sends a shortened daughter sideways, whose
            # crown then hangs below and beside the main crown as a
            # detached blob.)
            eligible = [j for j, s in enumerate(stem_data)
                        if (s[6][1] - s[6][0]) > 0.15*h]
            if not eligible:
                break                     # nothing left tall enough to split
            j = max(eligible,
                    key=lambda jj: stem_data[jj][6][1] - stem_data[jj][6][0])
            span = stem_data[j][6][1] - stem_data[j][6][0]
            if ref_span is None:
                ref_span = span
            sep_scale = min(1.0, span/ref_span)
            if k >= 1:
                sep_scale *= 0.6
            stem_data[j], d1, d2, z_split = _split_stem(stem_data[j], splhr,
                                                        sep_scale)
            split_log.append({'split_index': int(k), 'parent_polyline': int(j),
                              'height_m': float(z_split)})
            stem_data.append(d1)
            stem_data.append(d2)

    # Branch budget: the original single fork doubled the upper-crown band
    # (fork = 2 crowns). Keep that level for ANY number of cascading splits --
    # the point budget (Np) does not grow with the number of leaders, so
    # scaling the branch count up would spread the same points over more,
    # sparser ribbons and visibly fragment the canopy. More leaders therefore
    # share the forked-tree budget; raise Np for denser multi-leader crowns.
    fork_mult = 2 if n_splits >= 1 else 1
        
    # Generate stem mesh
    verts = []
    faces = []
    lv = 0
    for stem in stem_data:
        (xc,yc,zc,p,kx,ky,zrange) = stem
        Ndotsl = 12
        t = np.linspace(0,2*pi,Ndotsl+1)
        t = t[:-1]
        for i in range(len(zc)):
            rn = r*(1-(zc[i]/h))
            xn = rn*np.sin(t)+xc[i]
            yn = rn*np.cos(t)+yc[i]
            for j in range(xn.shape[0]):
                verts.append([xn[j],yn[j],zc[i]])
        for i in range(len(zc)-1):
            for j in range(t.shape[0]-1):
                faces.append([Ndotsl*i+j+1+lv,Ndotsl*i+j+lv,Ndotsl*(i+1)+j+lv,Ndotsl*(i+1)+j+1+lv])
            faces.append([Ndotsl*i+lv,Ndotsl*i+t.shape[0]-1+lv,Ndotsl*(i+1)+t.shape[0]-1+lv,Ndotsl*(i+1)+lv])
        lv = len(verts)
    
    mesh = trimesh.Trimesh(vertices=verts,faces=faces)
    
    # Sample points from stem meshes
    points = mesh.sample(Nstem)
    points_all = np.concatenate((points,class_id['stem']*np.ones((Nstem,1))),axis=1)
    
    # Generate branches and foliage
    h1r = [model_params['min_can_height'][0], model_params['min_can_height'][1]]
    h2r = [model_params['max_can_width_height'][0], model_params['max_can_width_height'][1]]
    h1 = h*(h1r[0]+(h1r[1]-h1r[0])*np.random.random())
    h2 = h*(h2r[0]+(h2r[1]-h2r[0])*np.random.random())
    
    branchz1 = h1+(h2-h1)*np.random.random(np.random.randint(0.166*model_params['num_branches'][0],0.3*model_params['num_branches'][1]))
    branchz2 = h2+(h-h2)*np.random.random(np.random.randint(fork_mult*0.833*model_params['num_branches'][0],fork_mult*0.7*model_params['num_branches'][1]))
    
    branchz = np.concatenate((branchz1,branchz2))
    xyc = []
    for stem in stem_data:
        (x,y,z,p_spline,kx,ky,zrange) = stem
        (xc,yc) = get_spline_points(p_spline,kx,ky,branchz)
        xyc.append((xc,yc))
    
    # Generate mesh data for branches and foliage
    verts = []
    faces = []
    branch_skeleton = []
    built_z = []      # base height of each branch that actually got a ribbon
    c = 0
    rolloff = 1.0
    rolloffh = h2
    rolloffh2 = h1
    for i in range(len(branchz)):
        
        # determine which stem to put on
        okstems = []
        for (j,stem) in enumerate(stem_data):
            (x,y,z,p_spline,kx,ky,zrange) = stem
            if branchz[i] >= zrange[0] and branchz[i] <= zrange[1]:
                okstems.append(j)
        if len(okstems) == 0:
            pass
        elif len(okstems) == 1:
            j = okstems[0]
            xc = xyc[j][0]
            yc = xyc[j][1]
        else:
            j = okstems[np.random.choice(len(okstems),1)[0]]
            xc = xyc[j][0]
            yc = xyc[j][1]
        
        if (i % 10) == 0:
            gap = np.random.random()
            gap2 = (2*pi-gap)*np.random.random()
        theta = 2*pi*np.random.random()
        if theta > gap2 and theta < (gap2+gap):
            continue
        rn = r*(1-(branchz[i]/h))
        rout = (model_params['max_can_width']/7.0)*sqrt(h-branchz[i])
        if branchz[i] < rolloffh2:
            rout = 0
        elif branchz[i] < rolloffh:
            rout = rout-(rolloff*rout*((rolloffh-branchz[i])/rolloffh))
        
        rout2 = rout+0.4*rout*np.random.random()-0.2*rout
        thick = 0.2*rout2
        verts.append([rn*sin(theta)+xc[i],rn*cos(theta)+yc[i],branchz[i]])
        verts.append([(rout2/3)*sin(theta)+thick*cos(theta)+xc[i],(rout2/3)*cos(theta)-thick*sin(theta)+yc[i],branchz[i]+np.random.random()])
        verts.append([(rout2/3)*sin(theta)-thick*cos(theta)+xc[i],(rout2/3)*cos(theta)+thick*sin(theta)+yc[i],branchz[i]+np.random.random()])
        verts.append([(2*rout2/3)*sin(theta)+thick*cos(theta)+xc[i],(2*rout2/3)*cos(theta)-thick*sin(theta)+yc[i],branchz[i]+np.random.random()])
        verts.append([(2*rout2/3)*sin(theta)-thick*cos(theta)+xc[i],(2*rout2/3)*cos(theta)+thick*sin(theta)+yc[i],branchz[i]+np.random.random()])
        verts.append([rout2*sin(theta)+xc[i],rout2*cos(theta)+yc[i],branchz[i]+np.random.random()])
        
        branch_skeleton.append({
            'base': [float(verts[-6][0]), float(verts[-6][1]), float(verts[-6][2])],
            'tip': [float(verts[-1][0]), float(verts[-1][1]), float(verts[-1][2])],
            'height_m': float(branchz[i]),
            'azimuth_rad': float(theta),
            'length_m': float(np.linalg.norm(np.array(verts[-1]) - np.array(verts[-6]))),
        })
        built_z.append(branchz[i])
        
        faces.append([6*c,6*c+1,6*c+2])
        faces.append([6*c+1,6*c+2,6*c+4])
        faces.append([6*c+1,6*c+4,6*c+3])
        faces.append([6*c+3,6*c+4,6*c+5])
        c += 1
    
    mesh2 = trimesh.Trimesh(vertices=verts,faces=faces)

    # Sample foliage points. The point budget is allocated across 1 m height
    # slices in proportion to the number of branches rooted in each slice
    # (within a slice, faces stay weighted by their area), so a section of
    # crown carrying more branches receives proportionally more points and
    # dense multi-leader sections sample as densely as sparse ones. Plain
    # area-weighted sampling would instead give every slice the same points
    # per unit of mesh area, leaving branch-dense sections starved.
    points = None
    if len(built_z):
        face_areas = np.asarray(mesh2.area_faces)
        if face_areas.sum() > 0:
            slice_id = np.floor(np.asarray(built_z)).astype(int)
            n_slices = int(slice_id.max()) + 1
            count_s = np.bincount(slice_id, minlength=n_slices).astype(float)
            area_b = face_areas.reshape(-1, 4).sum(axis=1)   # mesh area per branch
            area_s = np.bincount(slice_id, weights=area_b, minlength=n_slices)
            boost = np.zeros(len(built_z))
            ok = area_s[slice_id] > 0
            boost[ok] = count_s[slice_id[ok]] / area_s[slice_id[ok]]
            fw = face_areas * np.repeat(boost, 4)            # per-face weights
            if fw.sum() > 0:
                points = mesh2.sample(Nfol, face_weight=fw)
    if points is None:
        points = mesh2.sample(Nfol)
    points[:,2] += model_params['foliage_noise']*np.random.random(points.shape[0])
    points_class = np.concatenate((points,class_id['folliage']*np.ones((Nfol,1))),axis=1)
    
    points_all = np.concatenate((points_all,points_class))
    
    if not return_skeleton:
        return points_all
    
    # --- Skeleton: trunk polylines + arc lengths (from stem_data) ---
    trunk_polylines = []
    trunk_segments = []
    trunk_len = 0.0
    for _stem in stem_data:
        (_xs, _ys, _zs, _p, _kx, _ky, _zr) = _stem
        _poly = np.column_stack([_xs, _ys, _zs]).astype(float)
        _seglen = float(np.sum(np.linalg.norm(np.diff(_poly, axis=0), axis=1)))
        trunk_polylines.append(_poly.tolist())
        trunk_segments.append({'z_range': [float(_zr[0]), float(_zr[1])],
                               'length_m': _seglen})
        trunk_len += _seglen
    
    skeleton = {
        'base_xy': [float(tx), float(ty)],
        'height_m': float(h),
        'base_radius_m': float(r),
        'fork': int(fork),
        'fork_height_m': (float(split_log[0]['height_m']) if split_log else None),
        'n_splits': int(n_splits),
        'splits': split_log,
        'min_canopy_height_m': float(h1),
        'max_canopy_width_height_m': float(h2),
        'trunk': {
            'total_length_m': float(trunk_len),
            'n_segments': len(trunk_segments),
            'segments': trunk_segments,
            'polylines': trunk_polylines,
        },
        'branches': branch_skeleton,
        'n_branches': len(branch_skeleton),
    }
    return points_all, skeleton

# gen_simtree: Generates sample point cloud of synthetically-generated tree
def gen_simtree(Np=1024, stem_frac_lambda=0.3, model_params=None, tx=None, ty=None,
                return_skeleton=False):
    
    # create example default parameters if none specified
    if model_params is None:
        model_params = {}
        model_params['height_range'] = [30, 50]
        model_params['diam_range'] = [0.5, 1]
        model_params['split_height_range'] = [0.15, 0.5]
        model_params['split_prob'] = 0.5
        model_params['num_branches'] = [60, 100]
        model_params['min_can_height'] = [0.2, 0.5]
        model_params['max_can_width'] = 7.0
        model_params['max_can_width_height'] = [0.4, 0.8]
        model_params['tree_top_dist'] = 2.5
        model_params['tree_mid_dist'] = 0.5
        model_params['foliage_noise'] = 0.5
        # cascading trunk splits for forked trees (1 = the original single fork)
        model_params['num_splits'] = [1, 1]
    
    # stem_frac_lambda: fraction of stem hits relative to foliage
    if isinstance(stem_frac_lambda,list): # stem_frac within a specified range
        stem_frac_lambda = stem_frac_lambda[0] + np.random.random()*(stem_frac_lambda[1]-stem_frac_lambda[0])
    
    # Work out numbers of points
    Nground = 0
    Nstem = int(stem_frac_lambda*Np)
    Nfol = Np-(Nground+Nstem)
    
    # Generate initial point cloud
    points_all = np.zeros((0,4))
    
    # Generate stem and foliage
    result = gen_single_simtree(Nstem=Nstem,Nfol=Nfol,model_params=model_params,
                                tx=tx,ty=ty,return_skeleton=return_skeleton)
    if return_skeleton:
        (points_tree, skeleton) = result
        skeleton['n_points'] = int(Np)
        skeleton['n_stem_points'] = int(Nstem)
        skeleton['n_foliage_points'] = int(Nfol)
        skeleton['stem_frac_lambda'] = float(stem_frac_lambda)
    else:
        points_tree = result
    
    points_all = np.concatenate((points_all,points_tree))
    
    if return_skeleton:
        return points_all, skeleton
    return points_all


