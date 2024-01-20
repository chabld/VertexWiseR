"""Python function to plot folded and unfolded hippocampus
   Modified from jordandekraker/hippunfold_toolbox https://github.com/jordandekraker/hippunfold_toolbox/blob/main/hippunfold_toolbox/plotting.py
   -removed dentate gyrus label
   -removed third dimension
   -load template point and cell data from github
   -no longer need to import hippunfold_toolbox
"""
import numpy as np
import pickle 
import copy
import urllib.request

from brainspace.mesh.mesh_io import read_surface
from brainspace.plotting import  plot_surf
from brainspace.mesh import mesh_creation as mc

def surfplot_canonical_foldunfold(cdata, hemis=['L','R'],size=[350,400],**qwargs):
    '''
    Plots canonical folded and unfolded surfaces (hipp/dentate; folded/unfolded). This is good for cdata that isn't specific to one subject (eg. maybe it has been averaged across many subjects).
    
    cdata: array with the shape Vx2xF, where V is the number of vertices, 2 is the number of hemispheres (unless specified), and F is the number of rows/features
    kwargs: see https://brainspace.readthedocs.io/en/latest/generated/brainspace.plotting.surface_plotting.plot_surf.html#brainspace.plotting.surface_plotting.plot_surf
    '''
    hipdat = pickle.load(urllib.request.urlopen("https://raw.githubusercontent.com/CogBrainHealthLab/VertexWiseR/main/data/hip_points_cells.pkl"))

    rh = mc.build_polydata(hipdat[0], cells=hipdat[1])
    ru = mc.build_polydata(hipdat[2], cells=hipdat[1])
    ru.Points = ru.Points[:,[1,0,2]] # reorient unfolded
    ru.Points[:,1]=-ru.Points[:,1]
    # flip to get left hemisphere
    lh = mc.build_polydata(rh.Points.copy(), cells=rh.GetCells2D().copy())
    lh.Points[:,0] = -lh.Points[:,0]
    lu = mc.build_polydata(ru.Points.copy(), cells=ru.GetCells2D().copy())
    lu.Points[:,0] = -lu.Points[:,0]
    
    # do some cdata formatting
    cdata = np.reshape(cdata,[cdata.shape[0],len(hemis),-1])

    # set up layout
    surfDict = {'Lf':lh, 'Lu':lu, 'Rf':rh, 'Ru':ru}
    surfList = np.ones((cdata.shape[2],len(hemis)*2),dtype=object)
    arrName = np.ones((cdata.shape[2],len(hemis)*2),dtype=object)
    for h,hemi in enumerate(hemis):
        if hemi=='L':
            surfList[:,[h,h+1]] = np.array([f"{hemi}f",f"{hemi}u"])
            for f in range(cdata.shape[2]):
                lh.append_array(cdata[:,h,f], name=f'feature{f}', at='point')
                lu.append_array(cdata[:,h,f], name=f'feature{f}', at='point')
        elif hemi=='R':
            surfList[:,[h*2,h*2 +1]] = np.array([f"{hemi}u",f"{hemi}f"])
            for f in range(cdata.shape[2]):
                rh.append_array(cdata[:,h,f], name=f'feature{f}', at='point')
                ru.append_array(cdata[:,h,f], name=f'feature{f}', at='point')
        for f in range(cdata.shape[2]):
            arrName[f,:] = f'feature{f}'
    
    # extra parameters
    new_qwargs = dict(zoom=1.7, nan_color=(0,0,0,0))
    new_qwargs.update(qwargs)
    new_size=copy.deepcopy(size)
    new_size[0] = new_size[0]*len(hemis)
    new_size[1] = new_size[1]*cdata.shape[2]
    if 'color_bar' in qwargs:
        new_size[0] = new_size[0]+60
    p = plot_surf(surfDict,surfList, array_name=arrName, size=new_size,view=['ventral','dorsal','dorsal','ventral'],return_plotter=True, **new_qwargs)
    return p

