#=============================================================================
import numpy as np
import vtk as v
from numpy import loadtxt as pylabload
import numpy_support as ah
import sys, time
import struct
from vtk.util.numpy_support import vtk_to_numpy
import vtktonumpy as ah_vtk
import numpy_support as ah
import read

def my_get_pointdata(mag_config_name,offset,filenameout='data',type='pvtu',attribute_mode='cell',mirror=None):
    

    filename=''.join([filenameout,repr(offset).zfill(4),'.vtu'])
    reader = v.vtkXMLUnstructuredGridReader()
    reader.SetFileName(filename)
    reader.Update()

    print('=== Reading ',filename,' ===')

    nodes_vtk_array = reader.GetOutput().GetPoints().GetData()
    
    nodes_nummpy_array = vtk_to_numpy(nodes_vtk_array)
    
    return {'pointdata':reader.GetOutput(), 'points': nodes_nummpy_array}


def rot2D(mag_config_name,offset,nphi,filenameout='data',type='pvtu'):
    """
    creates a mock vtu datastructure with rotated values.
    No connectivity is set, vtu is just used as a container for pointdata     since
    the following routines expect to get a vtu unstructured datatype.
    """

    data=my_get_pointdata(mag_config_name,offset,filenameout=filenameout,type=type)

    points2D=data.get('points')
    
    data3D = v.vtkUnstructuredGrid()

    nelem2D = len(points2D)
    nelem3D = nelem2D*nphi

    print('preparing',nelem3D, 'points')

    points3D = np.empty((nelem3D,3))

    phis=np.linspace(0.,2.*np.pi*float(nphi)/float(nphi+1),nphi)

    i3D = 0
    for i2D in range(nelem2D):
# pick next 2D point:
        point = points2D[i2D]
        r = point[0]
        z = point[1]
        for phi in phis:
            sinphi = np.sin(phi)
            cosphi = np.cos(phi)
            points3D[i3D,0] = r*cosphi
            points3D[i3D,1] = r*sinphi
            points3D[i3D,2] = z
            i3D=i3D+1
    
    

    vtkpoints=ah_vtk.array2vtkPoints(points3D)
    #vtkpoints = numpy_to_vtk(points3D)
    data3D.SetPoints(vtkpoints)

# rotate variables:
    ivar = 0
    while ivar < data.get('pointdata').GetPointData().GetNumberOfArrays():
        var = data.get('pointdata').GetPointData().GetArrayName(ivar)
        #print 'rotating variable:',ivar,var
        i3D=0
        #print 'is ivaaaar',ivar,var
# Treat vectors:
        if var == 'b1' or var == 'v1':
            array2D_vec = np.empty((nelem2D,3))
            array3D_vec = np.empty((nelem3D,3))
            bx = np.empty(nelem3D)
            by = np.empty(nelem3D)
            bz = np.empty(nelem3D)
            array2D_vec[:,0] = read.extract(data.get('pointdata'),data.get('pointdata').GetPointData().GetArrayName(ivar),attribute_mode='point')
            array2D_vec[:,1] = read.extract(data.get('pointdata'),data.get('pointdata').GetPointData().GetArrayName(ivar+1),attribute_mode='point')
            array2D_vec[:,2] = read.extract(data.get('pointdata'),data.get('pointdata').GetPointData().GetArrayName(ivar+2),attribute_mode='point')
            for i2D in range(nelem2D):
                for phi in phis:
# Should have three components (makes no sense otherwise):
# Input in cylindrical coordinates, convention:
# b1==br; b2==bz; b3==bphi  and similar for velocities.
                    array3D_vec[i3D,0] = array2D_vec[i2D,0]
                    array3D_vec[i3D,1] = array2D_vec[i2D,1]
                    array3D_vec[i3D,2] = array2D_vec[i2D,2]
# Now create the rotated vectors:
                    sinphi = np.sin(phi)
                    cosphi = np.cos(phi)
                    bz[i3D] = array3D_vec[i3D,1]  
                    bx[i3D] = array3D_vec[i3D,0] * cosphi - array3D_vec[i3D,2] * sinphi
                    by[i3D] = array3D_vec[i3D,2] * cosphi + array3D_vec[i3D,0] * sinphi
                    i3D=i3D+1
            
            if var == 'v1':
                v1 = bz
                v2 = bx
                v3 = by
                
            if var == 'b1':
                b1 = bz
                b2 = bx
                b3 = by

            vtkarray3D = ah_vtk.array2vtk(bx)
            #vtkarray3D = numpy_to_vtk(bx)
            vtkarray3D.SetName(data.get('pointdata').GetPointData().GetArrayName(ivar))
            data3D.GetPointData().AddArray(vtkarray3D)

            vtkarray3D = ah_vtk.array2vtk(by)
            #vtkarray3D = numpy_to_vtk(by)
            vtkarray3D.SetName(data.get('pointdata').GetPointData().GetArrayName(ivar+1))
            data3D.GetPointData().AddArray(vtkarray3D)
 
            vtkarray3D = ah_vtk.array2vtk(bz)
            #vtkarray3D = numpy_to_vtk(bx)
            vtkarray3D.SetName(data.get('pointdata').GetPointData().GetArrayName(ivar+2))
            data3D.GetPointData().AddArray(vtkarray3D)
  
            ivar = ivar+3
        else:
# Treat scalars:
            array2D = np.empty(nelem2D)
            array3D = np.empty(nelem3D)
            array2D = read.extract(data.get('pointdata'),data.get('pointdata').GetPointData().GetArrayName(ivar),attribute_mode='point')

            for i2D in range(nelem2D):
                for phi in phis:
                    array3D[i3D] = array2D[i2D]
                    i3D=i3D+1
                    
            #vtkarray3D = numpy_to_vtk(array3D)
            #vtkarray3D.SetName(var)
            #data3D.GetPointData().AddArray(vtkarray3D)
            
            if ivar ==  0:
                rho = array3D

            if ivar == 1:
                flrho1 = array3D
            
            if ivar == 1:
                flrho2 = array3D
                
            if ivar == 3:
                flrho3 = array3D
            
            if ivar == 7:
                p = array3D
                
            if ivar == 12:
                lfac = array3D
            
            if mag_config_name == 'hydro':
                if ivar == 24:           
                    bturb = array3D
            

            ivar = ivar+1
    
    # bug with polo :: replace by rho in return !!! 

    if mag_config_name == 'hydro':

        return {'pointdata': data3D, 'points': points3D, 'rho': rho, 'flrho1': flrho1, 'flrho2': flrho2, 'flrho3': flrho3, 'v1': v1, 'v2': v2, 'v3': v3, 'b1': b1,
                'b2':b2, 'b3': b3,'bturb': bturb, 'p': p, 'lfac': lfac, 'fileid':str(offset)}
    
    if mag_config_name == 'toro' or mag_config_name == 'polo' or mag_config_name == 'helico':

        return {'pointdata': data3D, 'points': points3D, 'rho': rho, 'flrho1': flrho1, 'flrho2': flrho2, 'flrho3': flrho3, 'v1': v1, 'v2': v2, 'v3': v3, 'b1': b1,
                         'b2':b2, 'b3': b3, 'p': p, 'lfac': lfac, 'fileid':str(offset)}


