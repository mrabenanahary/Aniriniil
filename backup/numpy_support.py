"""This module adds support to easily import and export NumPy
(http://numpy.scipy.org) arrays into/out of VTK arrays.  The code is
loosely based on TVTK (https://svn.enthought.com/enthought/wiki/TVTK).

This code depends on an addition to the VTK data arrays made by Berk
Geveci to make it support Python's buffer protocol (on Feb. 15, 2008). 

The main functionality of this module is provided by the two functions:
    numpy_to_vtk, 
    vtk_to_numpy.


Caveats:
--------

 - Bit arrays in general do not have a numpy equivalent and are not
   supported.  Char arrays are also not easy to handle and might not
   work as you expect.  Patches welcome.

 - You need to make sure you hold a reference to a Numpy array you want
   to import into VTK.  If not you'll get a segfault (in the best case).
   The same holds in reverse when you convert a VTK array to a numpy
   array -- don't delete the VTK array.


Created by Prabhu Ramachandran in Feb. 2008.
"""

import vtk
import vtkConstants
import numpy

# Useful constants for VTK arrays.
VTK_ID_TYPE_SIZE = vtk.vtkIdTypeArray().GetDataTypeSize()
if VTK_ID_TYPE_SIZE == 4:
    ID_TYPE_CODE = numpy.int32
elif VTK_ID_TYPE_SIZE == 8:
    ID_TYPE_CODE = numpy.int64

VTK_LONG_TYPE_SIZE = vtk.vtkLongArray().GetDataTypeSize()
if VTK_LONG_TYPE_SIZE == 4:
    LONG_TYPE_CODE = numpy.int32
    ULONG_TYPE_CODE = numpy.uint32
elif VTK_LONG_TYPE_SIZE == 8:
    LONG_TYPE_CODE = numpy.int64
    ULONG_TYPE_CODE = numpy.uint64


def get_vtk_array_type(numpy_array_type):
    """Returns a VTK typecode given a numpy array."""
    # This is a Mapping from numpy array types to VTK array types.
    _np_vtk = {numpy.character:vtkConstants.VTK_UNSIGNED_CHAR,
                numpy.uint8:vtkConstants.VTK_UNSIGNED_CHAR,
                numpy.uint16:vtkConstants.VTK_UNSIGNED_SHORT,
                numpy.uint32:vtkConstants.VTK_UNSIGNED_INT,
                ULONG_TYPE_CODE:vtkConstants.VTK_UNSIGNED_LONG,
                numpy.int8:vtkConstants.VTK_CHAR,
                numpy.int16:vtkConstants.VTK_SHORT,
                numpy.int32:vtkConstants.VTK_INT,
                LONG_TYPE_CODE:vtkConstants.VTK_LONG,
                numpy.float32:vtkConstants.VTK_FLOAT,
                numpy.float64:vtkConstants.VTK_DOUBLE,
                numpy.complex64:vtkConstants.VTK_FLOAT,
                numpy.complex128:vtkConstants.VTK_DOUBLE}
    try:
        return _np_vtk[numpy_array_type]
    except KeyError:
        for key in _np_vtk:
            if numpy.issubdtype(numpy_array_type, key):
                return _np_vtk[key]

def get_vtk_to_numpy_typemap():
    """Returns the VTK array type to numpy array type mapping."""
    _vtk_np = {vtkConstants.VTK_BIT:numpy.bool,
                vtkConstants.VTK_CHAR:numpy.int8,
                vtkConstants.VTK_UNSIGNED_CHAR:numpy.uint8,
                vtkConstants.VTK_SHORT:numpy.int16,
                vtkConstants.VTK_UNSIGNED_SHORT:numpy.uint16,
                vtkConstants.VTK_INT:numpy.int32,
                vtkConstants.VTK_UNSIGNED_INT:numpy.uint32,
                vtkConstants.VTK_LONG:LONG_TYPE_CODE,
                vtkConstants.VTK_UNSIGNED_LONG:ULONG_TYPE_CODE,
                vtkConstants.VTK_ID_TYPE:ID_TYPE_CODE,
                vtkConstants.VTK_FLOAT:numpy.float32,
                vtkConstants.VTK_DOUBLE:numpy.float64}
    return _vtk_np


def get_numpy_array_type(vtk_array_type):
    """Returns a numpy array typecode given a VTK array type."""
    return get_vtk_to_numpy_typemap()[vtk_array_type]


def create_vtk_array(vtk_arr_type):
    """Internal function used to create a VTK data array from another
    VTK array given the VTK array type.
    """
    tmp = vtk.vtkDataArray.CreateDataArray(vtk_arr_type)
    # We need to manually dereference objects created with anything
    # but the constructor.
    tmp.UnRegister(None)
    return tmp


def numpy_to_vtk(num_array, deep=0):
    """Converts a contiguous real numpy Array to a VTK array object.

    This function only works for real arrays that are contiguous.
    Complex arrays are NOT handled.  It also works for multi-component
    arrays.  However, only 1, and 2 dimensional arrays are supported.
    This function is very efficient, so large arrays should not be a
    problem.

    If the second argument is set to 1, the array is deep-copied from
    from numpy. This is not as efficient as the default behavior
    (shallow copy) and uses more memory but detaches the two arrays
    such that the numpy array can be released.

    WARNING: You must maintain a reference to the passed numpy array, if
    the numpy data is gc'd and VTK will point to garbage which will in
    the best case give you a segfault.

    Parameters
    ----------

    - num_array :  a contiguous 1D or 2D, real numpy array.

    """

    z = numpy.asarray(num_array)

    shape = z.shape
    assert z.flags.contiguous, 'Only contiguous arrays are supported.'
    assert len(shape) < 3, \
           "Only arrays of dimensionality 2 or lower are allowed!"
    assert not numpy.issubdtype(z.dtype, complex), \
           "Complex numpy arrays cannot be converted to vtk arrays."\
           "Use real() or imag() to get a component of the array before"\
           " passing it to vtk."

    # First create an array of the right type by using the typecode.
    vtk_typecode = get_vtk_array_type(z.dtype)
    result_array = create_vtk_array(vtk_typecode)

    # Find the shape and set number of components.
    if len(shape) == 1:
        result_array.SetNumberOfComponents(1)
    else:
        result_array.SetNumberOfComponents(shape[1])

    result_array.SetNumberOfTuples(shape[0])

    # Ravel the array appropriately.
    arr_dtype = get_numpy_array_type(vtk_typecode)
    if numpy.issubdtype(z.dtype, arr_dtype):
        z_flat = numpy.ravel(z)
    else:
        z_flat = numpy.ravel(z).astype(arr_dtype)

    # Point the VTK array to the numpy data.  The last argument (1)
    # tells the array not to deallocate.
    result_array.SetVoidArray(z_flat, len(z_flat), 1)
    if deep:
        copy = result_array.NewInstance()
        # NewInstance sets the refcount to 3 and this causes a severe
        # memory leak.        
        copy.SetReferenceCount(2)
        copy.DeepCopy(result_array)
        result_array = copy
    return result_array

def vtk2array(vtk_array):
    return vtk_to_numpy(vtk_array)

def vtk_to_numpy(vtk_array):
    """Converts a VTK data array to a numpy array.

    Given a subclass of vtkDataArray, this function returns an
    appropriate numpy array containing the same data -- it actually
    points to the same data.  

    WARNING: This does not work for bit arrays.

    Parameters
    ----------

    - vtk_array : `vtkDataArray`

      The VTK data array to be converted.

    """
    print(vtk_array)
    
    typ = vtk_array.GetDataType()
    assert typ in get_vtk_to_numpy_typemap().keys(), \
           "Unsupported array type %s"%typ
    assert typ != vtkConstants.VTK_BIT, 'Bit arrays are not supported.'

    shape = vtk_array.GetNumberOfTuples(), \
            vtk_array.GetNumberOfComponents()

    # Get the data via the buffer interface
    dtype = get_numpy_array_type(typ)
    result = numpy.frombuffer(vtk_array, dtype=dtype)
    if shape[1] == 1:
        shape = (shape[0], )
    result.shape = shape
    return result
