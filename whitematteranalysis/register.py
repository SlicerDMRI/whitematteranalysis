import numpy

import vtk

import whitematteranalysis.fibers


class RegistrationInformation:
    def __init__(self):
        self._original_fibers = whitematteranalysis.fibers.FiberArray()
        self._moving_fibers = whitematteranalysis.fibers.FiberArray()
        # indices of moving fibers to compute the objective function
        self._moving_fiber_sample = []
        self.points_per_fiber = 5
        self.fiber_sample_size = 200

        # transformation matrices for internal use
        # (vtkTransform is returned by compute) 
        # rot x,y,z trans x,y,z scale x,y,z
        #self.transform = numpy.array([0, 0, 0, 0, 0, 0, 1, 1, 1]).astype(float)
        self.transform = numpy.array([0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0]).astype(float)
        self.modified = False
        self.x = []

    def initialize(self, polydata):
        # internal representation for fast similarity computation
        self._original_fibers.convert_from_polydata(polydata,
                                                self.points_per_fiber)
        #self._moving_fibers.convert_from_polydata(movingpd,
        #                                         self.points_per_fiber)
        self.modified = True
        self.apply_transform()

    def initialize_fiber_sample(self):
        # indices of moving fibers to compute the objective function
        self._moving_fiber_sample = numpy.random.random_integers(
            0, self._original_fibers.number_of_fibers - 1,
            self.fiber_sample_size)
        self.modified = True

    def apply_transform(self):
        # apply transform to moving fiber data IF the transform is modified
        if self.modified:
            self._moving_fibers = \
                self.transform_fiber_array(self._original_fibers,
                                           self.transform)

    def transform_fiber_array_NOT_USED(self, in_array, transform):
        """Transform in_array (of class FiberArray) by transform (9
        components, rotation about R,A,S, translation in R, A, S, and
        scale along R, A, S. Fibers are assumed to be in RAS.  Transformed
        fibers are returned."""
        
        out_array = whitematteranalysis.fibers.FiberArray()
        pd_in = in_array.convert_to_polydata()
        transformer = vtk.vtkTransformPolyDataFilter()
        if (vtk.vtkVersion().GetVTKMajorVersion() >= 6.0):
            transformer.SetInputData(pd_in)
        else:
            transformer.SetInput(pd_in) 
        vtktrans = self.convert_transform_to_vtk(transform)
        transformer.SetTransform(vtktrans)
        transformer.Update()
        pd_out = transformer.GetOutput()
        out_array.convert_from_polydata(pd_out, self.points_per_fiber)
        return out_array

    def transform_fiber_array(self, in_array, transform):
        """Transform in_array (of class FiberArray) by transform (9
        components, rotation about R,A,S, translation in R, A, S, and
        scale along R, A, S. Fibers are assumed to be in RAS.
        Transformed fibers are returned. """

        out_array = whitematteranalysis.fibers.FiberArray()
        out_array.number_of_fibers = in_array.number_of_fibers
        out_array.points_per_fiber = in_array.points_per_fiber
        # allocate array number of lines by line length
        out_array.fiber_array_r = numpy.zeros((in_array.number_of_fibers,
                                               in_array.points_per_fiber))
        out_array.fiber_array_a = numpy.zeros((in_array.number_of_fibers,
                                               in_array.points_per_fiber))
        out_array.fiber_array_s = numpy.zeros((in_array.number_of_fibers,
                                               in_array.points_per_fiber))

        vtktrans = self.convert_transform_to_vtk(transform)

        # for testing only
        #out_array_2 = self.transform_fiber_array_NOT_USED(in_array, transform)

        # Transform input array in place
        for lidx in range(0, in_array.number_of_fibers):
            for pidx in range(0, in_array.points_per_fiber):
                pt = vtktrans.TransformPoint(in_array.fiber_array_r[lidx, pidx],
                                             in_array.fiber_array_a[lidx, pidx], 
                                             in_array.fiber_array_s[lidx, pidx])
                out_array.fiber_array_r[lidx, pidx] = pt[0]
                out_array.fiber_array_a[lidx, pidx] = pt[1]
                out_array.fiber_array_s[lidx, pidx] = pt[2]

        # test. this confirmed results were equivalent to old method
        # with time consuming polydata conversion.
        #print "=========================**************====================="
        #print numpy.max(out_array.fiber_array_r - out_array_2.fiber_array_r)
        #print numpy.max(out_array.fiber_array_a - out_array_2.fiber_array_a)
        #print numpy.max(out_array.fiber_array_s - out_array_2.fiber_array_s)
        #print "=========================**************====================="

        # test. this confirmed in-place array modification results
        # were equivalent to old method with time consuming polydata
        # conversion.
        #print "=========================**************====================="
        #print numpy.max(in_array.fiber_array_r - out_array_2.fiber_array_r)
        #print numpy.max(in_array.fiber_array_a - out_array_2.fiber_array_a)
        #print numpy.max(in_array.fiber_array_s - out_array_2.fiber_array_s)
        #print "=========================**************====================="
        
        return out_array


    def set_transform(self, input_transform):
        input_transform = numpy.array(input_transform)
        # decide whether transform was modified
        if numpy.max(numpy.abs(self.transform - input_transform)) == 0.0:
            self.modified = False
        else:
            # directly set it. assume it has 9 components 
            self.transform = numpy.copy(input_transform)
            self.modified = True
        #print "Set_transform, modified:", self.modified, "T0:", self.transform, "T1:", input_transform

    def convert_transform_to_vtk(self, transform=None):
        """ Produce an output vtkTransform corresponding to the
        registration results. Optionally can input a 9-component
        transform vector."""
        
        if transform is None:
            transform = self.transform
        
        vtktrans = vtk.vtkTransform()

        vtktrans.RotateX(transform[0] * (180 / numpy.pi))
        vtktrans.RotateY(transform[1] * (180 / numpy.pi))
        vtktrans.RotateZ(transform[2] * (180 / numpy.pi))

        vtktrans.Translate(transform[3],
                           transform[4], transform[5])

        vtktrans.Scale(transform[6],
                       transform[7], transform[8])

        #// Update affine transformation: Add shearing
        #vtkMatrix4x4 *skewx= vtkMatrix4x4::New();  skewx->Identity();
        #skewx->SetElement(2, 1,tan(_szy*(pi/180.0)));   skewx->SetElement(1, 2,tan(_syz*(pi/180.0))); 
        #vtkMatrix4x4 *skewy= vtkMatrix4x4::New();   skewy->Identity();
        #skewy->SetElement(2, 0,tan(_szx*(pi/180.0)));   skewy->SetElement(0, 2,tan(_sxz*(pi/180.0)));
        #vtkMatrix4x4 *skewz= vtkMatrix4x4::New();   skewz->Identity();
        #skewz->SetElement(1, 0,tan(_sxy*(pi/180.0)));    skewz->SetElement(0, 1,tan(_syx*(pi/180.0))); 
        #tr->Concatenate(skewx);   tr->Concatenate(skewy);   tr->Concatenate(skewz);
        sxy = transform[9]
        sxz = transform[10]
        syx = transform[11]
        syz = transform[12]
        szx = transform[13]
        szy = transform[14]
        skewx = vtk.vtkMatrix4x4()
        skewy = vtk.vtkMatrix4x4()
        skewz = vtk.vtkMatrix4x4()
        skewx.SetElement(2, 1, numpy.tan(szy))
        skewx.SetElement(1, 2, numpy.tan(syz))
        skewy.SetElement(2, 0, numpy.tan(szx))
        skewy.SetElement(0, 2, numpy.tan(sxz))
        skewz.SetElement(1, 0, numpy.tan(sxy))
        skewz.SetElement(0, 1, numpy.tan(syx))
        vtktrans.Concatenate(skewx)
        vtktrans.Concatenate(skewy)
        vtktrans.Concatenate(skewz)
        #del skewx
        #del skewy
        #del skewz
        return vtktrans
