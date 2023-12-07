# -*- coding: utf-8 -*-

""" Class RenderPolyData for rendering tractography, and saving views. """

import os

import numpy as np
import vtk

from . import filter


def render(input_polydata, number_of_fibers=None, opacity=1, depth_peeling=False, scalar_bar=False, axes=False, scalar_range=None, data_mode="Cell", tube=True, colormap='jet', data_name=None, verbose=True):
    """ Function for easy matlab-like use of the rendering
    functionality.


    Returns RenderPolyData object which can be used directly for more
    functionality."""

    if verbose:
        print(f"<{os.path.basename(__file__)}> Initiating rendering.")
        
    if number_of_fibers is not None:
        if verbose:
            print(f"<{os.path.basename(__file__)}> Downsampling vtkPolyData: {number_of_fibers}")
        # downsample if requested
        input_polydata = filter.downsample(input_polydata, number_of_fibers, preserve_point_data=True, preserve_cell_data=True, verbose=verbose)

    if data_name is not None:
        if verbose:
            print(f"<{os.path.basename(__file__)}> Visualizing data: {data_name}")
        if data_mode == "Cell":
            input_polydata.GetCellData().SetActiveScalars(data_name)
        if data_mode == "Point":
            input_polydata.GetPointData().SetActiveScalars(data_name)

    ren = RenderPolyData()
    
    ren.render_polydata(input_polydata, opacity=opacity, depth_peeling=depth_peeling, scalar_bar=scalar_bar, axes=axes, scalar_range=scalar_range, data_mode=data_mode, tube=tube, colormap=colormap, verbose=verbose)

    if verbose:
        print(f"<{os.path.basename(__file__)}> Render pipeline created.")
    return ren

def save_views(render_object, directory=".", subjectID=None):
    render_object.save_views(directory, subjectID)


def get_jet_lookup_table():
    # create a jet colormap like matlab to make everything look better
    jet_r = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0625, 0.1250, 0.1875, 0.2500, 0.3125, 0.3750, 0.4375, 0.5000, 0.5625, 0.6250, 0.6875, 0.7500, 0.8125, 0.8750, 0.9375, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 0.9375, 0.8750, 0.8125, 0.7500, 0.6875, 0.6250, 0.5625, 0.5000]

    jet_g = [0,     0,      0,      0,      0,      0,      0,      0, 0.0625, 0.1250, 0.1875, 0.2500, 0.3125, 0.3750, 0.4375, 0.5000, 0.5625, 0.6250, 0.6875, 0.7500, 0.8125, 0.8750, 0.9375, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 0.9375, 0.8750, 0.8125, 0.7500, 0.6875, 0.6250, 0.5625, 0.5000, 0.4375, 0.3750, 0.3125, 0.2500, 0.1875, 0.1250, 0.0625,      0,      0,      0,      0,      0,      0,      0,      0,      0]

    jet_b = [0.5625, 0.6250, 0.6875, 0.7500, 0.8125, 0.8750, 0.9375, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 0.9375, 0.8750, 0.8125, 0.7500, 0.6875, 0.6250, 0.5625, 0.5000, 0.4375, 0.3750, 0.3125, 0.2500, 0.1875, 0.1250, 0.0625,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0]

    lut = vtk.vtkLookupTable()
    lut.SetNumberOfColors(len(jet_r))

    for i in range(len(jet_r)):
        lut.SetTableValue(i, jet_r[i], jet_g[i], jet_b[i])

    return(lut)

def get_hot_lookup_table():
    # create a hot colormap for p values.
    hot_r = [      1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1]

    hot_g = [    1.0000, 0.9841, 0.9683, 0.9524, 0.9365, 0.9206, 0.9048, 0.8889, 0.8730, 0.8571, 0.8413, 0.8254, 0.8095, 0.7937, 0.7778, 0.7619, 0.7460, 0.7302, 0.7143, 0.6984, 0.6825, 0.6667, 0.6508, 0.6349, 0.6190, 0.6032, 0.5873, 0.5714, 0.5556, 0.5397, 0.5238, 0.5079, 0.4921, 0.4762, 0.4603, 0.4444, 0.4286, 0.4127, 0.3968, 0.3810, 0.3651, 0.3492, 0.3333, 0.3175, 0.3016, 0.2857, 0.2698, 0.2540, 0.2381, 0.2222, 0.2063, 0.1905, 0.1746, 0.1587, 0.1429, 0.1270, 0.1111, 0.0952, 0.0794, 0.0635, 0.0476, 0.0317, 0.0159,      0]

    hot_b = [     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0]

    lut = vtk.vtkLookupTable()
    lut.SetNumberOfColors(len(hot_r))

    for i in range(len(hot_r)):
        lut.SetTableValue(i, hot_r[i], hot_g[i], hot_b[i])

    return(lut)

def argsort_by_jet_lookup_table(rgb_color):
    # create a jet colormap like matlab
    jet_r = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0625, 0.1250, 0.1875, 0.2500, 0.3125, 0.3750, 0.4375, 0.5000, 0.5625, 0.6250, 0.6875, 0.7500, 0.8125, 0.8750, 0.9375, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 0.9375, 0.8750, 0.8125, 0.7500, 0.6875, 0.6250, 0.5625, 0.5000]

    jet_g = [0,     0,      0,      0,      0,      0,      0,      0, 0.0625, 0.1250, 0.1875, 0.2500, 0.3125, 0.3750, 0.4375, 0.5000, 0.5625, 0.6250, 0.6875, 0.7500, 0.8125, 0.8750, 0.9375, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 0.9375, 0.8750, 0.8125, 0.7500, 0.6875, 0.6250, 0.5625, 0.5000, 0.4375, 0.3750, 0.3125, 0.2500, 0.1875, 0.1250, 0.0625,      0,      0,      0,      0,      0,      0,      0,      0,      0]

    jet_b = [0.5625, 0.6250, 0.6875, 0.7500, 0.8125, 0.8750, 0.9375, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 0.9375, 0.8750, 0.8125, 0.7500, 0.6875, 0.6250, 0.5625, 0.5000, 0.4375, 0.3750, 0.3125, 0.2500, 0.1875, 0.1250, 0.0625,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0]

    # add missing purples, cyans, yellows to this jet map.    
    jet_r = jet_r + [0.8, 0.8, 0.6, 0.2, 0.8, 0.8]
    jet_g = jet_g + [0.2, 0.2, 0.2, 0.8, 0.8, 0.8]
    jet_b = jet_b + [0.6, 0.8, 0.6, 0.8, 0.4, 0.2]

    # map from 0..255 to 0..1
    rgb_color = np.divide(rgb_color, 255)
    
    # sort input rgb into this colormap order
    match_jet = list()
    match_dist = list()
    for col_idx in range(len(rgb_color)):
        diff_r = rgb_color[col_idx,0]-jet_r
        diff_g = rgb_color[col_idx,1]-jet_g
        diff_b = rgb_color[col_idx,2]-jet_b

        mag = np.sqrt(np.multiply(diff_r,diff_r) + np.multiply(diff_g,diff_g) + np.multiply(diff_b,diff_b) )
        match_jet.append(np.argmin(mag))
        match_dist.append(np.min(mag))
        # uncomment for testing of worst color matches
        #if np.min(mag) > 0.3:
        #    print np.min(mag), "Color matched:", rgb_color[col_idx,:], "Idx:", match_jet[col_idx], ":", jet_r[match_jet[col_idx]], jet_g[match_jet[col_idx]], jet_b[match_jet[col_idx]], "\n"
    
    # Return indices that will sort these colors (centroids) in their match order
    return(np.argsort(np.array(match_jet)))
        
class RenderPolyData:

    """Makes a vtk render window to display polydata tracts"""

    def build_vtk_renderer(self):

        # offscreen rendering
        if (vtk.vtkVersion().GetVTKMajorVersion() == 5.0):
            graphics_factory = vtk.vtkGraphicsFactory()
            graphics_factory.SetOffScreenOnlyMode( 1);
            graphics_factory.SetUseMesaClasses( 1 );
            imaging_factory = vtk.vtkImagingFactory()
            imaging_factory.SetUseMesaClasses( 1 );
  
        self.renderer = vtk.vtkRenderer()
        self.renderer.SetBackground(1, 1, 1)
        self.render_window = vtk.vtkRenderWindow()
        # offscreen rendering
        self.render_window.SetOffScreenRendering(1)

        self.render_window.AddRenderer(self.renderer)
        # create a renderwindowinteractor
        #if self.interact:
        #self.iren = vtk.vtkRenderWindowInteractor()
        #self.iren.SetRenderWindow(self.render_window)

        # scalar bar        
        self.scalarbar = vtk.vtkScalarBarActor()
        # To avoid uninitialized warning in VTK 6
        self.scalarbar.SetTitle("") 
        # black text since background is white for printing
        self.scalarbar.GetLabelTextProperty().SetColor(0, 0, 0)
        self.scalarbar.GetTitleTextProperty().SetColor(0, 0, 0)
        
    def __init__(self):
        # parameters
        #self.range = [-1, 1]
        # performance options
        self.verbose = 0
        self.magnification = 4
        #self.interact = True
        self.number_of_tube_sides = 9
        self.verbose = False
        # permanent objects that persist with this one
        self.build_vtk_renderer()
        
    def __del__(self):
        try:
            if self.verbose:
                print(f"<{os.path.basename(__file__)}> in DELETE")
            #del self.renderer
            #del self.render_window
            #del self.iren
            #del self.scalarbar
        except Exception:
            print(f"<{os.path.basename(__file__)}> ERROR: deletion failed")

    def scalar_bar_on(self):
        self.renderer.AddActor2D(self.scalarbar)

    def scalar_bar_off(self):
        self.renderer.RemoveActor2D(self.scalarbar)

    def render_polydata(self, input_polydata, scalar_range=None,
                        tube=True, opacity=1, depth_peeling=False,
                        scalar_bar=False, scalar_bar_title=None,
                        axes=False, data_mode="Cell", colormap='jet', verbose=True):

        # re-do this every time because the user may close the window
        # and re-use the object
        #self.build_vtk()

        if verbose:
            print(f"<{os.path.basename(__file__)}> Rendering vtkPolyData.")
        
        # actor and mapper
        mapper = vtk.vtkPolyDataMapper()
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        # assign actor to the renderer
        self.renderer.AddActor(actor)


        # lookup table for everything
        #self.lut = vtk.vtkLookupTable()
        if colormap == 'hot':
            lut = get_hot_lookup_table()
        else:
            lut = get_jet_lookup_table()

        mapper.UseLookupTableScalarRangeOn()
        mapper.SetLookupTable(lut)

        # see if we have RGB data, if so the scalar bar is meaningless
        self.render_RGB = False
        if input_polydata.GetCellData().GetScalars():
            if input_polydata.GetCellData().GetScalars().GetNumberOfComponents() == 3:
                self.render_RGB = True
                self.renderer.RemoveActor2D(self.scalarbar)
        if verbose:
            print(f"<{os.path.basename(__file__)}> RGB: f{self.render_RGB}")

        if data_mode == "Cell":
            mapper.SetScalarModeToUseCellData()
        elif data_mode == "Point":
            mapper.SetScalarModeToUsePointData()

        #if input_polydata.GetNumberOfLines() > 0:
        #    mapper.SetScalarModeToUseCellData()
        
        if tube & (input_polydata.GetNumberOfLines() > 0):
            # make a tube filter
            tube_filter = vtk.vtkTubeFilter()
            tube_filter.SetNumberOfSides(self.number_of_tube_sides)
            if (vtk.vtkVersion().GetVTKMajorVersion() > 5.0):
                tube_filter.SetInputData(input_polydata)
            else:
                tube_filter.SetInput(input_polydata)
            if (vtk.vtkVersion().GetVTKMajorVersion() >= 6.0):
                mapper.SetInputConnection(tube_filter.GetOutputPort())
            else:
                mapper.SetInput(tube_filter.GetOutput())
        else:
            if (vtk.vtkVersion().GetVTKMajorVersion() >= 6.0):
                mapper.SetInputData(input_polydata)
            else:
                mapper.SetInput(input_polydata)            


        if scalar_range is None:
            # check polydata scalar range. 
            range_low, range_high = input_polydata.GetScalarRange()
            # if these numbers are between -1 and 1 assume laterality
            # data, so make sure colorbar is centered at 0 for clearer
            # interpretation
            if (range_low >= -1.0 ) and (range_high <= 1.0):
                range_final = max(abs(range_low), abs(range_high))
                lut.SetRange([-range_final, range_final])
            else:
                # otherwise just use the range directly
                lut.SetRange([range_low, range_high])

        else:
            lut.SetRange(scalar_range)

        if opacity != 1:
            actor.GetProperty().SetOpacity(opacity)

        if depth_peeling is not False:
            # 1. Use a render window with alpha bits (as initial value is 0 (false)):
            self.render_window.SetAlphaBitPlanes(1)
            # 2. Force to not pick a framebuffer with a multisample buffer
            # (as initial value is 8):
            self.render_window.SetMultiSamples(0)
            # 3. Choose to use depth peeling (if supported) (initial value is 0 (false)):
            self.renderer.SetUseDepthPeeling(1)
            # 4. Set depth peeling parameters
            # - Set the maximum number of rendering passes (initial value is 4):
            #renderer->SetMaximumNumberOfPeels(maxNoOfPeels);
            #// - Set the occlusion ratio (initial value is 0.0, exact image):
            #renderer->SetOcclusionRatio(occlusionRatio);
            
        if scalar_bar:
            self.scalarbar.SetLookupTable(lut)
            self.renderer.AddActor2D(self.scalarbar)
            if scalar_bar_title is not None:
                self.scalarbar.SetTitle(scalar_bar_title)
        else:
            self.renderer.RemoveActor2D(self.scalarbar)
        
        # Create the axes and the associated mapper and actor.
        self.axes = vtk.vtkAxes()
        self.axes.SetOrigin(0, 0, 0)
        self.axes.SetScaleFactor(100)
        self.axes.SymmetricOn()
        self.axes_mapper = vtk.vtkPolyDataMapper()
        if (vtk.vtkVersion().GetVTKMajorVersion() >= 6.0):
            self.axes_mapper.SetInputConnection(self.axes.GetOutputPort())
        else:
            self.axes_mapper.SetInput(self.axes.GetOutput())
        self.axes_actor = vtk.vtkActor()
        self.axes_actor.SetMapper(self.axes_mapper)
        self.axes_mapper.SetLookupTable(lut)
        if axes:
            self.renderer.AddActor(self.axes_actor)
        else:
            self.renderer.RemoveActor(self.axes_actor)
        
        try:
            self.render_window.Render()
        except Exception:
            print(f"<{os.path.basename(__file__)}> ERROR. Failed to render. Check display settings.")
            # do not re-raise the exception because we may be in the
            # middle of running this pipeline via ssh, and this is not
            # a fatal error for the computation

        # enable user interface interactor
        #self.iren.Start()

    def view_superior(self):
        cam = self.renderer.GetActiveCamera()
        cam.SetPosition([0, 0, 400])
        cam.SetViewUp(0, 1, 0)
        self.render_window.Render()
        self.renderer.ResetCameraClippingRange()

    def view_inferior(self):
        cam = self.renderer.GetActiveCamera()
        cam.SetPosition([0, 0, -400])
        cam.SetViewUp(0, 1, 0)
        self.render_window.Render()
        self.renderer.ResetCameraClippingRange()

    def view_left(self):
        cam = self.renderer.GetActiveCamera()
        cam.SetPosition([-400, 0, 0])
        #cam.SetRoll(90)
        cam.SetViewUp(0, 0, 1)
        self.render_window.Render()
        self.renderer.ResetCameraClippingRange()

    def view_right(self):
        cam = self.renderer.GetActiveCamera()
        cam.SetPosition([400, 0, 0])
        cam.SetViewUp(0, 0, 1)
        self.render_window.Render()
        self.renderer.ResetCameraClippingRange()

    def view_anterior(self):
        cam = self.renderer.GetActiveCamera()
        cam.SetPosition([0, 400, 0])
        cam.SetViewUp(0, 0, 1)
        #cam.SetRoll(90)
        self.render_window.Render()
        self.renderer.ResetCameraClippingRange()

    def view_posterior(self):
        cam = self.renderer.GetActiveCamera()
        cam.SetPosition([0, -400, 0])
        cam.SetViewUp(0, 0, 1)
        #cam.SetRoll(90)
        self.render_window.Render()
        self.renderer.ResetCameraClippingRange()

    def save_image(self, filename="test.jpg"):
        img = vtk.vtkWindowToImageFilter()
        img.SetInput(self.render_window)
        if (vtk.vtkVersion().GetVTKMajorVersion() >= 8.0):
            img.SetScale(self.magnification, self.magnification)
        else:
            img.SetMagnification(self.magnification)
        img.Update()
        #writer = vtk.vtkPNGWriter()
        writer = vtk.vtkJPEGWriter()
        if (vtk.vtkVersion().GetVTKMajorVersion() >= 6.0):
            writer.SetInputConnection(img.GetOutputPort())
        else:
            writer.SetInput(img.GetOutput())
        
        writer.SetFileName(filename)
        writer.Write()
        del writer
        del img

    def save_views(self, directory=".", subjectID=None,  verbose=True):

        if verbose:
            print(f"<{os.path.basename(__file__)}> Saving rendered views to disk: {directory}")
        
        if not os.path.isdir(directory):
            print(f"<{os.path.basename(__file__)}> ERROR: directory does not exist. {directory}")
            return

        #ext = ".png"
        ext = ".jpg"
        
        # Use subject ID as part of filename for easier visual 
        # identification of problem cases
        if subjectID is not None:
            fname_sup = f"view_sup_{subjectID}{ext}"
            fname_inf = f"view_inf_{subjectID}{ext}"
            fname_left = f"view_left_{subjectID}{ext}"
            fname_right = f"view_right_{subjectID}{ext}"
            fname_ant = f"view_ant_{subjectID}{ext}"
            fname_post = f"view_post_{subjectID}{ext}"
        else:
            fname_sup = f"view_sup{ext}"
            fname_inf = f"view_inf{ext}"
            fname_left = f"view_left{ext}"
            fname_right = f"view_right{ext}"
            fname_ant = f"view_ant{ext}"
            fname_post = f"view_post{ext}"
            
        # sometimes the model gets clipped, this mostly fixes it
        self.renderer.ResetCameraClippingRange()

        # init
        #self.renderer.RemoveActor2D(self.scalarbar)

        # superior
        self.view_superior()
        #if self.render_RGB == False:
        #    self.renderer.AddActor2D(self.scalarbar)
        #    self.save_image(os.path.join(directory, "view_sup_scalar_bar.png"))
        #self.renderer.RemoveActor2D(self.scalarbar)
        self.save_image(os.path.join(directory, fname_sup))

        # inferior
        self.view_inferior()
        #if self.render_RGB == False:
        #    self.renderer.AddActor2D(self.scalarbar)
        #    self.save_image(os.path.join(directory, "view_inf_scalar_bar.png"))
        #self.renderer.RemoveActor2D(self.scalarbar)
        self.save_image(os.path.join(directory, fname_inf))

        # left
        self.view_left()
        self.save_image(os.path.join(directory,  fname_left))

        # right
        self.view_right()
        self.save_image(os.path.join(directory,  fname_right))

        # anterior
        self.view_anterior()
        self.save_image(os.path.join(directory,  fname_ant))

        # posterior
        self.view_posterior()
        self.save_image(os.path.join(directory,  fname_post))

        self.render_window.Render()

def histeq(values,nbr_bins=256):

   #get image histogram
   imhist,bins = np.histogram(values,nbr_bins,normed=True)
   cdf = imhist.cumsum() #cumulative distribution function
   cdf = 255 * cdf / cdf[-1] #normalize

   #use linear interpolation of cdf to find new pixel values
   new_values = np.interp(values,bins[:-1],cdf)

   return new_values, cdf
