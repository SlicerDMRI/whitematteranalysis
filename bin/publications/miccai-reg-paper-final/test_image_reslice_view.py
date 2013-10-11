# This example shows how to load a 3D image into VTK and then reformat
# that image into a different orientation for viewing.  It uses
# vtkImageReslice for reformatting the image, and uses vtkImageActor
# and vtkInteractorStyleImage to display the image.  This InteractorStyle
# forces the camera to stay perpendicular to the XY plane.

import vtk

# these globals are needed for this to work. perhaps create a class
# later.  the below code was modified from a vtk example, to create a
# function to call for rendering.
reslice = vtk.vtkImageReslice()
window = vtk.vtkRenderWindow()
interactorStyle = vtk.vtkInteractorStyleImage()
interactor = vtk.vtkRenderWindowInteractor()
table = vtk.vtkLookupTable()
color = vtk.vtkImageMapToColors()
actor = vtk.vtkImageActor()
renderer = vtk.vtkRenderer()
actions = {}
    
def ButtonCallback(obj, event):
    #print "button"
    if event == "LeftButtonPressEvent":
        actions["Slicing"] = 1
    else:
        actions["Slicing"] = 0

def MouseMoveCallback(obj, event):
    #print "mouse move"
    (lastX, lastY) = interactor.GetLastEventPosition()
    (mouseX, mouseY) = interactor.GetEventPosition()
    if actions["Slicing"] == 1:
        #print "slicing"
        deltaY = mouseY - lastY
        reslice.GetOutput().UpdateInformation()
        sliceSpacing = reslice.GetOutput().GetSpacing()[2]
        matrix = reslice.GetResliceAxes()
        # move the center point that we are slicing through
        center = matrix.MultiplyPoint((0, 0, sliceSpacing*deltaY, 1))
        matrix.SetElement(0, 3, center[0])
        matrix.SetElement(1, 3, center[1])
        matrix.SetElement(2, 3, center[2])
        window.Render()
    else:
        interactorStyle.OnMouseMove()
  
def view_image(image):
    # Calculate the center of the volume
    image.UpdateInformation()
    (xMin, xMax, yMin, yMax, zMin, zMax) = image.GetWholeExtent()
    (xSpacing, ySpacing, zSpacing) = image.GetSpacing()
    (x0, y0, z0) = image.GetOrigin()

    center = [x0 + xSpacing * 0.5 * (xMin + xMax),
            y0 + ySpacing * 0.5 * (yMin + yMax),
            z0 + zSpacing * 0.5 * (zMin + zMax)]

    # Matrices for axial, coronal, sagittal, oblique view orientations
    axial = vtk.vtkMatrix4x4()
    axial.DeepCopy((1, 0, 0, center[0],
                    0, 1, 0, center[1],
                    0, 0, 1, center[2],
                    0, 0, 0, 1))
    
    coronal = vtk.vtkMatrix4x4()
    coronal.DeepCopy((1, 0, 0, center[0],
                      0, 0, 1, center[1],
                      0,-1, 0, center[2],
                      0, 0, 0, 1))

    sagittal = vtk.vtkMatrix4x4()
    sagittal.DeepCopy((0, 0,-1, center[0],
                       1, 0, 0, center[1],
                       0,-1, 0, center[2],
                       0, 0, 0, 1))

    oblique = vtk.vtkMatrix4x4()
    oblique.DeepCopy((1, 0, 0, center[0],
                      0, 0.866025, -0.5, center[1],
                      0, 0.5, 0.866025, center[2],
                      0, 0, 0, 1))

    # Extract a slice in the desired orientation
    #reslice = vtk.vtkImageReslice()
    #reslice.SetInputConnection(image)
    reslice.SetInput(image)
    reslice.SetOutputDimensionality(2)
    #reslice.SetResliceAxes(sagittal)
    reslice.SetResliceAxes(axial)
    reslice.SetInterpolationModeToLinear()

    # Create a greyscale lookup table
    #table = vtk.vtkLookupTable()
    table.SetRange(0, 2000) # image intensity range
    table.SetValueRange(0.0, 1.0) # from black to white
    table.SetSaturationRange(0.0, 0.0) # no color saturation
    table.SetRampToLinear()
    table.Build()

    # Map the image through the lookup table
    #color = vtk.vtkImageMapToColors()
    color.SetLookupTable(table)
    color.SetInputConnection(reslice.GetOutputPort())

    # Display the image
    #actor = vtk.vtkImageActor()
    actor.SetInput(color.GetOutput())

    #renderer = vtk.vtkRenderer()
    renderer.AddActor(actor)

    #window = vtk.vtkRenderWindow()
    window.AddRenderer(renderer)

    # Set up the interaction
    #interactorStyle = vtk.vtkInteractorStyleImage()
    #interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetInteractorStyle(interactorStyle)
    window.SetInteractor(interactor)
    window.Render()

    # Create callbacks for slicing the image
    #actions = {}
    actions["Slicing"] = 0

    interactorStyle.AddObserver("MouseMoveEvent", MouseMoveCallback)
    interactorStyle.AddObserver("LeftButtonPressEvent", ButtonCallback)
    interactorStyle.AddObserver("LeftButtonReleaseEvent", ButtonCallback)

    # Start interaction
    interactor.Start()

    return window, renderer, interactor, interactorStyle


