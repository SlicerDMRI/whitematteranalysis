
nodePattern = '*FiberBundleNode*'

# get fiber bundle nodes
nodes = slicer.util.getNodes(nodePattern)


# find the slicer objects needed for picking
layoutManager = slicer.app.layoutManager()
view = layoutManager.threeDWidget(0).threeDView()
interactor = view.interactorStyle().GetInteractor()

# add a mouse event for selecting a fiber bundle
interactor.AddObserver("KeyPressEvent", processEvent, 1.0)


# this picker here is not the same one we get in the event.
# this is strange
# it has no pick points and a different memory address
picker = interactor.GetPicker()
ren = picker.GetRenderer()
ren.GetPickResultProps()

# clear picking before we start
#picker.DeletePickList()
#ren.GetPickResultProps()



# display what we picked
renWin = view.renderWindow()
renWin.Render()

def processEvent(caller=None, event=None):
    print 'hello 13 ---'
    # QVTKInteractor
    #print caller
    print event
    # this gets a picker that works.
    picker = caller.GetPicker()
    print picker
    ren = picker.GetRenderer()
    prp = ren.GetPickResultProps()
    print "props picked", prp.GetNumberOfItems()
    #for idx in range(prp.GetNumberOfItems()):
        #print prp.GetItemAsObject(idx)
        
    print '==='
    
    #print picker
    #picker = interactor.GetPicker()
    if 0:
        print picker 
        print interactor.GetPicker()
        ren = picker.GetRenderer()
        prp = ren.GetPickResultProps()
        print "props picked", prp.GetNumberOfItems()
        print picker.GetProp3D()
        print '---'


def processEvent(caller=None, event=None):
    print 'hello 13 ---'
    lm = slicer.app.layoutManager()
    td = lm.threeDWidget(0)
    ms = vtk.vtkCollection()
    td.getDisplayableManagers(ms)
    for i in range(ms.GetNumberOfItems()): 
        m = ms.GetItemAsObject(i)
        if m.GetClassName() == "vtkMRMLModelDisplayableManager":
            print i, m.GetClassName()
            print m.GetPickedNodeID()
            print m.GetPickedRAS()
        
ask steve call vtkproppicker in model displayable manager
and it's not exposed to python
getpickednodeid
maybe getpicked RAS ,e tc

best practices
fbnode managed by model displayable manager

rely on model displayable manager
code from python you can get it
fb node set pickable true
end of markups module wiki page info for developers
markups 3d fiducial displayable manager

can pass this through to the picker from vtkmrmlmodeldisplayablemanager
interactor.GetLastEventPosition()
translation to ras to get pick

helper function on one of the displayable managers to say what does this belong to

key press events n modules loadable markups


        

picked node: vtkMRMLFiberBundleLineDisplayNode1
picked cell: 380
