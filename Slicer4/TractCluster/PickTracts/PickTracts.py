import os
import unittest
from __main__ import vtk, qt, ctk, slicer

#
# PickTracts
#

class ModuleEventsList:
  """A helper class to manage a list of observers for a module. Call removeObservers on module exit."""
    
  def __init__(self):
    # list of events we have added to the scene
    # keep list of pairs: [observee,tag] so they can be removed easily
    self.observerTags = []

  def observeEvent(self, observee, eventName, eventCallback):
    tag = observee.AddObserver(eventName, eventCallback, 1.0)
    self.observerTags.append([observee,tag])

  def removeObservers(self):
    # remove all observers and reset
    for observee,tag in self.observerTags:
      observee.RemoveObserver(tag)
    self.observerTags = []


class ThreeDNodePicker:
  """A helper class to pick nodes in the 3D View. Can pick fiber bundles or models."""

  def __init__(self):
    # find the slicer objects needed for picking in the 3D window
    layoutManager = slicer.app.layoutManager()
    view = layoutManager.threeDWidget(0).threeDView()
    self.interactor = view.interactorStyle().GetInteractor()
    self.renderWindow = self.interactor.GetRenderWindow()
    # this assumes the first renderer is the one we want.
    self.renderer = self.renderWindow.GetRenderers().GetItemAsObject(0)
    
    # create a vtk picker to select from the scene
    # this prop picker will pick entire actors
    self.picker = vtk.vtkPropPicker()
    # type of node we want to pick
    self.setNodeTypeToFiberBundle()

    # event manager
    self.eventsList = ModuleEventsList()

  def setNodeType(self, nodeType):
    self.nodeType = nodeType
    self.nodePattern = '*'+nodeType+'*'

  def setNodeTypeToFiberBundle(self):
      self.setNodeType('FiberBundleNode')

  def setNodeTypeToModel(self):
      self.setNodeType('ModelNode')
        
  def startPicking(self, eventName="KeyPressEvent", eventCallback = None):
    # use this class callback function as default
    if not eventCallback:
        eventCallback = self.processEvent
    self.eventsList.observeEvent(self.interactor, eventName, eventCallback)

  def stopPicking(self):
    self.eventsList.removeObservers()
    
  def processEvent(self, caller=None, event=None):
      """ example event processing code. when using the class the
       process event callback should be defined when calling
       startPicking.  This is an example."""
      print "Using example processEvent. Please create your own based on this example code."
      key = self.interactor.GetKeySym()
      print "key:", key
      ret = self.pick()
      if ret:
          selection = self.getPickedNode()
          print "picked:", selection.GetName()
        
  def pick(self):    
    pos = self.interactor.GetEventPosition()
    sz = self.renderWindow.GetSize()
    ret = 0  
    if (pos[0] < sz[0]) & (pos[1] < sz[1]):
        ret = self.picker.PickProp(pos[0], pos[1], self.renderer)
    return ret  

  def getPickedNode(self):
    picked_node = None
    # find the picked mapper and its polydata
    actor = self.picker.GetActor()
    if actor:
        mapper = self.picker.GetActor().GetMapper()
        if mapper:
            pd_picked =  self.picker.GetActor().GetMapper().GetInput()
                
            # dictionary of the kind of nodes we are picking
            nodes = slicer.util.getNodes(self.nodePattern)
                
            # check which one was chosen in the pick
            for key in nodes:
                m = nodes[key]
                for n in range(m.GetNumberOfDisplayNodes()):
                    dn = m.GetNthDisplayNode(n)
                    if dn.GetOutputPolyData() == pd_picked:
                        picked_node = m
    # return the node object
    return picked_node          
  

class ModelDisplayHelper:
  """A helper class to handle all low-level display functionality for one or all models or fiber bundles. """

  def __init__(self):
    self.highlightColor = [.5, .2, .9]
    self.originalColors = []
    self.displayState = 'Lines'
    print "Find out tubes or not, and visible or not, from the data at init!"
    
  def setColor(self, m, color):
    print "set color", color
    for n in range(m.GetNumberOfDisplayNodes()):
        dn = m.GetNthDisplayNode(n)
        oldColor = dn.GetColor()
        print "old color", oldColor, "setting color", color
        if not (oldColor[0] == color[0] and oldColor[1] == color[1] and oldColor[2] == color[2]):
            print "old color and new are different"
            dn.SetColor(color)
            if not (oldColor == self.highlightColor):
                print "saving prev color", oldColor
                self.originalColors.append([m, oldColor])
            
  def resetColor(self, m):
    oldColor = None
    for node,color in self.originalColors:
        print "test node, color", color
        if node is m:
            print "found node on color list", node.GetName()
            oldColor = color
    if oldColor:
        print "found old color", oldColor
        for n in range(m.GetNumberOfDisplayNodes()):
            dn = m.GetNthDisplayNode(n)
            dn.SetColor(oldColor)
        
  def resetAllColors(self):
    # remove all observers and reset
    for node,color in self.originalColors:
      self.setColor(node,color)
    self.originalColors = []
        
  def allVisibleOn(self, nodeType):
    # dictionary of the kind of nodes we are picking
    if nodeType == 'FiberBundleNode' or nodeType == 'ModelNode':
        nodePattern = '*'+nodeType+'*'
        nodes = slicer.util.getNodes(nodePattern)
        for key in nodes:
            m = nodes[key]
            print "show model", m.GetName()
            self.visibleOn(m)
    else:
        print "NodeType", nodeType, "is unsupported"

  def allVisibleOff(self, nodeType):
    # dictionary of the kind of nodes we are picking
    if nodeType == 'FiberBundleNode' or nodeType == 'ModelNode':
        nodePattern = '*'+nodeType+'*'
        nodes = slicer.util.getNodes(nodePattern)
        for key in nodes:
            m = nodes[key]
            print "hide model", m.GetName()
            self.visibleOff(m)
    else:
        print "NodeType", nodeType, "is unsupported"
        
  def visibleOn(self,m):
    print "show model", m.GetClassName()
    className =  m.GetClassName()
    if className == 'vtkMRMLFiberBundleNode':
        #for n in range(m.GetNumberOfDisplayNodes()):
        #    dn = m.GetNthDisplayNode(n)
        #    dn.SetVisibility(1)
        if self.displayState == 'Tubes':
            m.GetTubeDisplayNode().SetVisibility(1)
            m.GetLineDisplayNode().SetVisibility(0)
            m.GetGlyphDisplayNode().SetVisibility(0)
        else:
            m.GetTubeDisplayNode().SetVisibility(0)
            m.GetLineDisplayNode().SetVisibility(1)
            m.GetGlyphDisplayNode().SetVisibility(0)
    else:
        # for fiber bundles this sets everything
        # visible, all three display formats.
        # so only use it for models.
        m.SetDisplayVisibility(1)

  def updateAll(self, nodeType):
    # called if we switch from tubes to lines, this seems not to have
    # a switch in MRML, as far as what display is active.
    # so update any visible nodes. This is likely
    # a bug, the desired behavior would be to store what is displayed
    # in the MRML file so it opens up with tubes again.
    print "updateAll", nodeType
    if nodeType == 'FiberBundleNode':
        nodePattern = '*'+nodeType+'*'
        nodes = slicer.util.getNodes(nodePattern)
        for key in nodes:
            m = nodes[key]
            print "update model", m.GetName()
            # update visibility of correct display type, if it's visible now
            if m.GetTubeDisplayNode().GetVisibility() or m.GetLineDisplayNode().GetVisibility():
                self.visibleOn(m)
    else:
        print "NodeType", nodeType, "is not in need of update"

        
  def visibleOff(self,m):
    print "hide model", m.GetName()
    print "type", m.GetClassName()  
    # for fiber bundles this sets everything
    # invisible, all three display formats
    m.SetDisplayVisibility(0)
        
  def highlightOff(self, m):
    # undo special highlight of the current model 
    # put back the original color
    print "highlightOff", m.GetName()
    self.resetColor(m)

  def highlightOn(self, m):
    # highlight display of the current model 
    print "highlightOn", m.GetName()
    self.setColor(m, self.highlightColor)

  def render(self):
    print "Render 3D"
    
class IDClipboard:
  """A helper class to maintain one interactively edited list or clipboard of nodes, and its display at the individual clipboard level. """

  def __init__(self):
    self.name = 'Selection'
    self.ID = 0
    self.contents = list()
    # display at the individual clipboard level
    self.display = ModelDisplayHelper()
    
  def remove(self, m):
    # Check if this model is already on the clipboard
    if m in self.contents:
        # if so, remove
        self.contents.remove(m)

  def add(self, m):
    # Check if this model node is already on the clipboard
    if m not in self.contents:
        if not self.contents:
            # If this is the first model, make it the default highlight color
            print "default color?"
            #self.ModelInteraction(clipboard,$id,clipboardColor) [Model($m,node) GetColor]
        # add model to clipboard list
        self.contents.append(m)
        # sort list
        self.contents.sort()

  def clear(self):
    # remove all from clipboard
    self.contents = list()
    self.highlightOff()
        
  def select(self, m):
    self.display.highlightOn(m)
    self.add(m)
    self.display.render()

  def deselect(self, m):
    self.display.highlightOff(m)
    self.remove(m)
    self.display.render()

  def highlight(self):
    for m in self.contents: 
        self.display.highlightOn(m)
    self.display.render()
    self.highlightMode = 1

  def highlightOff(self):
    for m in self.contents: 
        self.display.highlightOff(m)
    self.display.render()
    self.highlightMode = 0

  def visibleOn(self):
    for m in self.contents: 
        self.display.visibleOn(m)
    self.display.render()
    self.visibilityMode = 1

  def visibleOff(self):
    for m in self.contents: 
        self.display.visibleOff(m)
    self.display.render()
    self.visibilityMode = 0

  def ModelInteractionToggleClipboardHighlight(self):

    if self.ModelInteractionClipboardHighlight == 1:
        self.ModelInteractionClipboardHighlightOff()
    elif self.ModelInteractionClipboardHighlight == 0:
        self.ModelInteractionClipboardHighlightOn()

  def ModelInteractionToggleClipboardVisibility(self):

    if self.ModelInteractionClipboardVisibility == 1:
        self.ModelInteractionClipboardInVisibility()
    elif self.ModelInteractionClipboardVisibility == 0:
        self.ModelInteractionClipboardVisible()

        
class IDClipboardManager:
  """A helper class to manage many editable clipboards of nodes and their appearance in 3D, as well as their GUI display."""

  def __init__(self):
    self.clipboards = list()
    c = IDClipboard()
    self.clipboards.append(c)
    self.active = c
    self.displayMode = 'ModelNames'
    self.clipboardTextBox = None
    self.nameTextBox = None
    # choices: 'AllVisible', 'AllClipboardsInvisible', 'AllInvisible'
    self.visibilityState = 'AllInvisible'
    # choices: 'Visible', 'Invisible'
    self.clipboardVisibilityState = 'Invisible'
    # choices: 'Tubes', 'Lines', 'None'
    self.displayState = 'None'
        
   # type of node we want to pick
    self.setNodeTypeToFiberBundle()
    # display at the highest level (all models on or off)
    self.display = ModelDisplayHelper()
 
  def setNodeType(self, nodeType):
    self.nodeType = nodeType

  def setNodeTypeToFiberBundle(self):
    self.setNodeType('FiberBundleNode')

  def setNodeTypeToModel(self):
    self.setNodeType('ModelNode')
      
  def toggleDisplayState(self, state):
    print "TOGGLE DISPLAY IMPLEMENT ME"
    if state == 'Tubes':
        self.displayState = 'Tubes'
        self.display.displayState = 'Tubes'
        self.display.updateAll(self.nodeType)
        for c in self.clipboards:
            c.display.displayState = 'Tubes'
    elif state == 'Lines':
        self.displayState = 'Lines'         
        self.display.displayState = 'Lines'
        self.display.updateAll(self.nodeType)
        for c in self.clipboards:
            c.display.displayState = 'Lines'

  def allModelsTubes(self):
    print "display tubes for all models of current type in slicer"
    # for all fiber tracts regardless if on clipboards.
    self.updateRender()

  def allModelsLines(self):
    print "show lines for all models of current type in slicer"
    # for all fiber tracts regardless if on clipboards.
    self.updateRender()

  def toggleClipboardVisibilityState(self):
    # display of current clipboard
    print "TOGGLE CLIP VISIBILITY IMPLEMENT ME"
    print "STATE:", self.clipboardVisibilityState
    if self.clipboardVisibilityState == 'Visible':
        self.clipboardVisibilityState = 'Invisible'
        self.checkActive()
        self.active.visibleOff()
    elif self.clipboardVisibilityState == 'Invisible':
        self.clipboardVisibilityState = 'Visible'
        self.checkActive()
        self.active.visibleOn()                
    print "NEW STATE:", self.clipboardVisibilityState

  def toggleVisibilityState(self):
    print "TOGGLE VISIBILITY IMPLEMENT ME"
    print "STATE:", self.visibilityState
    
    if self.visibilityState == 'AllVisible':
        # show just unselected ones ("to-do" list of models)
        self.visibilityState = 'AllClipboardsInvisible'
        self.allClipboardsInvisible()
    elif self.visibilityState == 'AllClipboardsInvisible':
        # turn them all off (now can toggle clipboard to see just it)
        self.visibilityState = 'AllInvisible'
        self.allModelsInvisible() 
    elif self.visibilityState == 'AllInvisible':
        # turn all on (now can see everything at once)
        self.visibilityState = 'AllVisible'
        self.allModelsVisible()

    print "NEW STATE:", self.visibilityState
    
  def allModelsVisible(self):
    print "display all models of current type in slicer"
    # for all fiber tracts regardless if on clipboards.
    self.display.allVisibleOn(self.nodeType)
    self.updateRender()

  def allModelsInvisible(self):
    print "hide all models of current type in slicer"
    # for all fiber tracts regardless if on clipboards.
    self.display.allVisibleOff(self.nodeType)
    self.updateRender()
    
  def allClipboardsVisible(self):
    print "show all clipboards"
    for c in self.clipboards:
        c.visibleOn()
    self.updateRender()

  def allClipboardsInvisible(self):
    print "hide all clipboards"
    for c in self.clipboards:
        c.visibleOff()
    #self.ModelInteractionClipboardVisibility = 0
    self.updateRender()

  def updateGUI(self):
    # make sure the GUI is up to date
    print "call something to update a GUI"
    # clear the text box to put current list there

    if not self.clipboardTextBox:
        print "no text box widget set"
        return

    if not self.nameTextBox:
        print "no name text box widget set"
        return

    guiString = ''
    
    # Either put numbers or model names, depending on
    # the mode
    if self.displayMode is "ModelIDs":
        # fix off by one error with regards to slicer internal IDs vs. model numbers
        # and model IDs in MRML file
        for m in self.active.contents:
            #guiString.append( [Model($m,node) GetModelID]
            guiString = guiString + m.GetID() + '\n'
            #print "get model id", m
    elif self.displayMode is "ModelNames":
        for m in self.active.contents:
            #lappend guiString [Model($m,node) GetName]
            guiString = guiString + m.GetName() + '\n'
            #print "get model name", m
    print "add model names to GUI"

    #self.clipboardTextBox.text = self.active.contents
    self.clipboardTextBox.setText(guiString)
    
    # make sure the menu of clipboards is up to date
    #$ModelInteraction(clipboardMenuButton) config -text \
    #    $ModelInteraction(clipboard,$id,name)

    #self.menu $ModelInteraction(clipboardMenu)
    #$menu delete 0 end

  def updateRender(self):
    #self.Render3D()
    print "Render 3D"
       
  def remove(self, m):
    self.checkActive()
    self.active.remove(m)
    self.updateGUI()

  def select(self, m):
    self.checkActive()
    self.active.select(m)
    self.updateGUI()    

  def deselect(self, m):
    self.checkActive()
    self.active.deselect(m)
    self.updateGUI()    

  def checkActive(self):
    # in case we have no clipboard (it was deleted or 
    # does not exist yet), then add one.
    if not self.active:
        self.addClipboard()
        # make sure the GUI is up to date
        self.updateGUI()
              
  def getActive(self):
    self.checkActive()
    return self.active

  def setActive(self, c):
    if self.active is not c:
        # turn off highlight of previous clipboard
        self.active.highlightOff()
        self.active = c
        # turn on highlight of new clipboard
        self.active.highlightOn()
        # display entries for this clipboard
        self.updateGUI()
        # make these models visible
        self.active.visibleOn()
        self.updateRender()

  def addToClipboard(self, m):
    self.checkActive()
    self.active.add(m)
    self.updateGUI()

  def clearClipboard(self):
    self.checkActive()
    self.active.clear()
    self.updateGUI()
    self.updateRender()

  def deleteEmptyClipboards(self):
    for c in self.clipboards:
        if not c.contents:
            self.deleteClipboard(c)
    self.updateGUI()

  def deleteClipboard(self, c=None):
    if c is None:
        c = self.getActive()
        #else:
        # really this one will already be active
        # if delete is called from GUI
        self.setActive(c)
    # turn off display of this clipboard
    self.active.highlightOff()
    # remove clipboard from list
    self.clipboards.remove(c)
    # delete the object
    del c
    # update active clipboard?
    # which one should it be?
    # probably first on the list
    print "what more to update here after deletion?"
    # update gui and render (again?)
    # perhaps do this in the calling function
    
  def deleteAllClipboards(self):
    for c in self.clipboards:
        self.deleteClipboard(c)
    # clear list of clipboards
    self.clipboards = list()
    # this removes clipboard pulldown menu contents
    self.updateGUI()
    self.updateRender()

  def addClipboard(self):
    c = IDClipboard()
    self.clipboards.append(c)
    self.setActive(c)
    self.updateGUI()
    self.updateRender()

  def nameClipboard(self, name):
    print "give a clipboard a name that the user has requested"
    c = self.getActive()
    c.name = name
    print name
    # the name came from the GUI so no need to update it until create menu
    print "update menus"

  def exportClipboardsToTextFile(self, filename):
    f = open(filename,'w')
    for c in self.clipboards: 
        output = "GROUP " + c.name + "\n"
        f.write(output)
        for m in c.contents:
            f.write(str(m))
            # perhaps convert to model ID or filename or something. cluster index
            #puts $fid [Model($m,node) GetModelID]
    f.close()


class PickTracts:
  def __init__(self, parent):
    parent.title = "PickTracts" # TODO make this more human readable by adding spaces
    parent.categories = ["Examples"]
    parent.dependencies = []
    parent.contributors = ["Jean-Christophe Fillion-Robin (Kitware), Steve Pieper (Isomics)"] # replace with "Firstname Lastname (Org)"
    parent.helpText = """
    This is an example of scripted loadable module bundled in an extension.
    """
    parent.acknowledgementText = """
    This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc. and Steve Pieper, Isomics, Inc.  and was partially funded by NIH grant 3P41RR013218-12S1.
""" # replace with organization, grant and thanks.
    self.parent = parent

    # Add this test to the SelfTest module's list for discovery when the module
    # is created.  Since this module may be discovered before SelfTests itself,
    # create the list if it doesn't already exist.
    try:
      slicer.selfTests
    except AttributeError:
      slicer.selfTests = {}
    slicer.selfTests['PickTracts'] = self.runTest

  def runTest(self):
    tester = PickTractsTest()
    tester.runTest()

#
# qPickTractsWidget
#

class PickTractsWidget:
  def __init__(self, parent = None):
    if not parent:
      self.parent = slicer.qMRMLWidget()
      self.parent.setLayout(qt.QVBoxLayout())
      self.parent.setMRMLScene(slicer.mrmlScene)
    else:
      self.parent = parent
    self.layout = self.parent.layout()

    # one persistent logic object for current settings
    self.logic = PickTractsLogic()
    self.manager = self.logic.manager
    self.picker = ThreeDNodePicker()    
    # tell the picker what node type we are interested in
    self.picker.nodePattern = '*FiberBundleNode*'

    if not parent:
      self.setup()
      self.parent.show()

  def onNameChange(self):
    widget = self.nameTextBox
    print "hello"
    self.manager.nameClipboard(widget.text)


  def setup(self):
    
    # Instantiate and connect widgets ...

    # buttons: new group, change name, delete group
    # checkbox: show more info about models (number of lines etc)
    # checkbox: show names vs IDs?
    # radio buttons: handle models or tracts
    # button/icon: tubes vs lines
    #---
    # export to model hierarchy
    # import from model hierarchy
    # as above for text file
    
    # menu for all of the groups
    
    # editable name of the current group
    # (should be a popup if change name button is pressed?)
    self.nameTextBox = qt.QLineEdit()
    self.nameTextBox.text = self.manager.getActive().name
    self.layout.addWidget(self.nameTextBox)
    reloadFormLayout = qt.QFormLayout(self.nameTextBox)
    self.nameTextBox.connect('editingFinished()', self.onNameChange)
    self.manager.nameTextBox = self.nameTextBox

    # text box for list of models in the current group
    self.clipboardTextBox = qt.QTextBrowser()
    self.clipboardTextBox.setText("")
    self.layout.addWidget(self.clipboardTextBox)
    reloadFormLayout = qt.QFormLayout(self.clipboardTextBox)
    self.manager.clipboardTextBox = self.clipboardTextBox
    #self.nameTextBox.connect('editingFinished()', self.onNameChange)
    
    #
    # Reload and Test area
    #
    reloadCollapsibleButton = ctk.ctkCollapsibleButton()
    reloadCollapsibleButton.text = "Reload && Test"
    self.layout.addWidget(reloadCollapsibleButton)
    reloadFormLayout = qt.QFormLayout(reloadCollapsibleButton)

    # reload button
    # (use this during development, but remove it when delivering
    #  your module to users)
    self.reloadButton = qt.QPushButton("Reload")
    self.reloadButton.toolTip = "Reload this module."
    self.reloadButton.name = "PickTracts Reload"
    reloadFormLayout.addWidget(self.reloadButton)
    self.reloadButton.connect('clicked()', self.onReload)

    # reload and test button
    # (use this during development, but remove it when delivering
    #  your module to users)
    self.reloadAndTestButton = qt.QPushButton("Reload and Test")
    self.reloadAndTestButton.toolTip = "Reload this module and then run the self tests."
    reloadFormLayout.addWidget(self.reloadAndTestButton)
    self.reloadAndTestButton.connect('clicked()', self.onReloadAndTest)

    #
    # Parameters Area
    #
    parametersCollapsibleButton = ctk.ctkCollapsibleButton()
    parametersCollapsibleButton.text = "Parameters"
    self.layout.addWidget(parametersCollapsibleButton)

    # Layout within the dummy collapsible button
    parametersFormLayout = qt.QFormLayout(parametersCollapsibleButton)

    #
    # input volume selector
    #
    self.inputSelector = slicer.qMRMLNodeComboBox()
    self.inputSelector.nodeTypes = ( ("vtkMRMLScalarVolumeNode"), "" )
    self.inputSelector.addAttribute( "vtkMRMLScalarVolumeNode", "LabelMap", 0 )
    self.inputSelector.selectNodeUponCreation = True
    self.inputSelector.addEnabled = False
    self.inputSelector.removeEnabled = False
    self.inputSelector.noneEnabled = False
    self.inputSelector.showHidden = False
    self.inputSelector.showChildNodeTypes = False
    self.inputSelector.setMRMLScene( slicer.mrmlScene )
    self.inputSelector.setToolTip( "Pick the input to the algorithm." )
    parametersFormLayout.addRow("Input Volume: ", self.inputSelector)

    #
    # output volume selector
    #
    self.outputSelector = slicer.qMRMLNodeComboBox()
    self.outputSelector.nodeTypes = ( ("vtkMRMLScalarVolumeNode"), "" )
    self.outputSelector.addAttribute( "vtkMRMLScalarVolumeNode", "LabelMap", 0 )
    self.outputSelector.selectNodeUponCreation = False
    self.outputSelector.addEnabled = True
    self.outputSelector.removeEnabled = True
    self.outputSelector.noneEnabled = False
    self.outputSelector.showHidden = False
    self.outputSelector.showChildNodeTypes = False
    self.outputSelector.setMRMLScene( slicer.mrmlScene )
    self.outputSelector.setToolTip( "Pick the output to the algorithm." )
    parametersFormLayout.addRow("Output Volume: ", self.outputSelector)

    #
    # check box to trigger taking screen shots for later use in tutorials
    #
    self.enableScreenshotsFlagCheckBox = qt.QCheckBox()
    self.enableScreenshotsFlagCheckBox.checked = 0
    self.enableScreenshotsFlagCheckBox.setToolTip("If checked, take screen shots for tutorials. Use Save Data to write them to disk.")
    parametersFormLayout.addRow("Enable Screenshots", self.enableScreenshotsFlagCheckBox)

    #
    # scale factor for screen shots
    #
    self.screenshotScaleFactorSliderWidget = ctk.ctkSliderWidget()
    self.screenshotScaleFactorSliderWidget.singleStep = 1.0
    self.screenshotScaleFactorSliderWidget.minimum = 1.0
    self.screenshotScaleFactorSliderWidget.maximum = 50.0
    self.screenshotScaleFactorSliderWidget.value = 1.0
    self.screenshotScaleFactorSliderWidget.setToolTip("Set scale factor for the screen shots.")
    parametersFormLayout.addRow("Screenshot scale factor", self.screenshotScaleFactorSliderWidget)

    #
    # Apply Button
    #
    self.applyButton = qt.QPushButton("Apply")
    self.applyButton.toolTip = "Run the algorithm."
    self.applyButton.enabled = False
    parametersFormLayout.addRow(self.applyButton)

    # connections
    self.applyButton.connect('clicked(bool)', self.onApplyButton)
    self.inputSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)
    self.outputSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)

    # Add vertical spacer
    self.layout.addStretch(1)

  def cleanup(self):
    pass

  def onSelect(self):
    self.applyButton.enabled = self.inputSelector.currentNode() and self.outputSelector.currentNode()

  def onApplyButton(self):
    logic = PickTractsLogic()
    enableScreenshotsFlag = self.enableScreenshotsFlagCheckBox.checked
    screenshotScaleFactor = int(self.screenshotScaleFactorSliderWidget.value)
    print("Run the algorithm")
    logic.run(self.inputSelector.currentNode(), self.outputSelector.currentNode(), enableScreenshotsFlag,screenshotScaleFactor)

  def onReload(self,moduleName="PickTracts"):
    """Generic reload method for any scripted module.
    ModuleWizard will subsitute correct default moduleName.
    """
    import imp, sys, os, slicer

    widgetName = moduleName + "Widget"

    # reload the source code
    # - set source file path
    # - load the module to the global space
    filePath = eval('slicer.modules.%s.path' % moduleName.lower())
    p = os.path.dirname(filePath)
    if not sys.path.__contains__(p):
      sys.path.insert(0,p)
    fp = open(filePath, "r")
    globals()[moduleName] = imp.load_module(
        moduleName, fp, filePath, ('.py', 'r', imp.PY_SOURCE))
    fp.close()

    # rebuild the widget
    # - find and hide the existing widget
    # - create a new widget in the existing parent
    parent = slicer.util.findChildren(name='%s Reload' % moduleName)[0].parent().parent()
    for child in parent.children():
      try:
        child.hide()
      except AttributeError:
        pass
    # Remove spacer items
    item = parent.layout().itemAt(0)
    while item:
      parent.layout().removeItem(item)
      item = parent.layout().itemAt(0)

    # delete the old widget instance
    if hasattr(globals()['slicer'].modules, widgetName):
      getattr(globals()['slicer'].modules, widgetName).cleanup()

    # create new widget inside existing parent
    globals()[widgetName.lower()] = eval(
        'globals()["%s"].%s(parent)' % (moduleName, widgetName))
    globals()[widgetName.lower()].setup()
    setattr(globals()['slicer'].modules, widgetName, globals()[widgetName.lower()])

  def onReloadAndTest(self,moduleName="PickTracts"):
    try:
      self.onReload()
      evalString = 'globals()["%s"].%sTest()' % (moduleName, moduleName)
      tester = eval(evalString)
      tester.runTest()
    except Exception, e:
      import traceback
      traceback.print_exc()
      qt.QMessageBox.warning(slicer.util.mainWindow(), 
          "Reload and Test", 'Exception!\n\n' + str(e) + "\n\nSee Python Console for Stack Trace")
         
  def processSelectionEvent(self, caller=None, event=None):
    # events from the interactor
    if event == "KeyPressEvent":
        key = self.picker.interactor.GetKeySym()
        if key == 's' or key == 'S' or key == 'd' or key == 'D':
            # try to pick something in the scene
            ret = self.picker.pick()
            if ret:
                node = self.picker.getPickedNode()
                print "Picked: ", node.GetName()
                if node:
                    if key == 's' or key == 'S':
                        self.manager.select(node)
                    elif key == 'd' or key == 'D':
                        self.manager.deselect(node)         

        elif key == 'v' or key == 'V':
            self.manager.toggleVisibilityState() 
        elif key == 'c' or key == 'C':
            self.manager.toggleClipboardVisibilityState() 
        elif key == 't' or key == 'T':
            self.manager.toggleDisplayState('Tubes')          
        elif key == 'l' or key == 'L':
            self.manager.toggleDisplayState('Lines')         
    else:
      pass

    
  def enter(self):
    #   So that this module's event bindings don't conflict with other 
    #   modules, use our bindings only when the user is in this module.
    print "set up event handling and initialize display"
    #self.ModelInteractionAllVisibility = 0
    #self.ModelInteractionClipboardVisibility = 0
    #self.ModelInteractionClipboardHighlight = 0
    # add a mouse event for selecting a fiber bundle
    self.picker.startPicking(eventName="KeyPressEvent", eventCallback = self.processSelectionEvent)

  def exit(self):
    print "remove event handling and reset display"
    # remove the event we added
    self.picker.stopPicking()


#
# PickTractsLogic
#

class PickTractsLogic:
  """This class should implement all the actual 
  computation done by your module.  The interface 
  should be such that other python code can import
  this class and make use of the functionality without
  requiring an instance of the Widget
  """
  def __init__(self):
    self.manager = IDClipboardManager()
    
  def hasImageData(self,volumeNode):
    """This is a dummy logic method that 
    returns true if the passed in volume
    node has valid image data
    """
    if not volumeNode:
      print('no volume node')
      return False
    if volumeNode.GetImageData() == None:
      print('no image data')
      return False
    return True

  def delayDisplay(self,message,msec=1000):
    #
    # logic version of delay display
    #
    print(message)
    self.info = qt.QDialog()
    self.infoLayout = qt.QVBoxLayout()
    self.info.setLayout(self.infoLayout)
    self.label = qt.QLabel(message,self.info)
    self.infoLayout.addWidget(self.label)
    qt.QTimer.singleShot(msec, self.info.close)
    self.info.exec_()

  def takeScreenshot(self,name,description,type=-1):
    # show the message even if not taking a screen shot
    self.delayDisplay(description)

    if self.enableScreenshots == 0:
      return

    lm = slicer.app.layoutManager()
    # switch on the type to get the requested window
    widget = 0
    if type == -1:
      # full window
      widget = slicer.util.mainWindow()
    elif type == slicer.qMRMLScreenShotDialog().FullLayout:
      # full layout
      widget = lm.viewport()
    elif type == slicer.qMRMLScreenShotDialog().ThreeD:
      # just the 3D window
      widget = lm.threeDWidget(0).threeDView()
    elif type == slicer.qMRMLScreenShotDialog().Red:
      # red slice window
      widget = lm.sliceWidget("Red")
    elif type == slicer.qMRMLScreenShotDialog().Yellow:
      # yellow slice window
      widget = lm.sliceWidget("Yellow")
    elif type == slicer.qMRMLScreenShotDialog().Green:
      # green slice window
      widget = lm.sliceWidget("Green")

    # grab and convert to vtk image data
    qpixMap = qt.QPixmap().grabWidget(widget)
    qimage = qpixMap.toImage()
    imageData = vtk.vtkImageData()
    slicer.qMRMLUtils().qImageToVtkImageData(qimage,imageData)

    annotationLogic = slicer.modules.annotations.logic()
    annotationLogic.CreateSnapShot(name, description, type, self.screenshotScaleFactor, imageData)

  def run(self,inputVolume,outputVolume,enableScreenshots=0,screenshotScaleFactor=1):
    """
    Run the actual algorithm
    """

    self.delayDisplay('Running the aglorithm')

    self.enableScreenshots = enableScreenshots
    self.screenshotScaleFactor = screenshotScaleFactor

    self.takeScreenshot('PickTracts-Start','Start',-1)

    return True



        
class PickTractsTest(unittest.TestCase):
  """
  This is the test case for your scripted module.
  """

  def delayDisplay(self,message,msec=1000):
    """This utility method displays a small dialog and waits.
    This does two things: 1) it lets the event loop catch up
    to the state of the test so that rendering and widget updates
    have all taken place before the test continues and 2) it
    shows the user/developer/tester the state of the test
    so that we'll know when it breaks.
    """
    print(message)
    self.info = qt.QDialog()
    self.infoLayout = qt.QVBoxLayout()
    self.info.setLayout(self.infoLayout)
    self.label = qt.QLabel(message,self.info)
    self.infoLayout.addWidget(self.label)
    qt.QTimer.singleShot(msec, self.info.close)
    self.info.exec_()

  def setUp(self):
    """ Do whatever is needed to reset the state - typically a scene clear will be enough.
    """
    slicer.mrmlScene.Clear(0)

  def runTest(self):
    """Run as few or as many tests as needed here.
    """
    self.setUp()
    self.test_PickTracts1()

  def test_PickTracts1(self):
    """ Ideally you should have several levels of tests.  At the lowest level
    tests sould exercise the functionality of the logic with different inputs
    (both valid and invalid).  At higher levels your tests should emulate the
    way the user would interact with your code and confirm that it still works
    the way you intended.
    One of the most important features of the tests is that it should alert other
    developers when their changes will have an impact on the behavior of your
    module.  For example, if a developer removes a feature that you depend on,
    your test should break so they know that the feature is needed.
    """

    self.delayDisplay("Starting the test")
    #
    # first, get some data
    #
    import urllib
    downloads = (
        ('http://slicer.kitware.com/midas3/download?items=5767', 'FA.nrrd', slicer.util.loadVolume),
        )

    for url,name,loader in downloads:
      filePath = slicer.app.temporaryPath + '/' + name
      if not os.path.exists(filePath) or os.stat(filePath).st_size == 0:
        print('Requesting download %s from %s...\n' % (name, url))
        urllib.urlretrieve(url, filePath)
      if loader:
        print('Loading %s...\n' % (name,))
        loader(filePath)
    self.delayDisplay('Finished with download and loading\n')

    volumeNode = slicer.util.getNode(pattern="FA")
    logic = PickTractsLogic()
    self.assertTrue( logic.hasImageData(volumeNode) )
    self.delayDisplay('Test passed!')
