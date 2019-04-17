from __future__ import print_function
import os
import unittest
import time
from __main__ import vtk, qt, ctk, slicer

#
# PickTracts
#

try:
    _fromUtf8 = qt.QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s

#class UI_Write_FileDialog()

#class Write_FileDialog(qt.QFileDialog):
  
#  def __init__(self):
#      qt.QFileDialog.__init__(self)
#      self.ui = UI_Write_FileDialog()
#      self.ui.setupui(self)

class UI_Newgroup_Dialog(object):
  
  def setupui(self,Dialog):
      Dialog.setObjectName(_fromUtf8("Dialog"))
      Dialog.resize(387, 153)
      self.QDialogLayout = qt.QFormLayout()    
      self.temLineEdit = qt.QLineEdit()
      self.temSaveButton = qt.QPushButton("Create New Group")
      Dialog.setLayout(self.QDialogLayout)
      self.QDialogLayout.addRow("Edit a Group Name:",self.temLineEdit)
      self.QDialogLayout.addWidget(self.temSaveButton)
      #self.temSaveButton.connect('clicked()',self.createnewgroup)

class Newgroup_Dialog(qt.QDialog):
  
  def __init__(self):
      qt.QDialog.__init__(self)
      self.ui = UI_Newgroup_Dialog()
      self.ui.setupui(self)

class UI_Changegroupname_Dialog(object):
  
  def setupui(self,Dialog):
      Dialog.setObjectName(_fromUtf8("Dialog"))
      Dialog.resize(387, 153)
      self.QDialogLayout = qt.QFormLayout()    
      self.temComboBox = qt.QComboBox()
      self.temLineEdit = qt.QLineEdit()
      self.temLineEdit.insert(self.temComboBox.itemText(0))
      self.temComboBox.connect('currentIndexChanged(int)',self.comboboxchanged)
      #self.temComboBox.editable = True
      self.temSaveButton = qt.QPushButton("Change Group Name")
      Dialog.setLayout(self.QDialogLayout)
      self.QDialogLayout.addRow("Choose a Group Name:",self.temComboBox)
      self.QDialogLayout.addRow("Edit a New Name:",self.temLineEdit)
      self.QDialogLayout.addWidget(self.temSaveButton)
      #self.temSaveButton.connect('clicked()',self.createnewgroup)
  
  def comboboxchanged(self):
      self.temLineEdit.clear()
      self.temLineEdit.insert(self.temComboBox.currentText)
class Changegroupname_Dialog(qt.QDialog):
  
  def __init__(self):
      qt.QDialog.__init__(self)
      self.ui = UI_Changegroupname_Dialog()
      self.ui.setupui(self)
  
  def additem(self,groupname):
      self.ui.temComboBox.addItem(groupname)

  def cleargroupnamelist(self):
      self.ui.temComboBox.clear()
  
class ModuleEventsList(object):
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


class ThreeDNodePicker(object):
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
    # example event processing code. when using the class the process event callback should be defined when calling startPicking.  This is an example.
    print("Using example processEvent. Please create your own based on this example code.")
    key = self.interactor.GetKeySym()
    print("key:", key)
    ret = self.pick()
    if ret:
        selection = self.getPickedNode()
        print("picked:", selection.GetName())
        
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
            slicer.mrmlScene.StartState(slicer.mrmlScene.BatchProcessState)
            for key in nodes:
                m = nodes[key]
                for n in range(m.GetNumberOfDisplayNodes()):
                    dn = m.GetNthDisplayNode(n)
                    if dn.GetOutputPolyData() == pd_picked:
                        picked_node = m
            slicer.mrmlScene.EndState(slicer.mrmlScene.BatchProcessState)
    # return the node object
    return picked_node   

  def findPickedNode(self,Name):
    picked_node = None
    # find the picked mapper and its polydata
                
    # dictionary of the kind of nodes we are picking
    nodes = slicer.util.getNodes(self.nodePattern)
                
    # check which one was chosen in the pick
    slicer.mrmlScene.StartState(slicer.mrmlScene.BatchProcessState)
    for key in nodes:
        m = nodes[key]
        if m.GetName() == Name:
            picked_node = m
    slicer.mrmlScene.EndState(slicer.mrmlScene.BatchProcessState)
    # return the node object
    return picked_node        
  

class ModelDisplayHelper(object):
  """A helper class to handle all low-level display functionality for one or all models or fiber bundles. """

  def __init__(self):
    self.highlightColor = [.5, .2, .9]
    self.originalColors = []
    self.displayState = 'Lines'
    print("Find out tubes or not, and visible or not, from the data at init!")
    
  def setColor(self, m, color):
    print("set color", color)
    for n in range(m.GetNumberOfDisplayNodes()):
        dn = m.GetNthDisplayNode(n)
        oldColor = dn.GetColor()
        print("old color", oldColor, "setting color", color)
        if not (oldColor[0] == color[0] and oldColor[1] == color[1] and oldColor[2] == color[2]):
            print("old color and new are different")
            dn.SetColor(color)
            if not (oldColor == self.highlightColor):
                print("saving prev color", oldColor)
                self.originalColors.append([m, oldColor])
            
  def resetColor(self, m):
    oldColor = None
    for node,color in self.originalColors:
        print("test node, color", color)
        if node is m:
            print("found node on color list", node.GetName())
            oldColor = color
    if oldColor:
        print("found old color", oldColor)
        slicer.mrmlScene.StartState(slicer.mrmlScene.BatchProcessState)
        for n in range(m.GetNumberOfDisplayNodes()):
            dn = m.GetNthDisplayNode(n)
            dn.SetColor(oldColor)
        slicer.mrmlScene.EndState(slicer.mrmlScene.BatchProcessState)
        
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
        slicer.mrmlScene.StartState(slicer.mrmlScene.BatchProcessState)
        for key in nodes:
            m = nodes[key]
            print("show model", m.GetName())
            self.visibleOn(m)
        slicer.mrmlScene.EndState(slicer.mrmlScene.BatchProcessState)
    else:
        print("NodeType", nodeType, "is unsupported")

  def allVisibleOff(self, nodeType):
    # dictionary of the kind of nodes we are picking
    if nodeType == 'FiberBundleNode' or nodeType == 'ModelNode':
        nodePattern = '*'+nodeType+'*'
        nodes = slicer.util.getNodes(nodePattern)
        slicer.mrmlScene.StartState(slicer.mrmlScene.BatchProcessState)
        for key in nodes:
            m = nodes[key]
            print("hide model", m.GetName())
            self.visibleOff(m)
        slicer.mrmlScene.EndState(slicer.mrmlScene.BatchProcessState)
    else:
        print("NodeType", nodeType, "is unsupported")
        
  def visibleOn(self,m):
    print("show model", m.GetClassName())
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
    print("updateAll", nodeType)
    if nodeType == 'FiberBundleNode':
        nodePattern = '*'+nodeType+'*'
        nodes = slicer.util.getNodes(nodePattern)
        slicer.mrmlScene.StartState(slicer.mrmlScene.BatchProcessState)
        for key in nodes:
            m = nodes[key]
            print("update model", m.GetName())
            # update visibility of correct display type, if it's visible now
            if m.GetTubeDisplayNode().GetVisibility() or m.GetLineDisplayNode().GetVisibility():
                self.visibleOn(m)
        slicer.mrmlScene.EndState(slicer.mrmlScene.BatchProcessState)
    else:
        print("NodeType", nodeType, "is not in need of update")

        
  def visibleOff(self,m):
    print("hide model", m.GetName())
    print("type", m.GetClassName())
    # for fiber bundles this sets everything
    # invisible, all three display formats
    m.SetDisplayVisibility(0)
        
  def highlightOff(self, m):
    # undo special highlight of the current model 
    # put back the original color
    print("highlightOff", m.GetName())
    self.resetColor(m)

  def highlightOn(self, m):
    # highlight display of the current model 
    print("highlightOn", m.GetName())
    self.setColor(m, self.highlightColor)

  def render(self):
    print("Render 3D")

class TractsGroup(object):
  """A helper class to maintain one saved tracts group. """
  
  def __init__(self):
    self.groupname = None
    self.contents = list()
    self.display = ModelDisplayHelper()

  def setgroupname(self, name):
    self.groupname = name

  def addcontents(self, tract):
    self.contents.append(tract)


class TractsGroupManager(object):
  """A helper class to manage saved tracts groups. """
  
  def __init__(self):
    self.groups = list()
    self.groupnames = list()
    #self.newTractsTextBox = None
    self.oldTractsTextBox = None
    self.groupnameSelector = None 
    self.group = TractsGroup()
    self.display = ModelDisplayHelper()
  
  def savetogroup(self,namestring,clipmanager):
    if namestring not in self.groupnames:
        self.groupnames.append(namestring)
        self.addgroup(namestring,clipmanager)
    else:
        self.addtogroup(namestring,clipmanager)
    #clipmanager.updateGUI()

  def addgroup(self,namestring,clipmanager):
    g = TractsGroup()
    g.setgroupname(namestring)
    for tract in clipmanager.newt.contents:
        print("1",tract)
        g.addcontents(tract)
    while len(clipmanager.newt.contents) > 0:
        clipmanager.newt.contents.pop()
    self.groups.append(g)
    print("2",g)

  def addtogroup(self,namestring,clipmanager):
    for g in self.groups:
        if g.groupname == namestring:
            for tract in clipmanager.newt.contents:
                print("3",tract)
                g.addcontents(tract)
            while len(clipmanager.newt.contents) > 0:
                clipmanager.newt.contents.pop()
            print("4",g)
            return
  
  def listgroup(self,namestring,clipmanager):
    while len(clipmanager.oldt.contents) > 0:
        clipmanager.oldt.contents.pop()
    for g in self.groups:
        if g.groupname == namestring:
            for tract in g.contents:
                print("5",tract)
                clipmanager.oldt.add(tract)
            #clipmanager.updateGUI()
            return

  def viewgroup(self,namestring,clipmanager):
    for g in self.groups:
        if g.groupname == namestring:
            self.display.allVisibleOff(clipmanager.nodeType)
            for tract in g.contents:
                self.display.visibleOn(tract)
            self.display.render()
            return
  
  def delete(self,tract):
    for g in self.groups:
        for t in g.contents:
            if t == tract:
                g.contents.remove(t)

  def changegroupname(self,groupname,newgroupname):
    while len(self.groupnames) > 0:
        self.groupnames.pop()
    for g in self.groups:
        if groupname == g.groupname:
            g.groupname = newgroupname
        self.groupnames.append(g.groupname)

          
class IDTracts(object):
  """A helper class to maintain ungrouped selected tracts. """

  def __init__(self):
    self.contents = list()
  
  def add(self, m):
    # Check if this model node is already on the list
    if m not in self.contents:
        # add model to ner selected tracts list
        self.contents.append(m)
        # sort list
        self.contents.sort()

  def remove(self, m):
    # Check if this model is already on the list
    if m in self.contents:
        # if so, remove
        self.contents.remove(m)
    
class IDClipboard(object):
  """A helper class to maintain one interactively edited list or clipboard of nodes, and its display at the individual clipboard level. """

  def __init__(self):
    self.name = 'All'
    self.ID = 0
    #self.group = list()
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
            print("default color?")
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

class IDClipboardManager(object):
  """A helper class to manage many editable clipboards of nodes and their appearance in 3D, as well as their GUI display."""

  def __init__(self):
    self.clipboards = list()
    c = IDClipboard()
    self.newt = IDTracts()
    self.oldt = IDTracts()
    self.clipboards.append(c)
    self.active = c
    self.displayMode = 'ModelNames'
    self.clipboardTextBox = None
    #self.newTractsTextBox = None
    self.oldTractsTextBox = None
    #self.nameTextBox = None
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
    print("TOGGLE DISPLAY IMPLEMENT ME")
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
    print("display tubes for all models of current type in slicer")
    # for all fiber tracts regardless if on clipboards.
    self.updateRender()

  def allModelsLines(self):
    print("show lines for all models of current type in slicer")
    # for all fiber tracts regardless if on clipboards.
    self.updateRender()

  def toggleClipboardVisibilityState(self):
    # display of current clipboard
    print("TOGGLE CLIP VISIBILITY IMPLEMENT ME")
    print("STATE:", self.clipboardVisibilityState)
    if self.clipboardVisibilityState == 'Visible':
        self.clipboardVisibilityState = 'Invisible'
        self.checkActive()
        self.active.visibleOff()
    elif self.clipboardVisibilityState == 'Invisible':
        self.clipboardVisibilityState = 'Visible'
        self.checkActive()
        self.active.visibleOn()                
    print("NEW STATE:", self.clipboardVisibilityState)

  def toggleVisibilityState(self):
    print("TOGGLE VISIBILITY IMPLEMENT ME")
    print("STATE:", self.visibilityState)
    
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

    print("NEW STATE:", self.visibilityState)
    
  def allModelsVisible(self):
    print("display all models of current type in slicer")
    # for all fiber tracts regardless if on clipboards.
    self.display.allVisibleOn(self.nodeType)
    self.updateRender()

  def allModelsInvisible(self):
    print("hide all models of current type in slicer")
    # for all fiber tracts regardless if on clipboards.
    self.display.allVisibleOff(self.nodeType)
    self.updateRender()
    
  def allClipboardsVisible(self):
    print("show all clipboards")
    for c in self.clipboards:
        c.visibleOn()
    self.updateRender()

  def allClipboardsInvisible(self):
    print("hide all clipboards")
    for c in self.clipboards:
        c.visibleOff()
    #self.ModelInteractionClipboardVisibility = 0
    self.updateRender()

  def updateGUI(self):
    # make sure the GUI is up to date
    print("call something to update a GUI")
    # clear the text box to put current list there

    #if not self.clipboardTextBox:
        #print "no text box widget set"
        #return

    #if not self.nameTextBox:
        #print "no name text box widget set"
        #return

    guiString = ''
    guiString2 = ''
    guiString3 = ''

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
    print("add model names to GUI")

    # Either put numbers or model names, depending on
    # the mode
    #if self.displayMode is "ModelIDs":
        # fix off by one error with regards to slicer internal IDs vs. model numbers
        # and model IDs in MRML file
        #for m in self.newt.contents:
            #guiString.append( [Model($m,node) GetModelID]
            #guiString2 = guiString2 + m.GetID() + '\n'
            #print "get model id", m
    #elif self.displayMode is "ModelNames":
        #for m in self.newt.contents:
            #lappend guiString [Model($m,node) GetName]
            #guiString2 = guiString2 + m.GetName() + '\n'
            #print "get model name", m

    # Either put numbers or model names, depending on
    # the mode
    if self.displayMode is "ModelIDs":
        # fix off by one error with regards to slicer internal IDs vs. model numbers
        # and model IDs in MRML file
        for m in self.oldt.contents:
            #guiString.append( [Model($m,node) GetModelID]
            guiString3 = guiString3 + m.GetID() + '\n'
            #print "get model id", m
    elif self.displayMode is "ModelNames":
        for m in self.oldt.contents:
            #lappend guiString [Model($m,node) GetName]
            guiString3 = guiString3 + m.GetName() + '\n'
            #print "get model name", m

    #self.clipboardTextBox.text = self.active.contents
    #self.clipboardTextBox.setText(guiString)
    #self.newTractsTextBox.setText(guiString2)
    self.oldTractsTextBox.setText(guiString3)
    
    # make sure the menu of clipboards is up to date
    #$ModelInteraction(clipboardMenuButton) config -text \
    #    $ModelInteraction(clipboard,$id,name)

    #self.menu $ModelInteraction(clipboardMenu)
    #$menu delete 0 end

  def updateRender(self):
    #self.Render3D()
    print("Render 3D")
       
  def remove(self, m):
    self.checkActive()
    self.active.remove(m)
    self.updateGUI()

  def select(self, m, TractsGroupManager, namestring):
    self.checkActive()
    self.active.select(m)
    self.newt.add(m)
    TractsGroupManager.savetogroup(namestring,self)
    TractsGroupManager.listgroup(namestring,self)
    self.updateGUI()    

  def deselect(self, m, tractsgroupmanager,picktractswidget):
    self.checkActive()
    self.active.deselect(m)
    self.newt.remove(m)
    tractsgroupmanager.delete(m)
    picktractswidget.ongroupList()
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
    print("what more to update here after deletion?")
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
    print("give a clipboard a name that the user has requested")
    c = self.getActive()
    c.name = name
    #c.group.append(name)
    print(name)
    # the name came from the GUI so no need to update it until create menu
    print("update menus")

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


class PickTracts(object):
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

class PickTractsWidget(object):
  def __init__(self, parent = None):
    if not parent:
      self.parent = slicer.qMRMLWidget()
      self.parent.setLayout(qt.QVBoxLayout())
      self.parent.setMRMLScene(slicer.mrmlScene)
    else:
      self.parent = parent
    self.layout = self.parent.layout()

    # whther or not to edit tracts group name first time
    self.flag = 1
    self.flag_existgroup = 0
    
    self.qmessage1 = qt.QMessageBox()
    self.qmessage2 = qt.QMessageBox()
    self.QDialog = Newgroup_Dialog()
    self.QDialog2 = Changegroupname_Dialog()
    #self.QLineEdit = qt.QLineEdit(self.QDialog)
    self.TractsGroupManager = TractsGroupManager()
    self.newgroupname = ""
    self.QDir = qt.QDir()    

    self.filedialog1 = qt.QFileDialog()
    self.filedialog2 = qt.QFileDialog()
    # group name list
    #self.group = list()

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
    print("hello",widget.text)
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
    # Select Tracts area
    #
    selectCollapsibleButton = ctk.ctkCollapsibleButton()
    selectCollapsibleButton.text = "Select Tracts"
    #self.layout.addWidget(selectCollapsibleButton)
    selectFormLayout = qt.QFormLayout(selectCollapsibleButton)

    # editable name of the current group
    # (should be a popup if change name button is pressed?)
    #self.nameTextBox = qt.QLineEdit()
    #self.nameTextBox.text = self.manager.getActive().name
    #self.nameTextBox.readOnly = True
    #selectFormLayout.addRow(self.nameTextBox)
    #reloadFormLayout = qt.QFormLayout(self.nameTextBox)
    #self.nameTextBox.connect('editingFinished()', self.onNameChange)
    #self.manager.nameTextBox = self.nameTextBox

    # text box for list of models in the current group
    #selectFormLayout.addRow("All Selected Tracts:",self.nameTextBox)
    #self.clipboardTextBox = qt.QTextBrowser()
    #self.clipboardTextBox.setText("")
    #selectFormLayout.addRow(self.clipboardTextBox)
    #reloadFormLayout = qt.QFormLayout(self.clipboardTextBox)
    #self.manager.clipboardTextBox = self.clipboardTextBox
    #self.nameTextBox.connect('editingFinished()', self.onNameChange)
    #self.TractsGroupManager.clipboardTextBox = self.clipboardTextBox

    # Build and Manage Tracts Groups 
    self.groupnameSelector = qt.QComboBox()
    self.groupnameSelector.setToolTip("Edit or choose a tracts group name.")
    #selectFormLayout.addRow("Edit or Choose a Tracts Group Name:",self.groupnameSelector)
    self.groupnameSelector.setEditable(1)
    self.groupnameSelector.setDuplicatesEnabled(0)
    #self.manager.nameTextBox = self.groupnameSelector
    self.TractsGroupManager.groupnameSelector = self.groupnameSelector
    #self.groupnameSelector.connect('currentIndexChanged(int)',self.groupnameChanged)
    #self.groupnameSelector.connect('editTextChanged(QString)',self.groupnameActivated)

    # New Selected Tracts Table
    #self.newTractsTextBox = qt.QTextBrowser()
    #self.newTractsTextBox.setToolTip("List new selected tracts.")
    #selectFormLayout.addRow("Ungrouped Selected Tracts:", self.newTractsTextBox)
    #self.manager.newTractsTextBox = self.newTractsTextBox
    #self.TractsGroupManager.newTractsTextBox = self.newTractsTextBox
 
    #Save Button
    self.saveGroupButton = qt.QPushButton("Save Selected Tracts in the Group")
    self.saveGroupButton.toolTip = "Save selected tracts in the group."
    #selectFormLayout.addRow(self.saveGroupButton)
    self.saveGroupButton.connect('clicked(bool)',self.onSaveTracts)  

    #
    # Manage Tracts area
    #
    manageCollapsibleButton = ctk.ctkCollapsibleButton()
    manageCollapsibleButton.text = "Pick and Group Tracts"
    self.layout.addWidget(manageCollapsibleButton)
    manageFormLayout = qt.QFormLayout(manageCollapsibleButton)


    # Group List
    self.groupChoose = qt.QComboBox()
    self.groupChoose.setToolTip("Choose a tracts group name.")
    self.groupChoose.editable = False   
    manageFormLayout.addRow("Edit or Choose a Tracts Group:", self.groupChoose)
    self.groupChoose.addItem("")
    self.groupChoose.addItem("Create a New Group")
    self.groupChoose.addItem("Change a Group Name")
    self.groupChoose.connect('currentIndexChanged(int)',self.ongroupList)
    x = int
    self.groupChoose.connect('activated(int)',self.onnewgroup)

    # List Button
    self.listGroupButton = qt.QPushButton("List Tracts in the Group")
    self.listGroupButton.toolTip = "List tracts in the group."
    #manageFormLayout.addRow(self.listGroupButton)
    self.listGroupButton.connect('clicked(bool)',self.onListTracts) 

    # Tracts in the Selected Group
    self.oldTractsTextBox = qt.QTextBrowser()
    self.oldTractsTextBox.setToolTip("List tracts in the selected group.")
    manageFormLayout.addRow(self.oldTractsTextBox)
    self.manager.oldTractsTextBox = self.oldTractsTextBox

    # View Button
    self.viewGroupButton = qt.QPushButton("View Tracts in the Group")
    self.viewGroupButton.toolTip = "View tracts in the group."
    manageFormLayout.addRow(self.viewGroupButton)
    self.viewGroupButton.connect('clicked(bool)',self.onViewTracts) 

    loadsaveCollapsibleButton = ctk.ctkCollapsibleButton()
    loadsaveCollapsibleButton.text = "Load and Save"
    self.layout.addWidget(loadsaveCollapsibleButton)
    loadsaveFormLayout = qt.QFormLayout(loadsaveCollapsibleButton)

    # Input Button
    self.inputGroupButton = qt.QPushButton("Load Tract Groups File")
    self.inputGroupButton.toolTip = "Input selected tracts groups from a txt file."
    loadsaveFormLayout.addRow(self.inputGroupButton)
    self.inputGroupButton.connect('clicked(bool)',self.onReadFile) 

    # Output Button
    self.outputGroupButton = qt.QPushButton("Save Tract Groups File")
    self.outputGroupButton.toolTip = "Output selected tracts groups as a txt file."
    loadsaveFormLayout.addRow(self.outputGroupButton)
    self.outputGroupButton.connect('clicked(bool)',self.onWriteFile) 

    #
    # Global Visibility State area
    #
    stateCollapsibleButton = ctk.ctkCollapsibleButton()
    stateCollapsibleButton.text = "Visibility Settings"
    self.layout.addWidget(stateCollapsibleButton)
    stateFormLayout = qt.QFormLayout(stateCollapsibleButton)

    #
    # Global Visilibility State Selector
    #
    self.globalviewSelector = qt.QComboBox()
    self.globalviewSelector.setToolTip("Choose a global visibility state.")
    self.globalviewSelector.addItem("All Tracts Visible")
    self.globalviewSelector.addItem("Only Selected Tracts Visible")
    self.globalviewSelector.addItem("Only Unselected Tracts Visible")
    self.globalviewSelector.addItem("All Tracts Unvisible")
    stateFormLayout.addRow("Visibility of All Tracts:",self.globalviewSelector)
    self.globalviewSelector.setCurrentIndex(0)

    # View Apply Button
    #
    self.viewapplyButton = qt.QPushButton("Apply")
    self.viewapplyButton.toolTip = "Apply the visibility state."
    self.viewapplyButton.enabled = True
    stateFormLayout.addRow(self.viewapplyButton)
    self.viewapplyButton.connect('clicked(bool)',self.onViewapply)

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
    globals()[moduleName] = imp.load_module(moduleName, fp, filePath, ('.py', 'r', imp.PY_SOURCE))
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
    globals()[widgetName.lower()] = eval('globals()["%s"].%s(parent)' % (moduleName, widgetName))
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
                print("Picked: ", node.GetName())
                if node:
                    if key == 's' or key == 'S':
                        namestring = self.groupChoose.currentText
                        if namestring == "":
                            self.qmessage1.about(self.parent,"Error!","No Group Name!")
                            return
                        self.manager.select(node, self.TractsGroupManager, namestring)
                    elif key == 'd' or key == 'D':
                        self.manager.deselect(node,self.TractsGroupManager,self)         

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
    print("set up event handling and initialize display")
    #self.ModelInteractionAllVisibility = 0
    #self.ModelInteractionClipboardVisibility = 0
    #self.ModelInteractionClipboardHighlight = 0
    # add a mouse event for selecting a fiber bundle
    self.picker.startPicking(eventName="KeyPressEvent", eventCallback = self.processSelectionEvent)

  def exit(self):
    print("remove event handling and reset display")
    # remove the event we added
    self.picker.stopPicking()

  def onViewapply(self):
    print("ViewApply")
    viewselect = self.globalviewSelector.currentIndex
    print(viewselect)
    # Visibility State: All Tracts Visible
    if viewselect == 0:
        self.manager.allModelsVisible()
        self.manager.visibilityState = 'AllVisible'
        print("NEW STATE:", self.globalviewSelector.currentText)
    # Visibility State: Only Selected Tracts Visible
    if viewselect == 1:
        self.manager.allModelsInvisible()
        self.manager.visibilityState = 'AllInvisible'
        self.manager.checkActive()
        self.manager.active.visibleOn()
        self.manager.clipboardVisibilityState = 'Visible'                
        print("NEW STATE:", self.globalviewSelector.currentText)
    # Visibility State: Only Unselected Tracts Visible
    if viewselect == 2:
        self.manager.allModelsVisible()
        self.manager.visibilityState = 'AllClipboardsInvisible'
        self.manager.checkActive()
        self.manager.active.visibleOff()
        self.manager.clipboardVisibilityState = 'Invisible'
        print("NEW STATE:", self.globalviewSelector.currentText)
    # Visibility State: All Tracts Unvisible
    if viewselect == 3:
        self.manager.allModelsInvisible()
        self.manager.visibilityState = 'AllInvisible'
        print("NEW STATE:", self.globalviewSelector.currentText)

  def onSaveSelectedTracts(self,m):
    NameString = self.groupChoose.currentText
    if NameString == "":
        self.qmessage1.about(self.parent,"Error!","No Group Name!")
        return
    self.TractsGroupManager.savetogroup(NameString,self.manager,m)

      

  def onSaveTracts(self):
    NameString = self.groupnameSelector.currentText
    if NameString == "":
        self.qmessage1.about(self.parent,"Error!","No Group Name!")
        return
    #if len(self.manager.newt.contents) == 0:
        #self.qmessage2.about(self.parent,"Error!","No Ungrouped Tract!")
        #return
    # if it is the first time to edit the groupname
    if self.flag == 1:
        self.groupnameSelector.addItem(NameString)
        tI = self.groupnameSelector.findText(NameString)           
        self.groupnameSelector.setCurrentIndex(tI)
        #self.TractsGroupManager.addgroup()
        self.flag = 0
        self.groupnameSelector.setCurrentIndex(0)
        self.groupnameSelector.clearEditText()
        self.TractsGroupManager.savetogroup(NameString,self.manager)
        self.groupChoose.addItem(NameString)
    else:
        self.groupnameSelector.clearEditText()
        if self.groupnameSelector.findText(NameString) == -1:
            self.groupnameSelector.addItem(NameString) 
            tI = self.groupnameSelector.findText(NameString)           
            self.groupnameSelector.setCurrentIndex(tI)
            #self.TractsGroupManager.addgroup()
            self.groupnameSelector.setCurrentIndex(0)
            self.groupnameSelector.clearEditText()
            self.groupChoose.addItem(NameString)
        self.TractsGroupManager.savetogroup(NameString,self.manager)

  def ongroupList(self):
    print("Here!")
    namestring = self.groupChoose.currentText
    print(namestring)
    self.TractsGroupManager.listgroup(namestring,self.manager)
    self.manager.updateGUI()    

  def onListTracts(self):
    namestring = self.groupChoose.currentText
    print(namestring)
    self.TractsGroupManager.listgroup(namestring,self.manager)            

  def onViewTracts(self):
    namestring = self.groupChoose.currentText
    self.TractsGroupManager.viewgroup(namestring,self.manager)  
    
  def groupnameChanged(self):
    print("HaHa")
    string = ""
    self.flag_existgroup = 0
    for i in self.TractsGroupManager.groups:
        if i.groupname == self.groupnameSelector.currentText:
            self.flag_existgroup = 1
            for s in i.contents:
                self.clipboardTextBox.setText(i.contents)
    if self.flag_existgroup == 0:
        self.clipboardTextBox.setText("")
  

  def groupnameActivated(self):
    print("activated")
    self.clipboardTextBox.setText("")

  def onnewgroup(self):
    print("CH")
    if self.groupChoose.currentText == "Create a New Group":
        #self.QDialog.show()
        #self.QDialogLayout = qt.QVBoxLayout()
        #self.temLineEdit = qt.QLineEdit()
        #self.temSaveButton = qt.QPushButton("Create")
        #self.QDialog.setLayout(self.QDialogLayout)
        #self.QDialogLayout.addWidget(self.temLineEdit)
        #self.QDialogLayout.addWidget(self.temSaveButton)
        #self.newgroupname = ""
        #self.newgroupname = self.QDialog.temLienEdit.text
        self.QDialog.ui.temSaveButton.connect('clicked()',self.createnewgroup)
        self.QDialog.exec_()
        print(self.newgroupname)
        print("New")
    if self.groupChoose.currentText == "Change a Group Name":
        print(self.groupChoose.count)
        for i in range(0,self.groupChoose.count-2):
            groupname = ""
            groupname = self.groupChoose.itemText(i)
            print("copy",groupname)
            self.QDialog2.additem(groupname)
        self.QDialog2.ui.temSaveButton.connect('clicked()',self.changegroupname)
        self.QDialog2.exec_()
        #self.QDialog2.show()
        self.QDialog2.cleargroupnamelist()
        self.oldTractsTextBox.setText("")
 
  def changegroupname(self):
    i = self.QDialog2.ui.temComboBox.currentIndex
    groupname = self.QDialog2.ui.temComboBox.currentText
    newgroupname = self.QDialog2.ui.temLineEdit.text
    self.groupChoose.setItemText(i,newgroupname)
    self.QDialog2.ui.temComboBox.setItemText(i,newgroupname)
    print("changegroupname",newgroupname)
    if groupname != newgroupname:
        self.TractsGroupManager.changegroupname(groupname,newgroupname)
   

  def createnewgroup(self):
    self.newgroupname = ""
    self.newgroupname = self.QDialog.ui.temLineEdit.text
    self.groupChoose.addItem(self.newgroupname)
    i = self.groupChoose.findText("")
    if i != -1:
        self.groupChoose.removeItem(i)
    i = self.groupChoose.findText("Create a New Group")
    self.groupChoose.removeItem(i)
    i = self.groupChoose.findText("Change a Group Name")
    self.groupChoose.removeItem(i)
    self.groupChoose.addItem("Create a New Group")
    self.groupChoose.addItem("Change a Group Name")
    self.groupChoose.setCurrentIndex(0)
    self.QDialog.ui.temLineEdit.clear()
    self.oldTractsTextBox.clear()

  def onReadFile(self):
    self.filedialog2.setFileMode(1)
    self.filedialog2.setAcceptMode(0)
    self.filedialog2.setNameFilter("*.txt")
    self.filedialog2.fileSelected.connect(self.openfilecall)
    self.filedialog2.show()

  def onWriteFile(self):
    self.filedialog1.setFileMode(0)
    self.filedialog1.setAcceptMode(1)
    self.filedialog1.setNameFilter("*.txt")
    self.filedialog1.setDefaultSuffix("txt")
    self.filedialog1.fileSelected.connect(self.savefilecall)
    self.filedialog1.show()
                
  def openfilecall(self):
    pathfile = self.filedialog2.selectedFiles()
    groupfile = file(pathfile[0],"r")
    groupfile.close()
    groupfile = file(pathfile[0],"r")
    groupnum = 0
    string = groupfile.readline()
    while string != "END":
        if string[0:10] == "GROUP NAME":
            groupnum = groupnum+1
            groupstring = groupfile.readline()
            n = len(groupstring) - 1
            self.newgroupname = ""
            self.newgroupname = groupstring[0:n]
            self.groupChoose.addItem(self.newgroupname)
            i = self.groupChoose.findText("")
            if i != -1:
            	self.groupChoose.removeItem(i)
    	    i = self.groupChoose.findText("Create a New Group")
            self.groupChoose.removeItem(i)
            i = self.groupChoose.findText("Change a Group Name")
            self.groupChoose.removeItem(i)
            self.groupChoose.addItem("Create a New Group")
            self.groupChoose.addItem("Change a Group Name")
            self.groupChoose.setCurrentIndex(0)
            self.QDialog.ui.temLineEdit.clear()
            self.oldTractsTextBox.clear()
        string = groupfile.readline()
    print(groupnum)
    groupfile.close()
    groupfile = file(pathfile[0],"r")
    string = groupfile.readline()
    while string != "END":
        if string[0:10] == "GROUP NAME":
            groupstring = groupfile.readline()
            n = len(groupstring) - 1
        if string[0:8] == "TRACT ID":
            tractID = groupfile.readline()
            n1 = len(tractID) - 1
        if string[0:10] == "TRACT NAME":
            tractNAME = groupfile.readline()
            n2 = len(tractNAME) - 1
	    node = self.picker.findPickedNode(tractNAME[0:n2])
            print(tractID[0:n1])
            print(tractNAME[0:n2])
            i = self.groupChoose.findText(groupstring[0:n])
            self.groupChoose.setCurrentIndex(i)
            self.manager.select(node, self.TractsGroupManager, groupstring)
        string = groupfile.readline()
    groupfile.close()

  def savefilecall(self):
    pathfile = self.filedialog1.selectedFiles()
    groupfile = file(pathfile[0],"w")
    groupfile.write("TRACTS GROUPING STATE\n")
    curTime = time.strftime("%Y-%m-%d %X", time.localtime(time.time()))
    groupfile.write("EDITED TIME ")
    groupfile.write(curTime)
    groupfile.write("\n")
    groupfile.write("\n")
    for g in self.TractsGroupManager.groups:
        groupfile.write("GROUP NAME")
        groupfile.write("\n")
        groupfile.write(g.groupname)
        groupfile.write("\n")
        groupfile.write("COLOR ")
        groupfile.write("\n")
        for t in g.contents:
            groupfile.write("TRACT ID ")
            groupfile.write("\n")
            groupfile.write(t.GetID())
            groupfile.write("\n")
            groupfile.write("TRACT NAME ")
            groupfile.write("\n")
            groupfile.write(t.GetName())
            groupfile.write("\n")
        groupfile.write("\n")  
    groupfile.write("END")

#
# PickTractsLogic
#

class PickTractsLogic(object):
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
