# -*- coding: utf-8 -*-

""" mrml.py

output in MRML scene file format (XML)

"""
import os


def write(pd_filenames, colors, filename, ratio=1.0):
    writer = WriteMRML()
    writer.write(pd_filenames, colors, filename, ratio)
    del writer
    
class WriteMRML:
    
    def __init__(self):
        self.header = '<MRML  version="Slicer4" userTags="">\n'
        self.footer = '</MRML>\n'
        self.indent = ' '
        self.node_id = 0
        self.props_id = 0
        
    def write(self, pd_filenames, colors, filename, ratio=1.0):
        #print "converting colors to strings"
        color_list = list()
        for cidx in range(len(colors)):
            col = str(colors[cidx,0]/256.0) + " " + str(colors[cidx,1]/256.0) + " " + str(colors[cidx,2]/256.0)
            color_list.append(col)
        #print color_list
        print(f"<{os.path.basename(__file__)}> Writing f{len(pd_filenames)} filenames in MRML scene: {filename}")
        f = open(filename, "w")
        f.write(self.header)
        for pidx in range(len(pd_filenames)):
            name = os.path.splitext(os.path.split(pd_filenames[pidx])[1])[0]
            self.write_node(pd_filenames[pidx], color_list[pidx], name, f, ratio)
        f.write(self.footer)
        f.close()
        
    def write_node(self, pd_fname, color, name, f, ratio):

        self.node_id += 1
        idx = self.node_id
        props = list()
        
        f.write(self.indent)
        f.write("<FiberBundleStorage\n")
        f.write(self.indent)
        f.write(self.indent)
        f.write(f"id=\"vtkMRMLFiberBundleStorageNode{str(idx)}\"  ")
        f.write("name=\"FiberBundleStorage\"  hideFromEditors=\"true\"  selectable=\"true\"  selected=\"false\"  ")
        f.write(f"fileName=\"{pd_fname}\"  ")
        f.write("useCompression=\"1\"  readState=\"0\"  writeState=\"0\" ></FiberBundleStorage>")
        f.write("\n")
                
        self.props_id += 1
        idxp = self.props_id
        props.append(idxp)

        f.write(self.indent)
        f.write("<FiberBundleLineDisplayNode\n")
        f.write(self.indent)
        f.write(self.indent)
        f.write(f"id=\"vtkMRMLFiberBundleLineDisplayNode{str(idx)}\"  ")
        f.write("name=\"FiberBundleLineDisplayNode\"  hideFromEditors=\"true\"  selectable=\"true\"  selected=\"false\"  ")
        f.write(f"color=\"{color}\"  ")
        f.write("edgeColor=\"0 0 0\"  selectedColor=\"1 0 0\"  selectedAmbient=\"0.4\"  ambient=\"0\"  diffuse=\"1\"  selectedSpecular=\"0.5\"  specular=\"0\"  power=\"1\"  opacity=\"1\"  pointSize=\"1\"  lineWidth=\"1\"  representation=\"2\"  lighting=\"true\"  interpolation=\"1\"  shading=\"true\"  visibility=\"true\"  edgeVisibility=\"false\"  clipping=\"false\"  sliceIntersectionVisibility=\"false\"  sliceIntersectionThickness=\"1\"  frontfaceCulling=\"false\"  backfaceCulling=\"false\"  scalarVisibility=\"false\"  vectorVisibility=\"false\"  tensorVisibility=\"false\"  interpolateTexture=\"false\"  autoScalarRange=\"true\"  scalarRange=\"0 1\"  colorNodeID=\"vtkMRMLColorTableNodeRainbow\"   colorMode =\"0\"  ")
        f.write(f"DiffusionTensorDisplayPropertiesNodeRef=\"vtkMRMLDiffusionTensorDisplayPropertiesNode{str(idxp)}\"  ")
        f.write("></FiberBundleLineDisplayNode>")
        f.write("\n")

        self.props_id += 1
        idxp = self.props_id
        props.append(idxp)
        
        f.write(self.indent)
        f.write("<FiberBundleTubeDisplayNode\n")
        f.write(self.indent)
        f.write(self.indent)
        f.write(f"id=\"vtkMRMLFiberBundleTubeDisplayNode{str(idx)}\"  ")
        f.write("name=\"FiberBundleTubeDisplayNode\"  hideFromEditors=\"true\"  selectable=\"true\"  selected=\"false\"  ")
        f.write(f"color=\"{color}\"  ")
        f.write("edgeColor=\"0 0 0\"  selectedColor=\"1 0 0\"  selectedAmbient=\"0.4\"  ambient=\"0.25\"  diffuse=\"0.8\"  selectedSpecular=\"0.5\"  specular=\"0.25\"  power=\"20\"  opacity=\"1\"  pointSize=\"1\"  lineWidth=\"1\"  representation=\"2\"  lighting=\"true\"  interpolation=\"1\"  shading=\"true\"  visibility=\"false\"  edgeVisibility=\"false\"  clipping=\"false\"  sliceIntersectionVisibility=\"false\"  sliceIntersectionThickness=\"1\"  frontfaceCulling=\"false\"  backfaceCulling=\"false\"  scalarVisibility=\"false\"  vectorVisibility=\"false\"  tensorVisibility=\"false\"  interpolateTexture=\"false\"  autoScalarRange=\"true\"  scalarRange=\"0 1\"  colorNodeID=\"vtkMRMLColorTableNodeRainbow\"   colorMode =\"0\"  ")
        f.write(f"DiffusionTensorDisplayPropertiesNodeRef=\"vtkMRMLDiffusionTensorDisplayPropertiesNode{str(idxp)}\"  ")
        f.write("tubeRadius =\"0.5\"  tubeNumberOfSides =\"6\" ></FiberBundleTubeDisplayNode>")
        f.write("\n")
        
        self.props_id += 1
        idxp = self.props_id
        props.append(idxp)

        f.write(self.indent)
        f.write("<FiberBundleGlyphDisplayNode\n")
        f.write(self.indent)
        f.write(self.indent)
        f.write(f"id=\"vtkMRMLFiberBundleGlyphDisplayNode{str(idx)}\"  ")
        f.write("name=\"FiberBundleGlyphDisplayNode\"  hideFromEditors=\"true\"  selectable=\"true\"  selected=\"false\"  ")
        f.write(f"color=\"{color}\"  ")
        f.write("edgeColor=\"0 0 0\"  selectedColor=\"1 0 0\"  selectedAmbient=\"0.4\"  ambient=\"0\"  diffuse=\"1\"  selectedSpecular=\"0.5\"  specular=\"0\"  power=\"1\"  opacity=\"1\"  pointSize=\"1\"  lineWidth=\"1\"  representation=\"2\"  lighting=\"true\"  interpolation=\"1\"  shading=\"true\"  visibility=\"false\"  edgeVisibility=\"false\"  clipping=\"false\"  sliceIntersectionVisibility=\"false\"  sliceIntersectionThickness=\"1\"  frontfaceCulling=\"false\"  backfaceCulling=\"false\"  scalarVisibility=\"false\"  vectorVisibility=\"false\"  tensorVisibility=\"false\"  interpolateTexture=\"false\"  autoScalarRange=\"true\"  scalarRange=\"0 1\"  colorNodeID=\"vtkMRMLColorTableNodeRainbow\"   colorMode =\"0\"  ")
        f.write(f"DiffusionTensorDisplayPropertiesNodeRef=\"vtkMRMLDiffusionTensorDisplayPropertiesNode{str(idxp)}\"  ")
        f.write("twoDimensionalVisibility=\"false\" ></FiberBundleGlyphDisplayNode>")
        f.write("\n")


        f.write(self.indent)
        f.write("<FiberBundle\n")
        f.write(self.indent)
        f.write(self.indent)
        f.write(f"id=\"vtkMRMLFiberBundleNode{str(idx)}\"  ")
        f.write(f"name=\"{name}\"  ")
        f.write("hideFromEditors=\"false\"  selectable=\"true\"  selected=\"false\"  ")
        f.write(f"displayNodeRef=\"vtkMRMLFiberBundleLineDisplayNode{str(idx)}  ")
        f.write(f"vtkMRMLFiberBundleTubeDisplayNode{str(idx)}  ")
        f.write(f"vtkMRMLFiberBundleGlyphDisplayNode{str(idx)}\"  ")
        f.write(f"storageNodeRef=\"vtkMRMLFiberBundleStorageNode{str(idx)}\"  ")
        f.write("references=\"display:")
        f.write(f"vtkMRMLFiberBundleLineDisplayNode{str(idx)}  ")
        f.write(f"vtkMRMLFiberBundleTubeDisplayNode{str(idx)}  ")
        f.write(f"vtkMRMLFiberBundleGlyphDisplayNode{str(idx)};")
        f.write("storage:")
        f.write(f"vtkMRMLFiberBundleStorageNode{str(idx)};\"  ")
        f.write("userTags=\"\"  SelectWithAnnotationNode=\"0\"  SelectionWithAnnotationNodeMode=\"0\"  ")
        f.write(f"SubsamplingRatio=\"{str(ratio)}\" ></FiberBundle>")
        f.write("\n")

        for prop in props:
            self.write_prop_node(prop, f)
        

    def write_prop_node(self, idx, f):

        f.write(self.indent)
        f.write("<DiffusionTensorDisplayProperties\n")
        f.write(self.indent)
        f.write(self.indent)
        f.write(f"id=\"vtkMRMLDiffusionTensorDisplayPropertiesNode{str(idx)}\"  ")
        f.write("name=\"DiffusionTensorDisplayPropertiesNode\"  description=\"A user defined colour table, use the editor to specify it\"  hideFromEditors=\"true\"  selectable=\"true\"  selected=\"false\"  userTags=\"\" type=\"13\" numcolors=\"0\"  glyphGeometry=\"2\"  colorGlyphBy=\"3\"  glyphScaleFactor=\"50\"  glyphEigenvector=\"1\"  glyphExtractEigenvalues=\"1\"  lineGlyphResolution=\"20\"  tubeGlyphRadius=\"0.1\"  tubeGlyphNumberOfSides=\"4\"  ellipsoidGlyphThetaResolution=\"9\"  ellipsoidGlyphPhiResolution=\"9\"  superquadricGlyphGamma=\"1\"  superquadricGlyphThetaResolution=\"6\"  superquadricGlyphPhiResolution=\"6\" ></DiffusionTensorDisplayProperties>")
        f.write("\n")
