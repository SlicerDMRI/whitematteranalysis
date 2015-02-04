#!/usr/bin/env python
"""
Converts a PickTracts group file into a blob of XML suitable for Slicer's
Model Hierarchy and inserts it into a scene mrml file. 

Usage: 
    picktracts_converter.py <groupfile> <scene.mrml> <output.mrml>
"""
from docopt import docopt

if __name__ == '__main__': 
    arguments = docopt(__doc__)
    groupfilepath = arguments['<groupfile>']
    scenepath     = arguments['<scene.mrml>']
    outputpath    = arguments['<output.mrml>']

    groupfile = open(groupfilepath) 

    modeltext = "" # blob of xml to insert into scene

    # build a dictionary that maps group name to a list of tuples: 
    # (tract id,tract name)
    groups = {}

    # read in a picktracts group file
    # modified from PickTracts.py
    string = groupfile.readline().strip()
    while string != "END":
        if string == "GROUP NAME":
            group_name = groupfile.readline().strip()
        if string == "TRACT ID":
            tract_id = groupfile.readline().strip()
        if string == "TRACT NAME":
            tract_name = groupfile.readline().strip()
            tracts = groups.get(group_name,[])
            tract = (tract_id,tract_name)
            # don't include duplicates
            if tract not in tracts:
                groups[group_name] = tracts + [tract]

        string = groupfile.readline().strip()
    groupfile.close()

    # generate XML for a MRML scene
    for i, group_name in enumerate(groups.keys()):
        modeldisplayid   = "PickTractsDisplayID{}".format(i)
        modeldisplayname = "PickTractsDisplayName{}".format(i)
        modelparentid    = "PickTractsModelHierarchyNode{}".format(i)

        modeldisplaytext="""<ModelDisplay 
            id="{}"
            name="{}"
            hideFromEditors="true" selectable="true" selected="false" color="0.5 0.5 0.5" edgeColor="0 0 0" selectedColor="1 0 0" selectedAmbient="0.4" ambient="0" diffuse="1" selectedSpecular="0.5" specular="0" power="1" opacity="1" pointSize="1" lineWidth="1" representation="2" lighting="true" interpolation="1" shading="true" visibility="true" edgeVisibility="false" clipping="false" sliceIntersectionVisibility="false" sliceIntersectionThickness="1" frontfaceCulling="false" backfaceCulling="true" scalarVisibility="false" vectorVisibility="false" tensorVisibility="false" interpolateTexture="false" scalarRangeFlag="2" autoScalarRange="true" scalarRange="0 100"/>""".format(modeldisplayid,modeldisplayname)
        model_parent="""<ModelHierarchy id="{}" name="{}" displayNodeID="{}" hideFromEditors="false" selectable="true" selected="false" sortingValue="0" allowMultipleChildren="true" expanded="true"/>""".format(modelparentid,group_name,modeldisplayid)

        tractnodes = []
        for j, tractinfo in enumerate(groups[group_name]):
            tract_id, tract_name = tractinfo
            node_name = "PickTractsModelHierarchy_{}_{}".format(i,j)
            node_id   = "PickTractsModelHierarchyNode{}_{}".format(i,j)

            tract_node="""<ModelHierarchy id="{}" name="{}" hideFromEditors="true" selectable="true" selected="false" parentNodeRef="{}" associatedNodeRef="{}" sortingValue="1" allowMultipleChildren="true" expanded="true"/>""".format(node_id,node_name,modelparentid,tract_id)
            tractnodes.append(tract_node)

        modeltext += "\n".join([modeldisplaytext,model_parent] + tractnodes)

    # append XML to scene and output
    scene = open(scenepath).read()
    scene = scene.replace('</MRML>',modeltext+'\n</MRML>') 
    open(outputpath,'w').write(scene)
