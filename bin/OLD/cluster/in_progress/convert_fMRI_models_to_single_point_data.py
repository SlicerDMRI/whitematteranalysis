import glob
import numpy
import whitematteranalysis as wma
import pickle
import os
import vtk

indir = '/Users/odonnell/Dropbox/Coding/OUTPUTS/March2013/input_control_atlas/fMRI_vtk'
input_mask = "{0}/*".format(indir)
input_dirs = glob.glob(input_mask)

#input_dirs = input_dirs[0:3]

landmarks = list()
for input_dir in input_dirs:
    if os.path.isdir(input_dir):
        input_mask = "{0}/*.vtk".format(input_dir)
        input_poly_datas = glob.glob(input_mask)
        #print input_poly_datas
        center_list = list()
        for fname in input_poly_datas:
            print fname
            pd = wma.io.read_polydata(fname)
            pts =  pd.GetPoints()
            points_list = list()
            for pidx in range(pts.GetNumberOfPoints()):
                points_list.append(pts.GetPoint(pidx))

            center_list.append(numpy.mean(numpy.array(points_list),0))

        landmarks.append(center_list)
        
        subj_id = os.path.split(input_dir)[1]
        fname = 'landmarks_' + subj_id + '.p'
        pickle.dump( center_list, open( fname, "wb" ) )
    
#print landmarks
fname = 'landmarks_all.p'
pickle.dump( landmarks, open( fname, "wb" ) )
        
pd_new = vtk.vtkPolyData()
pts = vtk.vtkPoints()
colors = vtk.vtkIntArray()
for sidx in range(len(landmarks)):
    for lm in landmarks[sidx]:
        pts.InsertNextPoint(lm)
        colors.InsertNextTuple1(sidx)
        
pd_new.SetPoints(pts)
colors.SetName('Subject')
pd_new.GetPointData().AddArray(colors)
pd_new.GetPointData().SetActiveScalars('Subject')

glypher = vtk.vtkGlyph3D()
sphere = vtk.vtkSphereSource()
glypher.SetSourceConnection(sphere.GetOutputPort())
glypher.SetInputData(pd_new)
glypher.SetScaleModeToDataScalingOff()
glypher.SetColorModeToColorByScalar()
glypher.Update()

ren = wma.render.render(glypher.GetOutput(), data_mode='Point')

wma.io.write_polydata(glypher.GetOutput(), 'point_cloud.vtp')
