num_points = list()
num_lines = list()
for fname in input_poly_datas:
    print fname
    pd = wma.io.read_polydata(fname)
    num_points.append(pd.GetNumberOfPoints())
    num_lines.append(pd.GetNumberOfLines())
    


      
num_points = numpy.array(num_points)
num_lines = numpy.array(num_lines)

scipy.stats.ttest_ind(num_points[g0], num_points[g1])
scipy.stats.ttest_ind(num_lines[g0], num_lines[g1])
