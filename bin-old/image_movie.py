
def image_movie(image, slice1, slice2, skip, orientation='ax'):
    limits = [numpy.min(image), numpy.max(image)]
    print "limits:", limits
    #plt.figure()
    for slice in range(slice1,slice2,skip):
        plt.imshow(image[:,:,slice])
        #plt.imshow(image[:,slice,:])
        #plt.imshow(image[slice,:,:])
        plt.clim(limits)
        plt.draw()
        fname = 'slice_{:04}.png'.format(slice)
        #plt.savefig('slice%d.png' % slice)
        plt.savefig(fname)
        time.sleep(0.2)



    


