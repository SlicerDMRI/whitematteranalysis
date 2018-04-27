#!/usr/bin/env python
import os
import argparse
import urllib2
import tarfile
import sys

#-----------------
# Parse arguments
#-----------------
parser = argparse.ArgumentParser(
    description="Download a pre-provided anatomically curated fiber clustering white matter atlas.",
    epilog="Written by Fan Zhang, fzhang@bwh.harvard.edu",
    version='1.0')

parser.add_argument(
    'outputDirectory',
    help='The output directory, where the atlas will be downloaded, should be a new empty directory. It will be created if needed.')
parser.add_argument(
    '-atlas', action="store", dest="requested_atlas", type=str, 
    help='Name of the atlas. Currently, \'ORG-800FC-100HCP\' and \'ORG-2000FC-100HCP\' are available to download.')

args = parser.parse_args()

def download_file(url, output_file):
    
    try:
        print("Downloading: {0}".format(url))
        url_file = urllib2.urlopen(url)
        meta = url_file.info()
        file_size = int(meta.getheaders("Content-Length")[0])
        print("File size: {0} MB".format(file_size/1024.0/1024.0))

        output = open(output_file, 'wb')
        file_size_dl = 0
        block_sz = 8192
        while True:
            buffer = url_file.read(block_sz)
            if not buffer:
                break
            file_size_dl += len(buffer)
            output.write(buffer)
            p = float(file_size_dl) / file_size
            status = r"{0} bytes [{1:.2%}]".format(file_size_dl, p)
            status = status + chr(8)*(len(status)+1)
            sys.stdout.write(status)

        output.close()  
        print ' '      
    except:
        print 'Download fail!! Please check your network connection and rerun the script.'
        print sys.exc_info()
        exit()

def extract_requested_folder(members, sub_folder):

    for tarinfo in members:
        if sub_folder in tarinfo.name:
            yield tarinfo

def extract_atlas_from_archive(archive_file, output_dir, sub_folder=None, remove_after_extraction=False):

    try:
        zip_ref = tarfile.open(archive_file)

        if sub_folder is not None:
            members = extract_requested_folder(zip_ref, sub_folder)
            zip_ref.extractall(output_dir, members=members)
        else:
            zip_ref.extractall(output_dir)

        zip_ref.close()

        if remove_after_extraction is True:
            os.remove(archive_file)

    except:
        print 'Extraction fail!!'
        print sys.exc_info()
        exit()

outdir = os.path.abspath(args.outputDirectory)
if not os.path.exists(args.outputDirectory):
    print "Error: Output directory", args.outputDirectory, "does not exist, creating it."
    os.makedirs(outdir)

requested_atlas = args.requested_atlas

repo = 'https://github.com/SlicerDMRI/ORG-Atlases'
version = 'v1.0'

org_atlases_version_folder_name = 'org-atlases-' + version[1:]
org_atlases_version_folder = os.path.join(outdir, org_atlases_version_folder_name)

output_atlas_folder = os.path.join(org_atlases_version_folder, requested_atlas)
if os.path.exists(output_atlas_folder):
    print '\n<wm_download_org_atlas> Found an existing atlas at \''+output_atlas_folder+'\'. Remove this folder before redownload.'
    print ''
    exit()

archive_url = repo + '/archive/' + version + '.tar.gz'
downloaded_archive_file = os.path.join(outdir, version + '.tar.gz')

if requested_atlas == 'ORG-800FC-100HCP' or requested_atlas == 'ORG-2000FC-100HCP':
    FC_atlas_url = repo + '/releases/download/' + version + '/' + requested_atlas + '.tar.gz'
    REG_atlas_url = repo + '/releases/download/' + version + '/ORG-RegAtlas-100HCP.tar.gz'
    
    downloaded_FC_atlas_file = os.path.join(org_atlases_version_folder, requested_atlas + '.tar.gz')
    downloaded_Reg_atlas_file = os.path.join(org_atlases_version_folder, 'ORG-RegAtlas-100HCP.tar.gz')
else:
    print '<wm_download_org_atlas> ' + requested_atlas + 'is not available yet. Please check input.'
    print ''
    exit()

print ""
print "===== <wm_download_org_atlas.py> "
print "===== Release version   : ", version
print "===== Download from     : ", repo
print "===== Output directory  : ", outdir
print "===== Requested atlas   : ", requested_atlas

if requested_atlas == 'ORG-800FC-100HCP':
    print '** The ' +requested_atlas+ ' atlas is an anatomically curated white matter atlas, with annotated anatomical label provided for each cluster. '\
           + 'MRML files that are used to organize fiber clusters that belong to the same anatomical fiber tract are provided in the atlas.' 
elif requested_atlas == 'ORG-2000FC-100HCP':
    print '** The ' +requested_atlas+ ' atlas is provided for applications that can be benefit from a fine scale white matter parcellation. This an uncurated white matter atlas.'


print '<wm_download_org_atlas> Start downloading requested files...'

if requested_atlas == 'ORG-800FC-100HCP' or requested_atlas == 'ORG-2000FC-100HCP':
    print ' '
    download_file(archive_url, downloaded_archive_file)
    extract_atlas_from_archive(downloaded_archive_file, outdir, sub_folder=requested_atlas, remove_after_extraction=True)

    if not os.path.exists(os.path.join(org_atlases_version_folder, 'ORG-RegAtlas-100HCP/registration_atlas.vtk')):
        print ' '
        download_file(REG_atlas_url, downloaded_Reg_atlas_file)
        extract_atlas_from_archive(downloaded_Reg_atlas_file, org_atlases_version_folder, remove_after_extraction=True)
    else:
        print ' '
        print 'Registration atlas already donwloaded: ORG-RegAtlas-100HCP/registration_atlas.vtk'

    print ' '
    download_file(FC_atlas_url, downloaded_FC_atlas_file)
    extract_atlas_from_archive(downloaded_FC_atlas_file, org_atlases_version_folder, remove_after_extraction=True)

print ''
print '<wm_download_org_atlas> Successfully downloaded to', outdir
