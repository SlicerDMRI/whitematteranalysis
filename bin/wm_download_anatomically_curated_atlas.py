#!/usr/bin/env python
import os
import argparse
import urllib.request, urllib.error, urllib.parse
import zipfile
import sys
import ssl

def main():
    if (not os.environ.get('PYTHONHTTPSVERIFY', '') and getattr(ssl, '_create_unverified_context', None)): 
        ssl._create_default_https_context = ssl._create_unverified_context
    
    #-----------------
    # Parse arguments
    #-----------------
    parser = argparse.ArgumentParser(
        description="Download a pre-provided anatomically curated fiber clustering white matter atlas.",
        epilog="Written by Fan Zhang, fzhang@bwh.harvard.edu")
    
    parser.add_argument(
        'outputDirectory',
        help='The output directory, where the atlas will be downloaded, should be a new empty directory. It will be created if needed.')
    parser.add_argument(
        '-atlas', action="store", dest="requested_atlas", type=str, 
        help='Name of the atlas. Currently, \'ORG-800FC-100HCP\' and \'ORG-2000FC-100HCP\' are available to download.')
    
    args = parser.parse_args()
    
    def download_file(url, output_file):
    
        try:
            print(f"* Downloading: {url}")
            url_file = urllib.request.urlopen(url)
            meta = url_file.info()
            file_size = int(meta.get_all("Content-Length")[0])
            print(f"File size: {file_size/1024.0/1024.0} MB")
    
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
                status = fr"{file_size_dl} bytes [{p:.2%}]"
                status = status + chr(8)*(len(status)+1)
                sys.stdout.write(status)
    
            output.close()  
            print('')      
        except:
            print('[ERRIR] Download fails! Please check your network connection and rerun the script.')
            print(sys.exc_info())
            exit()
    
    def extract_from_zip(archive_file, output_dir, remove_after_extraction=False):
    
        try:
            zip_ref = zipfile.ZipFile(archive_file, 'r')
            zip_ref.extractall(output_dir)
            zip_ref.close()
        except:
            print('[ERROR] Extraction fails!')
            print(sys.exc_info())
            exit()
    
        if remove_after_extraction:
            os.remove(archive_file)
    
    outdir = os.path.abspath(args.outputDirectory)
    if not os.path.exists(args.outputDirectory):
        print("<wm_download_org_atlas> Output directory", args.outputDirectory, "does not exist, creating it.")
        os.makedirs(outdir)
    
    requested_atlas = args.requested_atlas
    
    repo = 'https://zenodo.org/record/4156927/'
    version = 'v1.2'
    
    org_atlases_version_folder_name = 'ORG-Atlases-' + version[1:]
    org_atlases_version_folder = os.path.join(outdir, org_atlases_version_folder_name)
    if not os.path.exists(org_atlases_version_folder):
        os.makedirs(org_atlases_version_folder)
    
    if requested_atlas == 'ORG-800FC-100HCP' or requested_atlas == 'ORG-2000FC-100HCP':
        FC_atlas_url = repo + 'files/' + requested_atlas + '.zip?download=1'
        REG_atlas_url = repo + 'files/' + 'ORG-RegAtlas-100HCP.zip?download=1'
    
        population_T1_url = repo + 'files/' + '100HCP-population-mean-T1.nii.gz?download=1'
        population_T2_url = repo + 'files/' + '100HCP-population-mean-T2.nii.gz?download=1'
        population_b0_url = repo + 'files/' + '100HCP-population-mean-b0.nii.gz?download=1'
        
        downloaded_FC_atlas_file = os.path.join(org_atlases_version_folder, requested_atlas + '.zip')
        downloaded_Reg_atlas_file = os.path.join(org_atlases_version_folder, 'ORG-RegAtlas-100HCP.zip')
    
        downloaded_population_T1_file = os.path.join(org_atlases_version_folder, '100HCP-population-mean-T1.nii.gz')
        downloaded_population_T2_file = os.path.join(org_atlases_version_folder, '100HCP-population-mean-T2.nii.gz')
        downloaded_population_b0_file = os.path.join(org_atlases_version_folder, '100HCP-population-mean-b0.nii.gz')
    else:
        print('<wm_download_org_atlas> ' + requested_atlas + 'is not available. Please check input.')
        print('')
        exit()
    
    print("")
    print("===== <wm_download_org_atlas.py> ")
    print("===== Release version   : ", version)
    print("===== Download from     : ", repo)
    print("===== Output directory  : ", outdir)
    print("===== Requested atlas   : ", requested_atlas)
    print("")
    
    if requested_atlas == 'ORG-800FC-100HCP':
        print('<wm_download_org_atlas> The ' +requested_atlas+ ' atlas is an anatomically curated white matter atlas, with annotated anatomical label provided for each cluster. '\
               + 'MRML files that are used to organize fiber clusters that belong to the same anatomical fiber tract are provided in the atlas.') 
        print('')
    elif requested_atlas == 'ORG-2000FC-100HCP':
        print('<wm_download_org_atlas> The ' +requested_atlas+ ' atlas is provided for applications that can be benefit from a fine scale white matter parcellation. This an uncurated white matter atlas.')
        print('')
    
    if requested_atlas == 'ORG-800FC-100HCP' or requested_atlas == 'ORG-2000FC-100HCP':
        print('<wm_download_org_atlas> Tractography registration atlas.')
        if not os.path.exists(os.path.join(org_atlases_version_folder, 'ORG-RegAtlas-100HCP/registration_atlas.vtk')):
            download_file(REG_atlas_url, downloaded_Reg_atlas_file)
            extract_from_zip(downloaded_Reg_atlas_file, org_atlases_version_folder, remove_after_extraction=True)
        else:
            print('<wm_download_org_atlas> Skip downloading: There is an existing registration atlas at \'ORG-RegAtlas-100HCP/registration_atlas.vtk\' in the output folder.')
            print('')
    
        print('<wm_download_org_atlas> Fiber clustering atlas.')
        output_atlas_folder = os.path.join(org_atlases_version_folder, requested_atlas)
        if not os.path.exists(output_atlas_folder):
            download_file(FC_atlas_url, downloaded_FC_atlas_file)
            extract_from_zip(downloaded_FC_atlas_file, org_atlases_version_folder, remove_after_extraction=True)
        else:
            print('<wm_download_org_atlas> Skip downloading: There is an existing fiber clustering atlas at \''+requested_atlas+'\' in the output folder.')
            print('')
    
        print('<wm_download_org_atlas> Population mean T1/T2/b0 images.')
        if not os.path.exists(downloaded_population_T1_file):
            download_file(population_T1_url, downloaded_population_T1_file)
        else:
            print('<wm_download_org_atlas> Skip downloading: There is an existing \'100HCP-population-mean-T1.nii.gz\' in the output folder.', end=' ')
            print('')
        if not os.path.exists(downloaded_population_T2_file):
            download_file(population_T2_url, downloaded_population_T2_file)
        else:
            print('<wm_download_org_atlas> Skip downloading: There is an existing \'100HCP-population-mean-T2.nii.gz\' in the output folder.', end=' ')
            print('')
        if not os.path.exists(downloaded_population_b0_file):
            download_file(population_b0_url, downloaded_population_b0_file)
        else:
            print('<wm_download_org_atlas> Skip downloading: There is an existing \'100HCP-population-mean-b0.nii.gz\' in the output folder.', end=' ')
            print('')
    
    print('')
    print('<wm_download_org_atlas> Successfully downloaded to', outdir)

if __name__ == '__main__':
    main()
