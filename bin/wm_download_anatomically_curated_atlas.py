#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import enum
import os
import ssl
import sys
import urllib.error
import urllib.parse
import urllib.request
import zipfile

ZENODO_RECORD_ROOT_URL = "https://zenodo.org/record/"
ZENODO_API_RECORD_ROOT_URL = "https://zenodo.org/api/records/"

fname_sep = "."


class DataExchangeFormatFileExtension(enum.Enum):
    JSON = "json"


class ORGAtlasVersion(enum.Enum):
    V1_1 = ("v1.1", str(2648284))
    V1_1_1 = ("v1.1.1", str(2648292))
    V1_2 = ("v1.2", str(4156927))
    V1_3_A = ("v1.3a", str(5109770))
    V1_3_B = ("v1.3b", str(7784967))
    V1_4 = ("v1.4", str(8082481))

    @staticmethod
    def get_version(name):
        return ORGAtlasVersion(name).value[0]

    @staticmethod
    def get_record(name):
        return ORGAtlasVersion(name).value[1]


def build_download_url(root_url, version):
    return urllib.parse.urljoin(root_url, ORGAtlasVersion.get_record(ORGAtlasVersion.__getitem__(version)))


def build_suffix(extension):
    return fname_sep + extension.value


def _build_arg_parser():

    parser = argparse.ArgumentParser(
        description="Download a pre-provided anatomically curated fiber clustering white matter atlas.",
        epilog="Written by Fan Zhang, fzhang@bwh.harvard.edu")
    
    parser.add_argument(
        'outputDirectory',
        help='The output directory, where the atlas will be downloaded, should be a new empty directory. It will be created if needed.')
    parser.add_argument(
        '-atlas', action="store", dest="requested_atlas", type=str, 
        help='Name of the atlas. Currently, \'ORG-800FC-100HCP\' and \'ORG-2000FC-100HCP\' are available to download.')
    parser.add_argument(
        '--version', choices=ORGAtlasVersion._member_names_, default=ORGAtlasVersion.V1_2.name, help="Atlas version.")

    return parser


def _parse_args(parser):

    return parser.parse_args()


def main():

    parser = _build_arg_parser()
    args = _parse_args(parser)

    if (not os.environ.get('PYTHONHTTPSVERIFY', '') and getattr(ssl, '_create_unverified_context', None)):
        ssl._create_default_https_context = ssl._create_unverified_context

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
                status = f"{file_size_dl} bytes [{p:.2%}]"
                status = status + chr(8)*(len(status)+1)
                sys.stdout.write(status)
    
            output.close()  
            print('')      
        except:
            print('[ERROR] Download fails! Please check your network connection and rerun the script.')
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
        print(f"<{os.path.basename(__file__)}> Output directory {args.outputDirectory} does not exist, creating it.")
        os.makedirs(outdir)
    
    requested_atlas = args.requested_atlas

    repo = build_download_url(ZENODO_RECORD_ROOT_URL, args.version)
    repo = f"{repo}/"
    version = ORGAtlasVersion.get_version(ORGAtlasVersion.__getitem__(args.version))

    metadata_url = build_download_url(ZENODO_API_RECORD_ROOT_URL, args.version)
    metadata_local_file_rootname = ORGAtlasVersion.get_record(ORGAtlasVersion.__getitem__(args.version))
    metadata_local_file_basename = metadata_local_file_rootname + build_suffix(
        DataExchangeFormatFileExtension.JSON)

    org_atlases_version_folder_name = f'ORG-Atlases-{version[1:]}'
    org_atlases_version_folder = os.path.join(outdir, org_atlases_version_folder_name)
    if not os.path.exists(org_atlases_version_folder):
        os.makedirs(org_atlases_version_folder)

    metadata_local_fname = os.path.join(
        org_atlases_version_folder,
        metadata_local_file_basename,
    )

    if requested_atlas == 'ORG-800FC-100HCP' or requested_atlas == 'ORG-2000FC-100HCP':
        FC_atlas_url = f'{repo}files/{requested_atlas}.zip?download=1'
        REG_atlas_url = f'{repo}files/ORG-RegAtlas-100HCP.zip?download=1'
    
        population_T1_url = f'{repo}files/100HCP-population-mean-T1.nii.gz?download=1'
        population_T2_url = f'{repo}files/100HCP-population-mean-T2.nii.gz?download=1'
        population_b0_url = f'{repo}files/100HCP-population-mean-b0.nii.gz?download=1'
        
        downloaded_FC_atlas_file = os.path.join(org_atlases_version_folder, f'{requested_atlas}.zip')
        downloaded_Reg_atlas_file = os.path.join(org_atlases_version_folder, 'ORG-RegAtlas-100HCP.zip')
    
        downloaded_population_T1_file = os.path.join(org_atlases_version_folder, '100HCP-population-mean-T1.nii.gz')
        downloaded_population_T2_file = os.path.join(org_atlases_version_folder, '100HCP-population-mean-T2.nii.gz')
        downloaded_population_b0_file = os.path.join(org_atlases_version_folder, '100HCP-population-mean-b0.nii.gz')
    else:
        print(f'<{os.path.basename(__file__)}> {requested_atlas} is not available. Please check input.')
        print('')
        exit()
    
    print("")
    print("===== <wm_download_org_atlas.py> ")
    print(f"===== Release version   : {version}")
    print(f"===== Download from     : {repo}")
    print(f"===== Output directory  : {outdir}")
    print(f"===== Requested atlas   : {requested_atlas}")
    print("")
    
    if requested_atlas == 'ORG-800FC-100HCP':
        print(f'<{os.path.basename(__file__)}> The {requested_atlas} atlas is an anatomically curated white matter atlas, with annotated anatomical label provided for each cluster. '\
               + 'MRML files that are used to organize fiber clusters that belong to the same anatomical fiber tract are provided in the atlas.') 
        print('')
    elif requested_atlas == 'ORG-2000FC-100HCP':
        print(f'<{os.path.basename(__file__)}> The {requested_atlas} atlas is provided for applications that can be benefit from a fine scale white matter parcellation. This an uncurated white matter atlas.')
        print('')
    
    if requested_atlas == 'ORG-800FC-100HCP' or requested_atlas == 'ORG-2000FC-100HCP':
        print(f'<{os.path.basename(__file__)}> Tractography registration atlas.')
        if not os.path.exists(os.path.join(org_atlases_version_folder, 'ORG-RegAtlas-100HCP/registration_atlas.vtk')):
            download_file(REG_atlas_url, downloaded_Reg_atlas_file)
            extract_from_zip(downloaded_Reg_atlas_file, org_atlases_version_folder, remove_after_extraction=True)
        else:
            print(f'<{os.path.basename(__file__)}> Skip downloading: There is an existing registration atlas at \'ORG-RegAtlas-100HCP/registration_atlas.vtk\' in the output folder.')
            print('')
    
        print(f'<{os.path.basename(__file__)}> Fiber clustering atlas.')
        output_atlas_folder = os.path.join(org_atlases_version_folder, requested_atlas)
        if not os.path.exists(output_atlas_folder):
            download_file(FC_atlas_url, downloaded_FC_atlas_file)
            extract_from_zip(downloaded_FC_atlas_file, org_atlases_version_folder, remove_after_extraction=True)
        else:
            print(f'<{os.path.basename(__file__)}> Skip downloading: There is an existing fiber clustering atlas at \'{requested_atlas}\' in the output folder.')
            print('')
    
        print(f'<{os.path.basename(__file__)}> Population mean T1/T2/b0 images.')
        if not os.path.exists(downloaded_population_T1_file):
            download_file(population_T1_url, downloaded_population_T1_file)
        else:
            print(f'<{os.path.basename(__file__)}> Skip downloading: There is an existing \'100HCP-population-mean-T1.nii.gz\' in the output folder.', end=' ')
            print('')
        if not os.path.exists(downloaded_population_T2_file):
            download_file(population_T2_url, downloaded_population_T2_file)
        else:
            print(f'<{os.path.basename(__file__)}> Skip downloading: There is an existing \'100HCP-population-mean-T2.nii.gz\' in the output folder.', end=' ')
            print('')
        if not os.path.exists(downloaded_population_b0_file):
            download_file(population_b0_url, downloaded_population_b0_file)
        else:
            print(f'<{os.path.basename(__file__)}> Skip downloading: There is an existing \'100HCP-population-mean-b0.nii.gz\' in the output folder.', end=' ')
            print('')
        if not os.path.exists(metadata_local_fname):
            download_file(metadata_url, metadata_local_fname)
        else:
            print(f"<wm_download_org_atlas> Skip downloading: There is an existing \'{metadata_local_file_basename}\' in the output folder.", end=" ")
            print("")

    print('')
    print(f'<{os.path.basename(__file__)}> Successfully downloaded to {outdir}')

if __name__ == '__main__':
    main()
