# -*- coding: utf-8 -*-

import enum
import pathlib

fname_sep = "."


class SlicerFileExtension(enum.Enum):
    MRML = "mrml"


class VTKFileExtension(enum.Enum):
    VTP = "vtp"


def build_suffix(extension):

    return fname_sep + extension.value


def find_filenames(path, extension, stem=False):
    filenames = sorted(list(pathlib.Path(path).glob("*" + build_suffix(extension))))
    if stem:
        return [elem.stem for elem in filenames]
    else:
        return filenames


def save2txt(fname, sequence):

    with open(fname, "w") as f:
        f.write("\n".join(sequence))
