""" congeal_multisubject.py

re-implementation of fiber tractography registration (group)

class MultiSubjectRegistration


"""

import os
import time
import numpy
import vtk
from joblib import Parallel, delayed

HAVE_PLT = 1
try:
    import matplotlib
    # Force matplotlib to not use any Xwindows backend.
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except:
    print "<wm_congeal_multisubject.py> Error importing matplotlib.pyplot package, can't plot objectives.\n"
    HAVE_PLT = 0

import whitematteranalysis as wma

class MultiSubjectRegistration:

    def __init__(self):
        # parameters that can be set by user
        self.output_directory = "."
        self.input_directory = None
        self.random_seed = None
        
        # performance options set by user
        self.verbose = False
        self.parallel_jobs = 4
        self.parallel_verbose = 0
        #self.points_per_fiber = 5
        #self.points_per_fiber = 9
        self.points_per_fiber = 15
        self.render = True
        self.nonlinear = False
        
        # optimizer parameters set by user
        self.maxfun = 300
        self.sigma = 20
        self.initial_step = 10
        self.final_step = 5
        self.mean_brain_size = 4000
        # A different sample of 1000 is taken every iteration, so we really
        # see more fibers than this over the course of the registration
        self.subject_brain_size = 1000
        
        # internal stuff
        self.polydatas = list()
        self.transforms = list()
        self.transforms_as_array = list()
        self.objectives_before = list()
        self.objectives_after = list()
        self.total_iterations = 0
        self.subject_ids = list()

        self.target_landmarks = list()
        self.nonlinear_grid_resolution = 3
        self.nonlinear_grid_3 = [-120, 0, 120]
        self.nonlinear_grid_4 = [-120, -60, 60, 120]
        self.nonlinear_grid_5 = [-120, -60, 0, 60, 120]
        self.nonlinear_grid_6 = [-120, -60, -20, 20, 60, 120]
        self.nonlinear_grid_8 = [-120, -85, -51, -17, 17, 51, 85, 120]
        self.nonlinear_grid_10 = [-120, -91, -65, -39, -13, 13, 39, 65, 91, 120]
        # random order so that optimizer does not start in one corner every time
        # the order was computed as
        # numpy.random.permutation(range(0,27))
        # numpy.random.permutation(range(0,64))
        # numpy.random.permutation(range(0,125))
        # numpy.random.permutation(range(0,216))
        # numpy.random.permutation(range(0,512))
        # numpy.random.permutation(range(0,1000))
        self.grid_order_3 = [22,  8,  0, 14,  5, 19, 20,  9, 17, 15,  2, 11,  1, 12, 21,  6, 25, 7,  3, 24, 13, 18, 26, 16, 23,  4, 10]
        self.grid_order_4 = [58, 54,  4, 62, 51, 41, 52, 45, 59,  8, 29, 17, 61, 46, 18, 22, 34, 42, 21,  0,  3, 39, 27, 13, 60, 12,  2, 15,  5, 28,  7, 43, 31, 38, 33, 55,  1, 37, 47, 30, 24, 35, 14, 50, 20, 36, 44, 53, 16, 57, 56, 10, 48, 26,  6, 25,  9, 32, 19, 11, 63, 23, 49, 40]
        self.grid_order_5 = [ 75,  18,  54,  61,  64,  73,  95,  13, 111, 118,  43,   7,  46, 56,   4, 124,  77,  98,  72,  60,  38,  80,  36,  27, 120, 119, 51,  81,   0,  93,  11,  41,  69,  83, 107,  12, 106,  30,  53, 105,  33,  91,  28,  17,  58,  90,  45,  94,  14,  26,  84,   1, 92,  21,  47,  59, 100,   2,   3,  87,  65, 102,  68,  20,  85, 79,  82,  15,  32,  88, 115,  74,   6,  19,  35,  99, 104, 109, 70, 101,  96,  66,  52,  48,  49,  31,  97, 122,  78, 113,  55, 112,  76,  44,  23, 103,  16,  10, 123,  86,  39,   8,  62, 110, 42, 114,  40, 117,  63,   9,  25,  67,  71,  37,  24, 116,  57, 89, 121,  34,   5,  29, 108,  50,  22]
        self.grid_order_6 = [165,  63, 129, 170, 148, 131,   1, 205, 181,  69,  35, 100,   6, 13,  82, 193, 159, 130,  54, 164, 174,  62, 108, 101, 169,  93, 112,  42, 110, 213, 107,  47,  45, 140, 138, 199, 125, 117,   8, 41, 215,  74, 143, 155, 144, 187,  26,  60, 139,  58,  97,  10, 113,  34, 116,  33, 202, 210,  27, 135, 152, 206, 128,  16, 177, 25,  67, 192, 147, 132, 160, 158, 161,  90, 102,  49,  21, 191, 32,  18,  81, 157, 114, 175,  94, 172, 207, 186, 167, 163, 196, 118,  28,  43, 133, 171, 211,  77,  56, 195, 173,  57,  96,  29, 64, 180,  89, 190, 115,  20,  52,  50,   4, 141,  98, 134, 109, 149, 176, 212,  11,   0, 146,  65,  91,  23,  53,  44, 123,  87, 24, 178, 184,  68, 124,  46,  76, 151, 127, 204, 154, 150, 106, 70,  37,  84,  17,  12, 189,   2,  92,  36,  71,  39,  30,  75, 179, 168,  73, 121,  86, 214, 188,  59, 209,  22,  19, 153, 162, 99, 182,  14,  48, 119, 203,  66,  61, 103, 208, 145,  79,  85, 142,  72, 126, 194, 104, 122, 198, 120, 200, 183, 201,   3,  78, 40,  83, 137,  31, 111,  15,  51,   9, 185,  55,  38, 156, 136, 7,  95,  80, 105, 166,  88, 197,   5]
        self.grid_order_8 = [455, 240,  90,  94, 412,  27, 287, 423, 102,   2, 171, 220, 159, 56,   6, 233, 454,  33,  73, 180, 213, 205, 253, 294, 105, 202, 224, 474, 299, 460, 282, 432, 468, 121, 307, 112, 409, 420,  76,  7, 239, 217, 279, 103, 419,  58, 425, 152, 387, 111, 391, 293, 179, 471,  20, 310, 264,  63, 100,  65, 346, 510,  82, 257, 302, 47, 242, 490, 141, 232, 211, 356, 110, 284,  86, 336, 396, 352, 451,  30, 497, 192, 427, 361, 144, 417, 466, 331, 416, 193, 338, 128, 163, 101,   4, 116, 208, 458, 266, 183, 198, 169, 487, 295, 218, 340, 298,  66, 123, 348, 447, 433, 383, 379, 268, 214, 456, 439, 311, 139, 230, 498, 484, 478, 262, 453, 251, 229, 207, 108, 83, 403, 470, 443, 271,  88,  32, 367, 162, 469, 308, 407,  51, 485, 377, 489, 496, 113, 345, 448, 504, 135, 323, 368, 149, 131, 488, 155, 385, 511, 365, 156, 166, 176, 260, 473, 339,   5, 283, 280,  61,  95,  28, 349,  72, 509, 446, 315, 191,  37,  34,  21, 197, 415, 364, 476, 341, 358, 147, 153, 125, 170, 181, 194, 436, 394, 390,  54,  81, 237,  15, 160, 269,  62, 386, 209, 275, 273, 64, 499, 120, 457,  75, 330, 384, 317, 501, 281, 254, 145, 444, 203, 303, 475, 133, 347, 472, 223, 143, 137, 395,  18,  29, 236, 459,  87, 329, 129, 508, 332, 380, 161,  19,  96,  24,  80,  23, 91,  69, 263, 369, 371,  77, 228, 250,  46,  70,  67,  10, 327, 463, 127,  11,  42, 119, 334, 406,  35,  45, 235,  41, 174, 274, 238, 344, 316, 414, 265, 404, 418, 486, 115, 114, 355, 505, 304, 357, 256, 221, 219, 450, 267, 150, 408, 292,  26, 297, 168, 118, 319, 177, 337,   8, 204, 405, 309, 225,  12, 479, 397, 278,   9, 276,  31, 410, 461, 399, 434, 321, 477, 503, 370,  99,  89, 188, 78, 438, 353,  71, 167, 320, 172, 398, 360, 382, 291, 245, 378, 246, 138,  17,  74, 124, 252,   0, 422, 244, 157,  50,  14, 350, 393, 306, 122, 190, 288, 402, 142, 222,  16,  92, 206, 313, 411, 216, 342, 261, 400, 389,  84, 117, 502, 483, 189, 290, 431,  43, 201, 227, 210, 430, 148,  68, 241, 363, 493, 182, 366,  57, 435, 465, 212, 199,  60, 146,  22, 226, 421, 507, 376,  53, 354, 130, 324, 184, 318,  49,  93, 333,  39, 258, 247, 440, 413,  59, 424, 165, 270, 429, 445, 464, 480, 401,  38, 442, 322, 289, 164, 231, 97, 328, 500, 175, 375, 301, 255, 286, 305,  44, 215, 428, 343, 52, 151, 449, 506, 388, 136, 248,  98, 132, 134, 234, 326, 374, 158, 243, 249,  48, 173, 200, 452, 482,  79, 140, 104, 312, 272, 359, 437, 373,  13, 381,  55, 196, 351, 106, 296,  40, 300, 392, 186,   3, 495, 126, 325, 154, 462, 362, 259, 481, 491, 441, 467, 335, 372, 109, 494,   1, 285,  25, 314,  36, 185, 178, 426, 492, 107, 187, 195, 277,  85]

        self.grid_order_10 = [591, 415, 411, 905, 539, 368, 714, 506, 488, 228, 856, 482, 858, 407, 403, 323, 829, 626, 371, 77, 232, 653, 891, 92, 370, 813, 175, 427, 408, 915, 995, 521, 789, 433, 67, 528, 325, 588, 544, 516, 652, 160, 782, 589, 899, 412, 621, 203, 584, 581, 938, 737, 888, 47, 234, 393, 658, 820, 93, 372, 781, 775, 937, 139, 580, 85, 985, 803, 960, 980, 237, 784, 733, 917, 128, 462, 447, 876, 461, 753, 641, 341, 843, 907, 110, 105, 600, 265, 23, 156, 363, 367, 841, 297, 874, 419, 951, 41, 762, 165, 682, 527, 704, 627, 613, 149, 622, 252, 609, 898, 113, 394, 817, 281, 900, 366, 178, 881, 213, 642, 46, 740, 647, 216, 935, 29, 698, 22, 352, 19, 127, 962, 478, 610, 759, 707, 319, 457, 383, 897, 357, 805, 748, 134, 320, 918, 345, 109, 410, 715, 335, 353, 4, 513, 45, 657, 585, 853, 114, 824, 177, 822, 132, 422, 689, 438, 131, 844, 167, 181, 391, 697, 973, 454, 258, 361, 570, 97, 768, 28, 864, 212, 37, 179, 596, 162, 280, 878, 170, 184, 708, 2, 651, 683, 640, 979, 307, 52, 747, 517, 862, 205, 618, 529, 840, 531, 913, 202, 764, 518, 147, 78, 699, 469, 423, 857, 34, 709, 552, 927, 389, 358, 992, 562, 649, 1, 563, 453, 887, 314, 800, 43, 629, 68, 837, 637, 574, 866, 72, 250, 390, 794, 569, 150, 896, 920, 299, 117, 44, 432, 716, 176, 26, 282, 532, 399, 788, 164, 692, 9, 676, 186, 694, 695, 833, 210, 48, 15, 861, 129, 838, 91, 263, 991, 661, 934, 835, 233, 424, 998, 555, 473, 157, 943, 309, 182, 603, 452, 311, 910, 810, 266, 169, 444, 241, 327, 956, 188, 257, 807, 140, 542, 359, 686, 20, 793, 769, 804, 523, 215, 583, 787, 174, 942, 656, 73, 294, 706, 582, 111, 577, 460, 606, 463, 121, 329, 116, 253, 35, 751, 732, 392, 154, 537, 719, 287, 931, 564, 339, 290, 828, 534, 474, 437, 662, 540, 826, 251, 957, 193, 235, 440, 104, 288, 435, 55, 272, 961, 166, 796, 249, 912, 405, 152, 53, 725, 701, 381, 612, 530, 119, 211, 993, 590, 648, 285, 718, 655, 260, 214, 514, 936, 439, 286, 509, 133, 221, 308, 275, 559, 362, 443, 825, 550, 868, 666, 254, 291, 994, 779, 680, 413, 705, 143, 62, 839, 21, 515, 3, 155, 735, 384, 70, 755, 968, 332, 123, 120, 986, 373, 922, 355, 561, 284, 501, 632, 338, 304, 445, 883, 467, 790, 665, 982, 579, 799, 63, 360, 734, 639, 758, 450, 490, 118, 197, 939, 832, 500, 894, 449, 988, 711, 587, 244, 378, 103, 909, 557, 667, 659, 965, 135, 465, 192, 558, 535, 889, 38, 448, 400, 746, 811, 446, 953, 673, 66, 356, 279, 996, 455, 895, 125, 631, 739, 239, 599, 821, 397, 636, 766, 18, 736, 82, 549, 99, 401, 61, 690, 436, 638, 989, 108, 16, 713, 472, 245, 507, 687, 972, 185, 240, 566, 710, 141, 785, 90, 209, 451, 495, 404, 903, 315, 42, 981, 79, 347, 729, 330, 479, 39, 480, 302, 295, 884, 669, 545, 893, 761, 717, 242, 486, 49, 846, 382, 946, 144, 264, 977, 771, 767, 406, 777, 831, 269, 5, 904, 434, 57, 496, 869, 497, 693, 208, 932, 56, 101, 928, 873, 122, 100, 354, 560, 145, 337, 475, 595, 721, 8, 420, 538, 818, 575, 752, 191, 852, 255, 678, 880, 187, 670, 301, 567, 236, 380, 774, 634, 770, 547, 270, 396, 303, 124, 83, 646, 802, 624, 74, 130, 966, 819, 750, 218, 336, 441, 508, 671, 385, 617, 724, 426, 568, 0, 503, 978, 823, 773, 340, 51, 684, 225, 967, 344, 136, 541, 7, 376, 668, 493, 386, 402, 96, 262, 276, 84, 741, 76, 33, 318, 231, 801, 431, 224, 797, 834, 172, 700, 875, 949, 619, 351, 305, 24, 278, 342, 916, 200, 941, 194, 468, 554, 879, 754, 757, 466, 812, 425, 576, 923, 925, 98, 416, 10, 911, 688, 827, 95, 908, 565, 40, 744, 969, 106, 954, 553, 702, 664, 151, 102, 201, 863, 81, 792, 615, 749, 87, 248, 886, 115, 964, 388, 50, 267, 890, 316, 892, 206, 958, 350, 578, 227, 476, 13, 247, 872, 107, 608, 25, 926, 456, 974,  6, 30, 870, 94, 499, 760, 919, 983, 12, 720, 198, 246, 65, 644, 489, 959, 271, 594, 546, 786, 798, 510, 333, 519, 171, 663, 86, 328, 229, 851, 89, 847, 17, 283, 743, 620, 421, 395, 60, 322, 442, 313, 484, 616, 409, 614, 975, 158, 369, 997, 293, 226, 756, 195, 58, 32, 296, 524, 536, 703, 940, 630, 650, 836, 153, 806, 999, 551, 428, 348, 783, 54, 633, 458, 161, 660, 877, 491, 944, 321, 865, 947, 906, 173, 256, 830, 742, 765, 848, 685, 791, 712, 990, 871, 728, 429, 625, 414, 492, 924, 326, 730, 238, 223, 142, 138, 289, 597, 349, 64, 901, 814, 572, 525, 220, 526, 727, 675, 207, 745, 11, 180, 723, 204, 921, 645, 27, 365, 882, 398, 722, 643, 573, 159, 36, 306, 259, 628, 31, 504, 855, 533, 607, 485, 520, 459, 334, 292, 933, 976, 623, 815, 312, 331, 854, 190, 88, 375, 778, 14, 487, 317, 592, 930, 952, 963, 859, 418, 168, 343, 885, 955, 71, 548, 571, 948, 346, 377, 763, 502, 795, 498, 481, 691, 842, 945, 902, 816, 277, 146, 196, 189, 681, 374, 126, 845, 364, 217, 324, 222, 971, 867, 298, 183, 635, 430, 59, 512, 984, 471, 987, 808, 850, 75, 556, 199, 593, 274, 605, 80, 849, 950, 914, 929, 477, 112, 300, 602, 230, 163, 522, 780, 674, 69, 268, 726, 679, 860, 586, 494, 672, 604, 505, 776, 611, 970, 511, 738, 731, 417, 483, 273, 601, 772, 387, 543, 696, 470, 464, 809, 310, 261, 654, 598, 677, 137, 148, 243, 379, 219]

        self.initialize_nonlinear_grid()

    def initialize_nonlinear_grid(self):
        """This initializes the nonlinear grid. This should be called before adding

        subjects' polydata. Calling it in init handles this.
        """
        
        self.target_landmarks = list()
        if self.nonlinear_grid_resolution == 3:
            grid = self.nonlinear_grid_3
            grid_order = self.grid_order_3
        elif self.nonlinear_grid_resolution == 4:
            grid = self.nonlinear_grid_4
            grid_order = self.grid_order_4
        elif self.nonlinear_grid_resolution == 5:
            grid = self.nonlinear_grid_5
            grid_order = self.grid_order_5                        
        elif self.nonlinear_grid_resolution == 6:
            grid = self.nonlinear_grid_6
            grid_order = self.grid_order_6
        elif self.nonlinear_grid_resolution == 8:
            grid = self.nonlinear_grid_8
            grid_order = self.grid_order_8
        elif self.nonlinear_grid_resolution == 10:
            grid = self.nonlinear_grid_10
            grid_order = self.grid_order_10
        else:
            print "<congeal_multisubject.py> Error: Unknown nonlinear grid mode:", self.nonlinear_grid_resolution
        tmp = list()
        for r in grid:
            for a in grid:
                for s in grid:
                    tmp.append([r, a, s])
        # now shuffle the order of these points to avoid biases
        for idx in grid_order:
            self.target_landmarks.extend(tmp[idx])
        self.target_points = wma.register_two_subjects_nonlinear.convert_numpy_array_to_vtk_points(self.target_landmarks)

    def update_nonlinear_grid(self):
        """This updates the nonlinear grid. Subjects must be added first using

        add_polydata.
        """
        
        # Compute the new grid
        new_target_landmarks = list()
        if self.nonlinear_grid_resolution == 3:
            grid = self.nonlinear_grid_3
            grid_order = self.grid_order_3
        elif self.nonlinear_grid_resolution == 4:
            grid = self.nonlinear_grid_4
            grid_order = self.grid_order_4
        elif self.nonlinear_grid_resolution == 5:
            grid = self.nonlinear_grid_5
            grid_order = self.grid_order_5                        
        elif self.nonlinear_grid_resolution == 6:
            grid = self.nonlinear_grid_6
            grid_order = self.grid_order_6
        elif self.nonlinear_grid_resolution == 8:
            grid = self.nonlinear_grid_8
            grid_order = self.grid_order_8
        elif self.nonlinear_grid_resolution == 10:
            grid = self.nonlinear_grid_10
            grid_order = self.grid_order_10
        else:
            print "<congeal_multisubject.py> Error: Unknown nonlinear grid mode:", self.nonlinear_grid_resolution
        tmp = list()
        for r in grid:
            for a in grid:
                for s in grid:
                    tmp.append([r, a, s])
        # now shuffle the order of these points to avoid biases
        for idx in grid_order:
            new_target_landmarks.append(tmp[idx])

        # Apply the existing transform to the target landmarks to compute new source landmarks
        new_transforms = list()
        for trans in self.transforms_as_array:
                source_landmarks = wma.register_two_subjects_nonlinear.convert_numpy_array_to_vtk_points(trans)
                tps =  wma.register_two_subjects_nonlinear.compute_thin_plate_spline_transform(source_landmarks,self.target_points)
                tps.Inverse()
                new_source_landmarks = list()
                for pt in new_target_landmarks:
                        pt2 = tps.TransformPoint(pt[0], pt[1], pt[2])
                        new_source_landmarks.append(pt2)
                new_transforms.append(numpy.array(new_source_landmarks).flatten())
        new_target_landmarks = numpy.array(new_target_landmarks).flatten()
        # Update all the relevant variables (the spline transform does not change but all source and target points do)
        print "UPDATE NONLINEAR GRID: ", len(self.target_landmarks), len(trans), "==>", len(new_target_landmarks), len(new_transforms[-1]),
        self.transforms_as_array = new_transforms
        self.target_landmarks = new_target_landmarks
        self.target_points = wma.register_two_subjects_nonlinear.convert_numpy_array_to_vtk_points(self.target_landmarks)

    def add_polydata(self, polydata, subject_id):
        """Add a subject's data to the groupwise registration. self.nonlinear

        must be set before calling this, if nonlinear registration is desired.
        """
        
        self.polydatas.append(polydata)
        if self.nonlinear:
            # This sets up identity transform to initialize. This will
            # be re-calculated with current grid resolution in
            # update_nonlinear_grid.
            # Set source and target points equal for initial identity transform:
            trans = wma.register_two_subjects_nonlinear.compute_thin_plate_spline_transform(self.target_points,self.target_points)
            self.transforms.append(trans)
            self.transforms_as_array.append(self.target_landmarks)
        else:
            trans = vtk.vtkTransform()
            self.transforms.append(trans)
            self.transforms_as_array.append(numpy.array([0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0]).astype(float))
        self.subject_ids.append(subject_id)
        
    def remove_mean_from_transforms(self):
        """ Remove mean rotations and mean scaling and mean
         translations from transforms.  The mean rotation and
         translation should not affect the objective function anyway.
         A mean scaling will affect the objective, i.e. if the brain
         shrinks all distances become smaller and the similarity is
         higher. This is not the desired effect."""

        if self.nonlinear:
            # remove any average displacement of each source point.
            # this means the mean of the source points must equal the target point
            transforms_array = numpy.array(self.transforms_as_array)
            meansource = numpy.mean(transforms_array, 0)
            landmarks = numpy.array(self.target_landmarks)
            meandisp = meansource - landmarks
            if self.verbose:
                print "MEAN DISPLACEMENT:", meandisp

            if self.verbose:
                print "<congeal.py> TRANSFORMS"
                print numpy.round(transforms_array * 100) / 100
                print "<congeal.py> Removing current (accumulated) mean transform before computing objective:"
                print meandisp

            for transform in self.transforms_as_array:
                transform[:] = transform - meandisp

            transforms_array = numpy.array(self.transforms_as_array)
            meansource = numpy.mean(transforms_array, 0)
            landmarks = numpy.array(self.target_landmarks)
            meandisp = meansource - landmarks
            if self.verbose:
                print "MEAN DISPLACEMENT 2:", meandisp

            matrix_average = numpy.zeros((3,4))
            target_landmarks = wma.register_two_subjects_nonlinear.convert_numpy_array_to_vtk_points(self.target_landmarks)
            for trans in self.transforms_as_array:
                affine_part = vtk.vtkLandmarkTransform()
                source_landmarks = wma.register_two_subjects_nonlinear.convert_numpy_array_to_vtk_points(trans)
                affine_part.SetSourceLandmarks(source_landmarks)
                affine_part.SetTargetLandmarks(target_landmarks)
                affine_part.SetModeToAffine()
                affine_part.Update()
                #print affine_part.GetMatrix()
                for row in [0,1,2]:
                    # the fourth row must be 0, 0, 0, 1 for homogeneous transform so don't need to average it
                    for col in [0,1,2,3]:
                        matrix_average[row,col] += affine_part.GetMatrix().GetElement(row,col)
            matrix_average = numpy.divide(matrix_average, len(self.transforms_as_array))

            # remove this average from the transforms
            # want: txform*meantxform^-1 * sourcept' = targetpt
            # new source points:
            # sourcept' = meantxform * sourcept 
            
            matrix = vtk.vtkMatrix4x4()
            for row in [0,1,2]:
                    # the fourth row must be 0, 0, 0, 1 for homogeneous transform so don't need to average it
                    for col in [0,1,2,3]:
                        matrix.SetElement(row,col, matrix_average[row,col])
            #print "MEAN TRANS", matrix
            mean_trans = vtk.vtkMatrixToLinearTransform()
            mean_trans.SetInput(matrix)
            
            # modify all the source points using the tps inverse and the mean transform we will remove
            new_source_pts = list()
            for trans in self.transforms_as_array:
                trans2 = list()
                # loop over all points defining this transform
                for pt in zip(trans[::3], trans[1::3], trans[2::3]):
                        pt2 = mean_trans.TransformPoint(pt[0], pt[1], pt[2])
                        trans2.extend(pt2)
                new_source_pts.append(numpy.array(trans2))

            for (trans, newtrans) in zip(self.transforms_as_array, new_source_pts):
                #print "REMOVE MEAN FROM TXFORMS:", trans, "NEW:", newtrans
                if self.verbose:
                    print "Mean removed. Source points changed by:", numpy.min(trans-newtrans), numpy.max(trans-newtrans)

            self.transforms_as_array = new_source_pts

            # TEST ONLY (to ensure the above was correct)
            target_landmarks = wma.register_two_subjects_nonlinear.convert_numpy_array_to_vtk_points(self.target_landmarks)
            matrix_average = numpy.zeros((3,4))
            for trans in self.transforms_as_array:
                affine_part = vtk.vtkLandmarkTransform()
                source_landmarks = wma.register_two_subjects_nonlinear.convert_numpy_array_to_vtk_points(trans)
                affine_part.SetSourceLandmarks(source_landmarks)
                affine_part.SetTargetLandmarks(target_landmarks)
                affine_part.SetModeToAffine()
                affine_part.Update()
                #print affine_part.GetMatrix()
                for row in [0,1,2]:
                    # the fourth row must be 0, 0, 0, 1 for homogeneous transform so don't need to average it
                    for col in [0,1,2,3]:
                        matrix_average[row,col] += affine_part.GetMatrix().GetElement(row,col)
            matrix_average = numpy.divide(matrix_average, len(self.transforms_as_array))
            # this print tests the affine part averages to identity
            print "Mean transform should be 0 after removing it from the group:", numpy.mean(numpy.abs(matrix_average[0:3,0:3] - numpy.identity(3)))

        else:
            # Here we are in the affine case, which is simpler.
            transforms_array = numpy.array(self.transforms_as_array)
            meantrans = numpy.mean(transforms_array, 0)
            if self.verbose:
                print "<congeal.py> TRANSFORMS"
                print numpy.round(transforms_array * 100) / 100
                print "<congeal.py> Removing current (accumulated) mean transform before computing objective:"
                print numpy.round(meantrans * 1000) / 1000        

            for transform in self.transforms_as_array:
                transform[0:6] = transform[0:6] - meantrans[0:6]
                transform[6:9] = transform[6:9] / meantrans[6:9]
                transform[9:15] = transform[9:15] - meantrans[9:15]

    def iterate(self):
        """ Run a single iteration of optimization, multiprocessing over input subjects."""
        
        self.total_iterations += 1
        start_time = time.time()

        # set up progress information saving if this is the first iteration
        if self.total_iterations == 1:
            self.progress_filename = os.path.join(self.output_directory, 'registration_performance.txt')
            progress_file = open(self.progress_filename, 'w')
            print >> progress_file, 'iteration','\t', 'sigma', '\t', 'nonlinear', '\t', 'subject_brain_fibers', '\t', 'fibers_per_subject_in_mean_brain','\t', 'mean_brain_fibers','\t', 'maxfun','\t', 'grid_resolution_if_nonlinear','\t', 'initial_step','\t', 'final_step','\t', 'objective_before','\t', 'objective_after', '\t', 'objective_change', '\t', 'objective_percent_change', '\t', 'mean_function_calls_per_subject','\t', 'min_function_calls_per_subject','\t', 'max_function_calls_per_subject','\t', 'subjects_hitting_maxfun','\t', 'total_subjects','\t', 'subjects_decreased','\t', 'mean_subject_change', '\t', 'mean_subject_decrease_if_decreased', '\t', 'time'
            progress_file.close()
            
        # make a directory for the current iteration
        if self.nonlinear:
            dirname = "iteration_%05d_sigma_%03d_grid_%03d" % (self.total_iterations, self.sigma, self.nonlinear_grid_resolution)
        else:
            dirname = "iteration_%05d_sigma_%03d" % (self.total_iterations, self.sigma)

        outdir = os.path.join(self.output_directory, dirname)
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        # make a directory for any intermediate rendering by subprocesses
        outdir_render = os.path.join(outdir, 'progress_rendered')
        if not os.path.exists(outdir_render):
            os.makedirs(outdir_render)

        # Calculate how many fibers are needed to sample from each subject to compute the mean brain at the requested size
        fibers_per_subject = self.mean_brain_size / (len(self.polydatas) - 1)
        if self.verbose:
            print "Fibers per subject for computing mean brain:", fibers_per_subject, "=", self.mean_brain_size, "/",  len(self.polydatas) -1

        # Set up lists of data to pass to the per-subject processes
        mean_list = list()
        subject_list = list()
        mode_list = list()
        sigma_list = list()
        subj_idx_list = list()
        iteration_list = list()
        outdir_list = list()
        stepsize_list = list()
        maxfun_list = list()
        render_list = list()
        grid_resolution_list = list()
        
        # Each subject will be registered to the current model or "mean brain"
        # Sample fibers from each subject for use in the "mean brain"
        subject_sampled_fibers = list()
        for (input_pd, trans) in zip(self.polydatas, self.transforms):
            pd = wma.filter.downsample(input_pd, fibers_per_subject, verbose=False, random_seed=self.random_seed)
            # apply the current transform to this polydata for computation of mean brain
            transformer = vtk.vtkTransformPolyDataFilter()
            if (vtk.vtkVersion().GetVTKMajorVersion() >= 6.0):
                transformer.SetInputData(pd)
            else:
                transformer.SetInput(pd)
            transformer.SetTransform(trans)
            transformer.Update()
            subject_sampled_fibers.append(transformer.GetOutput())
            del transformer
            
        # Loop over all subjects and prepare lists of inputs for subprocesses
        subj_idx = 0
        for input_pd in self.polydatas:
            # Compute current atlas model "mean brain" in a leave-one out fashion.
            # Otherwise, the optimal transform may be identity.
            appender = vtk.vtkAppendPolyData()
            for subj_idx2 in range(len(subject_sampled_fibers)):
                if subj_idx2 != subj_idx:
                    pd = subject_sampled_fibers[subj_idx2]
                    if (vtk.vtkVersion().GetVTKMajorVersion() >= 6.0):
                        appender.AddInputData(pd)
                    else:
                        appender.AddInput(pd)

            # Convert "mean brain" from vtk to numpy format
            appender.Update()
            mean_brain = appender.GetOutput()
            mean_fibers = wma.fibers.FiberArray()
            mean_fibers.convert_from_polydata(mean_brain, self.points_per_fiber)
            mean_fibers = numpy.array([mean_fibers.fiber_array_r,mean_fibers.fiber_array_a,mean_fibers.fiber_array_s])
            #  R,A,S is the first index
            # then fiber number
            # then points along fiber
            mean_list.append(mean_fibers)

            # Now get the current sample of fibers from the subject for registration to the "mean brain"
            pd = wma.filter.downsample(input_pd, self.subject_brain_size, verbose=False, random_seed=self.random_seed)
            fibers = wma.fibers.FiberArray()
            fibers.convert_from_polydata(pd, self.points_per_fiber)
            fibers_array = numpy.array([fibers.fiber_array_r,fibers.fiber_array_a,fibers.fiber_array_s])
            subject_list.append(fibers_array)

            # Append parameter information to lists of parameters for subprocesses
            sigma_list.append(self.sigma)
            if self.nonlinear:
                mode_list.append('Nonlinear')
            else:
                mode_list.append('Linear')                    
            subj_idx_list.append(subj_idx)
            subj_idx += 1
            iteration_list.append(self.total_iterations)
            outdir_list.append(outdir_render)
            stepsize_list.append(numpy.array([self.initial_step, self.final_step]))
            maxfun_list.append(self.maxfun)
            render_list.append(self.render)
            grid_resolution_list.append(self.nonlinear_grid_resolution)
            
        # Multiprocess over subjects
        print "\nITERATION", self.total_iterations, "STARTING MULTIPROCESSING. NUMBER OF JOBS:", self.parallel_jobs, "\n"

        # note we can't pass vtk objects to subprocesses since they can't be pickled.
        ret = Parallel(
            n_jobs=self.parallel_jobs, verbose=self.parallel_verbose)(
                delayed(congeal_multisubject_inner_loop)(fixed, moving, initial_transform, mode, sigma, subject_idx, iteration_count, output_directory, step_size, maxfun, render, grid_resolution)
                for (fixed, moving, initial_transform, mode, sigma, subject_idx, iteration_count, output_directory, step_size, maxfun, render, grid_resolution) in zip(mean_list, subject_list, self.transforms_as_array, mode_list, sigma_list, subj_idx_list, iteration_list, outdir_list, stepsize_list, maxfun_list, render_list, grid_resolution_list))

            
        #print "RETURNED VALUES", ret
        
        # Progress reporting: loop over all registration outputs.
        # Get the current transform for each subject and report the
        # objective values to the user by printing, saving, and plotting.
        self.transforms_as_array = list()
        objective_total_before = 0.0
        objective_total_after = 0.0
        sidx = 0
        functions_per_subject = list()
        objective_changes_per_subject = list()
        decreases = list()
        if HAVE_PLT:
            plt.figure(0)
            plt.title('Iteration '+str(self.total_iterations)+' Objective Values for All Subjects')
            plt.xlabel('objective function computations')
            plt.ylabel('objective value')
            
        for (trans, objectives, diff) in ret:
            self.transforms_as_array.append(trans)
            print "Iteration:", self.total_iterations, "Subject:", sidx, "Objective function computations:", len(objectives), "change", diff
            functions_per_subject.append(len(objectives))
            # Compute total objective for progress reporting.
            objective_total_before += objectives[0]
            if diff < 0:
                objective_total_after += objectives[-1]
                decreases.append(diff)
            else:
                objective_total_after += objectives[0]
            objective_changes_per_subject.append(diff)
            sidx += 1
            if HAVE_PLT:
                plt.figure(0)
                plt.plot(objectives, 'o-', label=sidx)

        number_of_subjects = sidx
        functions_per_subject = numpy.array(functions_per_subject)
        objective_changes_per_subject = numpy.array(objective_changes_per_subject)
        decreases = numpy.array(decreases)

        self.objectives_before.append(objective_total_before)
        self.objectives_after.append(objective_total_after)
        total_change =  self.objectives_after[-1] - self.objectives_before[-1]
        percent_change = total_change / self.objectives_before[-1]
        print "Iteration:", self.total_iterations, "TOTAL objective change:",  total_change
        print "Iteration:", self.total_iterations, "PERCENT objective change:",  percent_change

        if HAVE_PLT:
            plt.figure(0)
            if self.nonlinear:
                fname_fig_base = "iteration_%05d_sigma_%03d_grid_%03d" % (self.total_iterations, self.sigma, self.nonlinear_grid_resolution)
            else:
                fname_fig_base = "iteration_%05d_sigma_%03d_" % (self.total_iterations, self.sigma)
            # Place the legend below the plot so it does not overlap it when there are many subjects
            #lgd = plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), fancybox=False, shadow=False, ncol=1)
            fname_fig = 'objectives_per_subject_' + fname_fig_base + '.pdf'
            # save everything even if the legend is long and goes off the plot
            #plt.savefig(os.path.join(outdir, fname_fig), bbox_extra_artists=(lgd,), bbox_inches='tight')
            plt.savefig(os.path.join(outdir, fname_fig))
            plt.close(0)

        progress_file = open(self.progress_filename, 'a')
        elapsed_time = time.time() - start_time
        if len(decreases) == 0:
            mean_decreases = 0.0
        else:
            mean_decreases = numpy.mean(decreases)
        print >> progress_file, self.total_iterations,'\t', self.sigma, '\t', self.nonlinear, '\t', self.subject_brain_size, '\t', fibers_per_subject,'\t', self.mean_brain_size,'\t', self.maxfun,'\t', self.nonlinear_grid_resolution,'\t', self.initial_step,'\t', self.final_step,'\t', self.objectives_before[-1],'\t', self.objectives_after[-1],'\t', total_change,'\t',  percent_change,'\t', numpy.mean(functions_per_subject), '\t', numpy.min(functions_per_subject), '\t', numpy.max(functions_per_subject), '\t', numpy.sum(functions_per_subject >= self.maxfun), '\t', number_of_subjects,'\t', len(decreases),'\t', numpy.mean(objective_changes_per_subject), '\t', mean_decreases, '\t', elapsed_time
        progress_file.close()

        # remove_mean_from_transforms
        self.remove_mean_from_transforms()

        # update our transforms list for the next iteration
        self.transforms = list()
        for trans in self.transforms_as_array:
            if self.nonlinear:
                vtktrans = wma.register_two_subjects_nonlinear.convert_transform_to_vtk(trans, self.target_points)
                #print vtktrans
            else:
                vtktrans = wma.register_two_subjects.convert_transform_to_vtk(trans)
                #print vtktrans.GetMatrix()
            self.transforms.append(vtktrans)

        # save the current transforms to disk
        wma.registration_functions.write_transforms_to_itk_format(self.transforms, outdir, self.subject_ids)

        
    def save_transformed_polydatas(self, intermediate_save=False, midsag_symmetric=False):
        """ Output polydatas for final or intermediate iterations. """
        
        # this can be slow so notify the user what is happening
        print "\nSAVING TRANSFORMED TRACTOGRAPHY FROM ITERATION", self.total_iterations, "\n"
        
        transform_list = self.transforms
        subject_id_list = self.subject_ids
        if midsag_symmetric:
            transform_list = transform_list[::2]
            subject_id_list = subject_id_list[::2]

        if intermediate_save:
            # Make a directory for the current iteration.
            # This directory name must match the one created above in the iteration.
            if self.nonlinear:
                dirname = "iteration_%05d_sigma_%03d_grid_%03d" % (self.total_iterations, self.sigma, self.nonlinear_grid_resolution)
            else:
                dirname = "iteration_%05d_sigma_%03d" % (self.total_iterations, self.sigma)
            outdir = os.path.join(self.output_directory, dirname)
            if not os.path.exists(outdir):
                os.makedirs(outdir)
            
            # make a directory for output polydatas
            outdir_pds = os.path.join(outdir, 'transformed_output')
            if not os.path.exists(outdir_pds):
                os.makedirs(outdir_pds)
                
            output_pds = wma.registration_functions.transform_polydatas_from_disk(self.input_directory, transform_list, outdir_pds)

        else:
            # make a directory for the final output
            outdir = os.path.join(self.output_directory, 'output_tractography')
            if not os.path.exists(outdir):
                os.makedirs(outdir)

            if self.verbose:
                for trans in  self.transforms:
                    if self.nonlinear:
                        print trans
                    else:
                        print trans.GetMatrix()

            output_pds = wma.registration_functions.transform_polydatas_from_disk(self.input_directory, transform_list, outdir)
        
            # Save the current atlas representation to disk.
            # Right now this is all the input fibers from all subjects.
            output_pds = list()
            for (pd, trans) in zip (self.polydatas, self.transforms):
                transformer = vtk.vtkTransformPolyDataFilter()
                if (vtk.vtkVersion().GetVTKMajorVersion() >= 6.0):
                    transformer.SetInputData(pd)
                else:
                    transformer.SetInput(pd)
                transformer.SetTransform(trans)
                transformer.Update()
                output_pds.append(transformer.GetOutput())

            appender = vtk.vtkAppendPolyData()
            for pd in output_pds:
                if (vtk.vtkVersion().GetVTKMajorVersion() >= 6.0):
                    appender.AddInputData(pd)
                else:
                    appender.AddInput(pd)
            appender.Update()
            wma.io.write_polydata(appender.GetOutput(), os.path.join(self.output_directory, 'registration_atlas.vtk'))
            del appender

            # save the transforms to text files
            wma.registration_functions.write_transforms_to_itk_format(transform_list, outdir, subject_id_list)
            
def congeal_multisubject_inner_loop(mean, subject, initial_transform, mode, sigma, subject_idx, iteration_count, output_directory, step_size, maxfun, render, grid_resolution):

    """This is the code executed by each subprocess that launches the

    registration of one subject to the current atlas model or mean brain.
    """
    
    #print "\n BEGIN ITERATION", iteration_count, "subject", subject_idx, "sigma:", sigma, "mean brain:", mean.shape, "subject:", subject.shape, "initial transform length:", len(initial_transform), "steps:", step_size[0], step_size[1], "maxfun:", maxfun, type(initial_transform), "Grid:", grid_resolution, "Mode:", mode, "initial transform:", initial_transform,
    
    # Set up registration objects and parameters that are specific to affine vs nonlinear
    if mode == 'Linear':
        register = wma.register_two_subjects.RegisterTractography()
        register.process_id_string = "_subject_%05d_iteration_%05d_sigma_%03d" % (subject_idx, iteration_count, sigma) 

    elif mode == "Nonlinear":
        register = wma.register_two_subjects_nonlinear.RegisterTractographyNonlinear()
        register.nonlinear_grid_resolution = grid_resolution
        register.initialize_nonlinear_grid()
        register.process_id_string = "_subject_%05d_iteration_%05d_sigma_%03d_grid_%03d" % (subject_idx, iteration_count, sigma, grid_resolution) 

    else:
        print "ERROR: Unknown registration mode"

    # Set up parameters that are used for both affine and nonlinear
    register.maxfun = maxfun
    register.sigma = sigma        
    register.parallel_jobs = 1
    register.threshold = 0
    register.fixed = mean
    register.moving = subject
    register.initial_transform = initial_transform
    register.verbose = False
    register.output_directory = output_directory
    register.initial_step = step_size[0]
    register.final_step = step_size[1]
    register.render = render

    # Run the current iteration of optimization.    
    register.compute()

    # Only update the transform if the objective function improved.
    # With affine registration, some subjects may have converged already to the current model.
    obj_diff = register.objective_function_values[-1] - register.objective_function_values[0]
    #print "\n END ITERATION", iteration_count, "subject", subject_idx, "OBJECTIVE CHANGE:", obj_diff
    if obj_diff < 0:
        #print "UPDATING MATRIX"
        return register.final_transform, register.objective_function_values, obj_diff
    else:
        return initial_transform, register.objective_function_values, obj_diff

