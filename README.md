# TomoStarPy

Some useful Python scripts to manipulate Relion tomo STAR files and much more.....


## Table of Contents
- [Installation](#installation)
- [Usage](#usage)
- [convmap_mod_get_particles.py](#convmap_mod_get_particlespy)
- [create_order_list.py](#create_order_listpy)
- [ctffind_tilts.py](#ctffind_tiltspy)
- [dose_correct_tomostar.py](#dose_correct_tomostarpy)
- [emc_csv_to_star.py](#emc_csv_to_starpy)
- [emc_ctf_tlt_to_xf.py](#emc_ctf_tlt_to_xfpy)
- [eulers_zyz_from_matrix.py](#eulers_zyz_from_matrixpy)
- [exclude_lines_by_range_file.py](#exclude_lines_by_range_filepy)
- [exclude_tilts.py](#exclude_tiltspy)
- [exclude_tilts_from_warp2_ts_xml.py](#exclude_tilts_from_warp2_ts_xmlpy)
- [filter_neighbors_only.py](#filter_neighbors_onlypy)
- [find_bad_ctf_fit_tilts.py](#find_bad_ctf_fit_tiltspy)
- [find_dark_tilts.py](#find_dark_tiltspy)
- [merge_tilts.py](#merge_tiltspy)
- [placeback_arrows.py](#placeback_arrowspy)
- [placeback_pdb.py](#placeback_pdbpy)
- [placeback_polygons.py](#placeback_polygonspy)
- [placeback_subvolume.py](#placeback_subvolumepy)
- [remove_duplicates_2_star.py](#remove_duplicates_2_starpy)
- [rotate_biomt_matrices.py](#rotate_biomt_matricespy)
- [rotate_coordinates_xy.py](#rotate_coordinates_xypy)
- [rotate_subtomo_along_filament.py](#rotate_subtomo_along_filamentpy)
- [rotate_xf.py](#rotate_xfpy)
- [select_cmm_mod_from_star.py](#select_cmm_mod_from_starpy)
- [shift_coordinates_by_tiltcom.py](#shift_coordinates_by_tiltcompy)
- [sort_by_distance.py](#sort_by_distancepy)
- [split_particles_star_to_single_tomos_star.py](#split_particles_star_to_single_tomos_starpy)
- [star_particle_editor.py](#star_particle_editorpy)
- [star_to_emc.py](#star_to_emcpy)
- [sub_subtomo.py](#sub_subtomopy)

## Installation

``` 
$ git clone https://github.com/fuzikt/tomostarpy.git
$ cd tomostarpy
$ ./INSTALL.sh
```

## Usage

Before usage don't forget to source the venv environment if it is not yet activated:
``` 
$ source <path_to_tomostarpy>/venv/bin/activate
``` 

## convmap_mod_get_particles.py

Get particle parameters (cc, coordinates, Euler angles) from a convmap file at positions defined by a mod file. Output is written in emClarity csv format.

```
--convmap    Input _convmap.mrc file.
--mod        Input mod file (Default: empty).
--o          Output emClarity csv file (Default: empty).
```

## create_order_list.py

Generates order_list.csv for RELION from *.tlt file according to dose-symmetric tilting scheme. (Default 2 pos, 2 neg generates: 0,1,-1,-2,2,3,-3,4...).

```
--i           Input *.tlt file.
--o           Output csv file.
--nr_pos      Number of consecutive positive tilts in output (Default: 2).
--nr_neg      Number of consecutive negative tilts in output (Default: 2).
--pacetomo    Pace-tomo style dose symmetric (0, ++, --, ++, --.....) (Default: False).
```

## ctffind_tilts.py

Performs constrained CTF fit on separate tilts using Ctffind4 or Gctf. First, the middle tilt ctf is estimated then the following tilts in each direction are constrained by the --defocus_window of the previous tilt. The max resolution of the fitting is limited by 1/cos(tilt_angle) to avoid overfitting.
```
--i                 Input tilt-series MRC file.
--o                 Output directory of results.
--tlt               *.tlt or *.rawtlt file of the tilt series.
--defocus_window    Defocus window in Angstroms used for the restrained ctf search of the subsequent tilts (Default: 10000).
--apix              Pixel-size in Angstroms (Default: 1.0).
--voltage           Acceleration voltage of the microscope in kV (Default: 300.0).
--cs_val            Cs value of the microscope (Default: 2.7).
--amp_cont          Amplitude contrast used for CTF fitting (Default: 0.07).
--spectr_size       Size of the Fourier spectrum used for CTF fitting (Default: 512).
--min_res           Minimum resolution in Angstroms used for zero-tilt CTF fitting (Default: 30.0).
--max_res           Maximum resolution in Angstroms used for zero-tilt CTF fitting (Default: 8.0).
--min_defoc         Minimum defocus in Angstroms used for zero-tilt CTF fitting (Default: 5000.0).
--max_defoc         Maximum defocus in Angstroms used for zero-tilt CTF fitting (Default: 75000.0).
--step_defoc        Defocus step in Angstroms used for CTF fitting (Default: 100.0).
--threads           Number of parallel threads used for calculation by Ctffind4 (Default: 10).
--gctf              Use gCTF instead of CtfFind4 (Default: False).
--gpu               GPU card ID to be used by gCTF (Default: 0).
--tmpDir            Temp directory for storage of the intermediate files (Default: ctffind_tmp).
--ctffind_cmd       Name of the Ctffind4 command (Default: ctffind).
--gctf_cmd          Name of the gCTF command (Default: gctf).
--path              Path to be included in env PATH for Ctffind4 or Gctf (Default: empty).
--verb              Verbosity level (0,1) (Default: 0).
```

## dose_correct_tomostar.py

Replace dose in warp .tomostar by value from .tltdose file.
```
--i        Input tomostar file.
--itlt     Input dosetlt file.
--o        Output tomostar file.
```

## emc_csv_to_star.py

Performs conversion from emClarity template matching csv file to Relion star file. It also recalculates the particle coordinates from emClarity partial tomogram to full tomogram.

```
--i        Input emClarity csv file.
--o        Output prefix. Prefix of the files generated by the script.
--reconsh  emClarity *_recon.sh file of the input csv file. (located in emc_project/recon/)
--tiltcom  tilt.com file (generated by IMOD) of the original non-splitted tomogram used by Relion
--mod      Modfile used for filtering particles from the input csv. Note: Set the "--modbin" binning factor used for template matching.
--modbin   Binning factor used for template matching. Needed only if modfile filtering is enabled. (Default 1)
--outbin   Binning factor used for output coordinates and mod files. NOT used for star file (always unbinned)! (Default 1)
--cs       Cs value of the microscope. Used in opticsgroup in the output star file. (Default 2.7)
--kv       Acceleration voltage of the microscope. Used in opticsgroup in the output star file. (Default 300.0)
--apix     Apix of the unbinned tomogram. Used in opticsgroup in the output star file.
--xtilt    Tomo X axis tilt in degrees. (Default: 0)
--ytilt    Tomo Y axis tilt in degrees. (Default: 0)
```

## emc_ctf_tlt_to_xf.py

Converts emClarity *_aliX_ctf.tlt from tomoCPR to IMOD format *.xf + *.tlt file and Ctffind4 format diag. output *.txt file.
```
--i        Input emClarity *_aliX_ctf.tlt file.
--o        Output prefix for IMOD style xf/tlt file and Ctffind4 style diag. output file.
```

## eulers_zyz_from_matrix.py

Calculates ZYZ convention Euler angles from a 3x3 rotation matrix.
```
--i        Rotation matrix members in format: m[1,1],m[1,2],m[1,3],m[2,1],m[2,2]
```

## exclude_lines_by_range_file.py

Removes lines defined in --exclude_file from input text file.
```
--i             Input file to exclude lines from.
--o             Output file.
--exclude_file  Textfile with comma separated list of excluded views (or range of views as in newstack).
```
## exclude_tilts.py

Removes tilts defined in --exclude_file from tilt series mrc-stack and corresponding *.rawtlt, *.mdoc.
```
--i             Input prefix mrc-stack.
--o             Output prefix.
--exclude_file  Textfile with comma separated list of excluded views (or range of views as in newstack).
```

## exclude_tilts_from_warp2_ts_xml.py
Removes tilts defined in --exclude_file from Warp2 tilt-series XML file and tomostar file.
```
    --i_xml         Input Warp2 tilt-series XML file.
    --i_tomostar    Input tomostar file.
    --o             Output prefix for XML and tomostar file.
    --exclude_file  Textfile with comma separated list of excluded views (or range of views as used in imod newstack).
```

Example 1: Exclude tilts from tomo1.xml and tomo1.tomostar files, the tilts to be excluded are defined in tomo1_excluded.txt.
```
mkdir excluded_tilts

exclude_tilts_from_warp2_ts_xml.py --i_xml warp_tiltseries/tomo1.xml --i_tomostar tomostar/tomo1.tomostar --o excluded_tilts/tomo1 --exclude_file tomo1_excluded.txt

# This will create excluded_tilts/tomo1.xml and excluded_tilts/tomo1.tomostar with the tilts defined in tomo1_excluded.txt removed.

Example of the content of tomo1_excluded.txt:
1,2,12-16,20,22-24
```

Example 2: A simple bash loop to exclude tilts from multiple tomograms (can be written as a script or "one-liner" in the terminal):
```
for i in `seq 1 10`; do
    exclude_tilts_from_warp2_ts_xml.py --i_xml warp_tiltseries/tomo${i}.xml --i_tomostar tomostar/tomo${i}.tomostar --o excluded_tilts/tomo${i} --exclude_file tomo${i}_excluded.txt
done
```

## filter_neighbors_only.py

Select only particles that have a certain amount of neighbors in a particular distance. Useful for filtering out template matched particles.
```
--i               Input STAR file name with particles.
--o               Output STAR file name.
--dist            Distance in Angstroms that consider particles as neighbors. (Default 1.0)
--min_neigh       Minimum number of neighbors at --dist particle has to be kept in selection! (Default 1)
--min_corr        Minimum cross-correlation value of the neighbour to be considered as a true neighbour! (Default 0)
--lb_corr'        Label of the cross-correlation value in the star file! (Default rlnLCCmax)
--max_ang_dist    Maximum angular distance in degrees between neighbors to be considered as a true neighbour! (Default 360)
```

## find_bad_ctf_fit_tilts.py

Analyze the Ctffind4 output txt file and write out text file with tilt-numbers that have CTF rings fitted over the threshold. Tilts are numbered from 1.
```
--i           Input Ctffind4 diagnostic txt file.
--o           Output text file with tilt numbers under the threshold.
--threshold   Threshold (in Angstroms) value for the resolution up to which CTF rings were fit successfully. (Default: 100)
```

## find_dark_tilts.py

Analyze the input tilt-series mrc stack and write out text file with tilt-numbers that have average signal under the threshold. Tilts are numbered from 1. Useful to remove dark tilts.
```
--i           Input tilt-series mrc-stack.
--o           Output text file with tilt numbers under the threshold.
--threshold   Threshold value for the average value of the tilt. (Default: 0.5)
```

## merge_tilts.py

Performs merging of separate tilt files from PACE-tomo or dose symmetric SerialEM script into tilt-series files.
```
--i                Input tilts directory.
--o                Output directory of merged tilt series.
--start_tilt       Angle of the "zero-tilt" - at which the tilt-series begin. Useful for pretilted tilt-series. (Default: 0)
--start_series_nr  The sequential number of the start nr. of tilt series (used for continuing the process). (Default: 1)
--prefix           Prefix of the output tilt series file. (Default: tomo)
--mdoc_dir         Directory where *.mdoc files are located. If empty, no *.rawtlt and merged *.mdoc are generated (Default: empty)
--mdoc_suffix      Suffix of the *.mdoc files after the mic name. (Default: .tif.mdoc)
--pre_dose         Dose applied before acquiring the first tilt. [e-/A^2]. (Default: 0.0)
--dose             Dose applied per tilt. [e-/A^2]. (Default: 0.0)
--apix             Pixel size of the merged tilt-series. [A/px]. (Default: 1.0)
--newstack_cmd     Name of the IMOD newstack command. (Default newstack)
--alterheader_cmd  Name of the IMOD alterheader command. (Default alterheader)
--path             Path to be included in env PATH for IMOD! (Default: empty)
--watch_int        Time interval in seconds for the watchdog. If set to -1 the watchdog is deactivated and runs a single-batch. (Default: 1)
--verb             Verbosity level (0,1,2). (Default: 1)
```

## placeback_arrows.py

Performs placeback of oriented arrows (*.bild format) according to the coordinates and euler angles defined in a star file.
```
--i           Input star file file.
--o           Output prefix
--length      Length of the arrow in Angstroms. (default: 200)
--thickness   Radius of the arrow base in Angstroms. (default: 1/20 of the arrow length)
--cmm         Generate cmm file as well.
--color_lb    Label from the star file that will be used for rainbow coloring of the arrows and cmm markers. (default: empty)
--invert      Invert the pointing direction of the arrow.
--tomo_name   Use only particles from tomogram equal in rlnTomoName. OPTIONAL: If not set all particles used in place-back
```

## placeback_pdb.py

Performs placeback of a PDB structure into full tomogram volume according to the coordinates and euler angles defined in a star file.
```
--i              Input star file file.
--i_pdb          Input star file file.
--o              Output prefix
--center_offset  PDB center offset. (default: 0,0,0)
--backbone       Store only backbone atoms in the output.
--tomo_name      Use only particles from tomogram equal in rlnTomoName. OPTIONAL: If not set all particles used in place-back
--separate_pdb   Write separate PDB file per particle
```

## placeback_polygons.py

Performs placeback of polygons (*.bild format) according to the coordinates and euler angles defined in a star file.
```
--i            Input star file file.
--o            Output prefix
--n            Order of the polygon (3 - triangle; 4 - square; 5 - pentagon etc. (default: 3))
--size         Distance between the center and the vertex in Angstroms. (default: 200)
--pre_rot      Pre rotation (in degrees) applied to the polygon before placing back. (default: 0)
--cmm          Generate cmm file as well.
--color_lb     Label from the star file that will be used for rainbow coloring of the arrows and cmm markers. (default: empty)
--tomo_name    Use only particles from tomogram equal in rlnTomoName. OPTIONAL: If not set all particles used in place-back
--single       Creates a single polygon with applied ---pre_rot and of size --size. Useful to find the --pre_rot parameter. No input star file needed.
```

## placeback_subvolume.py

Performs placeback of a subvolume into full tomogram volume according to the coordinates and euler angles defined in a star file.
```
--i                      Input star file file.
--isub                   Input file of sub-volume to be placed.
--itomo                  Input file of the stencil tomogram. Size and origin is taken from this file and applied on the output.
--o                      Output prefix
--cmm                    Create Chimera cmm file with the coordinates of the placed sub-volumes.
--tomo_name              Use only particles from tomogram equal in rlnTomoName. OPTIONAL: If not set all particles used in place-back
--bin                    Binning factor of the stencil tomogram. If not provided calculated from MRC header and star file
--no_partial             Do not place partial volumes. If the subvolume is partially out of the output volume, it is not placed at all. This avoids "half-cut" sub-volumes in the output.
--recenter               Recenter the particles by subtracting rlnOriginX/Y/ZAngst from X/Y/Z coordinates.
--color_lb               Label from the star file that will be used for rainbow coloring of the cmm markers.
--color_map              Create map with coloring values stored as pixel value.
--radial_color           Create map with coloring values storing the value of radial distance form the center of the box. Useful for radial coloring.
--color_map_threshold    Threshold value at which the contour of the --isub in the output color map should contain values
--color_map_extend       Extend the border around the threshold meeting value by this amount of pixels.
--xtilt                  Tomo X axis tilt in degrees. (Default: 0)
--ytilt                  Tomo Y axis tilt in degrees. (Default: 0)
--gpu                    Use GPU acceleration.
```

## remove_duplicates_2_star.py

Remove duplicate particles appearing both in --i1 --i2 from --i2 and write result into --o.
```
--i1    Input1 STAR filename.
--i2    Input2 STAR filename.
--o     Output STAR filename.
--tol   Tolerance in pixels of considering the coordinates of particles to be the same (Default: 3).
```

## rotate_biomt_matrices.py

Rotate and shift REMARK 350 BIOMT matrices from PDB file, using a given matrix or ZYZ convention euler angles.
```
--i      Input PDB file.
--o      Output file.
--i_mat  Input file with rotation and translation matrix in BIOMT format used for the rotation. If defined --rot, --tilt, and --psi are ignored.
--rot    Euler angle ROT in degrees. (Default 0.0)
--tilt   Euler angle TILT in degrees. (Default 0.0)
--psi    Euler angle PSI in degrees. (Default 0.0)
--x      X shift in Angstroms. (Default 0.0)
--y      Y shift in Angstroms. (Default 0.0)
--z      Z shift in Angstroms. (Default 0.0)
```

## rotate_coordinates_xy.py

Rotate the inplane coordinates X, Y of the subvolume.
```
--i     Input star file.
--o     Output star file.
--ang   Rotation angle in degrees. (Default 0.0)
--xmax  X-size of the unbinned tomogram. (Default 4092)
--ymax  Y-size of the unbinned tomogram. (Default 5760)
```

## rotate_subtomo_along_filament.py

Set the Euler angles (tilt, psi) of the particles along a filament according to the angles of the between the neighboring particle centers.
```
--i   Input STAR filename with particles.
--o   Output STAR filename.
```

## rotate_xf.py

Rotate the values in IMOD xf file to change the rotation of the finally aligned stack.
```
--i   Input IMOD xf file.
--o   Output IMOD xf file.
--ang Rotation angle in degrees. (Default 0.0)
```

## select_cmm_mod_from_star.py

Select particles from star file according to matching particle coordinates listed in Chimera cmm or IMOD mod file.
```
--i    Input STAR filename with particles.
--o    Output STAR filename.
--cmm  Chimera CMM file with desired coordinates.
--mod  IMOD mod file with desired coordinates.
--bin  Binning factor of the IMOD mod file.
```

## shift_coordinates_by_tiltcom.py

Shift particle coordinates in input STAR file, by recalculating the difference between the two reconstructions defined in --tiltcom1 and --tiltcom2.
```
--i         Input STAR file name with particles.
--o         Output STAR file name.
--tiltcom1  tilt.com file (generated by IMOD) of the tomogram of the input coordinates
--tiltcom2  tilt.com file (generated by IMOD) of the tomogram for the output coordinates
--xtilt     Tomo X axis tilt in degrees. (Default: 0)
--ytilt     Tomo Y axis tilt in degrees. (Default: 0)
```

## sort_by_distance.py

Sort particles in star file to order based on the minimal distance between consecutive particles. In general, this should form a filament path.
```
--i   Input STAR file name with particles.
--o   Output STAR file name.
```

## split_particles_star_to_single_tomos_star.py

Split a particle star file into separate files per tomoName containing the corresponding particles.
```
--i   Input star file.
--o   Output directory name.
```

## star_particle_editor.py

Visual editor to add, remove, combine particles from tomo STAR files.
```
--i          Input star file. Possible to open multiple star files delimited by comma.
--itomo      Input mrc file. Possible to open multiple mrc files delimited by comma.
--o          Output star file name. Possible to change later in the gui.
--tomo_name  Use only particles from tomogram equal in rlnTomoName. OPTIONAL: If not set all particles will be visualized.
--color_lb   Label from the star file that will be used for rainbow coloring of the markers.
--bin        User provided binning of the --itomo. If not set the apix of the first MRC file in --itomo is taken to calculate the binning.
--point_size Size of the points in pixels.
--angles_mrc emClarity *_angles.mrc from template matching, to be used for newly added particles orientations. Same named *_angles.list must be present in the directory.
```

## star_to_emc.py

Performs conversion from Relion star file to emClarity template matching csv file. It also generates *.recon.coords files for tomogram subareas.
```
--i       Input particle star file.
--o       Output prefix. Prefix of the files generated by the script.
--itomo   Input tomogram star file.
--splitX  Number of subareas to split the tomograms on X axis. Default: 2
--splitY  Number of subareas to split the tomograms on X axis. Default: 2
--bin     Binning factor used for the output mod files. Default: 1
```

## sub_subtomo.py

Creates a star file with sub-subtomo coordinates similar to sub-particle approach in SPA.
```
--i                     Input star file.
--o                     Output star file.
--cmm                   A CMM file defining the location(s) of the subparticle(s) (use instead of --vector). Coordinates should be in Angstrom.
--vector                Vector defining the location of the subparticle. (default: 0,0,1)
--length                Alternative length of the vector. Use to adjust the subparticle center (A). (default: length of the given vector)
--sym                   Symmetry of the particle. (default: C1)
--align_subparticles    Align subparticles to the standard orientation.
--randomize             Randomize the order of the symmetry matrices. Useful for preventing preferred orientations.
--unique                Keep only unique subparticles within angular distance (useful to remove overlapping subparticles on symmetry axis).
--mindist               Minimum distance between the subparticles in the image (all overlapping ones will be discarded; pixels).
--side                  Keep only particles within specified angular distance from side views (all others will be discarded; degrees).
--top                   Keep only particles within specified angular distance from top views (all others will be discarded; degrees).
--library_path          Define LD_LIBRARY_PATH used. Default: empty
```