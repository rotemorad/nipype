import os
import os.path as op

from nipype.interfaces.fsl.base import (isdefined)
from nipype.interfaces.fsl.utils import split_filename

from .base import FSLCommand, FSLCommandInputSpec, isdefined
from ..base import (
    TraitedSpec,
    File,
    Directory,
    traits, OutputMultiPath,
)
from ...utils.filemanip import split_filename


class FSLAnatInputSpec(FSLCommandInputSpec):
    input_img = File(
        exists=True,
        desc="filename of input image (for one image only)",
        argstr="-i %s",
        position=-1,
        mandatory=True,
        xor=['directory']
    )

    directory = Directory(
        exists=True,
        desc="directory name for existing .anat directory where this script will be run in place",
        argstr="-d %s",
        position=-1,
        mandatory=True,
        xor=['input_img']
    )

    output_directory = Directory(
        exists=True,
        desc='basename of directory for output (default is input image basename followed by .anat)',
        argstr='-o %s',
        mandatory=False,
    )

    weakbias = traits.Bool(desc="used for images with little and/or smooth bias fields",
                           argstr="--weakbias"
                           )

    clobber = traits.Bool(
        desc="if .anat directory exist (as specified by -o or default from -i) then delete it and make a new one",
        argstr="--clobber"
    )

    noreorient = traits.Int(desc="turn off step that does reorientation 2 standard (fslreorient2std)",
                            argstr="--noreorient %d")

    nocrop = traits.Bool(desc="turn off step that does automated cropping (robustfov)",
                         argstr="--nocrop"
                         )

    nobias = traits.Bool(desc="turn off steps that do bias field correction (via FAST)",
                         argstr="--nobias"
                         )

    noreg = traits.Bool(desc="turn off steps that do registration to standard (FLIRT and FNIRT)",
                        argstr="--noreg"
                        )

    nononlinreg = traits.Bool(desc="turn off step that does non-linear registration (FNIRT)",
                              argstr="--nononlinreg"
                              )

    noseg = traits.Bool(desc="turn off step that does tissue-type segmentation (FAST)",
                        argstr="--noseg"
                        )

    nosubcortseg = traits.Bool(desc="turn off step that does sub-cortical segmentation (FIRST)",
                               argstr="--nosubcortseg"
                               )

    image_type = traits.Str('image_type',
                            desc="specify the type of image (choose one of T1 T2 PD - default is T1)",
                            argstr="-t %s"
                            )

    s = traits.Int(desc="specify the value for bias field smoothing (the -l option in FAST)",
                   argstr="-s %d")

    nosearch = traits.Bool(desc="specify that linear registration uses the -nosearch option (FLIRT)",
                           argstr="--nosearch"
                           )

    betfparam = traits.Float(
        desc="specify f parameter for BET (only used if not running non-linear reg and also wanting brain extraction "
             "done)",
        argstr="--betfparam %.2f"
    )

    nocleanup = traits.Bool(desc="do not remove intermediate files",
                            argstr="--nocleanup"
                            )


class FSLAnatOutputSpec(TraitedSpec):
    Out = File(exists=False, extentions='.nii',
               desc=r"Contains either the T1, T2 or PD image (according to the -t input option) after cropping and\or "
                    r"orientation.")
    Out_orig = File(exists=False, extentions='.nii',
                    desc=r"The original image (exists if the image was cropped and\or reoriented)")
    Out_fullfov = File(exists=False, extentions='.nii',
                       desc=r"The image in full field-of-view (exists if the image was cropped and\or reoriented)")
    Out_orig2std = File(exists=False, extentions='.mat', desc="Transformation matrix to allow images to be moved "
                                                              "between original space and standard space")
    Out_std2orig = File(exists=False, extentions='.mat', desc="Transformation matrix to allow images to be moved "
                                                              "between standard space and original space ")
    Out_orig2roi = File(exists=False, extentions='.mat', desc="Transformation matrix to allow images to be moved "
                                                              "between original space and the region of interest")

    Out_roi2orig = File(exists=False, extentions='.mat', desc="Transformation matrix to allow images to be moved ")

    Out_roi2nonroi = File(exists=False, extentions='.mat', desc="Transformation matrix to allow images to be moved "
                                                                "between region of interest and non- region of "
                                                                "interest "
                                                                "(full field of view) in the standard space")
    Out_nonroi2roi = File(exists=False, extentions='.mat', desc="Transformation matrix to allow images to be moved "
                                                                "between non- region of interest (full field of view) "
                                                                "and region of interest in the standard space")
    Out_biascorr = File(exists=False, extentions='.nii', desc="The estimated restored input image after correction "
                                                              "for bias field, if tissue type segmentation occurs it "
                                                              "is refined again")
    Out_to_MNI_lin = File(exists=False, extentions='.nii', desc="Linear registration output")
    Out_to_MNI_nonlin = File(exists=False, extentions='.nii', desc="Non-linear registration output")
    Out_to_MNI_nonlin_field = File(exists=False, extentions='.nii', desc="Non-linear warp field")
    Out_to_MNI_nonlin_jac = File(exists=False, extentions='.nii', desc="Jacobian of the non-linear warp field")
    Out_vols = File(exists=False, extentions='.txt',
                    desc="A file containing a scaling factor and brain volumes, based on skull-contrained "
                         "registration, "
                         "suitable for head-size normalisation (as the scaling is based on the skull size, "
                         "not the brain size")
    Out_biascorr_brain = File(exists=False, extentions='.nii',
                              desc="The estimated restored input image after correction "
                                   "for bias field and brain extraction")
    Out_biascorr_brain_mask = File(exists=False, extentions='.nii', desc="The estimated restored input image after "
                                                                         "correction for bias field and extraction of "
                                                                         "brain mask")
    Out_fast_pve_0 = File(exists=False, extentions='.nii', desc="Cerebral spinal fluid segmentation")
    Out_fast_pve_1 = File(exists=False, extentions='.nii', desc="Gray matter segmentation")
    Out_fast_pve_2 = File(exists=False, extentions='.nii', desc="White matter segmentation")
    Out_fast_pveseg = File(exists=False, extentions='.nii', desc="A summary image showing the tissue with the "
                                                                 "greatest partial volume fraction per voxel")
    Out_subcort_seg = File(exists=False, extentions='.nii', desc="Summary image of all sub-cortical segmentations")
    Out_first_all_fast_firstseg = File(exists=False, extentions='.nii', desc="Summary image of all sub-cortical "
                                                                             "segmentations")
    Out_biascorr_to_std_sub = File(exists=False, extentions='.nii', desc="A transformation matrix of the sub-cortical "
                                                                         "optimised MNI registration")
    first_results = Directory(exists=False, desc="")  # TODO: add path?
    Out_vtk_surfaces = OutputMultiPath(File(exists=False), desc="VTK format meshes for each subcortical region")
    Out_bvars = OutputMultiPath(File(exists=False), desc="bvars for each subcortical region")


class FSLAnat(FSLCommand):
    """FSL fsl_anat wrapper for pipeline to processing anatomical images (e.g. T1-weighted scans).
    https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/fsl_anat
    Examples
    --------

    >>> from nipype.interfaces import fsl
    >>> fsl_anat = fsl.fsl_anat()
    >>> fsl_anat.inputs.in_file = 'structural.nii'
    >>> res = fsl_anat.run() #doctest: +SKIP
    """

    _cmd = "fsl_anat"
    input_spec = FSLAnatInputSpec
    output_spec = FSLAnatOutputSpec

    def _list_outputs(self):
        outputs = self.output_spec().get()
        basename = self.inputs.input_img
        cwd = self.inputs.output_directory
        if not isdefined(self.inputs.output_directory):
            cwd = os.path.abspath(self.inputs.input_img)
        kwargs = {'basename': basename, 'cwd': cwd}

        # Reorient the images to the standard orientation (fslreorient2std):
        if not self.inputs.no_orient:
            outputs["Out_fullfov"] = self._gen_fname(suffix="_fullfov", **kwargs)
            outputs["Out_orig"] = self._gen_fname(suffix="_orig", **kwargs)
            outputs["Out_orig2std"] = self._gen_fname(suffix="_orig2std", **kwargs)
            outputs["Out_std2orig"] = self._gen_fname(suffix="_std2orig", **kwargs)

        # Automatically crop the image (robustfov):
        if not self.inputs.no_crop:
            outputs["Out"] = self._gen_fname(**kwargs)
            outputs["Out_orig2roi"] = self._gen_fname(suffix="_orig2roi", **kwargs)
            outputs["Out_roi2orig"] = self._gen_fname(suffix="_roi2orig", **kwargs)

        if not self.inputs.no_crop and not self.inputs.no_orient:
            outputs["Out_roi2nonroi"] = self._gen_fname(suffix="_roi2nonroi", **kwargs)
            outputs["Out_nonroi2roi"] = self._gen_fname(suffix="_nonroi2roi", **kwargs)

        # Bias-field correction (FAST):
        if not self.inputs.no_bias:
            outputs["Out_biascorr"] = self._gen_fname(suffix="_biascorr", **kwargs)

        # Registration to standard space (FLIRT and FNIRT):
        if not self.inputs.no_reg:
            outputs["Out_to_MNI_lin"] = self._gen_fname(suffix="_to_MNI_lin", **kwargs)
            outputs["Out_to_MNI_nonlin"] = self._gen_fname(suffix="_to_MNI_nonlin", **kwargs)
            outputs["Out_to_MNI_nonlin_field"] = self._gen_fname(suffix="_to_MNI_nonlin_field", **kwargs)
            outputs["Out_to_MNI_nonlin_jac"] = self._gen_fname(suffix="_to_MNI_nonlin_jac", **kwargs)
            outputs["Out_vols"] = self._gen_fname(suffix="_vols", **kwargs)

        # Brain-extraction (FNIRT-based or BET):
        outputs["Out_biascorr_brain"] = self._gen_fname("biascorr_brain")
        outputs["Out_biascorr_brain_mask"] = self._gen_fname("biascorr_brain_mask")

        # If tissue-type segmentation (FAST) occurs:
        if not self.inputs.no_seg:
            outputs["Out_biascorr"] = self._gen_fname(suffix="_biascorr", **kwargs)  # Refined again in this stage
            outputs["Out_fast_pveseg"] = self._gen_fname(suffix="_fast_pveseg", **kwargs)
            outputs["Out_fast_pve_0"] = self._gen_fname(suffix="_fast_pve_0", **kwargs)
            outputs["Out_fast_pve_1"] = self._gen_fname(suffix="_fast_pve_1", **kwargs)
            outputs["Out_fast_pve_2"] = self._gen_fname(suffix="_fast_pve_2", **kwargs)

        # If sub-cortical segmentation (FIRST) occurs:
        if not self.inputs.no_subcort_seg:
            structures = [
                "L_Hipp",
                "R_Hipp",
                "L_Accu",
                "R_Accu",
                "L_Amyg",
                "R_Amyg",
                "L_Caud",
                "R_Caud",
                "L_Pall",
                "R_Pall",
                "L_Puta",
                "R_Puta",
                "L_Thal",
                "R_Thal",
                "BrStem",
            ]
            outputs["Out_subcort_seg"] = self._gen_fname(suffix="_subcort_seg", **kwargs)
            outputs["Out_biascorr_to_std_sub"] = self._gen_fname(suffix="_biascorr_to_std_sub", **kwargs)
            # The following are in a sub directory called first_results
            _first_gen_fname_opts = {"cwd": os.path.join(cwd, 'first_results')}
            outputs["Out_first_all_fast_firstseg"] = self._gen_fname(suffix="_first_all_fast_firstseg", **_first_gen_fname_opts)
            outputs["Out_vtk_surfaces"] = self._gen_mesh_names("vtk_surfaces", structures)
            outputs["Out_bvars"] = self._gen_mesh_names("bvars", structures)
            return outputs

    # TODO: note that everything but the mat and subcort_seg is inside a file called first_results

    def _gen_mesh_names(self, name, structures):
        path, prefix, ext = split_filename(self.inputs.out_file)
        if name == "vtk_surfaces":
            vtks = list()
            for struct in structures:
                vtk = prefix + "-" + struct + "_first.vtk"
                vtks.append(op.abspath(vtk))
            return vtks
        if name == "bvars":
            bvars = list()
            for struct in structures:
                bvar = prefix + "-" + struct + "_first.bvars"
                bvars.append(op.abspath(bvar))
            return bvars
        return None
