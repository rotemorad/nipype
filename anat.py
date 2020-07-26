import os
import os.path as op

from nipype.interfaces.fsl.base import (FSLCommand, isdefined)
from nipype.interfaces.fsl.utils import split_filename


class fsl_anat(FSLCommand):
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
    input_spec = None  # fsl_anatInputSpec
    output_spec = None  # fsl_anatOutputSpec

    def _list_outputs(self):
        outputs = self.output_spec().get()
        # Reorient the images to the standard orientation (fslreorient2std):
        if not self.inputs.no_orient:
            outputs["reorient2std"] = self._gen_fname("fullfov")
            outputs["original_input"] = self._gen_fname("orig")

        # Automatically crop the image (robustfov):
        if not self.inputs.no_crop:
            outputs["out_roi"] = self._gen_fname("")
            outputs["out_transform"] = self._gen_fname("out_transform")
            pass

        # Bias-field correction (FAST):
        if not self.inputs.no_bias:
            outputs["biascorr"] = self._gen_fname("biascorr")

        # Registration to standard space (FLIRT and FNIRT):
        if not self.inputs.no_reg:
            outputs["lin"] = self._gen_fname("to_MNI_lin")
            outputs["nonlin"] = self._gen_fname("to_MNI_nonlin")
            outputs["nonlin_field"] = self._gen_fname("to_MNI_nonlin_field")
            outputs["nonlin_jac"] = self._gen_fname("to_MNI_nonlin_jac")
            outputs["vols"] = self._gen_fname("vols.txt")

        # Brain-extraction (FNIRT-based or BET):
        outputs["biascorr_brain"] = self._gen_fname("biascorr_brain")
        outputs["biascorr_brain_mask"] = self._gen_fname("biascorr_brain_mask")

        # If tissue-type segmentation (FAST) occurs:
        if not self.inputs.no_seg:
            _fast_gen_fname_opts = {}
            if isdefined(self.inputs.out_basename):
                _fast_gen_fname_opts["basename"] = self.inputs.out_basename
                _fast_gen_fname_opts["cwd"] = os.getcwd()
            else:
                _fast_gen_fname_opts["basename"] = self.inputs.in_files[-1]
                _fast_gen_fname_opts["cwd"], _, _ = split_filename(_fast_gen_fname_opts["basename"])

            outputs["tissue_class_map"] = self._gen_fname(suffix="_seg", **_fast_gen_fname_opts)
            if isdefined(self.inputs.output_biascorrected):
                outputs["fast_bias"] = []
                if len(self.inputs.in_files) > 1:
                    # for multi-image segmentation there is one corrected image
                    # per input
                    for val, f in enumerate(self.inputs.in_files):
                        # image numbering is 1-based
                        outputs["fast_bias"].append(
                            self._gen_fname(**_fast_gen_fname_opts)
                        )
                else:
                    # single image segmentation has unnumbered output image
                    outputs["fast_bias"].append(
                        self._gen_fname(**_fast_gen_fname_opts)
                    )

            outputs["mixeltype"] = self._gen_fname(suffix="_mixeltype", **_fast_gen_fname_opts)
            outputs["partial_volume_map"] = self._gen_fname(suffix="_pveseg", **_fast_gen_fname_opts)
            outputs["partial_volume_files"] = []
            for i in range(3):
                outputs["partial_volume_files"].append(
                    self._gen_fname(suffix="_pve_%d" % i, **_fast_gen_fname_opts)
                )

            outputs["bias_field"] = []
            if len(self.inputs.in_files) > 1:
                # for multi-image segmentation there is one bias field image
                # per input
                for val, f in enumerate(self.inputs.in_files):
                    # image numbering is 1-based
                    outputs["bias_field"].append(
                        self._gen_fname(**_fast_gen_fname_opts)
                    )
            else:
                # single image segmentation has unnumbered output image
                outputs["bias_field"].append(
                    self._gen_fname(**_fast_gen_fname_opts)
                )
            return
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
            _first_gen_fname_opts = {}
            if isdefined(self.inputs.out_basename):
                _first_gen_fname_opts["basename"] = self.inputs.out_basename
                _first_gen_fname_opts["cwd"] = os.getcwd()
            else:
                _first_gen_fname_opts["basename"] = self.inputs.in_files[-1]
                _first_gen_fname_opts["cwd"], _, _ = split_filename(_first_gen_fname_opts["basename"])

            outputs["original_segmentations"] = self._gen_fname(suffix="_all_origsegs", **_first_gen_fname_opts)
            outputs["segmentation_file"] = self._gen_fname(suffix="_all_firstseg", **_first_gen_fname_opts)
            outputs["subcort_seg"] = self._gen_fname(suffix="_subcort_seg", **_first_gen_fname_opts)
            outputs["biascorr_to_std_sub"] = self._gen_fname(suffix="_biascorr_to_std_sub", **_first_gen_fname_opts)
            outputs["vtk_surfaces"] = self._gen_mesh_names("vtk_surfaces", structures)
            outputs["bvars"] = self._gen_mesh_names("bvars", structures)
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
