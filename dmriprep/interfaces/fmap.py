"""
Interfaces to deal with the various types of fieldmap sources
    .. testsetup::
        >>> tmpdir = getfixture('tmpdir')
        >>> tmp = tmpdir.chdir() # changing to a temporary directory
        >>> nb.Nifti1Image(np.zeros((90, 90, 60)), None, None).to_filename(
        ...     tmpdir.join('epi.nii.gz').strpath)
"""

import numpy as np
import nibabel as nb
from nipype import logging
from nipype.interfaces.base import (
    BaseInterfaceInputSpec, TraitedSpec, File, isdefined, traits,
    SimpleInterface)

LOGGER = logging.getLogger('nipype.interface')


class FieldToRadSInputSpec(BaseInterfaceInputSpec):
    in_file = File(exists=True, mandatory=True, desc='input fieldmap')
    fmap_range = traits.Float(desc='range of input field map')


class FieldToRadSOutputSpec(TraitedSpec):
    out_file = File(desc='the output fieldmap')
    fmap_range = traits.Float(desc='range of input field map')


class FieldToRadS(SimpleInterface):
    """
    The FieldToRadS converts from arbitrary units to rad/s
    """
    input_spec = FieldToRadSInputSpec
    output_spec = FieldToRadSOutputSpec

    def _run_interface(self, runtime):
        fmap_range = None
        if isdefined(self.inputs.fmap_range):
            fmap_range = self.inputs.fmap_range
        self._results['out_file'], self._results['fmap_range'] = _torads(
            self.inputs.in_file, fmap_range, newpath=runtime.cwd)
        return runtime


class FieldToHzInputSpec(BaseInterfaceInputSpec):
    in_file = File(exists=True, mandatory=True, desc='input fieldmap')
    range_hz = traits.Float(mandatory=True, desc='range of input field map')


class FieldToHzOutputSpec(TraitedSpec):
    out_file = File(desc='the output fieldmap')


class FieldToHz(SimpleInterface):
    """
    The FieldToHz converts from arbitrary units to Hz
    """
    input_spec = FieldToHzInputSpec
    output_spec = FieldToHzOutputSpec

    def _run_interface(self, runtime):
        self._results['out_file'] = _tohz(
            self.inputs.in_file, self.inputs.range_hz, newpath=runtime.cwd)
        return runtime


class Phasediff2FieldmapInputSpec(BaseInterfaceInputSpec):
    in_file = File(exists=True, mandatory=True, desc='input fieldmap')
    metadata = traits.Dict(mandatory=True, desc='BIDS metadata dictionary')


class Phasediff2FieldmapOutputSpec(TraitedSpec):
    out_file = File(desc='the output fieldmap')


class Phasediff2Fieldmap(SimpleInterface):
    """
    Convert a phase difference map into a fieldmap in Hz
    """
    input_spec = Phasediff2FieldmapInputSpec
    output_spec = Phasediff2FieldmapOutputSpec

    def _run_interface(self, runtime):
        self._results['out_file'] = phdiff2fmap(
            self.inputs.in_file,
            _delta_te(self.inputs.metadata),
            newpath=runtime.cwd)
        return runtime


def _despike2d(data, thres, neigh=None):
    """
    despiking as done in FSL fugue
    """

    if neigh is None:
        neigh = [-1, 0, 1]
    nslices = data.shape[-1]

    for k in range(nslices):
        data2d = data[..., k]

        for i in range(data2d.shape[0]):
            for j in range(data2d.shape[1]):
                vals = []
                thisval = data2d[i, j]
                for ii in neigh:
                    for jj in neigh:
                        try:
                            vals.append(data2d[i + ii, j + jj])
                        except IndexError:
                            pass
                vals = np.array(vals)
                patch_range = vals.max() - vals.min()
                patch_med = np.median(vals)

                if (patch_range > 1e-6 and
                        (abs(thisval - patch_med) / patch_range) > thres):
                    data[i, j, k] = patch_med
    return data


def _unwrap(fmap_data, mag_file, mask=None):
    from math import pi
    from nipype.interfaces.fsl import PRELUDE
    magnii = nb.load(mag_file)

    if mask is None:
        mask = np.ones_like(fmap_data, dtype=np.uint8)

    fmapmax = max(abs(fmap_data[mask > 0].min()), fmap_data[mask > 0].max())
    fmap_data *= pi / fmapmax

    nb.Nifti1Image(fmap_data, magnii.affine).to_filename('fmap_rad.nii.gz')
    nb.Nifti1Image(mask, magnii.affine).to_filename('fmap_mask.nii.gz')
    nb.Nifti1Image(magnii.get_data(), magnii.affine).to_filename('fmap_mag.nii.gz')

    # Run prelude
    res = PRELUDE(phase_file='fmap_rad.nii.gz',
                  magnitude_file='fmap_mag.nii.gz',
                  mask_file='fmap_mask.nii.gz').run()

    unwrapped = nb.load(res.outputs.unwrapped_phase_file).get_data() * (fmapmax / pi)
    return unwrapped


def get_ees(in_meta, in_file=None):
    """
    Calculate the *effective echo spacing* :math:`t_\\text{ees}`
    for an input :abbr:`EPI (echo-planar imaging)` scan.
    There are several procedures to calculate the effective
    echo spacing. The basic one is that an ``EffectiveEchoSpacing``
    field is set in the JSON sidecar. The following examples
    use an ``'epi.nii.gz'`` file-stub which has 90 pixels in the
    j-axis encoding direction.
    >>> meta = {'EffectiveEchoSpacing': 0.00059,
    ...         'PhaseEncodingDirection': 'j-'}
    >>> get_ees(meta)
    0.00059
    If the *total readout time* :math:`T_\\text{ro}` (``TotalReadoutTime``
    BIDS field) is provided, then the effective echo spacing can be
    calculated reading the number of voxels :math:`N_\\text{PE}` along the
    readout direction and the parallel acceleration
    factor of the EPI
      .. math ::
           =  T_\\text{ro} \\,  (N_\\text{PE} / f_\\text{acc} - 1)^{-1}
    where :math:`N_y` is the number of pixels along the phase-encoding direction
    :math:`y`, and :math:`f_\\text{acc}` is the parallel imaging acceleration factor
    (:abbr:`GRAPPA (GeneRalized Autocalibrating Partial Parallel Acquisition)`,
    :abbr:`ARC (Autocalibrating Reconstruction for Cartesian imaging)`, etc.).
    >>> meta = {'TotalReadoutTime': 0.02596,
    ...         'PhaseEncodingDirection': 'j-',
    ...         'ParallelReductionFactorInPlane': 2}
    >>> get_ees(meta, in_file='epi.nii.gz')
    0.00059
    Some vendors, like Philips, store different parameter names
    (see http://dbic.dartmouth.edu/pipermail/mrusers/attachments/\
20141112/eb1d20e6/attachment.pdf):
    >>> meta = {'WaterFatShift': 8.129,
    ...         'MagneticFieldStrength': 3,
    ...         'PhaseEncodingDirection': 'j-',
    ...         'ParallelReductionFactorInPlane': 2}
    >>> get_ees(meta, in_file='epi.nii.gz')
    0.00041602630141921826
    """

    import nibabel as nb
    from fmriprep.interfaces.fmap import _get_pe_index

    # Use case 1: EES is defined
    ees = in_meta.get('EffectiveEchoSpacing', None)
    if ees is not None:
        return ees

    # All other cases require the parallel acc and npe (N vox in PE dir)
    acc = float(in_meta.get('ParallelReductionFactorInPlane', 1.0))
    npe = nb.load(in_file).shape[_get_pe_index(in_meta)]
    etl = npe // acc

    # Use case 2: TRT is defined
    trt = in_meta.get('TotalReadoutTime', None)
    if trt is not None:
        return trt / (etl - 1)

    # Use case 3 (philips scans)
    wfs = in_meta.get('WaterFatShift', None)
    if wfs is not None:
        fstrength = in_meta['MagneticFieldStrength']
        wfd_ppm = 3.4  # water-fat diff in ppm
        g_ratio_mhz_t = 42.57  # gyromagnetic ratio for proton (1H) in MHz/T
        wfs_hz = fstrength * wfd_ppm * g_ratio_mhz_t
        return wfs / (wfs_hz * etl)

    raise ValueError('Unknown effective echo-spacing specification')


def get_trt(in_meta, in_file=None):
    """
    Calculate the *total readout time* for an input
    :abbr:`EPI (echo-planar imaging)` scan.
    There are several procedures to calculate the total
    readout time. The basic one is that a ``TotalReadoutTime``
    field is set in the JSON sidecar. The following examples
    use an ``'epi.nii.gz'`` file-stub which has 90 pixels in the
    j-axis encoding direction.
    >>> meta = {'TotalReadoutTime': 0.02596}
    >>> get_trt(meta)
    0.02596
    If the *effective echo spacing* :math:`t_\\text{ees}`
    (``EffectiveEchoSpacing`` BIDS field) is provided, then the
    total readout time can be calculated reading the number
    of voxels along the readout direction :math:`T_\\text{ro}`
    and the parallel acceleration factor of the EPI :math:`f_\\text{acc}`.
      .. math ::
          T_\\text{ro} = t_\\text{ees} \\, (N_\\text{PE} / f_\\text{acc} - 1)
    >>> meta = {'EffectiveEchoSpacing': 0.00059,
    ...         'PhaseEncodingDirection': 'j-',
    ...         'ParallelReductionFactorInPlane': 2}
    >>> get_trt(meta, in_file='epi.nii.gz')
    0.02596
    Some vendors, like Philips, store different parameter names:
    >>> meta = {'WaterFatShift': 8.129,
    ...         'MagneticFieldStrength': 3,
    ...         'PhaseEncodingDirection': 'j-',
    ...         'ParallelReductionFactorInPlane': 2}
    >>> get_trt(meta, in_file='epi.nii.gz')
    0.018721183563864822
    """

    # Use case 1: TRT is defined
    trt = in_meta.get('TotalReadoutTime', None)
    if trt is not None:
        return trt

    # All other cases require the parallel acc and npe (N vox in PE dir)
    acc = float(in_meta.get('ParallelReductionFactorInPlane', 1.0))
    npe = nb.load(in_file).shape[_get_pe_index(in_meta)]
    etl = npe // acc

    # Use case 2: TRT is defined
    ees = in_meta.get('EffectiveEchoSpacing', None)
    if ees is not None:
        return ees * (etl - 1)

    # Use case 3 (philips scans)
    wfs = in_meta.get('WaterFatShift', None)
    if wfs is not None:
        fstrength = in_meta['MagneticFieldStrength']
        wfd_ppm = 3.4  # water-fat diff in ppm
        g_ratio_mhz_t = 42.57  # gyromagnetic ratio for proton (1H) in MHz/T
        wfs_hz = fstrength * wfd_ppm * g_ratio_mhz_t
        return wfs / wfs_hz

    raise ValueError('Unknown total-readout time specification')


def _get_pe_index(meta):
    pe = meta['PhaseEncodingDirection']
    try:
        return {'i': 0, 'j': 1, 'k': 2}[pe[0]]
    except KeyError:
        raise RuntimeError('"%s" is an invalid PE string' % pe)


def _torads(in_file, fmap_range=None, newpath=None):
    """
    Convert a field map to rad/s units
    If fmap_range is None, the range of the fieldmap
    will be automatically calculated.
    Use fmap_range=0.5 to convert from Hz to rad/s
    """
    from math import pi
    import nibabel as nb
    from nipype.utils.filemanip import fname_presuffix

    out_file = fname_presuffix(in_file, suffix='_rad', newpath=newpath)
    fmapnii = nb.load(in_file)
    fmapdata = fmapnii.get_data()

    if fmap_range is None:
        fmap_range = max(abs(fmapdata.min()), fmapdata.max())
    fmapdata = fmapdata * (pi / fmap_range)
    out_img = nb.Nifti1Image(fmapdata, fmapnii.affine, fmapnii.header)
    out_img.set_data_dtype('float32')
    out_img.to_filename(out_file)
    return out_file, fmap_range


def _tohz(in_file, range_hz, newpath=None):
    """Convert a field map to Hz units"""
    from math import pi
    import nibabel as nb
    from nipype.utils.filemanip import fname_presuffix

    out_file = fname_presuffix(in_file, suffix='_hz', newpath=newpath)
    fmapnii = nb.load(in_file)
    fmapdata = fmapnii.get_data()
    fmapdata = fmapdata * (range_hz / pi)
    out_img = nb.Nifti1Image(fmapdata, fmapnii.affine, fmapnii.header)
    out_img.set_data_dtype('float32')
    out_img.to_filename(out_file)
    return out_file


def phdiff2fmap(in_file, delta_te, newpath=None):
    r"""
    Converts the input phase-difference map into a fieldmap in Hz,
    using the eq. (1) of [Hutton2002]_:
    .. math::
        \Delta B_0 (\text{T}^{-1}) = \frac{\Delta \Theta}{2\pi\gamma \Delta\text{TE}}
    In this case, we do not take into account the gyromagnetic ratio of the
    proton (:math:`\gamma`), since it will be applied inside TOPUP:
    .. math::
        \Delta B_0 (\text{Hz}) = \frac{\Delta \Theta}{2\pi \Delta\text{TE}}
    """
    import math
    import numpy as np
    import nibabel as nb
    from nipype.utils.filemanip import fname_presuffix
    #  GYROMAG_RATIO_H_PROTON_MHZ = 42.576

    out_file = fname_presuffix(in_file, suffix='_fmap', newpath=newpath)
    image = nb.load(in_file)
    data = (image.get_data().astype(np.float32) / (2. * math.pi * delta_te))
    nii = nb.Nifti1Image(data, image.affine, image.header)
    nii.set_data_dtype(np.float32)
    nii.to_filename(out_file)
    return out_file


def _delta_te(in_values, te1=None, te2=None):
    """Read :math:`\Delta_\text{TE}` from BIDS metadata dict"""
    if isinstance(in_values, float):
        te2 = in_values
        te1 = 0.

    if isinstance(in_values, dict):
        te1 = in_values.get('EchoTime1')
        te2 = in_values.get('EchoTime2')

        if not all((te1, te2)):
            te2 = in_values.get('EchoTimeDifference')
            te1 = 0

    if isinstance(in_values, list):
        te2, te1 = in_values
        if isinstance(te1, list):
            te1 = te1[1]
        if isinstance(te2, list):
            te2 = te2[1]

    # For convienience if both are missing we should give one error about them
    if te1 is None and te2 is None:
        raise RuntimeError('EchoTime1 and EchoTime2 metadata fields not found. '
                           'Please consult the BIDS specification.')
    if te1 is None:
        raise RuntimeError(
            'EchoTime1 metadata field not found. Please consult the BIDS specification.')
    if te2 is None:
        raise RuntimeError(
            'EchoTime2 metadata field not found. Please consult the BIDS specification.')

    return abs(float(te2) - float(te1))
