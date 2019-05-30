# -*- coding: utf-8 -*-
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""
Orchestrating the diffusion preprocessing workflow
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: init_dwi_preproc_wf
.. autofunction:: init_dwi_derivatives_wf
"""

from nipype.pipeline import engine as pe


def init_dwi_preproc_wf(dwi_file, ignore, use_syn, layout=None, num_dwi=1):
    ref_file = dwi_file

    metadata = layout.get_metadata(ref_file)

    # Find fieldmaps. Options: (epi|fieldmap|phasediff|phase1|phase2|syn)
    fmaps = []
    if 'fieldmaps' not in ignore:
        fmaps = layout.get_fieldmap(ref_file, return_list=True)
        for fmap in fmaps:
            fmap['metadata'] = layout.get_metadata(fmap[fmap['suffix']])

    # Run Syn if in the absence of fieldmap correction
    if (use_syn and not fmaps):
        fmaps.append({'suffix': 'syn'})

    # Build workflow
    workflow = pe.Workflow(name=wf_name)

    pe_direction=metadata.get("PhaseEncodingDirection")
    total_readout=metadata.get("TotalReadoutTime")
