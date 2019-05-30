# -*- coding: utf-8 -*-

import sys
import os
from copy import deepcopy

from nipype import __version__ as nipype_ver
from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from nipype.interfaces.io import BIDSDataGrabber

from niworkflows.interfaces.bids import BIDSInfo, BIDSFreeSurferDir
from niworkflows.utils.bids import collect_data


def init_dmriprep_wf(layout, subject_list, output_dir,
                     ignore_fmaps, eddy_niter, slice_outlier_threshold):

    dmriprep_wf = pe.Workflow(name='dmriprep_wf')
    dmriprep_wf.base_dir = output_dir

    for subject_id in subject_list:
        single_subject_wf = init_single_subject_wf(
            layout=layout,
            subject_id=subject_id,
            name="single_subject_" + subject_id + "_wf",
            output_dir=output_dir,
            ignore_fmaps=ignore_fmaps,
            eddy_niter=eddy_niter,
            slice_outlier_threshold=slice_outlier_threshold
        )

        if freesurfer:
            dmriprep_wf.connect(fsdir, 'subjects_dir',
                                single_subject_wf, 'inputnode.subjects_dir')
        else:
            dmriprep_wf.add_nodes([single_subject_wf])

    return dmriprep_wf


def init_single_subject_wf(layout, subject_id, output_dir,
                           ignore_fmaps, eddy_niter, slice_outlier_threshold):

    subject_data = collect_data(layout, subject_id, task_id, echo_idx)[0]

    # Make sure we always go through these two checks
    if subject_data['dwi'] == []:
        raise Exception("No diffusion images found for participant {}. "
                        "All workflows require diffusion images.".format(subject_id))

    if not subject_data['t1w']:
        raise Exception("No T1w images found for participant {}. "
                        "All workflows require T1w images.".format(subject_id))

    workflow = pe.Workflow(name=name)

    inputnode = pe.Node(niu.IdentityInterface(fields=['subjects_dir']),
                        name='inputnode')

    bidssrc = pe.Node(BIDSDataGrabber(subject_data=subject_data, anat_only=anat_only),
                      name='bidssrc')

    bids_info = pe.Node(BIDSInfo(
        bids_dir=layout.root, bids_validate=False), name='bids_info')

    summary = pe.Node(SubjectSummary(output_spaces=output_spaces, template=template),
                      name='summary', run_without_submitting=True)

    about = pe.Node(AboutSummary(version=__version__,
                                 command=' '.join(sys.argv)),
                    name='about', run_without_submitting=True)

    ds_report_summary = pe.Node(
        DerivativesDataSink(base_directory=reportlets_dir,
                            suffix='summary'),
        name='ds_report_summary', run_without_submitting=True)

    ds_report_about = pe.Node(
        DerivativesDataSink(base_directory=reportlets_dir,
                            suffix='about'),
        name='ds_report_about', run_without_submitting=True)

    workflow.connect([
        (inputnode, anat_preproc_wf, [('subjects_dir', 'inputnode.subjects_dir')]),
        (bidssrc, bids_info, [(('t1w', fix_multi_T1w_source_name), 'in_file')]),
        (inputnode, summary, [('subjects_dir', 'subjects_dir')]),
        (bidssrc, summary, [('t1w', 't1w'),
                            ('t2w', 't2w'),
                            ('bold', 'bold')]),
        (bids_info, summary, [('subject', 'subject_id')]),
        (bids_info, anat_preproc_wf, [(('subject', _prefix), 'inputnode.subject_id')]),
        (bidssrc, anat_preproc_wf, [('t1w', 'inputnode.t1w'),
                                    ('t2w', 'inputnode.t2w'),
                                    ('roi', 'inputnode.roi'),
                                    ('flair', 'inputnode.flair')]),
        (bidssrc, ds_report_summary, [(('t1w', fix_multi_T1w_source_name), 'source_file')]),
        (summary, ds_report_summary, [('out_report', 'in_file')]),
        (bidssrc, ds_report_about, [(('t1w', fix_multi_T1w_source_name), 'source_file')]),
        (about, ds_report_about, [('out_report', 'in_file')]),
    ])

    # Overwrite ``out_path_base`` of smriprep's DataSinks
    for node in workflow.list_node_names():
        if node.split('.')[-1].startswith('ds_'):
            workflow.get_node(node).interface.out_path_base = 'dmriprep'

    if anat_only:
        return workflow

    for bold_file in subject_data['bold']:
        func_preproc_wf = init_func_preproc_wf(bold_file=bold_file,
                                               layout=layout,
                                               ignore=ignore,
                                               freesurfer=freesurfer,
                                               use_bbr=use_bbr,
                                               t2s_coreg=t2s_coreg,
                                               bold2t1w_dof=bold2t1w_dof,
                                               reportlets_dir=reportlets_dir,
                                               output_spaces=output_spaces,
                                               template=template,
                                               medial_surface_nan=medial_surface_nan,
                                               cifti_output=cifti_output,
                                               output_dir=output_dir,
                                               omp_nthreads=omp_nthreads,
                                               low_mem=low_mem,
                                               fmap_bspline=fmap_bspline,
                                               fmap_demean=fmap_demean,
                                               use_syn=use_syn,
                                               force_syn=force_syn,
                                               debug=debug,
                                               template_out_grid=template_out_grid,
                                               use_aroma=use_aroma,
                                               aroma_melodic_dim=aroma_melodic_dim,
                                               err_on_aroma_warn=err_on_aroma_warn,
                                               num_bold=len(subject_data['bold']))

        workflow.connect([
            (anat_preproc_wf, func_preproc_wf,
             [('outputnode.t1_preproc', 'inputnode.t1_preproc'),
              ('outputnode.t1_brain', 'inputnode.t1_brain'),
              ('outputnode.t1_mask', 'inputnode.t1_mask'),
              ('outputnode.t1_seg', 'inputnode.t1_seg'),
              ('outputnode.t1_aseg', 'inputnode.t1_aseg'),
              ('outputnode.t1_aparc', 'inputnode.t1_aparc'),
              ('outputnode.t1_tpms', 'inputnode.t1_tpms'),
              ('outputnode.t1_2_mni_forward_transform', 'inputnode.t1_2_mni_forward_transform'),
              ('outputnode.t1_2_mni_reverse_transform', 'inputnode.t1_2_mni_reverse_transform'),
              # Undefined if --no-freesurfer, but this is safe
              ('outputnode.subjects_dir', 'inputnode.subjects_dir'),
              ('outputnode.subject_id', 'inputnode.subject_id'),
              ('outputnode.t1_2_fsnative_forward_transform',
               'inputnode.t1_2_fsnative_forward_transform'),
              ('outputnode.t1_2_fsnative_reverse_transform',
               'inputnode.t1_2_fsnative_reverse_transform')]),
        ])

    return workflow


def _prefix(subid):
    if subid.startswith('sub-'):
        return subid
    return '-'.join(('sub', subid))
