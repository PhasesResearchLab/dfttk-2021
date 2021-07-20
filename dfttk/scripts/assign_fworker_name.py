

from atomate.utils.utils import get_meta_from_structure, get_fws_and_tasks
def set_queue_options(
    original_wf,
    walltime=None,
    time_min=None,
    qos=None,
    pmem=None,
    fw_name_constraint=None,
    task_name_constraint=None,
):
    """
    Modify queue submission parameters of Fireworks in a Workflow.
    This powerup overrides paramters in the qadapter file by setting values in
    the 'queueadapter' key of a Firework spec. For example, the walltime
    requested from a queue can be modified on a per-workflow basis.
    Args:
        original_wf (Workflow):
        walltime (str): Total walltime to request for the job in HH:MM:SS
            format e.g., "00:10:00" for 10 minutes.
        time_min (str): Minimum walltime to request in HH:MM:SS format.
            Specifying both `walltime` and `time_min` can improve throughput on
            some queues.
        qos (str): QoS level to request. Typical examples include "regular",
            "flex", and "scavenger". For Cori KNL "flex" QoS, it is necessary
            to specify a `time_min` of no more than 2 hours.
        fw_name_constraint (str): name of the Fireworks to be tagged (all if
            None is passed)
        task_name_constraint (str): name of the Firetasks to be tagged (e.g.
            None or 'RunVasp')
    Returns:
        Workflow: workflow with modified queue options
    """
    qsettings = {}
    if walltime:
        qsettings.update({"walltime": walltime})
    if time_min:
        qsettings.update({"time_min": time_min})
    if qos:
        qsettings.update({"qos": qos})
    if pmem:
        qsettings.update({"pmem": pmem})

    idx_list = get_fws_and_tasks(
        original_wf,
        fw_name_constraint=fw_name_constraint,
        task_name_constraint=task_name_constraint,
    )

    for idx_fw, idx_t in idx_list:
        original_wf.fws[idx_fw].spec.update({"_queueadapter": qsettings})

    return original_wf


def set_execution_options(
    original_wf,
    fworker_name=None,
    category=None,
    fw_name_constraint=None,
    task_name_constraint=None,
):
    """
    set _fworker spec of Fireworker(s) of a Workflow. It can be used to specify
    a queue; e.g. run large-memory jobs on a separate queue.

    Args:
        original_wf (Workflow):
        fworker_name (str): user-defined tag to be added under fw.spec._fworker
            e.g. "large memory", "big", etc
        category (str): category of FWorker that should pul job
        fw_name_constraint (str): name of the Fireworks to be tagged (all if
            None is passed)
        task_name_constraint (str): name of the Firetasks to be tagged (e.g.
            None or 'RunVasp')

    Returns:
        Workflow: modified workflow with specified Fireworkers tagged
    """
    idx_list = get_fws_and_tasks(
        original_wf,
        fw_name_constraint=fw_name_constraint,
        task_name_constraint=task_name_constraint,
    )

    for idx_fw, idx_t in idx_list:
        if fworker_name:
            original_wf.fws[idx_fw].spec["_fworker"] = fworker_name
        if category:
            original_wf.fws[idx_fw].spec["_category"] = category
    return original_wf


import atomate.vasp.powerups as powerups
def Customizing_Workflows(original_wf, user_settings={}):
    """
    set _preserve_fworker spec of Fireworker(s) of a Workflow. Can be used to
    pin a workflow to the first fworker it is run with. Very useful when running
    on multiple machines that can't share files.
    """
    original_wf = powerups.preserve_fworker(original_wf)
    powerups_options = user_settings.get('powerups', {}) 
    if len(powerups_options) == 0:
        if 'user_incar_settings' in user_settings:
            powerups_options = user_settings['user_incar_settings'].get('powerups', {})

    if 'set_execution_options' in powerups_options:
        execution_options = powerups_options['set_execution_options']
        original_wf = set_execution_options(original_wf, 
            fworker_name=execution_options.get("fworker_name", None),
            category=execution_options.get("category", None),
            )

    if 'set_queue_options' in powerups_options:
        queue_options = powerups_options['set_queue_options']
        original_wf = set_queue_options(original_wf, pmem=queue_options.get("pmem", None))

    return original_wf
