from fireworks import Workflow, Firework
from atomate.utils.utils import get_meta_from_structure, get_fws_and_tasks
import sys


def add_priority(original_wf, root_priority, child_priority=None):
    """
    Adds priority to a workflow

    Args:
        original_wf (Workflow): original WF
        root_priority (int): priority of first (root) job(s)
        child_priority(int): priority of all child jobs. Defaults to
            root_priority

    Returns:
       Workflow: priority-decorated workflow
    """
    child_priority = child_priority or root_priority
    root_fw_ids = original_wf.root_fw_ids
    for fw in original_wf.fws:
        if fw.fw_id in root_fw_ids:
            fw.spec["_priority"] = root_priority
        else:
            fw.spec["_priority"] = child_priority
    return original_wf


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
def get_powerups(wfs):
    """
    get user powerups setting.
    """
    if isinstance(wfs, list) :
        for wflow in wfs:
            try:
                powerups_options = get_powerups_spec(wflow)
            except:
                try:
                    powerups_options = get_powerups_wf(wflow)
                except:
                    powerups_options = {}
            if len(powerups_options)>0: break   
    else:
        try:
            powerups_options = get_powerups_spec(wfs)
        except:
            try:
                powerups_options = get_powerups_wf(wfs)
            except:
                powerups_options = {}
    return powerups_options


debug = False
from collections.abc import Iterable
def get_powerups_wf(original_wf):
    """
    get user powerups setting.
    """
    idx_list = get_fws_and_tasks(original_wf)
    for idx_fw, idx_t in idx_list:
        f0 = original_wf.fws[idx_fw].tasks[idx_t]
        if not isinstance(f0, Iterable) or isinstance(f0, str) : continue
        for k0 in f0:
            if debug: print("level 0", k0, type(f0))
            if k0=='powerups' : 
                if debug: print("level 0", f0[k0])
                return f0[k0]
            else:
                try:
                    f1 = f0[k0]
                except:
                    f1 = k0
                if not isinstance(f1, Iterable) or isinstance(f1, str) : continue
                for k1 in f1:
                    if debug: print("level 1", k1, type(f1))
                    if str(k1)=='powerups' : 
                        if debug: print("level 1", f1[k1])
                        return f1[k1]
                    else:
                        try:
                            f2 = f1[k1]
                        except:
                            f2 = k1
                        if not isinstance(f2, Iterable) or isinstance(f2, str) : continue
                        for k2 in f2:
                            if debug: print("level 2", k2, type(f2))
                            if str(k2)=='powerups' : 
                                if debug: print("level 2", f2[k2])
                                return f2[k2]
                            else:
                                try:
                                    f3 = f2[k2]
                                except:
                                    f3=k2
                                if not isinstance(f3, Iterable) or isinstance(f3, str) : continue
                                for k3 in f3:
                                    if debug: print("level 3", k3, type(f3))
                                    if str(k3)=='powerups' : 
                                        if debug: print(type(f0),type(f1),type(f2),type(f3))
                                        if debug: print("level 3", f3[k3])
                                        return f3[k3]
                                    else:
                                        try:
                                            f4 = f3[k3]
                                        except:
                                            f4=k3
                                        if not isinstance(f4, Iterable) or isinstance(f4, str) : continue
                                        for k4 in f4:
                                            if debug: print("level 4", k4, type(f4))
                                            if str(k4)=='powerups' : 
                                                if debug: print("level 4", f4[k4])
                                                return f4[k4]                                        
    return {}


def get_powerups_spec(original_wf):
        """
        get user powerups setting.
        """
        if debug: print("level -2", type(original_wf), original_wf)
        f0 = original_wf.spec
        print ("Here 1", f0)
        if debug: print("level -1", f0, type(f0))
        for k0 in f0:
            print ("Here 2", k0)
            if debug: print("level 0", k0, type(f0))
            if k0=='powerups' : 
                if debug: print("level 0", f0[k0])
                return f0[k0]
            else:
                try:
                    f1 = f0[k0]
                except:
                    f1 = k0
                if not isinstance(f1, Iterable) or isinstance(f1, str) : continue
                for k1 in f1:
                    if debug: print("level 1", k1, type(f1))
                    if str(k1)=='powerups' : 
                        if debug: print("level 1", f1[k1])
                        return f1[k1]
                    else:
                        try:
                            f2 = f1[k1]
                        except:
                            f2 = k1
                        if not isinstance(f2, Iterable) or isinstance(f2, str) : continue
                        for k2 in f2:
                            if debug: print("level 2", k2, type(f2))
                            if str(k2)=='powerups' : 
                                if debug: print("level 2", f2[k2])
                                return f2[k2]
                            else:
                                try:
                                    f3 = f2[k2]
                                except:
                                    f3=k2
                                if not isinstance(f3, Iterable) or isinstance(f3, str) : continue
                                for k3 in f3:
                                    if debug: print("level 3", k3, type(f3))
                                    if str(k3)=='powerups' : 
                                        if debug: print(type(f0),type(f1),type(f2),type(f3))
                                        if debug: print("level 3", f3[k3])
                                        return f3[k3]
                                    else:
                                        try:
                                            f4 = f3[k3]
                                        except:
                                            f4=k3
                                        if not isinstance(f4, Iterable) or isinstance(f4, str) : continue
                                        for k4 in f4:
                                            if debug: print("level 4", k4, type(f4))
                                            if str(k4)=='powerups' : 
                                                if debug: print("level 4", f4[k4])
                                                return f4[k4]                                        
        return {}


def get_powerups_options(override_default_vasp_params):
    if 'user_incar_settings' in override_default_vasp_params:
        if 'powerups' in override_default_vasp_params['user_incar_settings']:
            return override_default_vasp_params['user_incar_settings']['powerups']
    return None


def Customizing_Workflows(wfs, powerups_options=None):
    if powerups_options is None: 
        powerups_options = get_powerups(wfs)
    if powerups_options is None: return wfs

    if isinstance(wfs, list) :
        _wfs = []
        for wflow in wfs:
            revised_wflow = Customizing_Workflows_wf(wflow,powerups_options=powerups_options)
            _wfs.append(revised_wflow)
        return _wfs
    else:
        try:
            revised_wflow = Customizing_Workflows_wf(wfs,powerups_options=powerups_options)
            return revised_wflow
        except:
            print ("powerups_options = ", powerups_options)
            print("***************WARNING! not a workflow",wfs)
            return wfs


def Customizing_Workflows_wf(original_wf, powerups_options=None):
    if powerups_options is None: 
        try:
            powerups_options = get_powerups_spec(original_wf)
        except:
            try:
                powerups_options = get_powerups_wf(original_wf)
            except:
                powerups_options = {}
    """
    set _preserve_fworker spec of Fireworker(s) of a Workflow. Can be used to
    pin a workflow to the first fworker it is run with. Very useful when running
    on multiple machines that can't share files.
    """

    if 'set_execution_options' in powerups_options:
        execution_options = powerups_options['set_execution_options']
        try:
            original_wf = set_execution_options(original_wf, 
                fworker_name=execution_options.get("fworker_name", None),
                category=execution_options.get("category", None),
                )
            if 'preserve_fworker' in execution_options:
                if execution_options['preserve_fworker']:
                    original_wf = powerups.preserve_fworker(original_wf)
        except:
            if execution_options.get("fworker_name", None):
                original_wf.spec["_fworker"] = execution_options.get("fworker_name", None)
            if execution_options.get("category", None):
                original_wf.spec["_category"] = execution_options.get("category", None)
            if execution_options.get("preserve_fworker", None):
                original_wf.spec["_preserve_fworker"] = execution_options.get("preserve_fworker", None)
        if 'priority' in execution_options:
            add_priority(original_wf, execution_options.get("priority", 0))

    if 'set_queue_options' in powerups_options:
        queue_options = powerups_options['set_queue_options']
        original_wf = set_queue_options(original_wf, pmem=queue_options.get("pmem", None))

    return original_wf
