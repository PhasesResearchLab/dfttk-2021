
import os
import datetime
import subprocess
from custodian.custodian import ErrorHandler, Job
from custodian.ansible.interpreter import Modder
from custodian.ansible.actions import FileActions


class ATATWalltimeHandler(ErrorHandler):
    """
    Check if a run is nearing the walltime. If so, write a file called ``stop``
    which will stop ATAT. You can specify the walltime either in the init,
    which is unfortunately necessary for SGE and SLURM systems. If you happen
    to be running on a PBS system and the PBS_WALLTIME variable is in the run
    environment, the wall time will be automatically determined if not set.
    """
    is_monitor = True

    # The WalltimeHandler should not terminate as we want ATAT to terminate
    # itself naturally with the ``stop`` file.
    is_terminating = False

    # This handler will be unrecoverable, but custodian shouldn't raise an
    # error
    raises_runtime_error = False

    def __init__(self, wall_time=None, buffer_time=3600):
        """
        Initializes the handler with a buffer time.
        Args:
            wall_time (int): Total walltime in seconds. If this is None and
                the job is running on a PBS system, the handler will attempt to
                determine the walltime from the PBS_WALLTIME environment
                variable. If the wall time cannot be determined or is not
                set, this handler will have no effect. Same for SLURM_TIMELIMIT,
                but this must be set by hand on some systems.
            buffer_time (int): The minimum amount of buffer time in secs at the
                end that the ``stop`` file will be written. The ``stop`` file is
                written when the time remaining is < the buffer time. Defaults to
                1 hour (to be safe).
        """
        if wall_time is not None:
            self.wall_time = wall_time
        elif "PBS_WALLTIME" in os.environ:
            self.wall_time = int(os.environ["PBS_WALLTIME"])
        elif "SBATCH_TIMELIMIT" in os.environ:
            self.wall_time = int(os.environ["SBATCH_TIMELIMIT"])
        else:
            self.wall_time = None
        self.buffer_time = buffer_time
        # Sets CUSTODIAN_WALLTIME_START as the start time to use for
        # future jobs in the same batch environment.  Can also be
        # set manually be the user in the batch environment.
        if "CUSTODIAN_WALLTIME_START" in os.environ:
            self.start_time = datetime.datetime.strptime(
                os.environ["CUSTODIAN_WALLTIME_START"], "%a %b %d %H:%M:%S %Z %Y")
        else:
            self.start_time = datetime.datetime.now()
            os.environ["CUSTODIAN_WALLTIME_START"] = datetime.datetime.strftime(
                self.start_time, "%a %b %d %H:%M:%S UTC %Y")

        self.prev_check_time = self.start_time

    def check(self):
        if self.wall_time:
            run_time = datetime.datetime.now() - self.start_time
            total_secs = run_time.total_seconds()

            time_left = self.wall_time - total_secs
            if time_left < self.buffer_time:
                return True

        return False

    def correct(self):
        # Write ``stop`` file
        actions = [{"file": "stop",
                    "action": {"_file_create": {'content': ''}}}]

        m = Modder(actions=[FileActions])
        for a in actions:
            m.modify(a["action"], a["file"])
        return {"errors": ["Walltime reached"], "actions": None}


class ATATInfDetJob(Job):
    def __init__(self, continuation=False, inflection=False):
        self.name = 'ATAT Inflection Detection Job'
        self.continuation = continuation
        self.inflection = inflection

    def run(self):
        inf_det_str = '-id' if self.inflection else ''
        run_args = ['robustrelax_vasp', inf_det_str, '-c', '0.05', '-rc', '""', '-vc', '""']
        if self.continuation:
            run_args.append(['-cip'])
        with open('ATAT.robustrelax.out', 'w') as f_std, \
                open('ATAT.robustrelax.err', "w", buffering=1) as f_err:
            # we have to use Popen and return in order for Custodian to be able to monitor the process in real time
            return subprocess.Popen(run_args, stdout=f_std, stderr=f_err)
