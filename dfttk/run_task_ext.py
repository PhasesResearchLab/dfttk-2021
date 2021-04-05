import os
import bson
import pickle
import gzip
import shutil
import datetime
from fireworks import explicit_serialize, FiretaskBase, FWAction
from atomate.vasp.database import VaspCalcDb
from atomate.utils.utils import load_class, env_chk
from atomate.vasp.firetasks.run_calc import RunVaspCustodian
from pymatgen.core import Structure
import dfttk.scripts.user_SETTINGS as user_SETTINGS

from monty.serialization import loadfn, dumpfn
if os.path.exists('SETTINGS.yaml'): #treat settings in 'SETTINGS.yaml' as globally accessible
    user_SETTINGS.user_settings=loadfn('SETTINGS.yaml')
    print("eeeeeeeeeeeee", user_SETTINGS.user_settings)
    global_user_SETTINGS = user_SETTINGS.user_settings


def run_task_ext(t,vasp_cmd,db_file,structure,tag):
    print(user_SETTINGS.user_settings)
    global global_user_SETTINGS
    print("lllllllll", global_user_SETTINGS)
    #if user_SETTINGS.user_settings.get('store_raw_vasprunxml', False):
    if True:
        t.append(nonscalc())
        t.append(RunVaspCustodian(vasp_cmd=vasp_cmd, auto_npar=">>auto_npar<<", gzip_output=False))
        t.append(InsertXMLToDb(db_file=db_file, structure=structure, 
            tag=tag, xml="vasprun.xml"))


@explicit_serialize
class nonscalc(FiretaskBase):
    '''
    nonselfconsistent calculation using denser k-mesh
    '''
    def run_task(self, fw_spec):
        shutil.copyfile("INCAR","INCAR.Static")
        with open("INCAR", "r") as f:
            lines = f.readlines()
        with open("INCAR", "w") as f:
            for line in lines:
                if line.lower().startswith("icharg"):
                    f.write('ICHARG=11\n')
                elif line.lower().startswith("lorbit"):
                    continue
                else:
                    f.write(line)

        shutil.copyfile("KPOINTS","KPOINTS.Static")
        with open("KPOINTS", "r") as f:
            lines = f.readlines()
        with open("KPOINTS", "w") as f:
            for i in range(0,3):
                f.write(lines[i])
            mesh = [int(x) for x in lines[3].split(" ") if x!=""]
            for i in range(0,3):
                f.write(' {}'.format(mesh[i]*2))
            f.write('\n')            



@explicit_serialize
class InsertXMLToDb(FiretaskBase):
    '''
    Store the CheckSymmetry result to MongoDB, the stored collection is named as 'relaxations'
    '''
    required_params = ["xml", "db_file", "tag"]
    optional_params = ['metadata','structure']

    def run_task(self, fw_spec):
        self.xml = self.get("xml", None)
        shutil.copyfile("INCAR","INCAR.nscf")
        shutil.copyfile("KPOINTS","KPOINTS.nscf")
        shutil.copyfile("INCAR.Static","INCAR")
        shutil.copyfile("KPOINTS.Static","KPOINTS")
        if self.xml is not None:
            with open (self.xml, 'rb') as f:
                xmldata = f.read()
            binxmldata = gzip.compress(bytes(xmldata))
            with open ("DOSCAR", 'rb') as f:
                doscar = f.read()
            bindoscar = gzip.compress(bytes(doscar))
            #with gzip.open("zen.txt.gz", "wb") as f:
            #f.write(bindata)
            self.db_file = env_chk(self.get("db_file"), fw_spec)
            self.vasp_db = VaspCalcDb.from_db_file(self.db_file, admin=True)

            structure = self.get('structure', Structure.from_file('POSCAR'))

            xml_data = {'metadata': {'tag': self.get('tag')},
                        self.xml+".gz": bson.Binary(pickle.dumps(binxmldata)),
                       'DOSCAR.gz': bson.Binary(pickle.dumps(bindoscar)),
                       'volume': structure.volume,
                       'last_updated':datetime.datetime.utcnow(),
                       'structure': structure.as_dict(),
                       'formula_pretty': structure.composition.reduced_formula}
            self.vasp_db.db['xmlgz'].insert_one(xml_data) 
