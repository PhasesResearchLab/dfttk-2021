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
import socket
import pytz
from dfttk.pythelec import get_code_version


def run_task_ext(t,vasp_cmd,db_file,structure,tag,override_default_vasp_params):
    try:
        store_raw_vasprunxml = override_default_vasp_params['user_incar_settings']['store_raw_vasprunxml']
    except:
        store_raw_vasprunxml = False
    if isinstance(store_raw_vasprunxml, int): kmesh_factor = store_raw_vasprunxml
    elif store_raw_vasprunxml: kmesh_factor = 2
    else: kmesh_factor = 0
    if kmesh_factor > 1:
        t.append(nonscalc())
        t.append(RunVaspCustodian(vasp_cmd=vasp_cmd, auto_npar=">>auto_npar<<", gzip_output=False))
        t.append(InsertXMLToDb(db_file=db_file, structure=structure, 
            tag=tag, xml="vasprun.xml", kmesh_factor = kmesh_factor))
    elif kmesh_factor == 1:
        t.append(InsertXMLToDb(db_file=db_file, structure=structure, 
            tag=tag, xml="vasprun.xml", kmesh_factor = 1))


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
            for i in range(0,2):
                f.write(lines[i])
            f.write("Gamma\n")
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
    optional_params = ['metadata','structure', 'kmesh_factor']

    def run_task(self, fw_spec):
        self.xml = self.get("xml", None)
        self.kmesh_factor = self.get("kmesh_factor", 1)
        with open("KPOINTS", "r") as f:
            kpoints = f.readlines()
        with open("INCAR", "r") as f:
            incar = f.readlines()
        if self.kmesh_factor>1:
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
                        'hostname': socket.gethostname(),
                        self.xml.replace(".","_")+"_gz": bson.Binary(pickle.dumps(binxmldata)),
                       'DOSCAR_gz': bson.Binary(pickle.dumps(bindoscar)),
                       'volume': structure.volume,
                       'INCAR': incar,
                       'KPOINTS': kpoints,
                       'VASP': get_code_version(),
                       'last_updated':datetime.datetime.utcnow(),
                       'US_eastern_time':datetime.datetime.now(pytz.timezone('US/Eastern')),
                       'structure': structure.as_dict(),
                       'formula_pretty': structure.composition.reduced_formula}
            self.vasp_db.db['xmlgz'].insert_one(xml_data) 
