import bson
import pickle
import gzip
import shutil
from fireworks import explicit_serialize, FiretaskBase, FWAction
from atomate.vasp.database import VaspCalcDb
from atomate.utils.utils import load_class, env_chk
from pymatgen.core import Structure


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
            bindata = gzip.compress(bytes(xmldata))
            #with gzip.open("zen.txt.gz", "wb") as f:
            #f.write(bindata)
            self.db_file = env_chk(self.get("db_file"), fw_spec)
            self.vasp_db = VaspCalcDb.from_db_file(self.db_file, admin=True)

            structure = self.get('structure', Structure.from_file('POSCAR'))

            xml_data = {'metadata': {'tag': self.get('tag')},
                       'type': self.xml+".gz",
                       'xmldata': bson.Binary(pickle.dumps(bindata)),
                       'volume': structure.volume,
                       'last_updated':datetime.datetime.utcnow(),
                       'structure': structure.as_dict(),
                       'formula_pretty': structure.composition.reduced_formula}
            self.vasp_db.db['xmlgz'].insert_one(xml_data) 
