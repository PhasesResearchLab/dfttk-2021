#!python
# This script is used to query the mongodb
#
from pymatgen import Structure
from atomate.vasp.database import VaspCalcDb
from monty.serialization import loadfn
from fireworks.fw_config import config_to_dict
from dfttk.utils import sort_x_by_y
from warnings import warn

def get_eq_structure_by_metadata(metadata, db_file=None):
    '''
    Get the equilibrium structure by metadata

    Parameters
    ----------
        metadata: dict
            The metadata use for searching the database
        db_file: filepath-like
            The file path of db.json(the settings for mongodb file)
            if it is None or >>db_file<<, then using the db_file of fireworks configurations
    Returns
    -------
        eq_structure: pymatgen.Strucutre object
            the equilibrium structure
    '''
    if (db_file is None) or (db_file == '>>db_file<<'):
        db_file = loadfn(config_to_dict()["FWORKER_LOC"])["env"]["db_file"]
    vasp_db = VaspCalcDb.from_db_file(db_file, admin=True)
    static_items = list(vasp_db.db['tasks'].find({'metadata.tag': metadata['tag']}).sort('_id', 1))
    structure_list = [itemi['output']['structure'] for itemi in static_items]
    if structure_list:
        #not empty
        eq_structure = Structure.from_dict(structure_list[0])
        return eq_structure

def get_static_structure_by_metadata(metadata, db_file=None):
    '''
    Get the static structure by metadata

    Parameters
    ----------
        metadata: dict
            The metadata use for searching the database
        db_file: filepath-like
            The file path of db.json(the settings for mongodb file)
            if it is None or >>db_file<<, then using the db_file of fireworks configurations
    Returns
    -------
        structure_list: list
            The list of different structures.
            The structures are sorted by volume
    '''
    if (db_file is None) or (db_file == '>>db_file<<'):
        db_file = loadfn(config_to_dict()["FWORKER_LOC"])["env"]["db_file"]
    vasp_db = VaspCalcDb.from_db_file(db_file, admin=True)
    static_items = list(vasp_db.collection.find({'$and':[ {'metadata.tag': metadata['tag']}, {'adopted': True} ]}))
    #static_items = list(vasp_db.db['tasks'].find({'metadata.tag': metadata}))
    structure_list = [Structure.from_dict(itemi['output']['structure']) for itemi in static_items]
    volumes = [itemi['output']['structure']['lattice']['volume'] for itemi in static_items]
    energies = [itemi['output']['energy_per_atom'] for itemi in static_items]
    band_gap = []

    static_settings = {}
    for itemi in static_items:
        if itemi['output']['is_gap_direct']:
            band_gap.append(itemi['output']['bandgap'])
        else:
            band_gap.append(itemi['output']['direct_gap'])
        if not static_settings:
            pot = itemi['input']['pseudo_potential']['functional'].upper()
            if pot=="":
                pot = itemi['orig_inputs']['potcar']['functional'].upper()
                if pot=='Perdew-Zunger81'.upper(): pot="LDA"
            static_settings['user_potcar_functional'] = pot
            static_settings['user_incar_settings'] = itemi['input']['incar']

    structure_list = sort_x_by_y(structure_list, volumes)
    band_gap = sort_x_by_y(band_gap, volumes)
    energies = sort_x_by_y(energies, volumes)
    volumes = sorted(volumes)
    return (structure_list, energies, band_gap, static_settings, volumes)


def is_property_exist_in_db(metadata, db_file=None, collection='tasks'):
    '''
    Search the MongoDB collection by metadata
    '''
    if (db_file is None) or (db_file == '>>db_file<<'):
        db_file = loadfn(config_to_dict()["FWORKER_LOC"])["env"]["db_file"]
    if collection == 'tasks':
        return get_static_structure_by_metadata(metadata=metadata, db_file=db_file)
    else:
        #try:
        if True:
            vasp_db = VaspCalcDb.from_db_file(db_file, admin=True)
            search_items = list(vasp_db.db[collection].find({'metadata.tag': metadata['tag']},{'_id':0,'volume':1}))
            try:
                volumes = [f['volume'] for f in search_items]
            except:
                volumes = []
            if not volumes:
                try:
                    search_items = list(vasp_db.db[collection].find({'metadata.tag': metadata['tag']},{'_id':0,'initial_structure':1}))
                    volumes = [f['initial_structure']['lattice']['volume'] for f in search_items]
                except:
                    volumes = []
        #except:
        #    volumes = []
        return volumes


def remove_data_by_metadata(tag, db_file=None, rem_mode='vol', forcedelete=False):
    '''
    rem_mode: str/list
        allvol: remove all volume collection
        vol: remove volume collection except bandstructure and dos
        property: other collections except volume
        all: all
        aeccar: remove all aeccar related
        chgcar: remove chgcar
        locpot: remove locpot

    '''
    if (db_file is None) or (db_file == '>>db_file<<'):
        db_file = loadfn(config_to_dict()["FWORKER_LOC"])["env"]["db_file"]
    vasp_db = VaspCalcDb.from_db_file(db_file, admin=True)
    metadata = {'tag': tag}

    VOL1_COLLECTION = ['aeccar0', 'aeccar1', 'aeccar2', 'chgcar', 'locpot', 'elfcar']
    VOL2_COLLECTION = ['bandstructure', 'dos']
    VOL_COLLECTION = VOL1_COLLECTION + VOL2_COLLECTION
    OTHER_COLLECTION = ['borncharge', 'phonon', 'qha', 'qha_phonon', 'relax',
                        'relax_scheme', 'relaxations', 'tasks']
    if isinstance(rem_mode, str):
        rem_mode = rem_mode.lower()
        if rem_mode == 'all':
            collections = VOL_COLLECTION + OTHER_COLLECTION
        elif rem_mode == 'allvol':
            collections = VOL_COLLECTION
        elif rem_mode == 'vol':
            collections = VOL1_COLLECTION
        elif rem_mode == 'property':
            collections = OTHER_COLLECTION
        elif rem_mode == 'aeccar':
            collections = ['aeccar0', 'aeccar1', 'aeccar2']
        else:
            collections = [rem_mode]
    elif isinstance(rem_mode, list):
        collections = rem_mode
    else:
        raise ValueError('Unsupported remove mode, please provide a str or list')

    flag_remove = False
    if forcedelete:
        flag_remove = True
    else:
        if tag:
            if input('Are you sure? This will remove all data in {} collection with metadata.tag={}. (Y/N)'.format(collections, tag))[0].upper() == 'Y':
                flag_remove = True
        else:
            #tag is None, which means remove the collection
            if input('Are you sure? This will remove the {} collections. (Y/N)'.format(collections))[0].upper() == 'Y':
                flag_remove = True
    if flag_remove:
        for collectioni in collections:
            if collectioni in VOL_COLLECTION:
                collectioni_file = collectioni + '_fs.files'
                collectioni_chunk = collectioni + '_fs.chunks'
                #It has files and chunks
                if tag:
                    static_items = list(vasp_db.db['tasks'].find({'metadata.tag': tag}))
                    for itemi in static_items:
                        task_id = itemi['task_id']
                        files_id = list(vasp_db.db[collectioni_file].find({'metadata.task_id': task_id}))
                        if files_id:
                            vasp_db.db[collectioni_chunk].remove({'files_id': files_id[0]['_id']})
                            vasp_db.db[collectioni_file].remove({'metadata.task_id': task_id})
                            print('The volume data with metadata.tag={} in {} collection is removed'.format(tag, collectioni))
                else:
                    vasp_db.db[collectioni_chunk].remove()
                    vasp_db.db[collectioni_file].remove()
                    print('The data in {} collection is removed'.format(collectioni))
            else:
                if tag:
                    vasp_db.db[collectioni].remove({'metadata.tag': tag})
                    print('The data with metadata.tag={} in {} collection is removed'.format(tag, collectioni))
                else:
                    vasp_db.db[collectioni].remove()
                    print('The data in {} collection is removed'.format(collectioni))
