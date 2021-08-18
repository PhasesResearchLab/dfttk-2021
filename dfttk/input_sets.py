import six
import os
import subprocess
import copy

from pymatgen.io.vasp.inputs import Kpoints, Incar
from pymatgen.io.vasp.sets import DictSet, get_vasprun_outcar, get_structure_from_prev_run, _load_yaml_config


'''
# We set potentials back to the VASP recommended, since we primarily don't work on oxides.
POTCAR_UPDATES = {
        'Be': 'Be',  # 2 electrons, default was Be_sv (4 electrons)
        'Cu': 'Cu',  # 11 electrons, default was Cu_pv (17 electrons)
        'Fe': 'Fe',  # 8 electrons, default was Fe_pv (14 electrons)
        'Mg': 'Mg',  # 2 electrons, default was Mg_pv (8 electrons)
        'Ni': 'Ni',  # 10 electrons, default was Ni_pv (16 electrons)
        'Mo': 'Mo_sv',  # 14 electrons, default was Mo_pv (12 electrons)
        'Nb': 'Nb_sv',  # 13 electrons, default was Nb_pv (11 electrons)
        'Os': 'Os',  # 8 electrons, default was Os_pv (14 electrons)
        'Re': 'Re',  # 7 electrons, default was Re_pv (13 electrons)
        'Ti': 'Ti_sv',  # 12 electrons, default Ti_pv (10 electrons)
        'V': 'V_sv',  # 13 electrons, default V_pv (11 electrons)
    }
'''


#Reset the POTCAR
POTCAR_UPDATES = {
        'Cu': 'Cu',  # 11 electrons, default was Cu_pv (17 electrons)
        'Mo': 'Mo_sv',  # 14 electrons, default was Mo_pv (12 electrons)
        'Nb': 'Nb_sv',  # 13 electrons, default was Nb_pv (11 electrons)
        'Ti': 'Ti_sv',  # 12 electrons, default Ti_pv (10 electrons)
        'V': 'V_sv',  # 13 electrons, default V_pv (11 electrons)
}


def magnetic_check(structure):
    '''
    If the structure contain any magnetic elements, return True, otherwise return False
    magnetic elements:
        V(23)-Ni(28), Ru(44)-Pd(46), Ce(58)-Lu(71), Pa(91)-118
    Suggested by Yi Wang, Aug. 24, 2020
    '''
    magnetic_elements = list(range(23, 29))
    magnetic_elements.extend(list(range(44, 47)))
    magnetic_elements.extend(list(range(58, 72)))
    magnetic_elements.extend(list(range(91, 119)))
    return any(ele.Z in magnetic_elements for ele in structure.species)


def metal_check(structure):
    '''
    If the structure contain any metal elements, return True, otherwise return False
    metal elements:
        Li(3)-Be(4), Na(11)-Al(13), K(19)-Ga(31), Rb(37)-Sn(50), Cs(55)-At(85), Ra(87)-118
    Suggested by Yi Wang, Aug. 24, 2020
    '''
    metal_elements = list(range(3, 5))
    metal_elements.extend(list(range(11, 14)))
    metal_elements.extend(list(range(19, 32)))
    metal_elements.extend(list(range(37, 51)))
    metal_elements.extend(list(range(55, 86)))
    metal_elements.extend(list(range(87, 119)))
    return all(ele.Z in metal_elements for ele in structure.species)


class RelaxSet(DictSet):
    """
    Set for performing relaxations.

    The smearing must be set to give the correct forces.
    The default is tuned for metal relaxations.
    Kpoints have a 8000 kpoints per reciprocal atom default.
    """
    CONFIG = _load_yaml_config("MPRelaxSet")
    # we never are comparing relaxations, only using them for optimizing structures.
    CONFIG['INCAR'].pop('ENCUT')  # use the ENCUT set by PREC
    CONFIG['KPOINTS'].update({
        'grid_density': 8000,
    })
    CONFIG['KPOINTS'].pop('reciprocal_density') # to be explicit
    CONFIG['INCAR'].update({
        'EDIFF_PER_ATOM': 1e-7,
        'ISMEAR': 1,
        'SIGMA': 0.2,
        'LREAL': False,
        'PREC': 'Accurate',
        'ALGO': 'NORMAL',
        'LWAVE': False,
        'LCHARG': False,
        'ISIF': 3,
        "ICHARG": 2,
        'ENCUT': 520,
    })
    # now we reset the potentials
    CONFIG['POTCAR_FUNCTIONAL'] = 'PBE'
    CONFIG['POTCAR'].update(POTCAR_UPDATES)

    def __init__(self, structure, volume_relax=False, isif=None, **kwargs):
        """If volume relax is True, will do volume only, ISIF 7"""
        self.kwargs = copy.deepcopy(kwargs)
        self.volume_relax = volume_relax
        self.isif = isif
        uis = copy.deepcopy(self.kwargs.get('user_incar_settings', {}))
        new_config = copy.deepcopy(RelaxSet.CONFIG)
        if self.volume_relax and self.isif is not None:
            raise ValueError("isif cannot have a value while volume_relax is True.")
        if self.volume_relax:
            uis.update({'ISIF': 7})
            #uis['ISIF'] = 7
        if self.isif is not None:
            uis.update({'ISIF': self.isif})

        if 'ISPIN' not in uis:
            if magnetic_check(structure):
                uis.update({'ISPIN': 2})
            else:
                uis.update({'ISPIN': 1})

        if 'magmom' in uis:
            if 'MAGMOM' in new_config['INCAR']:
                new_config['INCAR'].pop('MAGMOM')
        elif uis['ISPIN']==1:
            if 'MAGMON' in uis.keys(): uis.pop['MAGMOM']
            if 'MAGMON' in new_config['INCAR']: new_config['INCAR'].pop['MAGMOM']

        if 'Relax_settings' in uis:
            relax = uis['Relax_settings']
            for ff in relax:
                if ff.lower()=='prec':
                    if 'ENCUT' in new_config['INCAR']:
                        new_config['INCAR'].pop('ENCUT')
                    new_config['INCAR'].update({ff:relax.get(ff)})
                elif ff=='KPAR':
                    new_config['INCAR'].update({ff:relax.get(ff)})
                elif ff=='grid_density':
                    new_config['KPOINTS'].update({ff:relax.get(ff)})
        #print (new_config)

        new_config['INCAR'].update(uis)
        pot = self.kwargs.get('user_potcar_functional', None)
        if pot:
            new_config['POTCAR_FUNCTIONAL'] = pot
        super(RelaxSet, self).__init__(structure, new_config, sort_structure=False, **self.kwargs)


class PreStaticSet(DictSet):
    """Set tuned for metal relaxations (correct smearing).
    Add `isif` parameter to the set to easily allow for overriding ISIF setting.
    Kpoints have a 6000 reciprocal density default.
    """
    CONFIG = _load_yaml_config("MPRelaxSet")
    CONFIG['KPOINTS'].update({
        'grid_density': 1000,
    })
    CONFIG['KPOINTS'].pop('reciprocal_density')  # to be explicit
    CONFIG['INCAR'].update({
        'EDIFF_PER_ATOM': 1e-5,
        'ENCUT': 520,  # MP compatibility
        'ISMEAR': 1,
        "NSW": 0,
        "IBRION": -1,
        'LREAL': False,
        'ALGO': 'NORMAL',
        # other settings from MPStaticSet
        "LAECHG": True,
        "LCHARG": False,
        "LWAVE": False,
        "LORBIT": 11,
        "LVHAR": True,
        "ICHARG": 2,
        "NEDOS": 5001,
    })
    # now we reset the potentials
    CONFIG['POTCAR_FUNCTIONAL'] = 'PBE'
    CONFIG['POTCAR'].update(POTCAR_UPDATES)

    def __init__(self, structure, **kwargs):
        # pop the old kwargs, backwards compatibility from the complex StaticSet
        old_kwargs = ['prev_incar', 'prev_kpoints', 'grid_density', 'lepsilon', 'lcalcpol']
        for k in old_kwargs:
            try:
                kwargs.pop(k)
            except KeyError:
                pass
        self.kwargs = copy.deepcopy(kwargs)
        uis = copy.deepcopy(self.kwargs.get('user_incar_settings', {}))
        new_config = copy.deepcopy(PreStaticSet.CONFIG)
        if 'ISPIN' not in uis:
            if magnetic_check(structure):
                uis.update({'ISPIN': 2})
            else:
                uis.update({'ISPIN': 1})

        if 'magmom' in uis:
            if 'MAGMOM' in new_config['INCAR']:
                new_config['INCAR'].pop('MAGMOM')
        elif uis['ISPIN']==1:
            if 'MAGMON' in uis.keys(): uis.pop['MAGMOM']
            if 'MAGMON' in new_config['INCAR']: new_config['INCAR'].pop['MAGMOM']

        new_config['INCAR'].update(uis)
        pot = self.kwargs.get('user_potcar_functional', None)
        if pot:
            new_config['POTCAR_FUNCTIONAL'] = pot
        super(PreStaticSet, self).__init__(structure, new_config, sort_structure=False, **self.kwargs)


class ForceConstantsSet(DictSet):
    """
    Set for calculating force constants calculations.

    Force constants are calculated by the finite difference method with symmetry considered.

    The smearing must be set to give the correct forces.
    The default is tuned for metals.

    Kpoints have a 8000 kpoints per reciprocal atom default.
    """
    CONFIG = _load_yaml_config("MPRelaxSet")
    # we never are comparing relaxations, only using them for optimizing structures.
    CONFIG['KPOINTS'].update({
        'grid_density': 8000,
    })
    CONFIG['KPOINTS'].pop('reciprocal_density') # to be explicit
    CONFIG['INCAR'].pop('ENCUT')  # use the ENCUT set by PREC
    CONFIG['INCAR'].update({
        'EDIFF_PER_ATOM': 1e-6,
        'ISMEAR': 1,
        'SIGMA': 0.2,
        'LREAL': False,
        'ISIF': 0,  # only calculate the forces, stress tensor is not needed
        'IBRION': 6,  # calculate force constants by finite differences with symmetry
        'POTIM': 0.015,  # displacement distance
        'NFREE': 2,  # how many displacments to do. 2 gives +POTIM and -POTIM
        'NSW': 1,  # backwards compatibility setting
        'PREC': 'Accurate',
        'ALGO': 'NORMAL',
        'SYMPREC': 1e-4,  # some supercells seem to have issues with primcel VASP algorithm
        "ICHARG": 2,
    })
    # now we reset the potentials
    CONFIG['POTCAR_FUNCTIONAL'] = 'PBE'
    CONFIG['POTCAR'].update(POTCAR_UPDATES)

    def __init__(self, structure, **kwargs):
        self.kwargs = copy.deepcopy(kwargs)
        uis = copy.deepcopy(self.kwargs.get('user_incar_settings', {}))
        new_config = copy.deepcopy(ForceConstantsSet.CONFIG)
        if 'ISPIN' not in uis:
            if magnetic_check(structure):
                uis.update({'ISPIN': 2})
            else:
                uis.update({'ISPIN': 1})

        if 'magmom' in uis:
            if 'MAGMOM' in new_config['INCAR']:
                new_config['INCAR'].pop('MAGMOM')
        elif uis['ISPIN']==1:
            if 'MAGMON' in uis.keys(): uis.pop['MAGMOM']
            if 'MAGMON' in new_config['INCAR']: new_config['INCAR'].pop['MAGMOM']

        new_config['INCAR'].update(uis)
        from pymatgen.io.vasp.inputs import Kpoints
        #kpoints = Kpoints.automatic_gamma_density(structure, 4000)
        kpoints = Kpoints(kpts=[[3,3,3],])
        new_config['KPOINTS'] = kpoints
        pot = self.kwargs.get('user_potcar_functional', None)
        if pot:
            new_config['POTCAR_FUNCTIONAL'] = pot
        super(ForceConstantsSet, self).__init__(
            structure, new_config, sort_structure=False, **self.kwargs)


class StaticSet(DictSet):
    """Set tuned for metal relaxations (correct smearing).
    Add `isif` parameter to the set to easily allow for overriding ISIF setting.
    Kpoints have a 6000 reciprocal density default.
    """
    CONFIG = _load_yaml_config("MPRelaxSet")
    CONFIG['KPOINTS'].update({
        'grid_density': 8000,
    })
    CONFIG['KPOINTS'].pop('reciprocal_density')  # to be explicit
    CONFIG['INCAR'].update({
        'EDIFF_PER_ATOM': 1e-6,
        'ENCUT': 520,  # MP compatibility
        'ISMEAR': -5,
        "NSW": 0,
        "IBRION": -1,
        'LREAL': False,
        'ALGO': 'NORMAL',
        # other settings from MPStaticSet
        "LAECHG": True,
        "LCHARG": True,
        "LWAVE": False,
        "LORBIT": 11,
        "LVHAR": True,
        "ICHARG": 2,
        "NEDOS": 5001,
    })
    # now we reset the potentials
    CONFIG['POTCAR_FUNCTIONAL'] = 'PBE'
    CONFIG['POTCAR'].update(POTCAR_UPDATES)

    def __init__(self, structure, isif=2, **kwargs):
        # pop the old kwargs, backwards compatibility from the complex StaticSet
        self.isif = isif
        
        old_kwargs = ['prev_incar', 'prev_kpoints', 'grid_density', 'lepsilon', 'lcalcpol']
        for k in old_kwargs:
            try:
                kwargs.pop(k)
            except KeyError:
                pass
        self.kwargs = copy.deepcopy(kwargs)
        uis = copy.deepcopy(self.kwargs.get('user_incar_settings', {}))
        uis['ISIF'] = isif
        new_config = copy.deepcopy(StaticSet.CONFIG)

        if 'ISPIN' not in uis:
            if magnetic_check(structure):
                uis.update({'ISPIN': 2})
            else:
                uis.update({'ISPIN': 1})

        if 'magmom' in uis:
            if 'MAGMOM' in new_config['INCAR']:
                new_config['INCAR'].pop('MAGMOM')
        elif uis['ISPIN']==1:
            if 'MAGMON' in uis.keys(): uis.pop['MAGMOM']
            if 'MAGMON' in new_config['INCAR']: new_config['INCAR'].pop['MAGMOM']

        new_config['INCAR'].update(uis)
        pot = self.kwargs.get('user_potcar_functional', None)
        if pot:
            new_config['POTCAR_FUNCTIONAL'] = pot
        super(StaticSet, self).__init__(structure, new_config, sort_structure=False, **self.kwargs)


class ATATIDSet():
    """Set tuned for Inflection Detection runs using ATAT with correct smearing for metals.
    Kpoints have a 8000 reciprocal density default.

    Overrides write_input to write the INCAR, KPPRA, USEPOT and DOSTATIC to the vasp.wrap instead.
    """

    def __init__(self, structure, grid_density=8000):
        self.structure = structure
        self.grid_density = grid_density


    def write_input(self, output_dir):
        """Write vasp.wrap and generate the other commands with robustrelax_vasp -mk"""
        # TODO: handle magmoms
        EDIFF_PER_ATOM = 1e-6
        EDIFF = len(self.structure)*EDIFF_PER_ATOM
        vasp_wrap = """[INCAR]
        EDIFF = {0}
        PREC = Accurate
        ALGO = Fast
        ENCUT = 520
        ISMEAR = 1
        SIGMA = 0.2
        IBRION = -1
        NSW = 1
        ISPIN = 2
        NELMIN = 4
        ISIF = 2
        LREAL = False
        ISYM = 0
        ICHARG = 1
        ISTART = 2
        USEPOT = PAWPBE
        KPPRA = {1}
        """.format(EDIFF, self.grid_density)
        with open(os.path.join(output_dir, 'vaspid.wrap'), 'w') as fp:
            fp.write(vasp_wrap)

class ForcesSet(DictSet):
    """Set tuned for generic force calculations (Gaussian smearing).
    """
    CONFIG = _load_yaml_config("MPRelaxSet")
    CONFIG['KPOINTS'].update({
        'grid_density': 8000,
    })
    CONFIG['KPOINTS'].pop('reciprocal_density')  # to be explicit
    CONFIG['INCAR'].pop('ENCUT',None)
    CONFIG['INCAR'].update({
        'EDIFF_PER_ATOM': 1e-8,
        'ISMEAR': 0,
        'SIGMA': 0.05,
        "NSW": 0,
        "IBRION": -1,
        'LREAL': False,
        'ALGO': 'NORMAL',
        # other settings from MPStaticSet
        "LCHARG": False,
        "LORBIT": 11,
        "LVHAR": True,
        "LWAVE": False,
        "ICHARG": 2,
        "NEDOS": 5001,
    })
    # now we reset the potentials
    CONFIG['POTCAR_FUNCTIONAL'] = 'PBE'
    CONFIG['POTCAR'].update(POTCAR_UPDATES)

    def __init__(self, structure, **kwargs):
        # pop the old kwargs, backwards compatibility from the complex StaticSet
        old_kwargs = ['prev_incar', 'prev_kpoints', 'grid_density', 'lepsilon', 'lcalcpol']
        for k in old_kwargs:
            try:
                kwargs.pop(k)
            except KeyError:
                pass
        self.kwargs = copy.deepcopy(kwargs)
        super(ForcesSet, self).__init__(structure, ForcesSet.CONFIG, sort_structure=False, **self.kwargs)


class BornChargeSet(DictSet):
    """Set tuned for metal relaxations (correct smearing).
    Add `isif` parameter to the set to easily allow for overriding ISIF setting.
    Kpoints have a 6000 reciprocal density default.
    Provide by Yi Wang, Aug. 24, 2020
    """
    CONFIG = _load_yaml_config("MPRelaxSet")
    CONFIG['KPOINTS'].update({
        'grid_density': 8000,
    })
    CONFIG['KPOINTS'].pop('reciprocal_density')  # to be explicit

    CONFIG['INCAR'].update({
        "LEPSILON": True,  #For Born Effective Charge
        "NCORE": 1,
        #"NPAR": 1,
        #"LCALCEPS": True,  #For Born Effective Charge
        "LRPA": False,
        'EDIFF_PER_ATOM': 1e-6,
        'ENCUT': 520,  # MP compatibility
        'ISMEAR': 0,
        "NSW": 0,
        "IBRION": -1,
        'LREAL': False,
        'ALGO': 'NORMAL',
        # other settings from MPStaticSet
        "LCHARG": False,
        "LWAVE": False,
        "ICHARG": 2,
    })

    CONFIG['POTCAR_FUNCTIONAL'] = 'PBE'
    CONFIG['POTCAR'].update(POTCAR_UPDATES)

    def __init__(self, structure, isif=2, **kwargs):
        # pop the old kwargs, backwards compatibility from the complex StaticSet
        self.isif = isif

        uis = copy.deepcopy(kwargs.get('user_incar_settings', {}))
        new_config = copy.deepcopy(BornChargeSet.CONFIG)
        """
        old_kwargs = ['prev_incar', 'prev_kpoints', 'grid_density', 'lepsilon', 'lcalcpol', \
            'user_potcar_functional', 'user_incar_settings']
        """
        old_kwargs = ['prev_incar', 'prev_kpoints', 'grid_density', 'lepsilon', 'lcalcpol']
        for k in old_kwargs:
            try:
                kwargs.pop(k)
            except KeyError:
                pass
        self.kwargs = copy.deepcopy(kwargs)

        if 'ISPIN' not in uis:
            if magnetic_check(structure):
                uis.update({'ISPIN': 2})
            else:
                uis.update({'ISPIN': 1})

        if 'magmom' in uis:
            if 'MAGMOM' in new_config['INCAR']:
                new_config['INCAR'].pop('MAGMOM')
        elif uis['ISPIN']==1:
            if 'MAGMON' in uis.keys(): uis.pop['MAGMOM']
            if 'MAGMON' in new_config['INCAR']: new_config['INCAR'].pop['MAGMOM']

        for key in uis.keys():
            if key not in new_config['INCAR']:
                if key in {'NELM', 'EDIFF', 'NEDOS', 'KPOINT_BSE'} : continue
                new_config['INCAR'][key] = uis[key]
            elif key == 'ISPIN':
                new_config['INCAR'][key] = uis[key]
            elif key == 'ISMEAR':
                new_config['INCAR'][key] = uis[key]
            elif key == 'SIGMA':
                new_config['INCAR'][key] = uis[key]

        if 'SIGMA' in new_config['INCAR'] and 'ISMEAR' in new_config['INCAR'] :
            if new_config['INCAR']['ISMEAR'] == -5:
                new_config['INCAR'].pop('SIGMA')
        pot = self.kwargs.get('user_potcar_functional', None)
        if pot:
            new_config['POTCAR_FUNCTIONAL'] = pot
        super(BornChargeSet, self).__init__(structure, new_config, sort_structure=False, **self.kwargs)


class ElasticSet(DictSet):
    """Set tuned for metal relaxations (correct smearing).
    Add `isif` parameter to the set to easily allow for overriding ISIF setting.
    Kpoints have a 6000 reciprocal density default.
    """
    CONFIG = _load_yaml_config("MPRelaxSet")
    #    'EDIFF_PER_ATOM': 1e-6,
    CONFIG['INCAR'] = {
        'EDIFF': 1e-6,
        'ENCUT': 520,  # MP compatibility
        'ISMEAR': -5,
        "IBRION": 2,
        'LREAL': False,
        'ALGO': 'NORMAL',
        # other settings from MPStaticSet
        "LAECHG": True,
        "LCHARG": True,
        "LWAVE": False,
        #"LORBIT": 11,
        "LVHAR": True,
        "ICHARG": 0,
        "NSW": 99,
        #"MAGMOM": CONFIG['INCAR']['MAGMOM'],
        "ISPIN": 2,
        "ISIF": 2,
        "PREC": "High"
    }
    # now we reset the potentials
    CONFIG['POTCAR_FUNCTIONAL'] = 'PBE'
    CONFIG['POTCAR'].update(POTCAR_UPDATES)

    def __init__(self, structure, **kwargs):
        # pop the old kwargs, backwards compatibility from the complex StaticSet
        
        uis = copy.deepcopy(kwargs.get('user_incar_settings', {}))
        new_config = copy.deepcopy(ElasticSet.CONFIG)
        """
        old_kwargs = ['prev_incar', 'prev_kpoints', 'grid_density', 'lepsilon', 'lcalcpol', \
            'user_potcar_functional', 'user_incar_settings']
        """
        old_kwargs = ['prev_incar', 'prev_kpoints', 'grid_density', 'lepsilon', 'lcalcpol']
        for k in old_kwargs:
            try:
                kwargs.pop(k)
            except KeyError:
                pass
        self.kwargs = copy.deepcopy(kwargs)

        if 'ISPIN' not in uis:
            if magnetic_check(structure):
                uis.update({'ISPIN': 2})
            else:
                uis.update({'ISPIN': 1})

        if 'magmom' in uis:
            if 'MAGMOM' in new_config['INCAR']:
                new_config['INCAR'].pop('MAGMOM')
        elif uis['ISPIN']==1:
            if 'MAGMON' in uis.keys(): uis.pop['MAGMOM']
            if 'MAGMON' in new_config['INCAR']: new_config['INCAR'].pop['MAGMOM']

        for key in uis.keys():
            if key not in new_config['INCAR']:
                if key in {'NELM', 'EDIFF', 'NEDOS', 'KPOINT_BSE'} : continue
                new_config['INCAR'][key] = uis[key]
            elif key == 'ISPIN':
                new_config['INCAR'][key] = uis[key]
            elif key == 'ISMEAR':
                new_config['INCAR'][key] = uis[key]
            elif key == 'SIGMA':
                new_config['INCAR'][key] = uis[key]

        if 'SIGMA' in new_config['INCAR'] and 'ISMEAR' in new_config['INCAR'] :
            if new_config['INCAR']['ISMEAR'] == -5:
                new_config['INCAR'].pop('SIGMA')

        from pymatgen.io.vasp.inputs import Kpoints
        if metal_check(structure):
            grid_density = 15625
            #new_config['INCAR']['ISMEAR'] = 1
            #new_config['INCAR']['SIGMA'] = 0.2
        else:
            grid_density = 8000
        kpoints = Kpoints.automatic_gamma_density(structure, grid_density)
        new_config['KPOINTS'] = kpoints
        pot = self.kwargs.get('user_potcar_functional', None)
        if pot:
            new_config['POTCAR_FUNCTIONAL'] = pot
        super(ElasticSet, self).__init__(structure, new_config, sort_structure=False, **self.kwargs)
        self.config = new_config


