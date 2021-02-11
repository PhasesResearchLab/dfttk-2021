#!python

import pymatgen
import os
import json
from monty.serialization import loadfn

def build_str_from_prototype(prototype, prototype_type="s"):
    """
    Build the structure from prototype.
    The prototype database file(aflow_prototypes.json) of pymatgen is used
        (The original is from AFLOW[http://www.aflowlib.org/CrystalDatabase/])
    Note: If this function is used, please cite:
        Mehl, M. J., Hicks, D., Toher, C., Levy, O., Hanson, R. M., Hart, G., & Curtarolo, S. (2017).
        The AFLOW library of crystallographic prototypes: part 1.
        Computational Materials Science, 136, S1-S828.
        http://doi.org/10.1016/j.commatsci.2017.01.017

    Parameters
    ----------
        prototype: str
            The prototype
        prototype_type: str
            The type of the prototype, only the first letter is used
                "p" for "pearson"
                "a" for "aflow"
                "s" for "strukturbericht"
                "m" for "mineral"
    Returns
    -------
    """
    prototype_type_dict = {"p": "pearson", "a": "aflow", "s": "strukturbericht", "m": "mineral"}
    prototype_type = prototype_type_dict[prototype_type[0].lower()]
    #print(prototype_type)
    AFLOW_PROTOTYPE_LIBRARY = loadfn(os.path.join(os.path.dirname(pymatgen.__file__),
                                              "analysis/aflow_prototypes.json"))
    struct = []
    for proto_i in AFLOW_PROTOTYPE_LIBRARY:
        if prototype == proto_i['tags'][prototype_type]:
            #print(proto_i['tags']["mineral"])
            struct.append(proto_i['snl'].structure)
            #proto_i['snl'].structure.to(fmt="POSCAR", filename="POSCAR")
    return struct