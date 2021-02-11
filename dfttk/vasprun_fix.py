from pymatgen.io.vasp.outputs import Vasprun
import xml.etree.ElementTree as ET

def float_string(s, verbose=True):
    """Try parsing a string as float, fixing errors if the string is #.######-### (e.g. not #.######E-###, missing E). Return the correct string."""
    try:
        float(s)
    except ValueError:
        x = s.split('-')
        new_s = ''.join(x[:-1]) + 'E-' + x[-1]
        if verbose:
            print('Fixing vasprun: {} to {}'.format(s, new_s))
        s = new_s

    return s

def fix_vasprun(fn, outfn=None, validate=True):
    """
    Fix vasprun files with bad float values in the dynamical matrix eigenvectors and rewrite out to the same file.

    Parameters
    ----------
    fn : string
        Input filename of an XML document.
    outfn : string, optional
        Output filename of the XML document. If None (default), will choose the input filename.
    validate : bool
        If True, will validate the outfn

    Returns
    -------

    """
    # Set them to 1e-300, so approx 0, but not 0 (possible numerical stability issue).
    # Not sure if this is exactly right, but it's better than failing.
    tree = ET.parse(fn)
    root = tree.getroot()
    last_calculation = root.findall("calculation")[-1] # contains the dynmat
    dynmat = last_calculation.find("dynmat")
    eigenvectors = dynmat.findall("varray")[-1]
    for v in eigenvectors:
        temp_text = v.text
        v.text = ' '.join([float_string(s) for s in temp_text.split()])
    if outfn is None:
        outfn = fn
    tree.write(outfn)
    if validate:
        v = Vasprun(outfn)

