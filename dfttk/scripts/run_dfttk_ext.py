# -*- coding: utf-8 -*-
# The template for batch run of DFTTK
import argparse
from pymatgen.ext.matproj import MPRester, Structure
from pymatgen.io.vasp.inputs import Potcar
#from dfttk.wflows import get_wf_gibbs, get_wf_EV_bjb, get_wf_gibbs_robust
from dfttk.wflows import get_wf_gibbs
from dfttk.utils import recursive_glob
from dfttk.structure_builders.parse_anrl_prototype import multi_replace
from fireworks.fw_config import config_to_dict
from monty.serialization import loadfn, dumpfn
from atomate.vasp.database import VaspCalcDb
import dfttk.pythelec as pythelec
from dfttk.pythelec import thelecMDB
from dfttk.pythfind import thfindMDB
from dfttk.pyEVfind import EVfindMDB
import warnings
import copy
import os
import sys
import shutil
from datetime import datetime
from dfttk.analysis.ywplot import myjsonout, thermoplot
from dfttk.analysis.ywutils import get_melting_temperature, reduced_formula, get_expt

no_MongoDB = False

def ext_thelec(args, plotfiles=None, vasp_db=None):
    global no_MongoDB
    print ("Postprocess for thermodynamic properties, Seebeck, Lorenz number etc. Yi Wang\n")
    """
    Postprocess for thermodynamic properties, Seebeck, Lorenz number etc

    Parameters
        STR_FOLDER = args.STRUCTURE_FOLDER
            folder/file containing structures
        MATCH_PATTERN = args.MATCH_PATTERN
            Match patterns for structure file, e.g. *POSCAR
        RECURSIVE = args.RECURSIVE
            recursive or not
        WORKFLOW = args.WORKFLOW
            workflow, current only get_wf_gibbs
        LAUNCH = args.LAUNCH
            Launch to lpad or not
        MAX_JOB = args.MAX_JOB
            Max job to submit
        SETTINGS = args.SETTINGS
            Settings file
        WRITE_OUT_WF = args.WRITE_OUT_WF
            Write out wf file or not
    """
    t0 = args.t0
    t1 = args.t1
    td = args.td
    xdn = args.xdn
    xup = args.xup
    ndosmx = args.ndosmx
    natom = args.natom
    gaussian = args.gaussian
    dope = args.dope
    doscar = args.doscar
    outf = args.outf
    qhamode = args.qhamode
    eqmode = args.eqmode
    elmode = args.elmode
    metatag = args.metatag
    everyT = args.everyT
    noel = args.noel
    smooth = args.smooth
    expt = args.expt
    xlim = args.xlim
    if abs(dope)<5.e-9:
        ndosmx = max(100001, int(ndosmx))
        gaussian = max(10000., float(gaussian))

    formula = None
    if args.local != "":
        print("\nRun using local data\n")

        readme = {}
        record_cmd(readme)
        proc = thelecMDB(t0, t1, td, xdn, xup, dope, ndosmx, gaussian, natom, outf,
            noel=noel, qhamode=qhamode, eqmode=eqmode, elmode=elmode,
            smooth=smooth, debug=args.debug, phasename=args.local,
            pyphon=True, renew=args.renew, fitF=args.fitF, args=args)
        volumes, energies, thermofile, comments = proc.run_console()

        if comments!=None: readme.update(comments)
        else: return
        if "ERROR" in readme.keys():
            record_cmd_print(thermofile, readme, dir=args.phasename)
            return

        print("\nFull thermodynamic properties have outputed into:", thermofile)
        #print(args.plot, "eeeeeeeee", volumes, energies, thermofile, comments)
        if args.plot==None: print("\nSupply '-plot phasename' for plot\n")
        else:
            from dfttk.analysis.ywplot import plotAPI
            #print("xxxxxxx",proc.get_formula())
            if plotAPI(readme, thermofile, volumes, energies, expt=expt, xlim=xlim, _fitCp=args.SGTEfitCp,
                formula = proc.get_formula(), debug=args.debug,
                plotlabel=args.plot, local=args.local):
                #print ("xxxxxxx",proc.get_formula())
                vtof = proc.get_free_energy_for_plot(readme)
                if vtof is not None:
                    plotAPI(readme, thermofile, volumes, energies, expt=expt, xlim=xlim, _fitCp=args.SGTEfitCp,
                    formula = proc.get_formula(), vtof=vtof, plotlabel=args.plot)
            """
            """
        #print("xxxxxxxxxxxx",readme)
        record_cmd_print(thermofile, readme)
    elif no_MongoDB:
            print("\n*********WARNING: CANNOT get MongoDB service, so I will proceed using local data")
            print("*********WARNING: CANNOT get MongoDB service, so I will proceed using local data")
            print("*********WARNING: CANNOT get MongoDB service, so I will proceed using local data\n")
    elif not args.plotonly:
        #if True:
        try:
            if vasp_db is None:
                db_file = loadfn(config_to_dict()["FWORKER_LOC"])["env"]["db_file"]
                vasp_db = VaspCalcDb.from_db_file(db_file, admin=False)
            static_calculations = vasp_db.collection.\
                find({'$and':[ {'metadata.tag': metatag}, {'adopted': True} ]})
            structure = Structure.from_dict(static_calculations[0]['output']['structure'])
            formula = reduced_formula(structure.composition.alphabetical_formula)
        except:
            no_MongoDB = True
            print("\n*********WARNING: CANNOT get MongoDB service, so I will proceed using local data")
            print("*********WARNING: CANNOT get MongoDB service, so I will proceed using local data")
            print("*********WARNING: CANNOT get MongoDB service, so I will proceed using local data\n")
        """
        """

    #call API
    if args.plotonly and plotfiles!=None:
        metatag, thermofile, volumes, energies, dir, formula = plotfiles
        sys.stdout.write('Processing {}, dir: {}, formula: {}\n'.format(metatag, dir, formula))

        #print(thermofile, volumes, energies, formula)
        #print(thermofile, dir, formula)
        readme={}
        from dfttk.analysis.ywplot import plotAPI
        plotAPI(readme, thermofile, None, energies, expt=expt, xlim=xlim, _fitCp=args.SGTEfitCp,
            formula = formula, vtof=None, plotlabel=args.plot)
    elif vasp_db==None and plotfiles!=None:
        metatag, thermofile, volumes, energies, dir, formula = plotfiles
        sys.stdout.write('Processing {}, dir: {}, formula: {}\n'.format(metatag, dir, formula))
        #print("xxxxxxxxxx", plotfiles)
        if expt!=None:
            _t1 = get_melting_temperature(expt, formula)
            if _t1!=None: t1 = _t1

        readme = {}
        record_cmd(readme)
        proc = thelecMDB(t0, t1, td, xdn, xup, dope, ndosmx, gaussian, natom, outf, vasp_db=vasp_db,
            noel=noel, metatag=metatag, qhamode=qhamode, eqmode=eqmode, elmode=elmode, everyT=everyT,
            smooth=smooth, debug=args.debug,
            phasename=dir, pyphon=args.pyphon, renew=args.renew, fitF=args.fitF, args=args)
        volumes, energies, thermofile, comments = proc.run_console()
        #print ("xxxxxxx", comments)
        if comments!=None: readme.update(comments)
        else: return
        if "ERROR" in readme.keys():
            #record_cmd_print(thermofile, readme, dir=args.phasename)
            return

        print("\nFull thermodynamic properties have outputed into:", thermofile)
        #print(args.plot, "eeeeeeeee", volumes, energies, thermofile, comments)
        if args.plot==None: print("\nSupply '-plot phasename' for plot\n")
        else:
            from dfttk.analysis.ywplot import plotAPI
            if plotAPI(readme, thermofile, volumes, energies, expt=expt, xlim=xlim, _fitCp=args.SGTEfitCp,
                formula = proc.get_formula(), debug=args.debug,
                plotlabel=args.plot):
                vtof = proc.get_free_energy_for_plot(readme)
                if vtof is not None:
                    plotAPI(readme, thermofile, volumes, energies, expt=expt, xlim=xlim, _fitCp=args.SGTEfitCp,
                    formula = proc.get_formula(), vtof=vtof, plotlabel=args.plot)
            """
            """
        record_cmd_print(thermofile, readme)

    elif metatag != None:
        if expt!=None:
            _t1 = get_melting_temperature(expt, formula)
            if _t1!=None: t1 = _t1
        readme = {}
        record_cmd(readme)
        proc = thelecMDB(t0, t1, td, xdn, xup, dope, ndosmx, gaussian, natom, outf, vasp_db=vasp_db,
            noel=noel, metatag=metatag, qhamode=qhamode, eqmode=eqmode, elmode=elmode, everyT=everyT,
            smooth=smooth, debug=args.debug,
            phasename=args.phasename, pyphon=args.pyphon, renew=args.renew, fitF=args.fitF, args=args)
        volumes, energies, thermofile, comments = proc.run_console()

        if comments!=None: readme.update(comments)
        else: return
        if "ERROR" in readme.keys():
            record_cmd_print(thermofile, readme, dir=args.phasename)
            return

        print("\nFull thermodynamic properties have outputed into:", thermofile)
        #print(args.plot, "eeeeeeeee", volumes, energies, thermofile, comments)
        if args.plot==None: print("\nSupply '-plot phasename' for plot\n")
        else:
            from dfttk.analysis.ywplot import plotAPI
            if plotAPI(readme, thermofile, volumes, energies, expt=expt, xlim=xlim, _fitCp=args.SGTEfitCp,
                formula = proc.get_formula(), debug=args.debug,
                plotlabel=args.plot):
                vtof = proc.get_free_energy_for_plot(readme)
                if vtof is not None:
                    plotAPI(readme, thermofile, volumes, energies, expt=expt, xlim=xlim, _fitCp=args.SGTEfitCp,
                    formula = proc.get_formula(), vtof=vtof, plotlabel=args.plot)
            """
            """
        record_cmd_print(thermofile, readme)
    elif args.vdos is not None:
        readme = {}
        record_cmd(readme)
        proc = thelecMDB(t0, t1, td, xdn, xup, dope, ndosmx, gaussian, natom, outf, renew=True, args=args)
        thermofile, comments, natoms = proc.run_single()
        if thermofile is None: return
        readme.update(comments)
        #record_cmd(thermofile, readme)

        print("\nFull thermodynamic properties have outputed into:", thermofile)
        if args.plot!=None:
            from dfttk.analysis.ywplot import plotAPI
            if plotAPI(readme, thermofile, None, None, expt=expt, xlim=xlim, _fitCp=args.SGTEfitCp,
                formula = proc.get_formula(), debug=args.debug,
                poscar=args.poscar,vdos=args.vdos, doscar=args.doscar, natoms=natoms, plotlabel=args.plot):
                record_cmd_print(thermofile, readme)
    elif args.local == "":
        pythelec.thelecAPI(t0, t1, td, xdn, xup, dope, ndosmx, gaussian, natom, outf, doscar)


def record_cmd(readme):
    cmdline = copy.deepcopy(sys.argv)
    cmdline[0] = cmdline[0].split('/')[-1]
    readme['command']='{}'.format(' '.join(cmdline))
    readme['Postprocess time']='{}'.format(datetime.now())
    #readme['start at']='{}'.format(datetime.now())
    #fp.write('#These results are produced by the following command line on {}\n'.format(datetime.now()))
    #fp.write('{}\n'.format(' '.join(cmdline)))


def record_cmd_print(fdir, readme, dir=None):
    dir = fdir
    #readme['finished at'] = '{}'.format(datetime.now())
    if not os.path.isdir(dir):
        dir = '/'.join(fdir.split('/')[0:-1])
        if dir == "": dir = "./"
        with open (dir+"/readme", "w") as fp:
            myjsonout(readme, fp, indent="", comma="")

        if "ERROR" in readme.keys():
            error ="**********FATAL ERROR encountered, you may check readme and E-V plot in the folder "+dir+"/figures"
            volumes = readme['E-V']['volumes']
            energies = readme['E-V']['energies']
            folder = os.path.join(dir,'figures')
            if not os.path.exists(folder): os.mkdir(folder)
            thermoplot(folder,"0 K total energies (eV/atom)",volumes, energies, plottitle=dir, lp=True)
            with open (os.path.join(dir,"ERROR"), "w") as fp:
                fp.write('{}\n'.format(readme['ERROR']))
                if dir!=None:
                    fp.write('{}\n'.format(error))
                    print ("\n", error, "\n")
        else:
            if os.path.exists(os.path.join(dir,"ERROR")): os.remove(os.path.join(dir,"ERROR"))


        with open ("runs.log", "a") as fp:
            try:
                #fp.write('phonon quality={}, LTC zigzag={}'.format(readme['phonon quality'], readme['LTC quality']))
                fp.write('phonon quality={}'.format(readme['phonon quality']))
                if os.path.exists(dir+'/fitF'): fp.write(', for fitF is  on: {}\n'.format(dir))
                else: fp.write(', for fitF is off: {}\n'.format(dir))
            except:
                fp.write('FATAL ERROR in {}\n'.format(dir))
                pass

def shared_aguments(pthelec):
    pthelec.add_argument("-py", "--pyphon", dest="pyphon", action='store_true', default=False,
                      help="use Yphon to recalculate vibrational properties. \n"
                           "Default: False")
    pthelec.add_argument("-T0", "-t0", dest="t0", nargs="?", type=float, default=0.0,
                      help="Low temperature limit. \n"
                           "Default: 0")
    pthelec.add_argument("-dgx", "--debye_gruneisen_x", dest="debye_gruneisen_x", nargs="?", type=float, default=2./3,
                      help="Value of x for the Debye gruneisen model. \n"
                           "Default: 2/3")
    pthelec.add_argument("-db_repair", "--db_repair", dest="db_repair", action='store_true', default=False,
                      help="repair database for prebiously FIZZLED calculation. \n"
                           "Default: False")
    pthelec.add_argument("-db_renew", "--db_renew", dest="db_renew", action='store_true', default=False,
                      help="renew the database. \n"
                           "Default: False")
    pthelec.add_argument("-T1", "-t1", dest="t1", nargs="?", type=float, default=4000,
                      help="High temperature limit. \n"
                           "Default: 4000")
    pthelec.add_argument("-dT", "-td", dest="td", nargs="?", type=float, default=10,
                      help="Temperature increment. \n"
                           "Default: 10")
    pthelec.add_argument("-xdn", "--xdn", dest="xdn", nargs="?", type=float, default=-100,
                      help="Low band energy limit. \n"
                           "Default: -100 (eV)")
    pthelec.add_argument("-xup", "--xup", dest="xup", nargs="?", type=float, default=100,
                      help="High band energy limit. \n"
                           "Default: 100")
    pthelec.add_argument("-dope", "--dope", dest="dope", nargs="?", type=float, default=0,
                      help="dope level (electrons). \n"
                           "Default: 0, change it to -1e-5 can speed up the calculation while \n"
                           "it may destropy the numerical stability")
    pthelec.add_argument("-ne", "--ndosmx", dest="ndosmx", nargs="?", type=int, default=10001,
                      help="new eDOS mesh. It recommend increase it 100001 if numberical instability seen. \n"
                           "Default: 10001 if dope >5.e-9 or 100001")
    pthelec.add_argument("-gauss", "--gauss", dest="gaussian", nargs="?", type=float, default=1000.,
                      help="densing factor for eDOS mesh near the Fermi energy. \n"
                           "Default: 1000")
    pthelec.add_argument("-natom", "--natom", dest="natom", nargs="?", type=int, default=1,
                      help="number of atoms in the DOSCAR. \n"
                           "Default: 1")
    pthelec.add_argument("-l", "--local", dest="local", nargs="?", type=str, default="",
                      help="path holding all input for Yphon. \n"
                           "Default: blank")
    pthelec.add_argument("-nT", "--nT", dest="nT", nargs="?", type=int, default=257,
                      help="number of temperatures, used together with -td -50. \n"
                           "Default: 257")
    pthelec.add_argument("-e", "--everyT", dest="everyT", nargs="?", type=int, default=1,
                      help="number of temperature points skipped from QHA analysis from the qha/qha_phonon collection. \n"
                           "Default: 1")
    pthelec.add_argument("-o", "-outf", dest="outf", nargs="?", type=str, default="fvib_ele",
                      help="output filename for calculated thermoelectric properties. \n"
                           "Default: fvib_ele")
    pthelec.add_argument("-noel", "-noel", dest="noel", action='store_true', default=False,
                      help="do not consider the thermal electron contribution. \n"
                           "Default: False")
    pthelec.add_argument("-tag", "--metatag", dest="metatag", nargs="?", type=str, default=None,
                      help="metatag: MongoDB metadata tag field. \n"
                           "Default: None")
    pthelec.add_argument("-qhamode", "-qhamode", dest="qhamode", nargs="?", type=str, default=None,
                      help="quasiharmonic mode: debye, phonon, or yphon. \n"
                           "Default: None")
    pthelec.add_argument("-pn", "-phasename", dest="phasename", nargs="?", type=str, default=None,
                      help="assigan phase name. \n"
                           "Default: None")
    pthelec.add_argument("-jp", "-jobpath", dest="jobpath", nargs="?", type=str, default=None,
                      help="For debug/development purpoase. Parent path where jobs were submittedi to check settings. \n"
                           "Default: None")
    pthelec.add_argument("-kpm", "--k_ph_mode", dest="k_ph_mode", nargs="?", type=int, default=0,
                      help="Mode to calculate thermal conductivity using slack model.\n"
                           "    0: Original model, using the Debye T calc. at lowest T; \n"
                           "    1: using thermodynamic gruneisen gamma and T-dependent debye T;  \n"
                           "Default: 0")
    pthelec.add_argument("-eq", "--eqmode", dest="eqmode", nargs="?", type=int, default=4,
                      help="Mode to calculate equilibrium volume and LTC.\n"
                           "    0: Symmetrical Central differential if the data is excellent; \n"
                           "    4: 4-parameter BM fitting if the data is faitly good;  \n"
                           "    5: 5-parameter BM fitting if the data is very good.  \n"
                           "Default: 4")
    pthelec.add_argument("-el", "--elmode", dest="elmode", nargs="?", type=int, default=0,
                      help="Mode to interpolate thermal electronic contribution:"
                           "                       0: interp1d;  \n"
                           "                       1: UnivariateSpline.  \n"
                           "Default: 0")
    pthelec.add_argument("-s", "-smooth", dest="smooth", action='store_true', default=False,
                      help="smooth the LTC. \n"
                           "Default: False")
    pthelec.add_argument("-plot", "-plot", dest="plot", nargs="?", type=str, default=None,
                      help="plot the figures and mark the theoretial line with the given label. \n"
                           "Default: None")
    pthelec.add_argument("-renew", "-renew", dest="renew", action='store_true', default=False,
                      help="renew/plot the figure. Otherwise, calculation will be skipped if the file 'fvib_ele' is seen.\n"
                           "Default: False")
    pthelec.add_argument("-refresh", "-refresh", dest="refresh", action='store_true', default=False,
                      help="recalculate the phonon dos.\n"
                           "Default: False")
    pthelec.add_argument("-fitCp", "--SGTEfitCp", dest="SGTEfitCp", action='store_true', default=False,
                      help="report SGTE fitting through the order of Cp, S, and H. \n"
                           "Default: False, report SGTE fitting using Gibbs energy.")
    pthelec.add_argument("-fitF", "-fitF", dest="fitF", action='store_true', default=False,
                      help="USE with CARE! apply for the case of poor data quality. Enforce linear\n"
                           "fitting to the vibrational and electronic free energy together enforce 4-parameter\n"
                           "Birch–Murnaghan fitting to the 0 K static energies. A file named 'fitF' must be\n"
                           "set in that phase folder using 'touch fitF' to make the option take effect.\n"
                           "Note: The regualar run will skip the phase folder if a file name 'fitF' is seen there.\n"
                           "Default: False")
    pthelec.add_argument("-po", "--plotonly", dest="plotonly", action='store_true', default=False,
                      help="No calculations, only plot the calulcated results if any. \n"
                           "Default: False")
    pthelec.add_argument("-g", "--debug", dest="debug", action='store_true', default=False,
                      help="turn on debug mode by reducing the mesh. \n"
                           "Default: False")
    pthelec.add_argument("-expt", "-expt", dest="expt", nargs="?", type=str, default=None,
                      help="file path (json format, list of dictionary) for experimental thermodynamic properties to \n"
                           "be compared with. Default: None")
    pthelec.add_argument("-xlim", "-xlim", dest="xlim", nargs="?", type=float, default=None,
                      help="Up temperature limit for plot. \n"
                           "Default: None")
    pthelec.add_argument("-dos", "--doscar", dest="doscar", nargs="?", type=str, default=None,
                      help="file path to DOSCAR file. Run thelec in single volume shot only. \n"
                           "Default: None")
    pthelec.add_argument("-pos", "--poscar", dest="poscar", nargs="?", type=str, default=None,
                      help="file path to POSCAR file. Run thelec in single volume shot only. \n"
                           "Default: None")
    pthelec.add_argument("-vdos", "--vdos", dest="vdos", nargs="?", type=str, default=None,
                      help="file path to phonon DOS file produced by Yphon. Run thelec in single volume shot only. \n"
                           "Default: None")


def run_ext_thelec(subparsers):
    # begin process by Yi Wang, July 23, 2020
    #SUB-PROCESS: thelec
    pthelec = subparsers.add_parser("thelec", help="Postprocess DFTTK results after DFT job completed.")
    shared_aguments(pthelec)
    pthelec.set_defaults(func=ext_thelec)
    # end process by Yi Wang, July 23, 2020

    #further extension for finding phonon calculation
    run_ext_thfind(subparsers)
    run_ext_EVfind(subparsers)


def run_ext_thfind(subparsers):
    #SUB-PROCESS: thfind
    if subparsers is not None:
        pthfind = subparsers.add_parser("thfind", help="Check the dfttk DFT calculation results followed by calling the 'thelec' module to get thermodynamic properties when the option '-get' is given.")
    else:
        pthfind = argparse.ArgumentParser(description='get all default DFTTK args.')
    pthfind.add_argument("-w", "--within", dest="within", nargs="?", type=str, default=None,
                      help="find calculations within element list\n"
                           "Default: None")
    pthfind.add_argument("-all", "--containall", dest="containall", nargs="?", type=str, default=None,
                      help="find calculations must contain all elements in the list\n"
                           "Default: None")
    pthfind.add_argument("-xall", "--excludeall", dest="excludeall", nargs="?", type=str, default=None,
                      help="exclude calculations that conain all elements in the list\n"
                           "Default: None")
    pthfind.add_argument("-xany", "--excludeany", dest="excludeany", nargs="?", type=str, default=None,
                      help="exclude calculations that conain any elements in the list\n"
                           "Default: None")
    pthfind.add_argument("-any", "--containany", dest="containany", nargs="?", type=str, default=None,
                      help="find calculations contain any elements in the list\n"
                           "Default: None")
    pthfind.add_argument("-v", "--nV", dest="nV", nargs="?", type=int, default=6,
                      help="Return phonon calculations finished for number of volumes larger or equals to. \n"
                           "Default: 6")
    pthfind.add_argument("-ss", "--supercellsize", dest="supercellN", nargs="?", type=int, default=0,
                      help="only return phonon calculation with supercell size larger than. \n"
                           "Default: 0")
    pthfind.add_argument("-fg", "--findbandgap", dest="findbandgap", action='store_true', default=False,
                      help="report the entries with band gap. \n"
                           "Default: False")
    pthfind.add_argument("-get", "--get", dest="get", action='store_true', default=False,
                      help="call thelec module to get the thermodyamic data for all found entries. \n"
                           "Default: False")
    pthfind.add_argument("-db", "--db_file", dest="db_file", nargs="?", type=str, default=None,
                      help="alternative database other than default\n"
                           "Default: None")
    pthfind.add_argument("-check", "--check", dest="check", action='store_true', default=False,
                      help="check database. \n"
                           "Default: False")
    pthfind.add_argument("-remove", "--remove", dest="remove", action='store_true', default=False,
                      help="remove database document entries under given conditions (under development). \n"
                           "Default: False")
    shared_aguments(pthfind)

    if subparsers is None: return pthfind

    pthfind.set_defaults(func=ext_thfind)


def ext_thfind(args, vasp_db=None):
    """
    find the metadata tag that has finished.

    Parameters
        STR_FOLDER = args.STRUCTURE_FOLDER
            folder/file containing structures
        MATCH_PATTERN = args.MATCH_PATTERN
            Match patterns for structure file, e.g. *POSCAR
        RECURSIVE = args.RECURSIVE
            recursive or not
        WORKFLOW = args.WORKFLOW
            workflow, current only get_wf_gibbs
    """
    qha_renew = args.db_repair or args.db_renew
    if args.db_file is not None:
        vasp_db = VaspCalcDb.from_db_file(args.db_file, admin=qha_renew)
        db_file = args.db_file       
    elif vasp_db is None:
        db_file = loadfn(config_to_dict()["FWORKER_LOC"])["env"]["db_file"]
        vasp_db = VaspCalcDb.from_db_file(db_file, admin=qha_renew)
    proc=thfindMDB(args,vasp_db)
    tags = proc.run_console()


    if qha_renew:
        from dfttk.scripts.QHAAnalysis_renew import QHAAnalysis_renew
        for t in tags:
            if isinstance(t,dict):
                tag = t['tag']

                qha_phonon = list(vasp_db.db['qha_phonon'].find({'metadata.tag': tag})) 
                if len(qha_phonon) > 0 : 
                    vasp_db.db['qha_phonon'].remove({'metadata.tag': tag})
                    print('The data with metadata.tag={} in {} collection is removed'.format(tag, 'qha_phonon'))
                qha = list(vasp_db.db['qha'].find({'metadata.tag': tag})) 
                if len(qha) > 0 : 
                    vasp_db.db['qha'].remove({'metadata.tag': tag})
                    print('The data with metadata.tag={} in {} collection is removed'.format(tag, 'qha'))
                    
                try:
                    phonon_calculations = list(vasp_db.db['phonon'].find({'$and':[ {'metadata.tag': tag}, {'adopted': True} ]}))
                    T = phonon_calculations[0]['temperatures']
                    t_min = min(T)
                    t_max = max(T)
                    t_step = T[1]-T[0]
                except:
                    t_min = 5
                    t_max = 2000
                    t_step = 5
                print("Renew calculation of Debye-gruneisen model for metadata tag:", t)
  
                proc = QHAAnalysis_renew(phonon=True, t_min=t_min, t_max=t_max,
                t_step=t_step, db_file=db_file, test_failure=False, admin=True,
                metadata={'tag':tag}, bp2gru = args.debye_gruneisen_x, 
                tag=tag)
                proc.run_task()
                """
                try: 
                    proc.run_task()
                except:
                    continue
                """

    if args.get:
        with open("runs.log", "a") as fp:
            fp.write ('\nPostprocessing run at {}\n\n'.format(datetime.now()))
        for t in tags:
            if isinstance(t,dict):
                print("\nDownloading data by metadata tag:", t['tag'], "\n")
                args.metatag = t['tag']
                args.phasename = t['phasename']
                ext_thelec(args, vasp_db=vasp_db)
                #args.metatag = None
                #args.phasename = None
            else:
                ext_thelec(args,plotfiles=t, vasp_db=vasp_db)

def run_ext_EVfind(subparsers):
    #SUB-PROCESS: EVfind
    pEVfind = subparsers.add_parser("EVfind", help="Find the metadata tags that have 0 K static calculaton finished.")
    pEVfind.add_argument("-w", "--within", dest="within", nargs="?", type=str, default=None,
                      help="find calculations within element list\n"
                           "Default: None")
    pEVfind.add_argument("-all", "--containall", dest="containall", nargs="?", type=str, default=None,
                      help="find calculations must contain all elements in the list\n"
                           "Default: None")
    pEVfind.add_argument("-xall", "--excludeall", dest="excludeall", nargs="?", type=str, default=None,
                      help="exclude calculations that conain all elements in the list\n"
                           "Default: None")
    pEVfind.add_argument("-xany", "--excludeany", dest="excludeany", nargs="?", type=str, default=None,
                      help="exclude calculations that conain any elements in the list\n"
                           "Default: None")
    pEVfind.add_argument("-any", "--containany", dest="containany", nargs="?", type=str, default=None,
                      help="find calculations contain any elements in the list\n"
                           "Default: None")
    pEVfind.add_argument("-v", "--nV", dest="nV", nargs="?", type=int, default=5,
                      help="Return phonon calculations finished for number of volumes larger or equals to. \n"
                           "Default: 5")
    pEVfind.add_argument("-fg", "--findbandgap", dest="findbandgap", action='store_true', default=False,
                      help="report the entries with band gap. \n"
                           "Default: False")
    pEVfind.add_argument("-p", "--print", dest="print", action='store_true', default=False,
                      help="report the entries with band gap. \n"
                           "Default: False")
    pEVfind.add_argument("-plot", "-plot", dest="plot", action='store_true', default=False,
                      help="plot the EV. \n"
                           "Default: False")
    pEVfind.set_defaults(func=ext_EVfind)


def ext_EVfind(args, vasp_db=None):
    """
    find the metadata tag that has finished.

    Parameters
        STR_FOLDER = args.STRUCTURE_FOLDER
            folder/file containing structures
        MATCH_PATTERN = args.MATCH_PATTERN
            Match patterns for structure file, e.g. *POSCAR
        RECURSIVE = args.RECURSIVE
            recursive or not
        WORKFLOW = args.WORKFLOW
            workflow, current only get_wf_gibbs
    """
    if vasp_db is None:
        db_file = loadfn(config_to_dict()["FWORKER_LOC"])["env"]["db_file"]
        vasp_db = VaspCalcDb.from_db_file(db_file, admin=False)
    proc=EVfindMDB(args, vasp_db=vasp_db)
    tags = proc.run_console()
