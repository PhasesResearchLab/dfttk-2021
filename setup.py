from setuptools import setup, find_packages
import os
import versioneer

module_dir = os.path.dirname(os.path.abspath(__file__))


def readme():
    with open(os.path.join(module_dir, 'README.rst'), encoding='UTF-8') as f:
        return f.read()


setup(
    name='dfttk',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    packages=find_packages(exclude=["*.tests", "*.tests.*", "tests.*", "tests"]),
    package_data = {'dfttk.structure_builders' : ["prototype_anrl.all", "aflow_prototype_db.json"],
                    'dfttk.analysis' : ['database/*']},
    description='Density functional theory workflows for finite temperature thermodynamics based on atomate workflows. Created by the Phases Research Lab',
    long_description=readme(),
    install_requires=['atomate==1.0.3','phonopy==2.9.1', 'fireworks==2.0.3', 'tinydb==4.4.0', 'ase==3.22.1', 'numpy==1.22.4', 'pymatgen==2022.5.26', 'pymongo==4.1.1', 'maggma==0.47.2', 'monty==2022.4.26', 'custodian==2022.5.26'],
    extras_require={
        'dev': [
            'sphinx',
            'sphinx_rtd_theme',
            'pytest',
            'twine',
        ]
    },
    author='Brandon Bocklund',
    author_email='brandonbocklund@gmail.com',
    url='https://github.com/phasesresearchlab/dfttk',
    license='MIT',
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',

        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Physics',
        'Topic :: Scientific/Engineering :: Chemistry',

        'License :: OSI Approved :: MIT License',

        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8'
    ],
    entry_points={
        'console_scripts': [
            'dfttk = dfttk.scripts.run_dfttk:run_dfttk',
        ]
    }
)
