==================
Contribution Guide
==================

.. |fork| image:: _static/fork_sign.png

.. |pull| image:: _static/pull_request.png


Style Guidelines
----------------

In general, code style should follow PEP8_ and PEP20_. Specifics are summarized below:

- Code should be indented using spaces, not tabs. One indentation = 4 spaces
- Lines longer than 100 should be manually wrapped, but prefer readability
- Minimize blank lines: 2 around top level classes functions, 1 in nested functions
- Workflows, Fireworks, and Firetasks should follow the same naming scheme as in atomate
- Include docstrings for classes and functions (see code), add comments where needed
- Function and variable names should be descriptive (not 'x' or 'xx') and all lowercase_with_underscores
- Class names should be descriptive CapitalWords

.. _PEP8: https://www.python.org/dev/peps/pep-0008/
.. _PEP20: https://www.python.org/dev/peps/pep-0020/

How to Contribute 
-----------------

It is recommended to first install `anaconda <https://docs.anaconda.com/anaconda/install/>`_ if one does not have it installed yet. 

1. Register an account in www.github.com or www.gitlab.com, depending on where the package that you want to contribute to is resided. 
2. Sign into your web account
3. find the the package that you want to contribute to
4. make a fork by click the |fork| sign on top-right corner of the repositary (for instance, see the web site of `DFTTK <https://github.com/yiwang62/dfttk>`_)
5. go back to your web account, you will see the forked repositary shown in your account
6. login in your local machine
7. clone the repository to your local machine and go to/create a folder that you want to reside your contribution, then install the development version by

  .. code-block:: bash

    git clone https://github.com/PhasesResearchLab/dfttk.git
    cd dfttk
    pip install -e .
    dfttk config -mp -aci #a folder named "config" will be created where running environmental info saved

8. Create a new branch for your addition or changes (``git checkout -b mybranchname``), something like

  .. code-block:: bash

    git checkout -b yourhname # "yourhname" is a name that you prefered to use

9. Write codes or make changes and commit them to your branch. (#Note from me, try your best not changing the original codes, only focusing expanding codes will save a lot of troubles when merging your changes to the package)

  .. code-block:: bash

    git commit -am "your notation for the changes"

10. Push your branch to the repository

  .. code-block:: bash

    git push #push your changes to github/gitlab

11. ``*When and only when your are sure that your changes are completely correct*``, submit a pull request by going to your web account and click the sign of |pull| on the top-left side of your repository


After you submit a merge request, other members of the group are able to review your changes and give feedback. Someone with a rank of Master or higher in the project can merge your commits into the master branch.

