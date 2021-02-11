Releasing DFTTK
===============

When releasing a new version of DFTTK

1. ``git pull`` to make sure you haven't missed any last-minute commits. **After this point, nothing else is making it into this version.**
#. Ensure that all tests pass locally on develop.
#. ``git push`` and verify all tests pass on all CI services.
#. Generate a list of commits since the last version with ``git --no-pager log --oneline --no-decorate 0.1^..origin/master``
   Replace ``0.1`` with the tag of the last public version.
#. Condense the change list into something user-readable. Update and commit CHANGES.rst with the release date.``
#. ``git tag 0.2 master -m "0.2"`` Replace ``0.2`` with the new version. 
   ``git show 0.2`` to ensure the correct commit was tagged
   ``git push origin master --tags``
#. The new version is tagged in the repository. Now the public package must be built and distributed.

Uploading to PyPI
-----------------

1. ``rm -R dist/*`` on Linux/OSX or ``del dist/*`` on Windows
2. With the commit checked out which was tagged with the new version:
   ``python setup.py sdist``

   **Make sure that the script correctly detected the new version exactly and not a dirty / revised state of the repo.**

   Assuming a correctly configured .pypirc:

   ``twine upload -u bocklund dist/*``

Some useful commands are following
----------------------------------

.. code-block:: bash

    git checkout master
    #
    git pull
    git pull upstream master
    git --version
    git remote set-url --all upstream git@github.com:PhasesResearchLab/dfttk.git
    git remote set-url --add upstream git@github.com:PhasesResearchLab/dfttk.git
    git remote set-url --delete upstream git@github.com:PhasesResearchLab/dfttk.git
    git remote remove upstream
    git remote rm upstream
    git remote -v
    git remote add upstream git@github.com:PhasesResearchLab/dfttk.git
    #
    git tag -d 0.3.0
    git push --delete upstream 0.3.0
    git tag 0.3.0 master -m "0.3.0"
    git push upstream master --tags
    rm dist/*
    python setup.py sdist
    twine upload -u yiwang62 dist/*
    
