=======================================
documentation for cmepy's documentation
=======================================

cmepydoc requires:

 * the `sphinx <http://sphinx.pocoo.org/>`_ python documentation generator

To build the documentation locally, run::

    make html

This will build the html documentation to ``_build/html``.

Publishing the documentation to github
--------------------------------------
There is a (rather dubious) script to build the html documentation and
publish it to github.

Provided:

 * cmepydoc is currently checked out in a git repository
 * the 'origin' remote repo corresponds to some github repository
 * there are no local uncommited changes
    **WARNING : LOCAL UNCOMITTED CHANGES WILL BE DESTROYED**
 * there is no existing ``gh-pages`` branch

then running::

    python publish.py

will rebuild the html documentation and place the resulting output inside a new
branch named 'gh-pages' of the
current repository. This branch will then be pushed to the origin (ie github).

The ``publish.py`` script accepts a few optional command line arguments, to
change the default publish branch, build directory, and remote repository.
