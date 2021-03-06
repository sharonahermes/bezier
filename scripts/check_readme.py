# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Check that the current README.rst is built from the template."""


from __future__ import print_function

import difflib
import functools
import os
import re


_SCRIPTS_DIR = os.path.dirname(__file__)
_ROOT_DIR = os.path.abspath(os.path.join(_SCRIPTS_DIR, '..'))
TEMPLATE_FILE = os.path.join(_ROOT_DIR, 'README.rst.template')
RELEASE_TEMPLATE_FILE = os.path.join(
    _ROOT_DIR, 'RELEASE_README.rst.template')
INDEX_FILE = os.path.join(_ROOT_DIR, 'docs', 'index.rst')
README_FILE = os.path.join(_ROOT_DIR, 'README.rst')
RTD_VERSION = 'latest'
REVISION = 'master'

PLAIN_CODE_BLOCK = '.. code-block:: python'
SPHINX_CODE_BLOCK1 = """\
.. testsetup:: getting-started

   import sys

   import mock
   import numpy as np

   import bezier

   # Fake the matplotlib/seaborn imports.
   plt_mod = mock.Mock(spec=['figure', 'show'])
   plt_mod.show.return_value = None
   sys.modules['matplotlib.pyplot'] = plt_mod
   mpl_mod = mock.Mock(pyplot=plt_mod, spec=[])
   sys.modules['matplotlib'] = mpl_mod
   seaborn_mod = mock.Mock(spec=['set'])
   seaborn_mod.set.return_value = None
   sys.modules['seaborn'] = seaborn_mod

.. doctest:: getting-started"""
SPHINX_CODE_BLOCK2 = """\
.. doctest:: getting-started
   :options: +NORMALIZE_WHITESPACE"""
SPHINX_CODE_BLOCK3 = '.. doctest:: getting-started'
TEST_CLEANUP = """\
.. testcleanup:: getting-started

   sys.modules.pop('matplotlib')
   sys.modules.pop('matplotlib.pyplot')
   sys.modules.pop('seaborn')

"""
INLINE_MATH_EXPR = re.compile(r':math:`(?P<math>.*?)`')
MOD_EXPR = re.compile(r':mod:`(?P<value>.*) <(?P<module>.*)>`')
DOC_EXPR = re.compile(r':doc:`(?P<value>.*) <(?P<path>.*)>`')
TOCTREE = """\
.. toctree::
   :hidden:
   :maxdepth: 4

   Bezier Package <reference/bezier>
   curve-curve-intersection
   algorithm-helpers
   development

"""
IMG_PREFIX = 'https://cdn.rawgit.com/dhermes/bezier/{revision}/docs/'
EXTRA_LINKS = """\
.. _Curves: https://bezier.readthedocs.io/en/{rtd_version}/reference/bezier.curve.html
.. _Surfaces: https://bezier.readthedocs.io/en/{rtd_version}/reference/bezier.surface.html
.. _Package: https://bezier.readthedocs.io/en/{rtd_version}/reference/bezier.html
.. _DEVELOPMENT doc: https://github.com/dhermes/bezier/blob/{revision}/DEVELOPMENT.rst
"""
BERNSTEIN_BASIS_SPHINX = """\
.. math::

   b_{j, n} = \\binom{n}{j} s^j (1 - s)^{n - j}"""
BERNSTEIN_BASIS_PLAIN = """\
.. image:: {img_prefix}images/bernstein_basis.png
   :align: center"""
BEZIER_DEFN_SPHINX = """\
.. math::

   B(s) = \\sum_{j = 0}^n b_{j, n} \\cdot v_j."""
BEZIER_DEFN_PLAIN = """\
.. image:: {img_prefix}images/bezier_defn.png
   :align: center"""
SUM_TO_UNITY_SPHINX = """\
.. math::

   b_{0, n} + b_{1, n} + \\cdots + b_{n, n} =
       \\left(s + (1 - s)\\right)^n = 1."""
SUM_TO_UNITY_PLAIN = """\
.. image:: {img_prefix}images/sum_to_unity.png
   :align: center"""
DOCS_IMG = """\
.. |docs| image:: https://readthedocs.org/projects/bezier/badge/?version={rtd_version}
   :target: https://bezier.readthedocs.io/en/{rtd_version}/
   :alt: Documentation Status
"""
PYPI_IMG = """
.. |pypi| image:: https://img.shields.io/pypi/v/bezier.svg
   :target: https://pypi.python.org/pypi/bezier
   :alt: PyPI Latest"""
VERSIONS_IMG = """
.. |versions| image:: https://img.shields.io/pypi/pyversions/bezier.svg
   :target: https://pypi.python.org/pypi/bezier
   :alt: Package Versions"""
ZENODO_IMG = """
.. |zenodo| image:: https://zenodo.org/badge/73047402.svg
   :target: https://zenodo.org/badge/latestdoi/73047402
   :alt: Zenodo DOI for ``bezier``"""
JOSS_IMG = """
.. |JOSS| image:: http://joss.theoj.org/papers/10.21105/joss.00267/status.svg
   :target: https://dx.doi.org/10.21105/joss.00267
   :alt: "Journal of Open Source Science" DOI for ``bezier``"""


def inline_math(match):
    """Convert Sphinx inline math to plain reST literal.

    Args:
        match (_sre.SRE_Match): A match (from ``re``) to be used
            in substitution.

    Returns:
        str: The ``match`` converted to literal (i.e. code font).
    """
    return '``{}``'.format(match.group('math'))


def mod_replace(match, sphinx_modules):
    """Convert Sphinx ``:mod:`` to plain reST link.

    Args:
        match (_sre.SRE_Match): A match (from ``re``) to be used
            in substitution.
        sphinx_modules (list): List to be track the modules that have been
            encountered.

    Returns:
        str: The ``match`` converted to a link.
    """
    sphinx_modules.append(match.group('module'))
    return '`{}`_'.format(match.group('value'))


def doc_replace(match, sphinx_docs):
    """Convert Sphinx ``:doc:`` to plain reST link.

    Args:
        match (_sre.SRE_Match): A match (from ``re``) to be used
            in substitution.
        sphinx_docs (list): List to be track the documents that have been
            encountered.

    Returns:
        str: The ``match`` converted to a link.
    """
    sphinx_docs.append(match.group('path'))
    return '`{}`_'.format(match.group('value'))


def get_diff(value1, value2, name1, name2):
    """Get a diff between two strings.

    Args:
        value1 (str): First string to be compared.
        value2 (str): Second string to be compared.
        name1 (str): Name of the first string.
        name2 (str): Name of the second string.

    Returns:
        str: The full diff.
    """
    lines1 = [line + '\n' for line in value1.splitlines()]
    lines2 = [line + '\n' for line in value2.splitlines()]
    diff_lines = difflib.context_diff(
        lines1, lines2, fromfile=name1, tofile=name2)
    return ''.join(diff_lines)


def populate_readme(revision, rtd_version, **extra_kwargs):
    """Populate README template with values.

    Args:
        revision (str): The branch, commit, etc. being referred to (e.g.
            ``master``).
        rtd_version (str): The version to use for RTD (Read the Docs) links
            (e.g. ``latest``).
        extra_kwargs (Dict[str, str]): Over-ride for template arguments.

    Returns:
        str: The populated README contents.

    Raises:
        ValueError: If the ``sphinx_modules`` encountered are not as expected.
        ValueError: If the ``sphinx_docs`` encountered are not as expected.
    """
    with open(TEMPLATE_FILE, 'r') as file_obj:
        template = file_obj.read()

    img_prefix = IMG_PREFIX.format(revision=revision)
    extra_links = EXTRA_LINKS.format(
        rtd_version=rtd_version, revision=revision)
    docs_img = DOCS_IMG.format(rtd_version=rtd_version)
    bernstein_basis = BERNSTEIN_BASIS_PLAIN.format(img_prefix=img_prefix)
    bezier_defn = BEZIER_DEFN_PLAIN.format(img_prefix=img_prefix)
    sum_to_unity = SUM_TO_UNITY_PLAIN.format(img_prefix=img_prefix)

    template_kwargs = {
        'code_block1': PLAIN_CODE_BLOCK,
        'code_block2': PLAIN_CODE_BLOCK,
        'code_block3': PLAIN_CODE_BLOCK,
        'testcleanup': '',
        'toctree': '',
        'bernstein_basis': bernstein_basis,
        'bezier_defn': bezier_defn,
        'sum_to_unity': sum_to_unity,
        'img_prefix': img_prefix,
        'extra_links': extra_links,
        'docs': '|docs| ',
        'docs_img': docs_img,
        'pypi': '\n\n|pypi| ',
        'pypi_img': PYPI_IMG,
        'versions': '|versions|\n\n',
        'versions_img': VERSIONS_IMG,
        'rtd_version': rtd_version,
        'coveralls_branch': 'master',
        'revision': revision,
        'zenodo': '|zenodo|',
        'zenodo_img': ZENODO_IMG,
        'joss': ' |JOSS|',
        'joss_img': JOSS_IMG,
    }
    template_kwargs.update(**extra_kwargs)
    readme_contents = template.format(**template_kwargs)

    # Apply regular expressions to convert Sphinx "roles" to plain reST.
    readme_contents = INLINE_MATH_EXPR.sub(inline_math, readme_contents)

    sphinx_modules = []
    to_replace = functools.partial(
        mod_replace, sphinx_modules=sphinx_modules)
    readme_contents = MOD_EXPR.sub(to_replace, readme_contents)
    if sphinx_modules != ['bezier.curve', 'bezier.surface']:
        raise ValueError('Unexpected sphinx_modules', sphinx_modules)

    sphinx_docs = []
    to_replace = functools.partial(
        doc_replace, sphinx_docs=sphinx_docs)
    readme_contents = DOC_EXPR.sub(to_replace, readme_contents)
    if sphinx_docs != ['reference/bezier', 'development']:
        raise ValueError('Unexpected sphinx_docs', sphinx_docs)

    return readme_contents


def readme_verify():
    """Populate the template and compare to ``README``.

    Raises:
        ValueError: If the current README doesn't agree with the expected
            value computed from the template.
    """
    expected = populate_readme(REVISION, RTD_VERSION)
    # Actually get the stored contents.
    with open(README_FILE, 'r') as file_obj:
        contents = file_obj.read()

    if contents != expected:
        err_msg = '\n' + get_diff(
            contents, expected, 'README.rst.actual', 'README.rst.expected')
        raise ValueError(err_msg)
    else:
        print('README contents are as expected.')


def release_readme_verify():
    """Specialize the template to a PyPI release template.

    Once populated, compare to ``RELEASE_README.rst.template``.

    Raises:
        ValueError: If the current template doesn't agree with the expected
            value specialized from the template.
    """
    version = '{version}'
    expected = populate_readme(
        version,
        version,
        docs='|docs|',
        pypi='',
        pypi_img='',
        versions=' ',
        versions_img='',
        coveralls_branch=version,
        zenodo='',
        zenodo_img='',
        joss='',
        joss_img='',
    )

    with open(RELEASE_TEMPLATE_FILE, 'r') as file_obj:
        contents = file_obj.read()

    if contents != expected:
        err_msg = '\n' + get_diff(
            contents, expected,
            'RELEASE_README.rst.actual',
            'RELEASE_README.rst.expected')
        raise ValueError(err_msg)
    else:
        print('RELEASE_README.rst.template contents are as expected.')


def docs_index_verify():
    """Populate the template and compare to ``docs/index.rst``.

    Raises:
        ValueError: If the current ``index.rst`` doesn't agree with the
            expected value computed from the template.
    """
    with open(TEMPLATE_FILE, 'r') as file_obj:
        template = file_obj.read()

    img_prefix = ''
    extra_links = ''
    docs_img = ''
    expected = template.format(
        code_block1=SPHINX_CODE_BLOCK1,
        code_block2=SPHINX_CODE_BLOCK2,
        code_block3=SPHINX_CODE_BLOCK3,
        testcleanup=TEST_CLEANUP,
        toctree=TOCTREE,
        bernstein_basis=BERNSTEIN_BASIS_SPHINX,
        bezier_defn=BEZIER_DEFN_SPHINX,
        sum_to_unity=SUM_TO_UNITY_SPHINX,
        img_prefix=img_prefix,
        extra_links=extra_links,
        docs='',
        docs_img=docs_img,
        pypi='\n\n|pypi| ',
        pypi_img=PYPI_IMG,
        versions='|versions|\n\n',
        versions_img=VERSIONS_IMG,
        rtd_version=RTD_VERSION,
        coveralls_branch='master',
        revision=REVISION,
        zenodo='|zenodo|',
        zenodo_img=ZENODO_IMG,
        joss=' |JOSS|',
        joss_img=JOSS_IMG,
    )

    with open(INDEX_FILE, 'r') as file_obj:
        contents = file_obj.read()

    if contents != expected:
        err_msg = '\n' + get_diff(
            contents, expected, 'index.rst.actual', 'index.rst.expected')
        raise ValueError(err_msg)
    else:
        print('docs/index.rst contents are as expected.')


def main():
    """Verify specialized versions of ``README.rst.template``."""
    readme_verify()
    release_readme_verify()
    docs_index_verify()


if __name__ == '__main__':
    main()
