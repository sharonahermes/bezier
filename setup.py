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

"""Setup file for bezier."""

import os

import setuptools


VERSION = '0.4.0.dev1'  # Also in codemeta.json
PACKAGE_ROOT = os.path.abspath(os.path.dirname(__file__))

with open(os.path.join(PACKAGE_ROOT, 'README.rst.template')) as file_obj:
    TEMPLATE = file_obj.read()

README = TEMPLATE.format(
    pypi='',
    pypi_img='',
    versions='',
    versions_img='',
    rtd_version=VERSION,
    coveralls_branch=VERSION,
    revision=VERSION,
)

REQUIREMENTS = (
    'numpy >= 1.11.2',
    'six >= 1.9.0',
)
EXTRAS_REQUIRE = {
    ':python_version<"3.4"': ['enum34'],
}
DESCRIPTION = (
    u'Helper for B\u00e9zier Curves, Triangles, and Higher Order Objects')


def setup():
    setuptools.setup(
        name='bezier',
        version=VERSION,
        description=DESCRIPTION,
        author='Danny Hermes',
        author_email='daniel.j.hermes@gmail.com',
        long_description=README,
        scripts=(),
        url='https://github.com/dhermes/bezier',
        packages=setuptools.find_packages('src'),
        package_dir={'': 'src'},
        license='Apache 2.0',
        platforms='Posix; MacOS X; Windows',
        include_package_data=True,
        zip_safe=True,
        install_requires=REQUIREMENTS,
        extras_require=EXTRAS_REQUIRE,
        ext_modules=[],
        classifiers=(
            'Development Status :: 4 - Beta',
            'Intended Audience :: Developers',
            'Intended Audience :: Science/Research',
            'Topic :: Scientific/Engineering :: Mathematics',
            'License :: OSI Approved :: Apache Software License',
            'Operating System :: OS Independent',
            'Programming Language :: Python :: 2',
            'Programming Language :: Python :: 2.7',
            'Programming Language :: Python :: 3',
            'Programming Language :: Python :: 3.5',
            'Programming Language :: Python :: 3.6',
        ),
    )


if __name__ == '__main__':
    setup()
