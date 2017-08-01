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

import os

import distutils.core
import distutils.extension
import Cython.Distutils
import numpy


HERE = os.path.abspath(os.path.dirname(__file__))
SPEEDUP_DIR = os.path.join(HERE, 'src', 'bezier', '_speedup')
SOURCE = os.path.join(SPEEDUP_DIR, 'speedup.pyx')
LIBSPEEDUP = os.path.join(SPEEDUP_DIR, 'speedup.o')


def main():
    npy_include_dir = numpy.get_include()
    ext_modules = [
        distutils.extension.Extension(
            'speedup',
            [SOURCE],
            include_dirs=[npy_include_dir],
            libraries=['gfortran'],
            extra_objects=[
                LIBSPEEDUP,
            ],
        ),
    ]
    distutils.core.setup(
        name='Cython-wrapped Fortran speedups for Bezier.',
        cmdclass={'build_ext': Cython.Distutils.build_ext},
        ext_modules=ext_modules,
   )


if __name__ == '__main__':
    main()
