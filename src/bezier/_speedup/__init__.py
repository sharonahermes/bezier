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

"""Sub-package for performance speed ups."""


import ctypes
import os

import numpy as np


HERE = os.path.abspath(os.path.dirname(__file__))
SHARED_LIB = os.path.join(HERE, 'libspeedup.so')
if os.path.exists(SHARED_LIB):
    LIBSPEEDUP = ctypes.cdll.LoadLibrary(SHARED_LIB)
else:
    LIBSPEEDUP = None
    MSG = '{!r} is missing'.format(SHARED_LIB)
    raise ImportError(MSG)
