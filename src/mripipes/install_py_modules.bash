#! /bin/bash +i

make
source build_envs.bash
py_pfx=`python -c 'import sys; print(sys.prefix)'`
env CFLAGS="${CFLAGS}" LFLAGS="-L${py_pfx}/lib ${LFLAGS}" \
    SRCFILES="${SRCFILES}" \
    python mripipes_setup.py install

