#! /bin/bash +i

make
source build_envs.bash
echo 'from bash' ${CFLAGS}
echo 'from bash' ${LFLAGS}
py_pfx=`python -c 'import sys; print(sys.prefix)'`
env CFLAGS="${CFLAGS}" LFLAGS="-L${py_pfx}/lib ${LFLAGS}" \
    python quaternion_setup.py install
env CFLAGS="${CFLAGS}" LFLAGS="-L${py_pfx}/lib ${LFLAGS}" \
    python fiasco_numpy_setup.py install
