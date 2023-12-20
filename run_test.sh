#! /bin/sh
set -e
set -v
chmod +x ./quotate.py

pip install -r requirements.txt

./quotate.py tests/H3S-in.cell tests/H3S-out.cell --ibrav 3

./quotate.py tests/SiO2-in.cell tests/SiO2-out.cell --ibrav 2

./quotate.py tests/Cd6Yb-in.poscar tests/Cd6Yb-out.cell --ibrav 3

./quotate.py tests/SiO2-in.cell tests/SiO2.scf.in --ibrav 2 -F espresso-in

./quotate.py tests/S.cell tests/S-out.cell --ibrav 5 
