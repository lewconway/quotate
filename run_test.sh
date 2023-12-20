#! /bin/sh
set -e
set -v
chmod +x ./quotate.py

pip install -r requirements.txt

./quotate.py tests/H3S-in.cell tests/H3S-out.cell --ibrav 3

./quotate.py tests/SiO2-in.cell tests/SiO2-out.cell --ibrav 2

./quotate.py tests/Cd6Yb-in.poscar tests/Cd6Yb-out.cell --ibrav 3

./quotate.py tests/SiO2-in.cell tests/SiO2.scf.in --ibrav 2 -F espresso-in

./quotate.py tests/gnome-3-H6GdIr-R-3m-o.res tests/H6GdIr-out.cell --ibrav 5 
