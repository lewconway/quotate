#! /bin/sh
set -e
set -v
chmod +x ./quotate.py

pip install -r requirements.txt

./quotate.py tests/H3S-in.cell tests/H3S-out.cell --ibrav 3

./quotate.py tests/SiO2-in.cell tests/SiO2-out.cell --ibrav 2
