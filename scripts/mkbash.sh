#!/bin/bash

if [[ $# -ne 1 ]]; then
    echo "ERROR: script accepts one argument">&2
    exit 2
fi
if [[ -f $1.sh ]]; then
    echo "ERROR: file '$1.sh' already exists"
    exit 2
fi

echo "#!/bin/bash" > $1.sh
chmod u+x $1.sh
vim $1.sh

