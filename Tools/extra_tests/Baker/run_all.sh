#!/bin/bash

export MOLCAS_MAXITER=50
export MOLCAS_OUTPUT=NAME

export BASIS="STO-3G"
export OPTS=""
export SCFTYPE=""

for mol in {01..33}*.xyz ; do
  export Project=$(basename ${mol} .xyz)
  echo Running ${Project} ...
  export XYZ=${mol}
  export SCF=${SCFTYPE}

  case "${Project}" in
    *histamine_H+)
      export SCF="${SCF}"$'\n'"Charge = 1"
      ;;
  esac

  pymolcas opt.input -oe ${Project}.output -b 1
done
