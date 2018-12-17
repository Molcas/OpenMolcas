#!/bin/bash

export MOLCAS_MAXITER=50
export MOLCAS_OUTPUT=NAME

export BASIS="3-21G"
export OPTS=""
export SCFTYPE=""
export SLAPAF_OPT="FindTS"

for mol in {01..25}*.xyz ; do
  export Project=$(basename ${mol} .xyz)
  echo Running ${Project} ...
  export XYZ=${mol}
  export SCF=${SCFTYPE}
  export SLAPAF=${SLAPAF_OPT}

  case "${Project}" in
    *CH3O_CH2OH)
      export SCF="${SCF}"$'\n'"UHF"
      ;;
    *cyclopropyl)
      export SCF="${SCF}"$'\n'"UHF"
      ;;
    *beta-formyloxy_ethyl)
      export SCF="${SCF}"$'\n'"UHF"
      ;;
    *H2O_PO3-_H2PO4-)
      export SCF="${SCF}"$'\n'"Charge = -1"
      ;;
    *HCONH3+_NH4+_CO)
      export SCF="${SCF}"$'\n'"Charge = +1"
      ;;
  esac

  case "${Project}" in
    *CH3O_CH2OH)
      export SLAPAF="${SLAPAF}"$'\n'"Gnrm = 0.1"
      ;;
    *beta-formyloxy_ethyl)
      export SLAPAF="${SLAPAF}"$'\n'"Cartesian"
      ;;
    *s-tetrazine_2HCN_N2)
      export SLAPAF="${SLAPAF}"$'\n'"Cartesian"
      ;;
    *HCOCl_HCl_CO)
      export SLAPAF="${SLAPAF}"$'\n'"Gnrm = 0.3"
      ;;
    *H2O_PO3-_H2PO4-)
      export SLAPAF="${SLAPAF}"$'\n'"Gnrm = 0.3"$'\n'"MaxStep = 0.15"
      ;;
    *HCONHOH_HCOHNHO)
      export SLAPAF="${SLAPAF}"$'\n'"Gnrm = 0.1"$'\n'"MaxStep = 0.1"
      ;;
  esac

  if [ -f ${Project}.cons ] ; then
    export SLAPAF="${SLAPAF}"$'\n'"$(< ${Project}.cons)"
  fi

  pymolcas opt.input -oe ${Project}.output -b 1
done
