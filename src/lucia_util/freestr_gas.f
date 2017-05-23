************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SUBROUTINE FREESTR_GAS
* Deallocate the memory that was set up in MEMSTR_GAS

*
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "mxpdim.fh"
#include "orbinp.fh"
#include "strbas.fh"
#include "csm.fh"
#include "WrkSpc.fh"
#include "cgas.fh"
#include "gasstr.fh"
#include "stinf.fh"
#include "crun.fh"
* allocations during strinf_gas
#include "distsym.fh"
*
*.  Offsets for occupation and reorder array of strings
*
      DO IGRP = 1, NGRP
        NSTRIN = NSTFGP(IGRP)
        LSTRIN = NSTRIN*NELFGP(IGRP)
        CALL GETMEM('OCSTR ','FREE','INTE',KOCSTR(IGRP),LSTRIN)
        CALL GETMEM('STREO ','FREE','INTE',KSTREO(IGRP),NSTRIN)
      END DO
*
*. Number of strings per symmetry and offset for strings of given sym
*. for groups
*
      CALL GETMEM('NSTSGP','FREE','INTE',KNSTSGP(1),NSMST*NGRP)
      CALL GETMEM('ISTSGP','FREE','INTE',KISTSGP(1),NSMST*NGRP)
*
*. Number of strings per symmetry and offset for strings of given sym
*. for types
*
      DO  ITP  = 1, NSTTP
        CALL GETMEM('NSTSO ','FREE','INTE',
     &              KNSTSO(ITP),NSPGPFTP(ITP)*NSMST)
        CALL GETMEM('ISTSO ','FREE','INTE',
     &              KISTSO(ITP),NSPGPFTP(ITP)*NSMST)
      END DO
*
**. Lexical adressing of arrays : use array indeces for complete active space
*
*. Not in use so
      DO  IGRP = 1, NGRP
        CALL GETMEM('Zmat  ','FREE','INTE',KZ(IGRP),NACOB*NELFGP(IGRP))
      END DO
*
*. Mappings between different groups
*
      DO  IGRP = 1, NGRP
        IEL = NELFGP(IGRP)
        IGAS = IGSFGP(IGRP)
        IORB = NOBPT(IGAS)
        ISTRIN = NSTFGP(IGRP)
*. IF creation is involve : Use full orbital notation
*  If only annihilation is involved, compact form will be used
        IF(ISTAC(IGRP,2).NE.0) THEN
          LENGTH = IORB*ISTRIN
          CALL GETMEM('ORBMAP','FREE','INTE',KSTSTM(IGRP,1),LENGTH)
          CALL GETMEM('STRMAP','FREE','INTE',KSTSTM(IGRP,2),LENGTH)
        ELSE IF(ISTAC(IGRP,1).NE.0) THEN
*. Only annihilation map so
          LENGTH = IEL*ISTRIN
          CALL GETMEM('ORBMAP','FREE','INTE',KSTSTM(IGRP,1),LENGTH)
          CALL GETMEM('STRMAP','FREE','INTE',KSTSTM(IGRP,2),LENGTH)
        ELSE
*. Neither annihilation nor creation (?!)
          KSTSTM(IGRP,1) = -1
          KSTSTM(IGRP,2) = -1
        END IF
      END DO
*
*. Symmetry of conjugated orbitals and orbital excitations
*
*     KCOBSM,KNIFSJ,KIFSJ,KIFSJO
      CALL GETMEM('Cobsm ','FREE','INTE',KCOBSM,NACOB)
      CALL GETMEM('Nifsj ','FREE','INTE',KNIFSJ,NACOB*NSMSX)
      CALL GETMEM('Ifsj  ','FREE','INTE',KIFSJ,NACOB**2 )
      CALL GETMEM('Ifsjo ','FREE','INTE',KIFSJO,NACOB*NSMSX)
*
*. Symmetry of excitation connecting  strings of given symmetry
*
      CALL GETMEM('Ststx ','FREE','INTE',KSTSTX,NSMST*NSMST)
*
*. Occupation classes
*
      CALL GETMEM('IOCLS ','FREE','INTE',KIOCLS,NMXOCCLS*NGAS)
*. Annihilation/Creation map of supergroup types
      CALL GETMEM('SPGPAN','FREE','INTE',KSPGPAN,NTSPGP*NGAS)
      CALL GETMEM('SPGPCR','FREE','INTE',KSPGPCR,NTSPGP*NGAS)
*
* Allocated during strinf_gas call
      CALL GETMEM('ISMDFGP','FREE','INTE',ISMDFGP, NSMST*NGRP)
      CALL GETMEM('NACTSYM','FREE','INTE',NACTSYM, NGRP)
      CALL GETMEM('ISMSCR','FREE','INTE',ISMSCR, NGRP)
      RETURN
      END
