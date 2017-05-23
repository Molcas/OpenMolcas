************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1994, Jeppe Olsen                                      *
************************************************************************
      SUBROUTINE MEMSTR_GAS
*
*
* Construct pointers for saving information about strings and
* their mappings
*
* GAS version
*
*========
* Input :
*========
* Number and groups of strings defined by /GASSTR/
* Symmetry information stored in         /CSM/
* String information stored in           /STINF/
*=========
* Output
*=========
* Pointers stored in common block /STRBAS/
*
* Jeppe Olsen , Winter of 1994
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
* Some dummy initializtions
      NTEST = 0
*
*. Start of string information
*
*.  Offsets for occupation and reorder array of strings
*
      DO IGRP = 1, NGRP
        NSTRIN = NSTFGP(IGRP)
        LSTRIN = NSTRIN*NELFGP(IGRP)
        CALL GETMEM('OCSTR ','ALLO','INTE',KOCSTR(IGRP),LSTRIN)
        CALL GETMEM('STREO ','ALLO','INTE',KSTREO(IGRP),NSTRIN)
      END DO
*
*. Number of strings per symmetry and offset for strings of given sym
*. for groups
*
      CALL GETMEM('NSTSGP','ALLO','INTE',KNSTSGP(1),NSMST*NGRP)
      CALL GETMEM('ISTSGP','ALLO','INTE',KISTSGP(1),NSMST*NGRP)
*
*. Number of strings per symmetry and offset for strings of given sym
*. for types
*
      DO  ITP  = 1, NSTTP
        CALL GETMEM('NSTSO ','ALLO','INTE',
     &              KNSTSO(ITP),NSPGPFTP(ITP)*NSMST)
        CALL GETMEM('ISTSO ','ALLO','INTE',
     &              KISTSO(ITP),NSPGPFTP(ITP)*NSMST)
      END DO
*
**. Lexical adressing of arrays : use array indeces for complete active space
*
*. Not in use so
      DO  IGRP = 1, NGRP
        CALL GETMEM('Zmat  ','ALLO','INTE',KZ(IGRP),NACOB*NELFGP(IGRP))
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
          CALL GETMEM('ORBMAP','ALLO','INTE',KSTSTM(IGRP,1),LENGTH)
          CALL GETMEM('STRMAP','ALLO','INTE',KSTSTM(IGRP,2),LENGTH)
        ELSE IF(ISTAC(IGRP,1).NE.0) THEN
*. Only annihilation map so
          LENGTH = IEL*ISTRIN
          CALL GETMEM('ORBMAP','ALLO','INTE',KSTSTM(IGRP,1),LENGTH)
          CALL GETMEM('STRMAP','ALLO','INTE',KSTSTM(IGRP,2),LENGTH)
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
      CALL GETMEM('Cobsm ','ALLO','INTE',KCOBSM,NACOB)
      CALL GETMEM('Nifsj ','ALLO','INTE',KNIFSJ,NACOB*NSMSX)
      CALL GETMEM('Ifsj  ','ALLO','INTE',KIFSJ,NACOB**2 )
      CALL GETMEM('Ifsjo ','ALLO','INTE',KIFSJO,NACOB*NSMSX)
*
*. Symmetry of excitation connecting  strings of given symmetry
*
      CALL GETMEM('Ststx ','ALLO','INTE',KSTSTX,NSMST*NSMST)
*
*. Occupation classes
*
      CALL GETMEM('IOCLS ','ALLO','INTE',KIOCLS,NMXOCCLS*NGAS)
*. Annihilation/Creation map of supergroup types
      CALL GETMEM('SPGPAN','ALLO','INTE',KSPGPAN,NTSPGP*NGAS)
      CALL GETMEM('SPGPCR','ALLO','INTE',KSPGPCR,NTSPGP*NGAS)
*
*. Last word of string information
*
*
      RETURN
      END
