!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2011, Steven Vancoillie                                *
!***********************************************************************
!***********************************************************************
! Written by Steven Vancoillie, May 2011
! A set of subroutines that can handle RHS arrays in either a serial or
! parallel environment, depending on the situation.
!***********************************************************************
! --> when running serially, the RHS arrays are stored on LUSOLV and are
! loaded into the WORK array when needed.
! --> when running in parallel, the RHS arrays are stored on disk as
! disk resident arrays (DRAs) with filename RHS_XX_XX_XX, where XX is a
! number referring to the case, symmetry, and RHS vector respectively,
! and are loaded onto a global array when needed.
!***********************************************************************

      SUBROUTINE RHS_SAVE (NIN,NIS,lg_W,iCASE,iSYM,iVEC)
      use definitions, only: iwp
!SVC: this routine reads an RHS array in SR format from disk
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
      use definitions, only: u6
#endif
      use caspt2_global, only: LURHS
      use fake_GA, only: GA_Arrays
      use caspt2_module, only: IOFFRHS
      IMPLICIT None
      integer(kind=iwp), intent(in):: NIN,NIS,lg_W,iCASE,iSYM,iVEC

      integer(kind=iwp) IDISK, NW
#ifdef _MOLCAS_MPP_
      integer(kind=iwp) myRank,ISTA,IEND,JSTA,JEND,mpt_W,LDW,NWPROC
#include "global.fh"
#include "mafdecls.fh"

      IF (Is_Real_Par()) THEN
        CALL GA_Sync()
        myRank = GA_NodeID()
        CALL GA_Distribution (lg_W,myRank,ISTA,IEND,JSTA,JEND)
        IF (IEND-ISTA+1.EQ.NIN .AND. ISTA.GT.0) THEN
          CALL GA_Access (lg_W,ISTA,IEND,JSTA,JEND,mpt_W,LDW)
          IF (LDW.NE.NIN) THEN
            WRITE(u6,*) 'RHS_SAVE: Assumption NIN==LDW wrong'
            CALL AbEnd()
          END IF
          NWPROC=NIN*(JEND-JSTA+1)
          IDISK=IOFFRHS(ISYM,ICASE)
          CALL DDAFILE(LURHS(IVEC),1,DBL_MB(mpt_W),NWPROC,IDISK)
          CALL GA_Release (lg_W,ISTA,IEND,JSTA,JEND)
        END IF
        CALL GA_Sync()
      ELSE
#endif
        NW=NIN*NIS
        IDISK=IOFFRHS(ISYM,ICASE)
        CALL DDAFILE(LURHS(IVEC),1,GA_Arrays(lg_W)%A,NW,IDISK)
#ifdef _MOLCAS_MPP_
      END IF
#endif

      END SUBROUTINE RHS_SAVE
