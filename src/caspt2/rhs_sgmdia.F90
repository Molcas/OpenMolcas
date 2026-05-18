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

      SUBROUTINE RHS_SGMDIA(NIN,NIS,lg_W,DIN,DIS)
      use definitions, only: iwp, wp
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
#endif
      use fake_GA, only: GA_Arrays
      IMPLICIT None

      integer(kind=iwp), intent(in):: NIN,NIS,lg_W
      real(kind=wp), Intent(in):: DIN(NIN),DIS(NIS)

! Apply the resolvent of the diagonal part of H0 to an RHS array

#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
      integer(kind=iwp) myRank,iLo,iHi,jLo,jHi,NROW,NCOL,mW,LDW

      IF (Is_Real_Par()) THEN
        CALL GA_Sync()
        myRank = GA_NodeID()
!-SVC: get the local vertical stripes of the lg_W vector
        CALL GA_Distribution (lg_W,myRank,iLo,iHi,jLo,jHi)
        IF (iLo.NE.0.AND.jLo.NE.0) THEN
          NROW=iHi-iLo+1
          NCOL=jHi-jLo+1
          CALL GA_Access (lg_W,iLo,iHi,jLo,jHi,mW,LDW)
          CALL SGMDIA(NROW,NCOL,DBL_MB(mW),LDW,DIN(iLo),DIS(jLo))
          CALL GA_Release_Update (lg_W,iLo,iHi,jLo,jHi)
        END IF
        CALL GA_Sync()
      ELSE
#endif
        CALL SGMDIA(NIN,NIS,GA_Arrays(lg_W)%A,NIN,DIN,DIS)
#ifdef _MOLCAS_MPP_
      END IF
#endif

      END SUBROUTINE RHS_SGMDIA
