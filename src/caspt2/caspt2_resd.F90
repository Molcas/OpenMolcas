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
! Copyright (C) 2021, Yoshio Nishimoto                                 *
!***********************************************************************
      SUBROUTINE CASPT2_ResD(Mode,NIN,NIS,lg_W1,lg_W2,DIN,DIS)
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
#endif
      use fake_GA, only: GA_Arrays
      use definitions, only: wp, iwp

      implicit none

      integer(kind=iwp), intent(in) :: Mode, NIN, NIS, lg_W1, lg_W2
      real(kind=wp), intent(in) :: DIN(NIN), DIS(NIS)

#ifdef _MOLCAS_MPP_
      integer(kind=iwp) :: myRank, iLo1, iHi1, jLo1, jHi1, iLo2, iHi2,  &
     &                     jLo2, jHi2, NROW, NCOL, mW1, LDW1, mW2, LDW2
#endif

! Apply the resolvent of the diagonal part of H0 to an RHS array

#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif

#ifdef _MOLCAS_MPP_
      IF (Is_Real_Par()) THEN
        CALL GA_Sync()
        myRank = GA_NodeID()
!-SVC: get the local vertical stripes of the lg_W vector
        CALL GA_Distribution (lg_W1,myRank,iLo1,iHi1,jLo1,jHi1)
        CALL GA_Distribution (lg_W2,myRank,iLo2,iHi2,jLo2,jHi2)
        !! Well, assume the same dimension
        IF (iLo1 > 0 .AND. jLo1 > 0 .AND. iLo2 > 0 .AND. jLo2 > 0) THEN
          NROW=iHi1-iLo1+1
          NCOL=jHi1-jLo1+1
          CALL GA_Access (lg_W1,iLo1,iHi1,jLo1,jHi1,mW1,LDW1)
          CALL GA_Access (lg_W2,iLo2,iHi2,jLo2,jHi2,mW2,LDW2)
          CALL CASPT2_ResD2(MODE,NROW,NCOL,DBL_MB(mW1),DBL_MB(mW2),     &
     &                      LDW1,DIN(iLo1),DIS(jLo1))
          CALL GA_Release_Update (lg_W1,iLo1,iHi1,jLo1,jHi1)
          CALL GA_Release_Update (lg_W2,iLo2,iHi2,jLo2,jHi2)
        END IF
        CALL GA_Sync()
      ELSE
#endif
        CALL CASPT2_ResD2(MODE,NIN,NIS,GA_Arrays(lg_W1)%A,              &
     &                    GA_Arrays(lg_W2)%A,NIN,DIN,DIS)
#ifdef _MOLCAS_MPP_
      END IF
#endif

      END SUBROUTINE CASPT2_ResD
