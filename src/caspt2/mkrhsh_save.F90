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
! Copyright (C) 1998, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------
! 1998  PER-AAKE MALMQUIST
! DEPARTMENT OF THEORETICAL CHEMISTRY
! UNIVERSITY OF LUND
! SWEDEN
!--------------------------------------------
      SUBROUTINE MKRHS_SAVE(ICASE,ISYM,IVEC,LW)
!SVC: special routine to save the RHS array. MKRHS works in serial, so
! in case of a true parallel run we need to put the local array in a
! global array and then save that to disk in a distributed fashion.
      use definitions, only: iwp
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
      use fake_GA, only: GA_Arrays
#endif
      use caspt2_module, only: NASUP,NISUP
      IMPLICIT None

      integer(kind=iwp), intent(in):: ICASE,ISYM,IVEC,LW

      integer(kind=iwp) NAS, NIS, lg_w

      NAS=NASUP(ISYM,ICASE)
      NIS=NISUP(ISYM,ICASE)

#ifdef _MOLCAS_MPP_
      IF (IS_REAL_PAR()) THEN
        CALL RHS_ALLO(NAS,NIS,lg_W)
        CALL RHS_PUT(NAS,NIS,lg_W,GA_Arrays(LW)%A)
      ELSE
#endif
        lg_W=LW
#ifdef _MOLCAS_MPP_
      END IF
#endif

      CALL RHS_SAVE(NAS,NIS,lg_W,iCASE,iSYM,iVEC)

#ifdef _MOLCAS_MPP_
      IF (IS_REAL_PAR()) THEN
        CALL RHS_FREE(lg_W)
      END IF
#endif
      END SUBROUTINE MKRHS_SAVE
