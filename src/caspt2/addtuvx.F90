!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
      Subroutine ADDTUVX(NP,NI,NQ,NK,NASHT,iOffP,iOffI,iOffQ,iOffK,     &
     &                   TUVX,nTUVX,PIQK,nPIQK,                         &
     &                   NUMERR)
      use definitions, only: iwp, wp
#ifndef _DEBUGPRINT_
      use Constants, only: One
#endif
      Implicit None
      integer(kind=iwp), intent(in):: NP,NI,NQ,NK,NASHT,iOffP,iOffI,    &
     &                                iOffQ,iOffK
      integer(kind=iwp), intent(in):: nTUVX, nPIQK
      real(kind=wp), intent(in)::  PIQK(nPIQK)
      real(kind=wp), intent(inout):: TUVX(nTUVX)
      integer(kind=iwp), intent(inout):: NUMERR

      integer(kind=iwp) iU, iUVX1, iUVX2, iV, iVX1, iVX2, iX, iX1, iX2
#ifndef _DEBUGPRINT_
#include "warnings.h"
#include "macros.fh"
      unused_var(NUMERR)
#else
      integer(kind=iwp) iPIQK, iT, iTUVX
#endif
!
! Add into correct positions in TUVX:
!
      DO iX=0,NK-1
         iX1=NASHT*(iX+iOffK)
         iX2=NQ   * iX
         DO iV=0,NQ-1
            iVX1=NASHT*(iV+iOffQ+iX1)
            iVX2=NI   *(iV      +iX2)
            DO iU=0,NI-1
               iUVX1=NASHT*(iU+iOffI+iVX1)
               iUVX2=NP   *(iU+      iVX2)
#ifdef _DEBUGPRINT_
               DO iT=1,NP
                  iTUVX=iT+iOffP+iUVX1
                  iPIQK=iT      +iUVX2
! Temporary test statements -- remove after debug!
                  IF(ITUVX.LT.1 .or. ITUVX.gt.NTUVX+1) THEN
                     ITUVX=NTUVX+1
                     NUMERR=NUMERR+1
                     IF (NUMERR.GT.100) THEN
                        WRITE(6,*)' THIS IS TOO MUCH -- STOP.'
                        CALL QUIT(_RC_INTERNAL_ERROR_)
                     END IF
                  END IF
! End of temporary test statements
                  TUVX(iTUVX)=TUVX(iTUVX)+PIQK(iPIQK)
               END DO
#else
               Call DaXpY_(nP,One,PIQK(1+      iUVX2),1,                &
     &                             TUVX(1+iOffP+iUVX1),1)
#endif
            END DO
         END DO
      END DO
!
      End Subroutine ADDTUVX
