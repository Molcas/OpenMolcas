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
Module pt2_guga
#include "pt2_guga.fh"
#ifdef _NEW_
!use definitions, only: iwp, wp
      ! MXLEV should be taken from the gugx module
      integer(kind=iwp), PARAMETER :: MXLEV=100 ,MXL3=(MXLEV*(MXLEV+1))/2

      real(kind=wp) ETA(MXLEV),CITHR,PKPREC

      CHARACTER(Len=8)CLAB10(64)

      integer(kind=iwp)                  LCI,MXCI,                      &
     &     NG1,NG2,NG3,NG3TOT,LG1,LG2,LG3,LF1,LF2,LF3,                  &
     &     IADR10(64,2),IDTAB(MXL3)

      integer(kind=iwp) NPLBUF,IPLBUF,JPLBUF,ISYMA,NSGMA,ISYMB,NSGMB
#endif
End Module pt2_guga
