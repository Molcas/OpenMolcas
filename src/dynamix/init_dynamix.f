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
!   . |  1    .    2    .    3    .    4    .    5    .    6    .    7 |  .    8
      SUBROUTINE Init_Dynamix
      IMPLICIT REAL*8 (a-h,o-z)
#include "MD.fh"
#include "WrkSpc.fh"
!      CHARACTER*180 Lines(10)
!
!#ifdef _SKIP_
!      Lines(1 )=_MOLCAS_VERSION_
!#ifdef _DEMO_
!      Lines(2 )='DEMO VERSION '
!#else
!      Lines(2 )=' '
!#endif
!      Lines(3 )=' '
!      Lines(4 )='D Y N A M I X'
!      Lines(5 )=' '
!      Lines(6 )='Written by Igor Schapiro'
!      Lines(7 )=' '
!      Lines(8 )='A program for molecular dynamics simulations'
!      Lines(9 )=' '
!      Lines(10)='Code compiled :'//_BUILD_DATE_
!      lLine=Len(Lines(1))
!      CALL Banner(Lines,10,lLine)
!
!      WRITE (6,*)
!      WRITE (6,'(A,I8,A)')   ' The Dynamix program has ',
!     &      mxMem,' double precision words memory available.'
!      WRITE (6,'(A,I8,A)')   '                       ',
!     &      mxMem*8/(1024**2),' MB'
!      WRITE (6,*)
!#endif
!
!     Set the defalut values
!
      CALL Get_nAtoms_Full(natom)
      THERMO=0
      TEMP=298.15
      VELO=0
      POUT=0
      PIN=natom*3
      iPrint=2
      DT=1.0D1
      RESTART=0.0D0
      lH5Restart= .False.
!
!      Etot0=0.0D0
!
      RETURN
!
      END
