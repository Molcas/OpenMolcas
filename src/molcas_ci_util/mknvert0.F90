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
! Copyright (C) Per Ake Malmqvist                                      *
!               Markus P. Fuelscher                                    *
!***********************************************************************
!#define _DEBUGPRINT_
      Subroutine mknVert0(SGS)
      use struct, only: SGStruct
      use Definitions, only: LF => u6
      Implicit None
      Type(SGStruct) SGS

      Integer IAC

      Associate (iSpin=> SGS%iSpin, nActEl=>SGS%nActEl, nLev=>SGS%nLev,  &
     &           IA0=>SGS%IA0, IB0=>SGS%IB0, IC0=>SGS%IC0, nVert0=>SGS%nVert0)

      IB0=ISPIN-1
      IA0=(NACTEL-IB0)/2
      IC0=NLEV-IA0-IB0

      IF ( ((2*IA0+IB0).NE.NACTEL) .OR.                                 &
     &     (IA0.LT.0) .OR.                                              &
     &     (IB0.LT.0) .OR.                                              &
     &     (IC0.LT.0) ) then
        Write(LF,*)'mknVert0 Error: Impossible specifications.'
        Write(LF,'(1x,a,3I8)')'NACTEL,NLEV,ISPIN:',NACTEL,NLEV,ISPIN
        Write(LF,'(1x,a,3I8)')'IA0,IB0,IC0:      ',IA0,IB0,IC0
        Write(LF,*)' This is a severe internal error, or possibly'
        Write(LF,*)' indicates a strange input which should have been'
        Write(LF,*)' diagnosed earlier. Please submit a bug report.'
        Call Abend()
      End If
      IAC=MIN(IA0,IC0)
      NVERT0=((IA0+1)*(IC0+1)*(2*IB0+IAC+2))/2-(IAC*(IAC+1)*(IAC+2))/6

      End Associate

      End Subroutine mknVert0
