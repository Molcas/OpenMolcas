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
! Copyright (C) 1994, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------*
! 1994  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*
      SUBROUTINE MKMAW(SGS)
      use stdalloc, only: mma_allocate
      use struct, only: SGStruct
      IMPLICIT None
      Type (SGStruct) SGS

      Integer IC, ID, ISUM, IU, IV

      Call mma_allocate(SGS%MAW,[1,SGS%nVert],[0,3],Label='SGS%MAW')

      Associate (nVert=>SGS%nVert, MVSta=>SGS%MVSta, MVEnd=>SGS%MVEnd, &
                 iDown=>SGS%Down, iUp=>SGS%Up, iDaw=>SGS%DAW, iRaw=>SGS%Raw, &
                 iMAW=>SGS%MAW)

! COPY LOWER PART OF DIRECT ARC WEIGHT TABLE INTO IMAW:
      DO 210 IV=MVSta,NVERT
        DO 211 IC=0,3
          IMAW(IV,IC)=IDAW(IV,IC)
  211   CONTINUE
  210 CONTINUE
! COPY UPPER PART OF REVERSE ARC WEIGHT TABLE INTO IMAW. HOWEVER,
!    NOTE THAT THE IMAW TABLE IS ACCESSED BY THE UPPER VERTEX.
      DO 230 IU=1,MVSta-1
        DO 220 IC=0,3
          ID=IDOWN(IU,IC)
          IMAW(IU,IC)=0
          IF(ID.NE.0) IMAW(IU,IC)=IRAW(ID,IC)
  220   CONTINUE
  230 CONTINUE
! FINALLY, ADD AN OFFSET TO ARCS LEADING TO MIDLEVELS:
      ISUM=1
      DO IV=MVSta,MVEnd
        DO IC=0,3
          IU=IUP(IV,IC)
          IF(IU.EQ.0) Cycle
          IMAW(IU,IC)=ISUM+IMAW(IU,IC)
        END DO
        ISUM=ISUM+IRAW(IV,4)
      END DO
      DO IV=MVSta,MVEnd
        DO IC=0,3
          IF(IDOWN(IV,IC).EQ.0) Cycle
          IMAW(IV,IC)=ISUM+IMAW(IV,IC)
        END DO
        ISUM=ISUM+IDAW(IV,4)
      END DO
#ifdef _DEBUGPRINT_
 1010 FORMAT(1X,I4,5X,5(1X,I6))
        WRITE(6,*)
        WRITE(6,*)' THE MODIFIED ARC WEIGHT TABLE IN MKMAW:'
        DO 280 IV=1,NVERT
          WRITE(6,1010) IV,(IMAW(IV,IC),IC=0,3)
  280   CONTINUE
        WRITE(6,*)
#endif
      End Associate

      END SUBROUTINE MKMAW
