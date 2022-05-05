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
! Copyright (C) 2021, Jie J. Bao                                       *
!***********************************************************************
! ****************************************************************
! history:                                                       *
! Jie J. Bao, on Dec. 08, 2021, created this file.               *
! ****************************************************************
      Subroutine PackPot1(Packed,Full,nPack,Factor)
      use nq_Info

!*****Input
      Real*8 Factor
      Real*8,DIMENSION(NPot1)::Full
!*****Output
      Real*8,DIMENSION(nPack)::Packed
!*****Auxiliary
      INTEGER iIrrep,p,q,iOff1,IOff2,nOrbs
      DO iIrrep=0,mIrrep-1
       nOrbs=mOrb(iIrrep)
       IOff1=OffOrbTri(iIrrep)
       IOff2=OffOrb2(iIrrep)
       Do P=1,nOrbs
        do Q=1,P
      Packed(IOff1+(P-1)*P/2+Q)=                                        &
     &Full(IOff2+(P-1)*nOrbs+Q)+Full(IOff2+(Q-1)*nOrbs+P)
        end do
       End Do
      END DO
      CALL DScal_(nPack,Factor,Packed,1)
      RETURN
      End Subroutine
