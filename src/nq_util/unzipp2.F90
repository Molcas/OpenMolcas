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
      Subroutine UnzipP2(P2Unzip,P2MO,nP2Act)
      use nq_Info

!*****Input
      INTEGER nP2Act
      Real*8,DIMENSION(nP2Act)::P2MO
!*****Output
      Real*8,DIMENSION(NASHT4)::P2Unzip
!*****AUXILIARY
      INTEGER NASHT2,NASHT3,IOFF1,IOff2,IOff3,                          &
     &I,J,K,L,IAct,JAct,kAct,LAct,iIrrep,jIrrep,kIrrep,lIrrep,          &
     &IJ,KL,IJKL
      Real*8 Fact

!***********************************************************************
!                                                                      *
      iTri(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
!                                                                      *
!***********************************************************************

      IF(NASHT4.eq.0) RETURN

      NASHT2=NASHT**2
      NASHT3=NASHT2*NASHT

      DO IIrrep = 0, mIrrep-1
      DO I=1,NASH(iIrrep)
       IAct=iOff_Ash(iIrrep)+I
       IOff1=(IAct-1)*NASHT3
       Do jIrrep = 0, mIrrep-1
       Do J=1,NASH(JIrrep)
        JAct=iOff_Ash(JIrrep)+J
        IOff2=IOff1+(JAct-1)*NASHT2
        IJ=iTri(IAct,JAct)
        do kIrrep = 0, mIrrep-1
        do K=1,NASH(KIrrep)
         KAct=IOff_Ash(KIrrep)+K
         IOff3=IOff2+(KAct-1)*NASHT
         do lIrrep = 0, mIrrep-1
         do L=1,NASH(lIrrep)
          LAct=IOff_Ash(LIrrep)+L
          KL=iTri(KAct,LAct)
          IJKL=iTri(ij,kl)
          Fact=0.5d0
         if((ij.ge.kl).and.(kAct.eq.lAct)) Fact=1.0d0
         if((kl.ge.ij).and.(iAct.eq.jAct)) Fact=1.0d0
          P2Unzip(IOff3+LAct)=P2MO(ijkl)*Fact
         end do
         end do
        end do
        end do
       End Do
       End Do
      END DO
      END DO

      RETURN
      End Subroutine
