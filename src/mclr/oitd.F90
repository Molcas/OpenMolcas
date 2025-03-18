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
      Subroutine OITD(rK,isym,D,Dtmp,act)
!
      use Arrays, only: G1t
      use Constants, only: Zero, One, Two
      use MCLR_Data, only: ipCM, ipMat, nA, nDens2
      use input_mclr, only: nSym,nAsh,nIsh,nOrb
      Implicit None
      Integer iSym
      Real*8 rK(*),D(*),Dtmp(nDens2)

      Logical act
      integer iS, iB, jB, jS
      integer i, j, itri
      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
!
      DTmp(:)=Zero
!
!     Note: even with NAC we set the inactive block,
!     because this is the SA density, not the transition density
      Do iS=1,nSym
        Do iB=1,nIsh(iS)
          Dtmp(1+(ipCM(iS)+(ib-1)*nOrb(iS)+ib-1)-1) = Two
        End Do
      End Do
      If (act) Then
       Do iS=1,nSym
        Do iB=1,nAsh(iS)
         Do jB=1,nAsh(iS)
          Dtmp(1+(ipCM(iS)+ib+nIsh(is)+(jB+nIsh(is)-1)*nOrb(is)-1)-1)=  &
     &    G1t((itri((nA(is)+ib),(nA(is)+jb))))
         End Do
        End Do
       End Do
      End If
!
      Do iS=1,nsym
         jS=ieor(iS-1,isym-1)+1
         If (nOrb(iS)*nOrb(jS).ge.1) Then
            Call DGEMM_('N','T',nOrb(iS),nOrb(jS),nOrb(iS),One,         &
     &                 Dtmp(1+ipCM(iS)-1),nOrb(iS),                     &
     &                 rK(ipMat(jS,iS)),nOrb(jS),                       &
     &                 Zero,D(ipMat(iS,jS)),nOrb(iS))
            Call DGEMM_('T','N',nOrb(iS),nOrb(jS),nOrb(jS),-One,        &
     &                 rK(ipMat(jS,iS)),nOrb(jS),                       &
     &                 Dtmp(1+ipCM(jS)-1),nOrb(jS),                     &
     &                 One,D(ipMat(iS,jS)),nOrb(iS))
         End If
      End Do
      End Subroutine OITD
