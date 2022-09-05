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
      Subroutine Get_Auxiliary_Shells(iSO,nSO,jOff,iSO2Shl,nSO2Shl,     &
     &                                iPair,nPair)
      Integer iSO(2,nSO), iSO2Shl(nSO2Shl), iPair(nPair)
!
!      Write (6,*) 'iSO'
!      Write (6,*) '==='
!      Do i = 1, nSO
!         Write (6,*) iSO(1,i), iSO(2,i)
!      End Do
!
!      Write (6,*) 'iSO2Shl'
!      Write (6,*) '======='
!      Do i = 1, nSO2Shl
!         Write (6,*) i, iSO2Shl(i)
!      End Do
       Do i = 1, nSO
          k=iSO(1,i) + jOff
          l=iSO(2,i) + jOff
          kSh=iSO2Shl(k)
          lSh=iSO2Shl(l)
!         Write (6,*) 'k,kSh=',k,kSh
!         Write (6,*) 'l,lSh=',k,lSh
          kl = Max(kSh,lSh)*(Max(kSh,lSh)-1)/2 + Min(kSh,lSh)
          iPair(kl)=1
      End Do
!     Write (6,*) 'iPairs'
!     Write (6,*) '======'
!     Do i = 1, nPair
!        Write (6,*) iPair(i)
!     End Do
!
      Return
      End
