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
! Copyright (C) 1992, Roland Lindh                                     *
!***********************************************************************
      Subroutine Get_D1A(CMO,D1A_MO,D1A_AO,                             &
     &                    nsym,nbas,nish,nash,ndens)
      use Constants, only: Zero, One
      use stdalloc, only: mma_allocate, mma_deallocate
      Implicit None
      Integer nSym, nDens
      Real*8 CMO(*) , D1A_MO(*) , D1A_AO(*)
      Integer nbas(nsym),nish(nsym),nash(nsym)

      Integer i, j, iTri
      Real*8, Allocatable:: Scr1(:), Tmp1(:,:), Tmp2(:,:)
      Integer iOff2, iOff3, ii, iSym, iBas, iAsh, iIsh

      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)


      iOff2 = 1
      iOff3 = 1
      ii=0
      Call mma_allocate(Scr1,2*nDens,Label='Scr1')
      Do iSym = 1,nSym
        iBas = nBas(iSym)
        iAsh = nAsh(iSym)
        iIsh = nIsh(iSym)
        Call dCopy_(iBas*iBas,[Zero],0,Scr1(iOff3),1)
        If ( iAsh.ne.0 ) then
          Call mma_allocate(Tmp1,iAsh,iAsh,Label='Tmp1')
          Call mma_allocate(Tmp2,iBas,iAsh,Label='Tmp2')
          Do i=1,iAsh
           Do j=1,iAsh
            Tmp1(j,i)=D1A_MO(itri(i+ii,j+ii))
           End do
          End do
          ii=ii+iash
          Call DGEMM_('N','T',                                          &
     &                iBas,iAsh,iAsh,                                   &
     &                One,CMO(iOff2+iIsh*iBas),iBas,                    &
     &                    Tmp1,iAsh,                                    &
     &               Zero,Tmp2,iBas)
          Call DGEMM_('N','T',                                          &
     &                iBas,iBas,iAsh,                                   &
     &                One,Tmp2,iBas,                                    &
     &                    CMO(iOff2+iIsh*iBas),iBas,                    &
     &               Zero,Scr1(iOff3),iBas)
          Call mma_deallocate(Tmp2)
          Call mma_deallocate(Tmp1)
        End If
        iOff2 = iOff2 + iBas*iBas
        iOff3 = iOff3 + iBas*iBas
      End Do
      Call Fold2(nSym,nBas,Scr1,D1A_AO)
      Call mma_deallocate(Scr1)
      Return
      End Subroutine Get_D1A
