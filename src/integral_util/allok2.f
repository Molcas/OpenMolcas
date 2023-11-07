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
!               1995, Martin Schuetz                                   *
!***********************************************************************
!#define _DEBUGPRINT_
      SubRoutine AlloK2()
!***********************************************************************
!                                                                      *
!  Object: Allocate space for K2 entities.                             *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, Sweden. November '92                 *
!             Martin Schuetz, Dept. of Theoretical Chemistry,          *
!             University of Lund, Sweden. Jun '95                      *
!***********************************************************************
      use k2_arrays, only: nDeDe, MaxDe, DoGrad_, DoHess_
      use iSD_data, only: iSD
      use Basis_Info, only: Shells
      use Sizes_of_Seward, only: S
      use Symmetry_Info, only: nIrrep
      use stdalloc, only: mma_allocate
      use k2_structure, only: k2data, Allocate_k2data, ZZZ_r, ZZZ_i,
     &                        k2_Processed, nIndK2, IndK2
      Implicit None

      Integer i, ixyz, nElem, nabSz, Nr_of_Densities, iS, nSkal, iShll,
     &        iAng, iCmp, iBas, iPrim, iAO, iShell, jS, jShll, jAng,
     &        jCmp, jBas, jPrim, jAO, jShell, iDeSiz, iSMLbl, nSO,
     &        nZeta, ijCmp, nHm, iIrrep, ik2, j, iTri,
     &        ijS, nk2_real, nData, nk2_integer
      Integer, external:: MemSO1
!
!---- Statement function
!
      nElem(i)=(i+1)*(i+2)/2
      nabSz(ixyz) = (ixyz+1)*(ixyz+2)*(ixyz+3)/6  - 1
      iTri(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
!
!
#ifdef _DEBUGPRINT_
      If (Allocated(k2Data).and.k2_processed) Then
         Write (6,*) 'Enter Allok2, k2_Status=Processed'
      Else If (Allocated(k2Data)) Then
         Write (6,*) 'Enter Allok2, k2_Status=Active'
      Else
         Write (6,*) 'Enter Allok2, k2_Status=InActive'
      End If
#endif
      If (Allocated(k2Data) .or. k2_processed) Return
!
      Call Nr_Shells(nSkal)

      ik2 = 0
      nk2_real=0
      nk2_integer=0
      Do iS = 1, nSkal
         iShll  = iSD( 0,iS)
         If (Shells(iShll)%Aux .and. iS.ne.nSkal) Cycle
         iPrim  = iSD( 5,iS)
         iAng   = iSD( 1,iS)
         iCmp   = iSD( 2,iS)
         Do jS = 1, iS
            jShll  = iSD( 0,jS)
            If (Shells(iShll)%Aux.and..Not.Shells(jShll)%Aux) Cycle
            If (Shells(jShll)%Aux .and. jS.eq.nSkal) Cycle
            jPrim  = iSD( 5,jS)
            jAng   = iSD( 1,jS)
            jCmp   = iSD( 2,jS)

            nZeta=iPrim*jPrim
            ijCmp=0
            If (DoGrad_) ijCmp=nElem(iAng)*nElem(jAng)
            nHm=iCmp*jCmp*(nabSz(iAng+jAng)-nabSz(Max(iAng,jAng)-1))
            If (DoHess_) nHm=0

            ik2 = ik2 + 1

            nData = nZeta * (10 + ijCmp*2) + nHm*nIrrep
            nk2_real = nk2_real + nData*nIrrep
            nk2_integer = nk2_integer + (nZeta+1)*nIrrep

         End Do
      End Do

      Call mma_allocate(ZZZ_r,nk2_real   ,Label='ZZZ_r')
      Call mma_allocate(ZZZ_i,nk2_integer,Label='ZZZ_i')
      Allocate(k2Data(1:nIrrep,1:ik2))
!
!     determine memory size for K2 entities
!     for this, run dummy K2 loop...
      ik2 = 0
      nDeDe = 0
      MaxDe = 0
      nr_of_Densities=1 ! Hardwired option

      nIndk2=S%nShlls*(S%nShlls+1)/2
      call mma_allocate(Indk2,3,nIndk2,Label='Indk2')
      Indk2(:,:)=0
!                                                                      *
!***********************************************************************
!                                                                      *
!-----Double loop over shells. These loops decide the integral type
!
      Do iS = 1, nSkal
         iShll  = iSD( 0,iS)
         If (Shells(iShll)%Aux .and. iS.ne.nSkal) Cycle
         iAng   = iSD( 1,iS)
         iCmp   = iSD( 2,iS)
         iBas   = iSD( 3,iS)
         iPrim  = iSD( 5,iS)
         iAO    = iSD( 7,iS)
         iShell = iSD(11,iS)
!
         Do jS = 1, iS
            jShll  = iSD( 0,jS)
            If (Shells(iShll)%Aux.and..Not.Shells(jShll)%Aux) Cycle
            If (Shells(jShll)%Aux .and. jS.eq.nSkal) Cycle
            jAng   = iSD( 1,jS)
            jCmp   = iSD( 2,jS)
            jBas   = iSD( 3,jS)
            jPrim  = iSD( 5,jS)
            jAO    = iSD( 7,jS)
            jShell = iSD(11,jS)

!
            If (nIrrep==1) Then
               iDeSiz = 1 + iPrim*jPrim +               iCmp*jCmp
            Else
               iDeSiz = 1 + iPrim*jPrim + (iBas*jBas+1)*iCmp*jCmp
            End If

            MaxDe = Max(MaxDe,iDeSiz)
            iSmLbl = 1
            nSO = MemSO1(iSmLbl,iCmp,jCmp,iShell,jShell,iAO,jAO)
            If (nSO.gt.0) Then
               nDeDe = nDeDe + nr_of_Densities*iDeSiz*nIrrep
            End If
!
            nZeta=iPrim*jPrim
            ijCmp=0
            If (DoGrad_) ijCmp=nElem(iAng)*nElem(jAng)
            nHm=iCmp*jCmp*(nabSz(iAng+jAng)-nabSz(Max(iAng,jAng)-1))
            If (DoHess_) nHm=0

            ijS=iTri(iShell,jShell)
            ik2=ik2+1
            Indk2(3,ijS)=ik2
            Do iIrrep = 1, nIrrep
               Call Allocate_k2data(k2data(iIrrep,ik2),nZeta,ijCmp,nHm)
            End Do

         End Do
      End Do
!
      Return
      End SubRoutine AlloK2
