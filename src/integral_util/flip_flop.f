************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
       Subroutine Flip_Flop(Primitive)
       use Basis_Info
       use Sizes_of_Seward, only:S
       Implicit Real*8 (a-h,o-z)
       Logical Primitive
#include "stdalloc.fh"
#include "Basis_Mode_Parameters.fh"
#include "Basis_Mode.fh"
*
      S%MaxBas(:)=0
      S%MaxPrm(:)=0
*
      Do iCnttp = 1, nCnttp
         nTest = dbsc(iCnttp)%nVal-1
         If (Shells(iCnttp)%Aux .and.
     &       iCnttp.eq.iCnttp_Dummy) nTest=-1
*        Do iCnt = 1, dbsc(iCnttp)%nCntr
*
            Do iAng=0, S%iAngMx
               If (iAng.gt.nTest)  Cycle
               iShll = dbsc(iCnttp)%iVal + iAng
               nExpi=Shells(iShll)%nExp
               If (nExpi.eq.0)   Cycle
               If (Shells(iShll)%nBasis_C.eq.0) Cycle
*
*              Decontract only the ordinary basis sets!
*
               Call mma_deallocate(Shells(iShll)%pCff)
               If (Primitive.and..Not.Shells(iShll)%Aux
     &                      .and..Not.Shells(iShll)%Frag) Then
                  Shells(iShll)%nBasis=nExpi
                  Call mma_allocate(Shells(iShll)%pCff,nExpi,
     &                              Shells(iShll)%nBasis,Label='pCff')
                  Shells(iShll)%pCff(:,:) = Shells(iShll)%Cff_p(:,:,1)
               Else
                  Shells(iShll)%nBasis=Shells(iShll)%nBasis_C
                  Call mma_allocate(Shells(iShll)%pCff,nExpi,
     &                              Shells(iShll)%nBasis,Label='pCff')
                  Shells(iShll)%pCff(:,:) = Shells(iShll)%Cff_c(:,:,1)
               End If
               S%MaxPrm(iAng) = Max(S%MaxPrm(iAng),nExpi)
               S%MaxBas(iAng) = Max(S%MaxBas(iAng),Shells(iShll)%nBasis)
*
            End Do ! iAng
*        End Do    ! iCnt
      End Do       ! iCnttp
*
      Return
      End
!#define _USE_CHECK_CFF_
#ifdef  _USE_CHECK_CFF_
      Subroutine Check_cff(Label)
      Use Basis_Info
      Character*(*) Label
      Write(6,*) Label
      Do i = 1, 30
         Write (6,*) i,
     &               Shells(i)%nExp, Shells(i)%nBasis,
     &               Allocated(Shells(i)%pCff),
     &               Allocated(Shells(i)%Cff_c),
     &               Allocated(Shells(i)%Cff_p)
      End Do
      Return
      End Subroutine Check_cff
#endif

