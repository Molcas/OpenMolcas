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
       Implicit Real*8 (a-h,o-z)
       Logical Primitive
#include "itmax.fh"
#include "info.fh"
#include "stdalloc.fh"
#include "Basis_Mode_Parameters.fh"
#include "Basis_Mode.fh"
*
      MaxBas(:)=0
      MaxPrm(:)=0
*
      Do iCnttp = 1, nCnttp
         nTest = nVal_Shells(iCnttp)-1
         If (AuxShell(iCnttp) .and.
     &       iCnttp.eq.iCnttp_Dummy) nTest=-1
*        Do iCnt = 1, dbsc(iCnttp)%nCntr
*
            Do iAng=0, iAngMx
               If (iAng.gt.nTest)  Cycle
               iShll = ipVal(iCnttp) + iAng
               nExpi=Shells(iShll)%nExp
               If (nExpi.eq.0)   Cycle
               If (nBasis_Cntrct(iShll).eq.0) Cycle
*
*              Decontract only the ordinary basis sets!
*
               Call mma_deallocate(Shells(iShll)%pCff)
               If (Primitive.and..Not.AuxShell(iShll)
     &                      .and..Not.FragShell(iShll)) Then
                  nBasis(iShll)=nExpi
                  Call mma_allocate(Shells(iShll)%pCff,nExpi,
     &                              nBasis(iShll),Label='pCff')
                  Shells(iShll)%pCff(:,:) = Shells(iShll)%Cff_p(:,:,1)
               Else
                  nBasis(iShll)=nBasis_Cntrct(iShll)
                  Call mma_allocate(Shells(iShll)%pCff,nExpi,
     &                              nBasis(iShll),Label='pCff')
                  Shells(iShll)%pCff(:,:) = Shells(iShll)%Cff_c(:,:,1)
               End If
               MaxPrm(iAng) = Max(MaxPrm(iAng),nExpi)
               MaxBas(iAng) = Max(MaxBas(iAng),nBasis(iShll))
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

