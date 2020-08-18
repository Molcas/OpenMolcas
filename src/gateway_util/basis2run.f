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
      Subroutine basis2run()
      use Basis_Info
      Implicit None
#include "itmax.fh"
#include "info.fh"
#include "stdalloc.fh"
      integer :: nPrim
      integer :: kExp
      integer :: iPrim, iCnttp, icnt, mdc, jSh, iShSrt, iAng, iBasis
      real*8, allocatable :: primitives(:,:)
      integer, allocatable :: primitive_ids(:,:), IndC(:)
      integer :: iyy, iCo, iAtoms, index_center
*
************************************************************************
*
      iAtoms=0
************************************************************************
*     Generate list of primitive basis functions
*
*     Loop over distinct shell types
      nPrim=0
*     Loop over basis sets
      Do iCnttp = 1, nCnttp
        if (iCnttp.eq.iCnttp_Dummy) cycle
        mdc = mdciCnttp(iCnttp)
        iShSrt = ipVal(iCnttp)
*     Loop over distinct centers
        Do icnt = 1, dbsc(iCnttp)%nCntr
          mdc = mdc + 1
*     Loop over symmetry-related centers
          Do iCo = 0, nIrrep/nStab(mdc)-1
*     Loop over shells associated with this center
*     Start with s type shells
            jSh = iShSrt
            if (Shells(jSh)%Aux.or.
     &          Shells(jSh)%Frag) cycle
            Do iAng = 0, nVal_Shells(iCnttp)-1
              nPrim = nPrim + Shells(jSh)%nExp * Shells(jSh)%nBasis
              jSh = jSh + 1
            End Do
          End Do
        End Do
      End Do

      call put_iScalar('nPrim',nPrim)

      Call mma_allocate(IndC,2*mCentr,label='IndC')
      call mma_allocate(primitive_ids, 3, nPrim,label='primitive_ids')
      call mma_allocate(primitives, 2, nPrim,label='primitives')

*     Loop over distinct shell types
      iPrim  =0
*     Loop over basis sets
      Do iCnttp = 1, nCnttp
        if (iCnttp.eq.iCnttp_Dummy) cycle
        mdc = mdciCnttp(iCnttp)
        iShSrt = ipVal(iCnttp)
*     Loop over distinct centers
        Do icnt = 1, dbsc(iCnttp)%nCntr
          mdc = mdc + 1
*     Loop over symmetry-related centers
          Do iCo = 0, nIrrep/nStab(mdc)-1
*     Loop over shells associated with this center
*     Start with s type shells
            jSh = iShSrt
            if (Shells(jSh)%Aux.or.
     &          ShellS(jSh)%Frag) cycle
*     Get the flat, desymmetrized id of the center
            iyy=Index_Center(mdc,iCo,IndC,iAtoms,mCentr)
            Do iAng = 0, nVal_Shells(iCnttp)-1
*     Pointer to the untouched contraction matrix as after input.
              Do iBasis = 1, Shells(jSh)%nBasis
                Do kExp = 1, Shells(jSh)%nExp
                  iPrim  = iPrim  + 1
                  primitive_ids(1,iPrim) = iyy
                  primitive_ids(2,iPrim) = iAng
                  primitive_ids(3,iPrim) = iBasis
                  primitives(1,iPrim) = Shells(jSh)%Exp(kExp)
                  primitives(2,iPrim) = Shells(jSh)%pCff(kExp,iBasis)
                End Do
              End Do
              jSh = jSh + 1
            End Do
          End Do
        End Do
*
      End Do

      call put_iArray('primitive ids', primitive_ids, 3*nPrim)
      call put_dArray('primitives', primitives, 2*nPrim)

      call mma_deallocate(primitive_ids)
      call mma_deallocate(primitives)
      call mma_deallocate(IndC)
      Return
      End
