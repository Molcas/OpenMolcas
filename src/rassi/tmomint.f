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
      SubRoutine TMOMInt(wavevector,iOpt)
************************************************************************
*                                                                      *
* Object: driver for computation of TMOM integrals                     *
*                                                                      *
************************************************************************
      Use MpmC
      Implicit Real*8 (A-H,O-Z)
      External EMFInt, EMFMem
*     ipList: list of pointers to the integrals of each component
*             of the operator
*     OperI: list which irreps a particular component of the operator
*            belongs to
*     OperC: list the character of each component of the operator
*     CoorO: list of origins of the operator, one for each component
      Integer, Dimension(:), Allocatable :: ipList, OperI, OperC
      Real*8, Dimension(:), Allocatable :: CoorO, Nuc
      Real*8 wavevector(3)
!define _DEBUG_
#ifdef _DEBUG_
#include "stdalloc.fh"
      Real*8, Dimension(:), Allocatable :: Int_R, Int_I, Temp_Int
#endif
#include "itmax.fh"
#include "info.fh"
#include "print.fh"
#include "nq_info.fh"
#include "real.fh"
#include "wldata.fh"
#include "property_label.fh"
#include "oneswi.fh"
#include "warnings.fh"
      Character*8 Label
      Dimension dum(1),idum(1)
*
      Call Set_Basis_Mode('Valence')
      Call Setup_iSD()
*
************************************************************************
************************************************************************
*                                                                      *
*     Electromagnetic field radiation integrals.                       *
*                                                                      *
*     Note that the integral is not symmetric or antisymmetric!        *
*                                                                      *
************************************************************************
************************************************************************
      rHrmt=-One ! Not used
*
*     B*s Magnetic * Spin: not that this boils down to just integrals
*     over A.
*
      If (iOpt.eq.2) Then
      nOrdOp = 0
      Label='TMOM0'
      nComp = 2
      Call Allocate_Aux()
*     Here we put in the k-vector
      Call FZero(CoorO,3*nComp)
      Call dcopy_(3,wavevector,1,CoorO,1)
*
*     The electromagnetic field operator contributes to all
*     irreducible irreps, hence OperI=255. Since the operator
*     itself is not symmetry-adapted OperC is set to a dummy value.
*
      OperI(:) = 255
      OperC(:) = 0 ! Dummy
*
      Call dcopy_(nComp,[Zero],0,Nuc,1)
      Call OneEl(EMFInt,EMFMem,Label,ipList,OperI,nComp,
     &           CoorO,nOrdOp,Nuc,rHrmt,OperC,
     &           dum,1,dum,idum,0,0,
     &           dum,1,0)
#ifdef _DEBUG_
*
*     This section of the code is for pure debugging and will replace
*     exact operator with truncated expansions of the operator in
*     terms of multipole integrals
*
      iOpt0=0
      iOpt1=1
      iOpt2=2
      iRc=-1
      Label='TMOM0'
      iComp=1
*     Pick up the size and the symmetry label.
      Call iRdOne(iRc,iOpt1,Label,iComp,idum,iSyLbl)
      nInts=idum(1)
*
      Call mma_allocate(Int_R,nInts+4,Label='Int_R')
      Call mma_allocate(Int_I,nInts+4,Label='Int_I')
      Call mma_allocate(Temp_Int,nInts+4,Label='Temp_Int')
*
      Int_R=0.0D0
      Int_R(nInts+1:nInts+3)=CoorO
      Int_I=0.0D0
      Int_I(nInts+1:nInts+3)=CoorO
      Temp_Int=0.0D0
*
      nMltpl=2
      iCase=1
      Phase=1.0D0
      Do iMltpl= 0, nMltpl
         Write (Label,'(A,I2)') 'Mltpl ',iMltpl
         nComp=(iMltpl+1)*(iMltpl+2)/2
         iComp=0
         xyz=1.0D0
         Do ix = iMltpl, 0, -1
            xyz=xyz*CoorO(1)**ix
            Do iy = iMltpl-ix, 0, -1
               xyz=xyz*CoorO(2)**iy
               iz = iMltpl-ix-iy
               xyz=xyz*CoorO(3)**iz
*
               iComp=iComp+1
               Call RdOne(iRc,iOpt2,Label,iComp,Temp_Int,iSyLbl)
*
               Fact=Phase*Gamma(Dble(iMltpl)+1.0D0)/
     &                   (Gamma(Dble(ix)+1.0D0)
     &                   *Gamma(Dble(iy)+1.0D0)
     &                   *Gamma(Dble(iz)+1.0D0))
*
               If (iCase.eq.1) Then
*                 Contribution to the real part
                  Call DaXpY_(nInts,Fact,Temp_Int,1,Int_R,1)
               Else
*                 Contribution to the imaginary part
                  Call DaXpY_(nInts,Fact,Temp_Int,1,Int_I,1)
               End If
            End Do
         End Do
*
         If (iCase.eq.1) Then
            iCase=2
         Else
            iCase=1
            Phase=-Phase
         End If
      End Do
*
*     Overwrite the integrals with a truncated expansion,
*
      Label='TMOM0'
      iComp=1
      Call WrOne(iRc,iOpt0,Label,iComp,Int_R,iSyLbl)
      iComp=2
      Call WrOne(iRc,iOpt0,Label,iComp,Int_I,iSyLbl)
*
      Call mma_deallocate(Int_R)
      Call mma_deallocate(Int_I)
      Call mma_deallocate(Temp_Int)
*
#endif
*
      Call Deallocate_Aux()
      End If
*
*     A*nabla. Note that when used the numbers are multiplied with -i to
*     generate A*p.
*
      If (iOpt.le.2) Then
      nOrdOp = 1
      Label='TMOM'
      nComp = 12
      Call Allocate_Aux()
*     Here we put in the k-vector
      Call FZero(CoorO,3*nComp)
      Call dcopy_(3,wavevector,1,CoorO,1)
*
*     The electromagnetic field operator contributes to all
*     irreducible irreps, hence OperI=255. Since the operator
*     itself is not symmetry-adapted OperC is set to a dummy value.
*
      OperI(:) = 255
      OperC(:) = 0 ! Dummy
*
      Call dcopy_(nComp,[Zero],0,Nuc,1)
      Call OneEl(EMFInt,EMFMem,Label,ipList,OperI,nComp,
     &           CoorO,nOrdOp,Nuc,rHrmt,OperC,
     &           dum,1,dum,idum,0,0,
     &           dum,1,0)
*
      Call Deallocate_Aux()
      End If
*
*     The A^2 term
*
      If (iOpt.gt.2) Then
      nOrdOp = 0
      Label='TMOM2'
      nComp = 2
      Call Allocate_Aux()
*     Here we put in the k-vector
      Call FZero(CoorO,3*nComp)
      Call dcopy_(3,wavevector,1,CoorO,1)
*     Change the argument to 2xA
      Call dscal_(3,2.0D0,CoorO,1)
*
*     The electromagnetic field operator contributes to all
*     irreducible irreps, hence OperI=255. Since the operator
*     itself is not symmetry-adapted OperC is set to a dummy value.
*
      OperI(:) = 255
      OperC(:) = 0 ! Dummy
*
      Call dcopy_(nComp,[Zero],0,Nuc,1)
      Call OneEl(EMFInt,EMFMem,Label,ipList,OperI,nComp,
     &           CoorO,nOrdOp,Nuc,rHrmt,OperC,
     &           dum,1,dum,idum,0,0,
     &           dum,1,0)
*
      Call Deallocate_Aux()
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Call Free_iSD()
*                                                                      *
************************************************************************
*                                                                      *
      Return
*
      Contains
      Subroutine Allocate_Aux()
      Implicit None
#include "stdalloc.fh"
*
      Call mma_Allocate(ipList,nComp,Label='ipList')
      Call mma_Allocate(OperI,nComp,Label='OperI')
      Call mma_Allocate(OperC,nComp,Label='OperC')
      Call mma_Allocate(CoorO,3*nComp,Label='CoorO')
      Call mma_Allocate(Nuc,nComp,Label='Nuc')
*
      Return
      End Subroutine Allocate_Aux
      Subroutine Deallocate_Aux()
      Implicit None
#include "stdalloc.fh"
*
      Call mma_Deallocate(OperC)
      Call mma_Deallocate(OperI)
      Call mma_Deallocate(ipList)
      Call mma_Deallocate(CoorO)
      Call mma_Deallocate(Nuc)
*
      Return
      End Subroutine Deallocate_Aux
*
      End Subroutine TMOMInt
