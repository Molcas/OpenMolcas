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
#include "itmax.fh"
#include "info.fh"
#include "print.fh"
#include "nq_info.fh"
#include "real.fh"
#include "WrkSpc.fh"
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
*     B*s Magnetic * Spin
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
      OperI(1   ) = 255
      OperI(1+1 ) = 255
      OperC(1   ) = 0 ! Dummy
      OperC(1+1 ) = 0 ! Dummy
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
*     A*p
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
      OperI(1   ) = 255
      OperI(1+1 ) = 255
      OperI(1+2 ) = 255
      OperI(1+3 ) = 255
      OperI(1+4 ) = 255
      OperI(1+5 ) = 255
      OperI(1+6 ) = 255
      OperI(1+7 ) = 255
      OperI(1+8 ) = 255
      OperI(1+9 ) = 255
      OperI(1+10) = 255
      OperI(1+11) = 255
      OperC(1   ) = 0 ! Dummy
      OperC(1+1 ) = 0 ! Dummy
      OperC(1+2 ) = 0 ! Dummy
      OperC(1+3 ) = 0 ! Dummy
      OperC(1+4 ) = 0 ! Dummy
      OperC(1+5 ) = 0 ! Dummy
      OperC(1+6 ) = 0 ! Dummy
      OperC(1+7 ) = 0 ! Dummy
      OperC(1+8 ) = 0 ! Dummy
      OperC(1+9 ) = 0 ! Dummy
      OperC(1+10) = 0 ! Dummy
      OperC(1+11) = 0 ! Dummy
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
      OperI(1   ) = 255
      OperI(1+1 ) = 255
      OperC(1   ) = 0 ! Dummy
      OperC(1+1 ) = 0 ! Dummy
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
