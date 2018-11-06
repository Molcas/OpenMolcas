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
      SubRoutine TMOSInt(wavevector)
************************************************************************
*                                                                      *
* Object: driver for computation of TMOS integrals                     *
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
      Real*8, Dimension(:), Allocatable :: CoorO
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
*
      Call Set_Basis_Mode('Valence')
      Call Setup_iSD()
*     Write (*,*) 'wavevector=',wavevector
*
************************************************************************
************************************************************************
*                                                                      *
*     Electromagnetic field radation integrals.                        *
*                                                                      *
*     Note that the integral is not symmetric or antisymmetric!        *
*                                                                      *
************************************************************************
************************************************************************
      rHrmt=-One ! Note used
*
*     B*s Magnetic * Spin
*
      nOrdOp = 0
      Label='TMOS0'
      nComp = 2
      Call Allocate_Auxiliary()
      Call GetMem('Nuc   ','ALLO','REAL',ipNuc,nComp)
*     Here we put in the k-vector
      Call FZero(CoorO,3*nComp)
      Call dcopy_(3,wavevector,1,CoorO,1)
*
*     The electromagnetic field operator contributes to all
*     irreducible irreps, hence OperI=255. Since the operator
*     it self is not symmetry adopted OperC is set to a dummy value.
*
      OperI(1   ) = 255
      OperI(1+1 ) = 255
      OperC(1   ) = 0 ! Dummy
      OperC(1+1 ) = 0 ! Dummy
*
      Call dcopy_(nComp,Zero,0,Work(ipNuc),1)
*     Write (*,*) 'Here we go!'
      Call OneEl(EMFInt,EMFMem,Label,ipList,OperI,nComp,
     &           CoorO,nOrdOp,Work(ipNuc),rHrmt,OperC,
     &           dum,1,dum,idum,0,0,
     &           dum,1,0)
*
      Call GetMem('Nuc   ','FREE','REAL',ipNuc,nComp)
      Call Deallocate_Auxiliary()
*
*     A*p
*
      nOrdOp = 1
      Label='TMOS'
      nComp = 12
      Call Allocate_Auxiliary()
      Call GetMem('Nuc   ','ALLO','REAL',ipNuc,nComp)
*     Here we put in the k-vector
      Call FZero(CoorO,3*nComp)
      Call dcopy_(3,wavevector,1,CoorO,1)
*
*     The electromagnetic field operator contributes to all
*     irreducible irreps, hence OperI=255. Since the operator
*     it self is not symmetry adopted OperC is set to a dummy value.
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
      Call dcopy_(nComp,Zero,0,Work(ipNuc),1)
*     Write (*,*) 'Here we go!'
      Call OneEl(EMFInt,EMFMem,Label,ipList,OperI,nComp,
     &           CoorO,nOrdOp,Work(ipNuc),rHrmt,OperC,
     &           dum,1,dum,idum,0,0,
     &           dum,1,0)
*
      Call GetMem('Nuc   ','FREE','REAL',ipNuc,nComp)
      Call Deallocate_Auxiliary()
*
*     The A^2 term
*
      nOrdOp = 0
      Label='TMOS2'
      nComp = 2
      Call Allocate_Auxiliary()
      Call GetMem('Nuc   ','ALLO','REAL',ipNuc,nComp)
*     Here we put in the k-vector
      Call FZero(CoorO,3*nComp)
      Call dcopy_(3,wavevector,1,CoorO,1)
*     Change the argument to 2xA
      Call dscal_(3,2.0D0,CoorO,1)
*
*     The electromagnetic field operator contributes to all
*     irreducible irreps, hence OperI=255. Since the operator
*     it self is not symmetry adopted OperC is set to a dummy value.
*
      OperI(1   ) = 255
      OperI(1+1 ) = 255
      OperC(1   ) = 0 ! Dummy
      OperC(1+1 ) = 0 ! Dummy
*
      Call dcopy_(nComp,Zero,0,Work(ipNuc),1)
*     Write (*,*) 'Here we go!'
      Call OneEl(EMFInt,EMFMem,Label,ipList,OperI,nComp,
     &           CoorO,nOrdOp,Work(ipNuc),rHrmt,OperC,
     &           dum,1,dum,idum,0,0,
     &           dum,1,0)
*
      Call GetMem('Nuc   ','FREE','REAL',ipNuc,nComp)
      Call Deallocate_Auxiliary()
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
      Subroutine Allocate_Auxiliary()
      Implicit None
#include "stdalloc.fh"
*
      Call mma_Allocate(ipList,nComp)
      Call mma_Allocate(OperI,nComp)
      Call mma_Allocate(OperC,nComp)
      Call mma_Allocate(CoorO,3*nComp)
*
      Return
      End Subroutine Allocate_Auxiliary
      Subroutine Deallocate_Auxiliary()
      Implicit None
#include "stdalloc.fh"
*
      Call mma_Deallocate(OperC)
      Call mma_Deallocate(OperI)
      Call mma_Deallocate(ipList)
      Call mma_Deallocate(CoorO)
*
      Return
      End Subroutine Deallocate_Auxiliary
*
      End Subroutine TMOSInt
