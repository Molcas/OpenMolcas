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
#include "compiler_features.h"
#ifdef _IN_MODULE_
      SubRoutine OneEl_Integrals(Kernel,KrnlMm,Label,ip,lOper,nComp,
     &                           CCoor,nOrdOp,rHrmt,iChO,Integrals)
      use Symmetry_Info, only: nIrrep
      Implicit Real*8 (A-H,O-Z)
      External Kernel, KrnlMm
#include "stdalloc.fh"
#include "real.fh"
      Character Label*8
      Real*8 CCoor(3,nComp)
      Integer ip(nComp), lOper(nComp), iChO(nComp), iStabO(0:7)
      Integer iTwoj(0:7)
      Data iTwoj/1,2,4,8,16,32,64,128/
      Dimension dum(1),idum(1)
      Real*8, Allocatable:: Integrals(:)
*                                                                      *
************************************************************************
*                                                                      *
*#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
      Write (6,*) ' In OneEl: Label', Label
      Write (6,*) ' In OneEl: nComp'
      Write (6,'(1X,8I5)') nComp
      Write (6,*) ' In OneEl: lOper'
      Write (6,'(1X,8I5)') lOper
      Write (6,*) ' In OneEl: n2Tri'
      Do iComp = 1, nComp
         ip(iComp) = n2Tri(lOper(iComp))
      End Do
      Write (6,'(1X,8I5)') (ip(iComp),iComp=1,nComp)
      Call RecPrt(' CCoor',' ',CCoor,3,nComp)
#endif
*                                                                      *
************************************************************************
*                                                                      *
*-----Compute the number of blocks from each component of the operator
*     and the irreps it will span.
*
      nIC = 0
      llOper = 0
      Do iComp = 1, nComp
         llOper = iOr(llOper,lOper(iComp))
         Do iIrrep = 0, nIrrep-1
            If (iAnd(lOper(iComp),iTwoj(iIrrep)).ne.0) nIC = nIC + 1
         End Do
      End Do
#ifdef _DEBUGPRINT_
      Write (6,*) ' nIC =',nIC
#endif
      If (nIC.eq.0) Then
         Call WarningMessage(2,'OneEl_Integrals: nIC.eq.0')
         Call Abend()
      End If
      Call SOS(iStabO,nStabO,llOper)
*                                                                      *
************************************************************************
*                                                                      *
*     Allocate memory for symmetry adapted one electron integrals.
*     Will just store the unique elements, i.e. low triangular blocks
*     and lower triangular elements in the diagonal blocks.
*
      ip(:)=-1
      LenTot=0
      Do iComp = 1, nComp
         ip(iComp)=1+LenTot
         LenInt=n2Tri(lOper(iComp))
         LenTot=LenTot+LenInt+4
      End Do
      call mma_allocate(Integrals,LenTot,Label='Integrals')
      Integrals(:)=Zero
*                                                                      *
************************************************************************
*                                                                      *
*---- Compute all SO integrals for all components of the operator.
*
      Call OneEl_Inner
     &           (Kernel,KrnlMm,Label,ip,lOper,nComp,CCoor,
     &            nOrdOp,rHrmt,iChO,
     &            dum,dum,1,idum,0,0,
     &            iStabO,nStabO,nIC,
     &            Dum,1,0,Integrals,LenTot)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End

#elif !defined (EMPTY_FILES)

! Some compilers do not like empty files
#include "macros.fh"
      dummy_empty_procedure(OneEl_Integrals)

#endif
