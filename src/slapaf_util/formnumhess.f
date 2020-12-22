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
      Subroutine FormNumHess(nIter,nInter,Delta,Stop,nAtom,Cubic,iNeg,
     &                       DipM)
      use Slapaf_Info, only: qInt, Shift, dqInt, Degen, Smmtrc
      use Slapaf_Parameters, only: Curvilinear, mTROld
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "stdalloc.fh"
#include "print.fh"
      Real*8 DipM(3,nIter)
      Logical Stop, Cubic, Found
      Real*8 rDum(1)
      Real*8, Allocatable:: FEq(:), dDipM(:), KtB(:), HB(:), Hx(:),
     &                      Degen2(:), H(:), IRInt(:)
*                                                                      *
************************************************************************
*                                                                      *
      mTR = mTROld
*                                                                      *
************************************************************************
*                                                                      *
      iRout = 182
      iPrint = nPrint(iRout)
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_allocate(dDipM,3*(nInter+mTR),Label='dDipM')
      dDipM(:)=Zero
*                                                                      *
************************************************************************
*                                                                      *
*-----Form the Hessian matrix via a 2-point numerical differentiation.
*
      Stop = .False.
      Call mma_allocate(H,nInter**2,Label='H')
      If (Cubic) Then
         Call mma_allocate(FEq,nInter**3,Label='FEq')
      Else
         Call mma_allocate(FEq,1,Label='FEq')
      End If

      Call NmHess(Shift,nInter,dqInt,nIter,H,Delta,qInt,FEq,Cubic,DipM,
     &            dDipM)

      Write (6,*)
      Write (6,*) ' Numerical differentiation is finished!'
      If (iPrint.GE.98) Call RecPrt(' Numerical force constant matrix',
     &     ' ',H,nInter,nInter)
*
      Call Add_Info('Numerical Hessian',H,nInter**2,2)
      Call Put_dArray('Hss_Q',H,nInter**2)
      Call Put_dArray('Hss_upd',rDum,0)
*
*-----That is in internal coordinates, now transform it to Cartesians
*     d^2E/dx^2 = dQ/dx d^2E/dQ^2 dQ/dx + d^2Q/dx^2 dE/dQ
*
      Call Qpg_dArray('KtB',Found,nKtB)
      If (.Not.Found) Call Abend()
      nDim=nKtB/nInter
      Call mma_allocate(KtB,nDim*nInter,Label='KtB')
      Call mma_allocate(HB,nDim*nInter,Label='HB')
      Call mma_allocate(Hx,nDim**2,Label='Hx')
      Call mma_allocate(Degen2,nDim,Label='Degen2')
      Call Get_dArray('KtB',KtB,nKtB)
      Call DGeMM_('N','T',nInter,nDim,nInter,One,H,nInter,
     &                    KtB,nDim,Zero,HB,nInter)
      Call DGeMM_('T','T',nDim,nDim,nInter,One,HB,nInter,
     &                    KtB,nDim,Zero,Hx,nDim)
      i=0
      Do ii=1,nAtom
         Do ij=1,3
            If (Smmtrc(ij,ii)) Then
               i=i+1
               Degen2(i) = Degen(ij,ii)
            End If
         End Do
      End Do
*
      If (Curvilinear) Then
*        Compute and add the d^2Q/dx^2 dE/dQ part
         Call dBuu(Degen2,nInter,nDim,dqInt(1,1),Hx,.True.)
      End If
*
      Call Put_dArray('Hss_X',Hx,nDim**2)
      Call mma_deallocate(KtB)
      Call mma_deallocate(HB)
      Call mma_deallocate(Hx)
      Call mma_deallocate(Degen2)
      Call mma_deallocate(H)
*
      If (Cubic) Then
         Call RecPrt(' Numerical cubic force constant matrix',' ',
     &               FEq,nInter**2,nInter)
         Call Add_Info('Numerical anharm. cons.',FEq,nInter**3,2)
      End If
      Call mma_deallocate(FEq)
*                                                                      *
************************************************************************
*                                                                      *
*---- Do an harmonic frequency analysis
*
      Call mma_allocate(IRInt,nInter+mTR,Label='IRInt')

      Call HrmFrq(nAtom,nInter,iNeg,dDipM,mTR,DipM,IRInt)

      Call Add_Info('Numerical IR Intensities',IRInt,nInter,2)
      Call mma_deallocate(IRInt)
      Write (6,*)
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_deallocate(dDipM)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
