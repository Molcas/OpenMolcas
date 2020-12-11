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
      Subroutine LNM(Cart,nAtoms,Hess,Scrt1,Scrt2,Vctrs,
     &               mAtoms,nDim,iAnr,nIter,iOptH,
     &               Schlegel,Analytic_Hessian,
     &               iOptC,iTabBonds,iTabAtoms,nBonds,nMax,nHidden)
      use Symmetry_Info, only: nIrrep
      use Slapaf_Info, only: Degen, Smmtrc
      Implicit Real*8 (a-h,o-z)
#include "print.fh"
#include "real.fh"
#include "stdalloc.fh"
#include "angstr.fh"
      Real*8 Cart(3,nAtoms+nHidden), Hess(3*nAtoms*(3*nAtoms+1)/2),
     &       Scrt1((3*nAtoms)**2),
     *       Scrt2((3*nAtoms)**2), Vctrs(3*nAtoms,nDim)
      Integer   iANr(nAtoms+nHidden), iTabBonds(3,nBonds),
     &          iTabAtoms(2,0:nMax,nAtoms+nHidden)
      Logical Schlegel, Analytic_Hessian,
     &        Found, RunOld
      Real*8, Allocatable:: TanVec(:), HTanVec(:)

      iRout=120
      iPrint=nPrint(iRout)
*#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
      Call RecPrt('In LNM: Cart',' ',Cart,3,nAtoms)
      If (nHidden.ne.0) Call RecPrt('In LNM: Cart(hidden atoms)',' ',
     &                              Cart(1,nAtoms+1),3,nHidden)
      Call RecPrt('In LNM: Vctrs',' ',Vctrs,3*nAtoms,nDim)
      Write (6,*) 'iAnr=',iANr
      Write (6,*) 'Analytic Hessian=',Analytic_Hessian
#endif
*                                                                      *
************************************************************************
*                                                                      *
*     CARTESIAN HESSIAN is generated approximately here or the analytic*
*     is read from file.                                               *
*                                                                      *
*-----Retrieve the analytic Hessian or compute the Hessian in cartesians
*     from the Hessian model function.
*
*                                                                      *
************************************************************************
*                                                                      *
      If (Analytic_Hessian) Then
*                                                                      *
************************************************************************
*                                                                      *
*
*        The Hessian matrix is read in the basis of the symmetric
*        displacements.
*
         Len3=ndim
         Len3=Len3*(Len3+1)/2

         Call qpg_dArray('Analytic Hessian',Found,Len)
         If (Found) Then
            Call Get_AnalHess(Hess,Len3)
            RunOld=.False.
         Else
            Call NameRun('RUNOLD')
            Call qpg_dArray('Analytic Hessian',Found,Len)
            If (.NOT.Found) Then
               Call WarningMessage(2,
     &                      ' Error in LNM: Analytic Hessian not found')
               Call Abend()
            End If
            Call Get_AnalHess(Hess,Len3)
            Call NameRun('#Pop')
            RunOld=.True.
         End If
#ifdef _DEBUGPRINT_
         Call TriPrt('LNM: Analytic Hessian',' ',Hess,nDim)
         Call RecPrt('LNM: Degen',' ',Degen,SIZE(Degen,1),SIZE(Degen,2))
#endif
*
*        Modify matrix with degeneracy factors and square the
*        matrix.
*
         Call FZero(Scrt1,nDim**2)
         ii = 0
         Do i = 1, 3*mAtoms
            iAtom=(i+2)/3
            ixyz = i - (iAtom-1)*3
            If (Smmtrc(ixyz,iAtom)) Then
               ii = ii + 1
               jj = 0
               Do j = 1, i
                  jAtom=(j+2)/3
                  jxyz = j - (jAtom-1)*3
                  If (Smmtrc(jxyz,jAtom)) Then
                     jj = jj + 1
                     ijTri=ii*(ii-1)/2 + jj
                     ij = (jj-1)*ndim + ii
                     ji = (ii-1)*ndim + jj
                     Tmp= Hess(ijTri)*sqrt(degen(ixyz,iAtom)
     &                                    *degen(jxyz,jAtom))
                     Scrt1(ij) = Tmp
                     Scrt1(ji) = Tmp
                  End If
               End Do
            End If
         End Do
#ifdef _DEBUGPRINT_
         Call TriPrt('Hessian(anal.)',' ',Hess,nDim)
         Write (6,*) 'nDim,nAtoms,mAtoms=',
     &                nDim,nAtoms,mAtoms
         Call RecPrt('Hessian(anal.)',' ',Scrt1,nDim,nDim)
#endif
*
**       If the analytic Hessian corresponds to the current iteration,
**       disable Hessian updating
*
         If (RunOld) Call NameRun('RUNOLD')
         Call Get_iScalar('HessIter',IterHess)
         If (RunOld) Call NameRun('#Pop')
         If (IterHess.eq.nIter) iOptH = iOr(8,iAnd(iOptH,32))
*                                                                      *
************************************************************************
*                                                                      *
      Else   ! Use the Hessian Model Function
*                                                                      *
************************************************************************
*                                                                      *
         Call ddV(Cart,nAtoms,Hess,iANr,Schlegel,iOptC,
     &            iTabBonds,iTabAtoms,nBonds,nMax,nHidden)
#ifdef _DEBUGPRINT_
         If (iPrint.ge.19)
     &      Call TriPrt(' The Model Hessian','(12f9.5)',Hess,3*nAtoms)
#endif
*                                                                      *
************************************************************************
*                                                                      *
*----    Square the matrix
*
         Do i = 1, 3*nAtoms
            Do j = 1, i
               ijTri=i*(i-1)/2 + j
               ij = (j-1)*(3*nAtoms) + i
               ji = (i-1)*(3*nAtoms) + j
               Scrt1(ij) = Hess(ijTri)
               Scrt1(ji) = Hess(ijTri)
            End Do
         End Do
*        Call RecPrt('Scrt1',' ',Scrt1,3*nAtoms,3*nAtoms)
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUGPRINT_
         If (iPrint.ge.19)
     &      Call RecPrt(' Scrt1',' ',Scrt1,3*nAtoms,3*nAtoms)
#endif
         If (nIrrep.eq.1) Go To 99
*
*------  Now project out the total symmetric part of the Hessian
*
         Call DGEMM_('N','N',
     &               3*nAtoms,nDim,3*nAtoms,
     &               1.0d0,Scrt1,3*nAtoms,
     &               Vctrs,3*nAtoms,
     &               0.0d0,Scrt2,3*nAtoms)
#ifdef _DEBUGPRINT_
         If (iPrint.ge.19)
     &      Call RecPrt(' Scrt2',' ',Scrt2,3*nAtoms,nDim)
#endif
         Call DGEMM_('T','N',
     &               nDim,nDim,3*nAtoms,
     &               1.0d0,Vctrs,3*nAtoms,
     &               Scrt2,3*nAtoms,
     &               0.0d0,Scrt1,nDim)
#ifdef _DEBUGPRINT_
         If (iPrint.ge.19)
     &      Call RecPrt(' The Symmetrized Hessian',' ',Scrt1,nDim,nDim)
#endif
 99      Continue
*                                                                      *
************************************************************************
*                                                                      *
*        Now project hessian using reaction path vector if TS search
*
         Call qpg_darray('TanVec',Found,nRP)
         If (Found) Then
            If (nRP.ne.3*mAtoms) Then
               Call WarningMessage(2,' Error in LNM: nRP.ne.3*mAtoms')
               Write (6,*) 'nRP,3*mAtoms=',nRP,mAtoms
               Call Abend()
            End If
*
            Call mma_allocate(TanVec,nRP,Label='TanVec')
            Call mma_allocate(HTanVec,nRP**2,Label='HTanVec')
            TanVec(:)=Zero
            HTanVec(:)=Zero

            Call Get_dArray('TanVec',TanVec,nRP)
*
*           Call RecPrt('TanVec',' ',TanVec,nRP,1)
            i = 0
            Do ix = 1, nRP
               iAtom = (ix+2)/3
               ixyz = ix - (iAtom-1)*3
               If (Smmtrc(ixyz,iAtom)) Then
                  i = i + 1
                  TanVec(i)=TanVec(ix)
               End If
            End Do
            nRP = nDim
*           Call RecPrt('TanVec',' ',TanVec,nRP,1)
*
*           Compute H|i>
*
            Call dGeMV_('N',nRP,nRP,
     &                 1.0d0,scrt1,nRP,
     &                       TanVec,1,
     &                 0.0d0,HTanVec,1)
*
*           Compute <i|H|i>
*
            eigen=0.0d0
            Do i=1,nRP
               eigen=eigen+HTanVec(i)*TanVec(i)
            End Do
*           Write (*,*) 'Eigen=',eigen
            If (eigen.lt.Zero) Go To 98
*
*           Form H - H|i><i| - |i><i|H
*
            Do i=0,nRP-1
               Do j=0,nRP-1
                  scrt1(nRP*i+j+1) = scrt1(nRP*i+j+1)
     &                             - HTanVec(1+i)*TanVec(1+j)
     &                             - TanVec(1+i)*HTanVec(1+j)
               End Do
            End Do
*
*
 98         Continue
            Call mma_deallocate(HTanVec)
            Call mma_deallocate(TanVec)
         End If
*        Call RecPrt('Scrt1(final)',' ',Scrt1,nDim,nDim)
*                                                                      *
************************************************************************
*                                                                      *
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Do j = 1, nDim
         Do i = 1, j
            ij = (j-1)*nDim + i
            ijTri= i*(i-1)/2 + j
            Hess(ijTri)=Scrt1(ij)
         End Do
      End Do
*     Call RecPrt('Scrt1(final)',' ',Scrt1,nDim,nDim)
*                                                                      *
************************************************************************
*                                                                      *
      End
