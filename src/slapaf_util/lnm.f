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
     &               mAtoms,nDim,iAnr,Smmtrc,Cx,Gx,nIter,iOptH,
     &               Degen,Schlegel,Analytic_Hessian,
     &               iOptC,iTabBonds,iTabAtoms,nBonds,nMax,nHidden,
     &               nMDstep,MMkept)
      use Symmetry_Info, only: nIrrep
      Implicit Real*8 (a-h,o-z)
#include "print.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "angstr.fh"
      Real*8 Cart(3,nAtoms+nHidden), Hess(3*nAtoms*(3*nAtoms+1)/2),
     &       Cx(3*mAtoms,nIter),
     &       Gx(3*mAtoms,nIter), Scrt1((3*nAtoms)**2),
     *       Scrt2((3*nAtoms)**2), Vctrs(3*nAtoms,nDim),
     &       Degen(3*mAtoms)
      Integer   iANr(nAtoms+nHidden), iTabBonds(3,nBonds),
     &          iTabAtoms(2,0:nMax,nAtoms+nHidden)
      Logical Smmtrc(3*mAtoms), Schlegel, Analytic_Hessian,
     &        Found, RunOld
      Character*180 Line,Get_Ln
      External Get_Ln
      Dimension XYZ(3)
*                                                                      *
************************************************************************
*                                                                      *
*define _DEBUG_
*                                                                      *
************************************************************************
*                                                                      *
      iRout=120
      iPrint=nPrint(iRout)
#ifdef _DEBUG_
      Call QEnter('LNM')
      iPrint=99
      If (iPrint.ge.19) Then
         Call RecPrt('In LNM: Cart',' ',Cart,3,nAtoms)
         If (nHidden.ne.0) Call RecPrt('In LNM: Cart(hidden atoms)',' ',
     &                                 Cart(1,nAtoms+1),3,nHidden)
         Call RecPrt('In LNM: Vctrs',' ',Vctrs,3*nAtoms,nDim)
         Write (6,*) 'iAnr=',iANr
         Write (6,*) 'Analytic Hessian=',Analytic_Hessian
      End If
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
         Call Get_AnalHess(ipHess,Len)
         RunOld=.False.
         If (Len.eq.0) Then
            Call NameRun('RUNOLD')
            Call Get_AnalHess(ipHess,Len)
            Call NameRun('#Pop')
            RunOld=.True.
         End If
         Len3=ndim
         Len3=Len3*(Len3+1)/2
         If (Len.ne.Len3) Then
            Call WarningMessage(2,' Error on RELAX in LNM: Len.ne.Len3')
            write(6,*) len,len3
            Call Abend()
         End If
         call dcopy_(Len,Work(ipHess),1,Hess,1)
         Call Free_Work(ipHess)
*
*        Modify matrix with degeneracy factors and square the
*        matrix.
*
         Call FZero(Scrt1,nDim**2)
         ii = 0
         Do i = 1, 3*mAtoms
            If (Smmtrc(i)) Then
               ii = ii + 1
               jj = 0
               Do j = 1, i
                  If (Smmtrc(j)) Then
                     jj = jj + 1
                     ijTri=ii*(ii-1)/2 + jj
                     ij = (jj-1)*ndim + ii
                     ji = (ii-1)*ndim + jj
                     Tmp= Hess(ijTri)*sqrt(degen(i)*degen(j))
                     Scrt1(ij) = Tmp
                     Scrt1(ji) = Tmp
                  End If
               End Do
            End If
         End Do
#ifdef _DEBUG_
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
cnf
c   A mean Hessian using some structures extracted form the Tinker MD?
cnf
         If (nMDstep .gt. 0) Then
            NTT = 3*nAtoms*(3*nAtoms+1)/2
            Fact = One/Dble(nMDstep)
            Call Allocate_Work(ipHtmp,NTT)
            Call dCopy_(NTT,[Zero],0,Work(ipHtmp),1)
            iMDsnap = IsFreeUnit(1)
            Call Molcas_Open(iMDsnap,'snap.hess')
            Line = Get_Ln(iMDsnap)
            If (Line(1:8) .ne. 'NMM') Then
               Write(6,'(A)') 'LNM: Cannot find NMM in snap.hess'
               Call Quit_OnUserError()
            Else
               Call Get_I1(2,nMMtot)
            End If
            Do iSnap = 1, nMDstep
               iHidden = 0
               Do iMM = 0,nMMtot-1
                  Line = Get_Ln(iMDsnap)
                  If (iWork(MMkept+iMM) .eq. 1) Then
                     iHidden = iHidden + 1
                     Call Get_F(3,XYZ,3)
                     Cart(1,nAtoms+iHidden) = XYZ(1)/Angstr
                     Cart(2,nAtoms+iHidden) = XYZ(2)/Angstr
                     Cart(3,nAtoms+iHidden) = XYZ(3)/Angstr
                  End If
               End Do
               If (iHidden .ne. nHidden) Then
                  Write(6,'(A,2i3)') 'LNM: iHidden ne nHidden',iHidden,
     &            nHidden
                  Call Quit_OnUserError()
               End If
               Call ddV(Cart,nAtoms,Hess,iANr,Schlegel,iOptC,
     &            iTabBonds,iTabAtoms,nBonds,nMax,nHidden)
               Call daxpy_(NTT,Fact,Hess,1,Work(ipHtmp),1)
            End Do
            Close(iMDsnap)
            Call Free_iWork(MMkept)
            Call dCopy_(NTT,Work(ipHtmp),1,Hess,1)
            Call Free_Work(ipHtmp)
         Else
            Call ddV(Cart,nAtoms,Hess,iANr,Schlegel,iOptC,
     &            iTabBonds,iTabAtoms,nBonds,nMax,nHidden)
         End If
#ifdef _DEBUG_
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
#ifdef _DEBUG_
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
#ifdef _DEBUG_
         If (iPrint.ge.19)
     &      Call RecPrt(' Scrt2',' ',Scrt2,3*nAtoms,nDim)
#endif
         Call DGEMM_('T','N',
     &               nDim,nDim,3*nAtoms,
     &               1.0d0,Vctrs,3*nAtoms,
     &               Scrt2,3*nAtoms,
     &               0.0d0,Scrt1,nDim)
#ifdef _DEBUG_
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
            Call GetMem('TanVec','ALLO','REAL',ipTan,nRP+nRP**2)
            Call FZero(Work(ipTan),nRP+nRP**2)
            ipHTan = ipTan + nRP
            Call Get_dArray('TanVec',Work(ipTan),nRP)
*
*           Call RecPrt('TanVec',' ',Work(ipTan),nRP,1)
            i = 0
            Do ix = 1, nRP
               If (Smmtrc(ix)) Then
                  i = i + 1
                  Work(ipTan+i-1)=Work(ipTan+ix-1)
               End If
            End Do
            nRP = nDim
*           Call RecPrt('TanVec',' ',Work(ipTan),nRP,1)
*
*           Compute H|i>
*
            Call dGeMV_('N',nRP,nRP,
     &                 1.0d0,scrt1,nRP,
     &                       Work(ipTan),1,
     &                 0.0d0,Work(ipHTan),1)
*
*           Compute <i|H|i>
*
            eigen=0.0d0
            Do i=0,nRP-1
               eigen=eigen+Work(ipHTan+i)*Work(ipTan+i)
            End Do
*           Write (*,*) 'Eigen=',eigen
            If (eigen.lt.Zero) Go To 98
*
*           Form H - H|i><i| - |i><i|H
*
            Do i=0,nRP-1
               Do j=0,nRP-1
                  scrt1(nRP*i+j+1) = scrt1(nRP*i+j+1)
     &                             - Work(ipHTan+i)*Work(ipTan+j)
     &                             - Work(ipTan+i)*Work(ipHTan+j)
               End Do
            End Do
*
*
 98         Continue
            Call GetMem('TanVec','FREE','REAL',ipTan,nRP+nRP**2)
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
#ifdef _DEBUG_
      Call QExit('LNM')
#endif
*                                                                      *
************************************************************************
*                                                                      *
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(Cx)
         Call Unused_real_array(Gx)
      End If
      End
