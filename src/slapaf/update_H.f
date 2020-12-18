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
      Subroutine Update_H(nWndw,H,nInter,nIter,
     &                    iOptC,dq,g,iOptH,
     &                    jPrint,GNrm,
     &                    nsAtom,Store,
     &                    AllowFindTS)
      use Slapaf_Info, only: MF
      use Slapaf_Parameters, only: IRC, Mode
      Implicit Real*8 (a-h,o-z)
#include "stdalloc.fh"
#include "real.fh"
      Real*8 H(nInter,nInter), dq(nInter,nIter), g(nInter,nIter)
      Logical Old_MF, Store, AllowFindTS
      Real*8, Allocatable:: Tmp(:,:)
*                                                                      *
************************************************************************
*                                                                      *
*---- Update the Hessian. The procedure here loops over gradients
*     of the nWndw last iteration. The reason for this is that the
*     initial Hessian is reset to new values at each iteration.
*     The anharmonic constants used here is the most recently
*     updated version.
*
      Call DrvUpH(nWndw,nIter,H,nInter,dq,g,iOptH,jPrint,IterHess)

      Call Chk4NAN(nInter*nInter,H,ierr)

      If (ierr.ne.0) Call SysAbendMsg('Update_H','NaNs in Hessian','')
      If (Store) Call Put_dArray('Hss_upd',H,nInter**2)
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*     Check if we have an old reaction mode
*
      Old_MF=DDot_(3*nsAtom,MF,1,MF,1).ne.Zero
      Old_MF = Old_MF .and. Mode.ne.0 .and. IRC.eq.0
      Call mma_allocate(Tmp,3,nsAtom,Label='Tmp')
      If (Old_MF) Then
         If (jPrint.ge.6)
     &      Write (6,*) ' Reading old reaction mode from disk'
         Tmp(:,:) = MF(:,:)
         Mode=1  ! any number between 1 and nInter
         iOptC=iOr(iOptC,8192)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Massage the new Hessian
*
      Call FixHess(H,nInter,iOptC,Tmp,GNrm,nsAtom,(IterHess.eq.nIter),
     &             AllowFindTS)
*                                                                      *
************************************************************************
*                                                                      *
*     Store the new reaction mode if any!
*
      If (Mode.gt.0 .and. Mode.le.nInter) Then
         If (jPrint.ge.6) Write (6,*)
     %      ' Storing new reaction mode on disk'
         MF(:,:)=Tmp(:,:)
      End If
      Call mma_deallocate(Tmp)
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      If (jPrint.ge.99) Then
         Call RecPrt('Update_H: Updated Hessian',' ',H,nInter,nInter)
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
