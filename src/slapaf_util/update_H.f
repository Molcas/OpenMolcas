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
     &                    iOptC,Mode,ipMF,dq,g,iNeg,iOptH,
     &                    HUpMet,nRowH,jPrint,GNrm,
     &                    GNrm_Threshold,nAtoms,IRC,Store)
      Implicit Real*8 (a-h,o-z)
#include "WrkSpc.fh"
#include "real.fh"
      Real*8 H(nInter,nInter), dq(nInter,nIter),
     &       g(nInter,nIter)
      Integer iNeg(2)
      Logical Old_MF, Store
      Character*6 HUpMet
*                                                                      *
************************************************************************
*                                                                      *
*---- Update the Hessian. The procedure here loops over gradients
*     of the nWndw last iteration. The reason for this is that the
*     initial Hessian is reset to new values at each iteration.
*     The anharmonic constants used here is the most recently
*     updated version.
*
      Write (6,*) 'LU32',nWndw,nIter,H,nInter,dq,g,iOptH,HUpMet,nRowH,
     &            jPrint,IterHess
      Call DrvUpH(nWndw,nIter,H,nInter,dq,g,iOptH,HUpMet,nRowH,
     &            jPrint,IterHess)
      Write (6,*) 'LU35'
      Call Chk4NAN(nInter*nInter,H,ierr)
      Write (6,*) 'LU37',ierr,H
      If (ierr.ne.0) Call SysAbendMsg('Update_H','NaNs in Hessian','')
      If (Store) Call Put_dArray('Hss_upd',H,nInter**2)
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*     Check if we have an old reaction mode
*
      Old_MF=DDot_(3*nAtoms,Work(ipMF),1,Work(ipMF),1).ne.Zero
      Old_MF = Old_MF .and. Mode.ne.0 .and. IRC.eq.0
      Write (6,*) 'LU48'
      Call Allocate_Work(ipTmp,3*nAtoms)
      If (Old_MF) Then
         If (jPrint.ge.6)
     &      Write (6,*) ' Reading old reaction mode from disk'
         Write (6,*) 'LU50'
         call dcopy_(3*nAtoms,Work(ipMF),1,Work(ipTmp),1)
         Mode=1  ! any number between 1 and nInter
         iOptC=iOr(iOptC,8192)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Massage the new Hessian
*
      Call FixHess(H,nInter,iOptC,Mode,nIter,ipTmp,GNrm,
     &             GNrm_Threshold,iNeg,nAtoms,(IterHess.eq.nIter))
*                                                                      *
************************************************************************
*                                                                      *
*     Store the new reaction mode if any!
*
      If (Mode.gt.0 .and. Mode.le.nInter) Then
         If (jPrint.ge.6) Write (6,*)
     %      ' Storing new reaction mode on disk'
         call dcopy_(3*nAtoms,Work(ipTmp),1,Work(ipMF),1)
      End If
      Call Free_Work(ipTmp)
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
