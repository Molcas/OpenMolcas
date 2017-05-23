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
      Subroutine DelGHOST_MBPT(ipCMO,ipCMO_t,lthCMO,
     &                         ipEOrb,ipEOrb_t,lthEOr)
************************************************************************
************************************************************************
      Implicit Real*8 (A-H,O-Z)

*     declaration of calling arguments
      Integer ipCMO,ipCMO_t,ipEOrb,ipEOrb_t,lthCMO,lthEOr

#include "real.fh"
#include "Molcas.fh"
#include "mxdim.fh"
#include "namact.fh"
#include "corbinf.fh"
#include "orbinf2.fh"
#include "mbpt2aux.fh"
#include "WrkSpc.fh"
#include "files_mbpt2.fh"
#include "print_mbpt2.fh"

*     declaration of local variables...
      Integer nZero(8)
      Character*(LENIN4) Name(mxBas)
      Logical Debug
      Data Debug/.False./


#include "SysDef.fh"

      If (.not.DelGHOST) Return

      If (Debug) Then
         Write (6,'(A,8I5)') 'nSym:',nSym
         Write (6,'(A,8I5)') 'nBas:',(nBas(i),i=1,nSym)
         Write (6,'(A,8I5)') 'nOrb:',(nOrb(i),i=1,nSym)
         Write (6,'(A,8I5)') 'nOcc:',(nOcc(i),i=1,nSym)
         Write (6,'(A,8I5)') 'nFro:',(nFro(i),i=1,nSym)
      End If
      Do iSym=1,nSym
         nDel(iSym)=nBas(iSym)-nOrb(iSym)
         nExt(iSym)=nOrb(iSym)-nOcc(iSym)-nFro(iSym)
         nDsto(iSym)=nDel(iSym)
         nZero(iSym)=0
      End Do
*
      Call GetMem('CMO   ','Allo','Real',ipCMO,lthCMO)
*
      Call GetMem('EOrb  ','Allo','Real',ipEOrb,lthEOr)
*
      Write(6,'(A)')
     &    '-------------------------------------------------------'
      Write(6,'(A)') ' GHOST virtual space removal'
      Write(6,'(A)')
     &    '-------------------------------------------------------'
      Write(6,'(A,8I4)')
      Write(6,'(A,8I4)')
     &  ' Secondary orbitals before selection:',(nExt(i),i=1,nSym)
      Write(6,'(A,8I4)')
     &  ' Deleted orbitals before selection:  ',(nDel(i),i=1,nSym)

      Call Get_iScalar('Unique atoms',nUniqAt)
      Call Get_cArray('Unique Basis Names',Name,(LENIN4)*nnB)
      Call Delete_GHOSTS(irc,nSym,nBas,nFro,nOcc,nZero,nExt,nDel,
     &       NAME,nUniqAt,thr_ghs,.False.,Work(ipCMO_t),Work(ipEOrb_t))

      If (irc.ne.0) Then
         write(6,*) 'Delete_GHOSTS returned rc= ',irc
         Call abend()
      EndIf
      Write(6,'(A,8I4)')
      Write(6,'(A)')
     &    '-------------------------------------------------------'
      Write(6,'(A,8I4)')
      Write(6,'(A,8I4)')
*
*...  set MO coefficients of the deleted orbitals to zero
*     Observe that these are not included at all in the basis
      iStart   = ipCMO
      iStart_t = ipCMO_t
      Do iSym=1,nSym
         call dcopy_(nOrb(iSym)*nBas(iSym),Work(iStart_t),1,
     &                                    Work(iStart),1)
         iStart   = iStart   + nOrb(iSym)*nBas(iSym)
         iStart_t = iStart_t + nOrb(iSym)*nBas(iSym)
         call dcopy_((nBas(iSym)-nOrb(iSym))*nBas(iSym),
     &               Zero,0,Work(iStart),1)
         iStart   = iStart   +
     &              (nBas(iSym)-nOrb(iSym))*nBas(iSym)
      End Do
      Call GetMem('CMO_t ','Free','Real',ipCMO_t,lthCMO)
*
*...  set energies of the deleted orbitals to zero
*
      iStart   = ipEOrb
      iStart_t = ipEOrb_t
      Do iSym=1,nSym
         call dcopy_(nOrb(iSym),Work(iStart_t),1,Work(iStart),1)
         iStart  = iStart   + nOrb(iSym)
         iStart_t= iStart_t + nOrb(iSym)
*
         call dcopy_(nBas(iSym)-nOrb(iSym),Zero,0,Work(iStart),1)
         iStart  = iStart   + nBas(iSym)-nOrb(iSym)
      End Do
      Call GetMem('EOrb_t','Free','Real',ipEOrb_t,lthEOr)
*
      Return
      End
