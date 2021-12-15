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
      SubRoutine StPert
      Implicit Real*8(a-h,o-z)

#include "Input.fh"
#include "Pointers.fh"
#include "disp_mclr.fh"
#include "WrkSpc.fh"
#include "spin.fh"
#include "cstate_mclr.fh"
#include "Files_mclr.fh"
      Character*16 Label
      Character*8  MckLbl
      Character*288 Header
*
      nLen=0
      Do iS=1,nSym
         nLen=nLen+lDisp(is)*(lDisp(is)+1)/2
      End Do
      Call GetMem('CONN','Allo','Real',ipHss,nLen)
      call dcopy_(nLen,[0.0d0],0,Work(ipHss),1)
*
      If (.Not.Mckinley) Then
         irc=-1
         iopt=1
         Call OPNMCK(irc,iopt,FNMCK,LUMCK)
         If (irc.ne.0) Then
            Write (6,*) 'StPert: Error opening MCKINT'
            Call Abend()
         End If
         irc=-1
         iopt=0
         LABEL='SEWARD'
         If (PT2) LABEL='PT2LAG'
         MckLbl='PERT    '
         Call cWrMck(iRC,iOpt,MckLbl,1,LABEL,iDummer)
         If (irc.ne.0) Then
             Write (6,*) 'StPert: Error writing to MCKINT'
             Write (6,'(A,A)') 'MckLbl=',MckLbl
            Call Abend()
         End If
         irc=-1
         iopt=0
         MckLbl='NDISP   '
         Call WrMck(iRC,iOpt,MckLbl,1,[ndisp],iDummer)
         If (irc.ne.0) Then
             Write (6,*) 'StPert: Error writing to MCKINT'
             Write (6,'(A,A)') 'MckLbl=',MckLbl
            Call Abend()
         End If
         irc=-1
         iopt=0
         MckLbl='TDISP   '
         Call WrMck(iRC,iOpt,MckLbl,1,ntpert,iDummer)
         If (irc.ne.0) Then
             Write (6,*) 'StPert: Error writing to MCKINT'
             Write (6,'(A,A)') 'MckLbl=',MckLbl
            Call Abend()
         End If
         irc=-1
         iopt=0
         MckLbl='Title'
         Call cWrMck(iRC,iOpt,MckLbl,1,Header,iDummer)
         If (irc.ne.0) Then
             Write (6,*) 'StPert: Error writing to MCKINT'
             Write (6,'(A,A)') 'MckLbl=',MckLbl
            Call Abend()
         End If
         irc=-1
         iopt=0
         MckLbl='nSym'
         Call WrMck(iRC,iOpt,MckLbl,1,[nSym],iDummer)
         If (irc.ne.0) Then
             Write (6,*) 'StPert: Error writing to MCKINT'
             Write (6,'(A,A)') 'MckLbl=',MckLbl
            Call Abend()
         End If
         irc=-1
         iopt=0
         MckLbl='nBas'
         Call WrMck(iRC,iOpt,MckLbl,1,nBas,iDummer)
         If (irc.ne.0) Then
             Write (6,*) 'StPert: Error writing to MCKINT'
             Write (6,'(A,A)') 'MckLbl=',MckLbl
            Call Abend()
         End If
         irc=-1
         iopt=0
         MckLbl='ldisp'
         Call WrMck(iRC,iOpt,MckLbl,1,ldisp,iDummer)
         If (irc.ne.0) Then
             Write (6,*) 'StPert: Error writing to MCKINT'
             Write (6,'(A,A)') 'MckLbl=',MckLbl
            Call Abend()
         End If
         irc=-1
         iopt=0
         MckLbl='chdisp'
         Call cWrMck(iRC,iOpt,MckLbl,1,swlbl(1),iDummer)
         If (irc.ne.0) Then
             Write (6,*) 'StPert: Error writing to MCKINT'
             Write (6,'(A,A)') 'MckLbl=',MckLbl
            Call Abend()
         End If
         irc=-1
         iopt=0
         MckLbl='NISH'
         Call WrMck(iRC,iOpt,MckLbl,1,nish,iDummer)
         If (irc.ne.0) Then
             Write (6,*) 'StPert: Error writing to MCKINT'
             Write (6,'(A,A)') 'MckLbl=',MckLbl
            Call Abend()
         End If
         irc=-1
         iopt=0
         MckLbl='NASH'
         Call WrMck(iRC,iOpt,MckLbl,1,nash,iDummer)
         If (irc.ne.0) Then
             Write (6,*) 'StPert: Error writing to MCKINT'
             Write (6,'(A,A)') 'MckLbl=',MckLbl
            Call Abend()
         End If
      End If
*
      If (SPINPOL) Then
         call coeff(ralphas,rbetaa,rbetas)
         rms=DBLE(ms2)/2.0d0
         nAct  = 0
         Do iSym = 1, nSym
            nAct = nAct + nAsh(iSym)
         End Do
         nG=nAct**2
         nG2=nAct**4
         Call GETMEM('FAMO_SPINp','ALLO','REAL',ipfamo_spinp,ndens2)
         Call GETMEM('FAMO_SPINm','ALLO','REAL',ipfamo_spinm,ndens2)
         Call GetMem('ipg2tmm','ALLO','REAL',ipg2mp,nG2)
         Call GetMem('ipg2tmm','ALLO','REAL',ipg2pp,nG2)
         Call GetMem('ipg2tmm','ALLO','REAL',ipg2mm,nG2)
         Call GetMem('ipg2spin','ALLO','REAL',ipfm,nG2)
         Call GetMem('ipg2spin','ALLO','REAL',ipfp,nG2)
         Call GetMem('ipg1m','ALLO','REAL',ipg1p,nG)
         Call GetMem('ipg1p','ALLO','REAL',ipg1m,nG)
         itype=2
         Call SpinDens(Work(ipin(ipCI)),Work(ipin(ipCI)),
     &                 STATE_SYM,
     &                 STATE_SYM,Work(ipg2mm),
     &                 Work(ipg2mp),work(ipg2pp),
     &                 Work(ipfm),Work(ipfp),
     %                 Work(ipg1m),work(ipg1p),
     &                 itype)
         Call GETMEM('Tmp2','ALLO','REAL',ipTmp2,ndens2)
         Call GETMEM('Tmp1','MAX','REAL',iptmp1,nMax)
         Call GETMEM('Tmp1','ALLO','REAL',iptmp1,nMax/2)
         Call Ex_spin(Work(ipg1p),Work(ipFAMO_Spinp),Work(ipTmp1),
     &                nMax/2,Work(ipTmp2))
         Call Ex_spin(Work(ipg1m),Work(ipFAMO_Spinm),Work(ipTmp1),
     &              nMax/2,Work(ipTmp2))
         Call GetMem('Tmp1','FREE','REAL',iptmp1,nmax/2)
         Call GetMem('tmp2','FREE','REAL',ipTmp2,ndens2)
      End If
*
      Return
      End
