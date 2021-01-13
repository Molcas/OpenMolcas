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
      SubRoutine StPert()
      Use Arrays, only: Hss, FAMO_SpinP, FAMO_SpinM,
     &                  G2mm, G2mp, G2pp, Fp, Fm, G1p, G1m
      use ipPage, only: W
      Implicit Real*8(a-h,o-z)

#include "real.fh"
#include "Input.fh"
#include "Pointers.fh"
#include "disp_mclr.fh"
#include "stdalloc.fh"
#include "spin.fh"
#include "cstate_mclr.fh"
#include "Files_mclr.fh"
      Character(LEN=16) Label
      Character(LEN=8)  MckLbl
      Character(LEN=288) Header
      Real*8, Allocatable:: Tmp1(:), Tmp2(:)
*
      nHss=0
      Do iS=1,nSym
         nHss=nHss+lDisp(is)*(lDisp(is)+1)/2
      End Do
      Call mma_allocate(Hss,nHss,Label='Hss')
      Hss(:)=Zero
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
         Call mma_allocate(famo_spinp,ndens2,Label='famo_spinp')
         Call mma_allocate(famo_spinm,ndens2,Label='famo_spinm')
         Call mma_allocate(G2mp,nG2,Label='G2mp')
         Call mma_allocate(G2pp,nG2,Label='G2pp')
         Call mma_allocate(G2mm,nG2,Label='G2mm')
         Call mma_allocate(Fm,nG2,Label='Fm')
         Call mma_allocate(Fp,nG2,Label='Fp')
         Call mma_allocate(G1p,nG,Label='G1p')
         Call mma_allocate(G1m,nG,Label='G1m')
         itype=2
         irc=ipin(ipCI)
         Call SpinDens(W(ipCI)%Vec,W(ipCI)%Vec,
     &                 STATE_SYM,STATE_SYM,G2mm,G2mp,G2pp,
     &                 Fm,Fp,G1m,G1p,iType)

         Call mma_allocate(Tmp2,ndens2,Label='Tmp2')
         Call mma_MaxDBLE(nMax)
         Call mma_allocate(Tmp1,nMax/2,Label='Tmp1')

         Call Ex_spin(G1p,FAMO_Spinp,Tmp1,nMax/2,Tmp2)
         Call Ex_spin(G1m,FAMO_Spinm,Tmp1,nMax/2,Tmp2)

         Call mma_deallocate(Tmp1)
         Call mma_deallocate(Tmp2)
      End If
*
      Return
      End
