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
      SubRoutine CIDens_TD(iCI,iS,rP,rD)
      use ipPage, only: W
      Implicit Real*8(a-h,o-z)
#include "detdim.fh"
#include "cicisp_mclr.fh"
#include "WrkSpc.fh"
#include "crun_mclr.fh"
#include "Input.fh"
#include "Pointers.fh"
#include "spinfo_mclr.fh"
#include "cands.fh"
      Real*8 rP(*),rD(*)
      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
* LS = CI
*
*     Ok we once more want to hide Jeppes routines from
*     the eyes of the world, so everyone belives that I have done
*     all the work.
*     If we have spin dependent Hamiltonian we will work with
*     SD in all parts of the program, no CSF is necessary,
*     otherwise we will do the optimazation in CSF's to increase
*     convergence.
*
*     Input:
*
*     response: true if the density should be used for response
*               calculations
*     LS : CI Coeff for left state
*     RS : CI Coeff for right state
*     iL : Symmetry of left state
*     iR : Symmetry of right state
*
*
*               +       +
*     iS=1 E  =a  a  + a  a     ! Singlett operator
*           pq  ap aq   Bp Bq
*
*                +       +
*     iS=-1 T  =a  a  - a a     ! Triplett operator
*            pq  ap aq   Bp Bq
*
*     Output:
*
*      rP : Two Electron Density
*      rD : One Electron Density
*
*
*      irc=ipin(iCI)
*      write(*,*)'iCI*iCI', ddot_(2*nConf1,W(iCI)%Vec,1,W(iCI)%Vec,1)

      If (nconf1.eq.0) return

      Call GetMem('1Dens2','ALLO','Real',ipDe,2*n1dens)
      Call GetMem('2Dens2','ALLO','Real',ipP,3*n2dens)
      call dcopy_(n1dens,[0.0d0],0,rD,1)
      call dcopy_(n2dens,[0.0d0],0,rP,1)
      If (nocsf.eq.0) Then
        nConfL=Max(ncsf(iS),nint(xispsm(iS,1)))
        nConfR=Max(ncsf(State_SYM),nint(xispsm(STATE_SYM,1)))
        nC=Max(nconfL,nconfR)
        Call GetMem('CIL','ALLO','REAL',ipL,nC)
        Call GetMem('CIR','ALLO','REAL',ipR,nC)
c
c iCI is as long as ipcid for the timedep case!
c
        irc=ipin(iCI)
        irc=ipin(ipCI)
        Call CSF2SD(W(iCI)%Vec,Work(ipR),iS)
        Call CSF2SD(W(ipCI)%Vec,Work(ipL),State_SYM)
*
*        write(*,*)'ipL*ipL',
*     &   ddot_(nConfL,Work(ipL),1,Work(ipL),1)
*        write(*,*)'ipR*ipR',
*     &   ddot_(nConfR,Work(ipR),1,Work(ipR),1)
*        Call RecPrt('ipL',' ',Work(ipL),nConfL,1)
*        Call RecPrt('ipR',' ',Work(ipR),nConfR,1)
*
        irc=ipnout(-1)
        icsm=iS
        issm=STATE_SYM
*
*       <P|E_pq|0> & <P|e_pqrs|0> -> ipDe & ipP
*       ipL is the bra side vector
        call dcopy_(n1dens,[0.0d0],0,Work(ipDe),1)
        call dcopy_(n2dens,[0.0d0],0,Work(ipP),1)
        Call Densi2(2,Work(ipDe),Work(ipP),
     &               Work(ipL),Work(ipR),0,0,0,n1dens,n2dens)
*
*        write(*,*)'ipDe*ipDe',
*     &   ddot_(n1dens,Work(ipDe),1,Work(ipDe),1)
*        Call RecPrt('ipDe',' ',Work(ipDe),n1dens,1)
*
        call dcopy_(n2dens,Work(ipp),1,rp,1)
        call dcopy_(n1dens,Work(ipde),1,rD,1)
*
*        write(*,*)'rD*rD',
*     &   ddot_(n1dens,rD,1,rD,1)
*        Call RecPrt('iprD',' ',rD,n1dens,1)
*
        irc=ipin(iCI)
        irc=ipin(ipCI)
        Call CSF2SD(W(iCI)%Vec(1+nconf1),Work(ipL),iS)
        Call CSF2SD(W(ipci)%Vec,Work(ipR),State_SYM)
*
*        write(*,*)'ipL*ipL',
*     &   ddot_(nConfL,Work(ipL),1,Work(ipL),1)
*        write(*,*)'ipR*ipR',
*     &   ddot_(nConfR,Work(ipR),1,Work(ipR),1)
*        Call RecPrt('ipL',' ',Work(ipL),nConfL,1)
*        Call RecPrt('ipR',' ',Work(ipR),nConfR,1)
*
        irc=ipnout(-1)
        issm=iS
        icsm=STATE_SYM
        call dcopy_(n1dens,[0.0d0],0,Work(ipDe),1)
        call dcopy_(n2dens,[0.0d0],0,Work(ipP),1)
        Call Densi2(2,Work(ipDe),Work(ipP),Work(ipl),Work(ipr),
     &               0,0,0,n1dens,n2dens)
*
*        write(*,*)'ipDe*ipDe',
*     &   ddot_(n1dens,Work(ipDe),1,Work(ipDe),1)
*        Call RecPrt('ipDe',' ',Work(ipDe),n1dens,1)
*
C
        call daxpy_(n2Dens,-1.0d0,Work(ipP),1,rp,1)
        call daxpy_(n1Dens,-1.0d0,Work(ipDe),1,rD,1)
C
*        call dscal_(n2dens,-1.0d0,rP,1)
*        call dscal_(n1dens,-1.0d0,rD,1)
C
*
*        write(*,*)'rD*rD',
*     &   ddot_(n1dens,rD,1,rD,1)
*        Call RecPrt('iprD',' ',rD,n1dens,1)
*
        Call GetMem('CIL','FREE','REAL',ipL,nConfL)
        Call GetMem('CIR','FREE','REAL',ipR,nConfR)
        Do i=1,ntash
        Do j=1,ntash
        Do k=1,ntash
        Do l=1,ntash
        ijkl=itri(ntash*(j-1)+i,k+(l-1)*ntash)
        !Write(6,'(I1,I1,I1,I1,F12.6)') i,j,k,l,rp(ijkl)
        End DO
        End DO
        End DO
        End DO
      End If
      Call GetMem('1Dens2','Free','Real',ipDe,n1dens)
      Call GetMem('2Dens2','Free','Real',ipP,n2dens)
      Return
      End
