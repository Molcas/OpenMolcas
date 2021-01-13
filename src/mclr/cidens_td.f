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
#include "stdalloc.fh"
#include "crun_mclr.fh"
#include "Input.fh"
#include "Pointers.fh"
#include "spinfo_mclr.fh"
#include "cands.fh"
      Real*8 rP(*),rD(*)
      Real*8, Allocatable:: De(:), Pe(:), CIL(:), CIR(:)

#ifdef _DEBUGPRINT_
      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
#endif

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

      Call mma_allocate(De,2*n1dens,Label='De')
      Call mma_allocate(Pe,3*n2dens,Label='Pe')
      call dcopy_(n1dens,[0.0d0],0,rD,1)
      call dcopy_(n2dens,[0.0d0],0,rP,1)

      If (nocsf.eq.0) Then
        nConfL=Max(ncsf(iS),nint(xispsm(iS,1)))
        nConfR=Max(ncsf(State_SYM),nint(xispsm(STATE_SYM,1)))
        nC=Max(nconfL,nconfR)
        Call mma_allocate(CIL,nC,Label='CIL')
        Call mma_allocate(CIR,nC,Label='CIR')
c
c iCI is as long as ipcid for the timedep case!
c
        irc=ipin(iCI)
        irc=ipin(ipCI)
        Call CSF2SD(W(iCI)%Vec,CIR,iS)
        Call CSF2SD(W(ipCI)%Vec,CIL,State_SYM)
*
*        write(*,*)'ipL*ipL',ddot_(nConfL,CIL,1,CIL,1)
*        write(*,*)'ipR*ipR',ddot_(nConfR,CIR,1,CIR,1)
*        Call RecPrt('CIL',' ',CIL,nConfL,1)
*        Call RecPrt('CIR',' ',CIR,nConfR,1)
*
        irc=ipnout(-1)
        icsm=iS
        issm=STATE_SYM
*
*       <P|E_pq|0> & <P|e_pqrs|0> -> ipDe & ipP
*       ipL is the bra side vector
        call dcopy_(n1dens,[0.0d0],0,De,1)
        call dcopy_(n2dens,[0.0d0],0,Pe,1)
        Call Densi2(2,De,Pe,CIL,CIR,0,0,0,n1dens,n2dens)
*
*        write(*,*)'De*De',ddot_(n1dens,De,1,De,1)
*        Call RecPrt('De',' ',De,n1dens,1)
*
        call dcopy_(n2dens,Pe,1,rp,1)
        call dcopy_(n1dens,De,1,rD,1)
*
*        write(*,*)'rD*rD',ddot_(n1dens,rD,1,rD,1)
*        Call RecPrt('rD',' ',rD,n1dens,1)
*
        irc=ipin(iCI)
        irc=ipin(ipCI)
        Call CSF2SD(W(iCI)%Vec(1+nconf1),CIL,iS)
        Call CSF2SD(W(ipci)%Vec,CIR,State_SYM)
*
*        write(*,*)'CIL*CIL',ddot_(nConfL,CIL,1,CIL,1)
*        write(*,*)'CIR*CIR',ddot_(nConfR,CIR,1,CIR,1)
*        Call RecPrt('CIL',' ',CIL,nConfL,1)
*        Call RecPrt('CIR',' ',CIR,nConfR,1)
*
        irc=ipnout(-1)
        issm=iS
        icsm=STATE_SYM
        call dcopy_(n1dens,[0.0d0],0,De,1)
        call dcopy_(n2dens,[0.0d0],0,Pe,1)
        Call Densi2(2,De,Pe,CIL,CIR,0,0,0,n1dens,n2dens)
*
*        write(*,*)'De*De',&   ddot_(n1dens,De,1,De,1)
*        Call RecPrt('De',' ',De,n1dens,1)
*
C
        call daxpy_(n2Dens,-1.0d0,Pe,1,rp,1)
        call daxpy_(n1Dens,-1.0d0,De,1,rD,1)
C
*        call dscal_(n2dens,-1.0d0,rP,1)
*        call dscal_(n1dens,-1.0d0,rD,1)
C
*
*        write(*,*)'rD*rD',ddot_(n1dens,rD,1,rD,1)
*        Call RecPrt('rD',' ',rD,n1dens,1)
*
        Call mma_deallocate(CIL)
        Call mma_deallocate(CIR)

#ifdef _DEBUGPRINT_
        Do i=1,ntash
        Do j=1,ntash
        Do k=1,ntash
        Do l=1,ntash
        ijkl=itri(ntash*(j-1)+i,k+(l-1)*ntash)
        Write(6,'(I1,I1,I1,I1,F12.6)') i,j,k,l,rp(ijkl)
        End DO
        End DO
        End DO
        End DO
#endif
      End If

      Call mma_deallocate(Pe)
      Call mma_deallocate(De)

      Return
#ifdef _WARNING_WORKAROUND_
      If (.False.) Call Unused_integer(irc)
#endif
      End
