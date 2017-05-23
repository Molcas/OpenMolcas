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
      SubRoutine CIDens(response,iLS,iRS,iL,iR,iS,rP,rD)
      Implicit Real*8(a-h,o-z)

#include "detdim.fh"
#include "cicisp_mclr.fh"
#include "WrkSpc.fh"
#include "crun_mclr.fh"

#include "Input.fh"
#include "Pointers.fh"
#include "spinfo_mclr.fh"
#include "cands.fh"
#include "dmrginfo_mclr.fh"
      Real*8 rP(*),rD(*)
      integer opout
      Logical Response
      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
* LS = CI
*
*     Ok, we once more want to hide Jeppe's routines from
*     the eyes of the world, so everyone believes that I have done
*     all the work.
*     If we have spin dependent Hamiltonian we will work with
*     SD in all parts of the program, no CSF is necessary,
*     otherwise we will do the optimization in CSF's to increase
*     convergence.
*
*     Input:
*
*     response: true if the density should be used for response calculations
*     iLS: CI Coeff for left state
*     iRS: CI Coeff for right state
*     iL : Symmetry of left state
*     iR : Symmetry of right state
*
*
*               +       +
*     iS=1 E  =a  a  + a  a
*           pq  ap aq   Bp Bq
*
*                +       +
*     iS=-1 T  =a  a  - a a
*            pq  ap aq   Bp Bq
*
*     Output:
*
*      rP : Two Electron Density
*      rD : One Electron Density
*
*

      If (nconf1.eq.0) return

      if(doDMRG)then  !yma
        call dmrg_dim_change_mclr(LRras2(1:8),ndim,0)
        call dmrg_dim_change_mclr(LRras2(1:8),nna,0)
        n1dens=ndim**2
        n2dens=n1dens*(n1dens+1)/2
      end if

      Call GetMem('1Dens2','ALLO','Real',ipDe,n1dens)
      Call GetMem('2Dens2','ALLO','Real',ipP,n2dens)
      call dcopy_(n1dens,0.0d0,0,rD,1)
      call dcopy_(n2dens,0.0d0,0,rP,1)
      If (nocsf.eq.0) Then
        nConfL=Max(ncsf(il),nint(xispsm(il,1)))
        nConfR=Max(ncsf(iR),nint(xispsm(iR,1)))
        Call GetMem('CIL','ALLO','REAL',ipL,nConfL)
        Call CSF2SD(Work(ipin1(iLS,nconfL)),Work(ipL),iL)
        irc=opout(ils)
        Call GetMem('CIR','ALLO','REAL',ipR,nConfR)
        Call CSF2SD(Work(ipin1(iRS,nconfR)),Work(ipR),iR)
        irc=opout(irs)
        irc=ipnout(-1)
        icsm=iR
        issm=iL
        Call Densi2(2,Work(ipDe),Work(ipP),
     &               Work(ipL),Work(ipR),0,0,0,n1dens,n2dens)

        If (.not.timedep) Then
         If (response) Then
          Do iA=1,nnA
           Do jA=1,nnA
            Do kA=1,nnA
             Do la=1,nnA
              ij1=nnA*(iA-1)+ja
              ij2=nna*(ja-1)+ia
              kl1=nnA*(ka-1)+la
              kl2=nna*(la-1)+ka
              rp(itri(ij1,kl1))=
     &         Work(ipp-1+itri(ij1,kl1))+work(ipp-1+itri(ij2,kl2))
             End Do
            End Do
           End Do
          End Do
          Do iA=1,nnA
           Do jA=1,nnA
              ij1=nnA*(iA-1)+ja
              ij2=nna*(ja-1)+ia
              rD(ij1)=Work(ipDe-1+ij1)+work(ipDe-1+ij2)
           End Do
          End Do

         Else
          call dcopy_(n2dens,Work(ipp),1,rp,1)
          call dcopy_(n1dens,Work(ipde),1,rD,1)
         End If
        Else
          call dcopy_(n2dens,Work(ipp),1,rp,1)
          call dcopy_(n1dens,Work(ipde),1,rD,1)
          iCSM=iL
          iSSM=iR
          Call Densi2(2,Work(ipDe),Work(ipP),Work(ipR),Work(ipL),
     &               0,0,0,n1dens,n2dens)
          call daxpy_(n2Dens,-1.0d0,Work(ipP),1,rp,1)
          call daxpy_(n1Dens,-1.0d0,Work(ipDe),1,rD,1)
         End If
        Call GetMem('CIL','FREE','REAL',ipL,nConfL)
        Call GetMem('CIR','FREE','REAL',ipR,nConfR)
      else
        issm=iL
        icsm=iR
        Call Densi2(2,Work(ipDe),Work(ipP), Work(ipin(iLS)),
     &              Work(ipin(iRS)), 0,0,0,
     &              n1dens,n2dens)
        If (.not.timedep) Then
         If (response) Then
          Do iA=1,nnA
           Do jA=1,nnA
            Do kA=1,nnA
             Do la=1,nnA
              ij1=nnA*(iA-1)+ja
              ij2=nna*(ja-1)+ia
              kl1=nnA*(ka-1)+la
              kl2=nna*(la-1)+ka
              rp(itri(ij1,kl1))=
     &         Work(ipp-1+itri(ij1,kl1))+work(ipp-1+itri(ij2,kl2))
             End Do
            End Do
           End Do
          End Do
          Do iA=1,nnA
           Do jA=1,nnA
              ij1=nnA*(iA-1)+ja
              ij2=nna*(ja-1)+ia
              rD(ij1)=Work(ipDe-1+ij1)+work(ipDe-1+ij2)
           End Do
          End Do
         Else
          call dcopy_(n2dens,Work(ipp),1,rp,1)
          call dcopy_(n1dens,Work(ipde),1,rD,1)
         End If
        Else
          call dcopy_(n2dens,Work(ipp),1,rp,1)
          call dcopy_(n1dens,Work(ipde),1,rD,1)
          iCSM=iL
          iSSM=iR
          Call Densi2(2,Work(ipDe),Work(ipP),
     &                Work(ipin(iRS)),Work(ipin(ils)),
     &                0,0,0,n1dens,n2dens)
          call daxpy_(n2Dens,-1.0d0,Work(ipP),1,rp,1)
          call daxpy_(n1Dens,-1.0d0,Work(ipDe),1,rD,1)
         End If
      End If
      Call GetMem('1Dens2','Free','Real',ipDe,n1dens)
      Call GetMem('2Dens2','Free','Real',ipP,n2dens)

      if(doDMRG)then  ! yma
        call dmrg_dim_change_mclr(RGras2(1:8),ndim,0)
        call dmrg_dim_change_mclr(RGras2(1:8),nna,0)
        n1dens=ndim**2
        n2dens=n1dens*(n1dens+1)/2
      end if

      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(iS)
      End
