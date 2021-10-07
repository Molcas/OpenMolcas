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
      SubRoutine CIDens_sa(RSP,iLS,iRS,iL,iR,rP,rD)

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
      Logical RSP
      integer opout
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
*     RSP: true if the density should be used for response calculations
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

      if(doDMRG)then
        call dmrg_dim_change_mclr(LRras2(1:8),ndim,0)
        call dmrg_dim_change_mclr(LRras2(1:8),nna,0)
        n1dens=ndim**2
        n2dens=n1dens*(n1dens+1)/2
      end if

      Call GetMem('1Dens2','ALLO','Real',ipDe,n1dens)
      Call GetMem('2Dens2','ALLO','Real',ipP,n2dens)
      call dcopy_(n1dens,[0.0d0],0,rD,1)
      call dcopy_(n2dens,[0.0d0],0,rP,1)
*
      nConfL=Max(ncsf(il),nint(xispsm(il,1)))
      nConfR=Max(ncsf(iR),nint(xispsm(iR,1)))

      Call GetMem('CIL','ALLO','REAL',ipL,nConfL)
      Call GetMem('CIR','ALLO','REAL',ipR,nConfR)
*
      Do i=0,nroots-1
        Call CSF2SD(Work(ipin(iLS)+i*ncsf(il)),Work(ipL),iL)
        irc=opout(ils)
        Call CSF2SD(Work(ipin(iRS)+i*ncsf(ir)),Work(ipR),iR)
        irc=opout(irs)
        irc=ipnout(-1)
        icsm=iR
        issm=iL
        Call Densi2(2,Work(ipDe),Work(ipP),
     &               Work(ipL),Work(ipR),0,0,0,n1dens,n2dens)

        If (RSP) Then
           Do iA=1,nnA
             Do jA=1,nnA
               Do kA=1,nnA
                Do la=1,nnA
                 ij1=nnA*(iA-1)+ja
                 ij2=nna*(ja-1)+ia
                 kl1=nnA*(ka-1)+la
                 kl2=nna*(la-1)+ka
                 if (ij1.ge.kl1)
     &           rp(itri(ij1,kl1))=rp(itri(ij1,kl1))+weight(1+i)*(
     &            Work(ipp-1+itri(ij1,kl1))+work(ipp-1+itri(ij2,kl2)))
                End Do
               End Do
             End Do
           End Do
           Do iA=1,nnA
              Do jA=1,nnA
                 ij1=nnA*(iA-1)+ja
                 ij2=nna*(ja-1)+ia
                 rD(ij1)=rD(ij1)+
     *                weight(1+i)*(Work(ipDe-1+ij1)+work(ipDe-1+ij2))
              End Do
           End Do
        Else
           call daxpy_(n2dens,weight(i+1),Work(ipp),1,rp,1)
           call daxpy_(n1dens,Weight(i+1),Work(ipDe),1,rD,1)
        End If
      End Do
*
      Call GetMem('CIL','FREE','REAL',ipL,nConfL)
      Call GetMem('CIR','FREE','REAL',ipR,nConfR)
      Call GetMem('1Dens2','Free','Real',ipDe,n1dens)
      Call GetMem('2Dens2','Free','Real',ipP,n2dens)

      if(doDMRG)then
        call dmrg_dim_change_mclr(RGras2(1:8),ndim,0)
        call dmrg_dim_change_mclr(RGras2(1:8),nna,0)
        n1dens=ndim**2
        n2dens=n1dens*(n1dens+1)/2
      end if
*
      Return
      End
