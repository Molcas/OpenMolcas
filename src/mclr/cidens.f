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
      use ipPage, only: W
      Implicit Real*8(a-h,o-z)

#include "real.fh"
#include "detdim.fh"
#include "cicisp_mclr.fh"
#include "stdalloc.fh"
#include "crun_mclr.fh"

#include "Input.fh"
#include "Pointers.fh"
#include "spinfo_mclr.fh"
#include "cands.fh"
#include "dmrginfo_mclr.fh"
      Real*8 rP(*),rD(*)
      integer opout
      Logical Response
      Real*8, Allocatable:: De(:), Pe(:), CIL(:), CIR(:)

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

      Call mma_allocate(De,n1dens,Label='De')
      Call mma_allocate(Pe,n2dens,Label='Pe')
      De(:)=Zero
      Pe(:)=Zero

      If (nocsf.eq.0) Then
        nConfL=Max(ncsf(il),nint(xispsm(il,1)))
        nConfR=Max(ncsf(iR),nint(xispsm(iR,1)))
        Call mma_allocate(CIL,nConfL,Label='CIL')
        irc=ipin1(iLS,nconfL)
        Call CSF2SD(W(iLS)%Vec,CIL,iL)
        irc=opout(iLs)
        Call mma_allocate(CIR,nConfR,Label='CIR')
        irc=ipin1(iRS,nconfR)
        Call CSF2SD(W(iRS)%Vec,CIR,iR)
        irc=opout(irs)
        irc=ipnout(-1)
        icsm=iR
        issm=iL
        Call Densi2(2,De,Pe,CIL,CIR,0,0,0,n1dens,n2dens)

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
     &         Pe(itri(ij1,kl1))+Pe(itri(ij2,kl2))
             End Do
            End Do
           End Do
          End Do
          Do iA=1,nnA
           Do jA=1,nnA
              ij1=nnA*(iA-1)+ja
              ij2=nna*(ja-1)+ia
              rD(ij1)=De(ij1)+De(ij2)
           End Do
          End Do

         Else
          call dcopy_(n2dens,Pe,1,rp,1)
          call dcopy_(n1dens,De,1,rD,1)
         End If
        Else
          call dcopy_(n2dens,Pe,1,rp,1)
          call dcopy_(n1dens,De,1,rD,1)
          iCSM=iL
          iSSM=iR
          Call Densi2(2,De,Pe,CIR,CIL,0,0,0,n1dens,n2dens)
          call daxpy_(n2Dens,-One,Pe,1,rp,1)
          call daxpy_(n1Dens,-One,De,1,rD,1)
         End If
         Call mma_deallocate(CIL)
         Call mma_deallocate(CIR)
      else
        issm=iL
        icsm=iR
        irc=ipin(iLS)
        irc=ipin(iRS)
        Call Densi2(2,De,Pe,W(iLS)%Vec,W(iRS)%Vec,0,0,0,n1dens,n2dens)
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
     &         Pe(itri(ij1,kl1))+Pe(itri(ij2,kl2))
             End Do
            End Do
           End Do
          End Do
          Do iA=1,nnA
           Do jA=1,nnA
              ij1=nnA*(iA-1)+ja
              ij2=nna*(ja-1)+ia
              rD(ij1)=De(ij1)+De(ij2)
           End Do
          End Do
         Else
          call dcopy_(n2dens,Pe,1,rp,1)
          call dcopy_(n1dens,De,1,rD,1)
         End If
        Else
          call dcopy_(n2dens,Pe,1,rp,1)
          call dcopy_(n1dens,De,1,rD,1)
          iCSM=iL
          iSSM=iR
          irc=ipin(iRS)
          irc=ipin(iLS)
          Call Densi2(2,De,Pe,W(iRS)%Vec,W(iLS)%Vec,0,0,0,n1dens,n2dens)
          call daxpy_(n2Dens,-One,Pe,1,rp,1)
          call daxpy_(n1Dens,-One,De,1,rD,1)
         End If
      End If
      Call mma_deallocate(Pe)
      Call mma_deallocate(De)

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
