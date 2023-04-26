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
      Subroutine RHS_NAC(Fock,SLag)
      use ipPage, only: W
      Implicit None
#include "Input.fh"
#include "Pointers.fh"
#include "stdalloc.fh"
#include "real.fh"
#include "sa.fh"
#include "detdim.fh"
#include "cicisp_mclr.fh"
#include "cands.fh"
      Real*8 Fock(*)
      real*8, optional :: SLag(*)
      Integer ng1,ng2,i,j,k,l,ij,kl,ijkl,ij2,kl2,ijkl2
      Integer iTri
      Integer ipIn,opOut,ipnOut,nConfL,nConfR,iRC,LuDens
      Real*8 factor
      External ipIn,opOut,ipnOut
*
      Integer iSLag !,jR,kR
      Real*8, Allocatable:: G1q(:), G1m(:), G1r(:), G2q(:), G2r(:),
     &                      CIL(:), CIR(:), T(:), F(:)
*                                                                      *
************************************************************************
*                                                                      *
      iTri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
*                                                                      *
************************************************************************
*                                                                      *
      ng1=ntAsh*(ntAsh+1)/2
      ng2=ng1*(ng1+1)/2
      If (PT2) Then
        Call mma_allocate(G1q,n1dens,Label='G1q')
        Call mma_allocate(G2q,n2dens,Label='G2q')
      Else
        Call mma_allocate(G1q,nG1,Label='G1q')
        Call mma_allocate(G2q,nG2,Label='G2q')
      End If
      Call mma_allocate(G1m,nG1,Label='G1m')
      Call mma_allocate(G1r,n1dens,Label='G1r')
      G1r(:)=Zero
      Call mma_allocate(G2r,n2dens,Label='G2r')
      G2r(:)=Zero
*
**    Calculate one- and two-particle transition matrices
**    from the CI vectors of the two NAC states
**    (code copied from CIdens_SA, same symmetry)
*
      nConfL=Max(nconf1,nint(xispsm(State_sym,1)))
      nConfR=Max(nconf1,nint(xispsm(State_sym,1)))
      Call mma_allocate(CIL,nConfL,Label='CIL')
      Call mma_allocate(CIR,nConfR,Label='CIR')
      If (PT2) Then
        Call PT2_SLag()
      Else
        irc=ipIn(ipCI)
        Call CSF2SD(W(ipCI)%Vec(1+(NSSA(2)-1)*nconf1),CIL,State_sym)
        iRC=opout(ipCI)
        Call CSF2SD(W(ipCI)%Vec(1+(NSSA(1)-1)*nconf1),CIR,State_sym)
        iRC=opout(ipCI)
        iRC=ipnout(-1)
        icsm=1
        issm=1
        Call Densi2(2,G1r,G2r,
     &                CIL,CIR,0,0,0,n1dens,n2dens)
      End If
      Call mma_deallocate(CIL)
      Call mma_deallocate(CIR)
*
**    Symmetrize densities
**    For the one-particle density, save the antisymmetric part too
*
      ij=0
      Do i=0,ntAsh-1
        Do j=0,i-1
          ij=ij+1
          G1q(ij)=(G1r(1+i*ntAsh+j)+
     &                    G1r(1+j*ntAsh+i))*Half
*         Note that the order of subtraction depends on how the matrix
*         will be used when contracting with derivative integrals
*         This is found to give the correct results:
          G1m(ij)=(G1r(1+j*ntAsh+i)-
     &                    G1r(1+i*ntAsh+j))*Half
        End Do
        ij=ij+1
        G1q(ij)=G1r(1+i*ntAsh+i)
        G1m(ij)=Zero
      End Do
*
      !! The anti-symmetric RDM is contructed somewhere in the CASPT2
      !! module. It will be read from disk in out_pt2.f.
      If (PT2) Call DCopy_(ng1,[zero],0,G1m,1)
*
      Do i=1,ntAsh**2
        j=itri(i,i)
        G2r(j)=Half*G2r(j)
      End Do
      Do i=0,ntAsh-1
        Do j=0,i-1
          ij=i*(i+1)/2+j
          Do k=0,ntAsh-1
            Do l=0,k
              kl=k*(k+1)/2+l
              If (ij.ge.kl) Then
                factor=Quart
                If (ij.eq.kl) factor=Half
                ijkl=ij*(ij+1)/2+kl
                ij2=i*ntAsh+j
                kl2=k*ntAsh+l
                G2q(1+ijkl)=factor*G2r(1+ij2*(ij2+1)/2+kl2)
                ij2=Max(j*ntAsh+i,l*ntAsh+k)
                kl2=Min(j*ntAsh+i,l*ntAsh+k)
                G2q(1+ijkl)=G2q(1+ijkl)+
     &                           factor*G2r(1+ij2*(ij2+1)/2+kl2)
                If (k.ne.l) Then
                  ij2=i*ntAsh+j
                  kl2=l*ntAsh+k
                  G2q(1+ijkl)=G2q(1+ijkl)+
     &                             factor*G2r(1+ij2*(ij2+1)/2+kl2)
                  If (ij.ne.kl) Then
                    ij2=Max(j*ntAsh+i,k*ntAsh+l)
                    kl2=Min(j*ntAsh+i,k*ntAsh+l)
                    G2q(1+ijkl)=G2q(1+ijkl)+
     &                              factor*G2r(1+ij2*(ij2+1)/2+kl2)
                  End If
                End If
              End If
            End Do
          End Do
        End Do
        ij=i*(i+1)/2+i
        Do k=0,ntAsh-1
          Do l=0,k
            kl=k*(k+1)/2+l
            If (ij.ge.kl) Then
              factor=Half
              If (ij.eq.kl) factor=One
              ijkl=ij*(ij+1)/2+kl
              ij2=i*ntAsh+i
              kl2=k*ntAsh+l
              G2q(1+ijkl)=factor*G2r(1+ij2*(ij2+1)/2+kl2)
              If (k.ne.l) Then
                kl2=l*ntAsh+k
                G2q(1+ijkl)=G2q(1+ijkl)+
     &                           factor*G2r(1+ij2*(ij2+1)/2+kl2)
              End If
            End If
          End Do
        End Do
      End Do
*
**    Write the symmetric densities in the runfile and
**    the antisymmetric one-particle in MCLRDENS
*
      iRC=0
      LuDens=20
      Call DaName(LuDens,'MCLRDENS')
      Call dDaFile(LuDens,1,G1m,ng1,iRC)
      Call DaClos(LuDens)
      Call Put_dArray('D1mo',G1q,ng1)
      Call Put_dArray('P2mo',G2q,ng2)
*
**    Store transition Fock matrix
*
      Do i=1,ntAsh
        Do j=1,ntAsh
          G1r(ntAsh*(j-1)+i)=G1q(iTri(i,j))
        End Do
      End Do
      Do i=1,ntAsh
        Do j=1,ntAsh
          ij=iTri(i,j)
          ij2=ntAsh*(i-1)+j
          Do k=1,ntAsh
            Do l=1,ntAsh
              kl=iTri(k,l)
              kl2=ntAsh*(k-1)+l
              factor=One
              If (ij.ge.kl .and. k.eq.l) factor=Two
              If (ij.lt.kl .and. i.eq.j) factor=Two
              ijkl=iTri(ij,kl)
              ijkl2=iTri(ij2,kl2)
              G2r(ijkl2)=factor*G2q(ijkl)
            End Do
          End Do
        End Do
      End Do

* Note: 1st arg = zero for no inactive density (TDM)
      Call mma_allocate(T,nDens2,Label='T')
      Call mma_allocate(F,nDens2,Label='F')
      Call FockGen(Zero,G1r,G2r,T,Fock,1)
      Call TCMO(T,1,-2)
      ij=0
      Do k=1,nSym
        Do i=0,nBas(k)-1
          Do j=0,i-1
            ij=ij+1
            F(ij)=T(ipMat(k,k)+nBas(k)*j+i)+
     &                   T(ipMat(k,k)+nBas(k)*i+j)
          End Do
          ij=ij+1
          F(ij)=T(ipMat(k,k)+nBas(k)*i+i)
        End Do
      End Do
      Call Put_dArray('FockOcc',F,nDens2)
*
      Call mma_deallocate(T)
      Call mma_deallocate(F)
      Call mma_deallocate(G1r)
      Call mma_deallocate(G2r)
      Call mma_deallocate(G2q)
      Call mma_deallocate(G1m)
      Call mma_deallocate(G1q)
*
      Return

       Contains

      Subroutine PT2_SLag
C
C     Almost the same to the subroutine in rhs_sa.f,
C     but slightly modified
C
      Implicit Real*8 (A-H,O-Z)
      ! integer opout
      integer jR,kR
C
      !! iR = iRLXRoot
      Do jR = 1, nRoots
        Do kR = 1, jR
          vSLag = 0.0D+00
C         write (*,*) "jr,kr= ", jr,kr
C         write (*,*) vslag
          iSLag = jR + nRoots*(kR-1)
          vSLag = SLag(iSLag)
C         write (*,*) vslag
C
          Call CSF2SD(W(ipCI)%Vec(1+(jR-1)*nconf1),CIL,1)
          ! iRC=opout(ipCI)
          Call CSF2SD(W(ipCI)%Vec(1+(kR-1)*nconf1),CIR,1)
          ! iRC=opout(ipCI)
          ! iRC=ipnout(-1)
          ! icsm=1
          ! issm=1
C
          If (abs(vSLag).gt.1.0d-10) Then
            Call Densi2(2,G1q,G2q,CIL,CIR,0,0,0,n1dens,n2dens)
            Call DaXpY_(n1dens,vSLag,G1q,1,G1r,1)
            Call DaXpY_(n2dens,vSLag,G2q,1,G2r,1)
          End If
C
          If (kR.ne.jR) Then
            iSLag = kR + nRoots*(jR-1)
            vSLag = SLag(iSLag)
            If (abs(vSLag).gt.1.0d-10) Then
              Call Densi2(2,G1q,G2q,CIR,CIL,0,0,0,n1dens,n2dens)
              Call DaXpY_(n1dens,vSLag,G1q,1,G1r,1)
              Call DaXpY_(n2dens,vSLag,G2q,1,G2r,1)
            End If
          End If
        End Do
      End Do
      nConf=ncsf(1) !! nconf is overwritten somewhere in densi2
C
      Return
C
      End Subroutine PT2_SLag
      End
