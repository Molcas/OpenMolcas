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
       Subroutine rhs_sa(Fock,SLag)
       use Arrays, only: Int1
       use ipPage, only: W
       Implicit Real*8 (a-h,o-z)

#include "Input.fh"
#include "Pointers.fh"
#include "stdalloc.fh"
#include "SysDef.fh"
#include "Files_mclr.fh"
#include "real.fh"
#include "sa.fh"
#include "dmrginfo_mclr.fh"
#include "detdim.fh"
#include "cicisp_mclr.fh"

       Real*8 Fock(*)
       real*8, optional :: SLag(*)
       Dimension rdum(1)
       Real*8, Allocatable:: T(:), F(:), G1q(:), G2q(:), G1r(:), G2r(:)
*                                                                      *
************************************************************************
*                                                                      *
       itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
*                                                                      *
************************************************************************
*                                                                      *
*
       if(doDMRG)then ! yma
         call dmrg_dim_change_mclr(RGras2(1:8),ntash,0)
       end if

       ng1=itri(ntash,ntash)
       ng2=itri(ng1,ng1)
*
       Call mma_allocate(T,ndens2,Label='T')
       Call mma_allocate(F,ndens2,Label='F')
       Call mma_allocate(G1q,ng1,Label='G1q')
       Call mma_allocate(G2q,ng2,Label='G2q')
       Call mma_allocate(G1r,ntash**2,Label='G1r')
       Call mma_allocate(G2r,itri(ntash**2,ntash**2),Label='G2r')
*
**     Pick up densities from JobIph file
*
       iR=iroot(istate)
       jdisk=itoc(3)
       Do ii=1,iR-1
         Call dDaFile(LUJOB ,0,rdum,ng1,jDisk)
         Call dDaFile(LUJOB ,0,rdum,ng1,jDisk)
         Call dDaFile(LUJOB ,0,rdum,Ng2,jDisk)
         Call dDaFile(LUJOB ,0,rdum,Ng2,jDisk)
       End Do
       Call dDaFile(LUJOB ,2,G1q,ng1,jDisk)
       Call dDaFile(LUJOB ,0,rdum,ng1,jDisk)
       Call dDaFile(LUJOB ,2,G2q,Ng2,jDisk)
       Call dDaFile(LUJOB ,0,rdum,Ng2,jDisk)
*
       !! Add SLag (rotations of states) contributions from the partial
       !! derivative of the CASPT2 energy. G1q and G2q are modified.
       !! The modified density will be used in out_pt2.f and ptrans_sa.f
       If (PT2.and.nRoots.gt.1) Call PT2_SLag()
*
       Call Put_dArray('P2mo',G2q,ng2)
       Call Put_dArray('D1mo',G1q,ng1)
*
       Do iB=1,ntash
        Do jB=1,ntash
        G1r(ib+(jb-1)*ntash) = G1q(itri(ib,jb))
        End Do
       End Do

       Do iB=1,ntash
        Do jB=1,ntash
         iDij=iTri(ib,jB)
         iRij=jb+(ib-1)*ntash
         Do kB=1,ntash
          Do lB=1,ntash
           iDkl=iTri(kB,lB)
           iRkl=lb+(kb-1)*ntash
           fact=One
           if(iDij.ge.iDkl .and. kB.eq.lB) fact=Two
           if(iDij.lt.iDkl .and. iB.eq.jB) fact=Two
           iijkl=itri(iDij,iDkl)
           iRijkl=itri(iRij,iRkl)
           G2r(iRijkl)=Fact*G2q(iijkl)
          End Do
         End Do
        End Do
       End Do
*
       if(doDMRG)then ! yma
         call dmrg_dim_change_mclr(RGras2(1:8),nna,0)
       end if

       Call FockGen(One,G1r,G2r,T,Fock,1)
*       Do iS=1,nsym
*        Call RecPrt(' ',' ',fock(ipMat(is,is)),nbas(is),nbas(is))
*       End Do

      If (.not.debug) Then !yma debug ??
       renergy=Zero
       Do ii=1,nsym
        Do jj=1,nbas(ii)
         renergy = renergy + T(ipmat(ii,ii)+jj-1+nbas(ii)*(jj-1))
        End DO
       End DO

      rcorei=Zero
      rcorea=Zero
      Do iS=1,nSym
       Do iB=1,nIsh(is)
       rcorei=rcorei+Two*Int1(ipCM(is)+nOrb(iS)*(ib-1)+ib-1)
       End Do

       Do iB=1,nAsh(iS)
        Do jB=1,nAsh(iS)
         iiB=nA(iS)+ib
         ijB=nA(iS)+jb
         iij=iTri(iib,ijb)
         iiB=nIsh(iS)+ib
         ijB=nIsh(iS)+jb
         rcorea=rcorea+G1q(iij)*Int1(ipCM(is)+nOrb(is)*(iib-1)+ijB-1)
        End Do
       End Do
      End Do
!      rcore=rCorei+rcoreA
!      write(*,*) 'In rhs_sa'
!      Write(*,*) 'Checking energy',0.5d0*renergy+potnuc+half*rcore !yma
!      Write(*,*) 'Checking energy',0.5d0*renergy,potnuc,rcore      !yma
!      write(*,*)
      End if
!      Do iS=1,nsym
!       Call RecPrt(' ',' ',fock(ipMat(is,is)),nbas(is),nbas(is))
!      End Do
*
       Call mma_deallocate(G1q)
       Call mma_deallocate(G2q)
*
       Call TCMO(T,1,-2)
       ijb=0
       Do is=1,nsym
        Do ib=1,nbas(is)
         Do jb=1,ib-1
          ijb=ijb+1
          F(ijb)=T(ipmat(is,is)+nbas(is)*(JB-1)+IB-1)
     &          +T(ipmat(is,is)+nbas(is)*(IB-1)+JB-1)
         End Do
         ijb=ijb+1
         F(ijb)=T(ipmat(is,is)+nbas(is)*(iB-1)+IB-1)
        End Do
       End Do
       Call Put_dArray('FockOcc',F,nDens2)

!       call recprt('RHS',' ',fock,ndens2,1)
*
       Call mma_deallocate(G1r)
       Call mma_deallocate(G2r)
       Call mma_deallocate(T)
       Call mma_deallocate(F)

*
       Return

       Contains

      Subroutine PT2_SLag

      Implicit Real*8 (A-H,O-Z)
      ! integer opout
      Real*8, Allocatable:: CIL(:), CIR(:)
      integer :: i,j

!     At present, Molcas accepts equally-weighted MCSCF reference,
!     so all SLag values are employed in the following computation.
!     For unequally-weighted reference as in GAMESS-US, some more
!     operations are required, but the CP-MCSCF part has to be
!     modified, so this may not be realized easily.

      nConfL=Max(nconf1,nint(xispsm(1,1)))
      nConfR=Max(nconf1,nint(xispsm(1,1)))
      call mma_allocate(CIL, nConfL, Label='CIL')
      call mma_allocate(CIR, nConfR, Label='CIR')
      !! iR = iRLXRoot
      Do jR = 1, nRoots
        Call CSF2SD(W(ipCI)%Vec(1+(jR-1)*nconf1),CIL,1)
        Do kR = 1, jR !! jR-1
          iSLag = jR + nRoots*(kR-1)
          vSLag = SLag(iSLag)
          If (abs(vSLag).le.1.0d-10) Cycle
C
          Call CSF2SD(W(ipCI)%Vec(1+(jR-1)*nconf1),CIL,1)
          ! iRC=opout(ipCI)
          Call CSF2SD(W(ipCI)%Vec(1+(kR-1)*nconf1),CIR,1)
          ! iRC=opout(ipCI)
          ! iRC=ipnout(-1)
          ! icsm=1
          ! issm=1
          Call Densi2(2,G1r,G2r,CIL,CIR,0,0,0,n1dens,n2dens)
          !! For RDM1
          ij=0
          Do i=0,ntAsh-1
            Do j=0,i-1
              ij=ij+1
              G1q(ij)=G1q(ij)+
     *          (G1r(1+i*ntAsh+j)+
     *           G1r(1+j*ntAsh+i))*Half*vSLag
            End Do
            ij=ij+1
            G1q(ij)=G1q(ij) + G1r(1+i*ntAsh+i)*vSLag
          End Do
          !! For RDM2
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
                    factor=Quart*vSLag
                    If (ij.eq.kl) factor=Half*vSLag
                    ijkl=ij*(ij+1)/2+kl
                    ij2=i*ntAsh+j
                    kl2=k*ntAsh+l
                    G2q(1+ijkl)=G2q(1+ijkl)
     *                + factor*G2r(1+ij2*(ij2+1)/2+kl2)
                    ij2=Max(j*ntAsh+i,l*ntAsh+k)
                    kl2=Min(j*ntAsh+i,l*ntAsh+k)
                    G2q(1+ijkl)=G2q(1+ijkl)
     &                + factor*G2r(1+ij2*(ij2+1)/2+kl2)
                    If (k.ne.l) Then
                      ij2=i*ntAsh+j
                      kl2=l*ntAsh+k
                      G2q(1+ijkl)=G2q(1+ijkl)
     &                  + factor*G2r(1+ij2*(ij2+1)/2+kl2)
                      If (ij.ne.kl) Then
                        ij2=Max(j*ntAsh+i,k*ntAsh+l)
                        kl2=Min(j*ntAsh+i,k*ntAsh+l)
                        G2q(1+ijkl)=G2q(1+ijkl)
     &                    + factor*G2r(1+ij2*(ij2+1)/2+kl2)
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
                  factor=Half*vSLag
                  If (ij.eq.kl) factor=One*vSLag
                  ijkl=ij*(ij+1)/2+kl
                  ij2=i*ntAsh+i
                  kl2=k*ntAsh+l
                  G2q(1+ijkl)=G2q(1+ijkl)
     *              + factor*G2r(1+ij2*(ij2+1)/2+kl2)
                  If (k.ne.l) Then
                    kl2=l*ntAsh+k
                    G2q(1+ijkl)=G2q(1+ijkl)
     &                + factor*G2r(1+ij2*(ij2+1)/2+kl2)
                  End If
                End If
              End Do
            End Do
          End Do
        End Do
      End Do
      call mma_deallocate(CIL)
      call mma_deallocate(CIR)
      nConf=ncsf(1)
C
      Return
C
      End Subroutine PT2_SLag
C
       End
