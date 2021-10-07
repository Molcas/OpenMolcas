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
       Implicit Real*8 (a-h,o-z)

#include "Input.fh"
#include "Pointers.fh"
#include "WrkSpc.fh"

#include "SysDef.fh"
#include "glbbas_mclr.fh"
#include "Files_mclr.fh"
#include "real.fh"
#include "sa.fh"
#include "dmrginfo_mclr.fh"
#include "detdim.fh"
#include "cicisp_mclr.fh"

       Real*8 Fock(*),SLag(*)
       Dimension rdum(1)
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
       Call Getmem('TEMP','ALLO','REAL',ipT,ndens2)
       Call Getmem('TEMP','ALLO','REAL',ipF,ndens2)
       Call Getmem('ONED','ALLO','REAL',ipG1q,ng1)
       Call Getmem('TWOD','ALLO','REAL',ipG2q,ng2)
       Call Getmem('ONED','ALLO','REAL',ipG1r,ntash**2)
       Call Getmem('TWOD','ALLO','REAL',ipG2r,itri(ntash**2,ntash**2))
*
**     Pick up densities from JobIph file
*
       iR=iroot(istate)
       jdisk=itoc(3)
       Do i=1,iR-1
         Call dDaFile(LUJOB ,0,rdum,ng1,jDisk)
         Call dDaFile(LUJOB ,0,rdum,ng1,jDisk)
         Call dDaFile(LUJOB ,0,rdum,Ng2,jDisk)
         Call dDaFile(LUJOB ,0,rdum,Ng2,jDisk)
       End Do
       Call dDaFile(LUJOB ,2,Work(ipG1q),ng1,jDisk)
       Call dDaFile(LUJOB ,0,rdum,ng1,jDisk)
       Call dDaFile(LUJOB ,2,Work(ipG2q),Ng2,jDisk)
       Call dDaFile(LUJOB ,0,rdum,Ng2,jDisk)
*
       !! Add SLag (rotations of states) contributions from the partial
       !! derivative of the CASPT2 energy. Work(ipG1q) and Work(ipG2q)
       !! are modified. The modified density will be used in out_pt2.f
       !! and ptrans_sa.f etc.
       If (PT2.and.nRoots.gt.1) Call PT2_SLag
*
       Call Put_P2MO(Work(ipG2q),ng2)
       Call Put_D1MO(Work(ipG1q),ng1)
*
       Do iB=1,ntash
        Do jB=1,ntash
        Work(ipG1r+ib-1+(jb-1)*ntash)=
     &     Work(ipg1q+itri(ib,jb)-1)
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
           Work(ipG2r-1+iRijkl)=Fact*Work(ipG2q+iijkl-1)
          End Do
         End Do
        End Do
       End Do
*
       if(doDMRG)then ! yma
         call dmrg_dim_change_mclr(RGras2(1:8),nna,0)
       end if

       Call FockGen(One,Work(ipG1r),Work(ipG2r),Work(ipT),Fock,1)
*       Do iS=1,nsym
*        Call RecPrt(' ',' ',fock(ipMat(is,is)),nbas(is),nbas(is))
*       End Do

      If (.not.debug) Then !yma debug ??
       renergy=Zero
       Do i=1,nsym
        Do j=1,nbas(i)
         renergy=renergy+
     &  Work(ipT-1+ipmat(i,i)+j-1+nbas(i)*(j-1))
        End DO
       End DO

      rcora=Zero
      rcorei=Zero
      rcorea=Zero
      Do iS=1,nSym
       Do iB=1,nIsh(is)
       rcorei=rcorei+Two*Work(kint1+ipCM(is)-1+
     &             nOrb(iS)*(ib-1)+ib-1)
       End Do

       Do iB=1,nAsh(iS)
        Do jB=1,nAsh(iS)
         iiB=nA(iS)+ib
         ijB=nA(iS)+jb
         iij=iTri(iib,ijb)
         iiB=nIsh(iS)+ib
         ijB=nIsh(iS)+jb
         rcorea=rcorea+Work(ipG1q+iij-1)*
     &           Work(kint1+ipCM(is)-1+nOrb(is)*(iib-1)+ijB-1)
        End Do
       End Do
      End Do
      rcore=rCorei+rcoreA
!      write(*,*) 'In rhs_sa'
!      Write(*,*) 'Checking energy',0.5d0*renergy+potnuc+half*rcore !yma
!      Write(*,*) 'Checking energy',0.5d0*renergy,potnuc,rcore      !yma
!      write(*,*)
      End if
!      Do iS=1,nsym
!       Call RecPrt(' ',' ',fock(ipMat(is,is)),nbas(is),nbas(is))
!      End Do
*
       Call Getmem('ONED','FREE','REAL',ipG1q,ng1)
       Call Getmem('TWOD','FREE','REAL',ipG2q,ng2)
*
       Call TCMO(Work(ipT),1,-2)
       ijb=0
       Do is=1,nsym
        Do ib=1,nbas(is)
         Do jb=1,ib-1
          Work(ipF+ijb)=Work(ipT+ipmat(is,is)-1+nbas(is)*(JB-1)+IB-1)
     &                 +work(ipT+ipmat(is,is)-1+nbas(is)*(IB-1)+JB-1)
          ijb=ijb+1
         End Do
         Work(ipF+ijb)=Work(ipT+ipmat(is,is)-1+nbas(is)*(iB-1)+IB-1)
         ijb=ijb+1
        End Do
       End Do
       Call Put_Fock_Occ(Work(ipf),nDens2)

!       call recprt('RHS',' ',fock,ndens2,1)
*
       Call Getmem('ONED','FREE','REAL',ipG1r,ntash**2)
       Call Getmem('TWOD','FREE','REAL',ipG2r,itri(ntash**2,ntash**2))
       Call Getmem('TEMP','FREE','REAL',ipT,ndens2)
       Call Getmem('TEMP','FREE','REAL',ipF,ndens2)

*
       Return

       Contains

      Subroutine PT2_SLag
C
      Implicit Real*8 (A-H,O-Z)
      integer opout
C
C     At present, Molcas accepts equally-weighted MCSCF reference,
C     so all SLag values are employed in the following computation.
C     For unequally-weighted reference as in GAMESS-US, some more
C     operations are required, but the CP-MCSCF part has to be
C     modified, so this may not be realized easily.
C
      nConfL=Max(nconf1,nint(xispsm(1,1)))
      nConfR=Max(nconf1,nint(xispsm(1,1)))
      Call GetMem('CIL','ALLO','REAL',ipL,nConfL)
      Call GetMem('CIR','ALLO','REAL',ipR,nConfR)
      !! iR = iRLXRoot
      Do jR = 1, nRoots
        Call CSF2SD(Work(ipIn(ipCI)+(jR-1)*nconf1),Work(ipL),1)
        Do kR = 1, jR !! jR-1
          iSLag = jR + nRoots*(kR-1)
          vSLag = SLag(iSLag)
          If (abs(vSLag).le.1.0d-10) Cycle
C
          Call CSF2SD(Work(ipIn(ipCI)+(jR-1)*nconf1),Work(ipL),1)
          iRC=opout(ipCI)
          Call CSF2SD(Work(ipIn(ipCI)+(kR-1)*nconf1),Work(ipR),1)
          iRC=opout(ipCI)
          iRC=ipnout(-1)
          icsm=1
          issm=1
          Call Densi2(2,Work(ipG1r),Work(ipG2r),
     &                  Work(ipL),Work(ipR),0,0,0,n1dens,n2dens)
          !! For RDM1
          ij=0
          Do i=0,ntAsh-1
            Do j=0,i-1
              Work(ipG1q+ij)=Work(ipG1q+ij)+
     *          (Work(ipG1r+i*ntAsh+j)+
     *           Work(ipG1r+j*ntAsh+i))*Half*vSLag
              ij=ij+1
            End Do
            Work(ipG1q+ij)=Work(ipG1q+ij)
     *        + Work(ipG1r+i*ntAsh+i)*vSLag
            ij=ij+1
          End Do
          !! For RDM2
          Do i=1,ntAsh**2
            j=itri(i,i)-1
            Work(ipG2r+j)=Half*Work(ipG2r+j)
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
                    Work(ipG2q+ijkl)=Work(ipG2q+ijkl)
     *                + factor*Work(ipG2r+ij2*(ij2+1)/2+kl2)
                    ij2=Max(j*ntAsh+i,l*ntAsh+k)
                    kl2=Min(j*ntAsh+i,l*ntAsh+k)
                    Work(ipG2q+ijkl)=Work(ipG2q+ijkl)
     &                + factor*Work(ipG2r+ij2*(ij2+1)/2+kl2)
                    If (k.ne.l) Then
                      ij2=i*ntAsh+j
                      kl2=l*ntAsh+k
                      Work(ipG2q+ijkl)=Work(ipG2q+ijkl)
     &                  + factor*Work(ipG2r+ij2*(ij2+1)/2+kl2)
                      If (ij.ne.kl) Then
                        ij2=Max(j*ntAsh+i,k*ntAsh+l)
                        kl2=Min(j*ntAsh+i,k*ntAsh+l)
                        Work(ipG2q+ijkl)=Work(ipG2q+ijkl)
     &                    + factor*Work(ipG2r+ij2*(ij2+1)/2+kl2)
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
                  Work(ipG2q+ijkl)=Work(ipG2q+ijkl)
     *              + factor*Work(ipG2r+ij2*(ij2+1)/2+kl2)
                  If (k.ne.l) Then
                    kl2=l*ntAsh+k
                    Work(ipG2q+ijkl)=Work(ipG2q+ijkl)
     &                + factor*Work(ipG2r+ij2*(ij2+1)/2+kl2)
                  End If
                End If
              End Do
            End Do
          End Do
        End Do
      End Do
      Call GetMem('CIL','FREE','REAL',ipL,nConfL)
      Call GetMem('CIR','FREE','REAL',ipR,nConfR)
      nConf=ncsf(1)
C
      Return
C
      End Subroutine PT2_SLag
C
       End
