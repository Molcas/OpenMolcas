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
       Subroutine rhs_sa(Fock)
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

       Real*8 Fock(*)
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
       End
