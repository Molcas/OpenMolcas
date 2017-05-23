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
      Subroutine RHS_NAC(Fock)
      Implicit None
#include "Input.fh"
#include "Pointers.fh"
#include "WrkSpc.fh"
#include "real.fh"
#include "sa.fh"
#include "detdim.fh"
#include "cicisp_mclr.fh"
#include "cands.fh"
      Real*8 Fock(*)
      Integer ng1,ng2,i,j,k,l,ij,kl,ijkl,ij2,kl2,ijkl2
      Integer ipG1q,ipG2q,ipG1m,ipG1r,ipG2r,ipF,ipT,iTri
      Integer ipIn,opOut,ipnOut,nConfL,nConfR,ipL,ipR,iRC,LuDens
      Real*8 factor
      External ipIn,opOut,ipnOut
*                                                                      *
************************************************************************
*                                                                      *
      iTri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
*                                                                      *
************************************************************************
*                                                                      *
      ng1=ntAsh*(ntAsh+1)/2
      ng2=ng1*(ng1+1)/2
      Call Getmem('ONED','ALLO','REAL',ipG1q,ng1)
      Call Getmem('TWOD','ALLO','REAL',ipG2q,ng2)
      Call Getmem('ONED-','ALLO','REAL',ipG1m,ng1)
      Call Allocate_Work(ipG1r,n1dens)
      Call Allocate_Work(ipG2r,n2dens)
      Call dcopy_(n1dens,zero,0,work(ipG1r),1)
      Call dcopy_(n2dens,zero,0,work(ipG2r),1)
*
**    Calculate one- and two-particle transition matrices
**    from the CI vectors of the two NAC states
**    (code copied from CIdens_SA, no symmetry)
*
      nConfL=Max(nconf1,nint(xispsm(1,1)))
      nConfR=Max(nconf1,nint(xispsm(1,1)))
      Call GetMem('CIL','ALLO','REAL',ipL,nConfL)
      Call GetMem('CIR','ALLO','REAL',ipR,nConfR)
      Call CSF2SD(Work(ipIn(ipCI)+(NSSA(2)-1)*nconf1),Work(ipL),1)
      iRC=opout(ipCI)
      Call CSF2SD(Work(ipIn(ipCI)+(NSSA(1)-1)*nconf1),Work(ipR),1)
      iRC=opout(ipCI)
      iRC=ipnout(-1)
      icsm=1
      issm=1
      Call Densi2(2,Work(ipG1r),Work(ipG2r),
     &              Work(ipL),Work(ipR),0,0,0,n1dens,n2dens)
      Call GetMem('CIL','FREE','REAL',ipL,nConfL)
      Call GetMem('CIR','FREE','REAL',ipR,nConfR)
*
**    Symmetrize densities
**    For the one-particle density, save the antisymmetric part too
*
      ij=0
      Do i=0,ntAsh-1
        Do j=0,i-1
          Work(ipG1q+ij)=(Work(ipG1r+i*ntAsh+j)+
     &                    Work(ipG1r+j*ntAsh+i))*Half
          Work(ipG1m+ij)=(Work(ipG1r+i*ntAsh+j)-
     &                    Work(ipG1r+j*ntAsh+i))*Half
          ij=ij+1
        End Do
        Work(ipG1q+ij)=Work(ipG1r+i*ntAsh+i)
        Work(ipG1m+ij)=Zero
        ij=ij+1
      End Do
*
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
                factor=Quart
                If (ij.eq.kl) factor=Half
                ijkl=ij*(ij+1)/2+kl
                ij2=i*ntAsh+j
                kl2=k*ntAsh+l
                Work(ipG2q+ijkl)=factor*Work(ipG2r+ij2*(ij2+1)/2+kl2)
                ij2=Max(j*ntAsh+i,l*ntAsh+k)
                kl2=Min(j*ntAsh+i,l*ntAsh+k)
                Work(ipG2q+ijkl)=Work(ipG2q+ijkl)+
     &                           factor*Work(ipG2r+ij2*(ij2+1)/2+kl2)
                If (k.ne.l) Then
                  ij2=i*ntAsh+j
                  kl2=l*ntAsh+k
                  Work(ipG2q+ijkl)=Work(ipG2q+ijkl)+
     &                             factor*Work(ipG2r+ij2*(ij2+1)/2+kl2)
                  If (ij.ne.kl) Then
                    ij2=Max(j*ntAsh+i,k*ntAsh+l)
                    kl2=Min(j*ntAsh+i,k*ntAsh+l)
                    Work(ipG2q+ijkl)=Work(ipG2q+ijkl)+
     &                              factor*Work(ipG2r+ij2*(ij2+1)/2+kl2)
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
              Work(ipG2q+ijkl)=factor*Work(ipG2r+ij2*(ij2+1)/2+kl2)
              If (k.ne.l) Then
                kl2=l*ntAsh+k
                Work(ipG2q+ijkl)=Work(ipG2q+ijkl)+
     &                           factor*Work(ipG2r+ij2*(ij2+1)/2+kl2)
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
      Call dDaFile(LuDens,1,Work(ipG1m),ng1,iRC)
      Call DaClos(LuDens)
      Call Put_D1MO(Work(ipG1q),ng1)
      Call Put_P2MO(Work(ipG2q),ng2)
*
**    Store transition Fock matrix
*
      Call Allocate_Work(ipT,nDens2)
      Call Allocate_Work(ipF,nDens2)
      Do i=1,ntAsh
        Do j=1,ntAsh
          Work(ipG1r+ntAsh*(j-1)+i-1)=Work(ipG1q+iTri(i,j)-1)
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
              Work(ipG2r+ijkl2-1)=factor*Work(ipG2q+ijkl-1)
            End Do
          End Do
        End Do
      End Do
* Note: 1st arg = zero for no inactive density (TDM)
      Call FockGen(Zero,Work(ipG1r),Work(ipG2r),Work(ipT),Fock,1)
      Call TCMO(Work(ipT),1,-2)
      ij=0
      Do k=1,nSym
        Do i=0,nBas(k)-1
          Do j=0,i-1
            Work(ipF+ij)=Work(ipT+ipMat(k,k)-1+nBas(k)*j+i)+
     &                   Work(ipT+ipMat(k,k)-1+nBas(k)*i+j)
            ij=ij+1
          End Do
          Work(ipF+ij)=Work(ipT+ipMat(k,k)-1+nBas(k)*i+i)
          ij=ij+1
        End Do
      End Do
      Call Put_Fock_Occ(Work(ipF),nDens2)
*
      Call Free_Work(ipT)
      Call Free_Work(ipF)
      Call Free_Work(ipG1r)
      Call Free_Work(ipG2r)
      Call Getmem('ONED','FREE','REAL',ipG1q,ng1)
      Call Getmem('TWOD','FREE','REAL',ipG2q,ng2)
      Call Getmem('ONED-','FREE','REAL',ipG1m,ng1)
*
      Return
      End
