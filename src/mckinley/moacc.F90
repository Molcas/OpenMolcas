!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1996, Anders Bernhardsson                              *
!***********************************************************************

subroutine MOAcc(AOInt,nint,Temp1,Temp2,Temp3,nTemp,MOInt,nMO,ishell,Ci,nCi,Cj,nCj,Ck,nCk,Cl,nCl,moip,nACO,pert,nOp,iBasa,iCmpa, &
                 icar,icnt,indgrd,D,fact,iao,iaost,Buffer,Tempi,nij,nkl,nbasi,nbasj,icmp,jcmp)
!***********************************************************************
!                                                                      *
!     Transforms a batch of unsymmetrized integrals to                 *
!     active integral batches and FM                                   *
!     All MO combinations are constructed                              *
!     They will be needed with unsymmetric perurbations                *
!                                                                      *
!     Author: Anders Bernhardsson, Dept. of Theoretical Chemistry,     *
!             University of Lund, Sweden. Januar '96                   *
!***********************************************************************

use Symmetry_Info, only: nIrrep, iChTbl, iOper
use Gateway_Info, only: CutInt

implicit real*8(A-H,O-Z)
#include "real.fh"
#include "etwas.fh"
!#include "print.fh"
real*8 AOInt(nkl,nij), MOint(nMO), Temp1(nTemp), Temp2(naco,naco), Ck(nCk), Cl(nCl), D(*), &
       Buffer(nbasi,icmp,nbasj,jcmp,0:nirrep-1, nAco*(naco+1)/2,*)
integer moip(0:7), nOp(4), ishell(4), iao(4), iAOST(4), ibasa(4), icmpa(4), indgrd(3,4,0:7)
logical pert(0:7)
real*8 Prmt(0:7)
data Prmt/1.d0,-1.d0,-1.d0,1.d0,-1.d0,1.d0,1.d0,-1.d0/
! Statement Function
xPrmt(i,j) = Prmt(iand(i,j))

iCB = 2**(icar-1)
rFact = xPrmt(ioper(nOp(icnt)),icb)*fact

iBas = iBasa(1)
jBas = iBasa(2)
kBas = iBasa(3)
lBas = iBasa(4)

kCmp = iCmpa(3)
lCmp = iCmpa(4)
kk = 0
do kIrrep=0,nIrrep-1
  sfact = dble(ichtbl(kirrep,nop(3)))
  do kAsh=1,nAsh(kIrrep)
    do k=1,kcmp*kbas
      kk = kk+1
      Ck(kk) = Ck(kk)*sFact
    end do
  end do
end do
kk = 0
do kIrrep=0,nIrrep-1
  sfact = dble(ichtbl(kirrep,nop(4)))
  do kAsh=1,nAsh(kIrrep)
    do k=1,lcmp*lbas
      kk = kk+1
      Cl(kk) = Cl(kk)*sFact
    end do
  end do
end do
rk = DNrm2_(nck,ck,1)
rl = DNrm2_(ncl,cl,1)
ij = 0
nt = lBas*lCmp*kcmp*kbas
do jc=1,jcmp
  do jb=1+iaost(2),iaost(2)+jbas
    do ic=1,icmp
      do ib=1+iaost(1),iaost(1)+ibas

        ij = ij+1
        vij = DNrm2_(nt,AOInt(1,ij),1)
        if (abs(vij*rk*rl) < cutint) goto 1000
        ipC = 0
        do kAsh=1,nAco
          ipM = (kAsh-1)*lbas*lcmp
          il = 0
          do i=1,lbas*lcmp
            Temp1(i+ipM) = 0.0d0
            do k=1,kCmp*kBas
              il = il+1
              Temp1(i+ipM) = Temp1(ipm+i)+Ck(k+ipc)*AOINT(il,ij)
            end do
          end do
          ipC = ipC+kBas*kCmp
        end do
        ipC = 0
        do lAsh=1,naco
          il = 0
          do kash=1,naco
            Temp2(kash,lash) = 0.0d0
            do l=1,lbas*lcmp
              il = il+1
              Temp2(kash,lash) = Temp2(kash,lash)+Cl(ipc+l)*Temp1(il)
            end do
          end do
          ipC = ipC+lBas*lCmp
        end do

        if (iShell(3) /= ishell(4)) then

          do iSPert=0,nIrrep-1
            if (pert(isPert)) then
              rFact2 = rFact*dble(iChtbl(ispert,nop(icnt)))
              k = abs(indgrd(icar,icnt,ispert))
              j = 0
              do lIrr=0,nIrrep-1
                do lAsh=1,nAsh(lIrr)
                  do kIrr=0,lIrr
                    irest = ieor(ieor(ioper(ispert),ioper(kirr)),ioper(lirr))
                    kMax = nAsh(kIrr)
                    if (kIrr == lIrr) kMax = lAsh
                    do kAsh=1,kMax
                      kk = kash+moip(kirr)
                      ll = lash+moip(lirr)
                      j = j+1
                      do jIrr=0,nIrrep-1
                        rPj = dble(iChTbl(jIrr,nop(2)))
                        iirr = nropr(ieor(iOPER(jirr),irest))
                        rPij = rPj*dble(iChTbl(iIrr,nop(1)))*rfact2
                        buffer(ib,ic,jb,jc,iirr,j,k) = buffer(ib,ic,jb,jc,iirr,j,k)+rpij*Temp2(kk,ll)+rpij*Temp2(ll,kk)
                      end do
                    end do
                  end do
                end do
              end do
            end if
          end do
        else

          do iSPert=0,nIrrep-1
            if (pert(isPert)) then
              rFact2 = rFact*dble(iChtbl(ispert,nop(icnt)))
              k = abs(indgrd(icar,icnt,ispert))
              j = 0
              do lIrr=0,nIrrep-1
                do lAsh=1,nAsh(lIrr)
                  do kIrr=0,lIrr
                    irest = ieor(ieor(ioper(ispert),ioper(kirr)),ioper(lirr))
                    kMax = nAsh(kIrr)
                    if (kIrr == lIrr) kMax = lAsh
                    do kAsh=1,kMax
                      kk = kash+moip(kirr)
                      ll = lash+moip(lirr)
                      j = j+1
                      do jIrr=0,nIrrep-1
                        rPj = dble(iChTbl(jIrr,nop(2)))
                        iirr = nropr(ieor(iOPER(jirr),irest))
                        rPij = rPj*dble(iChTbl(iIrr,nop(1)))*rfact2
                        buffer(ib,ic,jb,jc,iirr,j,k) = buffer(ib,ic,jb,jc,iirr,j,k)+rpij*Temp2(kk,ll)
                      end do
                    end do
                  end do
                end do
              end do
            end if
          end do
        end if
1000    continue
      end do
    end do
  end do
end do

kk = 0
do kIrrep=0,nIrrep-1
  sfact = dble(ichtbl(kirrep,nop(3)))
  do kAsh=1,nAsh(kIrrep)
    do k=1,kcmp*kbas
      kk = kk+1
      Ck(kk) = Ck(kk)*sFact
    end do
  end do
end do
kk = 0
do kIrrep=0,nIrrep-1
  sfact = dble(ichtbl(kirrep,nop(4)))
  do kAsh=1,nAsh(kIrrep)
    do k=1,lcmp*lbas
      kk = kk+1
      Cl(kk) = Cl(kk)*sFact
    end do
  end do
end do

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_integer(nint)
  call Unused_real(Temp3)
  call Unused_real_array(MOInt)
  call Unused_real(Ci)
  call Unused_integer(nCi)
  call Unused_real(Cj)
  call Unused_integer(nCj)
  call Unused_real_array(D)
  call Unused_integer_array(iao)
  call Unused_real(Tempi)
end if

end subroutine MOAcc
