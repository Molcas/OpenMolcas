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
! Copyright (C) 1995, Anders Bernhardsson                              *
!***********************************************************************

subroutine ClrBuf(idcrr,idcrs,idcrt,ngr,Shijij,iAnga,iCmp,iCmpa,iShll,iShell,jShell,iBasi,jBasj,kBask,lBasl,Dij1,Dij2,mDij,nDij, &
                  Dkl1,Dkl2,mDkl,nDkl,Dik1,Dik2,mDik,nDik,Dil1,Dil2,mDil,nDil,Djk1,Djk2,mDjk,nDjk,Djl1,Djl2,mDjl,nDjl,rFinal, &
                  nFinal,FckTmp,nFT,Scrtch1,nS1,Scrtch2,nS2,Temp,nTemp,TwoHam,nTwo,IndGrd,Indx,iAO,iAOst,iuvwx,n8,ltri,moip,nAcO, &
                  rmoin,nmoin,ntemptot,Buffer,nop,din,dan,new_fock)
!***********************************************************************
!                                                                      *
!       Called from: Twoel                                             *
!       takes care of the integrals                                    *
!       integrals -> fckmatrix,MO                                      *
!       in the near feature a disk based version                       *
!                                                                      *
!       Calling:   CntrDens : Gets the indexes for d1                  *
!                  MkFck : Add up the integrals on the Fock Matrix     *
!                                                                      *
!       Author: Anders Bernhardsson, Theoretical Chemistry,            *
!               University of Lund, Sweden, June '95                   *
!***********************************************************************

use McKinley_global, only: CPUStat, ipDisp, ipDisp2, nFckAcc, nMethod, nMOTrans, RASSCF
use pso_stuff, only: ndens
use Symmetry_Info, only: nIrrep
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: idcrr, idcrs, idcrt, ngr, iAnga(4), iCmp(4), icmpa(4), iShll(4), iShell(4), jShell(4), iBasi, &
                                 jBasj, kBask, lBasl, mDij, nDij, mDkl, nDkl, mDik, nDik, mDil, nDil, mDjk, nDjk, mDjl, nDjl, &
                                 nFinal, nFT, nS1, nS2, nTemp, nTwo, IndGrd(3,4,0:7), Indx(3,4), iAO(4), iAOst(4), iuvwx(4), &
                                 moip(0:7), nAcO, nmoin, ntemptot, nop(4)
logical(kind=iwp), intent(in) :: Shijij, n8, ltri, new_fock
real(kind=wp), intent(in) :: Dij1(mDij,nDij), Dij2(mDij,nDij), Dkl1(mDkl,nDkl), Dkl2(mDkl,nDkl), Dik1(mDik,nDik), Dik2(mDik,nDik), &
                             Dil1(mDil,nDil), Dil2(mDil,nDil), Djk1(mDjk,nDjk), Djk2(mDjk,nDjk), Djl1(mDjl,nDjl), Djl2(mDjl,nDjl), &
                             rFinal(nFinal), din(*), dan(*)
real(kind=wp), intent(out) :: FckTmp(nFT), Scrtch1(nS1), Scrtch2(nS2), Temp(nTemp)
real(kind=wp), intent(inout) :: TwoHam(nTwo), rmoin(nMOIN), buffer(*)
integer(kind=iwp) :: iCar, iCent, iCnt, iGr, ii, iIrrep, ij1, ij2, ij3, ij4, ik1, ik2, ik3, ik4, il1, il2, il3, il4, ip, ipFin, &
                     jk1, jk2, jk3, jk4, jl1, jl2, jl3, jl4, jOp(6), kl1, kl2, kl3, kl4, nabcd, nao, nijkl
logical(kind=iwp) :: pert(0:7)
real(kind=wp) :: dum1, dum2, dum3, ExFac, Time

nijkl = iBasi*jBasj*kBask*lBasl
nabcd = iCmp(1)*iCmp(2)*iCmp(3)*iCmp(4)
call Timing(dum1,Time,dum2,dum3)

ExFac = One

if (ltri) then
  if (.not. new_fock) then
    !------------------------------------------------------------*
    !
    !   Get the size of the work area that should be contracted
    !
    !------------------------------------------------------------*
    if (jShell(1) >= jShell(2)) then
      ij1 = iBasi
      ij2 = jBasj
      ij3 = iCmp(1)
      ij4 = iCmp(2)
    else
      ij1 = jBasj
      ij2 = iBasi
      ij3 = iCmp(2)
      ij4 = iCmp(1)
    end if
    if (jShell(3) >= jShell(4)) then
      kl1 = kBask
      kl2 = lBasl
      kl3 = iCmp(3)
      kl4 = iCmp(4)
    else
      kl1 = lBasl
      kl2 = kBask
      kl3 = iCmp(4)
      kl4 = iCmp(3)
    end if
    if (jShell(1) >= jShell(3)) then
      ik1 = iBasi
      ik2 = kBask
      ik3 = iCmp(1)
      ik4 = iCmp(3)
    else
      ik1 = kBask
      ik2 = iBasi
      ik3 = iCmp(3)
      ik4 = iCmp(1)
    end if
    if (jShell(1) >= jShell(4)) then
      il1 = iBasi
      il2 = lBasl
      il3 = iCmp(1)
      il4 = iCmp(4)
    else
      il1 = lBasl
      il2 = iBasi
      il3 = iCmp(4)
      il4 = iCmp(1)
    end if
    if (jShell(2) >= jShell(3)) then
      jk1 = jBasj
      jk2 = kBask
      jk3 = iCmp(2)
      jk4 = iCmp(3)
    else
      jk1 = kBask
      jk2 = jBasj
      jk3 = iCmp(3)
      jk4 = iCmp(2)
    end if
    if (jShell(2) >= jShell(4)) then
      jl1 = jBasj
      jl2 = lBasl
      jl3 = iCmp(2)
      jl4 = iCmp(4)
    else
      jl1 = lBasl
      jl2 = jBasj
      jl3 = iCmp(4)
      jl4 = iCmp(2)
    end if

    ! Here we go

    nao = iBasi*jBasj*kBask*lBasl*iCmp(1)*iCmp(2)*iCmp(3)*iCmp(4)
    call CtlDns(iDCRR,iDCRS,iDCRT,jOp)
  end if

  ! Add this contribution to the (Inactive) Fock Matrix

  ! Out from  CD : jOp

  do iCent=1,4
    do iCar=1,3
      pert(:) = .false.

      ! To which irreps does this derivative contribute?

      do iIrrep=0,nIrrep-1
        if (indgrd(iCar,iCent,iIrrep) /= 0) pert(iIrrep) = .true.
      end do

      if (Indx(iCar,iCent) > 0) then
        iGr = Indx(iCar,iCent)-1
        ipFin = 1+iGr*nijkl*nabcd
        if (.not. new_fock) then
          call MkFck(iAnga,iCmp,Shijij,iShll,iShell,iBasi,jBasj,kBask,lBasl,iAO,iAOst,nop,jop,Dij1,mDij,nDij,ij1,ij2,ij3,ij4,Dkl1, &
                     mDkl,nDkl,kl1,kl2,kl3,kl4,Dik1,mDik,nDik,ik1,ik2,ik3,ik4,Dil1,mDil,nDil,il1,il2,il3,il4,Djk1,mDjk,nDjk,jk1, &
                     jk2,jk3,jk4,Djl1,mDjl,nDjl,jl1,jl2,jl3,jl4,rFinal(ipFin),nAO,TwoHam,nTwo,Scrtch2,nS2,FckTmp,nFT,pert, &
                     iuvwx(iCent),iCent,iCar,indgrd,ipDisp)
          if (nMethod == RASSCF) call MkFck(iAnga,iCmp,Shijij,iShll,iShell,iBasi,jBasj,kBask,lBasl,iAO,iAOst,nop,jop,Dij2,mDij, &
                                            nDij,ij1,ij2,ij3,ij4,Dkl2,mDkl,nDkl,kl1,kl2,kl3,kl4,Dik2,mDik,nDik,ik1,ik2,ik3,ik4, &
                                            Dil2,mDil,nDil,il1,il2,il3,il4,Djk2,mDjk,nDjk,jk1,jk2,jk3,jk4,Djl2,mDjl,nDjl,jl1,jl2, &
                                            jl3,jl4,rFinal(ipFin),nAO,TwoHam,nTwo,Scrtch2,nS2,FckTmp,nFT,pert,iuvwx(iCent),iCent, &
                                            iCar,indgrd,ipDisp2)

        else
          ip = ipDisp(abs(indgrd(iCar,iCent,0)))
          call FckAcc_NoSym(iCmpa(1),iCmpa(2),iCmpa(3),iCmpa(4),Shijij,iShell,nijkl,rFinal(ipFin),TwoHam(ip),dan,ndens,iAO,iAOst, &
                            iBasi,jBasj,kBask,lBasl,ExFac)
          if (nMethod == RASSCF) then
            ip = ipDisp2(abs(indgrd(iCar,iCent,0)))
            call FckAcc_NoSym(iCmpa(1),iCmpa(2),iCmpa(3),iCmpa(4),Shijij,iShell,nijkl,rFinal(ipFin),TwoHam(ip),din,nDens,iAO, &
                              iAOst,iBasi,jBasj,kBask,lBasl,ExFac)
          end if
        end if

      else if (Indx(iCar,iCent) < 0) then
        Temp(1:nabcd*nijkl) = Zero
        do iCnt=1,4
          iGr = Indx(iCar,iCnt)
          if (iGr > 0) then
            ipFin = 1+(iGr-1)*nijkl*nabcd
            do ii=1,nabcd*nijkl
              Temp(ii) = Temp(ii)-rFinal(ipFin-1+ii)
            end do
          end if
        end do

        if (.not. new_fock) then
          call MkFck(iAnga,iCmp,Shijij,iShll,iShell,iBasi,jBasj,kBask,lBasl,iAO,iAOst,nop,jop,Dij1,mDij,nDij,ij1,ij2,ij3,ij4,Dkl1, &
                     mDkl,nDkl,kl1,kl2,kl3,kl4,Dik1,mDik,nDik,ik1,ik2,ik3,ik4,Dil1,mDil,nDil,il1,il2,il3,il4,Djk1,mDjk,nDjk,jk1, &
                     jk2,jk3,jk4,Djl1,mDjl,nDjl,jl1,jl2,jl3,jl4,Temp,nAO,TwoHam,nTwo,Scrtch2,nS2,FckTmp,nFT,pert,iuvwx(iCent), &
                     icent,iCar,indgrd,ipDisp)
          if (nMethod == RASSCF) call MkFck(iAnga,iCmp,Shijij,iShll,iShell,iBasi,jBasj,kBask,lBasl,iAO,iAOst,nop,jop,Dij2,mDij, &
                                            nDij,ij1,ij2,ij3,ij4,Dkl2,mDkl,nDkl,kl1,kl2,kl3,kl4,Dik2,mDik,nDik,ik1,ik2,ik3,ik4, &
                                            Dil2,mDil,nDil,il1,il2,il3,il4,Djk2,mDjk,nDjk,jk1,jk2,jk3,jk4,Djl2,mDjl,nDjl,jl1,jl2, &
                                            jl3,jl4,Temp,nAO,TwoHam,nTwo,Scrtch2,nS2,FckTmp,nFT,pert,iuvwx(iCent),icent,iCar, &
                                            indgrd,ipDisp2)

        else
          ip = ipDisp(abs(indgrd(iCar,iCent,0)))
          call FckAcc_NoSym(iCmpa(1),iCmpa(2),iCmpa(3),iCmpa(4),Shijij,iShell,nijkl,Temp,TwoHam(ip),dan,nDens,iAO,iAOst,iBasi, &
                            jBasj,kBask,lBasl,ExFac)
          if (nMethod == RASSCF) then
            ip = ipDisp2(abs(indgrd(iCar,iCent,0)))
            call FckAcc_NoSym(iCmpa(1),iCmpa(2),iCmpa(3),iCmpa(4),Shijij,iShell,nijkl,Temp,TwoHam(ip),din,nDens,iAO,iAOst,iBasi, &
                              jBasj,kBask,lBasl,ExFac)
          end if
        end if

      end if
    end do
  end do
  call Timing(dum1,Time,dum2,dum3)
  CPUStat(nFckAcc) = CPUStat(nFckAcc)+Time
end if

if (n8 .and. (nmethod == RASSCF)) call MakeMO(rFinal,Scrtch1,nTempTot,nFinal,iCmp,iCmpa,iBasi,jBasj,kBask,lBasl,nGr,Indx,moip, &
                                              naco,nop,indgrd,ishll,ishell,rmoin,nMOIN,iuvwx,iaost,Buffer,ianga)

call Timing(dum1,Time,dum2,dum3)
CPUStat(nMOTrans) = CPUStat(nMOTrans)+Time

return

end subroutine ClrBuf
