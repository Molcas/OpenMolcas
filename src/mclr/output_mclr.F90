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

!#define _DEBUGPRINT_
subroutine OutPut_MCLR(iKapDisp,isigdisp,iCiDisp,iCiSigDisp,iRHSDisp,iRHSCIDisp,converged)
!***********************************************************************
!                                                                      *
! Contracts the response coefficient to the hessian                    *
!                                                                      *
! Input                                                                *
!       iKapDisp : Disk locations of solutions to respons equation     *
!       iSigDisp : Disk locations of RHS                               *
!       iCIDisp  : Disk locations of CI Soulutions to response         *
!       iCISigDisp : Disk locations of RHS                             *
!       nHess    : Length of hessian                                   *
!                                                                      *
! Output to disk                                                       *
!                                                                      *
!       RespHess : Hessian etc                                         *
!       Hess     : Hessian etc                                         *
!                                                                      *
! Author: Anders Bernhardsson, 1996                                    *
!         Theoretical Chemistry, University of Lund                    *
!***********************************************************************

use MckDat, only: sLength
use Arrays, only: Hss
use ipPage, only: W
use stdalloc, only: mma_allocate, mma_deallocate
use MCLR_Data, only: nConf1, nDensC
use MCLR_Data, only: nHess, lDisp
use MCLR_Data, only: LuTEMP
use MCLR_Data, only: XISPSM
use input_mclr, only: nDisp, Debug, nSym, State_Sym, iMethod, McKinley, Coor, lCalc, nCSF, nTPert

implicit none
integer iKapDisp(nDisp), isigdisp(nDisp), iCiDisp(nDisp), iCiSigDisp(nDisp), iRHSDisp(nDisp), iRHSCiDisp(nDisp)
logical converged(8)
character(len=8) Label
#ifdef _DEBUGPRINT_
character(len=20) Label2
logical elec_On
integer ip
#endif
integer Pstate_sym, ldisp2(8), ielec(3)
logical CI
real*8 Pola(6)
real*8, allocatable :: RHss(:)
real*8, allocatable :: Kap1(:), Kap2(:), sKap(:), rKap1(:), rKap2(:)
real*8, allocatable :: Hess(:), Hess2(:), Temp(:), ELEC(:), EG(:), ELOUT(:)
integer, allocatable :: NrDisp(:), DegDisp(:)
integer nHss, mSym, kSym, iDum, iDisp, iSym, nConfm, ipCIP1, ipCIP2, ipSP, ipRP1, ipRP2, jDisp, jSpin, iDisk, Len, i, iLen, iDis, &
        iRC, kDisp, kSpin, MaxI, MinI, Index, iOpt, Lu_10
real*8 rTempC1, rTempK1, Fact, rTempK2, rTempK3, rTempC2, rTempC3
real*8, external :: DDot_
integer, external :: ipGet, ipIn, ipClose
integer, external :: IsFreeUnit

!                                                                      *
!***********************************************************************
!                                                                      *

#ifdef _DEBUGPRINT_
debug = .true.
#else
debug = .false.
#endif
nHss = size(Hss)
nhess = nDisp*(nDisp+1)/2
call mma_allocate(RHss,nHss,Label='RHss')
RHss(:) = 0.0d0

!----------------------------------------------------------------------*
!
! Ok construct hessian
!
!----------------------------------------------------------------------*

mSym = 0
kSym = 0
idum = 1
idisp = 0
do iSym=1,nSym

  ! Calculate length of the density, Fock and Kappa matrix etc
  ! notice that this matrixes not necessary are symmetric.
  ! Store pointers.
  !
  ! Input:
  ! iSym : Symmetry of perturbation
  !
  ! Output: Commonblocks (Pointers.fh)

  call Setup_MCLR(iSym)
  PState_SYM = ieor(State_Sym-1,iSym-1)+1
  nconfM = max(ncsf(PState_Sym),nint(xispsm(Pstate_Sym,1)))
  nconf1 = ncsf(PState_Sym)
  CI = .false.
  if ((iMethod == 2) .and. (nconf1 > 0)) CI = .true.
  if (CI .and. (nconf1 == 1) .and. (isym == 1)) CI = .false.

  ! Allocate areas for scratch and state variables

  call mma_allocate(Kap1,nDensC,Label='Kap1')
  call mma_allocate(Kap2,nDensC,Label='Kap2')
  call mma_allocate(sKap,nDensC,Label='sKap')
  call mma_allocate(rkap1,nDensC,Label='rKap1')
  call mma_allocate(rkap2,nDensC,Label='rKap2')

  if (CI) then
    ipcip1 = ipget(nconfM)
    ipcip2 = ipget(nconfM)
    ipsp = ipget(nconfM)
    iprp2 = ipget(nconfM)
    iprp1 = ipget(nconfM)
  end if

  !                            [2]
  ! Calculate the diagonal of E    and store in core/disc

  do jDisp=1,lDisp(iSym)
    iDisp = iDisp+1
    jspin = 0
    if (iand(nTPert(idisp),1) == 1) jSpin = 1
    if (jspin == 0) then
      nconf1 = ncsf(Pstate_sym)
    else
      nconf1 = nint(xispsm(Pstate_Sym,1))
    end if

    !------------------------------------------------------------------*
    !                                                                  *
    !    Read response from disk                                       *
    !                                                                  *
    !------------------------------------------------------------------*

    iDisk = iKapDisp(iDisp)

    ! If disk address =/= -1 arrays on file.

    if (iDisk /= -1) then
      Len = nDensC
      call dDaFile(LuTemp,2,Kap1,Len,iDisk)
      iDisk = iSigDisp(iDisp)
      call dDaFile(LuTemp,2,SKap,Len,iDisk)
      iDisk = iRHSDisp(iDisp)
      call dDaFile(LuTemp,2,rKap1,Len,iDisk)
      do i=1,ndensC
        SKap(i) = -SKap(i)-rKap1(i)
      end do

      !call Recprt('ORB-RHS',' ',rKap1,nDensC,1)
      !write(6,*) 'ddot orb-resp',ddot_(ndensC,Kap1,1,Kap1,1)
      !write(6,*) 'ddot orb-sigma',ddot_(ndensC,SKap,1,SKap,1)
      !write(6,*) 'ddot orb-rhs',ddot_(ndensC,rKap1,1,rKap1,1)

      call GADSum(Kap1,Len)
      call GADSum(SKap,Len)
      call GADSum(rKap1,Len)

      if (CI) then
        ilen = nconf1
        idis = iCIDisp(iDisp)
        irc = ipin(ipCIp1)
        call dDaFile(LuTemp,2,W(ipCIp1)%Vec,iLen,iDis)
        idis = iCISigDisp(idisp)
        irc = ipin(ipSp)
        call dDaFile(LuTemp,2,W(ipSp)%Vec,iLen,iDis)
        idis = iRHSCIDisp(idisp)
        irc = ipin(iprp1)
        call dDaFile(LuTemp,2,W(iprp1)%Vec,iLen,iDis)
        irc = ipin(ipSp)
        irc = ipin(iprp1)
        do i=1,nConf1
          W(ipSp)%Vec(i) = -W(ipSp)%Vec(i)-W(iprp1)%Vec(i)
        end do

        !write(6,*) 'ddot ci-resp',ddot_(nConf1,W(ipcip1)%Vec,1,W(ipcip1)%Vec,1)
        !write(6,*) 'ddot ci-sigma',ddot_(nConf1,W(ipSp)%Vec,1,W(ipSp)%Vec,1)
        !write(6,*) 'ddot ci-rhs',ddot_(nConf1,W(iprp1)%Vec,1,W(iprp1)%Vec,1)

        call GADSum(W(ipCIp1)%Vec,iLen)
        call GADSum(W(ipSp)%Vec,iLen)
        call GADSum(W(ipRp1)%Vec,iLen)
      end if

    else

      Len = nDensC
      Kap1(1:Len) = 0.0d0
      call GADSum(Kap1,Len)
      sKap(1:Len) = 0.0d0
      call GADSum(SKap,Len)
      rKap1(1:Len) = 0.0d0
      call GADSum(rKap1,Len)

      if (CI) then

        ilen = nconf1
        irc = ipin(ipCIp1)
        call FZero(W(ipCIp1)%Vec,iLen)
        call GADSum(W(ipCIp1)%Vec,iLen)
        irc = ipin(ipSp)
        call FZero(W(ipSp)%Vec,iLen)
        call GADSum(W(ipSp)%Vec,iLen)
        irc = ipin(ipRp1)
        call FZero(W(ipRp1)%Vec,iLen)
        call GADSum(W(ipRp1)%Vec,iLen)

      end if

    end if

    !*******************************************************************

    do kDisp=1,jdisp

      !         (x) (2) (y)   (x) (y)    (y) (x)
      ! E    = k   E   k   + F   k   +  F   k
      !  Resp

      kspin = 0
      if (iand(nTPert(kdisp+ksym),1) == 1) kSpin = 1
      if (kspin == 0) then
        nconf1 = ncsf(PState_Sym)
      else
        nConf1 = nint(xispsm(Pstate_Sym,1))
      end if
      if (.not. lCalc(kDisp+ksym)) goto 120

      !write(6,*) 'kDisp+kSym',kDisp+kSym
      !write(6,*) 'iKapDisp(kdisp+ksym)',iKapDisp(kdisp+ksym)

      iDisk = iKapDisp(kDisp+kSym)
      if (iDisk /= -1) then
        Len = nDensC
        call dDaFile(LuTemp,2,Kap2,Len,iDisk)
        iDisk = iRHSDisp(kDisp+kSym)
        call dDaFile(LuTemp,2,rKap2,Len,iDisk)

        call GASync()
        call GADSum(Kap2,Len)
        call GADSum(rKap2,Len)

        if (CI) then
          ilen = nconf1
          idis = iCIDisp(kDisp+ksym)
          irc = ipin(ipCIp2)
          call dDaFile(LuTemp,2,W(ipCIp2)%Vec,iLen,iDis)
          idis = iRHSCIDisp(kdisp+ksym)
          irc = ipin(iprp2)
          call dDaFile(LuTemp,2,W(iprp2)%Vec,iLen,iDis)
          irc = ipin(ipsp)
          rTempc1 = DDot_(nConf1,W(ipCIp2)%Vec,1,W(ipsp)%Vec,1)

          call GASync() ! <----------------- NOTE!
          call GADSum(W(ipCIp2)%Vec,iLen)
          call GADSum(W(iprp2)%Vec,iLen)

        else
          rtempc1 = 0.0d0
        end if

      else

        call GASync()
        Len = nDensC
        call FZero(Kap2,Len)
        call GADSum(Kap2,Len)
        call FZero(rKap2,Len)
        call GADSum(rKap2,Len)
        if (CI) then
          ilen = nconf1
          call GASync()   ! <----------------- NOTE!
          irc = ipin(ipCIp2)
          call FZero(W(ipCIp2)%Vec,iLen)
          call GADSum(W(ipCIp2)%Vec,iLen)
          irc = ipin(iprp2)
          call FZero(W(iprp2)%Vec,iLen)
          call GADSum(W(iprp2)%Vec,iLen)
          irc = ipin(ipsp)
          rTempc1 = DDot_(nConf1,W(ipCIp2)%Vec,1,W(ipsp)%Vec,1)
        else
          rtempc1 = 0.0d0
        end if

      end if

      rTempk1 = DDot_(nDensC,Kap2,1,SKap,1)

      Fact = 1.0d0
      if (kdisp == jdisp) Fact = 2.0d0
      rTempk2 = Fact*DDot_(nDensC,Kap1,1,rKap2,1)
      if (kdisp /= jdisp) then
        rtempk3 = 1.0d0*DDot_(nDensC,rKap1,1,Kap2,1)
      else
        rTempk3 = 0.0d0
      end if
      if (CI) then
        Fact = 1.0d0
        if (kdisp == jdisp) Fact = 2.0d0
        irc = ipin(ipCip1)
        irc = ipin(iprp2)
        rTempc2 = Fact*DDot_(nConf1,W(ipCip1)%Vec,1,W(iprp2)%Vec,1)
        if (kdisp /= jdisp) then
          irc = ipin(iprp1)
          irc = ipin(ipCIp2)
          rtempc3 = DDot_(nConf1,W(iprp1)%Vec,1,W(ipCIp2)%Vec,1)
        else
          rTempc3 = 0.0d0
        end if
      else
        rtempc2 = 0.0d0
        rtempc3 = 0.0d0
      end if

      !write(6,*) kdisp,jdisp
      !write(6,*) rTempk1,rtempk2,rtempk3
      !write(6,*) rtempc1,rtempc2,rtempc3

      Maxi = max(kDisp,jDisp)
      Mini = min(kDisp,jDisp)
      index = mSym+Maxi*(Maxi-1)/2+Mini

      Rhss(Index) = Rhss(Index)+rTempk1+rtempk2+rtempk3+rtempc1+rtempc2+rtempc3

120   continue
    end do

    !*******************************************************************

  end do
  kSym = kSym+lDisp(iSym)
  mSym = mSym+lDisp(iSym)*(lDisp(iSym)+1)/2

  ! Free areas for scratch and state variables

  call mma_deallocate(rKap2)
  call mma_deallocate(rKap1)
  call mma_deallocate(sKap)
  call mma_deallocate(Kap2)
  call mma_deallocate(Kap1)
  if (CI) irc = ipclose(ipcip1)
end do

call mma_allocate(Hess,nHss,Label='Hess')
call mma_allocate(Hess2,nHss,Label='Hess2')
call mma_allocate(Temp,nHss,Label='Temp')
Temp(:) = 0.0d0
call mma_allocate(ELEC,3*ndisp,Label='ELEC')
call mma_allocate(EG,3*ndisp,Label='EG')
call mma_allocate(ELOUT,3*ndisp,Label='ELOUT')
irc = ipclose(-1)

!----------------------------------------------------------------------*
!
!     OK now when we have out Hessian, what should we do with it!
!
!----------------------------------------------------------------------*

! If a basis set is dependent on perturbation add terms
! constructed in mckinley.

call dcopy_(6,[0.0d0],0,pola,1)
idum = 1
iopt = ibset(0,sLength)
irc = 3*ndisp
Label = 'DOTELGR'
call drdMCk(irc,iopt,LaBeL,idum,EG,idum)
call GADsum(Hss,nHss)
call dcopy_(nHss,Hss,1,Hess2,1)
#ifdef _DEBUGPRINT_
elec_On = .true.
if (irc /= 0) elec_On = .false.
if (debug) then
  ip = 1
  do iSym=1,nSym
    write(label2,'(A,I2)') 'CHessian symmetry',iSym
    if (lDisp(iSym) /= 0) call TriPrt(label2,' ',Hess2(ip),lDisp(iSym))
    ip = ip+ldisp(isym)*(1+ldisp(isym))/2
  end do
end if

if (debug) then
  call MMSORT2(HESS2,ELEC,pola,ielec)
  call Recprt('CONN',' ',Elec,3*nDisp,1)
end if

!write(6,*) 'I am here 1'
call Recprt('Rhss','(5G20.10)',RHss,nhss,1)
call Recprt('Hss','(5G20.10)',Hss,nHss,1)
#endif

call DaXpY_(mSym,1.0d0,RHss,1,Hess2,1)

#ifdef _DEBUGPRINT_
if (debug) then
  call MMSORT2(RHSS,ELEC,pola,ielec)
  call Recprt('RESP',' ',Elec,3*nDisp,1)
end if
#endif

call MMSORT2(HESS2,ELEC,pola,ielec)

#ifdef _DEBUGPRINT_
if (debug) then
  call Recprt('R+C',' ',Elec,3*nDisp,1)
  ip = 1
  do iSym=1,nSym
    write(label2,'(A,I2)') 'Hessian symmetry',iSym
    if (lDisp(iSym) /= 0) call TriPrt(label2,' ',Hess2(ip),lDisp(iSym))
    ip = ip+ldisp(isym)*(1+ldisp(isym))/2
  end do
end if
#endif

call mmSort(Hess2,Hess,ldisp2)

if (McKinley) then

  iRC = -1
  iOpt = 0
  Label = 'StatHess'
  call dRdMck(iRC,iOpt,Label,idum,Temp,idum)
  if (iRC /= 0) then
    write(6,*) 'OutPut: Error reading MCKINT'
    write(6,'(A,A)') 'Label=',Label
    call Abend()
  end if

# ifdef _DEBUGPRINT_
  if (debug) then
    ip = 1
    do iSym=1,nSym
      write(label2,'(a,i2)') 'SHessian symmetry',iSym
      if (lDisp2(iSym) /= 0) call TriPrt(label2,' ',Temp(ip),lDisp2(iSym))
      ip = ip+ldisp2(isym)*(1+ldisp2(isym))/2
    end do
  end if
# endif
  call DaXpY_(mSym,1.0d0,Temp,1,Hess,1)
end if
#ifdef _DEBUGPRINT_
if (debug) then
  ip = 1
  do iSym=1,nSym
    write(label2,'(a,i2)') 'Hessian symmetry',iSym
    if (lDisp2(iSym) /= 0) call TriPrt(label2,' ',Hess(ip),lDisp2(iSym))
    ip = ip+ldisp2(isym)*(1+ldisp2(isym))/2
  end do
end if
#endif

if (McKinley) then
  iRC = -1
  iOpt = 0
  Label = 'Hess'
  call dWrMck(iRC,iOpt,Label,iDum,Hess,iDum)
  if (iRC /= 0) then
    write(6,*) 'OutPut: Error writing to MCKINT'
    write(6,'(A,A)') 'Label=',Label
    call Abend()
  end if
  call Put_iScalar('No of Internal coordinates',ldisp2(1))
  call Put_AnalHess(Hess,ldisp2(1)*(ldisp2(1)+1)/2)
end if

if (.true.) then
  iRC = -1
  iOpt = 0
  call mma_allocate(NrDisp,ndisp,Label='NrDisp')
  Label = 'NRCTDISP'
  call RdMck(irc,iopt,Label,idum,NrDisp,idum)
  iRC = -1
  iOpt = 0
  call mma_allocate(DegDisp,ndisp,Label='DegDisp')
  Label = 'DegDisp'
  call RdMck(irc,iopt,Label,idum,DegDisp,idum)
  if (iRC /= 0) then
    write(6,*) 'OutPut: Error reading RELAX'
    write(6,'(A,A)') 'Label=',Label
    call Abend()
  end if

  !if (debug) call HssPrt_MCLR(DegDisp,Hess,ldisp2)
  !call Recprt('hess',' ',Hess,nhss,1)

  call daxpy_(3*ndisp,-1.0d0,EG,1,ELEC,1)
# ifdef _DEBUGPRINT_
  if (debug .and. elec_On) call Recprt('ELEC-ST',' ',EG,3*nDisp,1)
  if (debug .and. elec_On) call Recprt('ELEC-TOT',' ',Elec,3*nDisp,1)
# endif

  Lu_10 = 10
  Lu_10 = IsFreeUnit(Lu_10)
  call molcas_open(lu_10,'UNSYM')
  !open(unit=Lu_10,file='UNSYM')

  if (Mckinley) then
    call FreqAnal(DegDisp,NrDisp,Hess,converged,ELEC,ielec,elout,ldisp2,Lu_10)
    call Niclas(Hess,coor,Lu_10)
  end if
  write(6,*)
  write(6,*)
  write(6,*) '************************************'
  write(6,*) '*                                  *'
  write(6,*) '*       Polarizabilities           *'
  write(6,*) '*                                  *'
  write(6,*) '************************************'
  write(6,*)
  write(6,*)
  call Add_Info('POLARIZABILITIES',Pola,6,2)

  ! Go from energy derivative to polarizability, there is a difference
  ! in the sign in the definition.

  call DScal_(6,-1.0d0,Pola,1)

  call TriPrt(' ',' ',Pola,3)
  close(Lu_10)
  call mma_deallocate(NrDisp)
  call mma_deallocate(DegDisp)
end if
!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*

call mma_deallocate(ELOUT)
call mma_deallocate(EG)
call mma_deallocate(ELEC)
call mma_deallocate(Temp)
call mma_deallocate(Hess2)
call mma_deallocate(Hess)
call mma_deallocate(RHss)

end subroutine OutPut_MCLR
