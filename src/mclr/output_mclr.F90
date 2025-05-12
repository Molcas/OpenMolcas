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
!       iKapDisp   : Disk locations of solutions to respons equation   *
!       iSigDisp   : Disk locations of RHS                             *
!       iCIDisp    : Disk locations of CI Soulutions to response       *
!       iCISigDisp : Disk locations of RHS                             *
!       nHess      : Length of hessian                                 *
!                                                                      *
! Output to disk                                                       *
!                                                                      *
!       RespHess : Hessian etc                                         *
!       Hess     : Hessian etc                                         *
!                                                                      *
! Author: Anders Bernhardsson, 1996                                    *
!         Theoretical Chemistry, University of Lund                    *
!***********************************************************************

use Index_Functions, only: iTri, nTri_Elem
use Symmetry_Info, only: Mul
use MckDat, only: sLength
use ipPage, only: ipclose, ipget, ipin, W
use MCLR_Data, only: Hss, lDisp, LuTEMP, nConf1, nDensC, nHess, XISPSM
use input_mclr, only: Coor, Debug, iMethod, lCalc, McKinley, nCSF, nDisp, nSym, nTPert, State_Sym
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iKapDisp(nDisp), isigdisp(nDisp), iCiDisp(nDisp), iCiSigDisp(nDisp), iRHSDisp(nDisp), &
                                 iRHSCiDisp(nDisp)
logical(kind=iwp), intent(in) :: converged(8)
integer(kind=iwp) :: iDis, iDisk, iDisp, iDum, ielec(3), iLen, Indx, iOpt, ipCIP1, ipCIP2, ipRP1, ipRP2, ipSP, iRC, iSym, jDisp, &
                     jSpin, kDisp, kSpin, kSym, ldisp2(8), Length, Lu_10, mSym, nConfm, nHss, Pstate_sym
logical(kind=iwp) :: CI
real(kind=wp) :: Fact, Pola(6), rTempC1, rTempC2, rTempC3, rTempK1, rTempK2, rTempK3
character(len=8) :: Label
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: ip
logical(kind=iwp) :: elec_On
character(len=20) :: Label2
#endif
integer(kind=iwp), allocatable :: DegDisp(:), NrDisp(:)
real(kind=wp), allocatable :: EG(:), ELEC(:), ELOUT(:), Hess(:), Hess2(:), Kap1(:), Kap2(:), RHss(:), rKap1(:), rKap2(:), sKap(:), &
                              Temp(:)
integer(kind=iwp), external :: IsFreeUnit
real(kind=wp), external :: DDot_

!                                                                      *
!***********************************************************************
!                                                                      *

#ifdef _DEBUGPRINT_
debug = .true.
#else
debug = .false.
#endif
nHss = size(Hss)
nhess = nTri_Elem(nDisp)
call mma_allocate(RHss,nHss,Label='RHss')
RHss(:) = Zero

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
  PState_SYM = Mul(State_Sym,iSym)
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
    if (btest(nTPert(idisp),0)) jSpin = 1
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
      Length = nDensC
      call dDaFile(LuTemp,2,Kap1,Length,iDisk)
      iDisk = iSigDisp(iDisp)
      call dDaFile(LuTemp,2,SKap,Length,iDisk)
      iDisk = iRHSDisp(iDisp)
      call dDaFile(LuTemp,2,rKap1,Length,iDisk)
      SKap(:) = -SKap(:)-rKap1(:)

      !call Recprt('ORB-RHS',' ',rKap1,nDensC,1)
      !write(u6,*) 'ddot orb-resp',ddot_(nDensC,Kap1,1,Kap1,1)
      !write(u6,*) 'ddot orb-sigma',ddot_(nDensC,SKap,1,SKap,1)
      !write(u6,*) 'ddot orb-rhs',ddot_(nDensC,rKap1,1,rKap1,1)

      call GADSum(Kap1,Length)
      call GADSum(SKap,Length)
      call GADSum(rKap1,Length)

      if (CI) then
        ilen = nconf1
        idis = iCIDisp(iDisp)
        call ipin(ipCIp1)
        call dDaFile(LuTemp,2,W(ipCIp1)%A,iLen,iDis)
        idis = iCISigDisp(idisp)
        call ipin(ipSp)
        call dDaFile(LuTemp,2,W(ipSp)%A,iLen,iDis)
        idis = iRHSCIDisp(idisp)
        call ipin(iprp1)
        call dDaFile(LuTemp,2,W(iprp1)%A,iLen,iDis)
        call ipin(ipSp)
        call ipin(iprp1)
        W(ipSp)%A(1:nConf1) = -W(ipSp)%A(1:nConf1)-W(iprp1)%A(1:nConf1)

        !write(u6,*) 'ddot ci-resp',ddot_(nConf1,W(ipcip1)%A,1,W(ipcip1)%A,1)
        !write(u6,*) 'ddot ci-sigma',ddot_(nConf1,W(ipSp)%A,1,W(ipSp)%A,1)
        !write(u6,*) 'ddot ci-rhs',ddot_(nConf1,W(iprp1)%A,1,W(iprp1)%A,1)

        call GADSum(W(ipCIp1)%A,iLen)
        call GADSum(W(ipSp)%A,iLen)
        call GADSum(W(ipRp1)%A,iLen)
      end if

    else

      Length = nDensC
      Kap1(1:Length) = Zero
      call GADSum(Kap1,Length)
      sKap(1:Length) = Zero
      call GADSum(SKap,Length)
      rKap1(1:Length) = Zero
      call GADSum(rKap1,Length)

      if (CI) then

        ilen = nconf1
        call ipin(ipCIp1)
        W(ipCIp1)%A(1:iLen) = Zero
        call GADSum(W(ipCIp1)%A,iLen)
        call ipin(ipSp)
        W(ipSp)%A(1:iLen) = Zero
        call GADSum(W(ipSp)%A,iLen)
        call ipin(ipRp1)
        W(ipRp1)%A(1:iLen) = Zero
        call GADSum(W(ipRp1)%A,iLen)

      end if

    end if

    !*******************************************************************

    do kDisp=1,jdisp

      !         (x) (2) (y)   (x) (y)    (y) (x)
      ! E    = k   E   k   + F   k   +  F   k
      !  Resp

      kspin = 0
      if (btest(nTPert(kdisp+ksym),0)) kSpin = 1
      if (kspin == 0) then
        nconf1 = ncsf(PState_Sym)
      else
        nConf1 = nint(xispsm(Pstate_Sym,1))
      end if
      if (.not. lCalc(kDisp+ksym)) cycle

      !write(u6,*) 'kDisp+kSym',kDisp+kSym
      !write(u6,*) 'iKapDisp(kdisp+ksym)',iKapDisp(kdisp+ksym)

      iDisk = iKapDisp(kDisp+kSym)
      if (iDisk /= -1) then
        Length = nDensC
        call dDaFile(LuTemp,2,Kap2,Length,iDisk)
        iDisk = iRHSDisp(kDisp+kSym)
        call dDaFile(LuTemp,2,rKap2,Length,iDisk)

        call GASync()
        call GADSum(Kap2,Length)
        call GADSum(rKap2,Length)

        if (CI) then
          ilen = nconf1
          idis = iCIDisp(kDisp+ksym)
          call ipin(ipCIp2)
          call dDaFile(LuTemp,2,W(ipCIp2)%A,iLen,iDis)
          idis = iRHSCIDisp(kdisp+ksym)
          call ipin(iprp2)
          call dDaFile(LuTemp,2,W(iprp2)%A,iLen,iDis)
          call ipin(ipsp)
          rTempc1 = DDot_(nConf1,W(ipCIp2)%A,1,W(ipsp)%A,1)

          call GASync() ! <----------------- NOTE!
          call GADSum(W(ipCIp2)%A,iLen)
          call GADSum(W(iprp2)%A,iLen)

        else
          rtempc1 = Zero
        end if

      else

        call GASync()
        Length = nDensC
        Kap2(1:Length) = Zero
        call GADSum(Kap2,Length)
        rKap2(1:Length) = Zero
        call GADSum(rKap2,Length)
        if (CI) then
          ilen = nconf1
          call GASync()   ! <----------------- NOTE!
          call ipin(ipCIp2)
          W(ipCIp2)%A(1:iLen) = Zero
          call GADSum(W(ipCIp2)%A,iLen)
          call ipin(iprp2)
          W(iprp2)%A(1:iLen) = Zero
          call GADSum(W(iprp2)%A,iLen)
          call ipin(ipsp)
          rTempc1 = DDot_(nConf1,W(ipCIp2)%A,1,W(ipsp)%A,1)
        else
          rtempc1 = Zero
        end if

      end if

      rTempk1 = DDot_(nDensC,Kap2,1,SKap,1)

      Fact = One
      if (kdisp == jdisp) Fact = Two
      rTempk2 = Fact*DDot_(nDensC,Kap1,1,rKap2,1)
      if (kdisp /= jdisp) then
        rtempk3 = DDot_(nDensC,rKap1,1,Kap2,1)
      else
        rTempk3 = Zero
      end if
      if (CI) then
        Fact = One
        if (kdisp == jdisp) Fact = Two
        call ipin(ipCip1)
        call ipin(iprp2)
        rTempc2 = Fact*DDot_(nConf1,W(ipCip1)%A,1,W(iprp2)%A,1)
        if (kdisp /= jdisp) then
          call ipin(iprp1)
          call ipin(ipCIp2)
          rtempc3 = DDot_(nConf1,W(iprp1)%A,1,W(ipCIp2)%A,1)
        else
          rTempc3 = Zero
        end if
      else
        rtempc2 = Zero
        rtempc3 = Zero
      end if

      !write(u6,*) kdisp,jdisp
      !write(u6,*) rTempk1,rtempk2,rtempk3
      !write(u6,*) rtempc1,rtempc2,rtempc3

      Indx = mSym+iTri(kDisp,jDisp)

      Rhss(Indx) = Rhss(Indx)+rTempk1+rtempk2+rtempk3+rtempc1+rtempc2+rtempc3

    end do

    !*******************************************************************

  end do
  kSym = kSym+lDisp(iSym)
  mSym = mSym+nTri_Elem(lDisp(iSym))

  ! Free areas for scratch and state variables

  call mma_deallocate(rKap2)
  call mma_deallocate(rKap1)
  call mma_deallocate(sKap)
  call mma_deallocate(Kap2)
  call mma_deallocate(Kap1)
  if (CI) call ipclose(ipcip1)
end do

call mma_allocate(Hess,nHss,Label='Hess')
Hess(:) = Zero
call mma_allocate(Hess2,nHss,Label='Hess2')
call mma_allocate(Temp,nHss,Label='Temp')
Temp(:) = Zero
call mma_allocate(ELEC,3*ndisp,Label='ELEC')
ELEC(:) = Zero
call mma_allocate(EG,3*ndisp,Label='EG')
call mma_allocate(ELOUT,3*ndisp,Label='ELOUT')
call ipclose(-1)

!----------------------------------------------------------------------*
!
!     OK now when we have out Hessian, what should we do with it!
!
!----------------------------------------------------------------------*

! If a basis set is dependent on perturbation add terms
! constructed in mckinley.

pola(:) = Zero
idum = 1
iopt = ibset(0,sLength)
irc = 3*ndisp
Label = 'DOTELGR'
call drdMCk(irc,iopt,LaBeL,idum,EG,idum)
call GADsum(Hss,nHss)
Hess2(:) = Hss(:)
#ifdef _DEBUGPRINT_
elec_On = .true.
if (irc /= 0) elec_On = .false.
if (debug) then
  ip = 1
  do iSym=1,nSym
    write(label2,'(A,I2)') 'CHessian symmetry',iSym
    if (lDisp(iSym) /= 0) call TriPrt(label2,' ',Hess2(ip),lDisp(iSym))
    ip = ip+nTri_Elem(ldisp(isym))
  end do
end if

if (debug) then
  call MMSORT2(HESS2,ELEC,pola,ielec)
  call Recprt('CONN',' ',Elec,3*nDisp,1)
end if

!write(u6,*) 'I am here 1'
call Recprt('Rhss','(5G20.10)',RHss,nhss,1)
call Recprt('Hss','(5G20.10)',Hss,nHss,1)
#endif

Hess2(1:mSym) = Hess2(1:mSym)+RHss(1:mSym)

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
    ip = ip+nTri_Elem(ldisp(isym))
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
    write(u6,*) 'OutPut: Error reading MCKINT'
    write(u6,'(A,A)') 'Label=',Label
    call Abend()
  end if

# ifdef _DEBUGPRINT_
  if (debug) then
    ip = 1
    do iSym=1,nSym
      write(label2,'(a,i2)') 'SHessian symmetry',iSym
      if (lDisp2(iSym) /= 0) call TriPrt(label2,' ',Temp(ip),lDisp2(iSym))
      ip = ip+nTri_Elem(ldisp2(isym))
    end do
  end if
# endif
  Hess(1:mSym) = Hess(1:mSym)+Temp(1:mSym)
end if
#ifdef _DEBUGPRINT_
if (debug) then
  ip = 1
  do iSym=1,nSym
    write(label2,'(a,i2)') 'Hessian symmetry',iSym
    if (lDisp2(iSym) /= 0) call TriPrt(label2,' ',Hess(ip),lDisp2(iSym))
    ip = ip+nTri_Elem(ldisp2(isym))
  end do
end if
#endif

if (McKinley) then
  iRC = -1
  iOpt = 0
  Label = 'Hess'
  call dWrMck(iRC,iOpt,Label,iDum,Hess,iDum)
  if (iRC /= 0) then
    write(u6,*) 'OutPut: Error writing to MCKINT'
    write(u6,'(A,A)') 'Label=',Label
    call Abend()
  end if
  call Put_iScalar('No of Internal coordinates',ldisp2(1))
  call Put_AnalHess(Hess,nTri_Elem(ldisp2(1)))
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
    write(u6,*) 'OutPut: Error reading RELAX'
    write(u6,'(A,A)') 'Label=',Label
    call Abend()
  end if

  !if (debug) call HssPrt_MCLR(DegDisp,Hess,ldisp2)
  !call Recprt('hess',' ',Hess,nhss,1)

  ELEC(:) = ELEC(:)-EG(:)
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
  write(u6,*)
  write(u6,*)
  write(u6,*) '************************************'
  write(u6,*) '*                                  *'
  write(u6,*) '*       Polarizabilities           *'
  write(u6,*) '*                                  *'
  write(u6,*) '************************************'
  write(u6,*)
  write(u6,*)
  call Add_Info('POLARIZABILITIES',Pola,6,2)

  ! Go from energy derivative to polarizability, there is a difference
  ! in the sign in the definition.

  Pola(:) = -Pola(:)

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
