!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine Read_Multipole_Int(lMax,ip_sq_mu,nBas,ip_mu,Ttot,Temp,Origin,rMPq,nElem,nBas1,nBas2,nBasMax,nTemp,nSym,ipP,Restart, &
                              Utility)

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: lMax, nElem, ip_sq_mu(0:nElem-1), nSym, nBas(0:nSym-1), ip_mu(0:nElem-1), nBas1, nBas2, nBasMax, nTemp, ipP
real(kind=wp) :: Ttot(nBas2), Temp(nTemp), Origin(3,0:lMax), rMPq(0:nElem-1)
logical(kind=iwp) :: Restart, Utility
integer(kind=iwp) :: iComp, idum(1), ijSym, iOff, iOffs, iOfft, iOpt0, iOpt1, ip_all_ints, ip_iSyLbl, ip_nComp, ip_Tmp, ipScr, &
                     iRc, iSyLbl, iSym, jSym, l, mu, nComp, nInts, nInts_Tot, nScr
character(len=16) :: RunFile_dLabel, RunFile_iLabel, RunFile_iLabel2
character(len=8) :: Label
logical(kind=iwp) :: Found
#include "WrkSpc.fh"

!                                                                      *
!***********************************************************************
!                                                                      *
! Note that we will store two sets of the integrals.
!
! Set 1: used to compute the localized properties.
! Set 2: used in the finite differenc calculation of the
!        polarizability.
!
! Set 1 will be desymmetrized in case of symmetry.
!
! Pointers to Set 1 are stored in ip_sq_mu.
! Pointers to Set 2 are stored in ip_mu.

RunFile_dLabel = 'LoProp Integrals'
RunFile_iLabel = 'LoProp nInts'
RunFile_iLabel2 = 'LoProp iSyLbl'
nInts_Tot = 0
call Allocate_iWork(ip_nComp,nElem)
call Allocate_iWork(ip_iSyLbl,nElem)
if (Restart) then
  call Qpg_dArray(RunFile_dLabel,Found,nInts_tot)
  if (.not. Found) then
    write(u6,*) 'LoProp Integrals not available on the RunFile.'
    call Abend()
  end if
  call Allocate_Work(ip_all_ints,nInts_Tot)
  call Get_dArray(RunFile_dLabel,Work(ip_all_ints),nInts_tot)
  call Get_iArray(RunFile_iLabel,iWork(ip_nComp),nElem)
  call Get_iArray(RunFile_iLabel2,iWork(ip_iSyLbl),nElem)
end if
nInts = 0
iOpt0 = 0
iOpt1 = 1
Label = 'Mltpl  X'
mu = -1
iOff = 0
do l=0,lMax
  nComp = (l+1)*(l+2)/2
  write(Label(8:8),'(I1)') l
  do iComp=1,nComp
#   ifdef _DEBUGPRINT_
    write(u6,*) 'l,iComp=',l,iComp
#   endif
    mu = mu+1
    if (Restart) then
      call Allocate_Work(ip_mu(mu),iWork(ip_nComp+mu))
      nInts = iWork(ip_nComp+mu)-4
      call dCopy_(iWork(ip_nComp+mu),Work(ip_all_ints+iOff),1,Work(ip_mu(mu)),1)
      iSyLbl = iWork(ip_iSyLbl+mu)
      iOff = iOff+iWork(ip_nComp+mu)
    else
      iRc = -1
      iSyLbl = 0
      call iRdOne(iRc,iOpt1,Label,iComp,idum,iSyLbl)
      if (iRc /= 0) then
        write(u6,*) 'Polar: error reading length of mu!'
        write(u6,*) 'Mu=',mu
        call Abend()
      end if
      nInts = idum(1)
      call Allocate_Work(ip_mu(mu),nInts+4)
      call RdOne(iRc,iOpt0,Label,iComp,Work(ip_mu(mu)),iSyLbl)
      if (iRc /= 0) then
        write(u6,*) 'Polar: error reading mu!'
        write(u6,*) 'Mu=',mu
        call Abend()
      end if
      iWork(ip_iSyLbl+mu) = iSyLbl
      iWork(ip_nComp+mu) = nInts+4
      nInts_Tot = nInts_Tot+iWork(ip_nComp+mu)
    end if

    ! Transform multipole moment integrals to new basis

    call Allocate_Work(ip_Tmp,nBas1**2)
    call Fzero(Work(ip_Tmp),nBas1**2)
    if (nSym == 1) then
      ip_sq_mu(mu) = ip_Tmp
    else
      call GetMem('SMatrix','Allo','Real',ip_sq_mu(mu),nBas1**2)
    end if

    ! Now I square the dipole moment integral matrix because it
    ! has been stored as a low triangle

    iOfft = ip_mu(mu)
    iOffs = ip_Tmp
    do iSym=0,nSym-1
      do jSym=0,iSym
        ijSym = ieor(iSym,jSym)
        if (iand(iSyLbl,2**ijSym) == 0) then
          Go To 20
        end if
        if (nBas(iSym)*nBas(jSym) == 0) then
          Go To 20
        end if
        if (iSym == jSym) then

          ! Now I square the overlap matrix because it has been
          ! stored as a lower triangle

          call Square(Work(iOfft),Work(iOffs),1,nBas(iSym),nBas(iSym))

          iOfft = iOfft+nBas(iSym)*(nBas(iSym)+1)/2
        else
          call dcopy_(nBas(iSym)*nBas(jSym),Work(iOfft),1,Work(iOffs),1)
          iOfft = iOfft+nBas(iSym)*nBas(jSym)
        end if
        iOffs = iOffs+nBas(iSym)*nBas(jSym)
20      continue
      end do
    end do

    if (nSym /= 1) then

      ! Desymmetrize

      nScr = nBasMax*nBas1
      call Allocate_Work(ipScr,nScr)
      call FZero(Work(ip_sq_mu(mu)),nBas1**2)
      call Desymmetrize(Work(ip_Tmp),nBas2,Work(ipScr),nScr,Work(ip_sq_mu(mu)),nBas,nBas1,Work(ipP),nSym,iSyLbl)
      call Free_Work(ipScr)
      call Free_Work(ip_Tmp)

    end if

    ! Now I transform with my Transformation Matrix Ttot the
    ! multipole moment integral matrices

#   ifdef _DEBUGPRINT_
    call RecPrt('Multipole Integrals in AO Basis',' ',Work(ip_sq_mu(mu)),nBas1,nBas1)
#   endif
    call TransMu(Work(ip_sq_mu(mu)),nBas1,Ttot,Temp)
#   ifdef _DEBUGPRINT_
    call RecPrt('Multipole Integrals in LoProp Basis',' ',Work(ip_sq_mu(mu)),nBas1,nBas1)
#   endif

    ! Pick up the nuclear value for this operator

    rMPq(mu) = Work(ip_mu(mu)+nInts+3)
  end do

  ! Pick up the origin of this operator

  call dcopy_(3,Work(ip_mu(mu)+nInts),1,Origin(1,l),1)

end do
if ((.not. Restart) .and. (.not. Utility)) then
  call Allocate_Work(ip_all_ints,nInts_Tot)
  mu = -1
  iOff = 0
  do l=0,lMax
    nComp = (l+1)*(l+2)/2
    do iComp=1,nComp
      mu = mu+1
      call dCopy_(iWork(ip_nComp+mu),Work(ip_mu(mu)),1,Work(ip_all_ints+iOff),1)
      iOff = iOff+iWork(ip_nComp+mu)
    end do
  end do
  call Put_dArray(RunFile_dLabel,Work(ip_all_ints),nInts_Tot)
  call Put_iArray(RunFile_iLabel,iWork(ip_nComp),nElem)
  call Put_iArray(RunFile_iLabel2,iWork(ip_iSyLbl),nElem)
end if
if (.not. Utility) then
  call Free_Work(ip_all_ints)
end if
call Free_iWork(ip_nComp)
call Free_iWork(ip_iSyLbl)
#ifdef _DEBUGPRINT_
call RecPrt('Origin',' ',Origin,3,lMax+1)
call RecPrt('rMPq',' ',rMPq,1,nElem)
call xSpot('Exit  Read_Multipole_Int')
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Read_Multipole_Int
