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
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

subroutine getci_cvb(civec)

use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: civec(*)
#include "main_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"
#include "WrkSpc.fh"
#include "casinfo_cvb.fh"
#include "io_cvb.fh"
integer(kind=iwp) :: ibf, icivec, istate, istsym_d, isyml, iwr, lcim, nci, ncix(mxirrep)
real(kind=wp) :: cnrm, fac
integer(kind=iwp), external :: igetcnt2_cvb, mstackr_cvb
real(kind=wp), external :: dnrm2_
logical(kind=iwp), external :: ifcasci_cvb, & ! ... Files/Hamiltonian available ...
                               valid_cvb

icivec = nint(civec(1))
if (igetcnt2_cvb(icivec) == 1) return
if (.not. ifcasci_cvb()) return
call setcnt2_cvb(icivec,1)
iwr = 0

if (iform_ci(icivec) /= 0) then
  write(u6,*) ' Unsupported format in GETCI :',iform_ci(icivec)
  call abend_cvb()
end if

if (iwr == 0) then
  if (ip(1) >= 1) then
    write(u6,'(a)') ' '
    call prtfid_cvb(' Restoring CI vector from ',strtci)
  end if
  call fzero(work(iaddr_ci(icivec)),ndet)
else if (iwr == 1) then
  if ((ip(5) >= 1) .and. valid_cvb(savvbci)) then
    write(u6,'(a)') ' '
    call prtfid_cvb(' Saving VB CI vector to ',savvbci)
  end if
end if

do istsym_d=1,nstsym_d
  isyml = istsy_d(istsym_d)
  call getnci_cvb(ncix,istnel_d(istsym_d),istms2_d(istsym_d),istsy_d(istsym_d))
  nci = ncix(1)
  lcim = mstackr_cvb(nci)
  if (iwr == 0) then
    do istate=1,nstats_d(istsym_d)
      if (abs(weight_d(istate,istsym_d)) > 1.0e-20_wp) then
        call mkfn_cvb(strtci,ibf)
        call rdcivec_cvb(work(lcim),filename(ibf),.true.)
        fac = sqrt(weight_d(istate,istsym_d))
        call mol2vbma_cvb(work(iaddr_ci(icivec)),work(lcim),isyml,fac)
      end if
    end do
  else if (iwr == 1) then
    do istate=1,nstats_d(istsym_d)
      if (abs(weight_d(istate,istsym_d)) > 1.0e-20_wp) then
        call vb2mol_cvb(work(iaddr_ci(icivec)),work(lcim),isyml)
        cnrm = one/dnrm2_(nci,work(lcim),1)
        call dscal_(nci,cnrm,work(lcim),1)
        call mkfn_cvb(savvbci,ibf)
        call wrcivec_cvb(work(lcim),filename(ibf),.not. variat)
      end if
    end do
  end if
  call mfreer_cvb(lcim)
end do

return

end subroutine getci_cvb
