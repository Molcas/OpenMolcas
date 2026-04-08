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
! Copyright (C) 1995, Jeppe Olsen                                      *
!***********************************************************************

subroutine GETSTRN_GASSM_SPGP(ISMFGS,ITPFGS,ISTROC,NSTR,NEL,NNSTSGP,IISTSGP)
! Obtain all superstrings containing  strings of given sym and type
!
! (Superstring :contains electrons belonging to all gasspaces
!       string :contains electrons belonging to a given GAS space
!  A super string is thus a product of NGAS strings)
!
! Jeppe Olsen, Summer of 95
!              Optimized version, october 1995
!
! In this subroutine the ordering of strings belonging to a given type
! is defined !!
! Currently we are using the order
! Loop over GAS 1 strings
!  Loop over GAS 2 strings
!   Loop over GAS 3 strings --
!     Loop over gas N strings

use lucia_data, only: MXPNGAS, MXPNSMST, NELFGP, NGAS, OCCSTR
use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: ISMFGS(*), ITPFGS(*), NEL, NNSTSGP(MXPNSMST,*), IISTSGP(MXPNSMST,*)
integer(kind=iwp), intent(inout) :: ISTROC(NEL,*)
integer(kind=iwp), intent(out) :: NSTR
integer(kind=iwp) :: IBSTFGS(MXPNGAS), IGAS, IGASL, NELB, NELI, NSTA, NSTB, NSTFGS(MXPNGAS), NSTI, NSTRTOT

! Number of strings per GAS space
!write(u6,*) ' entering problem child'
do IGAS=1,NGAS
  NSTFGS(IGAS) = NNSTSGP(ISMFGS(IGAS),IGAS)
  IBSTFGS(IGAS) = IISTSGP(ISMFGS(IGAS),IGAS)
end do

#ifdef _DEBUGPRINT_
write(u6,*) '  GETSTR_GASSM_SPGP speaking'
write(u6,*) '  =========================='
write(u6,*) ' ISMFGS,ITPFGS (input)'
call IWRTMA(ISMFGS,1,NGAS,1,NGAS)
call IWRTMA(ITPFGS,1,NGAS,1,NGAS)
write(u6,*)
write(u6,*) ' NSTFGS, IBSTFGS (intermediate results)'
call IWRTMA(NSTFGS,1,NGAS,1,NGAS)
call IWRTMA(IBSTFGS,1,NGAS,1,NGAS)
#endif
! Last gasspace with a nonvanishing number of electrons
IGASL = 0
do IGAS=1,NGAS
  if (NELFGP(ITPFGS(IGAS)) /= 0) IGASL = IGAS
end do

NSTRTOT = product(NSTFGS(1:NGAS))

if ((IGASL /= 0) .and. (NSTRTOT /= 0)) then
  ! Loop over GAS spaces
  NELB = 0
  do IGAS=1,IGASL
    ! Number of electrons in GAS = 1, IGAS - 1
    if (IGAS > 1) NELB = NELB+NELFGP(ITPFGS(IGAS-1))
    ! Number of electron in IGAS
    NELI = NELFGP(ITPFGS(IGAS))
    if (NELI > 0) then

      ! The order of strings corresponds to a matrix A(I(after),Igas,I(before))
      ! where I(after) loops over strings in IGAS+1 - IGASL and
      ! I(before) loop over strings in 1 - IGAS -1
      NSTA = product(NSTFGS(IGAS+1:IGASL))

      NSTB = product(NSTFGS(1:IGAS-1))

      NSTI = NSTFGS(IGAS)

#     ifdef _DEBUGPRINT_
      write(u6,*) ' NSTI,NSTB,NSTA,NELB,NELI,NEL ',NSTI,NSTB,NSTA,NELB,NELI,NEL
      write(u6,*) ' IBSTFGS(IGAS)',IBSTFGS(IGAS)
#     endif

      call ADD_STR_GROUP(NSTI,IBSTFGS(IGAS),OCCSTR(ITPFGS(IGAS))%A,NSTB,NSTA,ISTROC,NELB+1,NELI,NEL)

      ! Loop over strings in IGAS
    end if
  end do
end if
NSTR = NSTRTOT

#ifdef _DEBUGPRINT_
write(u6,*) ' Info from  GETSTR_GASSM_SPGP'
write(u6,*) ' ============================'
write(u6,*)
write(u6,*) ' Symmetry and type strings :'
write(u6,*)
write(u6,*) '   AS    Sym  Type'
write(u6,*) ' ================='
do IGAS=1,NGAS
  write(u6,'(3I6)') IGAS,ISMFGS(IGAS),ITPFGS(IGAS)
end do
write(u6,*)
write(u6,*) ' Number of strings generated : ',NSTR
write(u6,*) ' Strings generated'
call PRTSTR(ISTROC,NEL,NSTR)
#endif

end subroutine GETSTRN_GASSM_SPGP
