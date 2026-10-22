!u**********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2006, Per Ake Malmqvist                                *
!               2021, Yoshio Nishimoto                                 *
!***********************************************************************
!--------------------------------------------*
! 2006  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*

subroutine DERSPE(NLEV,NG3,DF1,DF2,DF3,idxG3,DEPSA)

#ifdef _MOLCAS_MPP_
use Para_Info, only: King
#endif
use caspt2_module, only: ISCF
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp, byte

implicit none
integer(kind=iwp), intent(in) :: NLEV, NG3
real(kind=wp), intent(in) :: DF1(NLEV,NLEV), DF2(NLEV,NLEV,NLEV,NLEV), DF3(NG3)
integer(kind=byte), intent(in) :: idxG3(6,NG3)
real(kind=wp), intent(inout) :: DEPSA(NLEV,NLEV)
integer(kind=iwp) :: iA, iB, iG3, iT, iV, iX
real(kind=wp) :: Diagonal, OCC

! The active orbitals are all occupied with OCC electrons, so E_pq|0> = OCC*delta(p,q)|0>.
!   dF1(t,u)/df(a,b)         = OCC**2*d(t,u)*d(a,b) - OCC*d(t,b)*d(a,u)
!   dF2(t,u,v,x)/df(a,b)     = OCC**3*d(t,u)*d(v,x)*d(a,b)
!                            - OCC**2*(d(t,x)*d(v,u)*d(a,b)+d(t,b)*d(v,x)*d(a,u)+d(t,u)*d(v,b)*d(a,x))
!                            + OCC*(d(t,x)*d(v,b)*d(a,u)+d(t,b)*d(v,u)*d(a,x))
!   dF3(t,u,v,x,y,z)/df(a,b) = OCC**4*d(t,u)*d(v,x)*d(y,z)*d(a,b) - OCC**3*d(t,u)*d(y,z)*d(a,x)*d(b,v)

OCC = Two
if (ISCF == 2) OCC = One

Diagonal = Zero
#ifdef _MOLCAS_MPP_
if (King()) then
#endif
  do iT=1,NLEV
    Diagonal = Diagonal+OCC**2*DF1(iT,iT)
    do iV=1,NLEV
      Diagonal = Diagonal+OCC**3*DF2(iT,iT,iV,iV)-OCC**2*DF2(iT,iV,iV,iT)
    end do
  end do
  do iB=1,NLEV
    do iA=1,NLEV
      DEPSA(iA,iB) = DEPSA(iA,iB)-OCC*DF1(iB,iA)
      do iT=1,NLEV
        DEPSA(iA,iB) = DEPSA(iA,iB)-OCC**2*(DF2(iB,iA,iT,iT)+DF2(iT,iT,iB,iA))+OCC*(DF2(iT,iA,iB,iT)+DF2(iB,iT,iT,iA))
      end do
    end do
  end do
#ifdef _MOLCAS_MPP_
end if
#endif

do iG3=1,NG3
  if ((idxG3(1,iG3) /= idxG3(2,iG3)) .or. (idxG3(5,iG3) /= idxG3(6,iG3))) cycle
  iV = idxG3(3,iG3)
  iX = idxG3(4,iG3)
  if (iV == iX) Diagonal = Diagonal+OCC**4*DF3(iG3)
  DEPSA(iX,iV) = DEPSA(iX,iV)-OCC**3*DF3(iG3)
end do

do iA=1,NLEV
  DEPSA(iA,iA) = DEPSA(iA,iA)+Diagonal
end do

end subroutine DERSPE
