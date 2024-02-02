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

subroutine MAGN(EXCH,N,X,Y,Z,H,W,zJ,THRS,dM,sM,nT,T,sopt,WZ,ZB,S,M,m_paranoid,DBG)
! this Subroutine is a wrapper for various MAGN subroutines

implicit none
integer, parameter :: wp = kind(0.d0)
integer, intent(in) :: EXCH, N, nT
real(kind=8), intent(in) :: X, Y, Z, H, zJ
real(kind=8), intent(in) :: W(EXCH), T(nT)
complex(kind=8), intent(in) :: dM(3,EXCH,EXCH)
complex(kind=8), intent(in) :: sM(3,EXCH,EXCH)
logical, intent(in) :: sopt
real(kind=8), intent(out) :: ZB(nT), WZ(N)
real(kind=8), intent(out) :: S(3,nT), M(3,nT)
real(kind=8), intent(in) :: THRS
logical, intent(in) :: m_paranoid
logical, intent(in) :: DBG

if (abs(zJ) < tiny(0.0_wp)) then

  if (DBG) write(6,*) 'Enter MAGN_NO_MF :'

  call MAGN_NO_MF(EXCH,N,X,Y,Z,H,W,dM,sM,nT,T,sopt,WZ,ZB,S,M,DBG)

  if (DBG) write(6,*) 'Exit MAGN_NO_MF :'

else ! zJ /= 0.0_wp

  if (DBG) write(6,*) 'Enter MAGN_ZJ_PAR :'

  call MAGN_ZJ_PAR(EXCH,N,X,Y,Z,H,W,zJ,dM,sM,nT,T,sopt,WZ,ZB,S,M,thrs,m_paranoid,DBG)

  if (DBG) write(6,*) 'Exit MAGN_ZJ_PAR :'

end if

return

end subroutine MAGN
