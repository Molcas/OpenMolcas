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

subroutine ECP_shells(iAtmNr,nCore,List)

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iAtmNr, nCore
integer(kind=iwp), intent(out) :: List(0:3)
integer(kind=iwp), parameter :: s = 0, p = 1, d = 2, f = 3

!                                                                      *
!***********************************************************************
!                                                                      *
! Try to guess which valence orbitals are explicitly represented
call OrbType(iAtmNr,List,31)
select case (nCore)
  case (0)
    ! No info on number of core electrons, pure guess
    List(:) = 0
    select case (iAtmNr)
      case (0) ! Dummy
      case (1:2) ! H-He
      case (3:4) ! Li-Be
        List(s) = 1
      case (5:10) ! B-Ne
        List(s) = 1
        List(p) = 1
      case (11:12) ! Na-Mg
        List(s) = 1
      case (13:18) ! Al-Ar
        List(s) = 1
        List(p) = 1
      case (19:20) ! K-Ca
        List(s) = 1
      case (21:30) ! Sc-Zn
        List(s) = 1
        List(d) = 1
      case (31:36) ! Ga-Kr
        List(s) = 1
        List(p) = 1
        List(d) = 1
      case (37:38) ! Rb-Sr
        List(s) = 1
      case (39:48) ! Y-Cd
        List(s) = 1
        List(d) = 1
      case (49:54) ! In-Xe
        List(s) = 1
        List(p) = 1
        List(d) = 1
      case (55:56) ! Cs-Ba
        List(s) = 1
      case (57:70) ! La-Yb
        List(s) = 1
        List(f) = 1
      case (71:80) ! Lu-Hg
        List(s) = 1
        List(d) = 1
        List(f) = 1
      case (81:86) ! Tl-Rn
        List(s) = 1
        List(p) = 1
        List(d) = 1
        List(f) = 1
      case (87:88) ! Fr-Ra
        List(s) = 1
      case (89:102) ! Ac-No
        List(s) = 1
        List(f) = 1
      case (103:112) ! Lr-Cn
        List(s) = 1
        List(d) = 1
        List(f) = 1
      case (113:118) ! Nh-Og
        List(s) = 1
        List(p) = 1
        List(d) = 1
        List(f) = 1
      case (119:120) ! Uue-Ubn
        List(s) = 1
      case default
        write(u6,*) 'ECP_shells cannot handle atom numbers beyond 112.'
        call Abend()
    end select
  case (2)
    List(s) = List(s)-1
  case (4)
    List(s) = List(s)-2
  case (10)
    List(s) = List(s)-2
    List(p) = List(p)-1
  case (12)
    List(s) = List(s)-3
    List(p) = List(p)-1
  case (18)
    List(s) = List(s)-3
    List(p) = List(p)-2
  case (28)
    List(s) = List(s)-3
    List(p) = List(p)-2
    List(d) = List(d)-1
  case (36)
    List(s) = List(s)-4
    List(p) = List(p)-3
    List(d) = List(d)-1
  case (46)
    List(s) = List(s)-4
    List(p) = List(p)-3
    List(d) = List(d)-2
  case (47:60)
    ! This includes lanthanides, with f shell in core
    List(s) = List(s)-4
    List(p) = List(p)-3
    List(d) = List(d)-2
    List(f) = List(f)-1
  case (68)
    List(s) = List(s)-5
    List(p) = List(p)-4
    List(d) = List(d)-2
    List(f) = List(f)-1
  case (78)
    List(s) = List(s)-5
    List(p) = List(p)-4
    List(d) = List(d)-3
    List(f) = List(f)-1
  case default
end select
List(:) = max(0,List(:))
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine ECP_Shells
