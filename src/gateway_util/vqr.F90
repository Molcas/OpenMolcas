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

subroutine Vqr(LUQRP,MPLBL,ISIM,ZETA,NZ,OPMAT)
! Matrix of the atomic quasi-relativistic potential in the
! 1-center basis set. Mass-velocity and Darwin potentials,
! which are read from the library in atomic units.

use Constants, only: Zero, Two, Three, Pi
use AMatrix, only: DFAC
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: LUQRP, ISIM, NZ
character(len=20), intent(in) :: MPLbl
real(kind=wp), intent(in) :: ZETA(NZ)
real(kind=wp), intent(out) :: OPMAT(NZ*(NZ+1)/2)
integer(kind=iwp), parameter :: MPoint = 250
integer(kind=iwp) :: I, IJ, istatus, J, k, N, N2P1, Npoint
character(len=80) :: Line
character(len=40) :: DumCha, Keyw
character(len=4) :: OrbLab
real(kind=wp) :: ERRE, ERRE2, FI(Mpoint), FIVFJ(Mpoint), PREN, R(Mpoint), RES, RNORI, RNORJ, V(Mpoint), ZI, ZJ
real(kind=wp), external :: SIMPLM
character(len=40), external :: RdName

N = ISIM

! Reads the Mass Velocity and Darwin potentials for the valence orbitals

!WRITE(u6,600) ISIM,LUQRP

! locate the beginning of the data
KEYW = MPLBL
DUMCHA = RDNAME(LUQRP,KEYW)
if (DUMCHA == ' ') then
  write(u6,'(1X,A20," MV & DW potentials not found in unit ",I3)') KEYW,LUQRP
  call Quit_OnUserError()
end if
!write(u6,*) 'Data found in libray. Start reading.'
read(LUQRP,*) Npoint
if (NPOINT > MPoint) then
  write(u6,*) 'VQR: nPoint',nPoint
  write(u6,*) 'Abend: Increade mPoint'
  call Quit_OnUserError()
end if
read(LUQRP,'(4es20.13)') (R(k),k=1,Npoint)

do
  read(LUQRP,'(a80)',iostat=istatus) Line
  if (istatus < 0) call Error()
  if (Line(1:8) == MPLbl(1:8)) then
    ! there is no potential for this symmetry
    !write(u6,603)
    OPMAT(:) = Zero
    return
  else
    read(Line,'(a4)') OrbLab
    read(LUQRP,'(4es20.13)') (V(k),k=1,Npoint)

    ! check the symmetry
    if (((index(ORBLAB,'S') == 0) .or. (ISIM == 1)) .and. ((index(ORBLAB,'P') == 0) .or. (ISIM == 2)) .and. &
        ((index(ORBLAB,'D') == 0) .or. (ISIM == 3)) .and. ((index(ORBLAB,'F') == 0) .or. (ISIM == 4))) exit
  end if
end do
!write(u6,601) ORBLAB,NPOINT,R(1),R(NPOINT)

! Matrix of a numerical local operator V(r) given in a logarithmic mesh

PREN = Two**(N+1)/(sqrt(DFAC(2*N)*sqrt(Two*PI)))
N2P1 = 2*N+1
IJ = 0

do I=1,NZ
  ZI = ZETA(I)
  RNORI = PREN*sqrt(sqrt(ZI**N2P1))
  do K=1,NPOINT
    ERRE = R(K)
    FI(K) = RNORI*ERRE**(N-1)*exp(-ZI*ERRE*ERRE)
  end do

  do J=1,I
    ZJ = ZETA(J)
    RNORJ = PREN*sqrt(sqrt(ZJ**N2P1))
    do K=1,NPOINT
      ERRE = R(K)
      ERRE2 = ERRE*ERRE
      FIVFJ(K) = RNORJ*ERRE**(N-1)*exp(-ZJ*ERRE2)*FI(K)*V(K)*ERRE2
    end do

    RES = SIMPLM(NPOINT,FIVFJ,R)

    ! add up the contribution from 0 to R(1)
    IJ = IJ+1
    OPMAT(IJ) = RES+FIVFJ(1)*R(1)/Three
  end do
end do

return

!600   FORMAT (/1X,'Symmetry ',I1,' : ','Mass Velocity and Darwin potentials are read from unit ',I2)
!601   FORMAT (/' Orbital ',A4,' : ',I5,' points in the logarithmic mesh from ',F10.6,' to ',F10.6)
!602   FORMAT (/1X,' The symmetry label of one of the relativistic potentials is not S, P, D, or F')
!603   FORMAT (/1X,' There is no potential for this symmetry.')

contains

subroutine Error()

  write(u6,*) ' Troubles reading MV and DW potentials from QRP library'
  call Quit_OnUserError()

end subroutine Error

end subroutine Vqr
