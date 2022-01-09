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
!...  Matrix of the atomic quasi-relativistic potential in the
!     1-center basis set. Mass-velocity and Darwin potentials,
!     which are read from the library in atomic units.

implicit real*8(A-H,O-Z)
character OrbLab*4, MPLbl*20, Keyw*40, RdName*40, DumCha*40, Line*80
real*8 ZETA(NZ), OPMAT(NZ*(NZ+1)/2)
! Internal
parameter(MPoint=250)
real*8 R(Mpoint), V(Mpoint), FI(Mpoint), FIVFJ(Mpoint)
! auxiliary constant pool:    only DFAC is used here.
#include "const.fh"
external RDNAME, SIMPLM
data PI/3.141592653589793d0/
! Flags: MV bit 2, DW bit 3
!data iMVPot/2/,iDWPot/4/

!LU6 = 6
N = ISIM

! Reads the Mass Velocity and Darwin potentials for the valence orbitals

!WRITE(LU6,600) ISIM,LUQRP

! locate the beginning of the data
KEYW = MPLBL
DUMCHA = RDNAME(LUQRP,KEYW)
if (DUMCHA == ' ') then
  write(6,'(1X,A20," MV & DW potentials not found in unit ",I3)') KEYW,LUQRP
  call Quit_OnUserError()
end if
!write(6,*) 'Data found in libray. Start reading.'
read(LUQRP,*) Npoint
if (NPOINT > MPoint) then
  write(6,*) 'VQR: nPoint',nPoint
  write(6,*) 'Abend: Increade mPoint'
  call Quit_OnUserError()
end if
read(LUQRP,'(4d20.13)') (R(k),k=1,Npoint)

1 read(LUQRP,'(a80)',end=999) Line
if (Line(1:8) == MPLbl(1:8)) then
  ! there is no potential for this symmetry
  !write(6,603)
  do IJ=1,NZ*(NZ+1)/2
    OPMAT(IJ) = 0d0
  end do
  return
else
  read(Line,'(a4)') OrbLab
  read(LUQRP,'(4d20.13)') (V(k),k=1,Npoint)

  ! check the symmetry
  if (((index(ORBLAB,'S') /= 0) .and. (ISIM /= 1)) .or. ((index(ORBLAB,'P') /= 0) .and. (ISIM /= 2)) .or. &
      ((index(ORBLAB,'D') /= 0) .and. (ISIM /= 3)) .or. ((index(ORBLAB,'F') /= 0) .and. (ISIM /= 4))) then
    GO TO 1
  end if
end if
!write(6,601) ORBLAB,NPOINT,R(1),R(NPOINT)

! Matrix of a numerical local operator V(r) given in a logarithmic mesh

PREN = 2d0**(N+1)/(sqrt(DFAC(2*N)*sqrt(2d0*PI)))
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
    OPMAT(IJ) = RES+FIVFJ(1)*R(1)/3d0
  end do
end do

return
999 continue
write(6,*) ' Troubles reading MV and DW potentials from QRP library'
call Quit_OnUserError()

!600   FORMAT (/1X,'Symmetry ',I1,' : ','Mass Velocity and Darwin potentials are read from unit ',I2)
!601   FORMAT (/' Orbital ',A4,' : ',I5,' points in the logarithmic mesh from ',F10.6,' to ',F10.6)
!602   FORMAT (/1X,' The symmetry label of one of the relativistic potentials is not S, P, D, or F')
!603   FORMAT (/1X,' There is no potential for this symmetry.')

end subroutine Vqr
