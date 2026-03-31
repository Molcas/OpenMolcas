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
! Copyright (C) Chad E. Hoyer                                          *
!***********************************************************************
!  DQVDiabat
!
!> @brief
!>   Compute diabats with DQV
!> @author C. E. Hoyer
!>
!> @details
!> This subroutine takes in properties that have been
!> computed by RASSI and uses them to compute diabatic
!> states and thus diabatic energies and couplings.
!> Currently, the user must compute x, y, z, xx,
!> yy, zz, and 1/r, and they must be computed in that
!> order.
!>
!> @param[in] PROP Properties computed in RASSI
!> @param[in] HAM
!***********************************************************************

subroutine DQVDiabat(PROP,HAM)

use Cntrl, only: AlphZ, BetaE, ICOMP, NPROP, NSTATE, PNAME
use Constants, only: Zero, One, Two, Half, Quart, Pi
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: PROP(NSTATE,NSTATE,NPROP), HAM(NSTATE,NSTATE)
integer :: I, IPROP, ISTA, JSTA, K, PNUM(7)
real(kind=wp) :: ATerm, BTerm, Chng, CosO, CTerm, HAMT(NSTATE,NSTATE), HDIA(NSTATE,NSTATE), HDIAI(NSTATE,NSTATE), RotAngF, &
                 RotAngO, SinO, T1, T2, ThrSch, TII, TIJ, TJJ, TROT(NSTATE,NSTATE), TROTT(NSTATE,NSTATE), TRQ(NSTATE,NSTATE)
integer(kind=iwp), parameter :: itMAX = 50
real(kind=wp), parameter :: MTE = 1.0e-8_wp, MTF = 1.0e-14_wp, THRS = 1.0e-8_wp

! Printing some stuff

write(u6,*)
write(u6,*)
write(u6,'(6X,A)') repeat('*',100)
write(u6,'(6X,A,98X,A)') '*','*'
write(u6,'(6X,A,33X,A,34X,A)') '*',' The DQV Diabatization Section ','*'
write(u6,'(6X,A,98X,A)') '*','*'
write(u6,'(6X,A)') repeat('*',100)
write(u6,*)
write(u6,*)

call CollapseOutput(1,'DQV Diabatization section')
write(u6,'(3X,A)') '-------------------------'

! Find the properties we need
do IPROP=1,7
  PNUM(IPROP) = 0
end do
do IPROP=1,NPROP
  if (PNAME(IPROP) == 'MLTPL  1') then
    if (ICOMP(IPROP) == 1) PNUM(1) = IPROP
    if (ICOMP(IPROP) == 2) PNUM(2) = IPROP
    if (ICOMP(IPROP) == 3) PNUM(3) = IPROP
  end if
  if (PNAME(IPROP) == 'MLTPL  2') then
    if (ICOMP(IPROP) == 1) PNUM(4) = IPROP
    if (ICOMP(IPROP) == 4) PNUM(5) = IPROP
    if (ICOMP(IPROP) == 6) PNUM(6) = IPROP
  end if
  if (PNAME(IPROP) == 'EF0    1') then
    if (ICOMP(IPROP) == 1) PNUM(7) = IPROP
  end if
end do
do IPROP=1,7
  if (PNUM(IPROP) == 0) then
    write(u6,*) 'DQVDiabat: Some required properties are not available'
    call Abend()
  end if
end do

! CEH Now we maximize f_DQV through a Jacobi sweep algorithm

! Initialize rotation matrix to unit matrix (theta=0)
call unitmat(TROT,NSTATE)

!Now, we have to form the trace of the quadrupole
!The user has to specify dipole then quadrupole
do ISTA=1,NSTATE
  do JSTA=1,NSTATE
    TRQ(ISTA,JSTA) = PROP(ISTA,JSTA,PNUM(4))+PROP(ISTA,JSTA,PNUM(5))+PROP(ISTA,JSTA,PNUM(6))
  end do
end do

call RecPrt('The TRQ matrix in DQVDiabat','',TRQ,NSTATE,NSTATE)
write(u6,*) ''

do i=1,itMAX
  THRSCH = Zero

  do ISTA=1,NSTATE
    do JSTA=ISTA+1,NSTATE

      ! Equations 9 and 10 of Kleir et al. with \mu written
      ! as a vector in addition to the quadrupole trace

      ATERM = PROP(ISTA,JSTA,PNUM(1))**2+PROP(ISTA,JSTA,PNUM(2))**2+PROP(ISTA,JSTA,PNUM(3))**2+ALPHZ*TRQ(ISTA,JSTA)**2+ &
              BETAE*PROP(ISTA,JSTA,PNUM(7))**2-Quart*(PROP(ISTA,ISTA,PNUM(1))**2+PROP(ISTA,ISTA,PNUM(2))**2+ &
              PROP(ISTA,ISTA,PNUM(3))**2+ALPHZ*TRQ(ISTA,ISTA)**2+BETAE*PROP(ISTA,ISTA,PNUM(7))**2)- &
              Quart*(PROP(JSTA,JSTA,PNUM(1))**2+PROP(JSTA,JSTA,PNUM(2))**2+PROP(JSTA,JSTA,PNUM(3))**2+ALPHZ*TRQ(JSTA,JSTA)**2+ &
              BETAE*PROP(JSTA,JSTA,PNUM(7))**2)+Half*(PROP(ISTA,ISTA,PNUM(1))*PROP(JSTA,JSTA,PNUM(1))+ &
              PROP(ISTA,ISTA,PNUM(2))*PROP(JSTA,JSTA,PNUM(2))+PROP(ISTA,ISTA,PNUM(3))*PROP(JSTA,JSTA,PNUM(3))+ &
              ALPHZ*TRQ(ISTA,ISTA)*TRQ(JSTA,JSTA)+BETAE*PROP(ISTA,ISTA,PNUM(7))*PROP(JSTA,JSTA,PNUM(7)))

      BTERM = PROP(ISTA,JSTA,PNUM(1))*(PROP(ISTA,ISTA,PNUM(1))-PROP(JSTA,JSTA,PNUM(1)))+ &
              PROP(ISTA,JSTA,PNUM(2))*(PROP(ISTA,ISTA,PNUM(2))-PROP(JSTA,JSTA,PNUM(2)))+ &
              PROP(ISTA,JSTA,PNUM(3))*(PROP(ISTA,ISTA,PNUM(3))-PROP(JSTA,JSTA,PNUM(3)))+ &
              ALPHZ*TRQ(ISTA,JSTA)*(TRQ(ISTA,ISTA)-TRQ(JSTA,JSTA))+BETAE*PROP(ISTA,JSTA,PNUM(7))* &
              (PROP(ISTA,ISTA,PNUM(7))-PROP(JSTA,JSTA,PNUM(7)))

      if (abs(BTERM) > MTF) then
        if ((ATERM > MTE) .or. (ATERM < -MTE)) then
          ! This part of the code computes the rotation angle.
          CTERM = -BTERM/ATERM
          ROTANGF = atan(CTERM)

          if ((ROTANGF < Zero) .and. (BTERM > Zero)) ROTANGF = ROTANGF+PI

          if ((ROTANGF > Zero) .and. (BTERM < Zero)) ROTANGF = ROTANGF-PI

          ROTANGO = ROTANGF*Quart

          if ((ROTANGO > PI*Quart) .or. (ROTANGO < -PI*Quart)) write(u6,*) 'Quadrants 2 and 3'

          COSO = cos(ROTANGO)
          SINO = sin(ROTANGO)

        else

          write(u6,*) 'A is pretty close to zero.  Something is fishy.'
          COSO = Zero
          SINO = Zero

        end if

        !Equation 13 of Kleir et al.
        CHNG = ATERM+sqrt(ATERM**2+BTERM**2)

        if (CHNG > THRSCH) THRSCH = CHNG

        if (CHNG > MTE) then

          !Update T rotation matrix
          do k=1,NSTATE
            T1 = TROT(k,ISTA)
            T2 = TROT(k,JSTA)
            TROT(k,ISTA) = COSO*T1+SINO*T2
            TROT(k,JSTA) = -SINO*T1+COSO*T2
          end do

          ! Rotates transition dipole matrices
          ! These equations can be found in a book on Jacobi sweeps

          TII = PROP(ISTA,ISTA,PNUM(1))
          TJJ = PROP(JSTA,JSTA,PNUM(1))
          TIJ = PROP(ISTA,JSTA,PNUM(1))
          do k=1,NSTATE
            T1 = PROP(ISTA,k,PNUM(1))
            T2 = PROP(JSTA,k,PNUM(1))
            PROP(ISTA,k,PNUM(1)) = COSO*T1+SINO*T2
            PROP(k,ISTA,PNUM(1)) = PROP(ISTA,k,PNUM(1))
            PROP(JSTA,k,PNUM(1)) = COSO*T2-SINO*T1
            PROP(k,JSTA,PNUM(1)) = PROP(JSTA,k,PNUM(1))
          end do
          PROP(ISTA,JSTA,PNUM(1)) = (COSO**2-SINO**2)*TIJ+COSO*SINO*(TJJ-TII)
          PROP(JSTA,ISTA,PNUM(1)) = PROP(ISTA,JSTA,PNUM(1))
          PROP(ISTA,ISTA,PNUM(1)) = COSO**2*TII+SINO**2*TJJ+Two*COSO*SINO*TIJ
          PROP(JSTA,JSTA,PNUM(1)) = SINO**2*TII+COSO**2*TJJ-Two*COSO*SINO*TIJ

          TII = PROP(ISTA,ISTA,PNUM(2))
          TJJ = PROP(JSTA,JSTA,PNUM(2))
          TIJ = PROP(ISTA,JSTA,PNUM(2))
          do k=1,NSTATE
            T1 = PROP(ISTA,k,PNUM(2))
            T2 = PROP(JSTA,k,PNUM(2))
            PROP(ISTA,k,PNUM(2)) = COSO*T1+SINO*T2
            PROP(k,ISTA,PNUM(2)) = PROP(ISTA,k,PNUM(2))
            PROP(JSTA,k,PNUM(2)) = COSO*T2-SINO*T1
            PROP(k,JSTA,PNUM(2)) = PROP(JSTA,k,PNUM(2))
          end do
          PROP(ISTA,JSTA,PNUM(2)) = (COSO**2-SINO**2)*TIJ+COSO*SINO*(TJJ-TII)
          PROP(JSTA,ISTA,PNUM(2)) = PROP(ISTA,JSTA,PNUM(2))
          PROP(ISTA,ISTA,PNUM(2)) = COSO**2*TII+SINO**2*TJJ+Two*COSO*SINO*TIJ
          PROP(JSTA,JSTA,PNUM(2)) = SINO**2*TII+COSO**2*TJJ-Two*COSO*SINO*TIJ

          TII = PROP(ISTA,ISTA,PNUM(3))
          TJJ = PROP(JSTA,JSTA,PNUM(3))
          TIJ = PROP(ISTA,JSTA,PNUM(3))
          do k=1,NSTATE
            T1 = PROP(ISTA,k,PNUM(3))
            T2 = PROP(JSTA,k,PNUM(3))
            PROP(ISTA,k,PNUM(3)) = COSO*T1+SINO*T2
            PROP(k,ISTA,PNUM(3)) = PROP(ISTA,k,PNUM(3))
            PROP(JSTA,k,PNUM(3)) = COSO*T2-SINO*T1
            PROP(k,JSTA,PNUM(3)) = PROP(JSTA,k,PNUM(3))
          end do
          PROP(ISTA,JSTA,PNUM(3)) = (COSO**2-SINO**2)*TIJ+COSO*SINO*(TJJ-TII)
          PROP(JSTA,ISTA,PNUM(3)) = PROP(ISTA,JSTA,PNUM(3))
          PROP(ISTA,ISTA,PNUM(3)) = COSO**2*TII+SINO**2*TJJ+Two*COSO*SINO*TIJ
          PROP(JSTA,JSTA,PNUM(3)) = SINO**2*TII+COSO**2*TJJ-Two*COSO*SINO*TIJ

          ! Rotates trace of quadrupole matrix

          TII = TRQ(ISTA,ISTA)
          TJJ = TRQ(JSTA,JSTA)
          TIJ = TRQ(ISTA,JSTA)
          do k=1,NSTATE
            T1 = TRQ(ISTA,k)
            T2 = TRQ(JSTA,k)
            TRQ(ISTA,k) = COSO*T1+SINO*T2
            TRQ(k,ISTA) = TRQ(ISTA,k)
            TRQ(JSTA,k) = COSO*T2-SINO*T1
            TRQ(k,JSTA) = TRQ(JSTA,k)
          end do
          TRQ(ISTA,JSTA) = (COSO**2-SINO**2)*TIJ+COSO*SINO*(TJJ-TII)
          TRQ(JSTA,ISTA) = TRQ(ISTA,JSTA)
          TRQ(ISTA,ISTA) = COSO**2*TII+SINO**2*TJJ+2*COSO*SINO*TIJ
          TRQ(JSTA,JSTA) = SINO**2*TII+COSO**2*TJJ-2*COSO*SINO*TIJ

          ! Rotates electric potential
          TII = PROP(ISTA,ISTA,PNUM(7))
          TJJ = PROP(JSTA,JSTA,PNUM(7))
          TIJ = PROP(ISTA,JSTA,PNUM(7))
          do k=1,NSTATE
            T1 = PROP(ISTA,k,PNUM(7))
            T2 = PROP(JSTA,k,PNUM(7))
            PROP(ISTA,k,PNUM(7)) = COSO*T1+SINO*T2
            PROP(k,ISTA,PNUM(7)) = PROP(ISTA,k,PNUM(7))
            PROP(JSTA,k,PNUM(7)) = COSO*T2-SINO*T1
            PROP(k,JSTA,PNUM(7)) = PROP(JSTA,k,PNUM(7))
          end do
          PROP(ISTA,JSTA,PNUM(7)) = (COSO**2-SINO**2)*TIJ+COSO*SINO*(TJJ-TII)
          PROP(JSTA,ISTA,PNUM(7)) = PROP(ISTA,JSTA,PNUM(7))
          PROP(ISTA,ISTA,PNUM(7)) = COSO**2*TII+SINO**2*TJJ+Two*COSO*SINO*TIJ
          PROP(JSTA,JSTA,PNUM(7)) = SINO**2*TII+COSO**2*TJJ-Two*COSO*SINO*TIJ

        end if
      end if
    end do
  end do

  if (THRSCH < THRS) then
    write(u6,*) 'Converged in this many iterations:'
    write(u6,*) i
    write(u6,*) ''
    exit
  end if

end do
if (i > itMAX) write(u6,*) 'Max iterations'

call RecPrt('Diabatic Coefficients','',TROT,NSTATE,NSTATE)

do ISTA=1,NSTATE
  do JSTA=1,NSTATE
    TROTT(JSTA,ISTA) = TROT(JSTA,ISTA)**2
  end do
end do
call RecPrt('Weights of adiabatic states','',TROTT,NSTATE,NSTATE)

do ISTA=1,NSTATE
  do JSTA=1,NSTATE
    HAMT(ISTA,JSTA) = HAM(ISTA,JSTA)
  end do
end do

TROTT = TROT

call DGEMM_('T','N',NSTATE,NSTATE,NSTATE,One,TROTT,NSTATE,HAMT,NSTATE,Zero,HDIAI,NSTATE)

call DGEMM_('N','N',NSTATE,NSTATE,NSTATE,One,HDIAI,NSTATE,TROTT,NSTATE,Zero,HDIA,NSTATE)

write(u6,*) ''
write(u6,*) 'The eigenvectors may no longer be in'
write(u6,*) 'energetic order'
write(u6,*) ''
call RecPrt('Diabatic Hamiltonian','',HDIA,NSTATE,NSTATE)

call CollapseOutput(0,'DQV Diabatization section')
write(u6,*) ''

! Molcas verify calls
call Add_Info('Adiabat1',HAM(1,1),1,4)
call Add_Info('Adiabat2',HAM(2,2),1,4)
call Add_Info('Adiabat3',HAM(3,3),1,4)
call Add_Info('DQVHam11',HDIA(1,1),1,4)
call Add_Info('DQVHam12',HDIA(1,2),1,4)
call Add_Info('DQVHam13',HDIA(1,3),1,4)
call Add_Info('DQVHam22',HDIA(2,2),1,4)
call Add_Info('DQVHam23',HDIA(2,3),1,4)
call Add_Info('DQVHam33',HDIA(3,3),1,4)
! End Molcas verify calls

end subroutine DQVDiabat
