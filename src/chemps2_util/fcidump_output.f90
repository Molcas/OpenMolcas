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
! Copyright (C) 2016, Sebastian Wouters                                *
!***********************************************************************
! Subroutine to write FCIDUMP file
! Written by Sebastian Wouters, Leuven, Aug 2016
SUBROUTINE FCIDUMP_OUTPUT( NACT, NELEC, TWOMS, ISYM, ORBSYM, ECONST, OEI, TEI, LINSIZE, NUM_TEI )

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NACT, NELEC, TWOMS, ISYM
  INTEGER, INTENT(IN) :: ORBSYM( NACT )
  REAL*8,  INTENT(IN) :: ECONST
  INTEGER, INTENT(IN) :: LINSIZE, NUM_TEI
  REAL*8,  INTENT(IN) :: OEI( LINSIZE )
  REAL*8,  INTENT(IN) :: TEI( NUM_TEI )
  INTEGER :: i, j, k, l, ij, kl, ijkl
  INTEGER :: writeout, isfreeunit

  EXTERNAL ISFREEUNIT


  writeout=isfreeunit(28)
!  open ( unit=writeout, file="FCIDUMP_CHEMPS2", action="write", status="replace" )
  call molcas_open(writeout,"FCIDUMP_CHEMPS2")
  write ( writeout, "(A11,I3,A7,I3,A5,I2,A1)" ) " &FCI NORB=", NACT, ",NELEC=", NELEC, ",MS2=", TWOMS, ","
  write ( writeout, "(A9)", ADVANCE = "NO" ) "  ORBSYM="
  do i=1,NACT
    write ( writeout, "(I1,A1)", ADVANCE = "NO" ) ORBSYM( i ), ","
  enddo
  write( writeout, * )
  write( writeout, "(A7,I1,A1)" ) "  ISYM=", ISYM, ","
  write( writeout, "(A2)" ) " /"

  do i=1,NACT
    do j=1,i
      ij = ( ( i - 1 ) * i ) / 2 + ( j - 1 )
      do k=1,NACT
        do l=1,k
          kl = ( ( k - 1 ) * k ) / 2 + ( l - 1 )
          if ( kl .LE. ij ) then
            ijkl = 1 + ( ij * ( ij + 1 ) ) / 2 + kl
            if ( ABS( TEI( ijkl ) ) .GE. 1.0e-16 ) then
              write( writeout, "(1X,ES23.16E2,I4,I4,I4,I4)" ) TEI( ijkl ), i, j, k, l
            end if
          end if
        enddo
      enddo
    enddo
  enddo

  do i=1,NACT
    do j=1,i
      ij = 1 + ( ( i - 1 ) * i ) / 2 + ( j - 1 )
      if ( ABS( OEI( ij ) ) .GE. 1.0e-16 ) then
         write( writeout, "(1X,ES23.16E2,I4,I4,I4,I4)" ) OEI( ij ), i, j, 0, 0
      end if
    enddo
  enddo

  write( writeout, "(1X,ES23.16E2,I4,I4,I4,I4)" ) ECONST, 0, 0, 0, 0

  close ( writeout )

END SUBROUTINE FCIDUMP_OUTPUT
