************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SUBROUTINE RMVERT(NLEV,NVERT,IDRT,IDOWN,NLIM,VER)
C Purpose: Remove vertices from a DRT table.
      use stdalloc, only: mma_allocate, mma_deallocate
      IMPLICIT None
      Integer :: NLEV, NVERT
      Integer :: IDRT(NVERT,5),IDOWN(NVERT,0:3)
      Integer :: VER(NVERT)

      Integer, PARAMETER :: LTAB=1,NTAB=2
      Integer ::  NLIM(NLEV)
      Integer, Allocatable:: CONN(:)
      Logical Test
      Integer IV, L, N, NCHANGES, IC, ID, NLD, NV

      Call mma_allocate(CONN,nVert,Label='CONN')

C KILL VERTICES THAT DO NOT OBEY RESTRICTIONS.
      DO IV=1,NVERT-1
        VER(IV)=1
        L=IDRT(IV,LTAB)
        N=IDRT(IV,NTAB)
        IF(N.LT.NLIM(L)) VER(IV)=0
      END DO
      VER(NVERT)=1

      NCHANGES=1 ! Initiate first loop
      Do While (NCHANGES>0)
C REMOVE ARCS HAVING A DEAD UPPER OR LOWER VERTEX.
C COUNT THE NUMBER OF ARCS REMOVED OR VERTICES KILLED.
      NCHANGES=0
      DO IV=1,NVERT-1
        IF(VER(IV).EQ.0) THEN
          DO IC=0,3
            ID=IDOWN(IV,IC)
            IF(ID.GT.0) THEN
              IDOWN(IV,IC)=0
              NCHANGES=NCHANGES+1
            END IF
          END DO
        ELSE
          NLD=0
          DO IC=0,3
            ID=IDOWN(IV,IC)
            IF(ID.GT.0) THEN
              IF(VER(ID).EQ.0) THEN
                IDOWN(IV,IC)=0
                NCHANGES=NCHANGES+1
              ELSE
                NLD=NLD+1
              END IF
            END IF
          END DO
          IF(NLD.EQ.0) THEN
            VER(IV)=0
            NCHANGES=NCHANGES+1
          END IF
        END IF
      END DO
C ALSO CHECK ON CONNECTIONS FROM ABOVE:
      CONN(:)=0
      CONN(1)=VER(1)
      DO IV=1,NVERT-1
        IF(VER(IV).EQ.1) THEN
          DO IC=0,3
            ID=IDOWN(IV,IC)
            Test = ID.GT.0
            If (Test) Then
               Test = VER(ID).EQ.1
            End If
            IF(Test) THEN
                CONN(ID)=1
            END IF
          END DO
        END IF
      END DO
      DO IV=1,NVERT
        IF(VER(IV).EQ.1 .AND. CONN(IV).EQ.0) THEN
          VER(IV)=0
          NCHANGES=NCHANGES+1
        END IF
      END DO

      End Do

      Call mma_deallocate(CONN)

C IF NO CHANGES, THE REMAINING GRAPH IS VALID.
C EVERY VERTEX OBEYS THE RESTRICTIONS. EVERY VERTEX IS
C CONNECTED ABOVE AND BELOW (EXCEPTING THE TOP AND BOTTOM)
C TO OTHER CONFORMING VERTICES.
C THE PROCEDURE IS GUARANTEED TO FIND A STABLE SOLUTIONS,
C SINCE EACH ITERATION REMOVES ARCS AND/OR VERTICES FROM THE
C FINITE NUMBER WE STARTED WITH.
      IF(VER(1).EQ.0) THEN
        WRITE(6,*)'RASSI/RMVERT: Too severe restrictions.'
        WRITE(6,*)'Not one single configuration is left.'
        CALL ABEND()
      END IF

      NV=0
      DO IV=1,NVERT
        IF(VER(IV).EQ.1) THEN
          NV=NV+1
          VER(IV)=NV
        END IF
      END DO
      NVERT=NV

      END
