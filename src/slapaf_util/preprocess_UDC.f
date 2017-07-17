************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2014,2015, Ignacio Fdez. Galvan                        *
************************************************************************
*  Preprocess_UDC
*
*> @brief
*>   Quickly read the constraints to take some decisions
*> @author Ignacio Fdez. Galv&aacute;n
*>
*> @details
*> Process the constraints in the \p Lu file, but just to find out whether
*> some constraints are present and their values.
*> This is used for:
*>   - Detecting if the "EDIFF" constraint is present and with a zero value,
*>     in which case a conical intersection algorithm could be activated
*>     (depending on the symmetry and spin in the runfile).
*>   - Detecting if any constraint is explicitly declared as "soft",
*>     this is recommended for some TS searches with numerical differentiation.
*>   - Detecting whether an explicit MEP/IRC constraint has been included,
*>     this will probably give problems.
*>
*> @param[in]  Lu     Unit number of the file with the constraints
*> @param[in]  iPrint Print level
************************************************************************
      SUBROUTINE Preprocess_UDC(Lu,iPrint)
      IMPLICIT NONE
      INTEGER :: Lu,iPrint,iPos,Error,nLines,i,j
      CHARACTER(LEN=180) :: Line1,Line2,EDiffName,Get_Ln
      REAL*8 :: EDiffValue
#include "info_slapaf.fh"
#include "real.fh"
#include "nadc.fh"

      EDiffName = ''
      EDiffZero = .FALSE.
      iState(1) = 0
      iState(2) = 0
*
*     An arbitrary initial value given to EDiffValue, this value is
*     irrelevant since the only thing that matters is whether it is 0.0
*
      EDiffValue = One
*
*     Ugly hack to be able to restore the file to its original position
*     (since Get_Ln may read more than one line)
*     Count the lines until the end of the file and backspace
*
      nLines = 0
      DO
        READ(Lu,'(A)',IOSTAT=Error) Line1
        nLines = nLines+1
        IF (Error.ne.0) EXIT
      END DO
      DO i=1,nLines
        BACKSPACE(Lu)
      END DO
*
*     Read the primitive constraints
*
*     First read until the "VALUES" line
*     Each line is split at the "=" sign by auto,
*     so read two lines if there is no "=" sign
*
      DO
        Line1 = Get_Ln(Lu)
        CALL UpCase(Line1)
        CALL LeftAd(Line1)
        IF (Line1(1:4).eq.'VALU') EXIT
        iPos = INDEX(Line1,'=')
        IF (iPos.gt.0) THEN
          Line2 = Line1(iPos+1:)
          Line1 = Line1(:iPos-1)
        ELSE
          Line2 = Get_Ln(Lu)
          CALL UpCase(Line2)
        END IF
        CALL LeftAd(Line2)
*       If a primitive is defined as "EDIFF", save its name
        IF (Line2(1:4).eq.'EDIF') THEN
          EDiffName = Line1
          iPos = INDEX(Line2,' ')
          Line2 = Line2(iPos+1:)
          CALL LeftAd(Line2)
          READ(Line2,*,IOSTAT=Error) i
          IF (Error.ne.0) i=0
          iPos = INDEX(Line2,' ')
          Line2 = Line2(iPos+1:)
          CALL LeftAd(Line2)
          READ(Line2,*,IOSTAT=Error) j
          IF (Error.ne.0) j=0
          iState(1)=MAX(i,j)
          iState(2)=MIN(i,j)
        END IF
*       If a primitive of the same type as that used by MEP/IRC
*       is defined, signal it.
        IF (Line2(1:4).eq.MEP_Type(1:4)) THEN
          MEPCons=.TRUE.
        END IF
      END DO
*
*     Now read the constraint values
*
*     Read until the "END" line
*     Each line is split at the "=" sign by auto,
*     so read two lines if there is no "=" sign
*
      DO
        Line1 = Get_Ln(Lu)
        CALL UpCase(Line1)
        CALL LeftAd(Line1)
*       Read continuation lines
        DO WHILE (INDEX(Line1,'&').ne.0)
          Line1 = Get_Ln(Lu)
          CALL UpCase(Line1)
          CALL LeftAd(Line1)
        END DO
        IF (Line1(1:4).eq.'END ') EXIT
        iPos = INDEX(Line1,'=')
        IF (iPos.gt.0) THEN
          Line2 = Line1(iPos+1:)
          Line1 = Line1(:iPos-1)
        ELSE
          Line2 = Get_Ln(Lu)
          CALL UpCase(Line2)
        END IF
*       If the name matches the "EDIFF" primitive, read the value
        IF (Line1.eq.EDiffName) THEN
          READ(Line2,*,IOSTAT=Error) EDiffValue
          IF (Error.ne.0) EDiffValue = One
        END IF
*       If the constraint is explicitly "soft"
        IF (INDEX(Line2,'SOFT').ne.0) lSoft = .TRUE.
      END DO
*
*     Return the file to its original position
*     (cannot use REWIND because this may be the full input for slapaf)
*     First advance to the end of file and then backspace
*
      DO
        READ(Lu,'(A)',IOSTAT=Error) Line1
        IF (Error.ne.0) EXIT
      END DO
      DO i=1,nLines
        BACKSPACE(Lu)
      END DO
*
*     If an "EDIFF" constraint is being used with a value exactly 0.0,
*     this may be a conical intersection search.
*     If there is no "EDIFF" at all, there's no use computing NACs
*
      IF (EDiffValue.eq.Zero) THEN
        EDiffZero = .TRUE.
        IF (iPrint.ge.6) THEN
          WRITE(6,*) 'Energy difference constraint with zero value.'
          WRITE(6,*) 'This may be a conical intersection search.'
        END IF
      ELSE IF (EDiffName(1:4).ne.'    ') THEN
        IF (iPrint.ge.6) THEN
          WRITE(6,*) 'Energy difference constraint with non-zero value.'
          WRITE(6,*) 'This will not be a conical intersection search.'
        END IF
      ELSE
        NADC=.FALSE.
      END IF

      END
