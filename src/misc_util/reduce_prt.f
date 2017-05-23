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
      LOGICAL FUNCTION Reduce_Prt()
*
      IMPLICIT NONE
      INTEGER :: i,Err
      CHARACTER(LEN=80) :: Word
      CHARACTER(LEN=100) :: SuperName
      CHARACTER(LEN=100), EXTERNAL :: Get_SuperName, Get_ProgName
*
      Reduce_Prt = .FALSE.
*
*     Do not reduce printing in last_energy
*
      SuperName = Get_SuperName()
      IF (SuperName .eq. 'last_energy') RETURN
*
*     Reduce printing if iter > 1
*
      CALL GetEnvF("MOLCAS_ITER",Word)
      READ (Word,*) i
      IF (i .gt. 1) Reduce_Prt = .TRUE.
*
*     ... but not if MOLCAS_REDUCE_PRT = NO
*
      IF (Reduce_Prt) THEN
        CALL GetEnvF("MOLCAS_REDUCE_PRT",Word)
        IF (Word(1:1) .eq. 'N') Reduce_Prt = .FALSE.
      END IF
*
*     ... or if we are not inside a loop (EMIL_InLoop < 1)
*
      IF (Reduce_Prt) THEN
        CALL GetEnvF("EMIL_InLoop",Word)
        i = 0
        READ (Word,*,IOSTAT=Err) i
        IF (i .lt. 1) Reduce_Prt = .FALSE.
      END IF
*
*     ... or if first iteration of a saddle branch (SADDLE_FIRST = 1)
*
      IF (Reduce_Prt) THEN
        CALL GetEnvF("SADDLE_FIRST",Word)
        i = 0
        READ (Word,*,IOSTAT=Err) i
        IF (i .eq. 1) Reduce_Prt = .FALSE.
      END IF
*
*     In any case, reduce printing inside numerical gradients,
*     unless specified otherwise (MOLCAS_REDUCE_NG_PRT = NO).
*
      IF (.NOT. Reduce_Prt) THEN
        IF ((SuperName .eq. 'numerical_gradient') .AND.
     &      (Get_ProgName() .ne. 'numerical_gradient')) THEN
          CALL GetEnvF("MOLCAS_REDUCE_NG_PRT",Word)
          IF (Word(1:1) .ne. 'N') Reduce_Prt = .TRUE.
        END IF
      END IF
*
      RETURN
*
      END FUNCTION
