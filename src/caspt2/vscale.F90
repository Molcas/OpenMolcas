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
! Copyright (C) 1994, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------*
! 1994  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*

#ifdef _MOLCAS_MPP_
      SUBROUTINE V_SCALE (EIG,SCA,V,nRows,NAS,LDV,NIN,COND)
      use caspt2_module, only: ThrShS
      use constants, only: Zero, One
      use definitions, only: iwp, wp, u6

      IMPLICIT None

      integer(kind=iwp), intent(in):: nRows, NAS, LDV, NIN
      real(kind=wp), intent(in):: EIG(NAS),SCA(NAS)
      real(kind=wp), Intent(inout):: V(LDV,*)
      real(kind=wp), Intent(out):: COND(NIN)

      integer(kind=iwp) jVEC, J, I, iVec
      real(kind=wp) EVAL, FACT, SZ

      jVEC=0
      DO J=1,NAS
        EVAL=EIG(J)
        IF(EVAL.GE.THRSHS) THEN
          jVEC=jVEC+1
          FACT=One/SQRT(EVAL)
          IF(jVEC.EQ.J) THEN
            CALL DSCAL_(nRows,FACT,V(1,J),1)
          ELSE
            CALL DYAX(nRows,FACT,V(1,J),1,V(1,jVEC),1)
          END IF
        END IF
      END DO
      IF (jVEC.NE.NIN) THEN
        WRITE(u6,*) 'V_SCALE: '//                                       &
     &      'inconsitency in linear dependence removal, ABORT'
        call AbEnd()
      END IF
! Addition, for the scaled symmetric ON.
      DO I=1,nRows
        CALL DSCAL_(NIN,SCA(I),V(I,1),LDV)
      END DO
! The condition number, after scaling, disregarding linear dep.
      IF(NIN.GE.2) THEN
        DO jVEC=1,NIN
          SZ=Zero
          DO iVEC=1,nRows
            SZ=SZ+V(iVEC,jVEC)**2
          END DO
          COND(jVEC)=SZ
        END DO
      END IF
      END SUBROUTINE V_SCALE

#elif defined (NAGFOR)
      ! Some compilers do not like empty files
      subroutine empty_vscale()
      end subroutine empty_vscale
#endif
