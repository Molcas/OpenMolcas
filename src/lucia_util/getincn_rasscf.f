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
      SUBROUTINE GETINCN_RASSCF(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,
     &                IXCHNG,IKSM,JLSM,
     &                IPNT2,NSMOB,INH1,ICOUL)
      use wadr, only: tuvx
* interface to RASSCF common blocks
      IMPLICIT None
      INTEGER ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,IXCHNG,IKSM,JLSM,NSMOB,
     & ICOUL
*. For Jesper and openMP
      INTEGER IPNT2(*),INH1(*)
      Real*8 XINT(*)
*
      CALL GETINCN_RASSCFS(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,
     &                IXCHNG,IKSM,JLSM,TUVX,NSMOB,ICOUL)
c Avoid unused argument warnings
      IF (.FALSE.) THEN
        CALL Unused_integer_array(IPNT2)
        CALL Unused_integer_array(INH1)
      END IF
      END SUBROUTINE GETINCN_RASSCF
