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
* Copyright (C) 1990, Jeppe Olsen                                      *
************************************************************************
      SUBROUTINE ICISPC(MNRS10,MXRS30,IPRNT)
*
* Obtain internal CI spaces relevant for MRSDCI
*       /STRINP/+/LUCINP/ = > /CICISP/
* Jeppe Olsen , Dec 1990
*
* Input
* =====
* Information in STRINP
*
* Output
* ======
* Common block CICISP
*
* Internal CI spaces
********************************************************************
*   *  Basic space  * Allowed internal excit * Delta NA * Delta Nb *
********************************************************************
* 1 *  Zero order   *           0            *    0     *    0     *
********************************************************************
* ====================
*. Input common blocks
* ====================
*./LUCINP : EXTSPC is used
#include "detdim.fh"
*./STRINP/
#include "strinp_mclr.fh"
*/ORBINP/
#include "orbinp_mclr.fh"
* ====================
*. Output common block
* ====================
#include "cicisp_mclr.fh"
* NICISP : Number of internal CI spaces constructed
* IASTFI : Alpha string type for internal CI space
* IBSTFI : Beta string type for internal CI space
* IACTI  : Given internal space is active
* MXR3IC : Max number of elecs in RAS 3 space for internal CI space
* MNR1IC : Min number of elecs in RAS 1 space for internal CI space
* IZCI   : Internal zero order space
* IRCI(IEX,DELTAA+5,DELTAB+5) : Number of zero order space
* NELCI : Number of electrons per CI space
* obtained by (IEX-1) fold internal excitation , with a NAEL + DELTAA
* alpha electrons and  NBEL + DELTAB beta electrons
*
*     Call qEnter('ICISPC')
*
      NTEST = 00000
      NTEST = MAX(NTEST,IPRNT)
*
      ICI = 1
      MNR1IC(ICI) = MNRS10
      MXR3IC(ICI) = MXRS30
      IASTFI(ICI) = IAZTP
      IBSTFI(ICI) = IBZTP
      NAELCI(ICI) = NELEC(IAZTP)
      NBELCI(ICI) = NELEC(IBZTP)
      NELCI(ICI)  = NAELCI(ICI)+NBELCI(ICI)
      IACTI(1) = 1
      NICISP = ICI
* EAW Just zero order
      Call iCopy(3*49,0,0,irci,1)
* EAW
*. Number and distribution of electrons in each space
      DO 100 IEX = 1, 3
      DO 100 IDA = -4,2
      DO 100 IDB = -4,2
        IF(IRCI(IEX,IDA+5,IDB+5).NE.0) THEN
           ICI = IRCI(IEX,IDA+5,IDB+5)
           NAELCI(ICI) = NELEC(IASTFI(ICI))
           NBELCI(ICI) = NELEC(IBSTFI(ICI))
           NELCI(ICI) = NAELCI(ICI)+NBELCI(ICI)
        END IF
100   CONTINUE
*
*. Default max in RAS1 and min in RAS3
      DO 150 ICI = 1, NICISP
        MXR1IC(ICI) = MIN(2*NORB1,NELCI(ICI))
        MNR3IC(ICI) = MAX(0,NELCI(ICI)-2*(NORB1+NORB2))
150   CONTINUE
*
      IF(NTEST .GE. 1 ) THEN
        WRITE(6,*) ' Number of internal CI spaces ', NICISP
        WRITE(6,*)
     &  ' Space a-type b-type nael nbel mnrs1 mxrs1 mnrs3 mxrs3 '
        WRITE(6,*)
     &  ' ===================================================== '
         DO 1020 ICI = 1, NICISP
          IF(IACTI(ICI).EQ.1)
     &    WRITE(6,'(I5,2I7,2I5,4I6)')
     &    ICI,IASTFI(ICI),IBSTFI(ICI),NAELCI(ICI),NBELCI(ICI),
     &    MNR1IC(ICI),MXR1IC(ICI),MNR3IC(ICI),MXR3IC(ICI)
1020    CONTINUE
      END IF
*
*     Call qExit('ICISPC')
*
      RETURN
      END
