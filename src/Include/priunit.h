/***********************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
***********************************************************************/
#if defined(__CVERSION__)
struct common_priunit {
  int lucmd, lupri, luerr, lustat, luw4, lupot, ninfo, nwarn, iprstatr;
};
extern struct common_priunit priunit_;
#else
!     FILE: priunit.h
      CHARACTER*80 SEPARATOR
      PARAMETER (SEPARATOR = '----------------------------------------' &
     &                     //'----------------------------------------')
      INTEGER LUCMD
      INTEGER LUPRI, LUERR, LUSTAT, LUW4, LUPOT, NINFO, NWARN, IPRSTAT
      COMMON /PRIUNIT/ LUCMD,                                           &
     &        LUPRI, LUERR, LUSTAT, LUW4, LUPOT, NINFO, NWARN, IPRSTAT
! ---  end of priunit.h ---
#endif
