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
#ifndef _CHOINTS_
 (file)  CHODIAG       "$WorkDir/$Project".ChDiag          rw
 (file)  CHOMAP        "$WorkDir/$Project".ChMap           rw
 (file)  CHORST        "$WorkDir/$Project".ChRst           rw
 (file)  CHRED         "$WorkDir/$Project".ChRed           rw
 (file)  CHOSEL        "$WorkDir/$Project".ChSel           rw*
 (file)  CHVEC         "$WorkDir/$Project".ChVec           rw*
 (file)  CHFV          "$WorkDir/$Project".ChFV            rw*
 (file)  CHRSTL        "$WorkDir/$Project".ChRstL          rwp
 (file)  CHRDL         "$WorkDir/$Project".ChRdL           rwp
 (file)  CHVCL         "$WorkDir/$Project".ChVcL           rw*p
 (file)  VTMP          "$WorkDir/$Project."Vtmp            rw*p
#ifdef _LDF_
 (file)  LDFC          "$WorkDir/$Project".LDFC            rw*
 (file)  LDFAP         "$WorkDir/$Project".LDFAP           rw*
#endif
#define _CHOINTS_
#endif
