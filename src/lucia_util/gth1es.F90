!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

!#define _DEBUGPRINT_
function GTH1ES(IREOTS,IPNT,H,IBSO,MXPNGAS,IBTSOB,NACOBS,IORB,ITP,ISM,JORB,JTP,JSM,IJSM)
!
! correct combination of row and column symmetry is assumed
! IJSM = 1 => Lower triangular packed
!      else=> No triangular packing
!
! Last Revision January 98 (IJSM added)

use Index_Functions, only: iTri
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
real(kind=wp) :: GTH1ES
integer(kind=iwp), intent(in) :: IREOTS(*), IPNT(*), IBSO(*), MXPNGAS, IBTSOB(MXPNGAS,*), NACOBS(*), IORB, ITP, ISM, JORB, JTP, &
                                 JSM, IJSM
real(kind=wp), intent(in) :: H(*)
integer(kind=iwp) :: I_ABS, IJ, IREO, JABS, JREO, NI, NJ

I_ABS = IORB+IBTSOB(ITP,ISM)-1
IREO = IREOTS(I_ABS)
JABS = JORB+IBTSOB(JTP,JSM)-1
JREO = IREOTS(JABS)
!write(u6,*) ' GTH1ES : IREO JREO ',IREO,JREO

!write(u6,*) ' GTH1ES : IBSO ',IBSO(ISM)
IJ = -2**30
if (IJSM == 1) then
  if (ISM > JSM) then
    NI = NACOBS(ISM)
    IJ = IPNT(ISM)-1+(JREO-IBSO(JSM))*NI+IREO-IBSO(ISM)+1
  else if (ISM == JSM) then
    IJ = IPNT(ISM)-1+iTri(IREO-IBSO(ISM)+1,JREO-IBSO(JSM)+1)
  else if (ISM < JSM) then
    NJ = NACOBS(JSM)
    IJ = IPNT(JSM)-1+(IREO-IBSO(ISM))*NJ+JREO-IBSO(JSM)+1
  end if
else
  NI = NACOBS(ISM)
  IJ = IPNT(ISM)-1+(JREO-IBSO(JSM))*NI+IREO-IBSO(ISM)+1
end if

GTH1ES = H(IJ)
#ifdef _DEBUGPRINT_
write(u6,*) ' one electron integral'
write(u6,*) ' IORB ITP ISM ',IORB,ITP,ISM
write(u6,*) ' JORB JTP JSM ',JORB,JTP,JSM
write(u6,*) ' IJ and H(IJ) ',IJ,H(IJ)
#endif

end function GTH1ES
