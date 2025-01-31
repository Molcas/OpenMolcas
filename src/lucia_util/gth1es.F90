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

function GTH1ES(IREOTS,IPNT,H,IBSO,MXPNGAS,IBTSOB,NACOBS,IORB,ITP,ISM,JORB,JTP,JSM,IJSM)
!
! correct combination of row and column symmetry is assumed
! IJSM = 1 => Lower triangular packed
!      else=> No triangular packing
!
! Last Revision January 98 (IJSM added )

implicit real*8(A-H,O-Z)
! Input
integer IREOTS(*), IPNT(*), IBTSOB(MXPNGAS,*), IBSO(*)
integer NACOBS(*)
dimension H(*)

IABS = IORB+IBTSOB(ITP,ISM)-1
IREO = IREOTS(IABS)
JABS = JORB+IBTSOB(JTP,JSM)-1
JREO = IREOTS(JABS)
!write(6,*) ' GTH1ES : IREO JREO ',IREO,JREO

!write(6,*) ' GTH1ES : IBSO ',IBSO(ISM)
IJ = -2**30
if (IJSM == 1) then
  if (ISM > JSM) then
    NI = NACOBS(ISM)
    IJ = IPNT(ISM)-1+(JREO-IBSO(JSM))*NI+IREO-IBSO(ISM)+1
  else if (ISM == JSM) then
    IJMAX = max(IREO-IBSO(ISM)+1,JREO-IBSO(JSM)+1)
    IJMIN = min(IREO-IBSO(ISM)+1,JREO-IBSO(JSM)+1)
    IJ = IPNT(ISM)-1+IJMAX*(IJMAX-1)/2+IJMIN
  else if (ISM < JSM) then
    NJ = NACOBS(JSM)
    IJ = IPNT(JSM)-1+(IREO-IBSO(ISM))*NJ+JREO-IBSO(JSM)+1
  end if
else
  NI = NACOBS(ISM)
  IJ = IPNT(ISM)-1+(JREO-IBSO(JSM))*NI+IREO-IBSO(ISM)+1
end if

GTH1ES = H(IJ)
NTEST = 0
if (NTEST /= 0) then
  write(6,*) ' one electron integral'
  write(6,*) ' IORB ITP ISM ',IORB,ITP,ISM
  write(6,*) ' JORB JTP JSM ',JORB,JTP,JSM
  write(6,*) ' IJ and H(IJ) ',IJ,H(IJ)
end if

end function GTH1ES
