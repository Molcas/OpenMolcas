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

subroutine Freeze_Default(iANr,nShell,lMax)

integer nShell(0:lMax)
parameter(nAtoms=106)
integer iDefaults(0:3,0:nAtoms)
!---- First row, H-He: no frozen orbitals
data((iDefaults(i,j),i=0,3),j=0,2)/0,0,0,0, &
                                   0,0,0,0, &
                                   0,0,0,0/
!---- Second row, Li-Ne: 1s
data((iDefaults(i,j),i=0,3),j=3,10)/1,0,0,0, &
                                    1,0,0,0, &
                                    1,0,0,0, &
                                    1,0,0,0, &
                                    1,0,0,0, &
                                    1,0,0,0, &
                                    1,0,0,0, &
                                    1,0,0,0/
!---- Third row, Na-Si: 1s2s
data((iDefaults(i,j),i=0,3),j=11,14)/2,0,0,0, &
                                     2,0,0,0, &
                                     2,0,0,0, &
                                     2,0,0,0/
!---- Third row, P-Ar: 1s2s2p
data((iDefaults(i,j),i=0,3),j=15,18)/2,1,0,0, &
                                     2,1,0,0, &
                                     2,1,0,0, &
                                     2,1,0,0/
!---- Fourth row, K-Cu: 1s2s2p3s
data((iDefaults(i,j),i=0,3),j=19,30)/3,1,0,0, &
                                     3,1,0,0, &
                                     3,1,0,0, &
                                     3,1,0,0, &
                                     3,1,0,0, &
                                     3,1,0,0, &
                                     3,1,0,0, &
                                     3,1,0,0, &
                                     3,1,0,0, &
                                     3,1,0,0, &
                                     3,1,0,0, &
                                     3,1,0,0/
!---- Fourth row, Zn-Kr: 1s2s2p3s3p3d
data((iDefaults(i,j),i=0,3),j=31,36)/3,2,1,0, &
                                     3,2,1,0, &
                                     3,2,1,0, &
                                     3,2,1,0, &
                                     3,2,1,0, &
                                     3,2,1,0/
!---- Fifth row, Rb-Ag: 1s2s2p3s3p4s3d
data((iDefaults(i,j),i=0,3),j=37,47)/4,2,1,0, &
                                     4,2,1,0, &
                                     4,2,1,0, &
                                     4,2,1,0, &
                                     4,2,1,0, &
                                     4,2,1,0, &
                                     4,2,1,0, &
                                     4,2,1,0, &
                                     4,2,1,0, &
                                     4,2,1,0, &
                                     4,2,1,0/
!---- Fifth row, Cd-Xe: 1s2s2p3s3p4s3d4p
data((iDefaults(i,j),i=0,3),j=48,54)/4,3,1,0, &
                                     4,3,1,0, &
                                     4,3,1,0, &
                                     4,3,1,0, &
                                     4,3,1,0, &
                                     4,3,1,0, &
                                     4,3,1,0/
!---- Sixth row, Cs-Yb: 1s2s2p3s3p4s3d4p5s4d5p
data((iDefaults(i,j),i=0,3),j=55,70)/4,3,2,0, &
                                     4,3,2,0, &
                                     4,3,2,0, &
                                     4,3,2,0, &
                                     4,3,2,0, &
                                     4,3,2,0, &
                                     4,3,2,0, &
                                     4,3,2,0, &
                                     4,3,2,0, &
                                     4,3,2,0, &
                                     4,3,2,0, &
                                     4,3,2,0, &
                                     4,3,2,0, &
                                     4,3,2,0, &
                                     4,3,2,0, &
                                     4,3,2,0/
!---- Sixth row, Lu-Au: 1s2s2p3s3p4s3d4p4d4f
data((iDefaults(i,j),i=0,3),j=71,79)/4,3,2,1, &
                                     4,3,2,1, &
                                     4,3,2,1, &
                                     4,3,2,1, &
                                     4,3,2,1, &
                                     4,3,2,1, &
                                     4,3,2,1, &
                                     4,3,2,1, &
                                     4,3,2,1/
!---- Sixth row, Hg-Rn: 1s2s2p3s3p4s3d4p4d4f5s5p
data((iDefaults(i,j),i=0,3),j=80,86)/5,4,2,1, &
                                     5,4,2,1, &
                                     5,4,2,1, &
                                     5,4,2,1, &
                                     5,4,2,1, &
                                     5,4,2,1, &
                                     5,4,2,1/
!---- Seventh row, Fr-Cm: 1s2s2p3s3p4s3d4p5s4d5p4f5d
data((iDefaults(i,j),i=0,3),j=87,96)/5,4,3,1, &
                                     5,4,3,1, &
                                     5,4,3,1, &
                                     5,4,3,1, &
                                     5,4,3,1, &
                                     5,4,3,1, &
                                     5,4,3,1, &
                                     5,4,3,1, &
                                     5,4,3,1, &
                                     5,4,3,1/
data((iDefaults(i,j),i=0,3),j=97,103)/6,5,3,1, &
                                      6,5,3,1, &
                                      6,5,3,1, &
                                      6,5,3,1, &
                                      6,5,3,1, &
                                      6,5,3,1, &
                                      6,5,3,1/
data((iDefaults(i,j),i=0,3),j=104,106)/6,5,3,2, &
                                       6,5,3,2, &
                                       6,5,3,2/

if ((iANr < 0) .or. (iANr > nAtoms)) then
  write(6,*) 'Freeze_Defaults: iAnr is out of range!'
  write(6,*) 'iANr=',iANr
  call Abend()
end if

call iCopy(lMax+1,[0],0,nShell,1)

do i=0,min(lMax,3)
  nShell(i) = iDefaults(i,iAnr)
end do

return

end subroutine Freeze_Default
