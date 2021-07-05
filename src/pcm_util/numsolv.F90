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

integer function NumSolv(Solvent)

implicit real*8(a-h,o-z)
character*30 Solvent, Solvent_

IdSolv = 0
Solvent_ = Solvent
call Upcase(Solvent_)
if (Solvent_ == 'WATER') IdSolv = 1
if (Solvent_ == 'ACETONITRILE') IdSolv = 2
if (Solvent_ == 'METHANOL') IdSolv = 3
if (Solvent_ == 'ETHANOL') IdSolv = 4
if (Solvent_ == 'ISOQUINOLINE') IdSolv = 5
if (Solvent_ == 'QUINOLINE') IdSolv = 6
if (Solvent_ == 'CHLOROFORM') IdSolv = 7
if (Solvent_ == 'ETHYLETHER') IdSolv = 8
if (Solvent_ == 'METHYLENECHLORIDE') IdSolv = 9
if (Solvent_ == 'DICHLOROETHANE') IdSolv = 10
if (Solvent_ == 'CARBONTETRACHLORIDE') IdSolv = 11
if (Solvent_ == 'BENZENE') IdSolv = 12
if (Solvent_ == 'TOLUENE') IdSolv = 13
if (Solvent_ == 'CHLOROBENZENE') IdSolv = 14
if (Solvent_ == 'NITROMETHANE') IdSolv = 15
if (Solvent_ == 'HEPTANE') IdSolv = 16
if (Solvent_ == 'CYCLOHEXANE') IdSolv = 17
if (Solvent_ == 'ANILINE') IdSolv = 18
if (Solvent_ == 'ACETONE') IdSolv = 19
if (Solvent_ == 'TETRAHYDROFURAN') IdSolv = 20
if (Solvent_ == 'DIMETHYLSULFOXIDE') IdSolv = 21
if (Solvent_ == 'ARGON') IdSolv = 22
if (Solvent_ == 'KRYPTON') IdSolv = 23
if (Solvent_ == 'XENON') IdSolv = 24

if (IdSolv == 0) then
  !call ErrTra()
  write(6,*) ' Unrecognized solvent: ',Solvent
  write(6,10) 'WATER','ACETONITRILE','METHANOL','ETHANOL','ISOQUINOLINE','QUINOLINE','CHLOROFORM','ETHYLETHER', &
              'METHYLENECHLORIDE','DICHLOROETHANE','CARBONTETRACHLORIDE','BENZENE','TOLUENE','CHLOROBENZENE','NITROMETHANE', &
              'HEPTANE','CYCLOHEXANE','ANILINE','ACETONE','TETRAHYDROFURAN','DIMETHYLSULFOXIDE','ARGON','KRYPTON','XENON'
  call Abend()
end if
NumSolv = IdSolv

return

10 format(' Allowed solvents are:',/,24(A30,/))

end function NumSolv
