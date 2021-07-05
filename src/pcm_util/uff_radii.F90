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

function UFF_radii(IA)

implicit real*8(A-H,O-Z)
parameter(MaxIA=104)
dimension Rii(0:MaxIA)
data Pt5/0.5d0/
! Warning: the following are atomic diameters.
data Rii/0.000D+00, &
         !  H        He
         2.886D+00,2.362D+00, &
         ! Li        Be         B         C         N         O         F        Ne
         2.451D+00,2.745D+00,4.083D+00,3.851D+00,3.660D+00,3.500D+00,3.364D+00,3.243D+00, &
         ! Na        Mg        Al        Si         P         S        Cl        Ar
         2.983D+00,3.021D+00,4.499D+00,4.295D+00,4.147D+00,4.035D+00,3.947D+00,3.868D+00, &
         !  K        Ca
         3.812D+00,3.399D+00, &
                   ! Sc        Ti         V        Cr        Mn        Fe        Co        Ni        Cu        Zn
                   3.295D+00,3.175D+00,3.144D+00,3.023D+00,2.961D+00,2.912D+00,2.872D+00,2.834D+00,3.495D+00,2.763D+00, &
                             ! Ga        Ge        As        Se        Br        Kr
                             4.383D+00,4.280D+00,4.230D+00,4.205D+00,4.189D+00,4.141D+00, &
         ! Rb        Sr (+2)
         4.114D+00,3.641D+00, &
                   !  Y  (+3)  Zr (+4)   Nb (+5)   Mo (+6)   Tc (+5)   Ru (+2)   Rh (+3)   Pd (+2)   Ag (+1)   Cd (+2)
                   3.345D+00,3.124D+00,3.165D+00,3.052D+00,2.998D+00,2.963D+00,2.929D+00,2.899D+00,3.148D+00,2.848D+00, &
                             !       In            Sn           Sb              Te          I       Xe
                             4.463D+00,4.392D+00,4.420D+00,4.470D+00,4.50D+00,4.404D+00, &
         ! Cs        Ba (+2)
         4.517D+00,3.703D+00, &
                   ! La (+3)   Ce (+3)   Pr (+3)   Nd (+3)   Pm (+3)   Sm (+3)   Eu (+3)
                   ! Gd (+3)   Tb (+3)   Dy (+3)   Ho (+3)   Er (+3)   Tm (+3)   Yb (+3)
                   3.522D+00,3.556D+00,3.606D+00,3.575D+00,3.547D+00,3.520D+00,3.493D+00, &
                   3.368D+00,3.451D+00,3.428D+00,3.409D+00,3.391D+00,3.374D+00,3.355D+00, &
                   ! Lu (+3)   Hf (+4)   Ta (+5)    W(+4,+6) Re(+5,+7) Os (+6)   Ir (+3)   Pt        Au        Hg
                   3.640D+00,3.141D+00,3.170D+00,3.069D+00,2.954D+00,3.120D+00,2.840D+00,2.754D+00,3.293D+00,2.705D+00, &
                             ! Tl        Pb        Bi (+3)   Po (+2)   At        Rn (+4)
                             4.347D+00,4.297D+00,4.370D+00,4.709D+00,4.750D+00,4.765D+00, &
         ! Fr        Ra (+2)
         4.900D-01,3.677D+00, &
                   ! Ac (+3)   Th (+4)   Pa (+4)    U (+4)   Np (+4)   Pu (+4)   Am (+4)
                   ! Cm (+3)   Bk (+3)   Cf (+3)   Es (+3)   Fm (+3)   Md (+3)   No (+3)
                   3.478D+00,3.396D+00,3.424D+00,3.395D+00,3.424D+00,3.424D+00,3.381D+00, &
                   3.326D+00,3.339D+00,3.313D+00,3.299D+00,3.286D+00,3.274D+00,3.248D+00, &
                   ! Lw(+3)
                   3.236D+00,3.500D+00/

UFF_radii = Rii(min(max(IA,0),MaxIA))*Pt5

return

end function UFF_radii
