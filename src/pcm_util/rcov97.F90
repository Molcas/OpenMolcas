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

function RCov97(IA,IB)
! This function returns an estimated covalent bond distance (in Ang)
! between atoms of atomic numbers IA and IB. Setting IB to 0 returns
! the covalent radius of IA. Parameters for atoms heavier than At are
! taken from the UFF force field (A.K.Rappe',C.J.Casewit,K.S.Colwell,
! W.A.Goddard III and W.M.Skiff, J.Am.Chem.Soc. 114,10024 (1992)).

implicit real*8(A-H,O-Z)
parameter(MxAtN=104)
dimension Rii(0:MxAtN)
save Rii
data Rii/0.0d0, &
         0.354d0,0.849d0, &
         1.336d0,1.010d0,0.838d0,0.757d0,0.700d0,0.658d0,0.668d0,0.920d0, &
         1.539d0,1.421d0,1.244d0,1.117d0,1.101d0,1.064d0,1.044d0,1.032d0, &
         1.953d0,1.761d0, &
                 1.513d0,1.412d0,1.402d0,1.345d0,1.382d0,1.270d0,1.241d0,1.164d0,1.302d0,1.193d0, &
                         1.260d0,1.197d0,1.211d0,1.190d0,1.192d0,1.147d0, &
         2.260d0,2.052d0, &
                 1.698d0,1.564d0,1.473d0,1.467d0,1.322d0,1.478d0,1.332d0,1.338d0,1.386d0,1.403d0, &
                         1.459d0,1.398d0,1.407d0,1.386d0,1.382d0,1.267d0, &
         2.570d0,2.277d0, &
                 1.943d0,1.841d0,1.823d0,1.816d0,1.801d0,1.780d0,1.771d0,1.735d0,1.732d0,1.710d0,1.696d0,1.673d0,1.660d0,1.637d0, &
                 1.671d0,1.611d0,1.511d0,1.392d0,1.372d0,1.372d0,1.371d0,1.364d0,1.262d0,1.340d0, &
                         1.518d0,1.459d0,1.512d0,1.500d0,1.545d0,1.420d0, &
         2.880d0,2.512d0, &
                 1.983d0,1.721d0,1.711d0,1.684d0,1.666d0,1.657d0,1.660d0,1.801d0,1.761d0,1.750d0,1.724d0,1.712d0,1.689d0,1.679d0, &
                 1.698d0,1.850d0/

RCov97 = Rii(min(max(IA,0),MxAtN))+Rii(min(max(IB,0),MxAtN))

return

end function RCov97
