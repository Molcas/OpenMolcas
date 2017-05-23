#***********************************************************************
# This file is part of OpenMolcas.                                     *
#                                                                      *
# OpenMolcas is free software; you can redistribute it and/or modify   *
# it under the terms of the GNU Lesser General Public License, v. 2.1. *
# OpenMolcas is distributed in the hope that it will be useful, but it *
# is provided "as is" and without any express or implied warranties.   *
# For more details see the full text of the license in the file        *
# LICENSE or in <http://www.gnu.org/licenses/>.                        *
#***********************************************************************
# Version 0.1
# used to build a docker image with a Molcas installation in it
FROM ubuntu:14.04
MAINTAINER Steven Vancoillie "molcas@steven.se"
RUN apt-get update && apt-get install -y gcc g++ gfortran perl bash curl cmake
ADD . /opt/molcas/src
