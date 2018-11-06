#!/bin/bash

# This is a small script to download QCMaquis and its components from the SCINE website
# with curl, i.e., without a graphical browser.

url=$1 # should be: https://scine.ethz.ch/download/ under normal circumstances (note the trailing slash!)
package=$2 # choose from: dmrg dmrg_interface_utils hdf5_qcm nevpt2
first_name="$3" # Your first name
last_name="$4" # Your last name
email="$5" # Your e-mail address

# Visit the download page in order to retrive a correct CSRF token
curl -s -L -X GET -c cookies.txt $url > /dev/null
csrf_token=`grep 'csrftoken' cookies.txt | awk '{print $7}'`

# Send a valid post request; note that is is also necessary to set the correct referer with "-e"
curl -s -L -X POST -b cookies.txt -e "$url" -F "first_name=${first_name}" -F "last_name=${last_name}" -F "email=${email}" -F 'eula_accepted=on' -F "package=qcmaquis_${package}" -F "csrfmiddlewaretoken=${csrf_token}" ${url} -o ${package}.tar.bz2 > /dev/null

if [ -f cookies.txt ]; then
  rm cookies.txt
fi

# Check the download integrity
bzip2 -t ${package}.tar.bz2
if [ $? -ne 0 ]; then
  echo "Download of ${package} failed." 1>&2
  exit 1
fi
