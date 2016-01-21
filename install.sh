#!/bin/bash

#If bash profile file doesn't exist create one
if ! [ -f ~/.profile ]
then
	touch ~/.profile
fi

#Delete any FAST installation lines
sed -i -e '/^#Adding FAST_v/{N;d;}' ~/.profile

#Write FAST installation lines
echo -e "\n#Adding FAST_v1.0/bin directory to PATH\nexport PATH=\"${PWD}/bin:\$PATH\"\n" >> ~/.profile
echo -e "\n#Adding FAST_v1.0 to PYTHONPATH\nexport PYTHONPATH=\"${PWD}:\$PYTHONPATH\"\n" >> ~/.profile

chmod +x bin/*