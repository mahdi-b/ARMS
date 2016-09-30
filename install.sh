
# python dependencies
pip install argparse
pip install enum34
pip install ete2

mkdir programs
cd programs

# spades
wget http://spades.bioinf.spbau.ru/release3.9.0/SPAdes-3.9.0-Linux.tar.gz
tar -xzf SPAdes-3.9.0-Linux.tar.gz
mv SPAdes-3.9.0-Linux spades
rm SPAdes-3.9.0-Linux.tar.gz


#pear
wget http://sco.h-its.org/exelixis/web/software/pear/files/pear-0.9.10-bin-64.tar.gz
tar -xzf pear-0.9.10-bin-64.tar.gz
cd pear-0.9.10-bin-64
mv pear-0.9.10-bin-64 pear
cd ..
mv pear-0.9.10-bin-64 pear
rm pear-0.9.10-bin-64.tar.gz


#fastx
wget http://hannonlab.cshl.edu/fastx_toolkit/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2
tar -xjf fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2
sudo cp ./bin/* /usr/local/bin
sudo ldconfig
rm fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2


#flexbar
sudo apt install flexbar


#vsearch
wget https://github.com/torognes/vsearch/archive/v1.11.1.tar.gz
tar -xzf v1.11.1.tar.gz
cd vsearch-1.11.1
./autogen.sh
./configure
make
sudo make install
mv bin/vsearch .
cd ..
mv vsearch-1.11.1 vsearch
rm v1.11.1.tar.gz

#crop
# Might need GSL developer library
sudo apt install libgsl-dev
wget https://github.com/tingchenlab/CROP/archive/master.zip
unzip master.zip
cd CROP-master
make
mv CROPLinux crop
cd ..
mv CROP-master crop
rm master.zip


#swarm
wget https://github.com/torognes/swarm/archive/master.zip
unzip master.zip
cd swarm-master/src
make
cd ../bin/
mv swarm ..
cd ../..
mv swarm-master swarm
rm master.zip

ls /usr/bin/fastx_barcode_splitter.pl
ls /usr/bin/flexbar
ls ~/ARMS/programs/pear/pear-0.9.5-bin-64
ls ~/ARMS/programs/spades/bin/spades.py
ls ~/ARMS/programs/vsearch/vsearch
ls ~/ARMS/programs/swarm/swarm-2.1.9-linux-x86_64
ls ~/ARMS/programs/crop/crop
