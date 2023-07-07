# ATLAS_HH_bbgg
## MadGraph5-Pythia8
F
llowing tutorial is assume using lxplus account. Since result files of these simulation are large, it is preferred to install the softwares in your CERN EOS storage space.
### MG5 installation
```bash
cd /eos/useri/<initial>/<user name>/
wget https://launchpad.net/mg5amcnlo/3.0/3.5.x/+download/MG5_aMC_v3.5.0.tar.gz
tar xvf MG5_aMC_v3.5.0.tar.gz
export mg5dir=$PWD/MG5_aMC_v3_5_0 
```
### pythia8 installation inside MG5
```bash
cd $mg5dir 
export PYTHIA8DATA=$PWD/HEPTools/pythia8/share/Pythia8/xmldoc/
./bin/mg5_aMC
```
Then, in MadGraph:
```MG5
install pythia8
```
During the process, it might ask you whether you want install other packages in support of pythia8, enter "y".

### MG simulation
```bash
cd $mg5dir
source /cvmfs/sft.cern.ch/lcg/views/LCG_101/x86_64-centos7-gcc11-opt/setup.sh 
export PYTHIA8DATA=$PWD/HEPTools/pythia8/share/Pythia8/xmldoc/
./bin/mg5_aMC
```
Here is a example for generating the signal for diHiggs channel, similar for background or other channel
```MG5
set auto_convert_model T
import model heft
generate p p > h h, (h > b b~), (h > g g)
output sig_bbgg_14TeV
launch
```
Then you will see a menu for MG5, type 1 and hit enter to open up the pythia8 for format transformation to .hepmc, then hit enter again to the next step--setting parameters:
```MG5
set mh 125.0
set wh 0.004
set ebeam1 7000
set ebeam2 7000
set nevents 100000
set iseed 1823211
set etaj 2.5
set etajmin -2.5
set etab 2.5
set etabmin -2.5
set ptj 20
set ptb 20
```
The process could take about several seconds or several hours depending on the number of events required. After generate the events, don't forget to gunzip the .hepmc.tar.gz for preparation for delphes simulation:
```bash
cd <your output directory>/Events/run_01/
gunzip tag_1_pythia8_events.hepmc.gz
```
## Detector simulation on Delphes
### Installation of Delphes
```bash
git clone git://github.com/delphes/delphes.git Delphes
cd Delphes
export delphes=$PWD
```
You could also get it from their official website:
```bash
wget http://cp3.irmp.ucl.ac.be/downloads/Delphes-3.5.0.tar.gz
tar xvf Delphes-3.5.0.tar.gz 
cd Delphes-3.5.0
export delphes=$PWD
```
Then compile
```bash
cd $delphes
make -j 8
```
### Simulation
```bash
./DelphesHepMC2 cards/delphes_card_HLLHC.tcl <name for output root file>.root $mg5dir/<your output directory>/Events/run_01/tag_1_pythia8_events.hepmc
```



