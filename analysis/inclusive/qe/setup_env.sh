echo "Loading CLAS12 environment and software for analysis"
echo " "

##source /group/clas12/packages/setup.csh ## Deprecated
module use /scigroup/cvmfs/hallb/clas12/sw/modulefiles

module purge
module load clas12

setenv HIPO /group/clas12/packages/hipo/dev


echo " "
echo "Setting up environment and software for CLAS12 analysis Complete"
