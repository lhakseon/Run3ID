#!/bin/sh
cat>Job_${7}.sh<<EOF
#!/bin/sh  
source /cvmfs/cms.cern.ch/cmsset_default.sh 
cd /cvmfs/cms.cern.ch/el9_amd64_gcc12/cms/cmssw/CMSSW_14_1_0_pre4/src
eval `scramv1 runtime -sh`
cd \${_CONDOR_SCRATCH_DIR}                                                                                                                                                                                                                   
./${1} ${2} ${3} ${4} ${5} ${6}                                                                                                                                                                                                             
EOF

chmod 775 Job_${7}.sh

cat>condor_${7}<<EOF
accounting_group=group_cms                                                                   
use_x509userproxy = True                                                                     
x509userproxy = /tmp/x509up_u41070
executable = ./Job_${7}.sh
notification         = never
whenToTransferOutput = On_Exit
shouldTransferFiles  = yes
getenv = true 
request_memory       = 1992
transfer_input_files = ${1}
output               = \$(Cluster)_\$(Process)_${7}.out 
error                = \$(Cluster)_\$(Process)_${7}.err
Log                  = \$(Cluster)_\$(Process)_${7}.log
Queue
EOF

condor_submit condor_${7}
