Prepare ROOT Tree

Modify ID.C and ID.h to include photon ID selections and input variables

Compile the macro:

./rootcom ID all
./all sample_location outputrootfile -1 10000 evenweight #from submit.sh


Prepare 3 sets of output ROOT files:
        GJets
        Flat pt samples
        Background samples



Merge outputs:
        hadd Merge.root gjet_40*.root
        hadd Flat_all.root flat*.root
        hadd QCD_all.root qcd*.root

-------------------------------
Compute Effective Area (EA)


cd EAS
vi ContourBuilderpv.C   # set input = "Merge.root"
root -l -b runCont.C
        # Output: EA_1_60_TruePV_test.txt

-------------------------------
Compute pt Scaling

cd Isopt
vi CutID.h   # set input = "Flat_all.root"

root -l run_ParseEAs.C
        # Output: Mergedrun3barrel.root

cd scaling
root -l -b Ex_hcal.C     # input: Mergedrun3barrel.root
root -l -b Ex_ecal.C     # input: Mergedrun3barrel.root
        # Output: scale_para.txt  (pt scaling recorded)



--------------------------------
TMVA running


Repeat run_ParseEAs.C with Mergy.root (GJet sample):

root -l run_ParseEAs.C


Go to TMVA directory:

cd  Trainner/TMVA

Edit input scaling values in:
        Scr99.C
        Reader_neu.C
        Reg.C



Run TMVA workflow:

root -l  -b Scr99.C
root -l  -b Reg.C
root -l  -b Reader_neu.C

Output: mycustom.txt (contains SEF)

Used working points:

Tight (SEF = 0.73)
Medium (SEF = 0.82)
Loose (SEF = 0.89)


------------------------------------
plotting

cd Trainner/TMVA/plots
root -l aPlotRunner.C
# Output: SEF*_Barrel_Plots_updat_SR_TruePVID.root

root -l plotter.C   # use above root files as input


------------------------------------
Background plotting

Rerun run_ParseEAs.C with:
input = QCD_all.root
