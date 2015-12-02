%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%======================= FOLDER STRUCTURE =======================%

                              ~
                            /   \
                         data    data2  ----------
                        /       /         \        \
   ProductionRuns_Storage    rundata   pulse-post  ProductionRuns
       /   |   \            /   |   \                /   |   \
      F_xx_xxx_xxx         F_xx_xxx_xxx             F_xx_xxx_xxx

%----------------------------------------------------------------%

Simulations are run in ~/data2/ProductionRuns/F_xx_xxx_xxx, where F_xx_xxx_xxx is the simulation ID.
Simulation output is in binary.

Post-processing for each simulation is done in ~/data2/rundata/F_xx_xxx_xxx
Post-processing includes binary-->ascii conversion, matlab calculations and figure creation.
Figures and text files containing calculation output are saved in ~/data2/rundata/F_xx_xxx_xxx

The suite of matlab post-processing scripts and queue submission scripts are in ~/data2/pulse-post

To manage storage and disk space, files that are no longer needed are zipped and stored in ~/data/ProductionRuns_Storage/F_xx_xxx_xxx

NOTE: This folder structure is particular to Taryn Black and is shown here for context.
The commands below use this folder structure as an example for clarity.
Always check path names and modify to fit your structure as appropriate.
It is recommended that you use absolute paths (~/full/path/) rather than relative paths (../partial/path/).

%================================================================%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%=================== POST-PROCESSING WORKFLOW ===================%

1. MAKE A POST-PROCESSING DIRECTORY
Create new post-processing directory and move binary output files to this directory to be converted to ascii and processed.
Files moved: EP_G -- T_G -- U_G -- V_G -- W_G -- RO_G -- ROP_S1 -- ROP_S2 -- ROP_S3 -- mfixconst -- mfix.dat

FILE: MoveRenameRun.sh
LOCATION: ~/data2/pulse-post

EDITS:
L3: Define your simulation ID, which should also be the name of the folder it was run in --- e.g.: runnum=('F_10P_996_9999')
L4: Define the path where the simulation was run --- e.g.: sourcedir=~/data2/ProductionRuns/$runnum/
L5: Define the path where the simulation will be post-processed --- e.g.: destdir=~/data2/rundata/$runnum/
%%% NOTE: Leave /$runnum/ in paths - the script fills it based on your defined runnum

>> ./MoveRenameRun.sh

%----------------------------------------------------------------%

2. CALCULATE CHOKED VELOCITY FOR MINIMUM GAS FRACTION
This value is needed for the binary-->ascii conversion.

FILE: PulseCalculations.xlsx
LOCATION: PULSE_SHARE -- shared OneDrive folder.

EDITS:
B4-B7: Change the steady/min/max gas volume fraction and pulse frequency to match the simulation you are processing.
(Boxes with thick red outline)

OUTPUT: choked velocity for minimum gas fraction (box with thick green outline - F13)
TAKE NOTE OF THIS VALUE

%----------------------------------------------------------------%

3. CONVERT BINARY FILES TO ASCII
The conversion files are stored in /pulse-post for tracking but need to be run in the post-processing directories for individual simulations.
The Matlab post-processing scripts read ascii files.
Convert the binary files down into one .txt file per timestep.
Files have multiple columns representing different variables.
Each column contains the entire 3D grid for that variable at that timestep.
Conversion scripts c. Mary Benage

FILES: convert_Pulse.f90; Pulse_Convert.sh
LOCATION: ~/data2/pulse-post

Copy convert_Pulse.f90 and Pulse_Convert.sh into the directory to be processed
>> cp convert_Pulse.f90 ~/data2/rundata/F_xx_xxx_xxx
>> cp Pulse_Convert.f90 ~/data2/rundata/F_xx_xxx_xxx
>> cd ~/data2/rundata/F_xx_xxx_xxx

EDITS:
@convert_Pulse.f90
L57: Change initial_vel to the choked velocity for minimum gas fraction (from Step 2, in PulseCalculations.xlsx).
L77: Set the number of timesteps to be converted (e.g. if simulation is 300 seconds written at 5-second intervals starting at t=0, then timesteps=61).
@Pulse_Convert.sh
L24: Change path to your current path (e.g. ~/data2/rundata/F_xx_xxx_xxx).

>> ifort convert_Pulse.f90 -convert big_endian -o converter.exe
>> qsub Pulse_Convert.sh

STRUCTURE OF FILES WRITTEN OUT THAT ARE USED IN MATLAB:

EP_t**.txt
1      2      3      4      5      6      7
EP_G   EP_S1  EP_S2  EP_S3  X      Y      Z

U_G_t**.txt
1      2      3      4      5      6
U      V      W      X      Y      Z

T_G_t**.txt
1      2      3      4
T_G    X      Y      Z

Current_Density_t**.txt
1      2      3                         4                    5      6      7
RO_C   RO_G   RO_C/avg_RO_air(height)   avg_RO_air(height)   X      Y      Z

%%% NOTE: the Matlab post-processing script does not use all of the files that are written out, or all of the columns in each file (e.g. X/Y/Z). However, these are used in OpenDX post-processing.

%----------------------------------------------------------------%

4. ZIP BINARY FILES (OPTIONAL)
After converting to ascii, zip binary files in rundata to reduce disk usage.
Also zip binary files that were not converted and are not currently used.
After zipping is complete, delete the binary files from their directories.
Move zip files and simulation scripts to a storage directory.

>> cd ~/data2/rundata/F_xx_xxx_xxx
>> zip F_xx_xxx_xxx_converted.zip *_G *_S1 *_S2 *_S3
>> rm *_G *_S1 *_S2 *_S3
>> mv F_xx_xxx_xxx_converted.zip ~/data/ProductionRuns_Storage/F_xx_xxx_xxx
>> cd ~/data2/ProductionRuns/F_xx_xxx_xxx
>> zip F_xx_xxx_xxx_unused.zip *_G *_S1 *_S2 *_S3 *_G2 *_FLUX*
>> rm *_G *_S1 *_S2 *_S3 *_G2 *_FLUX*
>> mv * ~/data/ProductionRuns_Storage/F_xx_xxx_xxx

%----------------------------------------------------------------%

5. MATLAB POST-PROCESSING
Process simulation files and create movies and figures of the plume boundary, entrainment, particle concentration, gas temperature, density, and buoyancy.

FILES: pulse_post3D.m; submit_post.sh
LOCATION: ~/data2/pulse-post

MAJOR EDITS:
@pulse_post3D.m 
L22: Define the ID of the simulation(s) to be processed --- e.g. allruns = {'F_1_996_9999'}; or allruns = {'F_1_996_9999';'F_1_997_999'};
L26: Define path where .txt files (from Step 3) are stored and where post-processing output will be saved --- e.g. runpath = '~/data2/rundata/';
L31: Define path where Matlab post-processing scripts are stored --- e.g. postpath = '~/data2/pulse-post';
L38: Set maximum simulation time (in seconds) to be processed and displayed in movies. Can be shorter than the full simulation.
%%% NOTE: these are the key changes that need to be made (and probably only L22 will change every time). However, the user can change several other parameters, mostly to modify figure appearances. See detailed description below.
@submit_post.sh
L4: Set walltime for submission. I recommend ~40 hours for a 300 second simulation (160x160x200 grid). So if you're processing three at once (three IDs in allruns), walltime=120:00:00.
L9: Define path where Matlab post-processing scripts are stored.

>> msub submit_post.sh

%================================================================%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
