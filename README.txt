Simulation for Windmill Si detectors at implantation sight, based
on Geant4 B4c example.
J. G. Cubiss (james.cubiss@york.ac.uk), 2022.

Simulation fires alpha particles, with the options to be followed
by cascadses of gamma rays and/or conversion electrons, and/or any
x rays (x ray energies and intensities taken from Firestone tables). 

The alpha particles are fired with energies and intensities
defined in first two columns of a_in.dat file. The columns after
that indicate the intensity of gamma transitions that follow
that alpha, with column number relating to line number of the
transition in g_in.dat.

In g_in.dat, the first column is the gamma energy, after that it
goes conversion electron energy, conversion electron coefficient,
alternating, as defined by the BrIcc calculator. Currently user
has to manually use BrIcc and input values to the g_in.dat file.

The user can use the config file to define the positions of the
annular and circular detector from the implantation foil (which 
is at the centre of the world). They can also select the Z of 
the alpha-decay daughter nucleus in order to automatically gen-
erate the correct x ray energies and intensities. User can also
select whether to include different decays (gammas, CEs, x rays)
using the options here.

User is recommended to compile the code somewhere, then use the
geant4_windmill_template directory and its contents to run the
simulation. The user is also advised to add the compiled simul-
ation directory to their PATH in .profile or .bashrc, for
instance:
if [ -d "$HOME/sfw/Windmill/Windmill_simulation/build" ] ; then
     PATH="$HOME/sfw/Windmill/Windmill_simulation/build:$PATH"
fi

User can then use the quickRun bash script which calls the exe-
cutable and runs the first argument as number of decays, for ins-
ance:
./quickRun 1000000
will run 1 million events, randomised over the intensities as
defined by user, and output data to a .root file. The energy
spectra will automatically be broadened into a crystal ball
function using the broaden.c macro. To stop this, please edit
the quickRun code accordingly.

Note that the code runs in multithreaded mode, currently set
to 4 threads. The verbosity is also set to minimise terminal
output and minimise run time.
