### Instructions to install python, pylada before using pylada_defects.py

1. Make a copy of your current $HOME/.bashrc
2. Remove all module calls from $HOME/.bashrc
3. Remove all additions to the PYTHONPATH from $HOME/.bashrc
4. Logout and logback in, do not source $HOME/.bashrc, as you need your PATH variable to be reset

#### installing python using conda 

5. Purge all modules currently loaded

```module purge```

6. If conda module exists on your HPC use

```module load conda/XXXX```

otherwise install conda in your $HOME and add conda path in your .bashrc

```export PATH="$HOME/anaconda2/bin:$PATH"```

7. Create the virtual environemnt in python 2 and install pip, ipython, git and cmake

```conda create -n nameofenv python=2 pip ipython cmake git```

8. Start the environment

```source activate nameofenv```

#### installing pylada

9. Install pylada using pip

```pip install -I git+https://github.com/pylada/pylada-light```

alternatively install in developers mode

```pip install -e git+https://github.com/pylada/pylada-light#egg=pylada```

Steps 5-9 should have pylada installed on your system!!

10. Add following function in your $HOME/.bashrc if you plan to use pylada environment 

```pyladaenv() {
   module purge
   "load-modules-that-you-may-need"
   source activate nameofenv
   }```

#### before using pylada_defects.py

11. Activate your pylada environment

12. Put the folder /pyspglib in somewhere equivalent

```$HOME/anaconda2/envs/venv/lib/python2.7/site-packages/pyspglib```

13. Install package tess

```pip install tess```

14. Put the dot_pylada_* (*-depending on your HPC scheduler SLURM or MOAB/TORQUE)
in your home directory as (.pylada)

```$HOME/.pylada```

15. Edit the .pylada and change **YOUR_USER_NAME** in line 21 with your actual username

16. Put ipython_config.py in your home

```$HOME/.ipython/profile_default/ipython_config.py```

17. Create somewhere a folder that will serve to put your custom_chains
and add it to the PYTHONPATH environment variable in the .bashrc

```export PYTHONPATH=$PYTHONPATH:/some/folder/somewhere```

18. **You are ready to run the examples and additional examples in pylada/tutorial**
