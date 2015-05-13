-----------------------------------------
 Django FlexDx Xpert Scale Up Web Version
-----------------------------------------

#FlexDx Xpert Scale Up

This README doc should help anybody interested in running the FlexDx Tuberculosis
transmission model using the web front-end. 

This model will run on any BSD based system with the GNU/Unix 'at' scheduler
command (OSX, Linux, etc.).  As long as the system runs python it should work.
You will need to install several python libraries, I recommend pip for that.
The file pip freeze in the root directory contains a list, but installing Django, 
matplotlib, numpy, scipy and Marmir should install all the other dependencies.

You should probably use a system like virtualenv to install the libraries
locally but it will work perfectly fine if you install the libraries globally.

This distribution does not contain the excel files of intermediate values generated per country for the country pre-sets.

#Installation

Make sure you have all the required python libraries Django, matplotlib, numpy
and scipy.  If you're having difficulty as an OSX user I've used a 'Superpack'
that contains matplotlib, numpy and scipy.  Search for it online.

Clone the repository and move into it ('git clone' or applicable); copy the 
file:

xpert_su/xpert_config.py_example
to
xpert_su/xpert_config.py

If you do not set the FLEXDX_HOME environment variable the logs from the
homebrew model will be put into the current working directory.  When the
homebrew is run from the web-model the logfiles will go into the repo
home directory (same directory as the manage.py file).  If you 
set FLEXDX_HOME then the logs will go to FLEXDX_HOME. 

then type:

`$ python manage.py runserver`  

`Performing system checks...`  

`System check identified no issues (0 silenced).`  
`April 29, 2015 - 01:31:16`  
`Django version 1.8, using settings 'xpert_su.settings'`  
`Starting development server at http://127.0.0.1:8000/`  
`Quit the server with CONTROL-C.`  

Then naviate your web-browser to http://localhost:8000 or 
http://127.0.0.1:8000/.  The model interface should be running there
now.

This model is CPU intensive and will require significantly more time
to execute on older machines.

Thanks,

The FlexDx Team
