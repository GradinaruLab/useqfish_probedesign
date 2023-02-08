# Probe Designer for uSeqFISH and HCR

This tool designs probes for a given gene. 

We are currently working on deploying an accessible web app for probe design, but in the mean time you can use it through the source code.

## Installation
Download `probe_design.py`, `requirements.txt`, `db` (available at https://doi.org/10.22002/b791z-fzd10), and the example file for the method you want (either `HCR3_probe_example.py` or `useqFISH_probe_example.py`). 

Open the terminal and cd to the directory with your files. Make sure you have all the requirements installed by running `pip install -r requirements.txt`.
Next, open the example files and edit the inputs to your needs. Finally, run the example file. This should create an excel file containing your probe designs.
