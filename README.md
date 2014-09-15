rhd2arf
=========
A Matlab utility to convert *.rhd files to an arf file.

Setup
-------
* Clone the repository into your Matlab environment eg: `cd ~/MATLAB; git clone https://github.com/kylerbrown/rhd2arf.git` .

* Add the folder rhd2arf to your Matlab path by navigating to the folder in Matlab's Current Folder window. Right click on rhd2arf and select Add Folder to Path


Usage
-------

* In Matlab, navigate to a folder containing ONLY rhd files from a single, continuous session and enter `rhd2arf` in the Matlab terminal.


Channel Naming Convention
------------------------------

Rhd2arf allows the channel name to contain arf metadata as follows (ctrl+R to rename channels):

`NAME_DATATYPE_CHANNEL`

where NAME is the name of the channel, DATATYPE is the arf datatype (see the [arf spec](https://github.com/melizalab/arf/blob/master/specification.org)) and CHANNEL is the channel number.

For example: `HP01_23_7`, will have the name HP01 the datatype 23 and the channel number 7.



