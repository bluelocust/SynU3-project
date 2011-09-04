*This contains an explanation of the code/files used for the SynU3 
project.
*This computational framework makes use of MEME Suite/FIMO, which 
requires a MEME file and a FASTA file.  Right now, these files need
to be plugged in at the MEME Suite submission website, but later
this manual process will be replaced with a web query or a local 
copy of MEME Suite software.
*Steps to visualize experimental sequences:
	1) Use make_fasta.py to process SEQ files
	2) Use FIMO and copy .txt output
	3) Use parse_MEME on FIMO output
	4) Use vis_a_vis files to visualize SynU3 sequences

FOLDER: Motifs

*all_motifs.py
	~Generates in silico SynU3 sequences using a Gillespie 
	algorithm
	~Converts strings of motifs in ACGT's and exports the 
	sequences in a FASTA format

*make_meme.py
	~Creates a MEME file for use in the MEME Suite of motif-
	finding tools
	~Which motifs do you want to look for?

***Ignore other folders.  Random .py files.***

FOLDER: vis_a_vis

*layers.py
	~Model in MVC architecture
	~Turns MEME output into matrices of motif data

*my_wx.py
	~View and Controller in MVC architecture
	~Uses wxpython and matplotlib to visualize sequences

FOLDER: synU3 clones

*make_fasta.py
	~turns SEQ files into FASTA file for MEME Suite

*parse_MEME.py
	~Simplifies output of MEME Suite for visualization 
	~Produces way-too-many intermediate .txt files that I used for
	debugging