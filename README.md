==================================================================
The GARNET package 
==================================================================
A tool designed to collect TF-DNA predictions from clustered position weight matrices using chromatin-accessibility data. 

Contact: Sara JC Gosline [sgosline@mit.edu]


Copyright (c) 2014 Sara JC Gosline


Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==================================================================
System Requirements:
==================================================================
1-Python 2.6.5 or higher: http://www.python.org

2-NetworkX (if using network option)

3-Python packages: numpy, scipy,matplotlib

==================================================================
Getting significant transcription factors from expression data
==================================================================
1-Get expression data in a tab-delimited format (gene\tfold change)

2-Collect chromatin accessibility data from related tissue/cell line as BED file

3-Download FASTA format for bed regions from http://usegalaxy.org

4-Enter all files into configuration file, run GARNET to determine transcription factors
