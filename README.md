# Create-Stimuli-and-Background-in-Matlab
This repository contains a series of Matlab scripts written by Prof. John A. Endler and modified by Cedric van den Berg

Instructions (by John A. Endler):

1.  Make new pattern outline or modify existing pattern (easier) to get desired relative areas.

2.       Save new outline as a png file with a unique name.  You can make several and
         save them all to the same directory.

3.       You will use the new programme   SetUPOneNudibranchPatternFromPNGfile.m
         Make sure that this programme is in the same directory as your png outline files.

4.  Run SetUPOneNudibranchPatternFromPNGfile.m

     A.  Select a PNG file from the list
     B.  A little window will pop up asking about the number of colour classes.
     C.  Enter the number of colours in the pattern as a single number, for example,
         for 2 colours enter 2, 3 colours 3, etc. Then click "ok".
     D.  Now the tricky bit.  You need to assign a colour code to each outline. 
         The outermost loop will be 1, if there are 2 colours the next one will
         be 2, and if there are 3 colours some will be 2 and some will be 3.
         You may have to experiment to "get the hang of it".  Enter a single
         number for each loop.  These are just codes for colour classes and
         the programme you already have will actually assing colours to them.
         If there are 3 colours (patch within a patch) then the body will be 
         coded 1, the larger patches will be 2, and the patches within larger
         patches will be 3.  They first appear as thin lines, and when you
         have entered their codes, the lines will become thicker and have 
         different colours for easy code checking.  Watch out because the 
         patch boundaries go from left to right, without regard to whether 
         they are patches within patches, etc. So, for example in the pattern
         "spots in spots" loops 2 (larger) and 3 (spots within spots) will 
         have to be enered as the programme goes from left to right, so somtimes
         a 2 will follow a 3 and vice versa.  You may have to do the same
         pattern several times until you "get the hang of it".  Each time
         you enter a code in the box the line with thicken and that code
         number will appear by that line.
    E.   After you have entered the code for the last patch the programme saves
         the outlines as a .mat  (matlab data file).  This will be read by the
         programme, as before.

5.       Use MakeSaveANudiPatternWithColours.m to add the colours to use in the programme.
         It will save it as another .mat file for reading by the other programmes.
        


