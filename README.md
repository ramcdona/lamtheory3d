This is a WIP implementation of Chou72 to calculate 3D composite laminate properties.  I first encountered this method in Bogetti95.  My attempts to implement Chou's theory as described by Bogetti revealed that Bogetti has several errors in the formulas in the paper.  I since obtained a copy of Chou and am primarially working from that reference.

Bogetti (pdf p. 50) includes some test cases I am using for comparison.  These are the same examples as from MIL-HDBK-17 (pdf p. 629).

Although the first example matches reasonably well, the second and third examples are very wrong.  The first example is a [0, 90] laminate, while the other two use plys of more diverse orientation.  This leads me to suspect the error is in the transformation matrices.  Unfortunately, Chou does not include these in his paper and I am relying on the equations from Bogetti that have already proven to have numerous errors. 

Obviously, the OpenVSP implementation of this is in C++.  However, I have converted the routines to Matlab and Python (with the help of AI) to make the test cases stand-alone and easy to work with.  All three implementations give the same (wrong) results.


![Chou PC, Carleone J, Hsu CM. "Elastic Constants of Layered Media". Journal of Composite Materials. 1972;6(1):80-93. doi:10.1177/002199837200600107](papers/chou-et-al-1972-elastic-constants-of-layered-media.pdf)
![Bogetti, Travis A.; Hoppel, Christopher P.; Drysdale, William H. "Three-Dimensional Effective Property and Strength Prediction of Thick Laminated Composite Media.", 1995](papers/ADA302236.pdf)
![MIL-HDBK-17-3F](papers/MIL-HDBK-17-3F.pdf)