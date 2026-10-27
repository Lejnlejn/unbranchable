### Obstructing veering triangulations

This repository contains the code used to support the analysis and results presented in my paper:

    Three-manifolds that do not admit veering triangulations

    Layne Hall

    October 2026

Of course, if any errors are found, please report them to me at layne.hall [at] warwick.ac.uk.
Other comments, suggestions and requests are welcome. 

Tool declaration: Publicly available LLM tools were used on occasion for debugging, documentation of
various external modules, and suggestions for speeding up written processes.

### unbranchable_tris:
The source for the computations made to assert several theorems and examples throughout the paper. Namely:
- The search unbranch_search.py is the script run to produce the examples for Theorem 5.1 and Theorem 5.2.
- In verify_unbranch_data.ipynb, we give a manual way to verify the census unbranch data without having to rereproduce it from the above.
- In mobius_surgery.ipynb, we demonstrate the example used in the proof of Theorem 5.7. 

### drawing_total_order:
The source code that implements the idea described in Remark 4.19 and to produce the resulting Figure 16.
This includes a more detailed explanation of the promised procedure.

### Requirements
This code is highly dependant on existing three-manifold software. We require current versions of:
- Python 3.9, Regina, SnapPy, Veering