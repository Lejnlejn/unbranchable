# Code using for building census of essential triangulations that do not
# admit branchings.

import regina, snappy

from tqdm import tqdm
from modules.branch_dfs import look_for_no_branching, has_z2_merge_structure_dfs

# We can get most of the examples without the protracted search:
output=[]
for M in tqdm(snappy.OrientableCuspedCensus()):
    tri=regina.Triangulation3(M)
    verdict = has_z2_merge_structure_dfs(tri,transverse=False)
    if verdict:
        # then verdict output is a sig that witnesses the non-branching
        output.append((str(M.name()), verdict)) 
        
with open("no_branch_up_to_10.txt", "w") as f:
    for name, sig in output:
        f.write(f"({name!r}, {sig!r})\n")

# To get the full census with more examples, we do random sweeping, though it takes much longer:
# output=[]
# for M in tqdm(snappy.OrientableCuspedCensus()):
#     tri=regina.Triangulation3(M)
#     verdict = look_for_no_branching(tri,transverse=False)
#     if verdict:
#         # then verdict output is a sig that witnesses the non-branching
#         output.append((str(M.name()), verdict))

# Note that since we invoke randomness, the results may be more or less than we found on our run.
# Still, they should likely only be off by a small factor. See the notebook verify_unbranch_data
# for a verification that our censuses donot admit branchings. This acts then acts as a proof in
# spite of the variation in output of the above. 

# The cooriented branchings, use the option transverse=True

# The  11 tetrahedra, one must import Shana Li's census from a local installation.
# Since they do not currently have their snappy name inside snappy itself (as of September 2026), our
# 11 tetrahedron census records two pieces of data:
# - the census triangulation sig
# - the sig for which we found an unbranching triangulation.